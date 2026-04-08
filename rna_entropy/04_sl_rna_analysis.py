"""Detect the trypanosomatid Spliced Leader (SL) RNA in *L. infantum* transcripts.

The SL RNA is a conserved 39-nucleotide sequence added to the 5' end of every
mRNA via *trans*-splicing in trypanosomatids.  Its presence in parasite
transcripts and absence in human mRNAs makes it a distinctive molecular marker
for RNA-level target identification.

This module:

1. Scans the *L. infantum* transcriptome for the SL sequence (allowing up to
   2 mismatches in the first 100 nt of each transcript).
2. Confirms SL absence in the human transcriptome via NCBI BLAST.
3. Calculates GC content and Shannon entropy of the SL itself.
4. Outputs per-transcript results and a summary JSON.

Usage:
    python -m rna_entropy.04_sl_rna_analysis
    python -m rna_entropy.04_sl_rna_analysis --force
    python -m rna_entropy.04_sl_rna_analysis --dry-run
"""

from __future__ import annotations

import csv
import json
import math
import re
import time
from collections import Counter
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SL_SEQUENCE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"
SL_LENGTH: Final[int] = len(SL_SEQUENCE)

TRANSCRIPTOME_FILE: Final[str] = "data/raw/transcriptome_linf.fasta"
OUTPUT_DIR: Final[str] = "results/rna"
OUTPUT_CSV: Final[str] = "results/rna/sl_rna_targets.csv"
OUTPUT_JSON: Final[str] = "results/rna/sl_rna_summary.json"

BLAST_URL: Final[str] = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
BLAST_PROGRAM: Final[str] = "blastn"
BLAST_DATABASE: Final[str] = "refseq_rna"
BLAST_ENTREZ_QUERY: Final[str] = "txid9606[ORGN]"
BLAST_EVALUE_THRESHOLD: Final[float] = 0.01

SEARCH_WINDOW: Final[int] = 100
MAX_MISMATCHES: Final[int] = 2

MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 120

BLAST_POLL_INTERVAL: Final[int] = 15
BLAST_MAX_WAIT: Final[int] = 300  # 5 minutes

USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (rna-entropy-analysis)"

CSV_COLUMNS: Final[list[str]] = [
    "gene_id",
    "gene_name",
    "has_sl_rna",
    "sl_position",
    "sl_mismatches",
]

logger = get_logger("sl_rna_analysis")


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------


def _request_with_retry(
    method: str,
    url: str,
    *,
    params: dict[str, str] | None = None,
    data: dict[str, str] | None = None,
    headers: dict[str, str] | None = None,
) -> requests.Response:
    """Perform an HTTP request with exponential backoff retry.

    Args:
        method: HTTP method ("GET" or "POST").
        url: The URL to request.
        params: Query-string parameters.
        data: Form-encoded body (for POST/PUT).
        headers: Extra HTTP headers.

    Returns:
        The successful ``requests.Response``.

    Raises:
        requests.RequestException: After all retries are exhausted.
    """
    merged_headers = {"User-Agent": USER_AGENT}
    if headers:
        merged_headers.update(headers)

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.request(
                method,
                url,
                params=params,
                data=data,
                headers=merged_headers,
                timeout=REQUEST_TIMEOUT,
            )
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            if attempt < MAX_RETRIES - 1:
                delay = RETRY_DELAYS[attempt]
                logger.warning(
                    "Request failed (attempt %d/%d): %s — retrying in %ds",
                    attempt + 1,
                    MAX_RETRIES,
                    exc,
                    delay,
                )
                time.sleep(delay)
            else:
                logger.error(
                    "All %d attempts failed for %s: %s",
                    MAX_RETRIES,
                    url,
                    exc,
                )
                raise

    raise requests.RequestException("Unexpected retry loop exit")


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def _parse_fasta(fasta_path: Path) -> list[tuple[str, str, str]]:
    """Parse a FASTA file into (gene_id, gene_name, sequence) tuples.

    The header line is split on whitespace; the first token becomes
    ``gene_id`` and the remainder becomes ``gene_name``.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        A list of (gene_id, gene_name, sequence) tuples.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    records: list[tuple[str, str, str]] = []
    current_id = ""
    current_name = ""
    current_seq: list[str] = []

    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                # Save previous record.
                if current_id:
                    records.append(
                        (current_id, current_name, "".join(current_seq))
                    )
                parts = line[1:].split(None, 1)
                current_id = parts[0] if parts else ""
                current_name = parts[1] if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line.upper())

    # Save last record.
    if current_id:
        records.append((current_id, current_name, "".join(current_seq)))

    return records


# ---------------------------------------------------------------------------
# SL RNA detection
# ---------------------------------------------------------------------------


def _count_mismatches(seq_a: str, seq_b: str) -> int:
    """Count positional mismatches between two equal-length strings.

    Args:
        seq_a: First nucleotide string.
        seq_b: Second nucleotide string.

    Returns:
        Number of positions where the characters differ.
    """
    return sum(a != b for a, b in zip(seq_a, seq_b))


def _find_sl_in_transcript(
    sequence: str,
) -> tuple[bool, int, int]:
    """Search for the SL sequence in the first SEARCH_WINDOW nucleotides.

    Slides the SL pattern across the search window and returns the
    best-matching position (fewest mismatches) if it falls within
    ``MAX_MISMATCHES``.

    Args:
        sequence: Full transcript nucleotide sequence.

    Returns:
        A tuple of (found, position, mismatches).  If not found,
        returns (False, -1, -1).
    """
    window = sequence[:SEARCH_WINDOW].upper()
    if len(window) < SL_LENGTH:
        return False, -1, -1

    best_pos = -1
    best_mismatches = SL_LENGTH  # worst case

    for start in range(len(window) - SL_LENGTH + 1):
        substr = window[start : start + SL_LENGTH]
        mismatches = _count_mismatches(SL_SEQUENCE, substr)
        if mismatches < best_mismatches:
            best_mismatches = mismatches
            best_pos = start

    if best_mismatches <= MAX_MISMATCHES:
        return True, best_pos, best_mismatches

    return False, -1, -1


# ---------------------------------------------------------------------------
# SL RNA properties
# ---------------------------------------------------------------------------


def _compute_gc_content(sequence: str) -> float:
    """Calculate GC content as a fraction in [0.0, 1.0].

    Args:
        sequence: Nucleotide sequence (case-insensitive).

    Returns:
        GC fraction.  Returns 0.0 for an empty sequence.
    """
    if not sequence:
        return 0.0
    upper = sequence.upper()
    gc_count = upper.count("G") + upper.count("C")
    return gc_count / len(upper)


def _compute_shannon_entropy(sequence: str) -> float:
    """Calculate Shannon entropy H(X) of a nucleotide sequence.

    Uses base-2 logarithm so the result is in bits.  For a perfectly
    uniform distribution over 4 bases the maximum is 2.0 bits; a
    highly conserved sequence will have lower entropy.

    Args:
        sequence: Nucleotide sequence (case-insensitive).

    Returns:
        Shannon entropy in bits.  Returns 0.0 for an empty sequence.
    """
    if not sequence:
        return 0.0

    upper = sequence.upper()
    length = len(upper)
    counts = Counter(upper)

    entropy = 0.0
    for count in counts.values():
        if count > 0:
            probability = count / length
            entropy -= probability * math.log2(probability)

    return entropy


# ---------------------------------------------------------------------------
# BLAST helpers
# ---------------------------------------------------------------------------


def _submit_blast_nucleotide(sequence: str) -> str:
    """Submit a BLASTn query against the human transcriptome.

    Args:
        sequence: Nucleotide sequence to search.

    Returns:
        The Request ID (RID) for polling results.

    Raises:
        RuntimeError: If the RID cannot be parsed from the response.
        requests.RequestException: On network failure.
    """
    form_data = {
        "PROGRAM": BLAST_PROGRAM,
        "DATABASE": BLAST_DATABASE,
        "QUERY": sequence,
        "ENTREZ_QUERY": BLAST_ENTREZ_QUERY,
        "CMD": "Put",
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": "5",
        "EXPECT": str(BLAST_EVALUE_THRESHOLD),
    }

    response = _request_with_retry("POST", BLAST_URL, data=form_data)

    rid_match = re.search(r"RID\s*=\s*(\S+)", response.text)
    if not rid_match:
        raise RuntimeError(
            f"Could not parse RID from BLAST submission response: "
            f"{response.text[:500]}"
        )

    rid = rid_match.group(1)
    logger.info("BLAST job submitted, RID=%s", rid)
    return rid


def _poll_blast(rid: str) -> str:
    """Poll NCBI BLAST until the job is ready and return the result body.

    Args:
        rid: The BLAST Request ID obtained from submission.

    Returns:
        The raw text/XML body of the completed BLAST results.

    Raises:
        TimeoutError: If the job does not complete within ``BLAST_MAX_WAIT`` seconds.
        RuntimeError: If the BLAST job reports a failure status.
        requests.RequestException: On network failure.
    """
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
    }

    elapsed = 0
    initial_wait = min(BLAST_POLL_INTERVAL, BLAST_MAX_WAIT)
    logger.info("Waiting %ds before first BLAST poll...", initial_wait)
    time.sleep(initial_wait)
    elapsed += initial_wait

    while elapsed < BLAST_MAX_WAIT:
        response = _request_with_retry("GET", BLAST_URL, params=params)
        body = response.text

        if "Status=WAITING" in body:
            logger.info(
                "BLAST RID=%s still running (%ds elapsed)...", rid, elapsed
            )
            time.sleep(BLAST_POLL_INTERVAL)
            elapsed += BLAST_POLL_INTERVAL
            continue

        if "Status=FAILED" in body:
            raise RuntimeError(f"BLAST job {rid} failed.")

        if "Status=UNKNOWN" in body:
            raise RuntimeError(f"BLAST job {rid} expired or is unknown.")

        return body

    raise TimeoutError(
        f"BLAST job {rid} did not complete within {BLAST_MAX_WAIT}s."
    )


def _parse_blast_evalue(xml_body: str) -> float | None:
    """Extract the best (lowest) e-value from BLAST XML output.

    Args:
        xml_body: Raw XML response body from NCBI BLAST.

    Returns:
        The best e-value as a float, or ``None`` if no hits were found.
    """
    if "<Hit>" not in xml_body:
        logger.info("No BLAST hits found.")
        return None

    evalue_match = re.search(r"<Hsp_evalue>([^<]+)</Hsp_evalue>", xml_body)
    if not evalue_match:
        return None

    return float(evalue_match.group(1))


def _check_sl_absent_in_human() -> bool:
    """BLAST the SL sequence against the human transcriptome.

    Returns:
        ``True`` if no significant hit is found (e-value > threshold),
        meaning the SL is absent in humans.  ``False`` if a significant
        hit is found.
    """
    try:
        rid = _submit_blast_nucleotide(SL_SEQUENCE)
        body = _poll_blast(rid)
        best_evalue = _parse_blast_evalue(body)
    except (RuntimeError, TimeoutError, requests.RequestException) as exc:
        logger.error("BLAST search for SL in human failed: %s", exc)
        # Conservative: cannot confirm absence.
        return False

    if best_evalue is None:
        logger.info("SL RNA: no hits in human transcriptome — confirmed absent.")
        return True

    if best_evalue > BLAST_EVALUE_THRESHOLD:
        logger.info(
            "SL RNA: best human hit e-value=%.2e (> %.2e) — confirmed absent.",
            best_evalue,
            BLAST_EVALUE_THRESHOLD,
        )
        return True

    logger.warning(
        "SL RNA: significant human hit found (e-value=%.2e).", best_evalue
    )
    return False


# ---------------------------------------------------------------------------
# CSV / JSON output
# ---------------------------------------------------------------------------


def _write_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    """Write SL RNA detection results to a CSV file.

    Args:
        rows: List of dicts with keys matching ``CSV_COLUMNS``.
        output_path: Path to the output CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Wrote %d row(s) to %s.", len(rows), output_path)


def _write_json(data: dict, output_path: Path) -> None:
    """Write a summary dict as pretty-printed JSON.

    Args:
        data: Dictionary to serialize.
        output_path: Path to the output JSON file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)

    logger.info("Wrote summary to %s.", output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyze_sl_rna(force: bool = False, dry_run: bool = False) -> str:
    """Detect the Spliced Leader RNA in *L. infantum* transcripts.

    Scans each transcript in the transcriptome FASTA for the conserved
    39 nt SL sequence (allowing up to 2 mismatches in the first 100 nt),
    then confirms its absence in humans via NCBI BLAST.

    Args:
        force: Re-run analysis even if the output CSV already exists.
        dry_run: Search for SL locally but skip the BLAST confirmation
                 against the human transcriptome.

    Returns:
        Path to the output CSV file.

    Raises:
        FileNotFoundError: If the transcriptome FASTA is missing.
    """
    csv_path = Path(OUTPUT_CSV)
    json_path = Path(OUTPUT_JSON)

    if csv_path.exists() and not force:
        logger.info(
            "SL RNA results already exist at %s. Use --force to re-run.",
            csv_path,
        )
        return str(csv_path)

    if dry_run:
        logger.info("[DRY RUN] Will search for SL RNA locally, skip BLAST.")

    # Step 1: Read transcripts.
    transcriptome_path = Path(TRANSCRIPTOME_FILE)
    records = _parse_fasta(transcriptome_path)
    total_transcripts = len(records)
    logger.info("Loaded %d transcripts from %s.", total_transcripts, transcriptome_path)

    # Step 2: Search for SL in each transcript.
    rows: list[dict[str, str]] = []
    sl_positive_count = 0

    for gene_id, gene_name, sequence in records:
        found, position, mismatches = _find_sl_in_transcript(sequence)
        if found:
            sl_positive_count += 1

        rows.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "has_sl_rna": str(found),
            "sl_position": str(position) if found else "",
            "sl_mismatches": str(mismatches) if found else "",
        })

    logger.info(
        "SL RNA detected in %d of %d transcripts.",
        sl_positive_count,
        total_transcripts,
    )

    # Step 3: BLAST SL against human transcriptome (skip in dry run).
    if dry_run:
        absent_in_human = True
        logger.info("[DRY RUN] Skipping BLAST — assuming SL absent in human.")
    else:
        logger.info("BLASTing SL sequence against human transcriptome...")
        absent_in_human = _check_sl_absent_in_human()

    # Step 4: Calculate SL RNA properties.
    sl_gc_content = _compute_gc_content(SL_SEQUENCE)
    sl_entropy = _compute_shannon_entropy(SL_SEQUENCE)

    logger.info(
        "SL RNA properties: GC=%.4f, Shannon entropy=%.4f bits.",
        sl_gc_content,
        sl_entropy,
    )

    # Step 5: Write outputs.
    _write_csv(rows, csv_path)

    summary = {
        "total_transcripts": total_transcripts,
        "sl_positive_count": sl_positive_count,
        "sl_gc_content": round(sl_gc_content, 6),
        "sl_entropy": round(sl_entropy, 6),
        "absent_in_human": absent_in_human,
        "sl_sequence": SL_SEQUENCE,
        "max_mismatches_allowed": MAX_MISMATCHES,
        "search_window_nt": SEARCH_WINDOW,
    }
    _write_json(summary, json_path)

    # Step 6: Final log.
    absent_label = "YES" if absent_in_human else "NO"
    logger.info(
        "SL RNA found in %d of %d transcripts. Absent in human: %s",
        sl_positive_count,
        total_transcripts,
        absent_label,
    )

    return str(csv_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Detect the trypanosomatid Spliced Leader RNA in "
            "L. infantum transcripts and confirm absence in humans."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run analysis even if results already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Search for SL locally but skip BLAST against human.",
    )
    args = parser.parse_args()

    result = analyze_sl_rna(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
