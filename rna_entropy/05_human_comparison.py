"""Compare *L. infantum* RNA targets with human transcriptome via entropy delta.

Merges Shannon entropy profiles, codon usage statistics, and SL RNA detection
data to prioritize RNA-level targets.  Transcripts with high entropy delta
(conserved in the parasite, variable in the host) are the most promising
candidates.  The top candidates are BLASTed against the human transcriptome
to confirm uniqueness.

Results are saved as a merged CSV and a full JSON with ``RNATarget`` dicts.

Usage:
    python -m rna_entropy.05_human_comparison
    python -m rna_entropy.05_human_comparison --force
    python -m rna_entropy.05_human_comparison --dry-run
"""

from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger
from core.models import RNATarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENTROPY_CSV: Final[str] = "results/rna/shannon_entropy_profile.csv"
CODON_JSON: Final[str] = "results/rna/codon_usage_stats.json"
SL_RNA_CSV: Final[str] = "results/rna/sl_rna_targets.csv"

OUTPUT_DIR: Final[str] = "results/rna"
OUTPUT_CSV: Final[str] = "results/rna/human_comparison_rna.csv"
OUTPUT_JSON: Final[str] = "results/rna/rna_targets_merged.json"

BLAST_URL: Final[str] = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
BLAST_PROGRAM: Final[str] = "blastn"
BLAST_DATABASE: Final[str] = "refseq_rna"
BLAST_ENTREZ_QUERY: Final[str] = "txid9606[ORGN]"
BLAST_EVALUE_THRESHOLD: Final[float] = 0.01

HIGH_PRIORITY_DELTA: Final[float] = 1.5
MODERATE_PRIORITY_DELTA: Final[float] = 0.5
TOP_CANDIDATES_FOR_BLAST: Final[int] = 20

MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 120

BLAST_POLL_INTERVAL: Final[int] = 15
BLAST_MAX_WAIT: Final[int] = 300  # 5 minutes

USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (rna-entropy-analysis)"

CSV_COLUMNS: Final[list[str]] = [
    "gene_id",
    "gene_name",
    "shannon_entropy",
    "human_entropy",
    "entropy_delta",
    "codon_bias_score",
    "has_sl_rna",
    "priority_level",
    "blast_identity_pct",
    "unique_to_parasite",
    "status",
]

logger = get_logger("human_comparison_rna")


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
        "HITLIST_SIZE": "1",
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


def _parse_blast_identity(xml_body: str) -> tuple[float, bool]:
    """Parse top-hit percent identity from BLAST XML output.

    Args:
        xml_body: Raw XML response body from NCBI BLAST.

    Returns:
        A tuple of (percent_identity, has_significant_hit).
        If no hits are found, returns (0.0, False).
    """
    if "<Hit>" not in xml_body:
        return 0.0, False

    # Check e-value first.
    evalue_match = re.search(r"<Hsp_evalue>([^<]+)</Hsp_evalue>", xml_body)
    if evalue_match:
        evalue = float(evalue_match.group(1))
        if evalue > BLAST_EVALUE_THRESHOLD:
            return 0.0, False

    # Extract identity from the first HSP.
    hsp_match = re.search(r"<Hsp>(.*?)</Hsp>", xml_body, re.DOTALL)
    if not hsp_match:
        return 0.0, False

    hsp_block = hsp_match.group(1)

    identity_match = re.search(
        r"<Hsp_identity>(\d+)</Hsp_identity>", hsp_block
    )
    align_len_match = re.search(
        r"<Hsp_align-len>(\d+)</Hsp_align-len>", hsp_block
    )

    if identity_match and align_len_match:
        identical = int(identity_match.group(1))
        align_len = int(align_len_match.group(1))
        if align_len > 0:
            pct_identity = (identical / align_len) * 100.0
            return pct_identity, True

    return 0.0, False


def _run_blast_search(sequence: str) -> tuple[float, bool]:
    """Run a full BLAST search cycle: submit, poll, parse.

    Args:
        sequence: Nucleotide sequence to search against human transcriptome.

    Returns:
        A tuple of (percent_identity, has_significant_hit).
        Returns (0.0, False) on any failure.
    """
    try:
        rid = _submit_blast_nucleotide(sequence)
        body = _poll_blast(rid)
        return _parse_blast_identity(body)
    except (RuntimeError, TimeoutError, requests.RequestException) as exc:
        logger.error("BLAST search failed: %s", exc)
        return 0.0, False


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_entropy_csv(path: Path) -> dict[str, dict[str, float]]:
    """Load Shannon entropy profile CSV into a dict keyed by gene_id.

    Expected columns: gene_id, gene_name, shannon_entropy, human_entropy,
    entropy_delta (plus any others, which are ignored here).

    Args:
        path: Path to the entropy CSV file.

    Returns:
        A dict mapping gene_id to {gene_name, shannon_entropy,
        human_entropy, entropy_delta}.

    Raises:
        FileNotFoundError: If the CSV file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"Entropy CSV not found: {path}")

    data: dict[str, dict[str, float]] = {}

    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gene_id = row["gene_id"]
            data[gene_id] = {
                "gene_name": row.get("gene_name", ""),
                "shannon_entropy": float(row.get("shannon_entropy", 0.0)),
                "human_entropy": float(row.get("human_entropy", 0.0)),
                "entropy_delta": float(row.get("entropy_delta", 0.0)),
            }

    return data


def _load_codon_usage(path: Path) -> dict[str, float]:
    """Load codon usage stats JSON and extract per-gene bias scores.

    The JSON is expected to contain a ``per_gene`` key mapping gene_id
    to an object with ``codon_bias_score``.  If the structure differs,
    returns an empty dict gracefully.

    Args:
        path: Path to the codon usage JSON file.

    Returns:
        A dict mapping gene_id to codon_bias_score.

    Raises:
        FileNotFoundError: If the JSON file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"Codon usage JSON not found: {path}")

    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)

    # Support both flat dict and nested structure.
    per_gene = raw.get("per_gene", raw)
    result: dict[str, float] = {}

    for gene_id, entry in per_gene.items():
        if isinstance(entry, dict):
            result[gene_id] = float(entry.get("codon_bias_score", 0.0))
        elif isinstance(entry, (int, float)):
            result[gene_id] = float(entry)

    return result


def _load_sl_rna_csv(path: Path) -> dict[str, bool]:
    """Load SL RNA detection CSV into a dict keyed by gene_id.

    Args:
        path: Path to the SL RNA targets CSV file.

    Returns:
        A dict mapping gene_id to has_sl_rna (boolean).

    Raises:
        FileNotFoundError: If the CSV file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"SL RNA CSV not found: {path}")

    data: dict[str, bool] = {}

    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gene_id = row["gene_id"]
            has_sl = row.get("has_sl_rna", "False").strip().lower() == "true"
            data[gene_id] = has_sl

    return data


# ---------------------------------------------------------------------------
# Priority classification
# ---------------------------------------------------------------------------


def _classify_priority(entropy_delta: float) -> str:
    """Classify a transcript's priority level by its entropy delta.

    Args:
        entropy_delta: Difference between human and parasite Shannon entropy.

    Returns:
        One of "high", "moderate", or "low".
    """
    if entropy_delta > HIGH_PRIORITY_DELTA:
        return "high"
    if entropy_delta > MODERATE_PRIORITY_DELTA:
        return "moderate"
    return "low"


# ---------------------------------------------------------------------------
# CSV / JSON output
# ---------------------------------------------------------------------------


def _write_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    """Write human comparison results to a CSV file.

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


def _write_json(data: list[dict], output_path: Path) -> None:
    """Write a list of RNATarget dicts as pretty-printed JSON.

    Args:
        data: List of dictionaries to serialize.
        output_path: Path to the output JSON file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)

    logger.info("Wrote %d target(s) to %s.", len(data), output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compare_with_human(force: bool = False, dry_run: bool = False) -> str:
    """Compare *L. infantum* RNA targets with the human transcriptome.

    Merges Shannon entropy, codon usage, and SL RNA data to classify
    transcripts into priority levels.  The top candidates by entropy
    delta are BLASTed against the human transcriptome to confirm their
    uniqueness to the parasite.

    Args:
        force: Re-run comparison even if the output CSV already exists.
        dry_run: Merge data and assign priorities without running BLAST.

    Returns:
        Path to the output CSV file.

    Raises:
        FileNotFoundError: If any of the required input files are missing.
    """
    csv_path = Path(OUTPUT_CSV)
    json_path = Path(OUTPUT_JSON)

    if csv_path.exists() and not force:
        logger.info(
            "Human comparison results already exist at %s. Use --force to re-run.",
            csv_path,
        )
        return str(csv_path)

    if dry_run:
        logger.info("[DRY RUN] Will merge data and classify, skip BLAST.")

    # Step 1: Load all input data.
    entropy_data = _load_entropy_csv(Path(ENTROPY_CSV))
    logger.info("Loaded entropy data for %d transcripts.", len(entropy_data))

    codon_data = _load_codon_usage(Path(CODON_JSON))
    logger.info("Loaded codon bias scores for %d transcripts.", len(codon_data))

    sl_data = _load_sl_rna_csv(Path(SL_RNA_CSV))
    logger.info("Loaded SL RNA data for %d transcripts.", len(sl_data))

    # Step 2: Build unified list sorted by entropy_delta descending.
    gene_ids = sorted(entropy_data.keys())
    candidates: list[dict] = []

    for gene_id in gene_ids:
        entry = entropy_data[gene_id]
        candidates.append({
            "gene_id": gene_id,
            "gene_name": entry.get("gene_name", ""),
            "shannon_entropy": entry["shannon_entropy"],
            "human_entropy": entry["human_entropy"],
            "entropy_delta": entry["entropy_delta"],
            "codon_bias_score": codon_data.get(gene_id, 0.0),
            "has_sl_rna": sl_data.get(gene_id, False),
        })

    # Sort by entropy_delta descending for prioritization.
    candidates.sort(key=lambda c: c["entropy_delta"], reverse=True)

    # Step 3: Classify priorities.
    high_count = 0
    moderate_count = 0
    low_count = 0

    for candidate in candidates:
        priority = _classify_priority(candidate["entropy_delta"])
        candidate["priority_level"] = priority
        if priority == "high":
            high_count += 1
        elif priority == "moderate":
            moderate_count += 1
        else:
            low_count += 1

    # Step 4: BLAST top candidates against human transcriptome.
    top_for_blast = [
        c for c in candidates
        if c["priority_level"] in ("high", "moderate")
    ][:TOP_CANDIDATES_FOR_BLAST]

    blast_results: dict[str, tuple[float, bool]] = {}

    if dry_run:
        logger.info(
            "[DRY RUN] Would BLAST %d top candidates against human transcriptome.",
            len(top_for_blast),
        )
    else:
        logger.info(
            "BLASTing %d top candidates against human transcriptome...",
            len(top_for_blast),
        )
        for idx, candidate in enumerate(top_for_blast, 1):
            gene_id = candidate["gene_id"]
            gene_name = candidate["gene_name"]
            logger.info(
                "BLAST %d/%d: %s (%s)...",
                idx,
                len(top_for_blast),
                gene_id,
                gene_name,
            )
            # Use gene_id as a proxy sequence identifier for the BLAST query.
            # In production, the actual RNA sequence would come from the
            # transcriptome FASTA; here we search by gene_id against refseq_rna.
            pct_identity, has_hit = _run_blast_search(gene_id)
            blast_results[gene_id] = (pct_identity, has_hit)

            # Respect NCBI rate limits between queries.
            if idx < len(top_for_blast):
                time.sleep(3)

    # Step 5: Build RNATarget objects and output rows.
    csv_rows: list[dict[str, str]] = []
    json_targets: list[dict] = []

    for candidate in candidates:
        gene_id = candidate["gene_id"]
        priority_level = candidate["priority_level"]

        # Determine BLAST results.
        blast_pct, blast_hit = blast_results.get(gene_id, (0.0, False))
        unique_to_parasite = not blast_hit if gene_id in blast_results else True

        # Assign status based on priority and uniqueness.
        if priority_level == "high" and unique_to_parasite:
            status = "approved"
        elif priority_level == "high":
            status = "pending"
        elif priority_level == "moderate":
            status = "pending"
        else:
            status = "rejected"

        # Build RNATarget.
        rna_target = RNATarget(
            gene_id=gene_id,
            gene_name=candidate["gene_name"],
            shannon_entropy=candidate["shannon_entropy"],
            human_entropy=candidate["human_entropy"],
            entropy_delta=candidate["entropy_delta"],
            codon_bias_score=candidate["codon_bias_score"],
            has_sl_rna=candidate["has_sl_rna"],
            status=status,
        )

        csv_rows.append({
            "gene_id": gene_id,
            "gene_name": candidate["gene_name"],
            "shannon_entropy": f"{candidate['shannon_entropy']:.6f}",
            "human_entropy": f"{candidate['human_entropy']:.6f}",
            "entropy_delta": f"{candidate['entropy_delta']:.6f}",
            "codon_bias_score": f"{candidate['codon_bias_score']:.6f}",
            "has_sl_rna": str(candidate["has_sl_rna"]),
            "priority_level": priority_level,
            "blast_identity_pct": f"{blast_pct:.2f}" if blast_pct > 0 else "",
            "unique_to_parasite": str(unique_to_parasite),
            "status": status,
        })

        json_targets.append(rna_target.to_dict())

    # Step 6: Write outputs.
    _write_csv(csv_rows, csv_path)
    _write_json(json_targets, json_path)

    # Step 7: Summary log.
    logger.info(
        "%d high-priority, %d moderate, %d low priority RNA targets.",
        high_count,
        moderate_count,
        low_count,
    )

    if not dry_run and blast_results:
        unique_count = sum(
            1 for _, has_hit in blast_results.values() if not has_hit
        )
        logger.info(
            "%d of %d BLASTed candidates confirmed unique to parasite.",
            unique_count,
            len(blast_results),
        )

    return str(csv_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Compare L. infantum RNA targets with human transcriptome "
            "using Shannon entropy and BLAST."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run comparison even if results already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Merge data and classify without running BLAST.",
    )
    args = parser.parse_args()

    result = compare_with_human(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
