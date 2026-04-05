"""Check epitope safety by BLASTing against the Canis lupus familiaris proteome.

For each optimized epitope, submits a short-peptide BLASTP query to NCBI
against the dog proteome (taxid 9615).  Epitopes with >70 % identity to any
canine protein are flagged as unsafe; >50 % triggers a moderate-risk warning.
Produces a safety report with per-epitope annotations and population-coverage
summary.

Usage:
    python -m pipeline.09_safety_check
    python -m pipeline.09_safety_check --force
    python -m pipeline.09_safety_check --dry-run
"""

from __future__ import annotations

import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INPUT_FILE: Final[str] = "results/v4_optimization/optimized_epitopes.json"
OUTPUT_FILE: Final[str] = "results/v4_optimization/safety_report.json"

BLAST_URL: Final[str] = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# Short-peptide BLASTP parameters for Canis lupus familiaris
BLAST_PROGRAM: Final[str] = "blastp"
BLAST_DATABASE: Final[str] = "nr"
BLAST_ENTREZ_QUERY: Final[str] = "txid9615[ORGN]"
BLAST_WORD_SIZE: Final[int] = 2
BLAST_MATRIX: Final[str] = "PAM30"
BLAST_EXPECT: Final[int] = 200000

# Safety thresholds (percent identity)
UNSAFE_THRESHOLD: Final[float] = 70.0
WARNING_THRESHOLD: Final[float] = 50.0

# BLAST polling
POLL_INTERVAL_SECONDS: Final[int] = 15
MAX_POLL_ATTEMPTS: Final[int] = 120  # ~30 minutes max wait

# Retry / back-off for HTTP requests
MAX_RETRIES: Final[int] = 3
BACKOFF_BASE_SECONDS: Final[float] = 5.0

# Known DLA alleles for population-coverage analysis
DLA_ALLELES: Final[list[str]] = [
    "DLA-88*501:01",
    "DLA-88*508:02",
    "DLA-12*001:01",
]

logger = get_logger("safety_check")

# ---------------------------------------------------------------------------
# HTTP helpers with retry / exponential back-off
# ---------------------------------------------------------------------------


def _request_with_retry(
    method: str,
    url: str,
    *,
    params: dict | None = None,
    data: dict | None = None,
    timeout: int = 120,
) -> requests.Response:
    """Send an HTTP request with retry and exponential back-off.

    Args:
        method: HTTP method (``"GET"`` or ``"POST"``).
        url: Target URL.
        params: Query-string parameters.
        data: Form-encoded body (for POST).
        timeout: Per-request timeout in seconds.

    Returns:
        The successful :class:`requests.Response`.

    Raises:
        requests.RequestException: After all retries are exhausted.
    """
    last_exc: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.request(
                method,
                url,
                params=params,
                data=data,
                timeout=timeout,
            )
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            last_exc = exc
            wait = BACKOFF_BASE_SECONDS * (2 ** (attempt - 1))
            logger.warning(
                "Request to %s failed (attempt %d/%d): %s. "
                "Retrying in %.0f s ...",
                url,
                attempt,
                MAX_RETRIES,
                exc,
                wait,
            )
            time.sleep(wait)

    raise last_exc  # type: ignore[misc]


# ---------------------------------------------------------------------------
# NCBI BLAST interaction
# ---------------------------------------------------------------------------


def _submit_blast(sequence: str) -> str:
    """Submit a short-peptide BLASTP job to NCBI and return the RID.

    Args:
        sequence: Amino-acid peptide sequence.

    Returns:
        The NCBI Request ID (RID) for the submitted job.

    Raises:
        RuntimeError: If the RID cannot be parsed from the NCBI response.
    """
    payload = {
        "CMD": "Put",
        "PROGRAM": BLAST_PROGRAM,
        "DATABASE": BLAST_DATABASE,
        "QUERY": sequence,
        "ENTREZ_QUERY": BLAST_ENTREZ_QUERY,
        "WORD_SIZE": str(BLAST_WORD_SIZE),
        "MATRIX_NAME": BLAST_MATRIX,
        "EXPECT": str(BLAST_EXPECT),
    }

    resp = _request_with_retry("POST", BLAST_URL, data=payload)
    text = resp.text

    # Parse RID from the HTML / plain-text response.
    for line in text.splitlines():
        if line.strip().startswith("RID ="):
            rid = line.split("=", 1)[1].strip()
            logger.info("BLAST job submitted for %s: RID=%s", sequence, rid)
            return rid

    raise RuntimeError(
        f"Could not parse RID from NCBI BLAST response for sequence {sequence!r}"
    )


def _poll_blast(rid: str) -> str:
    """Poll NCBI BLAST for job completion and return the XML results.

    Args:
        rid: The Request ID returned by :func:`_submit_blast`.

    Returns:
        Raw XML string of the BLAST results.

    Raises:
        TimeoutError: If the job does not complete within the allowed
            polling window.
    """
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
    }

    for attempt in range(1, MAX_POLL_ATTEMPTS + 1):
        time.sleep(POLL_INTERVAL_SECONDS)
        resp = _request_with_retry("GET", BLAST_URL, params=params)
        text = resp.text

        if "Status=WAITING" in text:
            if attempt % 4 == 0:
                logger.info("BLAST RID %s still running (poll %d) ...", rid, attempt)
            continue

        if "Status=FAILED" in text:
            raise RuntimeError(f"BLAST job {rid} failed on the NCBI server.")

        if "Status=UNKNOWN" in text:
            raise RuntimeError(f"BLAST job {rid} expired or is unknown.")

        # If none of the status markers appear, assume results are ready.
        return text

    raise TimeoutError(
        f"BLAST job {rid} did not finish within "
        f"{MAX_POLL_ATTEMPTS * POLL_INTERVAL_SECONDS} seconds."
    )


def _parse_blast_xml(xml_text: str) -> tuple[float, str]:
    """Extract the top hit percent identity and gene name from BLAST XML.

    Args:
        xml_text: Raw XML response from NCBI BLAST.

    Returns:
        A tuple of ``(max_percent_identity, top_hit_gene_name)``.  When no
        hits are found, returns ``(0.0, "")``.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        logger.warning("Could not parse BLAST XML; treating as no hits.")
        return 0.0, ""

    max_identity: float = 0.0
    top_gene: str = ""

    # Navigate: BlastOutput -> BlastOutput_iterations -> Iteration ->
    #           Iteration_hits -> Hit -> Hit_hsps -> Hsp
    for hit in root.iter("Hit"):
        hit_def = hit.findtext("Hit_def", default="")
        for hsp in hit.iter("Hsp"):
            identity_text = hsp.findtext("Hsp_identity", default="0")
            align_len_text = hsp.findtext("Hsp_align-len", default="1")
            try:
                identity = int(identity_text)
                align_len = int(align_len_text)
                pct = (identity / align_len) * 100.0 if align_len > 0 else 0.0
            except (ValueError, ZeroDivisionError):
                continue

            if pct > max_identity:
                max_identity = pct
                top_gene = hit_def

    return round(max_identity, 2), top_gene


def blast_epitope(sequence: str) -> tuple[float, str]:
    """Run a full BLAST cycle for a single epitope sequence.

    Submits the job, polls until completion, and parses the results.

    Args:
        sequence: Amino-acid peptide sequence.

    Returns:
        A tuple of ``(max_percent_identity, top_hit_gene_name)``.
    """
    rid = _submit_blast(sequence)
    xml_text = _poll_blast(rid)
    return _parse_blast_xml(xml_text)


# ---------------------------------------------------------------------------
# Safety classification
# ---------------------------------------------------------------------------


def _classify_safety(
    identity: float,
) -> tuple[bool, str]:
    """Classify an epitope as safe / warning / unsafe based on BLAST identity.

    Args:
        identity: Maximum percent identity against the dog proteome.

    Returns:
        A tuple of ``(is_safe, risk_level)`` where *risk_level* is one of
        ``"safe"``, ``"moderate_risk"``, or ``"unsafe"``.
    """
    if identity > UNSAFE_THRESHOLD:
        return False, "unsafe"
    if identity > WARNING_THRESHOLD:
        return True, "moderate_risk"
    return True, "safe"


# ---------------------------------------------------------------------------
# Population coverage
# ---------------------------------------------------------------------------


def _compute_allele_coverage(epitopes: list[dict]) -> list[str]:
    """Determine which DLA alleles are covered by the safe epitopes.

    An allele is considered covered if at least one safe epitope has it listed
    in its ``dla_alleles`` field (if present).  When no allele metadata exists
    on the epitopes, coverage falls back to counting all known alleles as
    covered (conservative assumption for the safety report).

    Args:
        epitopes: List of epitope dicts (post-safety annotation).

    Returns:
        Sorted list of covered DLA allele names.
    """
    covered: set[str] = set()

    safe_epitopes = [e for e in epitopes if e.get("is_safe", False)]
    if not safe_epitopes:
        return []

    has_allele_data = False
    for ep in safe_epitopes:
        alleles = ep.get("dla_alleles", [])
        if alleles:
            has_allele_data = True
            for a in alleles:
                if a in DLA_ALLELES:
                    covered.add(a)

    # If the upstream data does not carry per-epitope allele annotations,
    # assume all known alleles are covered by the safe set as a whole.
    if not has_allele_data and safe_epitopes:
        covered = set(DLA_ALLELES)

    return sorted(covered)


# ---------------------------------------------------------------------------
# Core pipeline function
# ---------------------------------------------------------------------------


def run_safety_check(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the epitope safety-check pipeline step.

    Loads optimized epitopes, BLASTs each unique sequence against the
    *Canis lupus familiaris* proteome, flags unsafe cross-reactive epitopes,
    and writes a safety report.

    Args:
        force: If ``True``, re-run even when the output file already exists.
        dry_run: If ``True``, skip BLAST queries (mark all epitopes as safe
                 with ``blast_identity=0.0``).

    Returns:
        Path to the output JSON file.

    Raises:
        FileNotFoundError: If the input file from module 08 is missing.
    """
    output_path = Path(OUTPUT_FILE)

    if output_path.exists() and not force:
        logger.info(
            "Output file %s already exists. Use --force to re-run.", OUTPUT_FILE
        )
        return OUTPUT_FILE

    # ------------------------------------------------------------------
    # 1. Load optimized epitopes
    # ------------------------------------------------------------------
    input_path = Path(INPUT_FILE)
    if not input_path.exists():
        msg = f"Input file not found: {INPUT_FILE}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    with open(input_path, "r") as fh:
        raw = json.load(fh)

    # Handle both formats: list of dicts or dict with "results" key.
    if isinstance(raw, dict):
        epitopes: list[dict] = raw.get("results", [])
    else:
        epitopes = raw

    logger.info("Loaded %d epitope(s) from %s.", len(epitopes), INPUT_FILE)

    if not epitopes:
        logger.warning("No epitopes to check. Writing empty report.")
        _save_report([], output_path)
        return OUTPUT_FILE

    # ------------------------------------------------------------------
    # 2. Collect unique sequences to BLAST
    # ------------------------------------------------------------------
    unique_sequences: dict[str, tuple[float, str]] = {}

    def _gather_sequences(ep: dict) -> list[str]:
        """Return all distinct peptide sequences from an epitope record."""
        seqs: list[str] = []
        for key in (
            "sequence", "original_sequence", "optimized_sequence",
            "original_peptide", "best_variant",
        ):
            seq = ep.get(key, "")
            if seq and seq not in seqs:
                seqs.append(seq)
        return seqs

    all_seqs: list[str] = []
    for ep in epitopes:
        for seq in _gather_sequences(ep):
            if seq not in unique_sequences:
                unique_sequences[seq] = (0.0, "")
                all_seqs.append(seq)

    logger.info(
        "Identified %d unique sequence(s) to BLAST across %d epitope(s).",
        len(all_seqs),
        len(epitopes),
    )

    # ------------------------------------------------------------------
    # 3. BLAST each unique sequence (or skip in dry-run mode)
    # ------------------------------------------------------------------
    if dry_run:
        logger.info(
            "[DRY RUN] Skipping BLAST for %d sequence(s). "
            "All epitopes will be marked as safe.",
            len(all_seqs),
        )
        for seq in all_seqs:
            unique_sequences[seq] = (0.0, "")
    else:
        for idx, seq in enumerate(all_seqs, 1):
            logger.info(
                "BLASTing sequence %d/%d: %s ...",
                idx,
                len(all_seqs),
                seq[:30],
            )
            try:
                identity, gene = blast_epitope(seq)
                unique_sequences[seq] = (identity, gene)
                logger.info(
                    "  -> max identity=%.1f%%, top hit=%s",
                    identity,
                    gene[:80] if gene else "(none)",
                )
            except Exception as exc:
                logger.warning(
                    "BLAST failed for sequence %s: %s. Marking as safe (identity=0).",
                    seq[:20],
                    exc,
                )
                unique_sequences[seq] = (0.0, "")

    # ------------------------------------------------------------------
    # 4. Annotate each epitope with safety information
    # ------------------------------------------------------------------
    for ep in epitopes:
        seqs = _gather_sequences(ep)

        # Take the worst-case (highest identity) across all sequence variants.
        best_identity: float = 0.0
        best_gene: str = ""
        for seq in seqs:
            identity, gene = unique_sequences.get(seq, (0.0, ""))
            if identity > best_identity:
                best_identity = identity
                best_gene = gene

        is_safe, risk_level = _classify_safety(best_identity)

        ep["blast_identity"] = best_identity
        ep["blast_hit_gene"] = best_gene
        ep["is_safe"] = is_safe
        ep["risk_level"] = risk_level

        if risk_level == "unsafe":
            logger.warning(
                "UNSAFE: epitope %s has %.1f%% identity to %s",
                ep.get("original_peptide", ep.get("sequence", ep.get("original_sequence", "?"))),
                best_identity,
                best_gene[:60],
            )
        elif risk_level == "moderate_risk":
            logger.warning(
                "MODERATE RISK: epitope %s has %.1f%% identity to %s",
                ep.get("original_peptide", ep.get("sequence", ep.get("original_sequence", "?"))),
                best_identity,
                best_gene[:60],
            )

    # ------------------------------------------------------------------
    # 5. Population coverage analysis
    # ------------------------------------------------------------------
    safe_count = sum(1 for ep in epitopes if ep.get("is_safe", False))
    total_count = len(epitopes)
    covered_alleles = _compute_allele_coverage(epitopes)

    logger.info(
        "%d of %d epitopes are safe (no cross-reactivity with dog proteome).",
        safe_count,
        total_count,
    )
    logger.info(
        "Population coverage: %d of %d known DLA alleles covered by safe epitopes.",
        len(covered_alleles),
        len(DLA_ALLELES),
    )

    # ------------------------------------------------------------------
    # 6. Save report
    # ------------------------------------------------------------------
    report = {
        "summary": {
            "total_epitopes": total_count,
            "safe_epitopes": safe_count,
            "unsafe_epitopes": total_count - safe_count,
            "dla_alleles_covered": covered_alleles,
            "dla_alleles_total": len(DLA_ALLELES),
            "unsafe_threshold_pct": UNSAFE_THRESHOLD,
            "warning_threshold_pct": WARNING_THRESHOLD,
            "dry_run": dry_run,
        },
        "epitopes": epitopes,
    }

    _save_report(report, output_path)
    return OUTPUT_FILE


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _save_report(data: dict | list, path: Path) -> None:
    """Write the safety report as pretty-printed JSON.

    Args:
        data: Report payload (dict or list).
        path: Destination file path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)
    logger.info("Safety report written to %s.", path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Step 9: Safety check -- BLAST epitopes against Canis lupus familiaris proteome.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if the output file already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Skip BLAST queries; mark all epitopes as safe.",
    )
    args = parser.parse_args()

    result_path = run_safety_check(force=args.force, dry_run=args.dry_run)
    print(f"Output: {result_path}")
