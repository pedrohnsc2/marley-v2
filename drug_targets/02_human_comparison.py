"""Compare *L. infantum* drug target enzymes with human homologs via NCBI BLAST.

For each enzyme in the ``drug_targets`` table, runs a BLASTp search against the
human proteome (taxid 9606) to measure sequence identity.  Enzymes with low
identity (< 60 %) also get their active-site residues compared via the UniProt
REST API.  High-identity enzymes (> 85 %) are rejected as unsuitable targets.

Results are saved to ``results/human_comparison.csv`` and persisted in the
Supabase ``drug_targets`` table.

Usage:
    python -m drug_targets.02_human_comparison
    python -m drug_targets.02_human_comparison --force
    python -m drug_targets.02_human_comparison --dry-run
"""

from __future__ import annotations

import csv
import re
import time
from pathlib import Path
from typing import Final

import requests

from core.db import get_all_drug_targets, update_drug_target
from core.logger import get_logger
from core.models import DrugTarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BLAST_URL: Final[str] = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
UNIPROT_BASE_URL: Final[str] = "https://rest.uniprot.org/uniprotkb"
LEISHMANIA_TAXID: Final[int] = 5671
HUMAN_TAXID: Final[int] = 9606

OUTPUT_DIR: Final[str] = "results"
OUTPUT_FILE: Final[str] = "results/human_comparison.csv"

MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 120

BLAST_POLL_INTERVAL: Final[int] = 15
BLAST_MAX_WAIT: Final[int] = 300  # 5 minutes

IDENTITY_REJECT_THRESHOLD: Final[float] = 0.85
IDENTITY_ACTIVE_SITE_THRESHOLD: Final[float] = 0.60
IDENTITY_PRIORITY_THRESHOLD: Final[float] = 0.40

USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (drug-target-discovery)"

CSV_COLUMNS: Final[list[str]] = [
    "gene_id",
    "gene_name",
    "pathway",
    "identity_score",
    "active_site_diff",
    "human_homolog_id",
    "classification",
]

logger = get_logger("human_comparison")


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


def _submit_blast(sequence: str) -> str:
    """Submit a BLASTp query against the human proteome.

    Args:
        sequence: Amino-acid sequence to search.

    Returns:
        The Request ID (RID) for polling results.

    Raises:
        RuntimeError: If the RID cannot be parsed from the response.
        requests.RequestException: On network failure.
    """
    form_data = {
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "QUERY": sequence,
        "ENTREZ_QUERY": "txid9606[ORGN]",
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
    # Initial wait recommended by NCBI before first poll.
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

        # If none of the status markers are present, the results are ready.
        return body

    raise TimeoutError(
        f"BLAST job {rid} did not complete within {BLAST_MAX_WAIT}s."
    )


def _parse_blast_result(xml_body: str) -> tuple[float, str]:
    """Parse top-hit identity and accession from BLAST XML output.

    Uses simple regex parsing to avoid requiring a full XML library for
    the small subset of fields we need.

    Args:
        xml_body: Raw XML response body from NCBI BLAST.

    Returns:
        A tuple of (percent_identity, hit_accession).  If no hits are found,
        returns (0.0, "").
    """
    # Check for "No hits found"
    if "<Hit>" not in xml_body:
        logger.info("No BLAST hits found.")
        return 0.0, ""

    # Extract the first Hit block.
    hit_match = re.search(r"<Hit>(.*?)</Hit>", xml_body, re.DOTALL)
    if not hit_match:
        return 0.0, ""

    hit_block = hit_match.group(1)

    # Extract accession from the first hit.
    acc_match = re.search(r"<Hit_accession>([^<]+)</Hit_accession>", hit_block)
    accession = acc_match.group(1) if acc_match else ""

    # Extract the first Hsp (highest-scoring segment pair) for identity.
    hsp_match = re.search(r"<Hsp>(.*?)</Hsp>", hit_block, re.DOTALL)
    if not hsp_match:
        return 0.0, accession

    hsp_block = hsp_match.group(1)

    identity_match = re.search(
        r"<Hsp_identity>(\d+)</Hsp_identity>", hsp_block
    )
    align_len_match = re.search(
        r"<Hsp_align-len>(\d+)</Hsp_align-len>", hsp_block
    )

    if identity_match and align_len_match:
        identical_residues = int(identity_match.group(1))
        alignment_length = int(align_len_match.group(1))
        if alignment_length > 0:
            percent_identity = (identical_residues / alignment_length) * 100.0
            return percent_identity, accession

    return 0.0, accession


def _run_blast_search(sequence: str) -> tuple[float, str]:
    """Run a full BLAST search cycle: submit, poll, parse.

    Args:
        sequence: Amino-acid sequence to search against human proteome.

    Returns:
        A tuple of (percent_identity, human_homolog_accession).
        Returns (0.0, "") on any failure.
    """
    try:
        rid = _submit_blast(sequence)
        body = _poll_blast(rid)
        return _parse_blast_result(body)
    except (RuntimeError, TimeoutError, requests.RequestException) as exc:
        logger.error("BLAST search failed: %s", exc)
        return 0.0, ""


# ---------------------------------------------------------------------------
# UniProt active-site helpers
# ---------------------------------------------------------------------------


def _fetch_active_sites_by_gene(
    gene_name: str, organism_id: int
) -> list[str]:
    """Fetch active-site residue annotations from UniProt for a gene.

    Args:
        gene_name: The gene or enzyme name to search for.
        organism_id: NCBI taxonomy ID for the organism.

    Returns:
        A sorted list of active-site annotation strings (e.g. "ACT_SITE 123").
        Returns an empty list if no annotations are found.
    """
    url = (
        f"{UNIPROT_BASE_URL}/search"
        f"?query=gene:{gene_name}+organism_id:{organism_id}"
        f"&fields=ft_act_site"
        f"&format=json"
        f"&size=1"
    )

    try:
        response = _request_with_retry("GET", url)
        data = response.json()
    except (requests.RequestException, ValueError) as exc:
        logger.warning(
            "Failed to fetch UniProt active sites for gene=%s org=%d: %s",
            gene_name,
            organism_id,
            exc,
        )
        return []

    return _extract_active_sites_from_json(data)


def _fetch_active_sites_by_accession(accession: str) -> list[str]:
    """Fetch active-site residue annotations from UniProt by accession ID.

    Args:
        accession: UniProt accession identifier.

    Returns:
        A sorted list of active-site annotation strings.
        Returns an empty list if no annotations are found.
    """
    url = f"{UNIPROT_BASE_URL}/{accession}?fields=ft_act_site&format=json"

    try:
        response = _request_with_retry("GET", url)
        data = response.json()
    except (requests.RequestException, ValueError) as exc:
        logger.warning(
            "Failed to fetch UniProt active sites for accession=%s: %s",
            accession,
            exc,
        )
        return []

    # Single-entry response wraps differently from search.
    features = data.get("features", [])
    sites: list[str] = []
    for feat in features:
        if feat.get("type") == "Active site":
            location = feat.get("location", {})
            start = location.get("start", {}).get("value", "")
            description = feat.get("description", "")
            sites.append(f"{start}:{description}")

    return sorted(sites)


def _extract_active_sites_from_json(data: dict) -> list[str]:
    """Extract active-site strings from a UniProt search JSON response.

    Args:
        data: Parsed JSON from the UniProt search endpoint.

    Returns:
        A sorted list of active-site annotation strings.
    """
    results = data.get("results", [])
    if not results:
        return []

    entry = results[0]
    features = entry.get("features", [])
    sites: list[str] = []
    for feat in features:
        if feat.get("type") == "Active site":
            location = feat.get("location", {})
            start = location.get("start", {}).get("value", "")
            description = feat.get("description", "")
            sites.append(f"{start}:{description}")

    return sorted(sites)


def _compute_active_site_diff(
    parasite_sites: list[str], human_sites: list[str]
) -> float:
    """Calculate the fraction of differing catalytic residues.

    Compares two lists of active-site annotations.  The diff score is
    the fraction of sites present in either organism but not both,
    divided by the total number of unique sites.

    Args:
        parasite_sites: Active-site annotations from the parasite enzyme.
        human_sites: Active-site annotations from the human homolog.

    Returns:
        A float in [0.0, 1.0] where 1.0 means completely different.
        Returns 1.0 when both lists are empty (no data implies maximal
        assumed divergence).
    """
    if not parasite_sites and not human_sites:
        return 1.0

    parasite_set = set(parasite_sites)
    human_set = set(human_sites)

    total = len(parasite_set | human_set)
    if total == 0:
        return 1.0

    differing = len(parasite_set ^ human_set)
    return differing / total


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def _classify_target(identity_score: float) -> str:
    """Classify a drug target by its human homolog identity score.

    Args:
        identity_score: Sequence identity as a fraction in [0.0, 1.0].

    Returns:
        One of "priority target", "moderate target", or "low priority".
    """
    if identity_score < IDENTITY_PRIORITY_THRESHOLD:
        return "priority target"
    if identity_score <= IDENTITY_ACTIVE_SITE_THRESHOLD:
        return "moderate target"
    return "low priority"


# ---------------------------------------------------------------------------
# Per-enzyme processing
# ---------------------------------------------------------------------------


def _process_enzyme(target: DrugTarget, dry_run: bool) -> dict[str, str]:
    """Run human comparison for a single drug target enzyme.

    Args:
        target: The ``DrugTarget`` to compare.
        dry_run: If True, skip API calls and return placeholder data.

    Returns:
        A dict matching ``CSV_COLUMNS`` with the comparison results.
    """
    gene_id = target.gene_id
    gene_name = target.gene_name
    pathway = target.pathway

    if dry_run:
        logger.info(
            "[DRY RUN] Would BLAST %s (%s) against human proteome.",
            gene_id,
            gene_name,
        )
        return {
            "gene_id": gene_id,
            "gene_name": gene_name,
            "pathway": pathway,
            "identity_score": "",
            "active_site_diff": "",
            "human_homolog_id": "",
            "classification": "",
        }

    # Step 1: BLAST against human proteome.
    logger.info("BLASTing %s (%s) against human proteome...", gene_id, gene_name)
    percent_identity, human_accession = _run_blast_search(target.sequence)
    identity_score = percent_identity / 100.0

    classification = _classify_target(identity_score)
    logger.info(
        "%s: %.1f%% identity with human (%s) — %s",
        gene_name,
        percent_identity,
        human_accession or "no hit",
        classification,
    )

    # Step 2: Active-site comparison for divergent enzymes.
    active_site_diff = 0.0
    if identity_score < IDENTITY_ACTIVE_SITE_THRESHOLD and human_accession:
        logger.info(
            "Fetching active-site annotations for %s (identity < 60%%)...",
            gene_name,
        )
        parasite_sites = _fetch_active_sites_by_gene(gene_name, LEISHMANIA_TAXID)
        human_sites = _fetch_active_sites_by_accession(human_accession)
        active_site_diff = _compute_active_site_diff(parasite_sites, human_sites)
        logger.info(
            "%s active-site diff: %.2f (parasite=%d sites, human=%d sites)",
            gene_name,
            active_site_diff,
            len(parasite_sites),
            len(human_sites),
        )

    # Step 3: Update database.
    update_fields: dict = {
        "identity_score": identity_score,
        "human_homolog_id": human_accession,
        "active_site_diff": active_site_diff,
    }

    if identity_score > IDENTITY_REJECT_THRESHOLD:
        update_fields["status"] = "rejected"
        logger.warning(
            "%s rejected: %.1f%% identity exceeds %.0f%% threshold.",
            gene_name,
            percent_identity,
            IDENTITY_REJECT_THRESHOLD * 100,
        )

    try:
        update_drug_target(gene_id, update_fields)
    except Exception:
        logger.error("Failed to update %s in database — continuing.", gene_id)

    return {
        "gene_id": gene_id,
        "gene_name": gene_name,
        "pathway": pathway,
        "identity_score": f"{identity_score:.4f}",
        "active_site_diff": f"{active_site_diff:.4f}",
        "human_homolog_id": human_accession,
        "classification": classification,
    }


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def _write_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    """Write comparison results to a CSV file.

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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compare_with_human(force: bool = False, dry_run: bool = False) -> str:
    """Compare all drug target enzymes with their human homologs.

    Fetches every enzyme from the ``drug_targets`` table and runs a BLASTp
    search against the human proteome to measure sequence identity.  Enzymes
    with low identity also get an active-site residue comparison via UniProt.
    High-identity enzymes are rejected.

    Args:
        force: Re-run comparison even if the output CSV already exists.
        dry_run: Log what would be done without making API calls or
                 writing to the database.

    Returns:
        Path to the output CSV file.

    Raises:
        RuntimeError: If no drug targets are found in the database.
    """
    output_path = Path(OUTPUT_FILE)

    if output_path.exists() and not force:
        logger.info(
            "Comparison file already exists at %s. Use --force to re-run.",
            output_path,
        )
        return str(output_path)

    if dry_run:
        logger.info("[DRY RUN] Would compare drug targets with human proteome.")

    # Fetch all drug targets from database.
    targets = get_all_drug_targets()

    if not targets:
        msg = "No drug targets found in database. Run 01_fetch_enzymes first."
        logger.error(msg)
        raise RuntimeError(msg)

    logger.info("Loaded %d drug target(s) from database.", len(targets))

    # Process each enzyme.
    rows: list[dict[str, str]] = []
    priority_count = 0

    for target in targets:
        row = _process_enzyme(target, dry_run=dry_run)
        rows.append(row)

        if not dry_run and row["classification"] == "priority target":
            priority_count += 1

    # Write output CSV (even in dry-run, write the placeholder file).
    if not dry_run:
        _write_csv(rows, output_path)
    else:
        logger.info(
            "[DRY RUN] Would write results to %s.", output_path
        )

    # Summary.
    total = len(targets)
    if not dry_run:
        logger.info(
            "%d of %d enzymes have < %.0f%% identity (priority targets).",
            priority_count,
            total,
            IDENTITY_PRIORITY_THRESHOLD * 100,
        )
    else:
        logger.info(
            "[DRY RUN] Would process %d enzyme(s) total.", total
        )

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compare L. infantum drug targets with human homologs via BLAST.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run comparison even if results CSV already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would be done without making API calls.",
    )
    args = parser.parse_args()

    result = compare_with_human(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
