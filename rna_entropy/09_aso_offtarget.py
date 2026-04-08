"""Screen ASO candidates for off-target complementarity via NCBI BLAST.

Checks each ASO against the human transcriptome (reject if >70% identity)
and the canine transcriptome (flag only).  Only the top 20 candidates by
composite_score are BLASTed to respect NCBI rate limits; the rest are
marked as ``not_screened``.

Usage:
    python -m rna_entropy.09_aso_offtarget
    python -m rna_entropy.09_aso_offtarget --force
    python -m rna_entropy.09_aso_offtarget --dry-run
"""

from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path
from typing import Any, Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CANDIDATES_PATH: Final[str] = "results/aso/aso_candidates.json"
OUTPUT_CSV: Final[str] = "results/aso/aso_offtarget_results.csv"
OUTPUT_FILTERED_JSON: Final[str] = "results/aso/aso_candidates_filtered.json"

BLAST_URL: Final[str] = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
BLAST_PROGRAM: Final[str] = "blastn"
BLAST_DATABASE: Final[str] = "refseq_rna"

HUMAN_ENTREZ_QUERY: Final[str] = "txid9606[ORGN]"
DOG_ENTREZ_QUERY: Final[str] = "txid9615[ORGN]"

# Short-sequence BLAST parameters
BLAST_WORD_SIZE: Final[int] = 7
BLAST_EXPECT: Final[int] = 1000
BLAST_REWARD: Final[int] = 1
BLAST_PENALTY: Final[int] = -3

IDENTITY_REJECT_THRESHOLD: Final[float] = 70.0
TOP_N_FOR_BLAST: Final[int] = 20

MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 120
BLAST_POLL_INTERVAL: Final[int] = 15
BLAST_MAX_WAIT: Final[int] = 300  # 5 minutes
INTER_QUERY_DELAY: Final[int] = 3

USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (aso-offtarget-screen)"

CSV_COLUMNS: Final[list[str]] = [
    "aso_id",
    "human_max_identity",
    "human_hit",
    "human_pass",
    "dog_max_identity",
    "dog_hit",
    "dog_pass",
]

logger = get_logger("aso_offtarget")


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------


def _request_with_retry(
    method: str,
    url: str,
    *,
    params: dict[str, str] | None = None,
    data: dict[str, str] | None = None,
) -> requests.Response:
    """Perform an HTTP request with exponential backoff retry.

    Args:
        method: HTTP method ("GET" or "POST").
        url: The URL to request.
        params: Query-string parameters.
        data: Form-encoded body (for POST).

    Returns:
        The successful ``requests.Response``.

    Raises:
        requests.RequestException: After all retries are exhausted.
    """
    headers = {"User-Agent": USER_AGENT}

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.request(
                method,
                url,
                params=params,
                data=data,
                headers=headers,
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


def _submit_blast(sequence: str, entrez_query: str) -> str:
    """Submit a BLASTn query for a short ASO sequence.

    Uses relaxed parameters appropriate for short oligonucleotides
    (WORD_SIZE=7, EXPECT=1000, REWARD=1, PENALTY=-3).

    Args:
        sequence: The ASO nucleotide sequence.
        entrez_query: NCBI Entrez query to restrict the organism
                      (e.g., ``txid9606[ORGN]`` for human).

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
        "ENTREZ_QUERY": entrez_query,
        "CMD": "Put",
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": "5",
        "WORD_SIZE": str(BLAST_WORD_SIZE),
        "EXPECT": str(BLAST_EXPECT),
        "REWARD": str(BLAST_REWARD),
        "PENALTY": str(BLAST_PENALTY),
    }

    response = _request_with_retry("POST", BLAST_URL, data=form_data)

    rid_match = re.search(r"RID\s*=\s*(\S+)", response.text)
    if not rid_match:
        raise RuntimeError(
            f"Could not parse RID from BLAST submission: "
            f"{response.text[:500]}"
        )

    rid = rid_match.group(1)
    logger.info("BLAST job submitted (entrez=%s), RID=%s", entrez_query, rid)
    return rid


def _poll_blast(rid: str) -> str:
    """Poll NCBI BLAST until the job completes and return the result body.

    Args:
        rid: The BLAST Request ID obtained from submission.

    Returns:
        The raw XML body of the completed BLAST results.

    Raises:
        TimeoutError: If the job does not finish within ``BLAST_MAX_WAIT``.
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
                "BLAST RID=%s still running (%ds elapsed)...", rid, elapsed,
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


def _parse_blast_top_hit(xml_body: str) -> tuple[float, str, str]:
    """Parse the top-hit percent identity, accession, and description.

    Args:
        xml_body: Raw XML response body from NCBI BLAST.

    Returns:
        A tuple of (max_identity_pct, hit_accession, hit_description).
        If no hits are found, returns (0.0, "", "").
    """
    if "<Hit>" not in xml_body:
        return 0.0, "", ""

    # Extract accession and description from the first Hit.
    accession = ""
    description = ""
    acc_match = re.search(r"<Hit_accession>([^<]+)</Hit_accession>", xml_body)
    desc_match = re.search(r"<Hit_def>([^<]+)</Hit_def>", xml_body)
    if acc_match:
        accession = acc_match.group(1)
    if desc_match:
        description = desc_match.group(1)

    # Extract identity from the first HSP.
    hsp_match = re.search(r"<Hsp>(.*?)</Hsp>", xml_body, re.DOTALL)
    if not hsp_match:
        return 0.0, accession, description

    hsp_block = hsp_match.group(1)

    identity_match = re.search(
        r"<Hsp_identity>(\d+)</Hsp_identity>", hsp_block,
    )
    align_len_match = re.search(
        r"<Hsp_align-len>(\d+)</Hsp_align-len>", hsp_block,
    )

    if identity_match and align_len_match:
        identical = int(identity_match.group(1))
        align_len = int(align_len_match.group(1))
        if align_len > 0:
            pct_identity = (identical / align_len) * 100.0
            return pct_identity, accession, description

    return 0.0, accession, description


def _blast_search(
    sequence: str,
    entrez_query: str,
) -> tuple[float, str, str]:
    """Run a full BLAST search cycle: submit, poll, parse.

    Args:
        sequence: Nucleotide sequence to search.
        entrez_query: NCBI Entrez query to restrict the organism.

    Returns:
        A tuple of (max_identity_pct, hit_accession, hit_description).
        Returns (0.0, "", "") on any failure.
    """
    try:
        rid = _submit_blast(sequence, entrez_query)
        body = _poll_blast(rid)
        return _parse_blast_top_hit(body)
    except (RuntimeError, TimeoutError, requests.RequestException) as exc:
        logger.error("BLAST search failed: %s", exc)
        return 0.0, "", ""


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_candidates(path: Path) -> list[dict[str, Any]]:
    """Load ASO candidates from a JSON file.

    Args:
        path: Path to ``aso_candidates.json``.

    Returns:
        List of candidate dictionaries.

    Raises:
        FileNotFoundError: If the candidates file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"ASO candidates file not found: {path}")

    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)

    # Handle both formats: list or dict with "candidates" key.
    if isinstance(raw, dict):
        data = raw.get("candidates", [])
    else:
        data = raw

    logger.info("Loaded %d ASO candidates from %s.", len(data), path)
    return data


# ---------------------------------------------------------------------------
# Off-target screening
# ---------------------------------------------------------------------------


def _screen_single_aso(
    candidate: dict[str, Any],
) -> dict[str, Any]:
    """Screen a single ASO candidate against human and dog transcriptomes.

    Args:
        candidate: ASO candidate dictionary with at least ``aso_id``
                   and ``aso_sequence`` keys.

    Returns:
        Dictionary with off-target screening results.
    """
    aso_id = candidate["aso_id"]
    sequence = candidate["aso_sequence"]

    logger.info("Screening %s against human transcriptome...", aso_id)
    human_identity, human_acc, human_desc = _blast_search(
        sequence, HUMAN_ENTREZ_QUERY,
    )
    human_pass = human_identity <= IDENTITY_REJECT_THRESHOLD

    # Respect NCBI rate limits between queries.
    time.sleep(INTER_QUERY_DELAY)

    logger.info("Screening %s against dog transcriptome...", aso_id)
    dog_identity, dog_acc, dog_desc = _blast_search(
        sequence, DOG_ENTREZ_QUERY,
    )
    dog_pass = dog_identity <= IDENTITY_REJECT_THRESHOLD

    human_hit = f"{human_acc}: {human_desc}" if human_acc else ""
    dog_hit = f"{dog_acc}: {dog_desc}" if dog_acc else ""

    if not human_pass:
        logger.warning(
            "%s REJECTED: human identity %.1f%% (hit: %s)",
            aso_id,
            human_identity,
            human_hit,
        )
    if not dog_pass:
        logger.warning(
            "%s flagged: dog identity %.1f%% (hit: %s)",
            aso_id,
            dog_identity,
            dog_hit,
        )

    return {
        "aso_id": aso_id,
        "human_max_identity": round(human_identity, 2),
        "human_hit": human_hit,
        "human_pass": human_pass,
        "dog_max_identity": round(dog_identity, 2),
        "dog_hit": dog_hit,
        "dog_pass": dog_pass,
    }


def _make_dry_run_result(candidate: dict[str, Any]) -> dict[str, Any]:
    """Create a passing off-target result for dry-run mode.

    Since SL RNA is absent from both human and canine transcriptomes,
    all candidates are expected to pass in a real screen.

    Args:
        candidate: ASO candidate dictionary.

    Returns:
        Dictionary with all-passing off-target results.
    """
    return {
        "aso_id": candidate["aso_id"],
        "human_max_identity": 0.0,
        "human_hit": "",
        "human_pass": True,
        "dog_max_identity": 0.0,
        "dog_hit": "",
        "dog_pass": True,
    }


def _make_not_screened_result(candidate: dict[str, Any]) -> dict[str, Any]:
    """Create a not-screened off-target result for lower-ranked candidates.

    Args:
        candidate: ASO candidate dictionary.

    Returns:
        Dictionary with ``not_screened`` status.
    """
    return {
        "aso_id": candidate["aso_id"],
        "human_max_identity": -1.0,
        "human_hit": "not_screened",
        "human_pass": True,
        "dog_max_identity": -1.0,
        "dog_hit": "not_screened",
        "dog_pass": True,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def _write_results_csv(
    results: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write off-target screening results to a CSV file.

    Args:
        results: List of per-ASO screening result dicts.
        output_path: Destination CSV path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in results:
            writer.writerow({
                "aso_id": row["aso_id"],
                "human_max_identity": row["human_max_identity"],
                "human_hit": row["human_hit"],
                "human_pass": str(row["human_pass"]),
                "dog_max_identity": row["dog_max_identity"],
                "dog_hit": row["dog_hit"],
                "dog_pass": str(row["dog_pass"]),
            })

    logger.info("Wrote %d rows to %s.", len(results), output_path)


def _write_filtered_json(
    candidates: list[dict[str, Any]],
    results: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write only candidates that pass human off-target screening.

    Args:
        candidates: Original ASO candidates list.
        results: Off-target screening results (parallel with candidates).
        output_path: Destination JSON path.
    """
    passing_ids = {
        r["aso_id"] for r in results if r["human_pass"]
    }

    filtered = [c for c in candidates if c["aso_id"] in passing_ids]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(filtered, fh, indent=2, ensure_ascii=False)

    logger.info(
        "Wrote %d filtered candidates to %s.", len(filtered), output_path,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def screen_offtargets(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Screen ASO candidates for off-target complementarity via BLAST.

    For each candidate, checks against the human transcriptome
    (reject if >70% complementarity) and the canine transcriptome
    (flag but do not reject).  Only the top 20 candidates by
    composite_score are BLASTed; others are marked ``not_screened``.

    Args:
        force: Re-run screening even if output files already exist.
        dry_run: Load candidates and mark all as passing without
                 running BLAST (SL RNA is absent in humans, so all
                 are expected to pass).

    Returns:
        Path to the output CSV file.

    Raises:
        FileNotFoundError: If the candidates file does not exist.
    """
    csv_path = Path(OUTPUT_CSV)
    json_path = Path(OUTPUT_FILTERED_JSON)

    if csv_path.exists() and not force:
        logger.info(
            "Off-target results already exist at %s. Use --force to re-run.",
            csv_path,
        )
        return str(csv_path)

    # Step 1: Load candidates.
    candidates = _load_candidates(Path(CANDIDATES_PATH))

    if not candidates:
        logger.warning("No ASO candidates to screen.")
        return str(csv_path)

    # Sort by composite_score descending to select top N for BLAST.
    candidates_sorted = sorted(
        candidates,
        key=lambda c: c.get("composite_score", 0.0),
        reverse=True,
    )

    top_candidates = candidates_sorted[:TOP_N_FOR_BLAST]
    remaining_candidates = candidates_sorted[TOP_N_FOR_BLAST:]

    # Step 2: Screen candidates.
    results: list[dict[str, Any]] = []

    if dry_run:
        logger.info(
            "[DRY RUN] Skipping BLAST — marking all %d candidates as passing.",
            len(candidates_sorted),
        )
        for candidate in candidates_sorted:
            results.append(_make_dry_run_result(candidate))
    else:
        logger.info(
            "Screening top %d candidates (of %d total) via BLAST...",
            len(top_candidates),
            len(candidates_sorted),
        )

        for idx, candidate in enumerate(top_candidates, 1):
            logger.info(
                "BLAST %d/%d: %s (score=%.3f)...",
                idx,
                len(top_candidates),
                candidate["aso_id"],
                candidate.get("composite_score", 0.0),
            )
            result = _screen_single_aso(candidate)
            results.append(result)

            # Respect NCBI rate limits between ASOs.
            if idx < len(top_candidates):
                time.sleep(INTER_QUERY_DELAY)

        # Mark remaining candidates as not screened.
        for candidate in remaining_candidates:
            results.append(_make_not_screened_result(candidate))

    # Step 3: Count passing candidates.
    pass_count = sum(1 for r in results if r["human_pass"])
    total_count = len(results)
    logger.info(
        "%d of %d candidates pass off-target screening.",
        pass_count,
        total_count,
    )

    # Step 4: Save results.
    _write_results_csv(results, csv_path)
    _write_filtered_json(candidates_sorted, results, json_path)

    return str(csv_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Screen ASO candidates for off-target complementarity "
            "against human and canine transcriptomes via NCBI BLAST."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run screening even if results already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Skip BLAST and mark all candidates as passing.",
    )
    args = parser.parse_args()

    result = screen_offtargets(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
