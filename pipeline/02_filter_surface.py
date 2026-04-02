"""Filter proteins for surface-exposed and secreted candidates using signal peptide predictions.

Uses the SignalP 6.0 web API to classify signal peptides, retaining only
proteins with a classic Sec/SPI signal peptide -- strong indicators of
surface exposure or secretion, both desirable traits for vaccine targets.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import requests
from Bio import SeqIO
from tqdm import tqdm

from core.db import upsert_candidate
from core.logger import get_logger
from core.models import Candidate

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SIGNALP_URL: str = "https://services.healthtech.dtu.dk/services/SignalP-6.0/"
INPUT_FILE: str = "data/raw/proteins.fasta"
OUTPUT_FILE: str = "data/raw/surface_proteins.fasta"
BATCH_SIZE: int = 50
POLL_INTERVAL: int = 10  # seconds
POLL_TIMEOUT: int = 300  # 5 minutes
SIGNAL_TYPE: str = "Sec/SPI"  # classic signal peptide

logger = get_logger("filter_surface")

# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def parse_fasta(fasta_path: str) -> list[tuple[str, str, str]]:
    """Parse a FASTA file into a list of (gene_id, gene_name, sequence) tuples.

    The gene_id is taken from the record id, the gene_name from the record
    description (falling back to the id when no description is available),
    and the sequence is the full amino-acid string.

    Args:
        fasta_path: Path to a FASTA file.

    Returns:
        A list of ``(gene_id, gene_name, sequence)`` tuples.

    Raises:
        FileNotFoundError: If *fasta_path* does not exist.
    """
    path = Path(fasta_path)
    if not path.exists():
        logger.error("FASTA file not found: %s", fasta_path)
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    results: list[tuple[str, str, str]] = []
    for record in SeqIO.parse(str(path), "fasta"):
        gene_id = record.id
        # Description often looks like "gene_id Some Gene Name OS=...";
        # strip the leading id to get the human-readable portion.
        description = record.description.replace(gene_id, "", 1).strip()
        gene_name = description if description else gene_id
        sequence = str(record.seq)
        results.append((gene_id, gene_name, sequence))

    logger.info("Parsed %d sequences from %s.", len(results), fasta_path)
    return results


# ---------------------------------------------------------------------------
# SignalP 6.0 API interaction
# ---------------------------------------------------------------------------


def _build_fasta_payload(sequences: list[tuple[str, str]]) -> str:
    """Convert a list of (header, sequence) pairs into FASTA-formatted text.

    Args:
        sequences: Each element is ``(header, sequence)``.

    Returns:
        A multi-record FASTA string.
    """
    lines: list[str] = []
    for header, seq in sequences:
        lines.append(f">{header}")
        lines.append(seq)
    return "\n".join(lines)


def submit_batch(
    sequences: list[tuple[str, str]], organism: str = "eukarya"
) -> dict:
    """Submit a batch of sequences to the SignalP 6.0 web API.

    The API expects multipart form data with a FASTA-formatted input and the
    target organism group.

    Args:
        sequences: List of ``(header, sequence)`` pairs.
        organism: Organism group for prediction. One of ``"eukarya"``,
                  ``"gram-"``, ``"gram+"``, ``"archaea"``.

    Returns:
        The parsed JSON response from the submission endpoint, which
        typically contains a ``job_id`` key.

    Raises:
        requests.HTTPError: On any non-2xx HTTP status.
    """
    fasta_text = _build_fasta_payload(sequences)

    form_data = {
        "fastaentry": fasta_text,
        "organism": organism,
        "format": "json",
    }

    try:
        response = requests.post(
            SIGNALP_URL,
            data=form_data,
            timeout=60,
        )
        response.raise_for_status()
        logger.info(
            "Submitted batch of %d sequences to SignalP.", len(sequences)
        )
        return response.json()
    except requests.ConnectionError:
        logger.error(
            "Cannot reach SignalP API at %s. The service may be unavailable.",
            SIGNALP_URL,
        )
        raise
    except requests.HTTPError as exc:
        logger.error(
            "SignalP API returned HTTP %s: %s",
            exc.response.status_code if exc.response is not None else "N/A",
            exc.response.text[:500] if exc.response is not None else "",
        )
        raise
    except requests.RequestException:
        logger.exception("Unexpected error submitting batch to SignalP.")
        raise


def poll_results(job_id: str) -> dict | None:
    """Poll the SignalP 6.0 API until results are ready or timeout is reached.

    Args:
        job_id: The job identifier returned by :func:`submit_batch`.

    Returns:
        The parsed JSON results, or ``None`` if the job did not complete
        within :data:`POLL_TIMEOUT` seconds.
    """
    results_url = f"{SIGNALP_URL.rstrip('/')}/results/{job_id}"
    elapsed = 0

    while elapsed < POLL_TIMEOUT:
        try:
            response = requests.get(results_url, timeout=30)
            if response.status_code == 200:
                data = response.json()
                # Some APIs return a "status" field while processing.
                if data.get("status") in ("complete", "done", None):
                    logger.info(
                        "Results ready for job %s after %ds.", job_id, elapsed
                    )
                    return data
            elif response.status_code == 202:
                # 202 Accepted -- still processing.
                pass
            else:
                logger.warning(
                    "Unexpected status %d while polling job %s.",
                    response.status_code,
                    job_id,
                )
        except requests.RequestException:
            logger.warning(
                "Network error polling job %s; will retry.", job_id
            )

        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL

    logger.error(
        "Timed out after %ds waiting for SignalP job %s.", POLL_TIMEOUT, job_id
    )
    return None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_signal_hits(results: dict) -> set[str]:
    """Extract gene_ids classified as Sec/SPI from a SignalP results dict.

    The function handles two common response shapes:
      1. ``{"results": [{"id": ..., "prediction": ...}, ...]}``
      2. ``{"<gene_id>": {"prediction": ...}, ...}``

    Args:
        results: Parsed JSON returned by :func:`poll_results`.

    Returns:
        A set of gene_id strings that have a Sec/SPI prediction.
    """
    hits: set[str] = set()

    # Shape 1: list-based results.
    if "results" in results and isinstance(results["results"], list):
        for entry in results["results"]:
            prediction = entry.get("prediction", "")
            if SIGNAL_TYPE in prediction:
                hits.add(entry.get("id", entry.get("gene_id", "")))
        return hits

    # Shape 2: dict keyed by gene_id.
    for gene_id, info in results.items():
        if isinstance(info, dict):
            prediction = info.get("prediction", "")
            if SIGNAL_TYPE in prediction:
                hits.add(gene_id)

    return hits


# ---------------------------------------------------------------------------
# Main filter function
# ---------------------------------------------------------------------------


def filter_surface_proteins(force: bool = False) -> str:
    """Run the surface-protein filter using SignalP 6.0.

    Reads protein sequences from :data:`INPUT_FILE`, submits them in batches
    to SignalP, retains those classified as *Sec/SPI*, persists each passing
    candidate to Supabase, and writes the filtered FASTA to :data:`OUTPUT_FILE`.

    Args:
        force: When ``False`` and *OUTPUT_FILE* already exists, skip
               processing and return the existing path.

    Returns:
        The path to the filtered FASTA file (:data:`OUTPUT_FILE`).
    """
    output_path = Path(OUTPUT_FILE)

    if not force and output_path.exists():
        logger.info(
            "Output file %s already exists; skipping (use force=True to "
            "rerun).",
            OUTPUT_FILE,
        )
        return OUTPUT_FILE

    # ------------------------------------------------------------------
    # 1. Parse input
    # ------------------------------------------------------------------
    proteins = parse_fasta(INPUT_FILE)
    if not proteins:
        logger.warning("No proteins found in %s. Nothing to filter.", INPUT_FILE)
        return OUTPUT_FILE

    total = len(proteins)
    logger.info("Starting surface filter for %d proteins.", total)

    # ------------------------------------------------------------------
    # 2. Process in batches via SignalP
    # ------------------------------------------------------------------
    surface_ids: set[str] = set()
    batches = [
        proteins[i : i + BATCH_SIZE] for i in range(0, total, BATCH_SIZE)
    ]

    for batch in tqdm(batches, desc="SignalP batches", unit="batch"):
        batch_seqs = [(gene_id, seq) for gene_id, _, seq in batch]
        try:
            submission = submit_batch(batch_seqs)
        except requests.RequestException:
            logger.error(
                "Failed to submit batch; skipping %d sequences. "
                "Ensure the SignalP 6.0 API is reachable at %s.",
                len(batch),
                SIGNALP_URL,
            )
            continue

        job_id = submission.get("job_id") or submission.get("id")
        if job_id:
            results = poll_results(job_id)
        else:
            # Some API versions return results directly (synchronous mode).
            results = submission

        if results is None:
            logger.warning(
                "No results for batch (job %s); skipping %d sequences.",
                job_id,
                len(batch),
            )
            continue

        hits = _extract_signal_hits(results)
        surface_ids.update(hits)

    # ------------------------------------------------------------------
    # 3. Build filtered set and persist
    # ------------------------------------------------------------------
    filtered: list[tuple[str, str, str]] = []
    for gene_id, gene_name, sequence in proteins:
        if gene_id in surface_ids:
            filtered.append((gene_id, gene_name, sequence))

            candidate = Candidate(
                gene_id=gene_id,
                gene_name=gene_name,
                sequence=sequence,
                has_signal_peptide=True,
                filters_passed=["surface_filter"],
            )
            try:
                upsert_candidate(candidate)
            except Exception:
                logger.warning(
                    "Could not upsert candidate %s to Supabase; "
                    "continuing with local output.",
                    gene_id,
                )

    # ------------------------------------------------------------------
    # 4. Write filtered FASTA
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as fh:
        for gene_id, gene_name, sequence in filtered:
            fh.write(f">{gene_id} {gene_name}\n{sequence}\n")

    passed = len(filtered)
    logger.info(
        "%d of %d proteins passed surface filter.", passed, total
    )

    return OUTPUT_FILE


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Filter surface-exposed proteins using SignalP 6.0.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if the output file already exists.",
    )
    args = parser.parse_args()

    result_path = filter_surface_proteins(force=args.force)
    logger.info("Surface filter complete. Output: %s", result_path)
