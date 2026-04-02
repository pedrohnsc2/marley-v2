"""Assess sequence conservation of candidate antigens across Leishmania species via BLAST alignment.

Runs remote BLAST (blastp) against the NCBI non-redundant database, restricted
to the Leishmania taxon (taxid 5658).  For each surface protein candidate the
script computes an average percent-identity score across comparison strains and
retains only those above a configurable threshold.

Usage:
    python -m pipeline.03_conservation          # normal run
    python -m pipeline.03_conservation --force  # re-run even if output exists
"""

from __future__ import annotations

import csv
import io
import os
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from tqdm import tqdm

from core.db import update_candidate
from core.logger import get_logger
from core.models import STATUS_APPROVED, STATUS_REJECTED

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INPUT_FILE: str = "data/raw/surface_proteins.fasta"
OUTPUT_FILE: str = "results/conserved_candidates.csv"
CONSERVATION_THRESHOLD: float = 0.80
REFERENCE_STRAIN: str = "Leishmania infantum JPCM5"
COMPARISON_STRAINS: list[str] = [
    "Leishmania infantum LEM3323",
    "Leishmania chagasi",
]
BLAST_DB: str = "nr"  # NCBI non-redundant database
BLAST_PROGRAM: str = "blastp"
BLAST_EVALUE: float = 1e-10
BLAST_MAX_HITS: int = 50

# Leishmania taxon id used to restrict BLAST searches
_LEISHMANIA_TAXID: int = 5658

# Retry / timeout settings for remote BLAST
_BLAST_RETRIES: int = 3
_BLAST_RETRY_DELAY_SECONDS: int = 30

logger = get_logger("conservation")

# ---------------------------------------------------------------------------
# BLAST helpers
# ---------------------------------------------------------------------------


def run_blast(sequence: str, gene_id: str) -> list[dict[str, Any]]:
    """Run a remote BLAST search and return parsed hit information.

    Uses Biopython's ``NCBIWWW.qblast`` to query the NCBI non-redundant
    database (``nr``) with blastp, restricted to the Leishmania taxon
    (taxid 5658).

    Args:
        sequence: Amino-acid sequence string to search.
        gene_id: Identifier for the query, used only in log messages.

    Returns:
        A list of dictionaries, each containing:
            - ``subject_id``       (str):   Accession of the matched subject.
            - ``identity_percent`` (float): Percent identity of the best HSP.
            - ``evalue``           (float): E-value of the best HSP.
            - ``organism``         (str):   Organism name parsed from the hit
                                            definition line.

    Raises:
        RuntimeError: If all retry attempts are exhausted.
    """
    hits: list[dict[str, Any]] = []
    entrez_query = f"txid{_LEISHMANIA_TAXID}[ORGN]"

    last_error: Exception | None = None

    for attempt in range(1, _BLAST_RETRIES + 1):
        try:
            logger.info(
                "BLAST search for %s (attempt %d/%d) ...",
                gene_id,
                attempt,
                _BLAST_RETRIES,
            )
            result_handle = NCBIWWW.qblast(
                program=BLAST_PROGRAM,
                database=BLAST_DB,
                sequence=sequence,
                entrez_query=entrez_query,
                expect=BLAST_EVALUE,
                hitlist_size=BLAST_MAX_HITS,
            )
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

            for alignment in blast_record.alignments:
                if not alignment.hsps:
                    continue

                best_hsp = alignment.hsps[0]
                identity_percent = (
                    best_hsp.identities / best_hsp.align_length * 100.0
                    if best_hsp.align_length > 0
                    else 0.0
                )

                # Attempt to extract organism from the hit definition.
                # NCBI definitions typically contain the organism in brackets:
                #   "hypothetical protein [Leishmania infantum JPCM5]"
                organism = _parse_organism(alignment.hit_def)

                hits.append(
                    {
                        "subject_id": alignment.accession,
                        "identity_percent": round(identity_percent, 2),
                        "evalue": best_hsp.expect,
                        "organism": organism,
                    }
                )

            logger.info(
                "BLAST for %s returned %d hit(s).", gene_id, len(hits)
            )
            return hits

        except Exception as exc:
            last_error = exc
            logger.warning(
                "BLAST attempt %d/%d for %s failed: %s",
                attempt,
                _BLAST_RETRIES,
                gene_id,
                exc,
            )
            if attempt < _BLAST_RETRIES:
                logger.info(
                    "Retrying in %d seconds ...", _BLAST_RETRY_DELAY_SECONDS
                )
                time.sleep(_BLAST_RETRY_DELAY_SECONDS)

    msg = (
        f"All {_BLAST_RETRIES} BLAST attempts failed for {gene_id}. "
        f"Last error: {last_error}"
    )
    logger.error(msg)
    raise RuntimeError(msg)


def _parse_organism(hit_def: str) -> str:
    """Extract the organism name from a BLAST hit definition line.

    NCBI hit definitions typically end with the organism in square brackets,
    e.g. ``"hypothetical protein [Leishmania infantum JPCM5]"``.

    Args:
        hit_def: The raw ``hit_def`` string from a BLAST alignment.

    Returns:
        The organism name, or ``"unknown"`` if parsing fails.
    """
    try:
        start = hit_def.rindex("[")
        end = hit_def.rindex("]")
        if start < end:
            return hit_def[start + 1 : end].strip()
    except ValueError:
        pass
    return "unknown"


# ---------------------------------------------------------------------------
# Conservation scoring
# ---------------------------------------------------------------------------


def calculate_conservation(blast_hits: list[dict[str, Any]]) -> float:
    """Compute an average conservation score from BLAST hits.

    Filters *blast_hits* for entries whose ``organism`` field matches any of
    the ``COMPARISON_STRAINS``.  For each matching strain the highest
    percent-identity hit is selected, and the final score is the average of
    those best-per-strain values, normalised to [0.0, 1.0].

    Args:
        blast_hits: Output of :func:`run_blast`.

    Returns:
        A conservation score in the range [0.0, 1.0], or ``0.0`` when no
        matching strains are found.
    """
    if not blast_hits:
        return 0.0

    # Collect best identity per comparison strain.
    best_per_strain: dict[str, float] = {}

    for hit in blast_hits:
        organism = hit.get("organism", "")
        for strain in COMPARISON_STRAINS:
            if strain.lower() in organism.lower():
                identity = hit.get("identity_percent", 0.0)
                if identity > best_per_strain.get(strain, 0.0):
                    best_per_strain[strain] = identity

    if not best_per_strain:
        return 0.0

    # Normalise from 0-100 percent to 0.0-1.0 score.
    average_identity = sum(best_per_strain.values()) / len(best_per_strain)
    score = round(average_identity / 100.0, 4)

    return min(score, 1.0)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------


def analyze_conservation(force: bool = False) -> str:
    """Run the conservation analysis pipeline step.

    Reads surface-protein candidates from *INPUT_FILE*, performs a BLAST
    search for each, computes conservation scores, and writes retained
    candidates to *OUTPUT_FILE*.

    Args:
        force: When ``True``, re-run even if *OUTPUT_FILE* already exists.

    Returns:
        The path to the output CSV file.

    Raises:
        FileNotFoundError: If *INPUT_FILE* does not exist.
    """
    # Early exit if output already present and force is not set.
    if not force and Path(OUTPUT_FILE).exists():
        logger.info(
            "Output file %s already exists. Use force=True to re-run.",
            OUTPUT_FILE,
        )
        return OUTPUT_FILE

    # Validate input.
    if not Path(INPUT_FILE).exists():
        msg = f"Input file not found: {INPUT_FILE}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    # Ensure output directory exists.
    Path(OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)

    # Parse input FASTA.
    records = list(SeqIO.parse(INPUT_FILE, "fasta"))
    total = len(records)
    logger.info("Loaded %d surface protein(s) from %s.", total, INPUT_FILE)

    if total == 0:
        logger.warning("No sequences found in %s. Nothing to do.", INPUT_FILE)
        # Write an empty CSV with headers only.
        _write_empty_csv()
        return OUTPUT_FILE

    results: list[dict[str, Any]] = []
    passed = 0

    for record in tqdm(records, desc="Conservation analysis", unit="seq"):
        gene_id = record.id
        gene_name = record.description.split(None, 1)[-1] if " " in record.description else gene_id
        sequence = str(record.seq)

        logger.info("Analysing %s (%s) ...", gene_id, gene_name)

        try:
            blast_hits = run_blast(sequence, gene_id)
        except RuntimeError:
            logger.error(
                "Skipping %s due to BLAST failure.", gene_id
            )
            score = 0.0
            blast_hits = []
        else:
            score = calculate_conservation(blast_hits)

        status = (
            STATUS_APPROVED
            if score >= CONSERVATION_THRESHOLD
            else STATUS_REJECTED
        )

        if score >= CONSERVATION_THRESHOLD:
            passed += 1

        logger.info(
            "%s conservation=%.4f status=%s", gene_id, score, status
        )

        # Persist score to Supabase.
        try:
            update_candidate(
                gene_id,
                {"conservation_score": score, "status": status},
            )
        except Exception:
            logger.warning(
                "Could not update Supabase for %s. Continuing.", gene_id
            )

        results.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "conservation_score": round(score, 4),
                "status": status,
            }
        )

    # Build DataFrame and write CSV.
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_FILE, index=False, quoting=csv.QUOTE_NONNUMERIC)

    logger.info(
        "%d of %d candidates passed conservation filter (threshold: %.2f)",
        passed,
        total,
        CONSERVATION_THRESHOLD,
    )
    logger.info("Results written to %s.", OUTPUT_FILE)

    return OUTPUT_FILE


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _write_empty_csv() -> None:
    """Write a CSV with headers only to *OUTPUT_FILE*."""
    df = pd.DataFrame(
        columns=["gene_id", "gene_name", "conservation_score", "status"]
    )
    df.to_csv(OUTPUT_FILE, index=False, quoting=csv.QUOTE_NONNUMERIC)
    logger.info("Wrote empty results to %s.", OUTPUT_FILE)


# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    force_flag = "--force" in sys.argv
    try:
        output_path = analyze_conservation(force=force_flag)
        logger.info("Conservation analysis complete. Output: %s", output_path)
    except Exception:
        logger.exception("Conservation analysis failed.")
        sys.exit(1)
