"""Filter proteins for surface-exposed and secreted candidates using signal peptide predictions.

Uses the SignalP 6.0 model via the BioLib Python SDK (pybiolib) to classify
signal peptides, retaining only proteins with a classic Sec/SPI signal
peptide -- strong indicators of surface exposure or secretion, both desirable
traits for vaccine targets.
"""

from __future__ import annotations

import csv
import io
import tempfile
from pathlib import Path

import biolib
from Bio import SeqIO
from tqdm import tqdm

from core.db import upsert_candidate
from core.logger import get_logger
from core.models import Candidate

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INPUT_FILE: str = "data/raw/proteins.fasta"
OUTPUT_FILE: str = "data/raw/surface_proteins.fasta"
CHUNK_SIZE: int = 500
BIOLIB_MODEL: str = "DTU/SignalP-6"

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
# SignalP 6.0 via BioLib SDK
# ---------------------------------------------------------------------------


def _write_temp_fasta(sequences: list[tuple[str, str]]) -> str:
    """Write sequences to a temporary FASTA file for BioLib input.

    Args:
        sequences: List of ``(header, sequence)`` pairs.

    Returns:
        The absolute path to the temporary FASTA file.
    """
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, prefix="signalp_chunk_"
    )
    try:
        for header, seq in sequences:
            tmp.write(f">{header}\n{seq}\n")
    finally:
        tmp.close()

    return tmp.name


def run_signalp(fasta_path: str) -> dict:
    """Run SignalP 6.0 on a FASTA file using the BioLib Python SDK.

    The BioLib SDK executes the DTU/SignalP-6 model synchronously and returns
    a tab-delimited text file with signal peptide predictions.

    Args:
        fasta_path: Path to a FASTA file containing protein sequences.

    Returns:
        A dictionary mapping gene_id to its prediction metadata::

            {"gene_id": {"prediction": "SP", "probability": 0.9998}, ...}

    Raises:
        RuntimeError: If the BioLib job fails or returns no parseable output.
    """
    logger.info("Running SignalP 6.0 via BioLib on %s ...", fasta_path)

    signalp = biolib.load(BIOLIB_MODEL)
    result = signalp.cli(
        args=(
            f"--fastafile {fasta_path} "
            f"--output_dir output "
            f"--organism eukarya "
            f"--format txt "
            f"--mode fast"
        )
    )

    # BioLib returns LazyLoadedFile objects; look for the TSV results file.
    output_files = result.list_output_files()
    logger.debug(
        "BioLib output files: %s",
        [f.path for f in output_files],
    )

    # SignalP 6.0 via BioLib outputs tab-delimited .txt files
    # (e.g. prediction_results.txt).  Columns:
    #   ID  Prediction  OTHER  SP(Sec/SPI)  CS Position
    for f in output_files:
        if not f.path.endswith(".txt"):
            continue

        content = f.get_data().decode("utf-8")
        reader = csv.DictReader(
            io.StringIO(content), delimiter="\t"
        )

        predictions: dict[str, dict] = {}
        for row in reader:
            gene_id = row.get("ID", "").strip()
            if not gene_id:
                continue
            prediction = row.get("Prediction", "").strip()
            # The SP(Sec/SPI) column holds the probability score.
            probability_str = row.get("SP(Sec/SPI)", "0").strip()
            try:
                probability = float(probability_str)
            except ValueError:
                probability = 0.0

            predictions[gene_id] = {
                "prediction": prediction,
                "probability": probability,
            }

        if predictions:
            logger.info(
                "SignalP returned predictions for %d sequences.",
                len(predictions),
            )
            return predictions

    raise RuntimeError(
        f"SignalP via BioLib produced no parseable .txt output. "
        f"Output files found: {[f.path for f in output_files]}"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _extract_signal_hits(results: dict) -> set[str]:
    """Extract gene_ids classified as having a signal peptide (SP).

    Expects the flat dict returned by :func:`run_signalp`::

        {"gene_id": {"prediction": "SP", "probability": 0.9998}, ...}

    A protein is considered a hit when its ``prediction`` value equals
    ``"SP"`` -- the SignalP 6.0 label for a classic Sec/SPI signal peptide.

    Args:
        results: Dict returned by :func:`run_signalp`.

    Returns:
        A set of gene_id strings that have a signal peptide prediction.
    """
    hits: set[str] = set()

    for gene_id, info in results.items():
        if not isinstance(info, dict):
            continue
        prediction = info.get("prediction", "")
        if prediction == "SP":
            hits.add(gene_id)

    return hits


# ---------------------------------------------------------------------------
# Main filter function
# ---------------------------------------------------------------------------


def filter_surface_proteins(force: bool = False) -> str:
    """Run the surface-protein filter using SignalP 6.0 via BioLib.

    Reads protein sequences from :data:`INPUT_FILE`, processes them in chunks
    through SignalP via the BioLib SDK, retains those classified as *Sec/SPI*,
    persists each passing candidate to Supabase, and writes the filtered FASTA
    to :data:`OUTPUT_FILE`.

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
    # 2. Process in chunks via SignalP (BioLib SDK)
    # ------------------------------------------------------------------
    surface_ids: set[str] = set()
    chunks = [
        proteins[i : i + CHUNK_SIZE] for i in range(0, total, CHUNK_SIZE)
    ]

    for chunk in tqdm(chunks, desc="SignalP chunks", unit="chunk"):
        chunk_seqs: list[tuple[str, str]] = [
            (gene_id, seq) for gene_id, _, seq in chunk
        ]

        tmp_fasta: str | None = None
        try:
            tmp_fasta = _write_temp_fasta(chunk_seqs)
            results = run_signalp(tmp_fasta)
            hits = _extract_signal_hits(results)
            surface_ids.update(hits)
            logger.info(
                "Chunk: %d/%d sequences had signal peptides.",
                len(hits),
                len(chunk),
            )
        except Exception:
            logger.exception(
                "Failed to process chunk of %d sequences via BioLib; "
                "skipping this chunk.",
                len(chunk),
            )
        finally:
            # Clean up the temporary FASTA file.
            if tmp_fasta is not None:
                try:
                    Path(tmp_fasta).unlink(missing_ok=True)
                except OSError:
                    pass

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
