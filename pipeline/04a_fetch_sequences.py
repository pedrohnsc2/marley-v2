"""
Fetch protein sequences for validated antigens from UniProt.

This pipeline step sits between immunogenicity scoring (step 04) and
downstream epitope prediction.  It resolves placeholder sequences for
experimentally validated antigens by querying the UniProt REST API,
caches results to a local FASTA file, and updates the Supabase
``candidates`` table with real amino-acid sequences.
"""

from __future__ import annotations

from pathlib import Path

from core.db import update_candidate
from core.logger import get_logger
from core.uniprot import fetch_all_validated_sequences

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTPUT_FASTA: str = "data/validated_sequences.fasta"

logger = get_logger("fetch_sequences")

# ---------------------------------------------------------------------------
# Main pipeline step
# ---------------------------------------------------------------------------


def fetch_validated_sequences(force: bool = False) -> str:
    """Fetch protein sequences for validated antigens from UniProt.

    Calls :func:`core.uniprot.fetch_all_validated_sequences` to resolve
    amino-acid sequences, then:

    1. Updates Supabase ``candidates`` rows with real sequences.
    2. Writes a cached multi-FASTA file to ``data/validated_sequences.fasta``.

    Args:
        force: If ``True``, re-fetch even when the output FASTA already
               exists.

    Returns:
        Path to the output FASTA file.
    """
    output_path = Path(OUTPUT_FASTA)

    if output_path.exists() and not force:
        logger.info(
            "Output file %s already exists. Use force=True to re-fetch.",
            OUTPUT_FASTA,
        )
        return OUTPUT_FASTA

    # Ensure output directory exists.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Fetch sequences from UniProt
    # ------------------------------------------------------------------
    logger.info("Fetching validated antigen sequences from UniProt ...")
    sequences = fetch_all_validated_sequences()

    if not sequences:
        logger.warning("No sequences retrieved. Nothing to write.")
        return OUTPUT_FASTA

    # ------------------------------------------------------------------
    # 2. Update Supabase candidates table
    # ------------------------------------------------------------------
    updated = 0
    for gene_id, sequence in sequences.items():
        try:
            update_candidate(gene_id, {"sequence": sequence})
            updated += 1
            logger.info("Updated sequence for %s in Supabase.", gene_id)
        except Exception:
            logger.warning(
                "Could not update Supabase for %s. Continuing.", gene_id,
            )

    logger.info("Updated %d/%d candidate(s) in Supabase.", updated, len(sequences))

    # ------------------------------------------------------------------
    # 3. Write cached FASTA
    # ------------------------------------------------------------------
    _write_fasta(sequences, output_path)

    return OUTPUT_FASTA


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Write gene_id -> sequence pairs to a multi-FASTA file.

    Args:
        sequences: Mapping of gene identifiers to amino-acid sequences.
        path: Destination file path.
    """
    with open(path, "w") as fh:
        for gene_id, sequence in sequences.items():
            fh.write(f">{gene_id}\n")
            # Wrap at 80 characters per FASTA convention.
            for i in range(0, len(sequence), 80):
                fh.write(sequence[i : i + 80] + "\n")

    logger.info("Wrote %d sequence(s) to %s.", len(sequences), path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Step 4a: Fetch validated antigen sequences from UniProt.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-fetch even if the output FASTA already exists.",
    )
    args = parser.parse_args()

    result_path = fetch_validated_sequences(force=args.force)
    print(f"Output: {result_path}")
