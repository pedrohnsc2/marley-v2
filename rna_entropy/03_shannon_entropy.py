"""Calculate Shannon entropy profiles for *L. infantum* and human transcripts.

Shannon entropy quantifies the uncertainty (information content) of
nucleotide composition at each position in a sequence.  For individual
sequences (as opposed to multiple-sequence alignments), we compute
entropy over a sliding window to capture local compositional variation.

Sequences with low entropy are highly biased in nucleotide composition
(e.g., poly-A tracts), while high-entropy regions approach maximum
randomness.  Comparing entropy between parasite and host can reveal
regions where the parasite genome deviates from host expectations --
relevant for mRNA-based vaccine insert design.

Usage:
    python -m rna_entropy.03_shannon_entropy
    python -m rna_entropy.03_shannon_entropy --force
    python -m rna_entropy.03_shannon_entropy --dry-run
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path
from typing import Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LINF_FASTA: Final[str] = "data/raw/transcriptome_linf.fasta"
HUMAN_FASTA: Final[str] = "data/raw/transcriptome_human.fasta"
OUTPUT_DIR: Final[str] = "results/rna"
ENTROPY_CSV: Final[str] = "results/rna/shannon_entropy_profile.csv"

# Sliding window size for per-sequence entropy calculation.
# 50 nucleotides provides a reasonable balance between resolution
# and statistical stability for base frequency estimation.
WINDOW_SIZE: Final[int] = 50

# The four canonical nucleotides used for entropy calculation.
# Any non-standard characters (N, ambiguity codes) are excluded
# from frequency counts within each window.
NUCLEOTIDES: Final[str] = "ACGT"

# Maximum possible Shannon entropy for 4 equally-probable symbols:
#   H_max = log2(4) = 2.0 bits
# This serves as the theoretical upper bound for validation.
MAX_ENTROPY: Final[float] = 2.0

logger = get_logger("shannon_entropy")

# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def _parse_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """Parse a FASTA file into a list of (header, sequence) tuples.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        List of ``(header, sequence)`` tuples where *header* is the full
        description line (without ``>``) and *sequence* is the concatenated
        residue string in uppercase.

    Raises:
        FileNotFoundError: If *fasta_path* does not exist.
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    sequences: list[tuple[str, str]] = []
    current_header = ""
    current_seq: list[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq).upper()))
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    # Don't forget the last sequence.
    if current_header:
        sequences.append((current_header, "".join(current_seq).upper()))

    return sequences


# ---------------------------------------------------------------------------
# Shannon entropy computation
# ---------------------------------------------------------------------------


def _window_entropy(window: str) -> float:
    """Calculate Shannon entropy for a single nucleotide window.

    Shannon entropy measures the average information content per symbol:

        H(X) = -sum( p(x) * log2(p(x)) )  for x in {A, C, G, T}

    where p(x) is the relative frequency of nucleotide x in the window.

    Properties:
        - H = 0.0 bits when the window is perfectly conserved (only one
          nucleotide type present, e.g., "AAAAAAA").
        - H = 2.0 bits when all four nucleotides are equally frequent
          (maximum uncertainty / randomness).
        - The convention 0 * log2(0) = 0 is applied when a nucleotide
          has zero frequency, because lim(p->0) p*log2(p) = 0.

    Args:
        window: A string of nucleotide characters to analyze.  Only
            canonical bases (A, C, G, T) are counted; ambiguity codes
            and other characters are ignored.

    Returns:
        Shannon entropy in bits, in the range [0.0, 2.0].
    """
    # Count only canonical nucleotides within this window.
    counts = Counter(nt for nt in window if nt in NUCLEOTIDES)
    total = sum(counts.values())

    if total == 0:
        return 0.0

    entropy = 0.0
    for nt in NUCLEOTIDES:
        count = counts.get(nt, 0)
        if count == 0:
            # Convention: 0 * log2(0) = 0
            # This follows from the limit: lim(p->0+) p * log2(p) = 0
            # and ensures that absent nucleotides contribute nothing
            # to the entropy sum.
            continue

        # p(x) = frequency of nucleotide x in the window.
        p = count / total

        # Shannon entropy contribution: -p * log2(p)
        # Since 0 < p <= 1, log2(p) <= 0, so -p*log2(p) >= 0.
        entropy -= p * math.log2(p)

    return entropy


def _sequence_entropy(sequence: str) -> float:
    """Calculate mean Shannon entropy across sliding windows of a sequence.

    Slides a window of ``WINDOW_SIZE`` nucleotides across the sequence
    with a step of 1 (fully overlapping windows), computes Shannon entropy
    for each window, and returns the arithmetic mean.

    For sequences shorter than ``WINDOW_SIZE``, the entire sequence is
    treated as a single window.

    Args:
        sequence: Nucleotide sequence in uppercase (DNA alphabet).

    Returns:
        Mean Shannon entropy in bits, in the range [0.0, 2.0].
        Returns 0.0 for empty sequences.
    """
    # Normalize: convert RNA U -> DNA T for consistent processing.
    seq = sequence.replace("U", "T")

    if len(seq) == 0:
        return 0.0

    # If the sequence is shorter than the window, use the whole thing.
    if len(seq) <= WINDOW_SIZE:
        return _window_entropy(seq)

    # Slide a window across the sequence and collect entropy values.
    window_entropies: list[float] = []
    num_windows = len(seq) - WINDOW_SIZE + 1

    for start in range(num_windows):
        window = seq[start : start + WINDOW_SIZE]
        h = _window_entropy(window)
        window_entropies.append(h)

    if not window_entropies:
        return 0.0

    # Arithmetic mean of all window entropies.
    return sum(window_entropies) / len(window_entropies)


def _gc_content(sequence: str) -> float:
    """Calculate GC content of a nucleotide sequence.

    Args:
        sequence: Nucleotide sequence in uppercase.

    Returns:
        GC content as a fraction in [0.0, 1.0].
        Returns 0.0 for empty sequences.
    """
    seq = sequence.replace("U", "T")
    if len(seq) == 0:
        return 0.0

    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def _extract_gene_id(header: str) -> str:
    """Extract the gene/accession ID from a FASTA header.

    Handles common FASTA header formats:
        - UniProt: ``sp|P12345|GENE_ORGANISM description``
        - Plain:   ``ACCESSION description``

    Args:
        header: The FASTA header line (without ``>``).

    Returns:
        The extracted gene ID or accession.
    """
    # UniProt format: sp|ACCESSION|NAME or tr|ACCESSION|NAME
    if "|" in header:
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1]

    # Plain format: first whitespace-delimited token.
    return header.split()[0] if header.strip() else "unknown"


def _extract_gene_name(header: str) -> str:
    """Extract a short gene name from a FASTA header.

    Attempts to find a meaningful name from UniProt-style headers.
    Falls back to the gene ID if no name field is found.

    Args:
        header: The FASTA header line (without ``>``).

    Returns:
        A short gene name string.
    """
    # UniProt format: sp|ACCESSION|NAME_ORGANISM
    if "|" in header:
        parts = header.split("|")
        if len(parts) >= 3:
            name_part = parts[2].split()[0]
            # Remove the _ORGANISM suffix (e.g., "TryR_LEIIN" -> "TryR").
            if "_" in name_part:
                return name_part.split("_")[0]
            return name_part

    return _extract_gene_id(header)


# ---------------------------------------------------------------------------
# Batch entropy calculation
# ---------------------------------------------------------------------------


def _compute_entropy_batch(
    sequences: list[tuple[str, str]],
) -> list[dict[str, str | float]]:
    """Calculate Shannon entropy and GC content for a batch of sequences.

    Args:
        sequences: List of ``(header, sequence)`` tuples.

    Returns:
        List of dicts with keys: gene_id, gene_name, length, gc_content,
        shannon_entropy.
    """
    results: list[dict[str, str | float]] = []

    for header, seq in sequences:
        gene_id = _extract_gene_id(header)
        gene_name = _extract_gene_name(header)
        length = len(seq)
        gc = _gc_content(seq)
        entropy = _sequence_entropy(seq)

        results.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "length": length,
            "gc_content": round(gc, 4),
            "shannon_entropy": round(entropy, 4),
        })

    return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def calculate_shannon_entropy(force: bool = False, dry_run: bool = False) -> str:
    """Calculate Shannon entropy profiles for L. infantum and human transcripts.

    For each gene/transcript, computes the mean Shannon entropy using a
    sliding window approach.  Also calculates GC content and the entropy
    difference (delta) between human and parasite sequences.

    The entropy delta highlights genes where the parasite's nucleotide
    composition deviates most from human patterns -- potentially useful
    for mRNA vaccine insert design where host immune recognition matters.

    Args:
        force: Re-analyze even if the output file already exists.
        dry_run: Process only the first 10 sequences per organism and
            log results without writing the full output.

    Returns:
        Path to the Shannon entropy profile CSV.

    Raises:
        FileNotFoundError: If input FASTA files are missing.
    """
    csv_path = Path(ENTROPY_CSV)

    # --- Check for existing output --------------------------------------------
    if csv_path.exists() and not force:
        logger.info(
            "Entropy profile already exists at %s. Use --force to re-analyze.",
            csv_path,
        )
        return str(csv_path)

    # --- Parse FASTA files ----------------------------------------------------
    logger.info("Reading L. infantum transcriptome: %s", LINF_FASTA)
    linf_seqs = _parse_fasta(LINF_FASTA)
    logger.info("Read %d L. infantum sequences.", len(linf_seqs))

    logger.info("Reading human transcriptome: %s", HUMAN_FASTA)
    human_seqs = _parse_fasta(HUMAN_FASTA)
    logger.info("Read %d human sequences.", len(human_seqs))

    # --- Dry run: limit to first 10 sequences ---------------------------------
    if dry_run:
        linf_seqs = linf_seqs[:10]
        human_seqs = human_seqs[:10]
        logger.info("[DRY RUN] Limited to first 10 sequences per organism.")

    # --- Calculate L. infantum entropy ----------------------------------------
    logger.info("Calculating Shannon entropy for L. infantum sequences...")
    linf_results = _compute_entropy_batch(linf_seqs)

    # --- Calculate human entropy (aggregate mean for comparison) ---------------
    logger.info("Calculating Shannon entropy for human sequences...")
    human_results = _compute_entropy_batch(human_seqs)

    # Compute the mean human entropy as a reference baseline.
    # Each L. infantum gene is compared against this aggregate human mean
    # because we don't have 1:1 ortholog mapping at this stage.
    human_entropies = [r["shannon_entropy"] for r in human_results]
    mean_human_entropy = (
        sum(human_entropies) / len(human_entropies) if human_entropies else 0.0
    )
    logger.info("Mean human Shannon entropy: %.4f bits", mean_human_entropy)

    # --- Compute entropy delta for each L. infantum gene ----------------------
    # entropy_delta = human_entropy - shannon_entropy
    # Positive delta: parasite gene has lower entropy than human average
    #   (more compositionally biased -- potentially interesting for vaccines).
    # Negative delta: parasite gene has higher entropy than human average.
    for record in linf_results:
        record["human_entropy"] = round(mean_human_entropy, 4)
        record["entropy_delta"] = round(
            mean_human_entropy - record["shannon_entropy"], 4
        )

    # --- Summary statistics ---------------------------------------------------
    linf_entropies = [r["shannon_entropy"] for r in linf_results]
    mean_linf_entropy = (
        sum(linf_entropies) / len(linf_entropies) if linf_entropies else 0.0
    )
    deltas = [r["entropy_delta"] for r in linf_results]
    mean_delta = sum(deltas) / len(deltas) if deltas else 0.0

    logger.info(
        "Calculated entropy for %d genes. Mean entropy: %.2f bits. Mean delta: %.2f bits",
        len(linf_results),
        mean_linf_entropy,
        mean_delta,
    )

    # --- Dry run: log results but skip file write -----------------------------
    if dry_run:
        logger.info("[DRY RUN] Top 5 genes by entropy delta (most biased):")
        sorted_by_delta = sorted(linf_results, key=lambda r: r["entropy_delta"], reverse=True)
        for rec in sorted_by_delta[:5]:
            logger.info(
                "[DRY RUN]   %s (%s): entropy=%.4f, delta=%.4f",
                rec["gene_id"],
                rec["gene_name"],
                rec["shannon_entropy"],
                rec["entropy_delta"],
            )
        return str(csv_path)

    # --- Write CSV output -----------------------------------------------------
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "gene_id", "gene_name", "length", "gc_content",
            "shannon_entropy", "human_entropy", "entropy_delta",
        ])
        for rec in linf_results:
            writer.writerow([
                rec["gene_id"],
                rec["gene_name"],
                rec["length"],
                rec["gc_content"],
                rec["shannon_entropy"],
                rec["human_entropy"],
                rec["entropy_delta"],
            ])

    logger.info("Wrote entropy profiles to %s", csv_path)

    return str(csv_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculate Shannon entropy profiles for L. infantum transcripts.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-analyze even if entropy profile already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Process only first 10 sequences and log results.",
    )
    args = parser.parse_args()

    result = calculate_shannon_entropy(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
