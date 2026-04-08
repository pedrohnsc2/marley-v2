"""Analyze codon usage bias between *L. infantum* and human transcriptomes.

Calculates Relative Synonymous Codon Usage (RSCU) for both organisms,
computes the euclidean distance between their RSCU vectors, and derives
a normalized codon bias score.  This score quantifies how different the
parasite's codon preferences are from the human host -- a high score
indicates strong codon bias that could be exploited for vaccine design.

Usage:
    python -m rna_entropy.02_codon_usage
    python -m rna_entropy.02_codon_usage --force
    python -m rna_entropy.02_codon_usage --dry-run
"""

from __future__ import annotations

import csv
import json
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
CODON_CSV: Final[str] = "results/rna/codon_usage_comparison.csv"
STATS_JSON: Final[str] = "results/rna/codon_usage_stats.json"

logger = get_logger("codon_usage")

# ---------------------------------------------------------------------------
# Genetic code (standard)
# ---------------------------------------------------------------------------
# Maps each of the 64 codons to its corresponding amino acid.
# Stop codons are represented as "*".

GENETIC_CODE: Final[dict[str, str]] = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Stop codons are excluded from RSCU analysis.
STOP_CODONS: Final[set[str]] = {"TAA", "TAG", "TGA"}

# Build reverse map: amino acid -> list of synonymous codons.
AMINO_ACID_CODONS: Final[dict[str, list[str]]] = {}
for _codon, _aa in GENETIC_CODE.items():
    if _aa != "*":
        AMINO_ACID_CODONS.setdefault(_aa, []).append(_codon)


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
# Coding region extraction
# ---------------------------------------------------------------------------


def _extract_coding_region(sequence: str) -> str | None:
    """Extract the coding region from ATG to the first in-frame stop codon.

    Scans the sequence for the first ATG start codon, then reads triplets
    until a stop codon (TAA, TAG, TGA) is encountered or the sequence ends.

    Args:
        sequence: Nucleotide sequence in uppercase (DNA alphabet: A, T, C, G).

    Returns:
        The coding region as a string of complete codons (including ATG but
        excluding the stop codon), or ``None`` if no ATG is found.
    """
    # Convert any RNA (U) to DNA (T) for uniform processing.
    sequence = sequence.replace("U", "T")

    # Find the first ATG.
    start_idx = sequence.find("ATG")
    if start_idx == -1:
        return None

    coding = []
    pos = start_idx
    while pos + 3 <= len(sequence):
        codon = sequence[pos : pos + 3]
        if codon in STOP_CODONS:
            break
        coding.append(codon)
        pos += 3

    # Require at least one codon beyond the start.
    if len(coding) < 2:
        return None

    return "".join(coding)


# ---------------------------------------------------------------------------
# Codon counting and RSCU
# ---------------------------------------------------------------------------


def _count_codons(sequences: list[tuple[str, str]]) -> Counter[str]:
    """Count all 64 codons across coding regions of the given sequences.

    Args:
        sequences: List of ``(header, sequence)`` tuples.

    Returns:
        A ``Counter`` mapping each codon to its observed count.
        Codons not found in any coding region have a count of 0.
    """
    codon_counts: Counter[str] = Counter()

    # Initialize all codons to ensure every codon is present.
    for codon in GENETIC_CODE:
        codon_counts[codon] = 0

    extracted = 0
    for _header, seq in sequences:
        cds = _extract_coding_region(seq)
        if cds is None:
            continue
        extracted += 1
        for i in range(0, len(cds), 3):
            codon = cds[i : i + 3]
            if len(codon) == 3 and codon in GENETIC_CODE:
                codon_counts[codon] += 1

    logger.info("Extracted coding regions from %d / %d sequences.", extracted, len(sequences))
    return codon_counts


def _calculate_rscu(codon_counts: Counter[str]) -> dict[str, float]:
    """Calculate Relative Synonymous Codon Usage (RSCU) for each codon.

    RSCU normalizes codon frequency by the number of synonymous codons
    encoding the same amino acid.  Under equal usage, RSCU = 1.0 for
    every codon.  Values > 1.0 indicate preference; values < 1.0
    indicate avoidance.

    Formula:
        RSCU(ij) = X(ij) / (1/n_i * sum(X(ik) for all synonyms k))

    where:
        - X(ij) is the observed count of codon j for amino acid i
        - n_i is the number of synonymous codons for amino acid i

    Args:
        codon_counts: Counter mapping codons to observed counts.

    Returns:
        Dict mapping each non-stop codon to its RSCU value.
    """
    rscu: dict[str, float] = {}

    for amino_acid, synonymous_codons in AMINO_ACID_CODONS.items():
        n_synonyms = len(synonymous_codons)

        # Total count of all synonymous codons for this amino acid.
        total_count = sum(codon_counts[c] for c in synonymous_codons)

        if total_count == 0:
            # No observations: assign RSCU = 1.0 (neutral) by convention.
            for codon in synonymous_codons:
                rscu[codon] = 1.0
            continue

        # Expected count per codon under equal usage.
        expected = total_count / n_synonyms

        for codon in synonymous_codons:
            rscu[codon] = codon_counts[codon] / expected

    return rscu


def _calculate_gc_content(sequences: list[tuple[str, str]]) -> float:
    """Calculate overall GC content across all sequences.

    Args:
        sequences: List of ``(header, sequence)`` tuples.

    Returns:
        GC content as a fraction in [0.0, 1.0].
    """
    total_gc = 0
    total_len = 0

    for _header, seq in sequences:
        seq_upper = seq.upper().replace("U", "T")
        total_gc += seq_upper.count("G") + seq_upper.count("C")
        total_len += len(seq_upper)

    if total_len == 0:
        return 0.0

    return total_gc / total_len


def _euclidean_distance(vec_a: dict[str, float], vec_b: dict[str, float]) -> float:
    """Calculate euclidean distance between two RSCU vectors.

    Both vectors must have the same set of keys (codons).

    Args:
        vec_a: First RSCU vector.
        vec_b: Second RSCU vector.

    Returns:
        Euclidean distance as a non-negative float.
    """
    common_keys = set(vec_a.keys()) & set(vec_b.keys())
    sum_sq = sum((vec_a[k] - vec_b[k]) ** 2 for k in common_keys)
    return math.sqrt(sum_sq)


def _normalize_distance(distance: float) -> float:
    """Normalize euclidean RSCU distance to a [0.0, 1.0] codon bias score.

    The theoretical maximum euclidean distance between two RSCU vectors
    depends on the number of synonymous codon families.  For the standard
    genetic code with 59 sense codons across 18 amino acid families (plus
    Met and Trp which have single codons), empirical maximum distance is
    approximately 12.0.  We use a sigmoid-like normalization:

        score = 1 - exp(-distance / scale)

    where scale = 4.0 produces a useful spread across biologically
    realistic distances (typically 1.0 -- 8.0).

    Args:
        distance: Raw euclidean distance between RSCU vectors.

    Returns:
        Normalized score in [0.0, 1.0].
    """
    scale: Final[float] = 4.0
    return 1.0 - math.exp(-distance / scale)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyze_codon_usage(force: bool = False, dry_run: bool = False) -> str:
    """Compare codon usage between L. infantum and human transcriptomes.

    Reads FASTA files produced by the fetch step, computes RSCU for both
    organisms, and quantifies the difference as a codon bias score.

    Args:
        force: Re-analyze even if output files already exist.
        dry_run: Log analysis plan using hardcoded reference values.

    Returns:
        Path to the codon usage comparison CSV.

    Raises:
        FileNotFoundError: If input FASTA files are missing.
        RuntimeError: If no coding regions are found in either organism.
    """
    csv_path = Path(CODON_CSV)
    stats_path = Path(STATS_JSON)

    # --- Check for existing output --------------------------------------------
    if csv_path.exists() and stats_path.exists() and not force:
        logger.info(
            "Codon usage files already exist. Use --force to re-analyze."
        )
        return str(csv_path)

    # --- Dry run mode ---------------------------------------------------------
    if dry_run:
        # Hardcoded reference values from literature.
        logger.info("[DRY RUN] Would analyze codon usage from:")
        logger.info("[DRY RUN]   L. infantum: %s", LINF_FASTA)
        logger.info("[DRY RUN]   Human:       %s", HUMAN_FASTA)
        logger.info("[DRY RUN] Known reference values:")
        logger.info("[DRY RUN]   L. infantum GC content: ~61%%")
        logger.info("[DRY RUN]   Human GC content:       ~50%%")
        logger.info("[DRY RUN]   Expected codon bias score: 0.60 -- 0.80")
        return str(csv_path)

    # --- Parse FASTA files ----------------------------------------------------
    logger.info("Reading L. infantum transcriptome: %s", LINF_FASTA)
    linf_seqs = _parse_fasta(LINF_FASTA)
    logger.info("Read %d L. infantum sequences.", len(linf_seqs))

    logger.info("Reading human transcriptome: %s", HUMAN_FASTA)
    human_seqs = _parse_fasta(HUMAN_FASTA)
    logger.info("Read %d human sequences.", len(human_seqs))

    # --- Count codons and compute RSCU ----------------------------------------
    logger.info("Counting codons in L. infantum coding regions...")
    linf_counts = _count_codons(linf_seqs)

    logger.info("Counting codons in human coding regions...")
    human_counts = _count_codons(human_seqs)

    linf_rscu = _calculate_rscu(linf_counts)
    human_rscu = _calculate_rscu(human_counts)

    # --- GC content -----------------------------------------------------------
    linf_gc = _calculate_gc_content(linf_seqs)
    human_gc = _calculate_gc_content(human_seqs)
    logger.info("GC content: L. infantum=%.1f%%, Human=%.1f%%", linf_gc * 100, human_gc * 100)

    # --- Euclidean distance and bias score ------------------------------------
    distance = _euclidean_distance(linf_rscu, human_rscu)
    bias_score = _normalize_distance(distance)
    logger.info(
        "Codon bias score: %.2f (euclidean distance: %.2f)",
        bias_score,
        distance,
    )

    # --- Write CSV output -----------------------------------------------------
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort codons alphabetically for reproducible output.
    sorted_codons = sorted(c for c in GENETIC_CODE if c not in STOP_CODONS)

    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "codon", "amino_acid",
            "linf_count", "linf_rscu",
            "human_count", "human_rscu",
        ])
        for codon in sorted_codons:
            amino_acid = GENETIC_CODE[codon]
            writer.writerow([
                codon,
                amino_acid,
                linf_counts[codon],
                f"{linf_rscu.get(codon, 0.0):.4f}",
                human_counts[codon],
                f"{human_rscu.get(codon, 0.0):.4f}",
            ])

    logger.info("Wrote codon comparison to %s", csv_path)

    # --- Write stats JSON -----------------------------------------------------
    stats = {
        "overall_distance": round(distance, 4),
        "codon_bias_score": round(bias_score, 4),
        "linf_gc_content": round(linf_gc, 4),
        "human_gc_content": round(human_gc, 4),
        "linf_sequence_count": len(linf_seqs),
        "human_sequence_count": len(human_seqs),
    }

    with open(stats_path, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2)

    logger.info("Wrote usage stats to %s", stats_path)

    return str(csv_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze codon usage bias: L. infantum vs. human.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-analyze even if output files already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log analysis plan without processing data.",
    )
    args = parser.parse_args()

    result = analyze_codon_usage(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
