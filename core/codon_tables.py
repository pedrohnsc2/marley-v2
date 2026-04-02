"""
Codon usage table for *Canis lupus familiaris* (taxonomy ID 9615).

Frequencies are derived from the Kazusa Codon Usage Database and represent
the fraction of usage for each codon within its synonymous group.  Values
are rounded to three decimal places and sum to ~1.0 per amino acid.

Public API
----------
CODON_TABLE : dict[str, list[tuple[str, float]]]
    Mapping from single-letter amino acid (plus ``"*"`` for stop) to a list
    of ``(codon, frequency)`` tuples.

get_optimal_codon(amino_acid) -> str
    Return the highest-frequency codon for the given amino acid.

reverse_translate(protein_seq) -> str
    Convert a protein sequence to a DNA sequence using optimal codons.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Canis lupus familiaris codon usage  (Kazusa, taxonomy ID 9615)
# ---------------------------------------------------------------------------

CODON_TABLE: dict[str, list[tuple[str, float]]] = {
    # Phenylalanine
    "F": [
        ("TTT", 0.427),
        ("TTC", 0.573),
    ],
    # Leucine
    "L": [
        ("TTA", 0.072),
        ("TTG", 0.126),
        ("CTT", 0.128),
        ("CTC", 0.202),
        ("CTA", 0.069),
        ("CTG", 0.403),
    ],
    # Isoleucine
    "I": [
        ("ATT", 0.342),
        ("ATC", 0.510),
        ("ATA", 0.148),
    ],
    # Methionine (start)
    "M": [
        ("ATG", 1.000),
    ],
    # Valine
    "V": [
        ("GTT", 0.177),
        ("GTC", 0.216),
        ("GTA", 0.108),
        ("GTG", 0.499),
    ],
    # Serine
    "S": [
        ("TCT", 0.182),
        ("TCC", 0.224),
        ("TCA", 0.144),
        ("TCG", 0.060),
        ("AGT", 0.147),
        ("AGC", 0.243),
    ],
    # Proline
    "P": [
        ("CCT", 0.289),
        ("CCC", 0.324),
        ("CCA", 0.272),
        ("CCG", 0.115),
    ],
    # Threonine
    "T": [
        ("ACT", 0.240),
        ("ACC", 0.367),
        ("ACA", 0.274),
        ("ACG", 0.119),
    ],
    # Alanine
    "A": [
        ("GCT", 0.264),
        ("GCC", 0.408),
        ("GCA", 0.222),
        ("GCG", 0.106),
    ],
    # Tyrosine
    "Y": [
        ("TAT", 0.422),
        ("TAC", 0.578),
    ],
    # Histidine
    "H": [
        ("CAT", 0.413),
        ("CAC", 0.587),
    ],
    # Glutamine
    "Q": [
        ("CAA", 0.260),
        ("CAG", 0.740),
    ],
    # Asparagine
    "N": [
        ("AAT", 0.451),
        ("AAC", 0.549),
    ],
    # Lysine
    "K": [
        ("AAA", 0.417),
        ("AAG", 0.583),
    ],
    # Aspartic acid
    "D": [
        ("GAT", 0.448),
        ("GAC", 0.552),
    ],
    # Glutamic acid
    "E": [
        ("GAA", 0.410),
        ("GAG", 0.590),
    ],
    # Cysteine
    "C": [
        ("TGT", 0.440),
        ("TGC", 0.560),
    ],
    # Tryptophan
    "W": [
        ("TGG", 1.000),
    ],
    # Arginine
    "R": [
        ("CGT", 0.083),
        ("CGC", 0.191),
        ("CGA", 0.108),
        ("CGG", 0.206),
        ("AGA", 0.202),
        ("AGG", 0.210),
    ],
    # Glycine
    "G": [
        ("GGT", 0.159),
        ("GGC", 0.345),
        ("GGA", 0.246),
        ("GGG", 0.250),
    ],
    # Stop codons
    "*": [
        ("TAA", 0.280),
        ("TAG", 0.200),
        ("TGA", 0.520),
    ],
}


# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------


def get_optimal_codon(amino_acid: str) -> str:
    """Return the most frequent codon for the given amino acid.

    Args:
        amino_acid: Single-letter amino acid code (upper-case), or ``"*"``
                    for stop.

    Returns:
        The codon (3-letter DNA string) with the highest usage frequency.

    Raises:
        KeyError: If *amino_acid* is not found in the codon table.
    """
    aa = amino_acid.upper()
    if aa not in CODON_TABLE:
        raise KeyError(f"Unknown amino acid: {amino_acid!r}")
    return max(CODON_TABLE[aa], key=lambda pair: pair[1])[0]


def reverse_translate(protein_seq: str) -> str:
    """Reverse-translate a protein sequence using optimal dog codons.

    Each amino acid is replaced by its highest-frequency codon from the
    *Canis lupus familiaris* codon table.  A stop codon (``*``) is **not**
    appended automatically -- include it in the input if desired.

    Args:
        protein_seq: Amino-acid sequence (single-letter codes, upper-case).

    Returns:
        A DNA string exactly 3x the length of *protein_seq*.

    Raises:
        KeyError: If any residue is not in the codon table.
    """
    return "".join(get_optimal_codon(aa) for aa in protein_seq.upper())
