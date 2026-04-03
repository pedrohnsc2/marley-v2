"""Design an mRNA vaccine construct from top-scored antigen candidates.

Selects CTL epitopes via IEDB MHC-I binding predictions, assembles a
multi-epitope protein construct with signal peptide, adjuvant, and linkers,
performs codon optimization for *Canis lupus familiaris*, and builds the
full mRNA cassette.  Outputs FASTA files, an identity card (JSON), and a
Markdown design-rationale report.

Usage:
    python -m pipeline.06_construct
    python -m pipeline.06_construct --force --signal-peptide IgK --adjuvant RS09
"""

from __future__ import annotations

import csv
import io
import json
import re
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import requests
from tqdm import tqdm

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Try importing Biopython's ProteinAnalysis (optional dependency)
# ---------------------------------------------------------------------------

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as _ProteinAnalysis

    HAS_BIOPYTHON = True
except ImportError:  # pragma: no cover
    HAS_BIOPYTHON = False

# ---------------------------------------------------------------------------
# Try importing the codon table created by the other agent
# ---------------------------------------------------------------------------

try:
    from core.codon_tables import CANIS_LUPUS_CODON_TABLE, get_optimal_codon
except ImportError:
    # Fallback: minimal Canis lupus familiaris codon preference table.
    # Each amino acid maps to a list of (codon, relative_frequency) tuples
    # sorted by descending frequency.  Only the top two codons per AA are
    # included; the first entry is used by default.
    CANIS_LUPUS_CODON_TABLE: dict[str, list[tuple[str, float]]] = {
        "A": [("GCC", 0.40), ("GCT", 0.26), ("GCA", 0.23), ("GCG", 0.11)],
        "R": [("CGG", 0.21), ("AGG", 0.20), ("AGA", 0.20), ("CGC", 0.19), ("CGT", 0.11), ("CGA", 0.09)],
        "N": [("AAC", 0.54), ("AAT", 0.46)],
        "D": [("GAC", 0.54), ("GAT", 0.46)],
        "C": [("TGC", 0.55), ("TGT", 0.45)],
        "Q": [("CAG", 0.73), ("CAA", 0.27)],
        "E": [("GAG", 0.58), ("GAA", 0.42)],
        "G": [("GGC", 0.34), ("GGG", 0.25), ("GGA", 0.25), ("GGT", 0.16)],
        "H": [("CAC", 0.58), ("CAT", 0.42)],
        "I": [("ATC", 0.48), ("ATT", 0.36), ("ATA", 0.16)],
        "L": [("CTG", 0.41), ("CTC", 0.20), ("CTT", 0.13), ("TTG", 0.13), ("TTA", 0.07), ("CTA", 0.07)],
        "K": [("AAG", 0.58), ("AAA", 0.42)],
        "M": [("ATG", 1.00)],
        "F": [("TTC", 0.55), ("TTT", 0.45)],
        "P": [("CCC", 0.33), ("CCT", 0.28), ("CCA", 0.27), ("CCG", 0.12)],
        "S": [("AGC", 0.24), ("TCC", 0.22), ("TCT", 0.18), ("TCA", 0.15), ("AGT", 0.12), ("TCG", 0.09)],
        "T": [("ACC", 0.36), ("ACA", 0.28), ("ACT", 0.24), ("ACG", 0.12)],
        "W": [("TGG", 1.00)],
        "Y": [("TAC", 0.56), ("TAT", 0.44)],
        "V": [("GTG", 0.46), ("GTC", 0.24), ("GTT", 0.18), ("GTA", 0.12)],
        "*": [("TGA", 0.47), ("TAA", 0.28), ("TAG", 0.25)],
    }

    def get_optimal_codon(amino_acid: str) -> str:
        """Return the most-preferred codon for *amino_acid* in *Canis lupus*."""
        codons = CANIS_LUPUS_CODON_TABLE.get(amino_acid)
        if not codons:
            raise ValueError(f"Unknown amino acid: {amino_acid!r}")
        return codons[0][0]


logger = get_logger("construct")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INPUT_FILE: str = "results/scored_candidates.csv"
FASTA_FILE: str = "data/raw/surface_proteins.fasta"
OUTPUT_DIR: str = "results/construct"

# Epitope selection
TOP_CANDIDATES: int = 10
TARGET_CTL_COUNT: int = 15
IC50_THRESHOLD: float = 500.0
OVERLAP_THRESHOLD: int = 5
MAX_EPITOPES_PER_GENE: int = 3

# IEDB MHC-I (CTL)
IEDB_API_URL: str = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHC_ALLELES: list[str] = ["DLA-8850101", "DLA-8850801", "DLA-8803401"]
PEPTIDE_LENGTH: int = 9
PREDICTION_METHOD: str = "netmhcpan_ba"
API_DELAY_SECONDS: float = 1.0

# IEDB MHC-II (HTL) — canine DLA-II not supported by IEDB, using HLA proxy
IEDB_MHCII_API_URL: str = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
MHCII_ALLELES: list[str] = ["HLA-DRB1*01:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01"]
MHCII_PEPTIDE_LENGTH: int = 15
MHCII_METHOD: str = "nn_align"
MHCII_IC50_THRESHOLD: float = 1000.0  # MHC-II uses relaxed threshold
TARGET_HTL_COUNT: int = 10

# Construct components
SIGNAL_PEPTIDES: dict[str, str] = {
    "tPA": "MDAMKRGLCCVLLLCGAVFVSAS",
    "IgK": "METDTLLLWVLLLWVPGSTGD",
}
ADJUVANTS: dict[str, str] = {
    "L7L12": (
        "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAEEQ"
        "SEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAKAKL"
        "EAAGATVTVK"
    ),
    "RS09": "APPHALS",
}
LINKERS: dict[str, str] = {
    "adjuvant_to_epitopes": "EAAAK",
    "ctl": "AAY",
    "htl": "GPGPG",
}

# mRNA components
FIVE_PRIME_UTR: str = (
    "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGACCCCGGCGCCGCCACC"
)
THREE_PRIME_UTR: str = (
    "GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAAC"
    "TACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAA"
    "AACATTTATTTTCATTGC"
)
POLY_A_LENGTH: int = 120
STOP_CODON: str = "TGA"

# Safety
VAXIJEN_URL: str = "http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html"
VAXIJEN_THRESHOLD: float = 0.4
ALLERTOP_URL: str = "https://www.ddg-pharmfac.net/AllerTOP/"

# Restriction sites to avoid in codon optimization
RESTRICTION_SITES: list[str] = [
    "GAATTC",  # EcoRI
    "GGATCC",  # BamHI
    "AAGCTT",  # HindIII
    "TCTAGA",  # XbaI
    "GCTAGC",  # NheI
    "GGTCTC",  # BsaI
]


# ---------------------------------------------------------------------------
# Lightweight dataclasses (not in core/models yet)
# ---------------------------------------------------------------------------


@dataclass
class SelectedEpitope:
    """A single epitope selected for the vaccine construct."""

    peptide: str
    gene_id: str
    gene_name: str
    allele: str
    ic50: float
    rank: float
    position: int = 0
    epitope_type: str = "CTL"  # "CTL" or "HTL"


@dataclass
class ConstructCard:
    """Identity card summarising the designed vaccine construct."""

    signal_peptide_name: str
    adjuvant_name: str
    epitope_count: int
    protein_length: int
    mrna_length: int
    gc_content: float
    molecular_weight: float = 0.0
    isoelectric_point: float = 0.0
    instability_index: float = 0.0
    gravy: float = 0.0
    aromaticity: float = 0.0
    vaxijen_score: float | None = None
    allertop_result: str | None = None
    restriction_sites_removed: int = 0
    epitopes: list[dict] = field(default_factory=list)


# ---------------------------------------------------------------------------
# FASTA parsing (mirrors module 04 helper)
# ---------------------------------------------------------------------------


def _load_sequences_from_fasta(fasta_path: str) -> dict[str, str]:
    """Parse a FASTA file into a mapping of gene_id to amino-acid sequence.

    The gene_id is extracted as the first whitespace-delimited token after
    the '>' character on each header line.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dictionary mapping gene identifiers to their full sequences.

    Raises:
        FileNotFoundError: If *fasta_path* does not exist.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(parts)

    return sequences


# ---------------------------------------------------------------------------
# IEDB API interaction (reuses pattern from module 04)
# ---------------------------------------------------------------------------


def _predict_binding(sequence: str, allele: str) -> list[dict]:
    """Submit a sequence to the IEDB MHC-I binding prediction API.

    Args:
        sequence: Amino-acid sequence to analyse.
        allele: MHC allele name (e.g. ``"DLA-8850101"``).

    Returns:
        List of dicts with keys ``peptide``, ``allele``, ``ic50``, ``rank``.
        Returns an empty list when the API call fails.
    """
    payload = {
        "method": PREDICTION_METHOD,
        "sequence_text": sequence,
        "allele": allele,
        "length": str(PEPTIDE_LENGTH),
    }

    try:
        response = requests.post(IEDB_API_URL, data=payload, timeout=120)
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("IEDB API request failed for allele %s: %s", allele, exc)
        return []

    return _parse_iedb_response(response.text, allele)


def _parse_iedb_response(text: str, allele: str) -> list[dict]:
    """Parse tab-separated IEDB MHC-I prediction output.

    Args:
        text: Raw response body from the IEDB API.
        allele: Allele string to attach to each record as fallback.

    Returns:
        Parsed list of prediction dicts.
    """
    results: list[dict] = []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")

    for row in reader:
        try:
            peptide = row.get("peptide", "").strip()
            ic50_raw = row.get("ic50", row.get("ic50_score", "")).strip()
            rank_raw = row.get("percentile_rank", row.get("rank", "0")).strip()

            if not peptide or not ic50_raw:
                continue

            results.append(
                {
                    "peptide": peptide,
                    "allele": row.get("allele", allele).strip(),
                    "ic50": float(ic50_raw),
                    "rank": float(rank_raw) if rank_raw else 0.0,
                }
            )
        except (ValueError, KeyError) as exc:
            logger.debug("Skipping unparsable IEDB row: %s (%s)", row, exc)
            continue

    return results


# ---------------------------------------------------------------------------
# 1. Epitope selection
# ---------------------------------------------------------------------------


def _peptides_overlap(pep_a: str, pep_b: str, seq: str, threshold: int = OVERLAP_THRESHOLD) -> bool:
    """Check whether two peptides overlap by at least *threshold* residues
    within the same parent *seq*.

    If either peptide is not found in *seq*, returns ``False``.
    """
    pos_a = seq.find(pep_a)
    pos_b = seq.find(pep_b)
    if pos_a == -1 or pos_b == -1:
        return False

    start = max(pos_a, pos_b)
    end = min(pos_a + len(pep_a), pos_b + len(pep_b))
    return (end - start) >= threshold


def select_epitopes(
    candidates: list[dict],
    sequences: dict[str, str],
    max_ctl: int = TARGET_CTL_COUNT,
) -> list[SelectedEpitope]:
    """Select CTL epitopes from top candidates via IEDB MHC-I predictions.

    For each of the top-ranked candidates (sorted by ``final_score``),
    queries the IEDB API for all canine DLA alleles, collects peptides
    with IC50 below the threshold, and greedily selects non-overlapping
    epitopes up to *max_ctl*.

    Args:
        candidates: Rows from ``scored_candidates.csv`` as dicts.
        sequences: Mapping of gene_id to protein sequence.
        max_ctl: Maximum number of CTL epitopes to select.

    Returns:
        List of :class:`SelectedEpitope` instances, sorted by IC50.
    """
    # Sort candidates by final_score descending, take top N.
    sorted_candidates = sorted(
        candidates,
        key=lambda c: float(c.get("final_score", 0)),
        reverse=True,
    )[:TOP_CANDIDATES]

    logger.info(
        "Selecting CTL epitopes from top %d candidates ...",
        len(sorted_candidates),
    )

    # Collect all strong binders across candidates and alleles.
    all_hits: list[dict] = []

    for cand in tqdm(sorted_candidates, desc="IEDB CTL prediction", unit="gene"):
        gene_id = cand.get("gene_id", "")
        gene_name = cand.get("gene_name", "")
        sequence = sequences.get(gene_id, cand.get("sequence", ""))

        if not sequence or len(sequence) < PEPTIDE_LENGTH:
            logger.warning("No usable sequence for %s. Skipping.", gene_id)
            continue

        for allele in MHC_ALLELES:
            predictions = _predict_binding(sequence, allele)
            time.sleep(API_DELAY_SECONDS)

            for pred in predictions:
                if pred["ic50"] < IC50_THRESHOLD:
                    all_hits.append(
                        {
                            "peptide": pred["peptide"],
                            "gene_id": gene_id,
                            "gene_name": gene_name,
                            "allele": pred["allele"],
                            "ic50": pred["ic50"],
                            "rank": pred["rank"],
                            "sequence": sequence,
                        }
                    )

    # Sort by IC50 ascending (strongest binders first).
    all_hits.sort(key=lambda h: h["ic50"])

    logger.info("Found %d peptides with IC50 < %.0f nM.", len(all_hits), IC50_THRESHOLD)

    # Greedy selection: non-overlapping, max per gene.
    selected: list[SelectedEpitope] = []
    gene_counts: dict[str, int] = {}

    for hit in all_hits:
        if len(selected) >= max_ctl:
            break

        gene_id = hit["gene_id"]
        peptide = hit["peptide"]
        sequence = hit["sequence"]

        # Enforce per-gene cap.
        if gene_counts.get(gene_id, 0) >= MAX_EPITOPES_PER_GENE:
            continue

        # Check overlap with already-selected epitopes from the same gene.
        dominated = False
        for existing in selected:
            if existing.gene_id == gene_id:
                if _peptides_overlap(peptide, existing.peptide, sequence):
                    dominated = True
                    break
        if dominated:
            continue

        position = sequence.find(peptide)
        selected.append(
            SelectedEpitope(
                peptide=peptide,
                gene_id=gene_id,
                gene_name=hit["gene_name"],
                allele=hit["allele"],
                ic50=hit["ic50"],
                rank=hit["rank"],
                position=position,
            )
        )
        gene_counts[gene_id] = gene_counts.get(gene_id, 0) + 1

    logger.info("Selected %d CTL epitopes.", len(selected))
    return selected


# ---------------------------------------------------------------------------
# 1b. MHC-II (HTL) epitope selection
# ---------------------------------------------------------------------------


def _predict_mhcii_binding(sequence: str, allele: str) -> list[dict]:
    """Submit a sequence to the IEDB MHC-II binding prediction API.

    Uses nn_align method against HLA-DRB1 alleles as cross-species proxy
    for canine DLA class II (which IEDB does not support directly).

    Args:
        sequence: Amino-acid sequence to analyse.
        allele: MHC-II allele name (e.g. ``"HLA-DRB1*01:01"``).

    Returns:
        List of dicts with keys ``peptide``, ``allele``, ``ic50``, ``rank``.
        Returns an empty list when the API call fails.
    """
    payload = {
        "method": MHCII_METHOD,
        "sequence_text": sequence,
        "allele": allele,
        "length": str(MHCII_PEPTIDE_LENGTH),
    }

    try:
        response = requests.post(IEDB_MHCII_API_URL, data=payload, timeout=120)
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("IEDB MHC-II API request failed for allele %s: %s", allele, exc)
        return []

    return _parse_iedb_response(response.text, allele)


def select_htl_epitopes(
    candidates: list[dict],
    sequences: dict[str, str],
    max_htl: int = TARGET_HTL_COUNT,
) -> list[SelectedEpitope]:
    """Select HTL epitopes from top candidates via IEDB MHC-II predictions.

    Mirrors :func:`select_epitopes` but targets the MHC-II endpoint with
    HLA-DRB1 alleles, 15-mer peptides, and a relaxed IC50 threshold of
    1000 nM.

    Args:
        candidates: Rows from ``scored_candidates.csv`` as dicts.
        sequences: Mapping of gene_id to protein sequence.
        max_htl: Maximum number of HTL epitopes to select.

    Returns:
        List of :class:`SelectedEpitope` instances with ``epitope_type="HTL"``,
        sorted by IC50.
    """
    sorted_candidates = sorted(
        candidates,
        key=lambda c: float(c.get("final_score", 0)),
        reverse=True,
    )[:TOP_CANDIDATES]

    logger.info(
        "Selecting HTL epitopes from top %d candidates ...",
        len(sorted_candidates),
    )

    all_hits: list[dict] = []

    for cand in tqdm(sorted_candidates, desc="IEDB HTL prediction", unit="gene"):
        gene_id = cand.get("gene_id", "")
        gene_name = cand.get("gene_name", "")
        sequence = sequences.get(gene_id, cand.get("sequence", ""))

        if not sequence or len(sequence) < MHCII_PEPTIDE_LENGTH:
            logger.warning("No usable sequence for %s. Skipping.", gene_id)
            continue

        for allele in MHCII_ALLELES:
            predictions = _predict_mhcii_binding(sequence, allele)
            time.sleep(API_DELAY_SECONDS)

            for pred in predictions:
                if pred["ic50"] < MHCII_IC50_THRESHOLD:
                    all_hits.append(
                        {
                            "peptide": pred["peptide"],
                            "gene_id": gene_id,
                            "gene_name": gene_name,
                            "allele": pred["allele"],
                            "ic50": pred["ic50"],
                            "rank": pred["rank"],
                            "sequence": sequence,
                        }
                    )

    all_hits.sort(key=lambda h: h["ic50"])

    logger.info(
        "Found %d HTL peptides with IC50 < %.0f nM.", len(all_hits), MHCII_IC50_THRESHOLD,
    )

    # Greedy selection: non-overlapping, max per gene.
    selected: list[SelectedEpitope] = []
    gene_counts: dict[str, int] = {}

    for hit in all_hits:
        if len(selected) >= max_htl:
            break

        gene_id = hit["gene_id"]
        peptide = hit["peptide"]
        sequence = hit["sequence"]

        if gene_counts.get(gene_id, 0) >= MAX_EPITOPES_PER_GENE:
            continue

        dominated = False
        for existing in selected:
            if existing.gene_id == gene_id:
                if _peptides_overlap(peptide, existing.peptide, sequence):
                    dominated = True
                    break
        if dominated:
            continue

        position = sequence.find(peptide)
        selected.append(
            SelectedEpitope(
                peptide=peptide,
                gene_id=gene_id,
                gene_name=hit["gene_name"],
                allele=hit["allele"],
                ic50=hit["ic50"],
                rank=hit["rank"],
                position=position,
                epitope_type="HTL",
            )
        )
        gene_counts[gene_id] = gene_counts.get(gene_id, 0) + 1

    logger.info("Selected %d HTL epitopes.", len(selected))
    return selected


# ---------------------------------------------------------------------------
# 2. Construct assembly
# ---------------------------------------------------------------------------


def assemble_construct(
    epitopes: list[SelectedEpitope],
    signal_peptide: str = "tPA",
    adjuvant: str = "L7L12",
    htl_epitopes: list[SelectedEpitope] | None = None,
) -> str:
    """Assemble the multi-epitope vaccine protein construct.

    Layout::

        [Signal peptide] + [Adjuvant] + EAAAK +
        [CTL1]-AAY-[CTL2]-AAY-...-[CTLn] +
        GPGPG + [HTL1]-GPGPG-[HTL2]-GPGPG-...-[HTLm]

    If *htl_epitopes* is ``None`` or empty, the HTL block is omitted and
    the construct is identical to the previous CTL-only layout (backwards
    compatible).

    Args:
        epitopes: Selected CTL epitopes.
        signal_peptide: Key into :data:`SIGNAL_PEPTIDES`.
        adjuvant: Key into :data:`ADJUVANTS`.
        htl_epitopes: Optional list of selected HTL epitopes.

    Returns:
        Full amino-acid sequence of the vaccine construct.

    Raises:
        ValueError: If *signal_peptide* or *adjuvant* is not recognised.
    """
    if signal_peptide not in SIGNAL_PEPTIDES:
        raise ValueError(
            f"Unknown signal peptide: {signal_peptide!r}. "
            f"Choose from {list(SIGNAL_PEPTIDES.keys())}."
        )
    if adjuvant not in ADJUVANTS:
        raise ValueError(
            f"Unknown adjuvant: {adjuvant!r}. "
            f"Choose from {list(ADJUVANTS.keys())}."
        )

    parts: list[str] = [
        SIGNAL_PEPTIDES[signal_peptide],
        ADJUVANTS[adjuvant],
        LINKERS["adjuvant_to_epitopes"],
    ]

    # CTL epitope cassette
    for i, ep in enumerate(epitopes):
        if i > 0:
            parts.append(LINKERS["ctl"])
        parts.append(ep.peptide)

    # HTL epitope cassette (appended after a GPGPG bridge)
    if htl_epitopes:
        parts.append(LINKERS["htl"])
        for i, ep in enumerate(htl_epitopes):
            if i > 0:
                parts.append(LINKERS["htl"])
            parts.append(ep.peptide)

    construct = "".join(parts)
    htl_count = len(htl_epitopes) if htl_epitopes else 0
    logger.info(
        "Assembled construct: %d aa (%s signal + %s adjuvant + %d CTL + %d HTL epitopes).",
        len(construct),
        signal_peptide,
        adjuvant,
        len(epitopes),
        htl_count,
    )
    return construct


# ---------------------------------------------------------------------------
# 3. Codon optimization
# ---------------------------------------------------------------------------


def _get_alternative_codon(amino_acid: str, avoid: str) -> str:
    """Return an alternative codon for *amino_acid* that is not *avoid*.

    Falls back to *avoid* itself if no alternative exists.
    """
    codons = CANIS_LUPUS_CODON_TABLE.get(amino_acid, [])
    for codon, _freq in codons:
        if codon != avoid:
            return codon
    return avoid


def reverse_translate_optimized(protein_seq: str) -> str:
    """Reverse-translate a protein sequence using Canis lupus codon preferences.

    Post-processing steps:
    1. Remove restriction enzyme recognition sites by swapping to the
       next-best synonymous codon.
    2. Break homopolymer runs longer than 4 nt to improve mRNA stability.

    Args:
        protein_seq: Amino-acid sequence (single-letter codes).

    Returns:
        Codon-optimised DNA coding sequence (without stop codon).
    """
    # --- Initial translation -----------------------------------------------
    codons: list[str] = []
    for aa in protein_seq:
        codons.append(get_optimal_codon(aa))

    # --- Remove restriction sites ------------------------------------------
    sites_removed = 0
    for _pass in range(3):  # iterate a few times in case swaps create new sites
        dna = "".join(codons)
        found_any = False
        for site in RESTRICTION_SITES:
            idx = dna.find(site)
            while idx != -1:
                found_any = True
                sites_removed += 1
                # Find which codon(s) overlap with the site and swap one.
                codon_idx = idx // 3
                if codon_idx < len(codons):
                    aa = protein_seq[codon_idx] if codon_idx < len(protein_seq) else None
                    if aa:
                        codons[codon_idx] = _get_alternative_codon(aa, codons[codon_idx])
                dna = "".join(codons)
                idx = dna.find(site, idx + 1)
        if not found_any:
            break

    # --- Break homopolymer runs (> 4 identical nucleotides) ----------------
    dna = "".join(codons)
    homopolymer_pattern = re.compile(r"([ATCG])\1{4,}")
    for match in homopolymer_pattern.finditer(dna):
        start = match.start()
        codon_idx = start // 3
        if codon_idx < len(codons) and codon_idx < len(protein_seq):
            aa = protein_seq[codon_idx]
            codons[codon_idx] = _get_alternative_codon(aa, codons[codon_idx])

    cds = "".join(codons)
    logger.info(
        "Codon optimisation complete: %d bp CDS, %d restriction site(s) removed.",
        len(cds),
        sites_removed,
    )
    return cds


# ---------------------------------------------------------------------------
# 4. mRNA assembly
# ---------------------------------------------------------------------------


def assemble_mrna(cds: str) -> str:
    """Assemble the full mRNA cassette around the codon-optimised CDS.

    Layout (BioNTech-style double 3'UTR)::

        5'UTR + CDS + stop codon + 3'UTR + 3'UTR + poly(A)

    Args:
        cds: Codon-optimised DNA coding sequence.

    Returns:
        Complete mRNA sequence as a DNA-alphabet string.
    """
    poly_a = "A" * POLY_A_LENGTH
    mrna = FIVE_PRIME_UTR + cds + STOP_CODON + THREE_PRIME_UTR + THREE_PRIME_UTR + poly_a
    logger.info("mRNA cassette assembled: %d nt total.", len(mrna))
    return mrna


# ---------------------------------------------------------------------------
# 5. GC content
# ---------------------------------------------------------------------------


def calculate_gc_content(dna_seq: str) -> float:
    """Calculate the GC content of a DNA sequence.

    Args:
        dna_seq: DNA sequence string (case-insensitive).

    Returns:
        GC fraction in [0.0, 1.0].  Returns 0.0 for empty input.
    """
    if not dna_seq:
        return 0.0
    dna_upper = dna_seq.upper()
    gc = sum(1 for nt in dna_upper if nt in ("G", "C"))
    return gc / len(dna_upper)


# ---------------------------------------------------------------------------
# 6. Physicochemical analysis
# ---------------------------------------------------------------------------


def compute_physicochemical(protein_seq: str) -> dict:
    """Compute physicochemical properties of the vaccine protein construct.

    Uses Biopython's ``ProteinAnalysis`` when available.  Returns a dict
    with zero-valued defaults if Biopython is not installed.

    Args:
        protein_seq: Full amino-acid sequence.

    Returns:
        Dictionary with keys: ``molecular_weight``, ``isoelectric_point``,
        ``instability_index``, ``gravy``, ``aromaticity``.
    """
    defaults: dict[str, float] = {
        "molecular_weight": 0.0,
        "isoelectric_point": 0.0,
        "instability_index": 0.0,
        "gravy": 0.0,
        "aromaticity": 0.0,
    }

    if not HAS_BIOPYTHON:
        logger.warning(
            "Biopython not available. Physicochemical properties will be zeroed."
        )
        return defaults

    try:
        # ProteinAnalysis cannot handle 'X' or other ambiguous residues.
        clean_seq = protein_seq.replace("X", "").replace("*", "")
        analysis = _ProteinAnalysis(clean_seq)
        return {
            "molecular_weight": round(analysis.molecular_weight(), 2),
            "isoelectric_point": round(analysis.isoelectric_point(), 2),
            "instability_index": round(analysis.instability_index(), 2),
            "gravy": round(analysis.gravy(), 4),
            "aromaticity": round(analysis.aromaticity(), 4),
        }
    except Exception as exc:
        logger.warning("Physicochemical analysis failed: %s", exc)
        return defaults


# ---------------------------------------------------------------------------
# 7. VaxiJen antigenicity prediction
# ---------------------------------------------------------------------------


def predict_antigenicity(protein_seq: str) -> float | None:
    """Predict antigenicity of the construct via the VaxiJen web server.

    Submits the protein sequence to the VaxiJen HTML form and parses the
    result page for the predicted score.

    Args:
        protein_seq: Amino-acid sequence.

    Returns:
        VaxiJen score (float), or ``None`` if the request fails.
    """
    try:
        payload = {
            "Ession": protein_seq,
            "Target": "0",  # bacteria (closest to parasite model)
            "thresh": str(VAXIJEN_THRESHOLD),
        }
        response = requests.post(
            VAXIJEN_URL,
            data=payload,
            timeout=30,
        )
        response.raise_for_status()

        # VaxiJen output contains a line like:
        #   "Overall Prediction for the Protective Antigen = 0.7531"
        match = re.search(
            r"Overall Prediction[^=]*=\s*([\d.]+)",
            response.text,
        )
        if match:
            score = float(match.group(1))
            logger.info("VaxiJen antigenicity score: %.4f", score)
            return score

        logger.warning("Could not parse VaxiJen score from response.")
        return None

    except requests.RequestException as exc:
        logger.warning("VaxiJen request failed: %s", exc)
        return None
    except (ValueError, AttributeError) as exc:
        logger.warning("VaxiJen parsing error: %s", exc)
        return None


# ---------------------------------------------------------------------------
# 8. AllerTOP allergenicity prediction
# ---------------------------------------------------------------------------


def predict_allergenicity(protein_seq: str) -> str | None:
    """Predict allergenicity of the construct via the AllerTOP web server.

    Args:
        protein_seq: Amino-acid sequence.

    Returns:
        Classification string (e.g. ``"PROBABLE NON-ALLERGEN"``), or
        ``None`` if the request fails.
    """
    try:
        payload = {"Ession": protein_seq}
        response = requests.post(
            ALLERTOP_URL,
            data=payload,
            timeout=30,
        )
        response.raise_for_status()

        # AllerTOP result page contains classification text like:
        #   "PROBABLE NON-ALLERGEN" or "PROBABLE ALLERGEN"
        for pattern in [
            r"(PROBABLE\s+NON-ALLERGEN)",
            r"(PROBABLE\s+ALLERGEN)",
        ]:
            match = re.search(pattern, response.text, re.IGNORECASE)
            if match:
                result = match.group(1).upper().strip()
                logger.info("AllerTOP result: %s", result)
                return result

        logger.warning("Could not parse AllerTOP result from response.")
        return None

    except requests.RequestException as exc:
        logger.warning("AllerTOP request failed: %s", exc)
        return None
    except (ValueError, AttributeError) as exc:
        logger.warning("AllerTOP parsing error: %s", exc)
        return None


# ---------------------------------------------------------------------------
# 9. Output generation
# ---------------------------------------------------------------------------


def write_outputs(
    construct: str,
    mrna: str,
    epitopes: list[SelectedEpitope],
    physicochemical: dict,
    card: ConstructCard,
    output_dir: str = OUTPUT_DIR,
) -> None:
    """Write all output artefacts to *output_dir*.

    Produces:
    - ``vaccine_construct.fasta`` -- protein sequence
    - ``vaccine_mrna.fasta`` -- mRNA cassette
    - ``construct_card.json`` -- machine-readable identity card
    - ``construct_report.md`` -- human-readable design rationale

    Args:
        construct: Protein sequence of the vaccine construct.
        mrna: Full mRNA cassette sequence.
        epitopes: Selected epitopes included in the construct.
        physicochemical: Dict of physicochemical properties.
        card: Populated :class:`ConstructCard`.
        output_dir: Destination directory.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # -- Protein FASTA ------------------------------------------------------
    fasta_protein = out / "vaccine_construct.fasta"
    with open(fasta_protein, "w") as fh:
        fh.write(
            f">Marley_vaccine_construct len={len(construct)} "
            f"signal={card.signal_peptide_name} adj={card.adjuvant_name}\n"
        )
        # Wrap at 70 characters.
        for i in range(0, len(construct), 70):
            fh.write(construct[i : i + 70] + "\n")
    logger.info("Protein FASTA written to %s.", fasta_protein)

    # -- mRNA FASTA ---------------------------------------------------------
    fasta_mrna = out / "vaccine_mrna.fasta"
    with open(fasta_mrna, "w") as fh:
        fh.write(
            f">Marley_vaccine_mRNA len={len(mrna)} "
            f"GC={card.gc_content:.1%}\n"
        )
        for i in range(0, len(mrna), 70):
            fh.write(mrna[i : i + 70] + "\n")
    logger.info("mRNA FASTA written to %s.", fasta_mrna)

    # -- Identity card (JSON) -----------------------------------------------
    card_path = out / "construct_card.json"
    card_dict = {
        "generated": datetime.now().isoformat(),
        "signal_peptide": card.signal_peptide_name,
        "adjuvant": card.adjuvant_name,
        "epitope_count": card.epitope_count,
        "protein_length_aa": card.protein_length,
        "mrna_length_nt": card.mrna_length,
        "gc_content": round(card.gc_content, 4),
        "molecular_weight_da": card.molecular_weight,
        "isoelectric_point": card.isoelectric_point,
        "instability_index": card.instability_index,
        "gravy": card.gravy,
        "aromaticity": card.aromaticity,
        "vaxijen_score": card.vaxijen_score,
        "vaxijen_threshold": VAXIJEN_THRESHOLD,
        "allertop_result": card.allertop_result,
        "restriction_sites_removed": card.restriction_sites_removed,
        "epitopes": card.epitopes,
    }
    with open(card_path, "w") as fh:
        json.dump(card_dict, fh, indent=2)
    logger.info("Identity card written to %s.", card_path)

    # -- Markdown report ----------------------------------------------------
    report_path = out / "construct_report.md"
    report = _build_markdown_report(construct, mrna, epitopes, physicochemical, card)
    with open(report_path, "w") as fh:
        fh.write(report)
    logger.info("Design report written to %s.", report_path)


def _build_markdown_report(
    construct: str,
    mrna: str,
    epitopes: list[SelectedEpitope],
    physicochemical: dict,
    card: ConstructCard,
    htl_epitopes: list[SelectedEpitope] | None = None,
) -> str:
    """Build a Markdown report describing the vaccine construct design.

    Returns:
        Full Markdown text.
    """
    today = datetime.now().strftime("%Y-%m-%d")
    lines: list[str] = []

    # Title
    lines.append("# Marley -- mRNA Vaccine Construct Design Report")
    lines.append("")
    lines.append(f"**Generated:** {today}")
    lines.append("")

    # Construct overview
    lines.append("## Construct Overview")
    lines.append("")
    lines.append(f"| Property | Value |")
    lines.append(f"|---|---|")
    lines.append(f"| Signal peptide | {card.signal_peptide_name} ({len(SIGNAL_PEPTIDES[card.signal_peptide_name])} aa) |")
    lines.append(f"| Adjuvant | {card.adjuvant_name} ({len(ADJUVANTS[card.adjuvant_name])} aa) |")
    htl_count = len(htl_epitopes) if htl_epitopes else 0
    lines.append(f"| CTL epitopes | {len(epitopes)} |")
    lines.append(f"| HTL epitopes | {htl_count} |")
    lines.append(f"| Total epitopes | {card.epitope_count} |")
    lines.append(f"| Linker (adjuvant) | EAAAK |")
    lines.append(f"| Linker (CTL) | AAY |")
    lines.append(f"| Linker (HTL) | GPGPG |")
    lines.append(f"| Total protein length | {card.protein_length} aa |")
    lines.append(f"| Total mRNA length | {card.mrna_length} nt |")
    lines.append("")

    # ASCII diagram
    lines.append("## Construct Architecture")
    lines.append("")
    lines.append("```")
    lines.append(
        "[Signal Peptide]-[Adjuvant]-EAAAK-[Epi1]-AAY-[Epi2]-AAY-...-AAY-[EpiN]"
    )
    lines.append("```")
    lines.append("")
    sig_len = len(SIGNAL_PEPTIDES[card.signal_peptide_name])
    adj_len = len(ADJUVANTS[card.adjuvant_name])
    lines.append(f"- Signal peptide ({card.signal_peptide_name}): residues 1--{sig_len}")
    lines.append(f"- Adjuvant ({card.adjuvant_name}): residues {sig_len + 1}--{sig_len + adj_len}")
    lines.append(f"- EAAAK linker: residues {sig_len + adj_len + 1}--{sig_len + adj_len + 5}")
    lines.append(f"- Epitope cassette: residues {sig_len + adj_len + 6}--{card.protein_length}")
    lines.append("")

    # Epitope table
    lines.append("## Selected CTL Epitopes")
    lines.append("")
    lines.append("| # | Peptide | Source Gene | Allele | IC50 (nM) | Rank |")
    lines.append("|---|---|---|---|---|---|")
    for i, ep in enumerate(epitopes, start=1):
        lines.append(
            f"| {i} | `{ep.peptide}` | {ep.gene_id} ({ep.gene_name}) "
            f"| {ep.allele} | {ep.ic50:.1f} | {ep.rank:.2f} |"
        )
    lines.append("")

    # Physicochemical properties
    lines.append("## Physicochemical Properties")
    lines.append("")
    lines.append("| Property | Value |")
    lines.append("|---|---|")
    lines.append(f"| Molecular weight | {physicochemical.get('molecular_weight', 0.0):.2f} Da |")
    lines.append(f"| Isoelectric point (pI) | {physicochemical.get('isoelectric_point', 0.0):.2f} |")
    lines.append(f"| Instability index | {physicochemical.get('instability_index', 0.0):.2f} |")
    stability = "stable" if physicochemical.get("instability_index", 100) < 40 else "unstable"
    lines.append(f"| Predicted stability | {stability} (threshold < 40) |")
    lines.append(f"| GRAVY | {physicochemical.get('gravy', 0.0):.4f} |")
    lines.append(f"| Aromaticity | {physicochemical.get('aromaticity', 0.0):.4f} |")
    lines.append("")

    # Safety assessment
    lines.append("## Safety Assessment")
    lines.append("")
    if card.vaxijen_score is not None:
        antigen_label = "probable antigen" if card.vaxijen_score >= VAXIJEN_THRESHOLD else "non-antigen"
        lines.append(
            f"- **VaxiJen score:** {card.vaxijen_score:.4f} "
            f"(threshold {VAXIJEN_THRESHOLD}) -- {antigen_label}"
        )
    else:
        lines.append("- **VaxiJen score:** not available (server unreachable)")

    if card.allertop_result is not None:
        lines.append(f"- **AllerTOP:** {card.allertop_result}")
    else:
        lines.append("- **AllerTOP:** not available (server unreachable)")
    lines.append("")

    # Codon optimisation stats
    lines.append("## Codon Optimisation")
    lines.append("")
    lines.append(f"- **Target organism:** *Canis lupus familiaris* (domestic dog)")
    lines.append(f"- **GC content:** {card.gc_content:.1%}")
    lines.append(f"- **Restriction sites removed:** {card.restriction_sites_removed}")
    lines.append(f"- **Avoided sites:** {', '.join(RESTRICTION_SITES)}")
    lines.append("")

    # mRNA structure
    lines.append("## mRNA Cassette Structure")
    lines.append("")
    lines.append("```")
    lines.append("5'cap -- 5'UTR -- CDS -- STOP -- 3'UTR -- 3'UTR -- poly(A)")
    lines.append("```")
    lines.append("")
    lines.append(f"- **5' UTR:** {len(FIVE_PRIME_UTR)} nt (from BNT162b2)")
    lines.append(f"- **CDS:** {card.protein_length * 3} nt")
    lines.append(f"- **Stop codon:** {STOP_CODON}")
    lines.append(f"- **3' UTR:** 2x {len(THREE_PRIME_UTR)} nt (double 3'UTR, BioNTech approach)")
    lines.append(f"- **Poly(A) tail:** {POLY_A_LENGTH} nt")
    lines.append(f"- **Total mRNA:** {card.mrna_length} nt")
    lines.append("")

    # References
    lines.append("## References")
    lines.append("")
    lines.append(
        "1. Patronov A, Doytchinova I. T-cell epitope vaccine design by "
        "immunoinformatics. *Open Biol.* 2013;3(1):120139."
    )
    lines.append(
        "2. Doytchinova IA, Flower DR. VaxiJen: a server for prediction of "
        "protective antigens, tumour antigens and subunit vaccines. "
        "*BMC Bioinformatics.* 2007;8:4."
    )
    lines.append(
        "3. Dimitrov I et al. AllerTOP v.2 -- a server for in silico "
        "prediction of allergens. *J Mol Model.* 2014;20(6):2278."
    )
    lines.append(
        "4. Reis AB et al. Immunity to *Leishmania* and the rational "
        "search for vaccines against canine leishmaniasis. "
        "*Trends Parasitol.* 2010;26(7):341-9."
    )
    lines.append(
        "5. Vogel AB et al. BNT162b vaccines protect rhesus macaques "
        "from SARS-CoV-2. *Nature.* 2021;592:283-289."
    )
    lines.append("")

    # Disclaimer
    lines.append("## Disclaimer")
    lines.append("")
    lines.append(
        "This is a **computational vaccine design** generated by the Marley "
        "reverse vaccinology pipeline. The construct has not been experimentally "
        "validated. All predictions (epitope binding, antigenicity, allergenicity, "
        "physicochemical properties) are *in silico* estimates and must be "
        "confirmed through laboratory testing, including expression, purification, "
        "immunogenicity assays, and animal challenge studies, before any "
        "translational application."
    )
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 10. Main entry point
# ---------------------------------------------------------------------------


def design_construct(
    force: bool = False,
    signal_peptide: str = "tPA",
    adjuvant: str = "L7L12",
) -> str:
    """Run the full mRNA vaccine construct design pipeline.

    Orchestrates epitope selection, construct assembly, codon optimisation,
    mRNA cassette generation, physicochemical analysis, safety prediction,
    and output generation.

    Args:
        force: Re-run even when output directory already contains results.
        signal_peptide: Key into :data:`SIGNAL_PEPTIDES` (``"tPA"`` or ``"IgK"``).
        adjuvant: Key into :data:`ADJUVANTS` (``"L7L12"`` or ``"RS09"``).

    Returns:
        Path to the output directory.
    """
    output_path = Path(OUTPUT_DIR)

    if (output_path / "construct_card.json").exists() and not force:
        logger.info(
            "Output already exists at %s. Use force=True to re-run.",
            OUTPUT_DIR,
        )
        return OUTPUT_DIR

    output_path.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load scored candidates
    # ------------------------------------------------------------------
    input_path = Path(INPUT_FILE)
    if not input_path.exists():
        msg = f"Input file not found: {INPUT_FILE}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    candidates: list[dict] = []
    with open(input_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            candidates.append(row)

    logger.info("Loaded %d scored candidate(s) from %s.", len(candidates), INPUT_FILE)

    if not candidates:
        logger.error("No candidates available. Cannot design construct.")
        raise ValueError("No candidates available for construct design.")

    # ------------------------------------------------------------------
    # 2. Load protein sequences
    # ------------------------------------------------------------------
    sequences: dict[str, str] = {}
    if Path(FASTA_FILE).exists():
        sequences = _load_sequences_from_fasta(FASTA_FILE)
        logger.info("Loaded %d sequence(s) from %s.", len(sequences), FASTA_FILE)
    else:
        logger.warning(
            "FASTA file %s not found. Epitope prediction will be limited.",
            FASTA_FILE,
        )

    # ------------------------------------------------------------------
    # 3. Select CTL epitopes
    # ------------------------------------------------------------------
    epitopes = select_epitopes(candidates, sequences)

    if not epitopes:
        logger.error(
            "No epitopes could be selected. Check IEDB connectivity and "
            "input sequences."
        )
        raise RuntimeError("Epitope selection returned no results.")

    # ------------------------------------------------------------------
    # 4. Assemble protein construct
    # ------------------------------------------------------------------
    construct = assemble_construct(epitopes, signal_peptide, adjuvant)

    # ------------------------------------------------------------------
    # 5. Physicochemical analysis
    # ------------------------------------------------------------------
    logger.info("Computing physicochemical properties ...")
    physicochemical = compute_physicochemical(construct)

    # ------------------------------------------------------------------
    # 6. Codon optimisation
    # ------------------------------------------------------------------
    logger.info("Running codon optimisation for Canis lupus familiaris ...")
    cds = reverse_translate_optimized(construct)
    gc = calculate_gc_content(cds)

    # Count restriction sites that remain (should be zero).
    remaining_sites = sum(
        1 for site in RESTRICTION_SITES if site in cds
    )
    if remaining_sites > 0:
        logger.warning(
            "%d restriction site(s) could not be removed from the CDS.",
            remaining_sites,
        )

    # ------------------------------------------------------------------
    # 7. mRNA assembly
    # ------------------------------------------------------------------
    mrna = assemble_mrna(cds)

    # ------------------------------------------------------------------
    # 8. Safety predictions (non-blocking)
    # ------------------------------------------------------------------
    logger.info("Predicting antigenicity (VaxiJen) ...")
    vaxijen_score = predict_antigenicity(construct)

    logger.info("Predicting allergenicity (AllerTOP) ...")
    allertop_result = predict_allergenicity(construct)

    # ------------------------------------------------------------------
    # 9. Build identity card
    # ------------------------------------------------------------------
    card = ConstructCard(
        signal_peptide_name=signal_peptide,
        adjuvant_name=adjuvant,
        epitope_count=len(epitopes),
        protein_length=len(construct),
        mrna_length=len(mrna),
        gc_content=gc,
        molecular_weight=physicochemical.get("molecular_weight", 0.0),
        isoelectric_point=physicochemical.get("isoelectric_point", 0.0),
        instability_index=physicochemical.get("instability_index", 0.0),
        gravy=physicochemical.get("gravy", 0.0),
        aromaticity=physicochemical.get("aromaticity", 0.0),
        vaxijen_score=vaxijen_score,
        allertop_result=allertop_result,
        restriction_sites_removed=0,  # updated below
        epitopes=[
            {
                "peptide": ep.peptide,
                "gene_id": ep.gene_id,
                "gene_name": ep.gene_name,
                "allele": ep.allele,
                "ic50": ep.ic50,
                "rank": ep.rank,
                "position": ep.position,
            }
            for ep in epitopes
        ],
    )

    # Count how many sites were in the naive translation vs. now.
    naive_cds = "".join(get_optimal_codon(aa) for aa in construct)
    naive_site_count = sum(
        naive_cds.count(site) for site in RESTRICTION_SITES
    )
    final_site_count = sum(cds.count(site) for site in RESTRICTION_SITES)
    card.restriction_sites_removed = max(0, naive_site_count - final_site_count)

    # ------------------------------------------------------------------
    # 10. Write outputs
    # ------------------------------------------------------------------
    write_outputs(construct, mrna, epitopes, physicochemical, card)

    # ------------------------------------------------------------------
    # 11. Summary
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("VACCINE CONSTRUCT DESIGN COMPLETE")
    logger.info("=" * 60)
    logger.info("  Protein length:  %d aa", len(construct))
    logger.info("  mRNA length:     %d nt", len(mrna))
    logger.info("  Epitopes:        %d", len(epitopes))
    logger.info("  GC content:      %.1f%%", gc * 100)
    logger.info("  MW:              %.2f Da", physicochemical.get("molecular_weight", 0.0))
    if vaxijen_score is not None:
        label = "PASS" if vaxijen_score >= VAXIJEN_THRESHOLD else "FAIL"
        logger.info("  VaxiJen:         %.4f (%s)", vaxijen_score, label)
    else:
        logger.info("  VaxiJen:         N/A")
    if allertop_result is not None:
        logger.info("  AllerTOP:        %s", allertop_result)
    else:
        logger.info("  AllerTOP:        N/A")
    logger.info("  Output dir:      %s", OUTPUT_DIR)
    logger.info("=" * 60)

    return OUTPUT_DIR


# ---------------------------------------------------------------------------
# 11. CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Step 6: Design an mRNA vaccine construct from top candidates.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if the output directory already contains results.",
    )
    parser.add_argument(
        "--signal-peptide",
        choices=list(SIGNAL_PEPTIDES.keys()),
        default="tPA",
        help="Signal peptide to use (default: tPA).",
    )
    parser.add_argument(
        "--adjuvant",
        choices=list(ADJUVANTS.keys()),
        default="L7L12",
        help="Adjuvant sequence to use (default: L7L12).",
    )
    args = parser.parse_args()

    result_dir = design_construct(
        force=args.force,
        signal_peptide=args.signal_peptide,
        adjuvant=args.adjuvant,
    )
    print(f"Output: {result_dir}")
