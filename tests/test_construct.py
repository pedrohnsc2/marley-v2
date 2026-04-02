"""Tests for module 06 -- mRNA vaccine construct designer."""

from __future__ import annotations

import importlib
from datetime import datetime
from unittest.mock import patch

import pytest

from core.models import Epitope, VaccineConstruct
from core.codon_tables import (
    CODON_TABLE,
    get_optimal_codon,
    reverse_translate,
)

# The module filename starts with a digit, so we use importlib.
construct_mod = importlib.import_module("pipeline.06_construct")

SelectedEpitope = construct_mod.SelectedEpitope
assemble_construct = construct_mod.assemble_construct
assemble_mrna = construct_mod.assemble_mrna
calculate_gc_content = construct_mod.calculate_gc_content
compute_physicochemical = construct_mod.compute_physicochemical
predict_antigenicity = construct_mod.predict_antigenicity
predict_allergenicity = construct_mod.predict_allergenicity


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _make_epitope(**overrides) -> Epitope:
    """Create an Epitope with sensible defaults, overridden by kwargs."""
    defaults = {
        "sequence": "KLFPGDEIFSV",
        "source_gene_id": "G001",
        "source_gene_name": "GP63",
        "epitope_type": "MHC-I",
        "allele": "DLA-88*001:01",
        "ic50": 42.5,
        "rank": 0.3,
        "start_position": 120,
    }
    defaults.update(overrides)
    return Epitope(**defaults)


def _make_selected_epitope(**overrides) -> SelectedEpitope:
    """Create a SelectedEpitope for construct assembly tests."""
    defaults = {
        "peptide": "KLFPGDEIFSV",
        "gene_id": "G001",
        "gene_name": "GP63",
        "allele": "DLA-88*001:01",
        "ic50": 42.5,
        "rank": 0.3,
        "position": 120,
    }
    defaults.update(overrides)
    return SelectedEpitope(**defaults)


# ---------------------------------------------------------------------------
# Model tests
# ---------------------------------------------------------------------------


def test_epitope_creation() -> None:
    """Create an Epitope and verify all fields are set correctly."""
    ep = _make_epitope()

    assert ep.sequence == "KLFPGDEIFSV"
    assert ep.source_gene_id == "G001"
    assert ep.source_gene_name == "GP63"
    assert ep.epitope_type == "MHC-I"
    assert ep.allele == "DLA-88*001:01"
    assert ep.ic50 == 42.5
    assert ep.rank == 0.3
    assert ep.start_position == 120


def test_epitope_to_dict_from_dict() -> None:
    """Round-trip serialization of Epitope via to_dict / from_dict."""
    original = _make_epitope()
    data = original.to_dict()

    assert isinstance(data, dict)
    assert data["sequence"] == original.sequence
    assert data["ic50"] == original.ic50

    restored = Epitope.from_dict(data)
    assert restored.sequence == original.sequence
    assert restored.source_gene_id == original.source_gene_id
    assert restored.allele == original.allele
    assert restored.ic50 == original.ic50
    assert restored.rank == original.rank


def test_vaccine_construct_creation() -> None:
    """Create a VaccineConstruct and verify defaults are populated."""
    ep = _make_epitope()
    vc = VaccineConstruct(
        construct_id="VC-001",
        protein_sequence="MFVFLVLLPLVSS",
        mrna_sequence="AUGUUUGUUUUU",
        signal_peptide_name="tPA",
        adjuvant_name="RS09",
        epitopes=[ep],
    )

    assert vc.construct_id == "VC-001"
    assert vc.signal_peptide_name == "tPA"
    assert vc.adjuvant_name == "RS09"
    assert len(vc.epitopes) == 1
    assert vc.protein_sequence == "MFVFLVLLPLVSS"
    # Defaults for numeric fields should be zero.
    assert vc.molecular_weight == 0.0
    assert vc.isoelectric_point == 0.0
    assert vc.instability_index == 0.0
    assert vc.gravy == 0.0
    assert vc.gc_content == 0.0
    assert vc.vaxijen_score is None
    assert vc.allergenicity is None


def test_vaccine_construct_to_dict() -> None:
    """Serialization should include an epitope_count key."""
    ep1 = _make_epitope(sequence="AAAAAAAAA")
    ep2 = _make_epitope(sequence="BBBBBBBBB")
    vc = VaccineConstruct(
        construct_id="VC-002",
        protein_sequence="MFVFLVLLPLVSS",
        mrna_sequence="AUGUUUGUUUUU",
        signal_peptide_name="tPA",
        adjuvant_name="RS09",
        epitopes=[ep1, ep2],
    )
    data = vc.to_dict()

    assert isinstance(data, dict)
    assert data["epitope_count"] == 2
    assert data["construct_id"] == "VC-002"


# ---------------------------------------------------------------------------
# Codon table tests
# ---------------------------------------------------------------------------


def test_codon_table_covers_all_amino_acids() -> None:
    """The codon table must cover all 20 standard amino acids plus stop."""
    covered_aas = set(CODON_TABLE.keys())
    assert STANDARD_AMINO_ACIDS.issubset(covered_aas), (
        f"Missing amino acids: {STANDARD_AMINO_ACIDS - covered_aas}"
    )
    assert "*" in covered_aas, "Stop codon entry '*' missing from codon table"


def test_codon_frequencies_sum() -> None:
    """For each amino acid, codon frequencies should sum to approximately 1.0."""
    for aa, codons in CODON_TABLE.items():
        total = sum(freq for _codon, freq in codons)
        assert abs(total - 1.0) < 0.05, (
            f"Frequencies for '{aa}' sum to {total}, expected ~1.0"
        )


def test_get_optimal_codon() -> None:
    """get_optimal_codon should return a valid 3-letter DNA codon."""
    codon = get_optimal_codon("M")
    assert len(codon) == 3
    assert codon == "ATG"  # Met has only one codon.
    assert all(c in "ATGC" for c in codon)

    # Leucine has multiple codons; result should still be valid DNA.
    leu_codon = get_optimal_codon("L")
    assert len(leu_codon) == 3
    assert all(c in "ATGC" for c in leu_codon)


def test_reverse_translate() -> None:
    """reverse_translate of 'MK' should return a 6-nt DNA string starting with ATG."""
    dna = reverse_translate("MK")
    assert len(dna) == 6
    assert dna[:3] == "ATG"
    assert all(c in "ATGC" for c in dna)


# ---------------------------------------------------------------------------
# Construct assembly tests
# ---------------------------------------------------------------------------


def test_assemble_construct_basic() -> None:
    """Given multiple epitopes, the protein should contain signal peptide,
    adjuvant, EAAAK linker, and epitopes joined by AAY linkers."""
    ep1 = _make_selected_epitope(peptide="KLFPGDEIFSV")
    ep2 = _make_selected_epitope(peptide="YMLDIFHEV")

    protein = assemble_construct([ep1, ep2])

    # The protein should contain the EAAAK rigid linker between
    # adjuvant and epitope region.
    assert "EAAAK" in protein
    # Epitopes should be joined by AAY flexible linkers.
    assert "AAY" in protein
    # Both epitope sequences must be present.
    assert "KLFPGDEIFSV" in protein
    assert "YMLDIFHEV" in protein


def test_assemble_construct_single_epitope() -> None:
    """With a single epitope, no AAY linker between epitopes is needed."""
    ep = _make_selected_epitope(peptide="KLFPGDEIFSV")

    protein = assemble_construct([ep])

    assert "KLFPGDEIFSV" in protein
    # The AAY linker should not appear between epitopes.
    epitope_region = protein[protein.index("KLFPGDEIFSV"):]
    assert "AAY" not in epitope_region.replace("KLFPGDEIFSV", "")


def test_assemble_construct_order() -> None:
    """Signal peptide must come first, followed by adjuvant."""
    ep = _make_selected_epitope(peptide="KLFPGDEIFSV")

    protein = assemble_construct([ep])

    epitope_pos = protein.index("KLFPGDEIFSV")
    eaaak_pos = protein.index("EAAAK")

    # EAAAK (between adjuvant and epitopes) must precede the epitope.
    assert eaaak_pos < epitope_pos, "EAAAK linker should precede epitopes"
    # The protein should start with the signal peptide, so EAAAK is not at 0.
    assert eaaak_pos > 0, "Signal peptide + adjuvant should precede EAAAK"


# ---------------------------------------------------------------------------
# mRNA assembly tests
# ---------------------------------------------------------------------------


def test_assemble_mrna() -> None:
    """The mRNA should start with a 5'UTR-like sequence and end with poly-A."""
    cds = reverse_translate("MK")
    mrna = assemble_mrna(cds)

    # Poly-A tail at the end.
    assert mrna.endswith("A" * 30), "mRNA should end with a poly-A tail"


def test_assemble_mrna_contains_orf() -> None:
    """The CDS must be present in the assembled mRNA output."""
    cds = reverse_translate("MKFV")
    mrna = assemble_mrna(cds)

    assert cds in mrna, "CDS should be embedded in the mRNA sequence"


def test_assemble_mrna_has_stop_codon() -> None:
    """A stop codon (TGA) must be present between the ORF and the 3'UTR."""
    cds = reverse_translate("MK")
    mrna = assemble_mrna(cds)

    # Find the CDS in the mRNA and check that a stop codon follows it.
    cds_end = mrna.index(cds) + len(cds)
    downstream = mrna[cds_end:]
    assert "TGA" in downstream or "TAA" in downstream or "TAG" in downstream, (
        "A stop codon should appear after the CDS"
    )


# ---------------------------------------------------------------------------
# GC content tests
# ---------------------------------------------------------------------------


def test_gc_content() -> None:
    """'ATGC' has 2/4 GC bases = 0.5."""
    assert calculate_gc_content("ATGC") == pytest.approx(0.5)


def test_gc_content_all_gc() -> None:
    """'GCGC' is 100% GC."""
    assert calculate_gc_content("GCGC") == pytest.approx(1.0)


def test_gc_content_no_gc() -> None:
    """'ATAT' is 0% GC."""
    assert calculate_gc_content("ATAT") == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Physicochemical tests
# ---------------------------------------------------------------------------


def test_compute_physicochemical() -> None:
    """A known short protein should return a molecular weight in a sane range."""
    # Insulin B-chain fragment (short, well-characterized).
    protein = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    result = compute_physicochemical(protein)

    mw = result["molecular_weight"]
    # Rough expected range for a ~30-residue peptide.
    assert 2000 < mw < 5000, f"Molecular weight {mw} outside expected range"


def test_compute_physicochemical_keys() -> None:
    """Result dict should contain all expected physicochemical keys."""
    result = compute_physicochemical("MFVFLVLLPLVSSQCVNL")
    expected_keys = {
        "molecular_weight",
        "isoelectric_point",
        "instability_index",
        "gravy",
    }
    assert expected_keys.issubset(result.keys()), (
        f"Missing keys: {expected_keys - result.keys()}"
    )


# ---------------------------------------------------------------------------
# Safety check tests (mock HTTP)
# ---------------------------------------------------------------------------


def test_predict_antigenicity_failure() -> None:
    """When the external HTTP call fails, predict_antigenicity returns None."""
    import requests as _requests

    with patch.object(
        construct_mod.requests,
        "post",
        side_effect=_requests.exceptions.ConnectionError("Network error"),
    ):
        result = predict_antigenicity("MFVFLVLLPLVSS")

    assert result is None


def test_predict_allergenicity_failure() -> None:
    """When the external HTTP call fails, predict_allergenicity returns None."""
    import requests as _requests

    with patch.object(
        construct_mod.requests,
        "post",
        side_effect=_requests.exceptions.ConnectionError("Network error"),
    ):
        result = predict_allergenicity("MFVFLVLLPLVSS")

    assert result is None
