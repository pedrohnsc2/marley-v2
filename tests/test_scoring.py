"""Tests for conservation and immunogenicity scoring."""

from __future__ import annotations

import importlib

import pytest

from core.models import (
    CONSERVATION_WEIGHT,
    IMMUNOGENICITY_WEIGHT,
    STATUS_PENDING,
    Candidate,
)

# The module filename starts with a digit, so we use importlib.
immuno_mod = importlib.import_module("pipeline.04_immunogenicity")
generate_peptides = immuno_mod.generate_peptides
calculate_immunogenicity = immuno_mod.calculate_immunogenicity
IC50_THRESHOLD = immuno_mod.IC50_THRESHOLD


# ---------------------------------------------------------------------------
# Candidate model tests
# ---------------------------------------------------------------------------


def test_candidate_creation() -> None:
    """A newly created Candidate should have sensible defaults."""
    c = Candidate(gene_id="G001", gene_name="Test Gene", sequence="ACDEF")

    assert c.gene_id == "G001"
    assert c.gene_name == "Test Gene"
    assert c.sequence == "ACDEF"
    assert c.has_signal_peptide is False
    assert c.conservation_score == 0.0
    assert c.immunogenicity_score == 0.0
    assert c.final_score == 0.0
    assert c.filters_passed == []
    assert c.status == STATUS_PENDING


def test_compute_final_score() -> None:
    """compute_final_score should produce the correct weighted average."""
    c = Candidate(gene_id="G002", gene_name="Scored Gene", sequence="ACDEF")
    c.conservation_score = 0.8
    c.immunogenicity_score = 0.6

    c.compute_final_score()

    expected = CONSERVATION_WEIGHT * 0.8 + IMMUNOGENICITY_WEIGHT * 0.6
    assert c.final_score == pytest.approx(expected)


def test_compute_final_score_zeros() -> None:
    """Zero sub-scores should yield a zero final score."""
    c = Candidate(gene_id="G003", gene_name="Zero", sequence="X")
    c.compute_final_score()
    assert c.final_score == 0.0


def test_to_dict_from_dict() -> None:
    """Round-trip through to_dict / from_dict should preserve all fields."""
    original = Candidate(
        gene_id="G004",
        gene_name="Roundtrip",
        sequence="MKTLLV",
        has_signal_peptide=True,
        conservation_score=0.75,
        immunogenicity_score=0.9,
        final_score=0.84,
        filters_passed=["surface_filter", "conservation"],
        status="approved",
    )

    data = original.to_dict()
    restored = Candidate.from_dict(data)

    assert restored.gene_id == original.gene_id
    assert restored.gene_name == original.gene_name
    assert restored.sequence == original.sequence
    assert restored.has_signal_peptide == original.has_signal_peptide
    assert restored.conservation_score == original.conservation_score
    assert restored.immunogenicity_score == original.immunogenicity_score
    assert restored.final_score == original.final_score
    assert restored.filters_passed == original.filters_passed
    assert restored.status == original.status


def test_from_dict_with_defaults() -> None:
    """from_dict should apply defaults for missing optional fields."""
    data = {"gene_id": "G005", "gene_name": "Minimal", "sequence": "AC"}
    c = Candidate.from_dict(data)

    assert c.has_signal_peptide is False
    assert c.conservation_score == 0.0
    assert c.final_score == 0.0
    assert c.filters_passed == []
    assert c.status == STATUS_PENDING


# ---------------------------------------------------------------------------
# Peptide generation tests
# ---------------------------------------------------------------------------


def test_generate_peptides() -> None:
    """generate_peptides should produce all overlapping 9-mers."""
    sequence = "ACDEFGHIKLMNP"  # length 13 -> 5 peptides of length 9
    peptides = generate_peptides(sequence, length=9)

    assert len(peptides) == 5
    assert peptides[0] == "ACDEFGHIK"
    # length 13, peptide length 9 -> indices 0..4, last peptide starts at index 4
    assert peptides[-1] == sequence[4:13]  # "FGHIKLMNP"
    # Verify each peptide is exactly 9 residues
    for p in peptides:
        assert len(p) == 9


def test_generate_peptides_short_sequence() -> None:
    """A sequence shorter than the peptide length should return an empty list."""
    assert generate_peptides("ACDE", length=9) == []


def test_generate_peptides_exact_length() -> None:
    """A sequence equal to peptide length should return exactly one peptide."""
    seq = "ACDEFGHIK"
    peptides = generate_peptides(seq, length=9)
    assert peptides == [seq]


# ---------------------------------------------------------------------------
# Immunogenicity scoring tests
# ---------------------------------------------------------------------------


def test_calculate_immunogenicity() -> None:
    """Predictions with mixed IC50 values should yield the correct fraction."""
    predictions = [
        {"peptide": "ACDEFGHIK", "allele": "DLA-88*50101", "ic50": 100.0, "rank": 1.0},
        {"peptide": "CDEFGHIKL", "allele": "DLA-88*50101", "ic50": 800.0, "rank": 5.0},
        {"peptide": "DEFGHIKLM", "allele": "DLA-88*50101", "ic50": 200.0, "rank": 2.0},
        {"peptide": "EFGHIKLMN", "allele": "DLA-88*50101", "ic50": 900.0, "rank": 8.0},
    ]

    score = calculate_immunogenicity(predictions)

    # 2 out of 4 unique peptides have IC50 < 500 (the threshold)
    assert score == pytest.approx(0.5)


def test_calculate_immunogenicity_all_binders() -> None:
    """When all peptides bind well, score should be 1.0."""
    predictions = [
        {"peptide": "AAA", "allele": "X", "ic50": 50.0, "rank": 0.1},
        {"peptide": "BBB", "allele": "X", "ic50": 100.0, "rank": 0.2},
    ]

    assert calculate_immunogenicity(predictions) == pytest.approx(1.0)


def test_calculate_immunogenicity_no_binders() -> None:
    """When no peptides bind, score should be 0.0."""
    predictions = [
        {"peptide": "AAA", "allele": "X", "ic50": 600.0, "rank": 10.0},
        {"peptide": "BBB", "allele": "X", "ic50": 999.0, "rank": 20.0},
    ]

    assert calculate_immunogenicity(predictions) == pytest.approx(0.0)


def test_calculate_immunogenicity_empty() -> None:
    """An empty prediction list should return 0.0."""
    assert calculate_immunogenicity([]) == 0.0


def test_calculate_immunogenicity_best_ic50_across_alleles() -> None:
    """The best IC50 per peptide across alleles should be used."""
    predictions = [
        # Same peptide, two alleles -- one above threshold, one below.
        {"peptide": "ACDEFGHIK", "allele": "DLA-88*50101", "ic50": 600.0, "rank": 5.0},
        {"peptide": "ACDEFGHIK", "allele": "DLA-88*50801", "ic50": 200.0, "rank": 1.0},
        # Another peptide, both above threshold.
        {"peptide": "CDEFGHIKL", "allele": "DLA-88*50101", "ic50": 700.0, "rank": 6.0},
        {"peptide": "CDEFGHIKL", "allele": "DLA-88*50801", "ic50": 800.0, "rank": 7.0},
    ]

    score = calculate_immunogenicity(predictions)

    # Peptide ACDEFGHIK: best IC50 = 200 < 500 -> binder
    # Peptide CDEFGHIKL: best IC50 = 700 >= 500 -> non-binder
    # 1 out of 2 unique peptides = 0.5
    assert score == pytest.approx(0.5)
