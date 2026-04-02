"""Tests for the conservation scoring logic."""

from __future__ import annotations

import importlib

import pytest

# The module filename starts with a digit, so we use importlib.
conservation_mod = importlib.import_module("pipeline.03_conservation")
calculate_conservation = conservation_mod.calculate_conservation
COMPARISON_STRAINS = conservation_mod.COMPARISON_STRAINS


def test_calculate_conservation_with_matches() -> None:
    """When BLAST hits match all comparison strains, score should reflect their identity."""
    blast_hits = [
        {
            "subject_id": "ACC001",
            "identity_percent": 95.0,
            "evalue": 1e-50,
            "organism": "Leishmania infantum LEM3323",
        },
        {
            "subject_id": "ACC002",
            "identity_percent": 90.0,
            "evalue": 1e-40,
            "organism": "Leishmania chagasi",
        },
    ]

    score = calculate_conservation(blast_hits)

    # Average identity = (95.0 + 90.0) / 2 = 92.5, normalised = 0.925
    assert score == pytest.approx(0.925, abs=0.001)


def test_calculate_conservation_no_matches() -> None:
    """Empty blast_hits should return 0.0."""
    assert calculate_conservation([]) == 0.0


def test_calculate_conservation_no_matching_strains() -> None:
    """Hits that don't match any comparison strain should return 0.0."""
    blast_hits = [
        {
            "subject_id": "ACC003",
            "identity_percent": 99.0,
            "evalue": 1e-80,
            "organism": "Trypanosoma cruzi",
        },
    ]

    assert calculate_conservation(blast_hits) == 0.0


def test_calculate_conservation_partial() -> None:
    """When only some comparison strains match, score averages only matched strains."""
    blast_hits = [
        {
            "subject_id": "ACC004",
            "identity_percent": 80.0,
            "evalue": 1e-30,
            "organism": "Leishmania infantum LEM3323",
        },
        # No hit for Leishmania chagasi
    ]

    score = calculate_conservation(blast_hits)

    # Only one strain matched at 80% -> 80.0 / 100.0 = 0.80
    assert score == pytest.approx(0.80, abs=0.001)


def test_calculate_conservation_best_per_strain() -> None:
    """When multiple hits match the same strain, the best identity should be used."""
    blast_hits = [
        {
            "subject_id": "ACC005",
            "identity_percent": 70.0,
            "evalue": 1e-20,
            "organism": "Leishmania chagasi",
        },
        {
            "subject_id": "ACC006",
            "identity_percent": 92.0,
            "evalue": 1e-50,
            "organism": "Leishmania chagasi",
        },
    ]

    score = calculate_conservation(blast_hits)

    # Best for L. chagasi = 92.0 -> 0.92
    assert score == pytest.approx(0.92, abs=0.001)
