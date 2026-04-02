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
            "organism": "Leishmania donovani donovani",
        },
        {
            "subject_id": "ACC002",
            "identity_percent": 90.0,
            "evalue": 1e-40,
            "organism": "Leishmania major strain Friedlin",
        },
        {
            "subject_id": "ACC003",
            "identity_percent": 88.0,
            "evalue": 1e-35,
            "organism": "Leishmania braziliensis M2904",
        },
    ]

    score = calculate_conservation(blast_hits)

    # Average identity = (95.0 + 90.0 + 88.0) / 3 = 91.0, normalised = 0.91
    assert score == pytest.approx(0.91, abs=0.001)


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
            "organism": "Leishmania donovani",
        },
        # No hit for L. major or L. braziliensis
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
            "organism": "Leishmania major",
        },
        {
            "subject_id": "ACC006",
            "identity_percent": 92.0,
            "evalue": 1e-50,
            "organism": "Leishmania major strain Friedlin",
        },
    ]

    score = calculate_conservation(blast_hits)

    # Best for L. major = 92.0 -> 0.92
    assert score == pytest.approx(0.92, abs=0.001)


def test_calculate_conservation_ignores_self_reference() -> None:
    """Hits from L. infantum (reference strain) should be excluded from scoring."""
    blast_hits = [
        {
            "subject_id": "ACC007",
            "identity_percent": 100.0,
            "evalue": 0.0,
            "organism": "Leishmania infantum JPCM5",
        },
        {
            "subject_id": "ACC008",
            "identity_percent": 85.0,
            "evalue": 1e-30,
            "organism": "Leishmania donovani",
        },
    ]

    score = calculate_conservation(blast_hits)

    # L. infantum hit should be ignored, only L. donovani counts
    assert score == pytest.approx(0.85, abs=0.001)
