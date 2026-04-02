"""Tests for the validated antigens data and loading logic."""

from __future__ import annotations

import importlib
from unittest.mock import MagicMock, patch

import pytest

# The module filename starts with a digit, so we use importlib.
immuno_mod = importlib.import_module("pipeline.04_immunogenicity")
VALIDATED_ANTIGENS = immuno_mod.VALIDATED_ANTIGENS
load_validated_antigens = immuno_mod.load_validated_antigens

REQUIRED_FIELDS = {"gene_id", "gene_name", "source", "evidence", "final_score"}


def test_validated_antigens_count() -> None:
    """VALIDATED_ANTIGENS should contain exactly 8 entries."""
    assert len(VALIDATED_ANTIGENS) == 8


def test_validated_antigens_scores() -> None:
    """All validated antigen scores should be between 0.0 and 1.0."""
    for antigen in VALIDATED_ANTIGENS:
        score = antigen["final_score"]
        assert 0.0 <= score <= 1.0, (
            f"{antigen['gene_id']} has score {score} outside [0.0, 1.0]"
        )


def test_validated_antigens_required_fields() -> None:
    """Each validated antigen dict should contain all required fields."""
    for antigen in VALIDATED_ANTIGENS:
        missing = REQUIRED_FIELDS - set(antigen.keys())
        assert not missing, (
            f"{antigen.get('gene_id', '???')} is missing fields: {missing}"
        )


@patch.object(immuno_mod, "upsert_candidate")
def test_load_validated_antigens_creates_candidates(mock_upsert: MagicMock) -> None:
    """load_validated_antigens should call upsert_candidate once per antigen."""
    count = load_validated_antigens()

    assert count == 8
    assert mock_upsert.call_count == 8

    # Verify the first call received a Candidate with the right gene_id.
    first_call_candidate = mock_upsert.call_args_list[0][0][0]
    assert first_call_candidate.gene_id == VALIDATED_ANTIGENS[0]["gene_id"]
    assert first_call_candidate.priority is True
