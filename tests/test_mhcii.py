"""Tests for MHC-II prediction in module 04 -- immunogenicity scoring."""

from __future__ import annotations

import importlib
from unittest.mock import MagicMock, patch

import pytest
import requests

# The module filename starts with a digit, so we use importlib.
immuno_mod = importlib.import_module("pipeline.04_immunogenicity")

predict_mhcii_binding = immuno_mod.predict_mhcii_binding
MHCII_ALLELES = immuno_mod.MHCII_ALLELES
MHCII_IC50_THRESHOLD = immuno_mod.MHCII_IC50_THRESHOLD


# ---------------------------------------------------------------------------
# MHC-II binding prediction
# ---------------------------------------------------------------------------


def test_predict_mhcii_binding() -> None:
    """Mock POST, return mock TSV with ic50 columns, verify parsing."""
    mock_tsv = (
        "allele\tseq_num\tstart\tend\tlength\tpeptide\tic50\trank\n"
        "HLA-DRB1*01:01\t1\t1\t15\t15\tACDEFGHIKLMNPQR\t125.5\t0.8\n"
        "HLA-DRB1*01:01\t1\t2\t16\t15\tCDEFGHIKLMNPQRS\t800.0\t5.2\n"
    )

    mock_response = MagicMock()
    mock_response.text = mock_tsv
    mock_response.raise_for_status = MagicMock()

    with patch.object(
        immuno_mod.requests, "post", return_value=mock_response
    ) as mock_post:
        result = predict_mhcii_binding("ACDEFGHIKLMNPQRS", "HLA-DRB1*01:01")

    mock_post.assert_called_once()
    assert len(result) == 2

    # Verify parsed fields.
    first = result[0]
    assert first["peptide"] == "ACDEFGHIKLMNPQR"
    assert first["ic50"] == pytest.approx(125.5)
    assert first["rank"] == pytest.approx(0.8)
    assert "allele" in first


def test_predict_mhcii_binding_failure() -> None:
    """Mock connection error, verify returns empty list."""
    with patch.object(
        immuno_mod.requests,
        "post",
        side_effect=requests.exceptions.ConnectionError("Network error"),
    ):
        result = predict_mhcii_binding("ACDEFGHIKLMNPQRS", "HLA-DRB1*01:01")

    assert result == []


# ---------------------------------------------------------------------------
# MHC-II constants
# ---------------------------------------------------------------------------


def test_mhcii_alleles_are_hla() -> None:
    """Verify MHCII_ALLELES contains HLA-DRB1 (not DLA).

    Canine DLA class II is not supported by IEDB, so the pipeline uses
    HLA-DRB1 alleles as a cross-species proxy.
    """
    assert len(MHCII_ALLELES) >= 1
    for allele in MHCII_ALLELES:
        assert allele.startswith("HLA-DRB1"), (
            f"Expected HLA-DRB1 allele, got {allele}"
        )
        assert not allele.startswith("DLA"), (
            f"MHC-II alleles should be HLA (proxy), not DLA: {allele}"
        )


def test_mhcii_threshold() -> None:
    """Verify MHCII_IC50_THRESHOLD == 1000.0 (relaxed vs MHC-I)."""
    assert MHCII_IC50_THRESHOLD == 1000.0
