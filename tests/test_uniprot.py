"""Tests for the UniProt REST API client (core/uniprot.py)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from core import uniprot
from core.uniprot import (
    SKIP_ANTIGENS,
    VALIDATED_ACCESSIONS,
    fetch_all_validated_sequences,
    fetch_fasta_by_accession,
    search_uniprot,
)


# ---------------------------------------------------------------------------
# fetch_fasta_by_accession
# ---------------------------------------------------------------------------


def test_fetch_fasta_by_accession() -> None:
    """Mock requests.get, return mock FASTA, verify sequence extraction."""
    mock_fasta = ">sp|A4HZU7|A2_LEIIN Amastigote surface protein 2\nMKTLLLTLVVV\nACDEFGHIKLM\n"

    mock_response = MagicMock()
    mock_response.text = mock_fasta
    mock_response.raise_for_status = MagicMock()

    with patch.object(uniprot.requests, "get", return_value=mock_response) as mock_get:
        result = fetch_fasta_by_accession("A4HZU7")

    mock_get.assert_called_once()
    assert result == "MKTLLLTLVVVACDEFGHIKLM"


def test_fetch_fasta_invalid_accession() -> None:
    """Mock 400 response, verify returns None."""
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
        "400 Bad Request"
    )

    with patch.object(uniprot.requests, "get", return_value=mock_response):
        result = fetch_fasta_by_accession("INVALID_ACC")

    assert result is None


# ---------------------------------------------------------------------------
# search_uniprot
# ---------------------------------------------------------------------------


def test_search_uniprot() -> None:
    """Mock search response JSON, verify returns sequence."""
    mock_json = {
        "results": [
            {
                "primaryAccession": "Q9NJC8",
                "sequence": {"value": "MKVLFAALLGFLTG"},
            }
        ]
    }

    mock_response = MagicMock()
    mock_response.json.return_value = mock_json
    mock_response.raise_for_status = MagicMock()

    with patch.object(uniprot.requests, "get", return_value=mock_response):
        result = search_uniprot("KMP-11")

    assert result == "MKVLFAALLGFLTG"


def test_search_uniprot_no_results() -> None:
    """Mock empty response, verify returns None."""
    mock_json = {"results": []}

    mock_response = MagicMock()
    mock_response.json.return_value = mock_json
    mock_response.raise_for_status = MagicMock()

    with patch.object(uniprot.requests, "get", return_value=mock_response):
        result = search_uniprot("nonexistent_protein_xyz")

    assert result is None


# ---------------------------------------------------------------------------
# fetch_all_validated_sequences
# ---------------------------------------------------------------------------


def test_fetch_all_validated_sequences() -> None:
    """Mock all API calls, verify dict has expected keys."""
    # All non-skipped antigens with known accessions use fetch_fasta_by_accession;
    # those with None accession use search_uniprot.
    expected_fetchable = {
        gene_id
        for gene_id in VALIDATED_ACCESSIONS
        if gene_id not in SKIP_ANTIGENS
    }

    with (
        patch.object(
            uniprot,
            "fetch_fasta_by_accession",
            return_value="MKFAKESEQUENCE",
        ) as mock_fasta,
        patch.object(
            uniprot,
            "search_uniprot",
            return_value="MKFAKESEARCH",
        ) as mock_search,
        patch.object(uniprot.time, "sleep"),
    ):
        result = fetch_all_validated_sequences()

    assert isinstance(result, dict)
    assert set(result.keys()) == expected_fetchable

    # Antigens with known accessions should use fetch_fasta_by_accession.
    accession_calls = {
        gene_id
        for gene_id, acc in VALIDATED_ACCESSIONS.items()
        if acc is not None and gene_id not in SKIP_ANTIGENS
    }
    assert mock_fasta.call_count == len(accession_calls)

    # Antigens needing search should use search_uniprot.
    search_calls = {
        gene_id
        for gene_id, acc in VALIDATED_ACCESSIONS.items()
        if acc is None and gene_id not in SKIP_ANTIGENS
    }
    assert mock_search.call_count == len(search_calls)


def test_fetch_skips_composite_antigens() -> None:
    """Verify LBSap and Lutzomyia are skipped (not fetched)."""
    assert "LBSap_antigens" in SKIP_ANTIGENS
    assert "Lutzomyia_longipalpis_proteins" in SKIP_ANTIGENS

    with (
        patch.object(
            uniprot,
            "fetch_fasta_by_accession",
            return_value="MKFAKESEQUENCE",
        ) as mock_fasta,
        patch.object(
            uniprot,
            "search_uniprot",
            return_value="MKFAKESEARCH",
        ) as mock_search,
        patch.object(uniprot.time, "sleep"),
    ):
        result = fetch_all_validated_sequences()

    # Composite antigens must not appear in results.
    assert "LBSap_antigens" not in result
    assert "Lutzomyia_longipalpis_proteins" not in result

    # Verify they were never queried.
    all_fasta_args = [call.args[0] for call in mock_fasta.call_args_list]
    all_search_args = [call.args[0] for call in mock_search.call_args_list]
    for skip_id in SKIP_ANTIGENS:
        assert skip_id not in all_fasta_args
        assert skip_id not in all_search_args
