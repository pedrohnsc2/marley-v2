"""Tests for genome fetching and raw data retrieval."""

from __future__ import annotations

import importlib
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

# The module filename starts with a digit, so we use importlib.
fetch_mod = importlib.import_module("pipeline.01_fetch_genome")
count_sequences = fetch_mod.count_sequences
fetch_genome = fetch_mod.fetch_genome


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    """Create a temporary FASTA file with three sequences."""
    content = textwrap.dedent("""\
        >seq1 First protein
        MKTLLLTLVVV
        >seq2 Second protein
        AACDEFGHIKLM
        >seq3 Third protein
        NPRSTVWY
    """)
    fasta = tmp_path / "test.fasta"
    fasta.write_text(content)
    return fasta


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_count_sequences(fasta_file: Path) -> None:
    """count_sequences should return the number of '>' header lines."""
    result = count_sequences(str(fasta_file))
    assert result == 3


def test_count_sequences_empty(tmp_path: Path) -> None:
    """An empty FASTA file should yield zero sequences."""
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    assert count_sequences(str(empty)) == 0


def test_count_sequences_missing_file() -> None:
    """count_sequences should raise FileNotFoundError for a missing path."""
    with pytest.raises(FileNotFoundError):
        count_sequences("/nonexistent/path/missing.fasta")


@patch.object(fetch_mod, "count_sequences", return_value=42)
@patch("pathlib.Path.exists", return_value=True)
def test_fetch_genome_existing_file(
    mock_exists: MagicMock,
    mock_count: MagicMock,
) -> None:
    """When the output file already exists, fetch_genome should skip download."""
    result = fetch_genome(force=False)
    assert result == fetch_mod.OUTPUT_FILE
    # count_sequences is called to log how many sequences are present
    mock_count.assert_called_once()


@patch.object(fetch_mod, "time")
@patch.object(fetch_mod, "_try_direct_download", return_value=False)
@patch.object(fetch_mod, "_try_api_download")
@patch("pathlib.Path.exists", return_value=False)
@patch("pathlib.Path.mkdir")
def test_fetch_genome_retry(
    mock_mkdir: MagicMock,
    mock_exists: MagicMock,
    mock_api: MagicMock,
    mock_direct: MagicMock,
    mock_time: MagicMock,
) -> None:
    """fetch_genome should retry on failure and succeed on the third attempt."""
    # Fail twice, succeed on third attempt
    mock_api.side_effect = [None, None, ">seq1\nACDEF"]

    with patch("pathlib.Path.write_text"):
        with patch.object(fetch_mod, "count_sequences", return_value=1):
            result = fetch_genome(force=True)

    assert result == fetch_mod.OUTPUT_FILE
    assert mock_api.call_count == 3
