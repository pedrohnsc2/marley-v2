"""Tests for surface protein filtering logic."""

from __future__ import annotations

import importlib
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# The module filename starts with a digit, so we use importlib.
filter_mod = importlib.import_module("pipeline.02_filter_surface")
parse_fasta = filter_mod.parse_fasta
filter_surface_proteins = filter_mod.filter_surface_proteins
_write_temp_fasta = filter_mod._write_temp_fasta


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def fasta_file(tmp_path: Path) -> Path:
    """Create a temporary FASTA file with two annotated sequences."""
    content = textwrap.dedent("""\
        >LINF_010010 hypothetical protein
        MKTLLLTLVVVAALQ
        >LINF_010020 surface antigen
        AACDEFGHIKLMNPR
    """)
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(content)
    return fasta


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_parse_fasta(fasta_file: Path) -> None:
    """parse_fasta should return (gene_id, gene_name, sequence) tuples."""
    results = parse_fasta(str(fasta_file))

    assert len(results) == 2

    gene_id, gene_name, sequence = results[0]
    assert gene_id == "LINF_010010"
    assert gene_name == "hypothetical protein"
    assert sequence == "MKTLLLTLVVVAALQ"

    gene_id_2, gene_name_2, sequence_2 = results[1]
    assert gene_id_2 == "LINF_010020"
    assert gene_name_2 == "surface antigen"
    assert sequence_2 == "AACDEFGHIKLMNPR"


def test_parse_fasta_missing_file() -> None:
    """parse_fasta should raise FileNotFoundError for a missing path."""
    with pytest.raises(FileNotFoundError):
        parse_fasta("/nonexistent/missing.fasta")


@patch("pathlib.Path.exists", return_value=True)
def test_filter_skips_existing(mock_exists: MagicMock) -> None:
    """When the output file already exists, filter_surface_proteins should short-circuit."""
    result = filter_surface_proteins(force=False)
    assert result == filter_mod.OUTPUT_FILE


def test_write_temp_fasta() -> None:
    """_write_temp_fasta should produce a valid FASTA file and return its path."""
    sequences = [
        ("gene_A", "ACDEFG"),
        ("gene_B", "HIKLMN"),
    ]

    path = _write_temp_fasta(sequences)
    content = Path(path).read_text()

    assert ">gene_A" in content
    assert "ACDEFG" in content
    assert ">gene_B" in content
    assert "HIKLMN" in content

    # Cleanup
    Path(path).unlink(missing_ok=True)


def test_write_temp_fasta_empty() -> None:
    """An empty input list should produce an empty file."""
    path = _write_temp_fasta([])
    content = Path(path).read_text()

    assert content.strip() == ""

    Path(path).unlink(missing_ok=True)
