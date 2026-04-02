"""Tests for the report generation module."""

from __future__ import annotations

import csv
import importlib
from pathlib import Path
from unittest.mock import patch

import pytest

from core.models import STATUS_APPROVED, STATUS_PENDING, Candidate

# The module filename starts with a digit, so we use importlib.
report_mod = importlib.import_module("pipeline.05_report")
_build_summary_stats = report_mod._build_summary_stats
_print_top_candidates = report_mod._print_top_candidates
_generate_csv = report_mod._generate_csv
_generate_markdown = report_mod._generate_markdown


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_candidate(**overrides) -> Candidate:
    """Create a Candidate with sensible defaults, overridden by kwargs."""
    defaults = {
        "gene_id": "G001",
        "gene_name": "Test Gene",
        "sequence": "ACDEF",
        "has_signal_peptide": True,
        "conservation_score": 0.85,
        "immunogenicity_score": 0.70,
        "final_score": 0.76,
        "status": STATUS_APPROVED,
    }
    defaults.update(overrides)
    return Candidate(**defaults)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_build_summary_stats() -> None:
    """_build_summary_stats should return correct counts for each category."""
    all_candidates = [
        _make_candidate(gene_id="G001", has_signal_peptide=True, conservation_score=0.9, status=STATUS_APPROVED),
        _make_candidate(gene_id="G002", has_signal_peptide=False, conservation_score=0.5, status=STATUS_APPROVED),
        _make_candidate(gene_id="G003", has_signal_peptide=True, conservation_score=0.0, status=STATUS_PENDING),
    ]
    scored = [all_candidates[0], all_candidates[1]]

    stats = _build_summary_stats(all_candidates, scored)

    assert stats["total"] == 3
    assert stats["surface_filter"] == 2  # G001, G003 have signal peptide
    assert stats["conservation_filter"] == 2  # G001, G002 have score > 0 and approved
    assert stats["final_candidates"] == 2


def test_print_top_candidates_empty(capsys) -> None:
    """_print_top_candidates with an empty list should not crash."""
    _print_top_candidates([])

    captured = capsys.readouterr()
    assert "No scored candidates to display" in captured.out


def test_generate_csv(tmp_path: Path) -> None:
    """_generate_csv should write a valid CSV with candidate data."""
    candidates = [
        _make_candidate(gene_id="G001", gene_name="Alpha", final_score=0.90),
        _make_candidate(gene_id="G002", gene_name="Beta", final_score=0.85),
    ]

    output_csv = str(tmp_path / "top_candidates.csv")

    with patch.object(report_mod, "OUTPUT_CSV", output_csv):
        result = _generate_csv(candidates)

    assert result == output_csv

    # Read back and verify contents.
    with open(output_csv, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    assert len(rows) == 2
    assert rows[0]["gene_id"] == "G001"
    assert rows[1]["gene_id"] == "G002"
    assert "final_score" in rows[0]


def test_generate_markdown(tmp_path: Path) -> None:
    """_generate_markdown should produce a report with expected sections."""
    candidates = [
        _make_candidate(gene_id="G001", gene_name="Alpha", final_score=0.90),
    ]
    stats = {
        "total": 100,
        "surface_filter": 50,
        "conservation_filter": 20,
        "final_candidates": 1,
    }

    output_md = str(tmp_path / "report.md")

    with patch.object(report_mod, "OUTPUT_REPORT", output_md):
        result = _generate_markdown(candidates, stats)

    assert result == output_md

    content = Path(output_md).read_text()

    # Verify expected sections are present.
    assert "# Marley" in content
    assert "## Summary Statistics" in content
    assert "## Top" in content
    assert "## Methodology" in content
    assert "## Next Steps" in content
    # Verify candidate data appears in the table.
    assert "G001" in content
    assert "Alpha" in content
