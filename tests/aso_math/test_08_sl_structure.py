"""Testes para o modulo 08 -- Estrutura secundaria do SL RNA."""

from __future__ import annotations

import importlib

import pytest

from aso_math.config import ASO_TARGET_END, ASO_TARGET_START, SL_SEQUENCE


def _import(name: str):
    return importlib.import_module(f"aso_math.08_sl_structure.{name}")


_structure = _import("structure")
_accessibility = _import("accessibility")

StructureResult = _structure.StructureResult
can_pair = _structure.can_pair
pair_score = _structure.pair_score
predict_structure = _structure.predict_structure

AccessibilityReport = _accessibility.AccessibilityReport
analyze_accessibility = _accessibility.analyze_accessibility


class TestBasePairing:
    def test_watson_crick_au(self) -> None:
        assert can_pair("A", "U") is True
        assert can_pair("U", "A") is True
    def test_watson_crick_gc(self) -> None:
        assert can_pair("G", "C") is True
        assert can_pair("C", "G") is True
    def test_wobble_gu(self) -> None:
        assert can_pair("G", "U") is True
        assert can_pair("U", "G") is True
    def test_invalid_pairs(self) -> None:
        assert can_pair("A", "A") is False
        assert can_pair("A", "C") is False
        assert can_pair("C", "U") is False
    def test_pair_scores(self) -> None:
        assert pair_score("A", "U") == 1.0
        assert pair_score("G", "C") == 1.0
        assert pair_score("G", "U") == 0.5
        assert pair_score("A", "A") == 0.0


class TestNussinov:
    def test_valid_dot_bracket_characters(self) -> None:
        result = predict_structure("GGGGAAAACCCC")
        assert all(c in set(".()") for c in result.dot_bracket)
    def test_dot_bracket_length_matches_sequence(self) -> None:
        seq = "AACUAACGCUAUAUAAG"
        result = predict_structure(seq)
        assert len(result.dot_bracket) == len(seq)
    def test_balanced_brackets(self) -> None:
        result = predict_structure("GGGGAAAACCCC")
        count = 0
        for ch in result.dot_bracket:
            if ch == "(": count += 1
            elif ch == ")": count -= 1
            assert count >= 0
        assert count == 0
    def test_simple_hairpin(self) -> None:
        result = predict_structure("GGGAAAACCC")
        assert result.n_pairs >= 2
        assert result.score >= 2.0
    def test_no_pairs_short_sequence(self) -> None:
        result = predict_structure("ACGU")
        assert result.n_pairs == 0
        assert result.dot_bracket == "...."
    def test_empty_sequence(self) -> None:
        result = predict_structure("")
        assert result.n_pairs == 0
        assert result.dot_bracket == ""
        assert result.score == 0.0
        assert result.mfe_estimate == 0.0
    def test_dna_to_rna_conversion(self) -> None:
        result = predict_structure("GGGTTTTCCC")
        assert "T" not in result.sequence
        assert "U" in result.sequence
    def test_mfe_estimate_negative(self) -> None:
        result = predict_structure("GGGGAAAACCCC")
        assert result.mfe_estimate < 0.0
    def test_paired_fraction_range(self) -> None:
        result = predict_structure("GGGGAAAACCCC")
        assert 0.0 <= result.paired_fraction <= 1.0
    def test_structure_result_is_dataclass(self) -> None:
        result = predict_structure("GGGGAAAACCCC")
        assert isinstance(result, StructureResult)
    def test_pairs_consistent_with_dot_bracket(self) -> None:
        result = predict_structure("GGGGAAAACCCC")
        for i, j in result.pairs:
            assert result.dot_bracket[i] == "("
            assert result.dot_bracket[j] == ")"


class TestAccessibility:
    def test_all_unpaired_returns_one(self) -> None:
        report = analyze_accessibility("..........", 0, 10)
        assert report.mean_accessibility == 1.0
        assert report.classification == "accessible"
        assert report.n_accessible == 10
        assert report.n_blocked == 0
    def test_all_paired_returns_zero(self) -> None:
        report = analyze_accessibility("((((()))))", 0, 10)
        assert report.mean_accessibility == 0.0
        assert report.classification == "blocked"
        assert report.n_accessible == 0
        assert report.n_blocked == 10
    def test_mixed_region(self) -> None:
        report = analyze_accessibility("..((..))..", 0, 10)
        assert 0.0 < report.mean_accessibility < 1.0
        assert report.n_accessible == 6
        assert report.n_blocked == 4
    def test_classification_accessible(self) -> None:
        report = analyze_accessibility("........()", 0, 10)
        assert report.classification == "accessible"
    def test_classification_partially_blocked(self) -> None:
        report = analyze_accessibility("......(())", 0, 10)
        assert report.classification == "partially_blocked"
    def test_classification_blocked(self) -> None:
        report = analyze_accessibility("(((..)))..", 2, 8)
        assert report.classification == "blocked"
    def test_subregion_analysis(self) -> None:
        report = analyze_accessibility("....((((....)))).......", 0, 4)
        assert report.mean_accessibility == 1.0
        assert report.classification == "accessible"
    def test_invalid_region_raises(self) -> None:
        with pytest.raises(ValueError):
            analyze_accessibility(".....", 3, 2)
        with pytest.raises(ValueError):
            analyze_accessibility(".....", 0, 10)
    def test_report_is_dataclass(self) -> None:
        report = analyze_accessibility("..........", 0, 5)
        assert isinstance(report, AccessibilityReport)


class TestSLRNAPrediction:
    def test_sl_rna_produces_structure(self) -> None:
        result = predict_structure(SL_SEQUENCE)
        assert len(result.dot_bracket) == len(SL_SEQUENCE)
        assert result.n_pairs >= 0
        assert isinstance(result, StructureResult)
    def test_sl_rna_has_pairs(self) -> None:
        result = predict_structure(SL_SEQUENCE)
        assert result.n_pairs > 0
        assert result.score > 0.0
    def test_sl_rna_mfe_negative(self) -> None:
        result = predict_structure(SL_SEQUENCE)
        assert result.mfe_estimate < 0.0
    def test_sl_rna_blocked_has_fewer_pairs(self) -> None:
        free = predict_structure(SL_SEQUENCE)
        blocked = set(range(ASO_TARGET_START, ASO_TARGET_END))
        bound = predict_structure(SL_SEQUENCE, blocked_positions=blocked)
        assert bound.n_pairs <= free.n_pairs
    def test_sl_rna_accessibility(self) -> None:
        result = predict_structure(SL_SEQUENCE)
        report = analyze_accessibility(
            result.dot_bracket, ASO_TARGET_START, ASO_TARGET_END,
        )
        assert 0.0 <= report.mean_accessibility <= 1.0
        assert report.classification in ("accessible", "partially_blocked", "blocked")


_MAIN_CACHE: dict[str, dict] = {}


def _run_main() -> dict:
    if "result" not in _MAIN_CACHE:
        mod = importlib.import_module("aso_math.08_sl_structure.run")
        _MAIN_CACHE["result"] = mod.main()
    return _MAIN_CACHE["result"]


class TestMain:
    def test_main_completes(self) -> None:
        envelope = _run_main()
        assert envelope["status"] == "success"
        assert envelope["module"] == "08_sl_structure"
        assert "data" in envelope
        assert "summary" in envelope
    def test_main_has_required_data(self) -> None:
        envelope = _run_main()
        data = envelope["data"]
        assert "free_structure" in data
        assert "bound_structure" in data
        assert "accessibility" in data
        assert "structural_impact" in data
    def test_main_metrics(self) -> None:
        envelope = _run_main()
        metrics = envelope["summary"]["key_metrics"]
        assert "free_n_pairs" in metrics
        assert "target_accessibility" in metrics
        assert "target_classification" in metrics
        assert "pairs_disrupted" in metrics
