"""Testes para o Modulo 07 -- Gap Window Optimization."""

from __future__ import annotations

import importlib
import pytest


def _load_enumerator():
    return importlib.import_module('aso_math.07_gap_window.gap_enumerator')


def _load_rnase_model():
    return importlib.import_module('aso_math.07_gap_window.rnase_h_model')


def _load_run():
    return importlib.import_module('aso_math.07_gap_window.run')


_CACHE = {}

def _run_main():
    if 'result' not in _CACHE:
        mod = _load_run()
        _CACHE['result'] = mod.main()
    return _CACHE['result']


class TestGapEnumerator:

    def test_all_configs_satisfy_constraints(self):
        mod = _load_enumerator()
        configs = mod.enumerate_gap_configs(aso_length=25)
        for cfg in configs:
            assert 2 <= cfg.lna_5p <= 5
            assert 2 <= cfg.lna_3p <= 5
            assert cfg.dna_gap >= 8
            assert cfg.total_lna <= 10
            assert cfg.total_length == 25

    def test_configs_not_empty(self):
        mod = _load_enumerator()
        configs = mod.enumerate_gap_configs(aso_length=25)
        assert len(configs) > 0

    def test_current_design_included(self):
        mod = _load_enumerator()
        configs = mod.enumerate_gap_configs(aso_length=25)
        found = any(
            c.lna_5p == 4 and c.dna_gap == 17 and c.lna_3p == 4
            for c in configs
        )
        assert found

    def test_notation_format(self):
        mod = _load_enumerator()
        cfg = mod.GapWindowConfig(lna_5p=3, dna_gap=19, lna_3p=3)
        assert cfg.notation() == "3-19-3"

    def test_total_lna_property(self):
        mod = _load_enumerator()
        cfg = mod.GapWindowConfig(lna_5p=4, dna_gap=17, lna_3p=4)
        assert cfg.total_lna == 8

    def test_no_configs_for_short_aso(self):
        mod = _load_enumerator()
        configs = mod.enumerate_gap_configs(aso_length=10)
        assert len(configs) == 0

    def test_sorted_by_gap_descending(self):
        mod = _load_enumerator()
        configs = mod.enumerate_gap_configs(aso_length=25)
        gaps = [c.dna_gap for c in configs]
        assert gaps == sorted(gaps, reverse=True)


class TestRNaseHModel:

    def test_gap_below_8_zero(self):
        mod = _load_rnase_model()
        assert mod.rnase_h_efficiency(7) == 0.0
        assert mod.rnase_h_efficiency(5) == 0.0

    def test_gap_8_9_efficiency(self):
        mod = _load_rnase_model()
        assert mod.rnase_h_efficiency(8) == 0.6
        assert mod.rnase_h_efficiency(9) == 0.6

    def test_gap_10_12_efficiency(self):
        mod = _load_rnase_model()
        assert mod.rnase_h_efficiency(10) == 0.9
        assert mod.rnase_h_efficiency(11) == 0.9
        assert mod.rnase_h_efficiency(12) == 0.9

    def test_gap_13_15_efficiency(self):
        mod = _load_rnase_model()
        assert mod.rnase_h_efficiency(13) == 1.0
        assert mod.rnase_h_efficiency(14) == 1.0
        assert mod.rnase_h_efficiency(15) == 1.0

    def test_gap_16_plus_efficiency(self):
        mod = _load_rnase_model()
        assert mod.rnase_h_efficiency(16) == 0.95
        assert mod.rnase_h_efficiency(20) == 0.95

    def test_nuclease_resistance_proportional(self):
        mod = _load_rnase_model()
        assert mod.nuclease_resistance(0) == 0.0
        assert mod.nuclease_resistance(5) == 0.5
        assert mod.nuclease_resistance(10) == 1.0
        assert mod.nuclease_resistance(12) == 1.0


class TestCompositeScoring:

    def test_score_in_valid_range(self):
        mod_e = _load_enumerator()
        mod_s = _load_rnase_model()
        aso = 'ACAGAAACTGATACTTATATAGCGT'
        sl = 'AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG'
        configs = mod_e.enumerate_gap_configs(aso_length=25)
        for cfg in configs:
            result = mod_s.score_config(cfg, aso, sl)
            assert 0.0 <= result.composite_score <= 1.0

    def test_optimal_rnase_h_gap_ranks_well(self):
        """Gap 13-15 (optimal RNase H) should rank in top half for 25 nt."""
        mod_e = _load_enumerator()
        mod_s = _load_rnase_model()
        aso = 'ACAGAAACTGATACTTATATAGCGT'
        sl = 'AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG'
        configs = mod_e.enumerate_gap_configs(aso_length=25)
        scored = []
        for cfg in configs:
            scored.append(mod_s.score_config(cfg, aso, sl))
        scored.sort(key=lambda s: -s.composite_score)
        top_half = scored[:len(scored) // 2 + 1]
        has_optimal = any(13 <= s.config.dna_gap <= 15 for s in top_half)
        assert has_optimal

    def test_scored_config_fields(self):
        mod_e = _load_enumerator()
        mod_s = _load_rnase_model()
        aso = 'ACAGAAACTGATACTTATATAGCGT'
        sl = 'AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG'
        cfg = mod_e.GapWindowConfig(lna_5p=4, dna_gap=17, lna_3p=4)
        result = mod_s.score_config(cfg, aso, sl)
        assert result.config == cfg
        assert isinstance(result.dg_binding, float)
        assert isinstance(result.tm_base, float)
        assert isinstance(result.tm_adjusted, float)
        assert isinstance(result.rnase_h_efficiency, float)
        assert isinstance(result.nuclease_resistance, float)
        assert isinstance(result.composite_score, float)

    def test_tm_adjusted_includes_lna_boost(self):
        mod_e = _load_enumerator()
        mod_s = _load_rnase_model()
        aso = 'ACAGAAACTGATACTTATATAGCGT'
        sl = 'AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG'
        cfg = mod_e.GapWindowConfig(lna_5p=4, dna_gap=17, lna_3p=4)
        result = mod_s.score_config(cfg, aso, sl)
        expected_boost = 8 * 3.0
        assert abs(result.tm_adjusted - result.tm_base - expected_boost) < 0.01


class TestEndToEnd:

    def test_main_completes_without_error(self):
        result = _run_main()
        assert result['status'] == 'success'
        assert result['module'] == '07_gap_window'
        assert 'data' in result
        assert 'ranking' in result['data']

    def test_main_has_current_design_comparison(self):
        result = _run_main()
        comparison = result['data']['current_design_comparison']
        assert comparison['found'] is True
        assert comparison['rank'] >= 1

    def test_ranking_is_sorted(self):
        result = _run_main()
        ranking = result['data']['ranking']
        scores = [entry['composite_score'] for entry in ranking]
        assert scores == sorted(scores, reverse=True)

    def test_key_metrics_present(self):
        result = _run_main()
        metrics = result['summary']['key_metrics']
        assert 'total_configs_evaluated' in metrics
        assert 'best_config' in metrics
        assert 'best_score' in metrics
        assert 'current_design_rank' in metrics
        assert 'current_design_score' in metrics
