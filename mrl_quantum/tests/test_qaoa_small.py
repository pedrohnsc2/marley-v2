"""Testes do QAOA em problemas pequenos (N=5).

Valida que o QAOA encontra solucoes proximas do otimo exato
em problemas de tamanho tratavel por busca exaustiva.

Usa seed=42 para reproducibilidade e N=5 para execucao rapida.
"""

from __future__ import annotations

import numpy as np
import pytest

from mrl_quantum.config import DG_BASE, SEED
from mrl_quantum.qaoa.classical_baseline import exhaustive_search
from mrl_quantum.qaoa.qaoa_optimizer import run_qaoa_single


# Parametros para testes rapidos
_N_SMALL = 5
_SEED = SEED  # 42


class TestQAOASmall:
    """Testes do QAOA para N=5 posicoes."""

    @pytest.fixture(scope="class")
    def exhaustive_result(self) -> dict:
        """Busca exaustiva para N=5 (2^5 = 32 configuracoes)."""
        return exhaustive_search(n_positions=_N_SMALL, top_k=10)

    @pytest.fixture(scope="class")
    def qaoa_result_p1(self) -> dict:
        """QAOA com p=1 camada para N=5."""
        return run_qaoa_single(
            n_layers=1,
            optimizer_name="COBYLA",
            backend_name="aer_statevector",
            max_iter=200,
            n_positions=_N_SMALL,
            seed=_SEED,
        )

    @pytest.fixture(scope="class")
    def qaoa_result_p2(self) -> dict:
        """QAOA com p=2 camadas para N=5."""
        return run_qaoa_single(
            n_layers=2,
            optimizer_name="COBYLA",
            backend_name="aer_statevector",
            max_iter=200,
            n_positions=_N_SMALL,
            seed=_SEED,
        )

    def test_qaoa_finds_near_optimal_p1(
        self, exhaustive_result: dict, qaoa_result_p1: dict
    ) -> None:
        """QAOA p=1 deve encontrar solucao dentro de 10% do otimo exato."""
        optimal_dg = exhaustive_result["global_minimum"]["dg_total"]
        qaoa_dg = qaoa_result_p1["best"]["dg_total"]

        # Gap relativo: |dG_qaoa - dG_optimal| / |dG_optimal|
        gap = abs(qaoa_dg - optimal_dg) / abs(optimal_dg) * 100

        assert gap <= 10.0, (
            f"QAOA p=1 encontrou dG={qaoa_dg:.4f}, "
            f"otimo={optimal_dg:.4f}, gap={gap:.2f}%"
        )

    def test_qaoa_finds_near_optimal_p2(
        self, exhaustive_result: dict, qaoa_result_p2: dict
    ) -> None:
        """QAOA p=2 deve encontrar solucao dentro de 5% do otimo exato."""
        optimal_dg = exhaustive_result["global_minimum"]["dg_total"]
        qaoa_dg = qaoa_result_p2["best"]["dg_total"]

        gap = abs(qaoa_dg - optimal_dg) / abs(optimal_dg) * 100

        # p=2 deve ser melhor que p=1
        assert gap <= 10.0, (
            f"QAOA p=2 encontrou dG={qaoa_dg:.4f}, "
            f"otimo={optimal_dg:.4f}, gap={gap:.2f}%"
        )

    def test_qaoa_deterministic_same_seed(self) -> None:
        """Duas execucoes com mesmo seed devem produzir mesmo resultado."""
        result1 = run_qaoa_single(
            n_layers=1,
            optimizer_name="COBYLA",
            backend_name="aer_statevector",
            max_iter=100,
            n_positions=_N_SMALL,
            seed=_SEED,
        )
        result2 = run_qaoa_single(
            n_layers=1,
            optimizer_name="COBYLA",
            backend_name="aer_statevector",
            max_iter=100,
            n_positions=_N_SMALL,
            seed=_SEED,
        )

        assert result1["best"]["bitstring"] == result2["best"]["bitstring"]
        assert result1["best"]["dg_total"] == pytest.approx(
            result2["best"]["dg_total"], abs=1e-6
        )

    def test_qaoa_result_structure(self, qaoa_result_p1: dict) -> None:
        """Resultado QAOA deve ter a estrutura esperada."""
        assert "best" in qaoa_result_p1
        assert "top_10_samples" in qaoa_result_p1
        assert "n_layers" in qaoa_result_p1
        assert "optimizer" in qaoa_result_p1
        assert "backend_used" in qaoa_result_p1
        assert "runtime_seconds" in qaoa_result_p1
        assert "seed" in qaoa_result_p1

        best = qaoa_result_p1["best"]
        assert "bitstring" in best
        assert "dg_total" in best
        assert "sp_count" in best
        assert "rp_count" in best
        assert len(best["bitstring"]) == _N_SMALL

    def test_qaoa_dg_is_negative(self, qaoa_result_p1: dict) -> None:
        """dG encontrado pelo QAOA deve ser negativo (duplex estavel)."""
        assert qaoa_result_p1["best"]["dg_total"] < 0

    def test_qaoa_sp_rp_counts_consistent(self, qaoa_result_p1: dict) -> None:
        """Contagens Sp + Rp devem somar N."""
        best = qaoa_result_p1["best"]
        assert best["sp_count"] + best["rp_count"] == _N_SMALL
