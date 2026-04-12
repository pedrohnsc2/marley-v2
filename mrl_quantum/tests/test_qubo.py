"""Testes da formulacao QUBO para estereoisomeros PS.

Valida a construcao do QuadraticProgram, a avaliacao de bitstrings,
e a consistencia entre a formulacao QUBO e o calculo manual de dG.
"""

from __future__ import annotations

import numpy as np
import pytest

from mrl_quantum.config import (
    DG_BASE,
    J_COUPLING,
    N_POSITIONS,
    RP_CONSECUTIVE_PENALTY,
    RP_CORRECTION,
    SP_CORRECTION,
)
from mrl_quantum.qaoa.qubo_formulation import (
    bitstring_to_config,
    build_qubo,
    evaluate_bitstring,
    get_aso_position_label,
)


class TestEvaluateBitstring:
    """Testes para avaliacao de bitstrings individuais."""

    def test_all_sp_gives_maximum_stability(self) -> None:
        """Configuracao all-Sp deve dar o menor dG (mais estavel).

        Todas as posicoes Sp: cada uma contribui SP_CORRECTION = -0.2 kcal/mol.
        Acoplamento J: (N-1) * J_COUPLING * 1 * 1 (todos Sp adjacentes).
        Sem penalidade Rp (nenhum Rp consecutivo).
        """
        n = 10
        all_sp = [1] * n

        result = evaluate_bitstring(all_sp, n_positions=n)

        # Calculo manual
        correction = n * SP_CORRECTION  # 10 * (-0.2) = -2.0
        coupling = (n - 1) * J_COUPLING * 1 * 1  # 9 * (-0.05) = -0.45
        penalty = 0.0  # Nenhum par Rp-Rp
        expected_dg = DG_BASE + correction + coupling + penalty

        assert result["dg_total"] == pytest.approx(expected_dg, abs=1e-3)
        assert result["sp_count"] == n
        assert result["rp_count"] == 0
        assert result["rp_triples"] == 0

    def test_all_rp_gives_minimum_stability(self) -> None:
        """Configuracao all-Rp deve dar o maior dG (menos estavel).

        Todas as posicoes Rp: cada uma contribui RP_CORRECTION = +0.1 kcal/mol.
        Sem acoplamento J (J * 0 * 0 = 0).
        Penalidade maxima: (N-1) pares Rp-Rp.
        """
        n = 10
        all_rp = [0] * n

        result = evaluate_bitstring(all_rp, n_positions=n)

        # Calculo manual
        correction = n * RP_CORRECTION  # 10 * 0.1 = 1.0
        coupling = 0.0  # J * 0 * 0 = 0
        pen_per_pair = RP_CONSECUTIVE_PENALTY / 2.0
        penalty = (n - 1) * pen_per_pair  # 9 * 0.25 = 2.25
        expected_dg = DG_BASE + correction + coupling + penalty

        assert result["dg_total"] == pytest.approx(expected_dg, abs=1e-3)
        assert result["sp_count"] == 0
        assert result["rp_count"] == n
        assert result["rp_triples"] == n - 2  # 8 tripletas

    def test_all_sp_more_stable_than_all_rp(self) -> None:
        """All-Sp deve ter dG mais negativo que all-Rp."""
        n = 10
        sp_result = evaluate_bitstring([1] * n, n_positions=n)
        rp_result = evaluate_bitstring([0] * n, n_positions=n)

        assert sp_result["dg_total"] < rp_result["dg_total"]

    def test_three_consecutive_rp_penalty(self) -> None:
        """Tres Rp consecutivos devem gerar penalidade maior que alternados."""
        n = 6
        # 3 Rp consecutivos: [0, 0, 0, 1, 1, 1]
        consecutive = [0, 0, 0, 1, 1, 1]
        # Rp alternados com Sp: [0, 1, 0, 1, 0, 1]
        alternating = [0, 1, 0, 1, 0, 1]

        result_consec = evaluate_bitstring(consecutive, n_positions=n)
        result_altern = evaluate_bitstring(alternating, n_positions=n)

        # Consecutivos Rp devem ter penalidade
        assert result_consec["penalty_term"] > result_altern["penalty_term"]
        # E portanto dG mais positivo (menos estavel)
        assert result_consec["rp_triples"] == 1
        assert result_altern["rp_triples"] == 0

    def test_j_coupling_effect(self) -> None:
        """Acoplamento J entre posicoes adjacentes Sp deve estabilizar.

        J_COUPLING = -0.05 (negativo), entao Sp-Sp adjacentes contribuem
        negativamente ao dG (estabilizam).
        """
        n = 4
        # Dois Sp adjacentes: [0, 1, 1, 0]
        adjacent_sp = [0, 1, 1, 0]
        # Dois Sp separados: [1, 0, 0, 1]
        separated_sp = [1, 0, 0, 1]

        result_adj = evaluate_bitstring(adjacent_sp, n_positions=n)
        result_sep = evaluate_bitstring(separated_sp, n_positions=n)

        # Mesma contagem de Sp/Rp, mas acoplamento diferente
        assert result_adj["sp_count"] == result_sep["sp_count"]
        assert result_adj["rp_count"] == result_sep["rp_count"]

        # Adjacentes Sp tem acoplamento J negativo (estabiliza)
        assert result_adj["coupling_term"] < 0
        # Separados nao tem pares Sp-Sp adjacentes
        assert result_sep["coupling_term"] == 0.0

    def test_bitstring_evaluation_manual_calculation(self) -> None:
        """Verifica calculo manual passo a passo para bitstring especifica."""
        n = 5
        bitstring = [1, 0, 1, 1, 0]  # Sp, Rp, Sp, Sp, Rp

        result = evaluate_bitstring(bitstring, n_positions=n)

        # Correcao individual: 3*Sp + 2*Rp
        correction = 3 * SP_CORRECTION + 2 * RP_CORRECTION
        # = 3 * (-0.2) + 2 * (0.1) = -0.6 + 0.2 = -0.4

        # Acoplamento J: x0*x1 + x1*x2 + x2*x3 + x3*x4
        # = 1*0 + 0*1 + 1*1 + 1*0 = 0 + 0 + 1 + 0 = 1
        coupling = 1 * J_COUPLING  # = -0.05

        # Penalidade Rp-Rp pares: (1-x0)(1-x1) + ... so pares (0,0)
        # Posicoes Rp: 1, 4 -> nenhum par adjacente Rp-Rp
        pen_per_pair = RP_CONSECUTIVE_PENALTY / 2.0
        penalty = 0.0  # Nenhum par Rp adjacente

        expected_dg = DG_BASE + correction + coupling + penalty

        assert result["dg_total"] == pytest.approx(expected_dg, abs=1e-3)
        assert result["correction_term"] == pytest.approx(correction, abs=1e-3)
        assert result["coupling_term"] == pytest.approx(coupling, abs=1e-3)
        assert result["penalty_term"] == pytest.approx(penalty, abs=1e-3)
        assert result["sp_count"] == 3
        assert result["rp_count"] == 2

    def test_bitstring_length_validation(self) -> None:
        """Bitstring com tamanho errado deve gerar ValueError."""
        with pytest.raises(ValueError, match="esperado"):
            evaluate_bitstring([1, 0, 1], n_positions=5)

    def test_full_25_positions(self) -> None:
        """Avaliacao com 25 posicoes (tamanho real do MRL-ASO-001)."""
        # Configuracao mista: primeiros 15 Sp, ultimos 10 Rp
        bitstring = [1] * 15 + [0] * 10

        result = evaluate_bitstring(bitstring, n_positions=25)

        assert result["sp_count"] == 15
        assert result["rp_count"] == 10
        assert result["dg_total"] < 0  # dG deve ser negativo
        assert len(result["configuration"]) == 50  # 25 posicoes * 2 chars each


class TestBuildQubo:
    """Testes para a construcao do QuadraticProgram."""

    def test_qubo_has_correct_variables(self) -> None:
        """QUBO deve ter n variaveis binarias."""
        n = 5
        qp = build_qubo(n_positions=n)
        assert qp.get_num_vars() == n
        assert qp.get_num_binary_vars() == n

    def test_qubo_is_minimization(self) -> None:
        """QUBO deve ser um problema de minimizacao."""
        qp = build_qubo(n_positions=5)
        from qiskit_optimization.problems.quadratic_objective import ObjSense
        assert qp.objective.sense == ObjSense.MINIMIZE

    def test_qubo_fval_matches_evaluate(self) -> None:
        """O fval do QUBO deve coincidir com evaluate_bitstring."""
        n = 5
        qp = build_qubo(n_positions=n)

        # Testar varias configuracoes
        for bits in [[0] * n, [1] * n, [1, 0, 1, 0, 1], [0, 1, 0, 1, 0]]:
            x = np.array(bits, dtype=float)
            fval = qp.objective.evaluate(x)
            eval_result = evaluate_bitstring(bits, n_positions=n)

            assert fval == pytest.approx(eval_result["dg_total"], abs=1e-3), (
                f"Divergencia para bitstring {bits}: "
                f"QUBO fval={fval}, evaluate dG={eval_result['dg_total']}"
            )


class TestBitstringToConfig:
    """Testes para conversao bitstring -> configuracao Rp/Sp."""

    def test_basic_conversion(self) -> None:
        """Conversao basica de bitstring para lista de strings."""
        config = bitstring_to_config([1, 0, 1, 0])
        assert config == ["Sp", "Rp", "Sp", "Rp"]

    def test_all_sp(self) -> None:
        """Todas as posicoes Sp."""
        config = bitstring_to_config([1, 1, 1])
        assert all(c == "Sp" for c in config)

    def test_all_rp(self) -> None:
        """Todas as posicoes Rp."""
        config = bitstring_to_config([0, 0, 0])
        assert all(c == "Rp" for c in config)


class TestPositionLabel:
    """Testes para labels de posicao."""

    def test_first_position(self) -> None:
        """Primeira posicao do ASO."""
        from mrl_quantum.config import ASO_SEQUENCE
        label = get_aso_position_label(0)
        assert label == f"{ASO_SEQUENCE[0]}0"

    def test_out_of_range(self) -> None:
        """Posicao fora do range retorna '?'."""
        label = get_aso_position_label(100)
        assert label.startswith("?")
