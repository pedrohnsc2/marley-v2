"""Formulacao QUBO para otimizacao de estereoisomeros fosforotioato.

Cada posicao i do ASO MRL-ASO-001 (25 nt) possui uma variavel binaria:
    x_i = 0 -> configuracao Rp
    x_i = 1 -> configuracao Sp

Funcao custo (minimizacao):
    C(x) = dG_BASE
           + sum_i correction(x_i)
           + sum_{i,i+1} J * x_i * x_{i+1}
           + sum_k penalty_consecutive_Rp(k)

Onde:
    - correction(x_i=1) = SP_CORRECTION = -0.2 kcal/mol (Sp estabiliza)
    - correction(x_i=0) = RP_CORRECTION = +0.1 kcal/mol (Rp desestabiliza)
    - J = -0.05 kcal/mol (cooperatividade PS entre posicoes adjacentes)
    - Penalidade: +0.5 kcal/mol se 3+ Rp consecutivos (reduz atividade RNase H)

Referencia: Eckstein F (2014) Nucleic Acids Ther 24(6):374-387
"""

from __future__ import annotations

import logging
from typing import Final

import numpy as np
from qiskit_optimization import QuadraticProgram

from mrl_quantum.config import (
    ASO_SEQUENCE,
    DG_BASE,
    J_COUPLING,
    N_POSITIONS,
    RP_CONSECUTIVE_PENALTY,
    RP_CORRECTION,
    SP_CORRECTION,
)

logger = logging.getLogger("mrl_quantum.qaoa.qubo")

# ---------------------------------------------------------------------------
# Constantes derivadas para a formulacao QUBO
# ---------------------------------------------------------------------------

# Correcao linear: correction(x_i) = SP_CORRECTION * x_i + RP_CORRECTION * (1 - x_i)
# = RP_CORRECTION + (SP_CORRECTION - RP_CORRECTION) * x_i
# Parte constante (RP_CORRECTION para todas as posicoes) vai no offset
_LINEAR_COEFF: Final[float] = SP_CORRECTION - RP_CORRECTION  # -0.3


def build_qubo(
    n_positions: int | None = None,
    dg_base: float | None = None,
    sp_correction: float | None = None,
    rp_correction: float | None = None,
    j_coupling: float | None = None,
    rp_consecutive_penalty: float | None = None,
) -> QuadraticProgram:
    """Constroi o QuadraticProgram para o problema de estereoisomeros PS.

    Minimiza dG total (quanto mais negativo, mais estavel o duplex).

    Args:
        n_positions: Numero de posicoes (default: N_POSITIONS = 25).
        dg_base: dG base do duplex sem correcao PS (default: DG_BASE).
        sp_correction: Correcao energetica Sp (default: SP_CORRECTION).
        rp_correction: Correcao energetica Rp (default: RP_CORRECTION).
        j_coupling: Acoplamento J entre vizinhos (default: J_COUPLING).
        rp_consecutive_penalty: Penalidade por 3+ Rp consecutivos
            (default: RP_CONSECUTIVE_PENALTY).

    Returns:
        QuadraticProgram configurado para minimizacao.
    """
    n = n_positions if n_positions is not None else N_POSITIONS
    dg = dg_base if dg_base is not None else DG_BASE
    sp_corr = sp_correction if sp_correction is not None else SP_CORRECTION
    rp_corr = rp_correction if rp_correction is not None else RP_CORRECTION
    j_coup = j_coupling if j_coupling is not None else J_COUPLING
    rp_pen = rp_consecutive_penalty if rp_consecutive_penalty is not None else RP_CONSECUTIVE_PENALTY

    logger.info(
        "Construindo QUBO: n=%d, dG_base=%.2f, Sp=%.2f, Rp=%.2f, J=%.3f, pen=%.2f",
        n, dg, sp_corr, rp_corr, j_coup, rp_pen,
    )

    qp = QuadraticProgram(f"ps_stereoisomer_{n}pos")

    # Criar variaveis binarias: x_i in {0=Rp, 1=Sp}
    var_names = [f"x{i}" for i in range(n)]
    for name in var_names:
        qp.binary_var(name)

    # -----------------------------------------------------------------------
    # Montar coeficientes lineares e quadraticos
    # -----------------------------------------------------------------------
    # correction(x_i) = rp_corr + (sp_corr - rp_corr) * x_i
    # O offset constante (dg + n * rp_corr) nao afeta a otimizacao,
    # mas incluimos via constant para que fval reflita dG total.

    linear_coeff = sp_corr - rp_corr  # contribuicao linear por variavel

    linear: dict[str, float] = {}
    for name in var_names:
        linear[name] = linear_coeff

    # Acoplamento J entre posicoes adjacentes: J * x_i * x_{i+1}
    quadratic: dict[tuple[str, str], float] = {}
    for i in range(n - 1):
        quadratic[(var_names[i], var_names[i + 1])] = j_coup

    # -----------------------------------------------------------------------
    # Penalidade por 3+ Rp consecutivos
    # -----------------------------------------------------------------------
    # Rp = x_i = 0, entao "3 Rp consecutivos" = (1-x_i)(1-x_{i+1})(1-x_{i+2}) >= 1
    # Expandindo: 1 - x_i - x_{i+1} - x_{i+2} + x_i*x_{i+1} + x_i*x_{i+2} + x_{i+1}*x_{i+2} - x_i*x_{i+1}*x_{i+2}
    # Para QUBO (quadratico), usamos aproximacao:
    # Penalizamos quando x_i + x_{i+1} + x_{i+2} <= 0 (todos Rp)
    # Termo de penalidade: pen * (1 - x_i) * (1 - x_{i+1}) nao captura triplas.
    #
    # Abordagem: introduzir penalidade suave para pares Rp-Rp adjacentes,
    # escalonada para que 3+ consecutivos sejam mais penalizados.
    # pen_pair = rp_pen * (1-x_i)(1-x_{i+1}) = rp_pen * (1 - x_i - x_{i+1} + x_i*x_{i+1})
    # Isso penaliza cada par Rp-Rp; 3 consecutivos geram 2 pares = 2 * rp_pen.
    #
    # Para manter o problema estritamente QUBO (sem termos cubicos),
    # usamos penalidade por par Rp-Rp com peso escalonado.
    pen_per_pair = rp_pen / 2.0  # Cada par contribui metade; 3 consec = 2 pares = rp_pen

    for i in range(n - 1):
        # pen_per_pair * (1 - x_i)(1 - x_{i+1})
        # = pen_per_pair * (1 - x_i - x_{i+1} + x_i*x_{i+1})
        # Constante: pen_per_pair (vai no offset)
        # Linear: -pen_per_pair * x_i, -pen_per_pair * x_{i+1}
        # Quadratico: +pen_per_pair * x_i * x_{i+1}
        linear[var_names[i]] = linear.get(var_names[i], 0.0) - pen_per_pair
        linear[var_names[i + 1]] = linear.get(var_names[i + 1], 0.0) - pen_per_pair
        key = (var_names[i], var_names[i + 1])
        quadratic[key] = quadratic.get(key, 0.0) + pen_per_pair

    # -----------------------------------------------------------------------
    # Offset constante: dG_base + n * rp_corr + (n-1) * pen_per_pair
    # -----------------------------------------------------------------------
    constant = dg + n * rp_corr + (n - 1) * pen_per_pair

    qp.minimize(linear=linear, quadratic=quadratic, constant=constant)

    logger.info(
        "QUBO construido: %d variaveis, %d termos lineares, %d termos quadraticos",
        n, len(linear), len(quadratic),
    )

    return qp


def evaluate_bitstring(
    bitstring: list[int] | np.ndarray,
    n_positions: int | None = None,
    dg_base: float | None = None,
    sp_correction: float | None = None,
    rp_correction: float | None = None,
    j_coupling: float | None = None,
    rp_consecutive_penalty: float | None = None,
) -> dict[str, object]:
    """Avalia uma configuracao de estereoisomeros e calcula o dG total.

    Args:
        bitstring: Vetor binario de tamanho n_positions.
            x_i = 1 -> Sp, x_i = 0 -> Rp.
        n_positions: Numero de posicoes (default: N_POSITIONS).
        dg_base: dG base (default: DG_BASE).
        sp_correction: Correcao Sp (default: SP_CORRECTION).
        rp_correction: Correcao Rp (default: RP_CORRECTION).
        j_coupling: Acoplamento J (default: J_COUPLING).
        rp_consecutive_penalty: Penalidade Rp (default: RP_CONSECUTIVE_PENALTY).

    Returns:
        Dicionario com dG total, configuracao Rp/Sp, e decomposicao de termos.
    """
    n = n_positions if n_positions is not None else N_POSITIONS
    dg = dg_base if dg_base is not None else DG_BASE
    sp_corr = sp_correction if sp_correction is not None else SP_CORRECTION
    rp_corr = rp_correction if rp_correction is not None else RP_CORRECTION
    j_coup = j_coupling if j_coupling is not None else J_COUPLING
    rp_pen = rp_consecutive_penalty if rp_consecutive_penalty is not None else RP_CONSECUTIVE_PENALTY

    x = np.asarray(bitstring, dtype=int)
    if len(x) != n:
        raise ValueError(
            f"Bitstring tem {len(x)} posicoes, esperado {n}"
        )

    # Termo de correcao individual
    correction_term = 0.0
    for i in range(n):
        if x[i] == 1:
            correction_term += sp_corr
        else:
            correction_term += rp_corr

    # Termo de acoplamento J entre posicoes adjacentes
    coupling_term = 0.0
    for i in range(n - 1):
        coupling_term += j_coup * x[i] * x[i + 1]

    # Penalidade por Rp consecutivos (usando mesma formula do QUBO)
    pen_per_pair = rp_pen / 2.0
    penalty_term = 0.0
    for i in range(n - 1):
        # (1 - x_i)(1 - x_{i+1}) = 1 se ambos Rp
        penalty_term += pen_per_pair * (1 - x[i]) * (1 - x[i + 1])

    # Contar tripletas Rp consecutivas (metrica informativa)
    rp_triples = 0
    for i in range(n - 2):
        if x[i] == 0 and x[i + 1] == 0 and x[i + 2] == 0:
            rp_triples += 1

    dg_total = dg + correction_term + coupling_term + penalty_term

    # Montar configuracao Rp/Sp
    config_str = "".join("Sp" if xi == 1 else "Rp" for xi in x)
    sp_count = int(np.sum(x))
    rp_count = n - sp_count

    return {
        "dg_total": round(float(dg_total), 4),
        "dg_base": dg,
        "correction_term": round(float(correction_term), 4),
        "coupling_term": round(float(coupling_term), 4),
        "penalty_term": round(float(penalty_term), 4),
        "configuration": config_str,
        "bitstring": x.tolist(),
        "sp_count": sp_count,
        "rp_count": rp_count,
        "rp_triples": rp_triples,
    }


def bitstring_to_config(bitstring: list[int] | np.ndarray) -> list[str]:
    """Converte bitstring para lista de configuracoes Rp/Sp.

    Args:
        bitstring: Vetor binario (0=Rp, 1=Sp).

    Returns:
        Lista de strings 'Rp' ou 'Sp' para cada posicao.
    """
    return ["Sp" if int(xi) == 1 else "Rp" for xi in bitstring]


def get_aso_position_label(index: int) -> str:
    """Retorna o label da posicao com nucleotideo e indice.

    Args:
        index: Indice 0-based da posicao no ASO.

    Returns:
        String no formato "A0", "C1", etc.
    """
    if 0 <= index < len(ASO_SEQUENCE):
        return f"{ASO_SEQUENCE[index]}{index}"
    return f"?{index}"
