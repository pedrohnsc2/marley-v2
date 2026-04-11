"""Cadeia de Markov para transicoes mutacionais no alvo de 25 nt.

Modela cada posicao do alvo SL RNA como um processo de Markov independente
com 4 estados (A, C, G, T). As taxas de transicao refletem a biologia
de trypanosomatideos:
  - Transicoes (purina<->purina, pirimidina<->pirimidina): taxa base * 2
  - Transversoes (purina<->pirimidina): taxa base * 1

O espaco de estados completo do alvo de 25 nt tem 4^25 ~ 10^15 estados,
impossivel de enumerar. Modelamos como 25 processos independentes de 4 estados.

Ref: Rogers MB et al. (2011) PLoS Genetics 7(8):e1002237 — taxa de mutacao
Ref: Kimura M (1980) J Mol Evol 16:111-120 — modelo transicao/transversao
"""

from __future__ import annotations

from typing import Any

import numpy as np

from aso_math.config import BASES, MUTATION_RATE


# ---------------------------------------------------------------------------
# Classificacao de substituicoes (transicao vs transversao)
# ---------------------------------------------------------------------------

# Purinas e pirimidinas
PURINES: frozenset[str] = frozenset({"A", "G"})
PYRIMIDINES: frozenset[str] = frozenset({"C", "T"})

# Razao transicao/transversao (ti/tv)
# Ref: Kimura (1980); tipico para genomas eucarioticos: 2.0-3.0
# Trypanosomatideos: ~2.0 (Rogers et al. 2011)
TI_TV_RATIO: float = 2.0


def is_transition(base_from: str, base_to: str) -> bool:
    """Verifica se a substituicao e uma transicao (purina<->purina ou piri<->piri).

    Transicoes (A<->G, C<->T) sao ~2x mais frequentes que transversoes.
    """
    both_purines = base_from in PURINES and base_to in PURINES
    both_pyrimidines = base_from in PYRIMIDINES and base_to in PYRIMIDINES
    return both_purines or both_pyrimidines


# ---------------------------------------------------------------------------
# Matriz de taxas de transicao Q para uma posicao
# ---------------------------------------------------------------------------


def build_rate_matrix(
    mutation_rate: float = MUTATION_RATE,
    ti_tv_ratio: float = TI_TV_RATIO,
) -> np.ndarray:
    """Constroi a matriz de taxas Q (4x4) para uma posicao.

    A matriz Q para um processo de Markov em tempo continuo:
      Q[i,j] = taxa de transicao do estado i para o estado j (i != j)
      Q[i,i] = -sum(Q[i,j] para j != i)  (linhas somam zero)

    As taxas sao normalizadas para que a taxa total de saida de qualquer
    estado reflita a taxa de mutacao observada.

    Modelo de substituicao:
      - Transicao (A<->G, C<->T): mu * alpha
      - Transversao (A<->C, A<->T, G<->C, G<->T): mu * beta
      - Onde alpha/beta = ti_tv_ratio

    Para normalizar: cada base tem 1 transicao e 2 transversoes possiveis
      taxa_total = mu * (alpha + 2*beta)
    Queremos que taxa_total = mutation_rate (taxa total de saida por base)
      => beta = mutation_rate / (ti_tv_ratio + 2)
      => alpha = ti_tv_ratio * beta

    Args:
        mutation_rate: Taxa de mutacao total por base por geracao.
        ti_tv_ratio: Razao transicao/transversao.

    Returns:
        Matriz Q (4x4) com taxas por geracao. Linhas: A=0, C=1, G=2, T=3.
    """
    bases = list(BASES)  # ["A", "C", "G", "T"]
    n = len(bases)

    # Normalizar taxas: total de saida = mutation_rate
    # Cada base tem 1 transicao + 2 transversoes possiveis
    beta = mutation_rate / (ti_tv_ratio + 2.0)
    alpha = ti_tv_ratio * beta

    q_matrix = np.zeros((n, n), dtype=np.float64)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if is_transition(bases[i], bases[j]):
                q_matrix[i, j] = alpha
            else:
                q_matrix[i, j] = beta

    # Elementos diagonais: soma das linhas = 0
    for i in range(n):
        q_matrix[i, i] = -np.sum(q_matrix[i, :])

    return q_matrix


def transition_probability_matrix(
    q_matrix: np.ndarray,
    n_generations: float,
) -> np.ndarray:
    """Calcula P(t) = exp(Q*t) — probabilidades de transicao apos t geracoes.

    Usa a decomposicao espectral para calcular a exponencial de matriz:
      P(t) = exp(Q*t)

    Para tempos pequenos (Q*t << 1), P(t) ~ I + Q*t (aproximacao linear).
    Para tempos grandes, usamos scipy-free eigendecomposition.

    Args:
        q_matrix: Matriz de taxas Q (4x4).
        n_generations: Numero de geracoes.

    Returns:
        Matriz de probabilidades de transicao P(t) (4x4).
    """
    # Decomposicao espectral: Q = V * diag(eigenvalues) * V^-1
    # P(t) = V * diag(exp(eigenvalues * t)) * V^-1
    eigenvalues, eigenvectors = np.linalg.eig(q_matrix)

    # Exponencial dos autovalores multiplicados pelo tempo
    exp_eigenvalues = np.diag(np.exp(eigenvalues * n_generations))

    # P(t) = V * exp(Lambda*t) * V^-1
    p_matrix = eigenvectors @ exp_eigenvalues @ np.linalg.inv(eigenvectors)

    # Garantir que valores sao reais (artefato numerico de eigendecomposition)
    p_matrix = np.real(p_matrix)

    # Clipar para [0, 1] — artefatos numericos podem dar valores ligeiramente negativos
    p_matrix = np.clip(p_matrix, 0.0, 1.0)

    # Renormalizar linhas para somar 1.0
    row_sums = p_matrix.sum(axis=1, keepdims=True)
    p_matrix = p_matrix / row_sums

    return p_matrix


# ---------------------------------------------------------------------------
# Taxa de mutacao efetiva por posicao, incluindo transicao vs transversao
# ---------------------------------------------------------------------------


def compute_position_mutation_rates(
    target_sequence: str,
    mutation_rate: float = MUTATION_RATE,
    ti_tv_ratio: float = TI_TV_RATIO,
) -> list[dict[str, Any]]:
    """Calcula a taxa de cada mutacao especifica em cada posicao do alvo.

    Para cada posicao (25) e cada base alternativa (3), calcula a taxa
    considerando se e transicao ou transversao.

    Args:
        target_sequence: Sequencia alvo de 25 nt no SL RNA.
        mutation_rate: Taxa de mutacao total por base por geracao.
        ti_tv_ratio: Razao transicao/transversao.

    Returns:
        Lista de dicionarios com taxa por mutacao especifica.
    """
    bases = list(BASES)

    # Normalizar: taxa_total_saida = mutation_rate por posicao
    beta = mutation_rate / (ti_tv_ratio + 2.0)
    alpha = ti_tv_ratio * beta

    position_rates: list[dict[str, Any]] = []

    for pos, wt_base in enumerate(target_sequence.upper()):
        for alt_base in bases:
            if alt_base == wt_base:
                continue

            is_ti = is_transition(wt_base, alt_base)
            rate = alpha if is_ti else beta

            position_rates.append({
                "position": pos,
                "wildtype_base": wt_base,
                "mutant_base": alt_base,
                "substitution_type": "transition" if is_ti else "transversion",
                "rate_per_generation": rate,
            })

    return position_rates


def compute_stationary_distribution(q_matrix: np.ndarray) -> np.ndarray:
    """Calcula a distribuicao estacionaria pi da cadeia de Markov.

    A distribuicao estacionaria satisfaz pi * Q = 0 e sum(pi) = 1.
    Para o modelo HKY/Kimura com taxas simetricas, pi = [0.25, 0.25, 0.25, 0.25].

    Args:
        q_matrix: Matriz de taxas Q (4x4).

    Returns:
        Vetor de distribuicao estacionaria (4,).
    """
    n = q_matrix.shape[0]

    # Resolver pi * Q = 0 com restricao sum(pi) = 1
    # Substituir a ultima equacao por sum(pi) = 1
    a_matrix = q_matrix.T.copy()
    a_matrix[-1, :] = 1.0

    b_vector = np.zeros(n)
    b_vector[-1] = 1.0

    pi = np.linalg.solve(a_matrix, b_vector)

    # Garantir que nao ha valores negativos (artefato numerico)
    pi = np.clip(pi, 0.0, 1.0)
    pi /= pi.sum()

    return pi
