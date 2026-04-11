"""Custos de fitness baseados em conservacao evolutiva do SL RNA.

O SL RNA de Leishmania e conservado ha ~500 milhoes de anos entre todos
os trypanosomatideos. TODA posicao na regiao alvo do ASO e 100% conservada
entre especies divergentes (L. infantum, L. major, L. donovani, T. brucei,
T. cruzi). Essa conservacao extrema implica selecao purificadora fortissima:
qualquer mutacao nessas posicoes e letal ou quasi-letal.

Modelo de fitness:
  w(mutation) = exp(-conservation_score * penalty_factor)

Onde:
  conservation_score: 1.0 para posicoes 100% conservadas (todas no SL RNA)
  penalty_factor: calibrado pela conservacao de ~500 Ma

A taxa de mutacao efetiva incorpora o custo de fitness:
  mu_eff(pos, base) = mu_raw(pos, base) * w(pos)

Ref: Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
Ref: Michaeli S (2011) Parasitology 138(12):1473-1491
Ref: Kimura M (1983) The Neutral Theory of Molecular Evolution, cap. 3
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from aso_math.config import (
    ASO_TARGET_END,
    ASO_TARGET_START,
    SL_SEQUENCE,
)


# ---------------------------------------------------------------------------
# Constantes de conservacao
# ---------------------------------------------------------------------------

# Score de conservacao por posicao do SL RNA alvo
# Para Leishmania: TODAS as posicoes sao 100% conservadas entre especies
# (entropia de Shannon = 0 em alinhamento multi-species).
# Score de 1.0 = totalmente conservada, 0.0 = sem restricao.
DEFAULT_CONSERVATION_SCORE: float = 1.0

# Fator de penalidade: calibrado pela observacao de que ZERO mutacoes
# foram fixadas em ~500 Ma de evolucao divergente.
# Se assumirmos que a populacao efetiva ancestral era ~10^6 e a taxa
# de mutacao e ~2e-9 por base por geracao, o numero esperado de mutacoes
# neutras fixadas em 500 Ma seria:
#   E[fixacoes] = 2 * Ne * mu * t_geracoes ~ 2 * 10^6 * 2e-9 * 4e11 ~ 1600
# Observamos ZERO. Isso implica selecao purificadora com s >> 1/(2*Ne).
# Conservadoramente, usamos penalty_factor = 10 (moderado).
# Com score=1.0: w = exp(-10) ~ 4.5e-5 (reducao de 99.995% na fitness)
PENALTY_FACTOR_DEFAULT: float = 10.0

# Cenario generoso (para sensibilidade): penalidade menor
PENALTY_FACTOR_GENEROUS: float = 3.0

# Cenario ultra-generoso: sem penalidade de fitness (fisicamente impossivel
# dado que SL RNA e essencial, mas usado como limite superior absoluto)
PENALTY_FACTOR_NONE: float = 0.0


# ---------------------------------------------------------------------------
# Modelo de fitness por posicao
# ---------------------------------------------------------------------------


def compute_conservation_scores(
    sl_sequence: str = SL_SEQUENCE,
    target_start: int = ASO_TARGET_START,
    target_end: int = ASO_TARGET_END,
    related_sequences: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Calcula o score de conservacao para cada posicao do alvo.

    Se sequencias de especies relacionadas forem fornecidas, calcula a
    entropia de Shannon em cada posicao. Caso contrario, usa o score
    default de 1.0 (100% conservado) baseado na literatura.

    A entropia de Shannon H em cada posicao e:
      H = -sum(f_i * log2(f_i)) para cada base i com f_i > 0

    O score de conservacao e: 1 - H/2 (normalizado para [0, 1])
    H_max = log2(4) = 2.0 para 4 bases

    Args:
        sl_sequence: Sequencia do SL RNA completo.
        target_start: Posicao inicial do alvo (0-indexed).
        target_end: Posicao final do alvo (0-indexed, exclusivo).
        related_sequences: Dict de {especie: sequencia_SL} para alinhamento.

    Returns:
        Lista com score de conservacao por posicao.
    """
    target = sl_sequence[target_start:target_end].upper()
    n_positions = len(target)

    scores: list[dict[str, Any]] = []

    if related_sequences and len(related_sequences) >= 2:
        # Calcular entropia de Shannon a partir do alinhamento
        # Assumimos que todas as sequencias estao alinhadas na mesma regiao
        all_seqs = [target] + [
            seq[target_start:target_end].upper()
            for seq in related_sequences.values()
        ]
        n_seqs = len(all_seqs)

        for pos in range(n_positions):
            # Contar frequencia de cada base nesta posicao
            base_counts: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0}
            for seq in all_seqs:
                if pos < len(seq):
                    base = seq[pos]
                    if base in base_counts:
                        base_counts[base] += 1

            # Entropia de Shannon
            total = sum(base_counts.values())
            entropy = 0.0
            for count in base_counts.values():
                if count > 0:
                    freq = count / total
                    entropy -= freq * math.log2(freq)

            # Score de conservacao: 1 - H/H_max
            h_max = 2.0  # log2(4)
            conservation = 1.0 - (entropy / h_max)

            scores.append({
                "position": pos,
                "sl_position": target_start + pos,
                "base": target[pos],
                "shannon_entropy": round(entropy, 6),
                "conservation_score": round(conservation, 6),
                "n_sequences_aligned": n_seqs,
            })
    else:
        # Sem sequencias de referencia: usar valor default da literatura
        # TODAS as posicoes do SL RNA de Leishmania sao conservadas
        for pos in range(n_positions):
            scores.append({
                "position": pos,
                "sl_position": target_start + pos,
                "base": target[pos],
                "shannon_entropy": 0.0,
                "conservation_score": DEFAULT_CONSERVATION_SCORE,
                "n_sequences_aligned": 1,
                "note": "Default: 100% conservado entre Leishmania spp. (literatura)",
            })

    return scores


def compute_fitness_cost(
    conservation_score: float,
    penalty_factor: float = PENALTY_FACTOR_DEFAULT,
) -> float:
    """Calcula o custo de fitness de uma mutacao baseado na conservacao.

    Modelo exponencial:
      w = exp(-conservation_score * penalty_factor)

    Interpretacao:
      - w proximo de 1.0: mutacao e quase neutra (posicao nao conservada)
      - w proximo de 0.0: mutacao e letal (posicao conservada)

    Para o SL RNA (conservation_score = 1.0, penalty = 10.0):
      w = exp(-10) = 4.54e-5

    Args:
        conservation_score: Score de conservacao da posicao [0, 1].
        penalty_factor: Fator de penalidade.

    Returns:
        Fitness relativa w em [0, 1]. Valores menores = mais deleteria.
    """
    return math.exp(-conservation_score * penalty_factor)


def compute_effective_mutation_rates(
    position_mutation_rates: list[dict[str, Any]],
    conservation_scores: list[dict[str, Any]],
    penalty_factor: float = PENALTY_FACTOR_DEFAULT,
) -> list[dict[str, Any]]:
    """Calcula taxas de mutacao efetivas (corrigidas por fitness).

    mu_eff(pos, base) = mu_raw(pos, base) * fitness(pos)

    A fitness reduz a taxa efetiva porque mutacoes em posicoes conservadas
    sao removidas por selecao purificadora antes de atingir frequencia
    apreciavel na populacao.

    Args:
        position_mutation_rates: Taxas brutas por posicao e base (de markov_chain).
        conservation_scores: Scores de conservacao por posicao.
        penalty_factor: Fator de penalidade para calculo de fitness.

    Returns:
        Lista com taxas efetivas por mutacao.
    """
    # Indexar scores por posicao para lookup rapido
    score_by_pos: dict[int, float] = {
        s["position"]: s["conservation_score"]
        for s in conservation_scores
    }

    effective_rates: list[dict[str, Any]] = []

    for rate_entry in position_mutation_rates:
        pos = rate_entry["position"]
        raw_rate = rate_entry["rate_per_generation"]

        # Score de conservacao da posicao
        cons_score = score_by_pos.get(pos, DEFAULT_CONSERVATION_SCORE)

        # Fitness da mutacao nesta posicao
        fitness = compute_fitness_cost(cons_score, penalty_factor)

        # Taxa efetiva
        effective_rate = raw_rate * fitness

        effective_rates.append({
            **rate_entry,
            "conservation_score": cons_score,
            "fitness_cost": round(1.0 - fitness, 10),
            "fitness_retained": fitness,
            "penalty_factor": penalty_factor,
            "effective_rate_per_generation": effective_rate,
        })

    return effective_rates


def summarize_fitness_landscape(
    conservation_scores: list[dict[str, Any]],
    penalty_factor: float = PENALTY_FACTOR_DEFAULT,
) -> dict[str, Any]:
    """Resumo estatistico do custo de fitness no alvo inteiro.

    Calcula a fitness media, minima e maxima para o alvo de 25 nt,
    e o fator de reducao total na taxa de mutacao efetiva.

    Args:
        conservation_scores: Scores de conservacao por posicao.
        penalty_factor: Fator de penalidade.

    Returns:
        Dicionario com resumo estatistico.
    """
    fitness_values = [
        compute_fitness_cost(s["conservation_score"], penalty_factor)
        for s in conservation_scores
    ]

    mean_fitness = float(np.mean(fitness_values))
    min_fitness = float(np.min(fitness_values))
    max_fitness = float(np.max(fitness_values))

    # Fator de reducao: quanto a fitness reduz a taxa de mutacao
    # Se todas as posicoes tem fitness = f, reducao = f
    reduction_factor = mean_fitness

    return {
        "n_positions": len(conservation_scores),
        "penalty_factor": penalty_factor,
        "mean_fitness": mean_fitness,
        "min_fitness": min_fitness,
        "max_fitness": max_fitness,
        "mean_conservation_score": float(np.mean([
            s["conservation_score"] for s in conservation_scores
        ])),
        "mutation_rate_reduction_factor": reduction_factor,
        "interpretation": (
            f"Mutacoes no alvo tem fitness media de {mean_fitness:.2e} "
            f"(reducao de {1.0/reduction_factor:.0e}x na taxa efetiva). "
            f"Posicoes 100% conservadas tornam virtualmente toda mutacao letal."
            if reduction_factor < 0.01
            else f"Fitness media de {mean_fitness:.4f} — restricao moderada."
        ),
    }
