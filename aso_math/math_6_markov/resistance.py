"""Modelos de Poisson e Kimura para emergencia e fixacao de resistencia.

Combina:
  1. Analise de escape: quais mutacoes rompem binding do ASO?
  2. Poisson: probabilidade de emergencia em populacao finita
  3. Kimura: probabilidade de fixacao dada selecao negativa
  4. Multi-mutante: caminhos de escape com 2+ mutacoes simultaneas
  5. Tempo esperado de resistencia

Ref: Kimura M (1962) Genetics 47:713-719 — fixacao
Ref: Kimura M (1980) J Mol Evol 16:111-120 — modelo neutro
Ref: Ross SM (2014) Introduction to Probability Models — Poisson
Ref: Crooke ST et al. (2017) Nucleic Acids Res — limiar de binding ASO
"""

from __future__ import annotations

import math
from itertools import combinations
from typing import Any

import numpy as np

from aso_math.config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    BASES,
    DG_FUNCTIONAL_THRESHOLD,
    GENERATION_TIME_HOURS,
    MUTATION_RATE,
    SL_RNA_COPY_NUMBER,
    SL_SEQUENCE,
)
from aso_math.thermo import compute_dg

# Importar funcoes do modelo existente (05) para reusar calculo de dG disruptido
from aso_math.config import NN_INIT, T_PHYSIOLOGICAL
from aso_math.thermo import _get_nn_params


# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

HOURS_PER_YEAR: float = 365.25 * 24.0

# Tamanhos de populacao para modelagem
POPULATION_SIZES: list[float] = [1e3, 1e6, 1e8, 1e10]

# Duracoes de tratamento em semanas
TREATMENT_DURATIONS_WEEKS: list[int] = [1, 2, 4, 8, 12, 26, 52]


# ---------------------------------------------------------------------------
# 1. Calculo de dG disruptido (reutiliza logica do modulo 05)
# ---------------------------------------------------------------------------


def compute_disrupted_dg(aso_seq: str, mismatch_position: int) -> float:
    """Calcula dG do ASO com mismatch na posicao indicada.

    Reimplementacao local da funcao do modulo 05 para manter independencia.
    Zera a contribuicao dos dinucleotideos adjacentes ao mismatch.
    Conservador: mismatches reais tem dG positivo (desestabilizante).

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        mismatch_position: Indice 0-based da posicao com mismatch no ASO.

    Returns:
        dG em kcal/mol com o mismatch.
    """
    seq = aso_seq.upper()
    n = len(seq)

    # Posicoes de dinucleotideos afetados
    disrupted: set[int] = set()
    if mismatch_position > 0:
        disrupted.add(mismatch_position - 1)
    if mismatch_position < n - 1:
        disrupted.add(mismatch_position)

    total_dh = NN_INIT[0]
    total_ds = NN_INIT[1]

    for i in range(n - 1):
        if i in disrupted:
            continue
        dh, ds = _get_nn_params(seq[i:i + 2])
        total_dh += dh
        total_ds += ds

    delta_g = total_dh - T_PHYSIOLOGICAL * (total_ds / 1000.0)
    return round(delta_g, 2)


def compute_multi_mismatch_dg(
    aso_seq: str,
    mismatch_positions: list[int],
) -> float:
    """Calcula dG do ASO com mismatches em MULTIPLAS posicoes.

    Extensao do modelo de mismatch unico: zera a contribuicao de todos
    os dinucleotideos adjacentes a qualquer posicao de mismatch.

    Para k mismatches, ate 2k dinucleotideos sao afetados (menos se
    mismatches sao adjacentes — overlap).

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        mismatch_positions: Lista de indices 0-based com mismatch.

    Returns:
        dG em kcal/mol com todos os mismatches.
    """
    seq = aso_seq.upper()
    n = len(seq)

    # Coletar todos os dinucleotideos afetados
    disrupted: set[int] = set()
    for mm_pos in mismatch_positions:
        if mm_pos > 0:
            disrupted.add(mm_pos - 1)
        if mm_pos < n - 1:
            disrupted.add(mm_pos)

    total_dh = NN_INIT[0]
    total_ds = NN_INIT[1]

    for i in range(n - 1):
        if i in disrupted:
            continue
        dh, ds = _get_nn_params(seq[i:i + 2])
        total_dh += dh
        total_ds += ds

    delta_g = total_dh - T_PHYSIOLOGICAL * (total_ds / 1000.0)
    return round(delta_g, 2)


# ---------------------------------------------------------------------------
# 2. Escape mutation analysis (posicao a posicao)
# ---------------------------------------------------------------------------


def analyze_escape_mutations(
    aso_seq: str = ASO_SEQUENCE,
    sl_seq: str = SL_SEQUENCE,
    target_start: int = ASO_TARGET_START,
    target_end: int = ASO_TARGET_END,
    dg_threshold: float = DG_FUNCTIONAL_THRESHOLD,
) -> dict[str, Any]:
    """Analisa mutacoes de escape posicao a posicao com DDG detalhado.

    Para cada uma das 75 mutacoes possiveis (25 pos x 3 bases alt):
      - Calcula dG com mismatch (disrupted_dg)
      - Calcula DDG = disrupted_dg - wildtype_dg
      - Classifica como escape se disrupted_dg > threshold

    Retorna analise detalhada incluindo numero minimo de mutacoes
    para escape em cada posicao.

    Args:
        aso_seq: Sequencia do ASO.
        sl_seq: Sequencia do SL RNA completo.
        target_start: Posicao inicial do alvo.
        target_end: Posicao final do alvo.
        dg_threshold: Limiar de dG para escape.

    Returns:
        Dicionario com analise completa de escape.
    """
    target = sl_seq[target_start:target_end].upper()
    aso = aso_seq.upper()
    aso_len = len(aso)
    wt_dg = compute_dg(aso)

    mutations: list[dict[str, Any]] = []
    escape_by_position: dict[int, list[dict[str, Any]]] = {}

    for target_pos in range(len(target)):
        original_base = target[target_pos]
        aso_mm_pos = aso_len - 1 - target_pos
        position_escapes: list[dict[str, Any]] = []

        for alt_base in BASES:
            if alt_base == original_base:
                continue

            disrupted_dg = compute_disrupted_dg(aso, aso_mm_pos)
            ddg = round(disrupted_dg - wt_dg, 4)
            is_escape = disrupted_dg > dg_threshold

            entry = {
                "target_position": target_pos,
                "sl_position": target_start + target_pos,
                "aso_mismatch_position": aso_mm_pos,
                "original": original_base,
                "mutant": alt_base,
                "disrupted_dg_kcal": disrupted_dg,
                "wildtype_dg_kcal": wt_dg,
                "ddg_kcal": ddg,
                "is_escape": is_escape,
            }
            mutations.append(entry)
            if is_escape:
                position_escapes.append(entry)

        escape_by_position[target_pos] = position_escapes

    # Contagem de posicoes que permitem escape com 1 mutacao
    positions_with_escape = sum(
        1 for pos_escapes in escape_by_position.values()
        if len(pos_escapes) > 0
    )

    n_escape = sum(1 for m in mutations if m["is_escape"])

    return {
        "wildtype_dg_kcal": wt_dg,
        "dg_threshold_kcal": dg_threshold,
        "total_mutations": len(mutations),
        "n_escape_single": n_escape,
        "n_positions_with_escape": positions_with_escape,
        "positions_total": len(target),
        "escape_fraction": round(n_escape / len(mutations), 4),
        "mutations": sorted(mutations, key=lambda m: m["disrupted_dg_kcal"], reverse=True),
        "escape_by_position": {
            pos: [
                f"{e['original']}{e['target_position']+1}{e['mutant']} (dG={e['disrupted_dg_kcal']})"
                for e in escapes
            ]
            for pos, escapes in escape_by_position.items()
        },
        "min_mutations_for_escape": 1 if n_escape > 0 else "undetermined (no single-mutation escape)",
    }


# ---------------------------------------------------------------------------
# 3. Kimura fixation probability
# ---------------------------------------------------------------------------


def kimura_fixation_probability(
    selection_coefficient: float,
    effective_population_size: float,
) -> float:
    """Calcula a probabilidade de fixacao pela formula de Kimura (1962).

    P_fix = (1 - exp(-2*s)) / (1 - exp(-2*N_e*s))

    Onde:
      s = coeficiente de selecao (positivo = vantajoso, negativo = deletrio)
      N_e = tamanho efetivo da populacao

    Casos especiais:
      - s = 0 (neutro): P_fix = 1/N_e
      - s < 0 (deletrio): P_fix << 1/N_e
      - s > 0 (vantajoso): P_fix >> 1/N_e

    Para mutacoes de resistencia ao ASO:
      - s < 0 porque a mutacao reduz funcao do SL RNA (custo de fitness)
      - Mesmo que o ASO esteja presente (pressao seletiva contra wildtype),
        o parasita precisa de SL RNA funcional para sobreviver

    Args:
        selection_coefficient: s (negativo = deletrio).
        effective_population_size: N_e.

    Returns:
        Probabilidade de fixacao P_fix.
    """
    s = selection_coefficient
    n_e = effective_population_size

    # Caso neutro: s muito proximo de zero
    if abs(s) < 1e-15:
        return 1.0 / n_e

    # Caso altamente deletrio: 2*N_e*s muito negativo
    # P_fix ~ (1 - exp(-2s)) * exp(2*N_e*s) para N_e*s << -1
    # Numericamente, exp(-2*N_e*s) domina quando N_e*|s| > 300
    two_s = 2.0 * s
    two_ns = 2.0 * n_e * s

    # Protecao contra overflow
    if two_ns < -700:
        # Denominador dominado por exp(-2*N_e*s) que e enorme
        # P_fix ~ (1 - exp(-2s)) / exp(-2*N_e*s) ~ 0
        return 0.0
    if two_ns > 700:
        # Mutacao altamente vantajosa: P_fix ~ 2s
        return min(1.0, two_s)

    numerator = 1.0 - math.exp(-two_s)
    denominator = 1.0 - math.exp(-two_ns)

    if abs(denominator) < 1e-300:
        return 1.0 / n_e

    return max(0.0, min(1.0, numerator / denominator))


def compute_fixation_analysis(
    fitness_values: list[float],
    population_sizes: list[float] | None = None,
    n_copies: int = SL_RNA_COPY_NUMBER,
) -> dict[str, Any]:
    """Analisa probabilidade de fixacao para cada cenario.

    Combina:
      - Fixacao na POPULACAO (Kimura): mutante precisa dominar
      - Fixacao no TANDEM ARRAY: mutacao precisa fixar em ~150 copias

    O coeficiente de selecao s e derivado do custo de fitness:
      s = fitness_mutant - fitness_wildtype = (w - 1)
    Para mutacoes deleterias: w < 1, logo s < 0.

    A probabilidade total de fixacao e:
      P_fix_total = P_fix_population * P_fix_array

    Args:
        fitness_values: Lista de fitness relativa para cada mutacao.
        population_sizes: Tamanhos de populacao a modelar.
        n_copies: Numero de copias no tandem array.

    Returns:
        Dicionario com probabilidades de fixacao por cenario.
    """
    if population_sizes is None:
        population_sizes = POPULATION_SIZES

    # Fitness media dos mutantes
    mean_fitness = float(np.mean(fitness_values)) if fitness_values else 0.0

    # Coeficiente de selecao (negativo = deletrio)
    s = mean_fitness - 1.0  # wildtype tem fitness 1.0

    results: dict[str, Any] = {
        "mean_mutant_fitness": mean_fitness,
        "selection_coefficient_s": round(s, 10),
        "tandem_array_copies": n_copies,
        "fixation_by_population": {},
    }

    for n_pop in population_sizes:
        # Fixacao na populacao (Kimura)
        p_fix_pop = kimura_fixation_probability(s, n_pop)

        # Fixacao no tandem array (deriva neutra como limite superior)
        # Na realidade, gene conversion no array e biased contra mutantes
        # Limite superior conservador: 1/n_copies (drift neutro)
        p_fix_array = 1.0 / n_copies

        # Probabilidade conjunta
        p_fix_total = p_fix_pop * p_fix_array

        pop_key = f"N_{n_pop:.0e}"
        results["fixation_by_population"][pop_key] = {
            "population_size": n_pop,
            "p_fixation_population_kimura": p_fix_pop,
            "p_fixation_array_neutral": p_fix_array,
            "p_fixation_total": p_fix_total,
            "interpretation": (
                f"Kimura P_fix = {p_fix_pop:.2e} na populacao, "
                f"1/{n_copies} = {p_fix_array:.4e} no array. "
                f"Total: {p_fix_total:.2e}"
            ),
        }

    return results


# ---------------------------------------------------------------------------
# 4. Poisson model for resistance emergence
# ---------------------------------------------------------------------------


def poisson_resistance_model(
    effective_mutation_rates: list[dict[str, Any]],
    escape_mutations: dict[str, Any],
    fixation_analysis: dict[str, Any],
    population_sizes: list[float] | None = None,
    gen_time: float = GENERATION_TIME_HOURS,
    treatment_durations_weeks: list[int] | None = None,
) -> dict[str, Any]:
    """Modelo de Poisson para emergencia de resistencia em populacao.

    E[R] = N * mu_eff * T * P_fix
    P(R >= 1) = 1 - exp(-E[R])

    Onde:
      N = tamanho da populacao
      mu_eff = taxa de mutacao efetiva (corrigida por fitness)
      T = duracao do tratamento em geracoes
      P_fix = probabilidade de fixacao (Kimura * array)

    Args:
        effective_mutation_rates: Taxas efetivas por mutacao.
        escape_mutations: Resultado de analyze_escape_mutations().
        fixation_analysis: Resultado de compute_fixation_analysis().
        population_sizes: Tamanhos de populacao.
        gen_time: Tempo de geracao em horas.
        treatment_durations_weeks: Duracoes de tratamento em semanas.

    Returns:
        Dicionario com probabilidades de resistencia por cenario.
    """
    if population_sizes is None:
        population_sizes = POPULATION_SIZES
    if treatment_durations_weeks is None:
        treatment_durations_weeks = TREATMENT_DURATIONS_WEEKS

    # Somar taxa efetiva de TODAS as mutacoes de escape
    escape_mutation_set = {
        (m["target_position"], m["mutant"])
        for m in escape_mutations["mutations"]
        if m["is_escape"]
    }

    total_escape_rate = 0.0
    for rate_entry in effective_mutation_rates:
        key = (rate_entry["position"], rate_entry["mutant_base"])
        if key in escape_mutation_set:
            total_escape_rate += rate_entry["effective_rate_per_generation"]

    results_by_scenario: list[dict[str, Any]] = []

    for n_pop in population_sizes:
        pop_key = f"N_{n_pop:.0e}"
        fix_data = fixation_analysis["fixation_by_population"].get(pop_key, {})
        p_fix_total = fix_data.get("p_fixation_total", 0.0)

        # Taxa combinada: mutacao + escape + fixacao
        lambda_per_generation = total_escape_rate * p_fix_total * n_pop

        for duration_weeks in treatment_durations_weeks:
            # Converter semanas para geracoes
            hours = duration_weeks * 7 * 24
            n_generations = hours / gen_time

            # E[R] = lambda * t
            expected_events = lambda_per_generation * n_generations

            # P(R >= 1) = 1 - exp(-E[R])
            if expected_events > 700:
                p_resistance = 1.0
            elif expected_events < 1e-300:
                p_resistance = 0.0
            else:
                p_resistance = 1.0 - math.exp(-expected_events)

            # Tempo esperado ate primeiro evento: E[T] = 1/lambda
            if lambda_per_generation > 0:
                expected_gen_to_first = 1.0 / lambda_per_generation
                expected_years = (expected_gen_to_first * gen_time) / HOURS_PER_YEAR
            else:
                expected_gen_to_first = math.inf
                expected_years = math.inf

            results_by_scenario.append({
                "population_size": n_pop,
                "treatment_weeks": duration_weeks,
                "treatment_generations": round(n_generations, 1),
                "lambda_per_generation": lambda_per_generation,
                "expected_resistance_events": expected_events,
                "p_resistance_during_treatment": p_resistance,
                "expected_years_to_first_event": (
                    expected_years if not math.isinf(expected_years)
                    else "infinity"
                ),
            })

    return {
        "total_escape_rate_effective": total_escape_rate,
        "n_escape_mutations": len(escape_mutation_set),
        "scenarios": results_by_scenario,
    }


# ---------------------------------------------------------------------------
# 5. Expected time to resistance (combinando tudo)
# ---------------------------------------------------------------------------


def compute_expected_resistance_time(
    total_escape_rate: float,
    p_fixation: float,
    population_size: float,
    gen_time: float = GENERATION_TIME_HOURS,
    p_functional: float = 1.0,
) -> dict[str, Any]:
    """Calcula o tempo esperado para resistencia combinando todos os fatores.

    E[T] = 1 / (N * mu_eff * P_fix * P_functional)

    Onde:
      N = tamanho da populacao
      mu_eff = taxa de mutacao efetiva (ja corrigida por fitness)
      P_fix = probabilidade de fixacao
      P_functional = probabilidade de reter funcao (ja incorporada em mu_eff, mas
                     incluida aqui para analises onde queremos variar separadamente)

    Args:
        total_escape_rate: Taxa total de escape efetiva por parasita por geracao.
        p_fixation: Probabilidade de fixacao (Kimura * array).
        population_size: Tamanho da populacao.
        gen_time: Tempo de geracao em horas.
        p_functional: Probabilidade de reter funcao (default 1.0 se ja em mu_eff).

    Returns:
        Dicionario com tempo esperado em varias unidades.
    """
    # Taxa combinada por geracao na populacao inteira
    combined_rate = population_size * total_escape_rate * p_fixation * p_functional

    if combined_rate <= 0:
        return {
            "combined_rate_per_generation": 0.0,
            "expected_generations": "infinity",
            "expected_hours": "infinity",
            "expected_days": "infinity",
            "expected_years": "infinity",
            "is_finite": False,
        }

    expected_gen = 1.0 / combined_rate
    expected_hours = expected_gen * gen_time
    expected_days = expected_hours / 24.0
    expected_years = expected_hours / HOURS_PER_YEAR

    return {
        "combined_rate_per_generation": combined_rate,
        "expected_generations": expected_gen,
        "expected_hours": expected_hours,
        "expected_days": expected_days,
        "expected_years": expected_years,
        "is_finite": True,
    }


# ---------------------------------------------------------------------------
# 6. Multi-position escape analysis (2+ mutations)
# ---------------------------------------------------------------------------


def analyze_multi_mutation_escape(
    aso_seq: str = ASO_SEQUENCE,
    sl_seq: str = SL_SEQUENCE,
    target_start: int = ASO_TARGET_START,
    target_end: int = ASO_TARGET_END,
    dg_threshold: float = DG_FUNCTIONAL_THRESHOLD,
    effective_rates: list[dict[str, Any]] | None = None,
    max_k: int = 4,
    sample_limit_per_k: int = 500,
) -> dict[str, Any]:
    """Analisa caminhos de escape requerendo k mutacoes simultaneas.

    Para k = 2, 3, 4 mutacoes:
      - Quais combinacoes de posicoes rompem binding?
      - Qual e a taxa combinada? (taxa^k para mutacoes independentes)
      - Qual e o tempo esperado?

    NOTA: Para k >= 2, o numero de combinacoes cresce rapidamente.
    Amostramos ate sample_limit_per_k combinacoes por k.

    A probabilidade de k mutacoes independentes ocorrerem simultaneamente:
      P(k mutacoes) = produto(mu_i para cada mutacao i)
    O tempo esperado cresce EXPONENCIALMENTE com k.

    Args:
        aso_seq: Sequencia do ASO.
        sl_seq: Sequencia do SL RNA.
        target_start: Posicao inicial do alvo.
        target_end: Posicao final do alvo.
        dg_threshold: Limiar de dG para escape.
        effective_rates: Taxas efetivas por mutacao (de fitness.py).
        max_k: Numero maximo de mutacoes simultaneas a analisar.
        sample_limit_per_k: Limite de combinacoes amostradas por k.

    Returns:
        Dicionario com analise de escape multi-mutante.
    """
    target = sl_seq[target_start:target_end].upper()
    aso = aso_seq.upper()
    aso_len = len(aso)
    n_positions = len(target)

    # Indexar taxas efetivas por posicao (media das 3 bases alternativas)
    rate_by_position: dict[int, float] = {}
    if effective_rates:
        for entry in effective_rates:
            pos = entry["position"]
            rate = entry["effective_rate_per_generation"]
            if pos not in rate_by_position:
                rate_by_position[pos] = 0.0
            rate_by_position[pos] += rate  # soma das 3 bases alternativas
    else:
        # Sem taxas efetivas: usar taxa bruta (limite superior)
        from aso_math.config import MUTATION_RATE
        for pos in range(n_positions):
            rate_by_position[pos] = MUTATION_RATE * 3  # 3 bases alternativas

    results_by_k: dict[str, Any] = {}

    for k in range(2, max_k + 1):
        # Todas as combinacoes de k posicoes
        all_combos = list(combinations(range(n_positions), k))
        n_total_combos = len(all_combos)

        # Amostrar se necessario
        if n_total_combos > sample_limit_per_k:
            rng = np.random.default_rng(seed=42)
            sampled_indices = rng.choice(n_total_combos, size=sample_limit_per_k, replace=False)
            sampled_combos = [all_combos[i] for i in sorted(sampled_indices)]
            is_sampled = True
        else:
            sampled_combos = all_combos
            is_sampled = False

        n_escape = 0
        escape_examples: list[dict[str, Any]] = []
        combined_rates: list[float] = []

        for combo in sampled_combos:
            # Converter posicoes do alvo para posicoes no ASO
            aso_mm_positions = [aso_len - 1 - pos for pos in combo]
            multi_dg = compute_multi_mismatch_dg(aso, aso_mm_positions)

            is_escape = multi_dg > dg_threshold

            if is_escape:
                n_escape += 1

                # Taxa combinada: produto das taxas em cada posicao
                combined_rate = 1.0
                for pos in combo:
                    combined_rate *= rate_by_position.get(pos, MUTATION_RATE * 3)

                combined_rates.append(combined_rate)

                if len(escape_examples) < 10:
                    escape_examples.append({
                        "positions": list(combo),
                        "aso_positions": aso_mm_positions,
                        "disrupted_dg_kcal": multi_dg,
                        "combined_rate": combined_rate,
                    })

        # Estatisticas
        escape_fraction = n_escape / max(len(sampled_combos), 1)

        # Extrapolacao: numero estimado de escapes se avaliassemos todos
        estimated_total_escapes = (
            round(escape_fraction * n_total_combos)
            if is_sampled else n_escape
        )

        # Taxa media e maxima dos escapes encontrados
        mean_combined_rate = float(np.mean(combined_rates)) if combined_rates else 0.0
        max_combined_rate = float(np.max(combined_rates)) if combined_rates else 0.0

        results_by_k[f"k={k}"] = {
            "n_mutations": k,
            "total_combinations": n_total_combos,
            "sampled_combinations": len(sampled_combos),
            "is_sampled": is_sampled,
            "n_escape_in_sample": n_escape,
            "escape_fraction": round(escape_fraction, 6),
            "estimated_total_escapes": estimated_total_escapes,
            "mean_combined_rate": mean_combined_rate,
            "max_combined_rate": max_combined_rate,
            "escape_examples": escape_examples,
            "interpretation": (
                f"{k} mutacoes simultaneas: {n_escape}/{len(sampled_combos)} "
                f"combinacoes escapam (taxa media = {mean_combined_rate:.2e}/gen). "
                f"Tempo esperado crece exponencialmente com k."
            ),
        }

    return {
        "max_k_analyzed": max_k,
        "results_by_k": results_by_k,
    }


# ---------------------------------------------------------------------------
# 7. Time comparison table
# ---------------------------------------------------------------------------


def build_resistance_time_table(
    poisson_results: dict[str, Any],
    reference_population: float = 1e8,
) -> dict[str, Any]:
    """Constroi tabela de tempos de resistencia para diferentes cenarios.

    Filtra os resultados do modelo de Poisson para a populacao de referencia
    e organiza por duracao de tratamento.

    Args:
        poisson_results: Resultado de poisson_resistance_model().
        reference_population: Populacao de referencia (default: caso clinico).

    Returns:
        Tabela formatada para o relatorio.
    """
    rows: list[dict[str, Any]] = []

    for scenario in poisson_results["scenarios"]:
        if scenario["population_size"] == reference_population:
            rows.append({
                "treatment_weeks": scenario["treatment_weeks"],
                "treatment_months": round(scenario["treatment_weeks"] / 4.33, 1),
                "p_resistance": scenario["p_resistance_during_treatment"],
                "expected_events": scenario["expected_resistance_events"],
                "expected_years": scenario["expected_years_to_first_event"],
            })

    return {
        "reference_population": reference_population,
        "table": rows,
    }
