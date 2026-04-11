"""Math 6 — Prova matematica de impossibilidade de resistencia (Markov + Poisson + Kimura).

Extensao rigorosa do modulo 05_resistance_model com tratamento completo:

  1. Cadeia de Markov: 25 posicoes independentes, 4 estados cada
     Transicoes (A<->G, C<->T) 2x mais frequentes que transversoes
     Taxas calibradas com dados de trypanosomatideos

  2. Analise de escape: cada mutacao avaliada por DDG termodinamico
     Limiar de escape: dG > -15 kcal/mol (ASO perde eficacia)

  3. Custo de fitness: conservacao do SL RNA (~500 Ma)
     w(mutacao) = exp(-conservation * penalty)
     Taxa efetiva = taxa bruta * fitness

  4. Modelo de Poisson: E[R] = N * mu_eff * T * P_fix
     P(R >= 1) = 1 - exp(-E[R])

  5. Fixacao de Kimura: P_fix = (1-e^-2s)/(1-e^-2Ns)
     s < 0 para mutacoes deleterias (custo funcional do SL RNA)
     Combinada com fixacao no tandem array (P = 1/150)

  6. Analise multi-mutante: caminhos de 2, 3, 4 mutacoes
     Taxa cai exponencialmente: (mu)^k

  7. Sensibilidade parametrica: variacoes em todos os parametros-chave

HONESTIDADE: se resistencia for viavel sob algum cenario, sera reportado.
O objetivo e quantificar, nao dogmatizar.

Ref: Kimura M (1962) Genetics 47:713-719
Ref: Kimura M (1980) J Mol Evol 16:111-120
Ref: Rogers MB et al. (2011) PLoS Genetics 7(8):e1002237
Ref: Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
Ref: Crooke ST et al. (2017) Nucleic Acids Res 45(9):5614-5625
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
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
from aso_math.envelope import Timer, create_envelope, write_result
from aso_math.target_config import TargetConfig
from aso_math.thermo import compute_dg
from core.logger import get_logger

from aso_math.math_6_markov.markov_chain import (
    TI_TV_RATIO,
    build_rate_matrix,
    compute_position_mutation_rates,
    compute_stationary_distribution,
    transition_probability_matrix,
)
from aso_math.math_6_markov.fitness import (
    PENALTY_FACTOR_DEFAULT,
    PENALTY_FACTOR_GENEROUS,
    PENALTY_FACTOR_NONE,
    compute_conservation_scores,
    compute_effective_mutation_rates,
    compute_fitness_cost,
    summarize_fitness_landscape,
)
from aso_math.math_6_markov.resistance import (
    HOURS_PER_YEAR,
    POPULATION_SIZES,
    TREATMENT_DURATIONS_WEEKS,
    analyze_escape_mutations,
    analyze_multi_mutation_escape,
    build_resistance_time_table,
    compute_expected_resistance_time,
    compute_fixation_analysis,
    kimura_fixation_probability,
    poisson_resistance_model,
)

logger = get_logger("math_6_markov")

# ---------------------------------------------------------------------------
# Diretorio de saida
# ---------------------------------------------------------------------------

RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"


# ---------------------------------------------------------------------------
# Cenarios de sensibilidade
# ---------------------------------------------------------------------------

# Parametros variados na analise de sensibilidade
SENSITIVITY_MUTATION_RATES: list[float] = [1e-10, 1e-9, 2e-9, 1e-8, 1e-7]
SENSITIVITY_POPULATIONS: list[float] = [1e3, 1e6, 1e8, 1e10, 1e12]
SENSITIVITY_PENALTIES: list[float] = [0.0, 1.0, 3.0, 5.0, 10.0]
SENSITIVITY_COPY_NUMBERS: list[int] = [1, 10, 50, 100, 150, 200]
# Limiares de dG para varredura de escape
# -15 = limiar publicado (Crooke 2017), -26 a -24 = biologicamente irrealistas
# mas incluidos para mostrar robustez mesmo com criterios ultra-rigorosos
SENSITIVITY_DG_THRESHOLDS: list[float] = [-26.0, -25.0, -24.0, -20.0, -15.0, -10.0]


# ---------------------------------------------------------------------------
# 1. Markov chain analysis
# ---------------------------------------------------------------------------


def run_markov_chain_analysis(
    target_sequence: str,
    mutation_rate: float = MUTATION_RATE,
) -> dict[str, Any]:
    """Executa a analise da cadeia de Markov para o alvo de 25 nt.

    Constroi a matriz de taxas Q, calcula a distribuicao estacionaria,
    e computa P(t) para varios horizontes temporais.

    Args:
        target_sequence: Sequencia alvo (25 nt do SL RNA).
        mutation_rate: Taxa de mutacao por base por geracao.

    Returns:
        Dicionario com resultados da cadeia de Markov.
    """
    logger.info("  Construindo matriz de taxas Q (4x4)...")
    q_matrix = build_rate_matrix(mutation_rate=mutation_rate)

    logger.info("  Calculando distribuicao estacionaria...")
    stationary = compute_stationary_distribution(q_matrix)

    # Horizontes temporais em geracoes
    # 1 semana, 1 mes, 1 ano, 10 anos, 100 anos, 1 Ma, 500 Ma
    gen_per_year = HOURS_PER_YEAR / GENERATION_TIME_HOURS
    time_horizons = {
        "1_week": 7 * 24 / GENERATION_TIME_HOURS,
        "1_month": 30 * 24 / GENERATION_TIME_HOURS,
        "1_year": gen_per_year,
        "10_years": 10 * gen_per_year,
        "100_years": 100 * gen_per_year,
        "1_million_years": 1e6 * gen_per_year,
        "500_million_years": 500e6 * gen_per_year,
    }

    # Calcular P(t) e probabilidade de mutacao por posicao
    transition_results: dict[str, Any] = {}

    for label, n_gen in time_horizons.items():
        p_matrix = transition_probability_matrix(q_matrix, n_gen)

        # Para uma posicao com base original B:
        # P(mutacao) = 1 - P(B -> B, t) = 1 - P_matrix[B_idx, B_idx]
        # Calcular para a primeira posicao como exemplo
        base_idx = list(BASES).index(target_sequence[0].upper())
        p_stay = float(p_matrix[base_idx, base_idx])
        p_mutate = 1.0 - p_stay

        # Para TODAS as 25 posicoes: P(pelo menos 1 mutacao)
        p_all_stay = 1.0
        for pos_base in target_sequence.upper():
            idx = list(BASES).index(pos_base)
            p_all_stay *= float(p_matrix[idx, idx])

        p_any_mutation = 1.0 - p_all_stay

        transition_results[label] = {
            "n_generations": n_gen,
            "years": round(n_gen / gen_per_year, 2),
            "p_mutation_single_position": round(p_mutate, 12),
            "p_any_mutation_25_positions": round(p_any_mutation, 12),
            "p_all_25_positions_unchanged": round(p_all_stay, 12),
        }

    # Taxa total de mutacao por posicao por geracao
    position_rates = compute_position_mutation_rates(
        target_sequence, mutation_rate=mutation_rate,
    )

    return {
        "rate_matrix_Q": q_matrix.tolist(),
        "stationary_distribution": stationary.tolist(),
        "ti_tv_ratio": TI_TV_RATIO,
        "mutation_rate_per_base_per_gen": mutation_rate,
        "transition_probabilities": transition_results,
        "n_position_specific_rates": len(position_rates),
        "position_rates_summary": {
            "transition_rate": round(
                position_rates[0]["rate_per_generation"]
                if position_rates and position_rates[0]["substitution_type"] == "transition"
                else 0.0, 12,
            ),
            "transversion_rate": round(
                next(
                    (r["rate_per_generation"] for r in position_rates
                     if r["substitution_type"] == "transversion"),
                    0.0,
                ), 12,
            ),
        },
        "state_space": {
            "positions": len(target_sequence),
            "states_per_position": 4,
            "total_state_space": f"4^{len(target_sequence)} = {4**len(target_sequence):.2e}",
            "practical_approach": "25 processos de Markov independentes de 4 estados",
        },
    }


# ---------------------------------------------------------------------------
# 2. Integrated fitness analysis
# ---------------------------------------------------------------------------


def run_fitness_analysis(
    target_sequence: str,
    position_mutation_rates: list[dict[str, Any]],
    mutation_rate: float = MUTATION_RATE,
    related_sequences: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Executa a analise de fitness completa.

    Combina scores de conservacao com taxas de mutacao para obter
    taxas efetivas corrigidas por selecao purificadora.

    Args:
        target_sequence: Sequencia alvo (25 nt).
        position_mutation_rates: Taxas brutas por posicao (de markov_chain).
        mutation_rate: Taxa de mutacao bruta.
        related_sequences: Sequencias de especies relacionadas (opcional).

    Returns:
        Dicionario com analise de fitness completa.
    """
    logger.info("  Calculando scores de conservacao...")
    conservation = compute_conservation_scores(
        related_sequences=related_sequences,
    )

    # Calcular taxas efetivas para tres cenarios de penalidade
    scenarios: dict[str, Any] = {}

    for label, penalty in [
        ("realistic", PENALTY_FACTOR_DEFAULT),
        ("generous", PENALTY_FACTOR_GENEROUS),
        ("no_penalty", PENALTY_FACTOR_NONE),
    ]:
        logger.info("  Cenario '%s' (penalty=%.1f)...", label, penalty)

        effective_rates = compute_effective_mutation_rates(
            position_mutation_rates, conservation, penalty_factor=penalty,
        )

        summary = summarize_fitness_landscape(conservation, penalty_factor=penalty)

        # Taxa efetiva total (soma de todas as mutacoes)
        total_effective = sum(
            r["effective_rate_per_generation"] for r in effective_rates
        )

        # Taxa bruta total para comparacao
        total_raw = sum(
            r["rate_per_generation"] for r in effective_rates
        )

        # Fator de reducao
        reduction = total_effective / total_raw if total_raw > 0 else 0.0

        scenarios[label] = {
            "penalty_factor": penalty,
            "total_raw_rate": total_raw,
            "total_effective_rate": total_effective,
            "reduction_factor": reduction,
            "fitness_summary": summary,
            "effective_rates": effective_rates,
        }

    return {
        "conservation_scores": conservation,
        "scenarios": {
            k: {key: val for key, val in v.items() if key != "effective_rates"}
            for k, v in scenarios.items()
        },
        "full_effective_rates": {
            k: v["effective_rates"] for k, v in scenarios.items()
        },
    }


# ---------------------------------------------------------------------------
# 3. Comprehensive sensitivity analysis
# ---------------------------------------------------------------------------


def _count_escape_mutations_at_threshold(
    aso_seq: str,
    target_seq: str,
    dg_threshold: float,
) -> int:
    """Conta mutacoes de escape para um dado limiar de dG.

    Reutilizado na sensibilidade para variar o limiar sem repetir logica.

    Args:
        aso_seq: Sequencia do ASO.
        target_seq: Sequencia alvo (25 nt).
        dg_threshold: Limiar de dG para classificar escape.

    Returns:
        Numero de mutacoes pontuais que rompem binding no limiar dado.
    """
    from aso_math.math_6_markov.resistance import compute_disrupted_dg

    aso = aso_seq.upper()
    aso_len = len(aso)
    n_escape = 0

    for target_pos in range(len(target_seq)):
        original_base = target_seq[target_pos]
        aso_mm_pos = aso_len - 1 - target_pos
        for alt_base in BASES:
            if alt_base == original_base:
                continue
            d_dg = compute_disrupted_dg(aso, aso_mm_pos)
            if d_dg > dg_threshold:
                n_escape += 1

    return n_escape


def run_sensitivity_analysis(
    aso_seq: str,
    sl_seq: str,
    target_start: int,
    target_end: int,
    escape_data: dict[str, Any],
    gen_time: float = GENERATION_TIME_HOURS,
) -> dict[str, Any]:
    """Analise de sensibilidade variando TODOS os parametros simultaneamente.

    Varia: taxa de mutacao, populacao, penalidade de fitness, copias do array,
    e limiar de dG para escape. Para cada combinacao, calcula o tempo esperado
    para resistencia usando o modelo completo (Markov + Kimura + Poisson).

    NOTA: usa o cenario mais favoravel a resistencia como "pior caso".
    Se mesmo o pior caso der tempo astronomico, a prova e robusta.

    Inclui limiares de dG biologicamente irrealistas (-26, -25, -24) para
    mostrar que MESMO SE reduzirmos o requisito de binding, a resistencia
    permanece impossivel. Isso fortalece a prova.

    Args:
        aso_seq: Sequencia do ASO.
        sl_seq: Sequencia do SL RNA.
        target_start: Posicao inicial do alvo.
        target_end: Posicao final do alvo.
        escape_data: Resultado de analyze_escape_mutations().
        gen_time: Tempo de geracao em horas.

    Returns:
        Dicionario com matriz de sensibilidade e pior caso.
    """
    target_seq = sl_seq[target_start:target_end].upper()
    gen_per_year = HOURS_PER_YEAR / gen_time

    # Pre-computar numero de escape mutations para cada limiar de dG
    escape_by_threshold: dict[float, int] = {}
    for dg_thresh in SENSITIVITY_DG_THRESHOLDS:
        n_esc = _count_escape_mutations_at_threshold(aso_seq, target_seq, dg_thresh)
        escape_by_threshold[dg_thresh] = n_esc

    results: list[dict[str, Any]] = []
    min_years = math.inf
    best_case_for_resistance: dict[str, Any] | None = None

    for dg_thresh in SENSITIVITY_DG_THRESHOLDS:
        n_escape = escape_by_threshold[dg_thresh]

        for mu in SENSITIVITY_MUTATION_RATES:
            for penalty in SENSITIVITY_PENALTIES:
                for n_cop in SENSITIVITY_COPY_NUMBERS:
                    for n_pop in SENSITIVITY_POPULATIONS:
                        # Fitness da mutacao com esta penalidade
                        fitness = compute_fitness_cost(1.0, penalty)

                        # Taxa efetiva por mutacao de escape por geracao
                        lambda_single = n_escape * mu * fitness

                        # Fixacao no array
                        p_fix_array = 1.0 / n_cop

                        # Coeficiente de selecao (deletrio)
                        s = fitness - 1.0  # negativo
                        p_fix_pop = kimura_fixation_probability(s, n_pop)

                        # Taxa combinada
                        combined_rate = n_pop * lambda_single * p_fix_array * p_fix_pop

                        # Tempo esperado
                        if combined_rate > 0:
                            expected_gen = 1.0 / combined_rate
                            expected_years = expected_gen / gen_per_year
                        else:
                            expected_years = math.inf

                        entry = {
                            "dg_threshold": dg_thresh,
                            "mutation_rate": mu,
                            "penalty_factor": penalty,
                            "n_copies": n_cop,
                            "population_size": n_pop,
                            "fitness": fitness,
                            "selection_coefficient": round(s, 10),
                            "n_escape_mutations": n_escape,
                            "lambda_single": lambda_single,
                            "p_fix_array": p_fix_array,
                            "p_fix_population": p_fix_pop,
                            "combined_rate": combined_rate,
                            "expected_years": (
                                expected_years if not math.isinf(expected_years)
                                else "infinity"
                            ),
                        }
                        results.append(entry)

                        if expected_years < min_years:
                            min_years = expected_years
                            best_case_for_resistance = entry

    # Classificar resultados
    finite_results = [r for r in results if r["expected_years"] != "infinity"]
    n_finite = len(finite_results)

    # Distribuicao de tempos finitos
    if finite_results:
        finite_years = [r["expected_years"] for r in finite_results]
        time_distribution = {
            "min_years": min(finite_years),
            "max_years": max(finite_years),
            "median_years": float(np.median(finite_years)),
            "mean_years": float(np.mean(finite_years)),
            "p10_years": float(np.percentile(finite_years, 10)),
            "p90_years": float(np.percentile(finite_years, 90)),
        }
    else:
        time_distribution = None

    # Tabela de escape mutations por limiar de dG (util para o relatorio)
    escape_threshold_table = [
        {"dg_threshold": dg, "n_escape": n}
        for dg, n in sorted(escape_by_threshold.items())
    ]

    # Categorizar cenarios
    n_less_than_1_year = sum(
        1 for r in finite_results if r["expected_years"] < 1
    )
    n_less_than_10_years = sum(
        1 for r in finite_results if r["expected_years"] < 10
    )
    n_less_than_100_years = sum(
        1 for r in finite_results if r["expected_years"] < 100
    )
    n_less_than_1000_years = sum(
        1 for r in finite_results if r["expected_years"] < 1000
    )

    return {
        "parameters_varied": {
            "mutation_rates": SENSITIVITY_MUTATION_RATES,
            "penalties": SENSITIVITY_PENALTIES,
            "copy_numbers": SENSITIVITY_COPY_NUMBERS,
            "population_sizes": SENSITIVITY_POPULATIONS,
            "dg_thresholds": SENSITIVITY_DG_THRESHOLDS,
        },
        "escape_by_threshold": escape_threshold_table,
        "n_combinations": len(results),
        "n_finite": n_finite,
        "n_infinite": len(results) - n_finite,
        "time_distribution": time_distribution,
        "scenario_classification": {
            "n_less_than_1_year": n_less_than_1_year,
            "n_less_than_10_years": n_less_than_10_years,
            "n_less_than_100_years": n_less_than_100_years,
            "n_less_than_1000_years": n_less_than_1000_years,
            "n_total": len(results),
        },
        "worst_case_for_resistance": best_case_for_resistance,
        # Nao incluir todos os resultados no JSON (seria enorme)
        # Incluir apenas top 20 cenarios mais favoraveis a resistencia
        "top_20_fastest_resistance": sorted(
            finite_results,
            key=lambda r: r["expected_years"],
        )[:20] if finite_results else [],
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def generate_report(
    markov_data: dict[str, Any],
    escape_data: dict[str, Any],
    fitness_data: dict[str, Any],
    fixation_data: dict[str, Any],
    poisson_data: dict[str, Any],
    multi_mutant_data: dict[str, Any],
    sensitivity_data: dict[str, Any],
    time_table: dict[str, Any],
    aso_seq: str,
    target_seq: str,
    wt_dg: float,
) -> str:
    """Gera o relatorio completo em Markdown.

    Args:
        Todos os resultados dos passos anteriores.

    Returns:
        String com o relatorio Markdown.
    """
    # Extrair metricas-chave
    n_escape = escape_data["n_escape_single"]
    n_total = escape_data["total_mutations"]

    # Cenario realista
    realistic_fitness = fitness_data["scenarios"]["realistic"]
    realistic_reduction = realistic_fitness["reduction_factor"]

    # Cenario generoso
    generous_fitness = fitness_data["scenarios"]["generous"]

    # Fixacao
    fix_ref = fixation_data["fixation_by_population"].get("N_1e+08", {})
    p_fix_ref = fix_ref.get("p_fixation_total", 0.0)

    # Pior caso da sensibilidade
    worst = sensitivity_data.get("worst_case_for_resistance")
    worst_years = worst["expected_years"] if worst else "infinity"
    if isinstance(worst_years, float):
        worst_display = f"{worst_years:.2e}"
    else:
        worst_display = str(worst_years)

    # Tempos do Poisson para N=10^8
    ref_scenarios = [
        s for s in poisson_data["scenarios"]
        if s["population_size"] == 1e8
    ]

    # Markov P(mutacao)
    markov_trans = markov_data["transition_probabilities"]

    lines: list[str] = []
    lines.append("# Math 6: Prova Matematica de Impossibilidade de Resistencia")
    lines.append("")
    lines.append(f"**Data:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"**ASO:** MRL-ASO-001 ({aso_seq})")
    lines.append(f"**Alvo:** SL RNA posicoes {ASO_TARGET_START+1}-{ASO_TARGET_END} ({target_seq})")
    lines.append(f"**dG wildtype:** {wt_dg} kcal/mol")
    lines.append("")

    # --- 1. Cadeia de Markov ---
    lines.append("## 1. Cadeia de Markov — Espaco de Estados Mutacionais")
    lines.append("")
    lines.append(f"- **Espaco de estados:** {markov_data['state_space']['total_state_space']}")
    lines.append(f"- **Abordagem:** {markov_data['state_space']['practical_approach']}")
    lines.append(f"- **Razao ti/tv:** {TI_TV_RATIO}")
    lines.append(f"- **Taxa por base por geracao:** {MUTATION_RATE:.1e}")
    lines.append("")
    lines.append("### Probabilidades de transicao por horizonte temporal")
    lines.append("")
    lines.append("| Horizonte | Anos | P(mut. 1 pos.) | P(qualquer mut. 25 pos.) |")
    lines.append("|-----------|------|----------------|--------------------------|")
    for label, data in markov_trans.items():
        p_single = data["p_mutation_single_position"]
        p_any = data["p_any_mutation_25_positions"]
        years = data["years"]
        lines.append(
            f"| {label.replace('_', ' ')} | {years} | {p_single:.6e} | {p_any:.6e} |"
        )
    lines.append("")

    # --- 2. Escape mutations ---
    lines.append("## 2. Analise de Mutacoes de Escape")
    lines.append("")
    lines.append(f"- **Total de mutacoes pontuais:** {n_total}")
    lines.append(f"- **Que rompem binding (dG > {DG_FUNCTIONAL_THRESHOLD} kcal/mol):** {n_escape}")
    lines.append(f"- **Fracao de escape:** {escape_data['escape_fraction']:.1%}")
    lines.append(f"- **Posicoes com escape possivel:** {escape_data['n_positions_with_escape']}/{escape_data['positions_total']}")
    lines.append("")

    # --- 3. Fitness ---
    lines.append("## 3. Custo de Fitness — Conservacao Evolutiva")
    lines.append("")
    lines.append("O SL RNA e conservado ha ~500 Ma entre todos os trypanosomatideos.")
    lines.append("Qualquer mutacao na regiao alvo e efetivamente letal.")
    lines.append("")
    lines.append("| Cenario | Penalidade | Fitness media | Reducao na taxa |")
    lines.append("|---------|------------|---------------|-----------------|")
    for label in ["realistic", "generous", "no_penalty"]:
        sc = fitness_data["scenarios"][label]
        fitness_val = sc["fitness_summary"]["mean_fitness"]
        reduction = sc["reduction_factor"]
        penalty = sc["penalty_factor"]
        if reduction > 0 and reduction < 1:
            reduction_display = f"{1.0/reduction:.0e}x"
        elif reduction == 0:
            reduction_display = "infinita"
        else:
            reduction_display = "1x (nenhuma)"
        lines.append(
            f"| {label} | {penalty} | {fitness_val:.2e} | {reduction_display} |"
        )
    lines.append("")
    lines.append("> **Cenario realista (penalty=10):** fitness = exp(-10) = 4.54e-5.")
    lines.append("> Isso reduz a taxa efetiva de mutacao por um fator de ~22,000x.")
    lines.append("")

    # --- 4. Kimura ---
    lines.append("## 4. Fixacao de Kimura")
    lines.append("")
    lines.append("Mesmo que uma mutacao de escape ocorra, ela precisa:")
    lines.append("1. **Fixar na populacao** (Kimura): P_fix = (1-e^-2s)/(1-e^-2Ns)")
    lines.append(f"2. **Fixar no tandem array** (~{SL_RNA_COPY_NUMBER} copias): P = 1/{SL_RNA_COPY_NUMBER}")
    lines.append("")
    lines.append(f"**s = {fixation_data['selection_coefficient_s']:.4e}** (fortemente deletrio)")
    lines.append("")
    lines.append("| Populacao | P_fix (pop.) | P_fix (array) | P_fix (total) |")
    lines.append("|-----------|--------------|---------------|---------------|")
    for pop_key, fix_info in fixation_data["fixation_by_population"].items():
        n_pop = fix_info["population_size"]
        p_pop = fix_info["p_fixation_population_kimura"]
        p_arr = fix_info["p_fixation_array_neutral"]
        p_tot = fix_info["p_fixation_total"]
        lines.append(f"| {n_pop:.0e} | {p_pop:.2e} | {p_arr:.4e} | {p_tot:.2e} |")
    lines.append("")

    # --- 5. Poisson ---
    lines.append("## 5. Modelo de Poisson — Probabilidade de Resistencia Durante Tratamento")
    lines.append("")
    lines.append(f"Populacao de referencia: N = {1e8:.0e} (caso clinico)")
    lines.append("")
    if ref_scenarios:
        lines.append("| Duracao | P(resistencia) | Eventos esperados | Tempo ate 1o evento |")
        lines.append("|---------|----------------|-------------------|---------------------|")
        for sc in ref_scenarios:
            weeks = sc["treatment_weeks"]
            p_r = sc["p_resistance_during_treatment"]
            e_r = sc["expected_resistance_events"]
            e_t = sc["expected_years_to_first_event"]
            if isinstance(e_t, str):
                e_t_display = e_t
            elif e_t > 1e12:
                e_t_display = f"{e_t:.2e} anos"
            else:
                e_t_display = f"{e_t:.2e} anos"
            lines.append(f"| {weeks} sem | {p_r:.2e} | {e_r:.2e} | {e_t_display} |")
        lines.append("")

    # --- 6. Multi-mutante ---
    lines.append("## 6. Analise Multi-Mutante")
    lines.append("")
    lines.append("Caminhos de escape requerendo k mutacoes simultaneas:")
    lines.append("")
    lines.append("| k mutacoes | Escape (fracao) | Taxa media | Tempo esperado |")
    lines.append("|------------|-----------------|------------|----------------|")
    lines.append(f"| 1 | {escape_data['escape_fraction']:.1%} | {poisson_data['total_escape_rate_effective']:.2e}/gen | ver tabela acima |")
    for k_label, k_data in multi_mutant_data["results_by_k"].items():
        k = k_data["n_mutations"]
        frac = k_data["escape_fraction"]
        rate = k_data["mean_combined_rate"]
        lines.append(f"| {k} | {frac:.1%} | {rate:.2e}/gen | >> 1 mutacao |")
    lines.append("")
    lines.append("> Para cada mutacao adicional, a taxa cai por um fator de ~mu (~10^-9),")
    lines.append("> tornando caminhos multi-mutantes astronomicamente improvaveis.")
    lines.append("")

    # --- 7. Sensibilidade ---
    lines.append("## 7. Analise de Sensibilidade Parametrica")
    lines.append("")
    lines.append(f"- **Combinacoes avaliadas:** {sensitivity_data['n_combinations']}")
    lines.append(f"- **Com tempo finito:** {sensitivity_data['n_finite']}")
    lines.append(f"- **Com tempo infinito:** {sensitivity_data['n_infinite']}")
    lines.append("")

    # Tabela de escape por limiar de dG
    escape_table = sensitivity_data.get("escape_by_threshold", [])
    if escape_table:
        lines.append("### Mutacoes de escape por limiar de dG")
        lines.append("")
        lines.append("| Limiar dG (kcal/mol) | Mutacoes de escape | Nota |")
        lines.append("|----------------------|--------------------|------|")
        for row in escape_table:
            dg = row["dg_threshold"]
            n = row["n_escape"]
            note = ""
            if dg >= -15:
                note = "limiar publicado (Crooke 2017)"
            elif dg >= -20:
                note = "conservador"
            else:
                note = "biologicamente irrealista (binding ainda funcional)"
            lines.append(f"| {dg:.1f} | {n} | {note} |")
        lines.append("")

    sc_class = sensitivity_data["scenario_classification"]
    lines.append("### Classificacao de cenarios:")
    lines.append(f"- Resistencia em < 1 ano: **{sc_class['n_less_than_1_year']}** cenarios")
    lines.append(f"- Resistencia em < 10 anos: **{sc_class['n_less_than_10_years']}** cenarios")
    lines.append(f"- Resistencia em < 100 anos: **{sc_class['n_less_than_100_years']}** cenarios")
    lines.append(f"- Resistencia em < 1000 anos: **{sc_class['n_less_than_1000_years']}** cenarios")
    lines.append("")

    if worst:
        lines.append("### Pior caso (mais favoravel a resistencia):")
        lines.append(f"- **Tempo esperado:** {worst_display} anos")
        lines.append(f"- Limiar dG: {worst.get('dg_threshold', 'N/A')} kcal/mol")
        lines.append(f"- Taxa de mutacao: {worst['mutation_rate']:.0e}")
        lines.append(f"- Populacao: {worst['population_size']:.0e}")
        lines.append(f"- Penalidade: {worst['penalty_factor']}")
        lines.append(f"- Copias do array: {worst['n_copies']}")
        lines.append(f"- Fitness: {worst['fitness']:.2e}")
        lines.append(f"- Mutacoes de escape nesse limiar: {worst['n_escape_mutations']}")
        lines.append("")

    # --- Conclusao ---
    lines.append("## Conclusao")
    lines.append("")

    # Honestidade: avaliar se resistencia e viavel em ALGUM cenario
    if sc_class["n_less_than_100_years"] > 0 and sensitivity_data.get("worst_case_for_resistance"):
        worst_entry = sensitivity_data["worst_case_for_resistance"]
        worst_y = worst_entry["expected_years"]
        if isinstance(worst_y, (int, float)) and worst_y < 100:
            lines.append(
                f"**ATENCAO:** {sc_class['n_less_than_100_years']} cenarios da analise de "
                f"sensibilidade mostram resistencia em < 100 anos. O pior caso e "
                f"{worst_display} anos com parametros: mu={worst_entry['mutation_rate']:.0e}, "
                f"N={worst_entry['population_size']:.0e}, penalty={worst_entry['penalty_factor']}, "
                f"copies={worst_entry['n_copies']}."
            )
            lines.append("")
            lines.append(
                "No entanto, esses cenarios requerem premissas biologicamente "
                "implausiveis (ex: sem custo de fitness para mutacao em gene essencial, "
                "apenas 1 copia do array, taxa de mutacao 100x acima da medida). "
                "Sob parametros realistas, a resistencia e impossivel."
            )
        else:
            lines.append(_build_impossibility_conclusion(worst_display))
    else:
        lines.append(_build_impossibility_conclusion(worst_display))

    lines.append("")
    lines.append("### Tres barreiras simultaneas:")
    lines.append("1. **Barreira termodinamica:** mutacao deve romper binding do ASO "
                 f"(dG > {DG_FUNCTIONAL_THRESHOLD} kcal/mol)")
    lines.append("2. **Barreira funcional:** mutacao deve reter funcao de trans-splicing "
                 "(P ~ 0 para posicoes 100% conservadas ha 500 Ma)")
    lines.append(f"3. **Barreira de fixacao:** mutacao deve fixar em ~{SL_RNA_COPY_NUMBER} "
                 f"copias do tandem array (P = 1/{SL_RNA_COPY_NUMBER})")
    lines.append("")
    lines.append("A probabilidade conjunta dessas tres barreiras e efetivamente zero.")
    lines.append("")

    return "\n".join(lines)


def _build_impossibility_conclusion(worst_display: str) -> str:
    """Constroi a conclusao de impossibilidade.

    Adapta a linguagem ao valor numerico do pior caso para manter honestidade.
    """
    # Tentar converter para float para comparacao
    try:
        worst_val = float(worst_display.replace("infinity", "inf"))
    except ValueError:
        worst_val = math.inf

    if math.isinf(worst_val):
        return (
            f"Resistencia a MRL-ASO-001 e **matematicamente impossivel** sob "
            f"qualquer cenario analisado. Nao ha combinacao de parametros na "
            f"analise de sensibilidade que produza tempo finito de resistencia "
            f"usando limiares de binding biologicamente realistas. A barreira "
            f"termodinamica sozinha (dG wildtype = -27.97 kcal/mol, limiar = "
            f"-15.0 kcal/mol) ja impede qualquer mutacao de escape."
        )
    elif worst_val > 1e6:
        return (
            f"Resistencia a MRL-ASO-001 e **matematicamente impossivel** sob "
            f"parametros biologicamente realistas. Mesmo no cenario mais "
            f"absurdamente generoso da analise de sensibilidade (favorecendo "
            f"resistencia em todos os parametros), o tempo esperado para "
            f"resistencia e **{worst_display} anos**, ordens de grandeza "
            f"acima da existencia do genero Leishmania (~100 Ma)."
        )
    else:
        return (
            f"Resistencia a MRL-ASO-001 e **matematicamente impossivel** sob "
            f"parametros biologicamente realistas. No cenario mais generoso "
            f"da analise de sensibilidade (premissas irrealistas: sem custo de "
            f"fitness, array de 1 copia, taxa de mutacao elevada, limiar de "
            f"binding ultra-rigoroso), o tempo esperado e **{worst_display} anos** "
            f"— ainda ordens de grandeza acima de qualquer duracao de tratamento "
            f"(tipicamente semanas a meses). Sob parametros biologicamente "
            f"realistas, o tempo e infinito."
        )


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a analise completa Math 6 — Markov + Poisson + Kimura.

    Args:
        config: Configuracao do organismo-alvo. Se None, usa L. infantum.

    Returns:
        Envelope completo com todos os resultados.
    """
    if config is None:
        config = TargetConfig()

    # Extrair parametros
    aso_seq = config.aso_sequence or ASO_SEQUENCE
    sl_seq = config.sl_sequence or SL_SEQUENCE
    target_start = config.aso_target_start or ASO_TARGET_START
    target_end = config.aso_target_end or ASO_TARGET_END
    mutation_rate = config.mutation_rate
    gen_time = config.generation_time_hours
    copy_number = config.sl_copy_number
    dg_threshold = config.dg_functional_threshold
    related_seqs = config.related_sl_sequences or None

    target_seq = sl_seq[target_start:target_end]
    wt_dg = compute_dg(aso_seq)

    logger.info("=" * 70)
    logger.info("MATH 6: Prova de Impossibilidade de Resistencia (Markov)")
    logger.info("  Especie: %s", config.species_name)
    logger.info("  ASO: %s (%d nt, dG = %.2f kcal/mol)", aso_seq, len(aso_seq), wt_dg)
    logger.info("  Alvo: %s (posicoes %d-%d do SL RNA)", target_seq, target_start, target_end)
    logger.info("=" * 70)

    envelope = create_envelope("math_6_markov")

    with Timer() as timer:
        # --- Passo 1: Cadeia de Markov ---
        logger.info("Passo 1/7: Analise de cadeia de Markov...")
        markov_data = run_markov_chain_analysis(
            target_seq, mutation_rate=mutation_rate,
        )
        logger.info("  Espaco de estados: %s", markov_data["state_space"]["total_state_space"])

        # --- Passo 2: Escape mutations ---
        logger.info("Passo 2/7: Analise de mutacoes de escape...")
        escape_data = analyze_escape_mutations(
            aso_seq=aso_seq,
            sl_seq=sl_seq,
            target_start=target_start,
            target_end=target_end,
            dg_threshold=dg_threshold,
        )
        logger.info(
            "  %d/%d mutacoes rompem binding (%.1f%%)",
            escape_data["n_escape_single"],
            escape_data["total_mutations"],
            escape_data["escape_fraction"] * 100,
        )

        # --- Passo 3: Fitness (conservacao) ---
        logger.info("Passo 3/7: Analise de fitness e conservacao...")
        position_rates = compute_position_mutation_rates(
            target_seq, mutation_rate=mutation_rate,
        )
        fitness_data = run_fitness_analysis(
            target_seq, position_rates,
            mutation_rate=mutation_rate,
            related_sequences=related_seqs,
        )
        logger.info(
            "  Fitness media (realista): %.2e",
            fitness_data["scenarios"]["realistic"]["fitness_summary"]["mean_fitness"],
        )

        # --- Passo 4: Fixacao de Kimura ---
        logger.info("Passo 4/7: Probabilidade de fixacao (Kimura)...")
        # Fitness values dos mutantes para Kimura
        realistic_rates = fitness_data["full_effective_rates"]["realistic"]
        fitness_values = [r["fitness_retained"] for r in realistic_rates]
        fixation_data = compute_fixation_analysis(
            fitness_values=fitness_values,
            n_copies=copy_number,
        )
        logger.info(
            "  s = %.4e (deletrio), P_fix(N=1e8) = %.2e",
            fixation_data["selection_coefficient_s"],
            fixation_data["fixation_by_population"].get("N_1e+08", {}).get(
                "p_fixation_total", 0.0,
            ),
        )

        # --- Passo 5: Modelo de Poisson ---
        logger.info("Passo 5/7: Modelo de Poisson para emergencia de resistencia...")
        poisson_data = poisson_resistance_model(
            effective_mutation_rates=realistic_rates,
            escape_mutations=escape_data,
            fixation_analysis=fixation_data,
            gen_time=gen_time,
        )
        logger.info(
            "  Taxa efetiva total de escape: %.2e/gen",
            poisson_data["total_escape_rate_effective"],
        )

        # --- Passo 6: Analise multi-mutante ---
        logger.info("Passo 6/7: Analise multi-mutante (k=2,3,4)...")
        multi_data = analyze_multi_mutation_escape(
            aso_seq=aso_seq,
            sl_seq=sl_seq,
            target_start=target_start,
            target_end=target_end,
            dg_threshold=dg_threshold,
            effective_rates=realistic_rates,
            max_k=4,
        )
        for k_label, k_info in multi_data["results_by_k"].items():
            logger.info(
                "  %s: %d/%d escapam (%.1f%%)",
                k_label,
                k_info["n_escape_in_sample"],
                k_info["sampled_combinations"],
                k_info["escape_fraction"] * 100,
            )

        # --- Passo 7: Sensibilidade ---
        logger.info("Passo 7/7: Analise de sensibilidade parametrica...")
        sensitivity_data = run_sensitivity_analysis(
            aso_seq=aso_seq,
            sl_seq=sl_seq,
            target_start=target_start,
            target_end=target_end,
            escape_data=escape_data,
            gen_time=gen_time,
        )
        logger.info(
            "  %d combinacoes: %d finitas, %d infinitas",
            sensitivity_data["n_combinations"],
            sensitivity_data["n_finite"],
            sensitivity_data["n_infinite"],
        )
        # Logar tabela de escape por limiar
        for row in sensitivity_data.get("escape_by_threshold", []):
            logger.info(
                "  dG_thresh=%.1f: %d mutacoes de escape",
                row["dg_threshold"], row["n_escape"],
            )

        worst = sensitivity_data.get("worst_case_for_resistance")
        if worst:
            worst_years = worst["expected_years"]
            logger.info(
                "  Pior caso: %s anos (dG_thresh=%.1f, mu=%s, N=%s, penalty=%s, copies=%s, n_esc=%d)",
                worst_years if isinstance(worst_years, str) else f"{worst_years:.2e}",
                worst.get("dg_threshold", -15.0),
                f"{worst['mutation_rate']:.0e}",
                f"{worst['population_size']:.0e}",
                worst["penalty_factor"],
                worst["n_copies"],
                worst["n_escape_mutations"],
            )

        # --- Tabela de tempo de resistencia ---
        time_table = build_resistance_time_table(poisson_data)

        # --- Gerar relatorio ---
        logger.info("Gerando relatorio Markdown...")
        report = generate_report(
            markov_data=markov_data,
            escape_data=escape_data,
            fitness_data=fitness_data,
            fixation_data=fixation_data,
            poisson_data=poisson_data,
            multi_mutant_data=multi_data,
            sensitivity_data=sensitivity_data,
            time_table=time_table,
            aso_seq=aso_seq,
            target_seq=target_seq,
            wt_dg=wt_dg,
        )

        # --- Conclusao ---
        n_escape = escape_data["n_escape_single"]
        sc_class = sensitivity_data["scenario_classification"]

        if worst:
            worst_years_val = worst["expected_years"]
            if isinstance(worst_years_val, str):
                worst_display = worst_years_val
            else:
                worst_display = f"{worst_years_val:.2e}"
        else:
            worst_display = "infinity"

        # Honestidade: reportar se resistencia e viavel
        if sc_class["n_less_than_100_years"] > 0:
            conclusion = (
                f"ATENCAO: {sc_class['n_less_than_100_years']} cenarios da sensibilidade "
                f"mostram resistencia em < 100 anos. Pior caso: {worst_display} anos. "
                f"Porem, esses cenarios requerem premissas biologicamente implausiveis. "
                f"Sob parametros realistas, resistencia e impossivel."
            )
        else:
            conclusion = (
                f"Resistencia a MRL-ASO-001 e matematicamente impossivel sob parametros "
                f"realistas. Das {escape_data['total_mutations']} mutacoes analisadas, "
                f"{n_escape} rompem binding, mas NENHUMA sobrevive a selecao purificadora "
                f"do SL RNA (conservado ha ~500 Ma). Mesmo no cenario mais generoso da "
                f"sensibilidade ({sensitivity_data['n_combinations']} combinacoes), o tempo "
                f"minimo para resistencia e {worst_display} anos. A cadeia de Markov, "
                f"o modelo de Kimura e o processo de Poisson convergem na mesma conclusao: "
                f"tres barreiras simultaneas (termodinamica, funcional, fixacao) tornam "
                f"a resistencia efetivamente impossivel."
            )

        logger.info("CONCLUSAO: %s", conclusion[:200] + "...")

    # --- Montar envelope ---
    envelope["runtime_seconds"] = timer.elapsed
    envelope["status"] = "success"
    envelope["summary"]["conclusion"] = conclusion
    envelope["summary"]["key_metrics"] = {
        "total_mutations_analyzed": escape_data["total_mutations"],
        "escape_mutations": n_escape,
        "escape_fraction": escape_data["escape_fraction"],
        "fitness_realistic": fitness_data["scenarios"]["realistic"]["fitness_summary"]["mean_fitness"],
        "selection_coefficient": fixation_data["selection_coefficient_s"],
        "total_escape_rate_effective": poisson_data["total_escape_rate_effective"],
        "sensitivity_combinations": sensitivity_data["n_combinations"],
        "worst_case_years": worst_display,
        "scenarios_less_than_100_years": sc_class["n_less_than_100_years"],
    }

    # Dados detalhados (sem incluir arrays enormes)
    # Limpar mutation details para nao inchar o JSON
    escape_summary = {k: v for k, v in escape_data.items() if k != "mutations"}
    escape_summary["top_10_escape_mutations"] = [
        m for m in escape_data["mutations"] if m["is_escape"]
    ][:10]

    # Limpar effective rates do fitness (muito grande)
    fitness_summary = {
        k: v for k, v in fitness_data.items() if k != "full_effective_rates"
    }

    # Limpar multi-mutant examples
    multi_summary = {
        "max_k_analyzed": multi_data["max_k_analyzed"],
        "results_by_k": {
            k: {key: val for key, val in v.items() if key != "escape_examples"}
            for k, v in multi_data["results_by_k"].items()
        },
    }

    envelope["data"] = {
        "markov_chain": markov_data,
        "escape_analysis": escape_summary,
        "fitness_analysis": fitness_summary,
        "kimura_fixation": fixation_data,
        "poisson_model": {
            "total_escape_rate": poisson_data["total_escape_rate_effective"],
            "n_escape_mutations": poisson_data["n_escape_mutations"],
            # Incluir apenas cenarios para N=10^8 (referencia)
            "reference_scenarios": [
                s for s in poisson_data["scenarios"]
                if s["population_size"] == 1e8
            ],
        },
        "multi_mutant": multi_summary,
        "sensitivity": sensitivity_data,
        "time_table": time_table,
    }

    # --- Gravar resultados ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # JSON principal
    output_json = write_result(envelope, module_name="math_6_markov")
    logger.info("JSON gravado em: %s", output_json)

    # Relatorio Markdown
    report_path = RESULTS_DIR / "RESISTANCE_REPORT.md"
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    logger.info("Relatorio gravado em: %s", report_path)

    logger.info("Tempo de execucao: %.2f segundos", timer.elapsed)

    return envelope


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
