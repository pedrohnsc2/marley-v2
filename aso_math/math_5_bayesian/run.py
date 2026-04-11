"""Math 5 — Otimizacao bayesiana de parametros de design de ASO.

Busca no espaco de design gapmer LNA-DNA-LNA por variantes potencialmente
superiores ao MRL-ASO-001 usando um Gaussian Process (GP) como modelo
surrogado e Expected Improvement (EI) como funcao de aquisicao.

Espaco de busca (4 dimensoes continuas mapeadas para design discreto):
    - Comprimento do ASO: 20-30 nt
    - Tamanho do gap de DNA: 10-20 posicoes
    - Posicao de inicio no SL RNA: 0-14 (para 25-mer)
    - Assimetria do flanco LNA: distribuicao 5' vs 3'

Objetivos (combinados em score composto normalizado):
    f1: dG de binding (minimizar = ligacao mais forte)
    f2: dG de hairpin (maximizar = menos auto-estrutura)
    f3: dG de self-dimer (maximizar = menos auto-dimerizacao)
    f4: Tm (penalidade fora do range 55-75 C)
    f5: Score de off-target (minimizar = mais seletivo)

Protocolo:
    1. 50 designs aleatorios como semente (Latin hypercube-like)
    2. 200 iteracoes de otimizacao bayesiana
    3. Refit do GP a cada 25 iteracoes (hiperparametros otimizados)
    4. Analise de Pareto e comparacao com MRL-ASO-001
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from aso_math.config import (
    ASO_KNOWN_DG,
    ASO_KNOWN_TM,
    ASO_SEQUENCE,
    ASO_TARGET_START,
    ASO_TARGET_END,
    SL_SEQUENCE,
)
from aso_math.envelope import Timer
from aso_math.math_5_bayesian.acquisition import expected_improvement, select_next_point
from aso_math.math_5_bayesian.design_space import (
    ASODesign,
    mrl_aso_001_params,
    params_to_design,
    sample_random_designs,
)
from aso_math.math_5_bayesian.gaussian_process import (
    GaussianProcess,
    optimize_hyperparameters,
)
from aso_math.target_config import TargetConfig
from aso_math.thermo import (
    compute_dg,
    compute_hairpin_dg,
    compute_self_dimer_dg,
    compute_tm,
    gc_content,
)
from core.logger import get_logger

logger = get_logger("math_5_bayesian")

# ---------------------------------------------------------------------------
# Diretorio de saida
# ---------------------------------------------------------------------------

RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Constantes da otimizacao
# ---------------------------------------------------------------------------

# Numero de designs aleatorios iniciais (semente)
N_INITIAL_SAMPLES: int = 50

# Numero de iteracoes de otimizacao bayesiana
N_OPTIMIZATION_ITERS: int = 200

# Frequencia de refit dos hiperparametros do GP
REFIT_EVERY: int = 25

# Numero de candidatos aleatorios para maximizar EI
N_EI_CANDIDATES: int = 5000

# Semente aleatoria para reprodutibilidade
RANDOM_SEED: int = 42

# Pesos dos objetivos no score composto
# Pesos maiores = objetivo mais importante na composicao
WEIGHTS: dict[str, float] = {
    "dg_binding": 0.35,     # ligacao forte e o objetivo principal
    "hairpin": 0.15,        # evitar auto-estrutura
    "self_dimer": 0.15,     # evitar auto-dimerizacao
    "tm_penalty": 0.20,     # Tm no range pratico (55-75 C)
    "off_target": 0.15,     # seletividade contra genoma hospedeiro
}

# Range ideal de Tm (penalidade zero dentro deste range)
TM_TARGET_MIN: float = 55.0
TM_TARGET_MAX: float = 75.0

# Bounds do espaco de busca (todas as dimensoes em [0, 1])
PARAM_BOUNDS: np.ndarray = np.array([
    [0.0, 1.0],  # comprimento
    [0.0, 1.0],  # gap fraction
    [0.0, 1.0],  # posicao de inicio
    [0.0, 1.0],  # assimetria LNA
])


# ---------------------------------------------------------------------------
# Funcao objetivo multi-criterio
# ---------------------------------------------------------------------------


def evaluate_design(design: ASODesign) -> dict[str, float]:
    """Avalia um design ASO em todos os 5 objetivos.

    Calcula as metricas termodinamicas e de seletividade do design,
    e combina em um score composto normalizado. O score composto e
    o que o GP modela — valor mais negativo = design melhor.

    Args:
        design: Design ASO valido com sequencia e arquitetura.

    Returns:
        Dicionario com todos os objetivos individuais e o score composto.
    """
    seq = design.sequence

    # --- Objetivos brutos ---
    dg_binding = compute_dg(seq)              # kcal/mol, mais negativo = melhor
    hairpin_dg = compute_hairpin_dg(seq)      # kcal/mol, menos negativo = melhor
    self_dimer_dg = compute_self_dimer_dg(seq)  # kcal/mol, menos negativo = melhor
    tm = compute_tm(seq)                       # Celsius
    gc = gc_content(seq)

    # --- Normalizacao e transformacao ---

    # f1: dG binding — normalizar para [-1, 0] onde -1 = muito forte (~-40 kcal/mol)
    # Range tipico: -15 a -40 kcal/mol
    f1_norm = max(-1.0, min(0.0, dg_binding / 40.0))

    # f2: Hairpin dG — normalizar para [-1, 0] onde 0 = sem hairpin
    # Range tipico: -5 a 0 kcal/mol
    f2_norm = max(-1.0, min(0.0, hairpin_dg / 5.0))
    # Invertemos: queremos maximizar (menos negativo = melhor)
    # Penalidade: mais negativo = pior, entao usamos o valor direto
    f2_penalty = -f2_norm  # 0 = sem hairpin (bom), 1 = hairpin forte (ruim)

    # f3: Self-dimer dG — mesmo tratamento que hairpin
    f3_norm = max(-1.0, min(0.0, self_dimer_dg / 10.0))
    f3_penalty = -f3_norm  # 0 = sem dimer (bom), 1 = dimer forte (ruim)

    # f4: Tm penalty — zero dentro do range ideal, cresce fora
    if TM_TARGET_MIN <= tm <= TM_TARGET_MAX:
        f4_penalty = 0.0
    elif tm < TM_TARGET_MIN:
        f4_penalty = min(1.0, (TM_TARGET_MIN - tm) / 20.0)
    else:
        f4_penalty = min(1.0, (tm - TM_TARGET_MAX) / 20.0)

    # f5: Off-target score — proxy simplificado baseado em GC content e comprimento
    # GC mais proximo de 50% = mais promiscuo; comprimento menor = mais off-targets
    # Score heuristico: P(match) ~ genome_size / 4^length
    host_genome_size = 2_400_000_000  # genoma canino ~2.4 Gb
    p_match = host_genome_size / (4 ** design.length)
    # Normalizar em escala log
    # Para L=20: p_match ~ 2.2e-3; para L=30: p_match ~ 2.2e-9
    f5_raw = min(1.0, max(0.0, -np.log10(max(p_match, 1e-15)) / 15.0))
    # Invertemos: f5_raw alto = poucos off-targets = bom
    f5_penalty = 1.0 - f5_raw

    # --- Score composto (minimizar) ---
    # Todos os componentes ja estao em [0, ~1], onde 0 = ideal
    # f1_norm ja e negativo, e queremos minimizar, entao contribui diretamente
    composite = (
        WEIGHTS["dg_binding"] * f1_norm          # mais negativo = melhor
        + WEIGHTS["hairpin"] * f2_penalty         # menor = melhor
        + WEIGHTS["self_dimer"] * f3_penalty      # menor = melhor
        + WEIGHTS["tm_penalty"] * f4_penalty      # menor = melhor
        + WEIGHTS["off_target"] * f5_penalty      # menor = melhor
    )

    return {
        "dg_binding_kcal": round(dg_binding, 4),
        "hairpin_dg_kcal": round(hairpin_dg, 4),
        "self_dimer_dg_kcal": round(self_dimer_dg, 4),
        "tm_celsius": round(tm, 2),
        "gc_content": round(gc, 4),
        "f1_binding_norm": round(f1_norm, 6),
        "f2_hairpin_penalty": round(f2_penalty, 6),
        "f3_dimer_penalty": round(f3_penalty, 6),
        "f4_tm_penalty": round(f4_penalty, 6),
        "f5_offtarget_penalty": round(f5_penalty, 6),
        "composite_score": round(composite, 6),
    }


# ---------------------------------------------------------------------------
# Analise de Pareto
# ---------------------------------------------------------------------------


def is_dominated(obj_a: np.ndarray, obj_b: np.ndarray) -> bool:
    """Verifica se obj_a e dominado por obj_b (minimizacao em todos os eixos).

    obj_a e dominado se obj_b e melhor ou igual em TODOS os objetivos
    e estritamente melhor em pelo menos um.
    """
    return bool(np.all(obj_b <= obj_a) and np.any(obj_b < obj_a))


def compute_pareto_front(
    objectives: np.ndarray,
) -> np.ndarray:
    """Identifica os indices dos pontos na frente de Pareto.

    Para um problema de minimizacao multi-objetivo, um ponto esta na
    frente de Pareto se nenhum outro ponto o domina em todos os objetivos.

    Args:
        objectives: Matriz (n, k) de valores objetivo (k objetivos).

    Returns:
        Array de indices dos pontos Pareto-otimos.
    """
    n = objectives.shape[0]
    is_pareto = np.ones(n, dtype=bool)

    for i in range(n):
        if not is_pareto[i]:
            continue
        for j in range(n):
            if i == j or not is_pareto[j]:
                continue
            if is_dominated(objectives[i], objectives[j]):
                is_pareto[i] = False
                break

    return np.where(is_pareto)[0]


def distance_to_pareto_front(
    point: np.ndarray,
    pareto_points: np.ndarray,
) -> float:
    """Calcula a distancia euclidiana minima de um ponto a frente de Pareto.

    Args:
        point: Vetor de objetivos do ponto de interesse, shape (k,).
        pareto_points: Matriz dos pontos Pareto, shape (m, k).

    Returns:
        Distancia minima. Zero se o ponto esta na frente.
    """
    if pareto_points.shape[0] == 0:
        return float("inf")
    distances = np.sqrt(np.sum((pareto_points - point) ** 2, axis=1))
    return float(np.min(distances))


# ---------------------------------------------------------------------------
# Loop principal de otimizacao bayesiana
# ---------------------------------------------------------------------------


def run_bayesian_optimization(
    sl_sequence: str,
    n_initial: int = N_INITIAL_SAMPLES,
    n_iterations: int = N_OPTIMIZATION_ITERS,
    seed: int = RANDOM_SEED,
) -> dict[str, Any]:
    """Executa a otimizacao bayesiana completa.

    Protocolo:
        1. Gerar n_initial designs aleatorios e avaliar
        2. Incluir MRL-ASO-001 como ponto de referencia
        3. Ajustar GP ao score composto
        4. Para cada iteracao: selecionar proximo ponto via EI, avaliar, atualizar GP
        5. Refit hiperparametros a cada REFIT_EVERY iteracoes

    Args:
        sl_sequence: Sequencia do SL RNA alvo.
        n_initial: Numero de designs iniciais aleatorios.
        n_iterations: Numero de iteracoes de otimizacao.
        seed: Semente aleatoria.

    Returns:
        Dicionario com todos os resultados da otimizacao.
    """
    rng = np.random.default_rng(seed)

    # --- Fase 1: Amostragem inicial ---
    logger.info("Fase 1: Gerando %d designs aleatorios iniciais...", n_initial)

    initial_designs = sample_random_designs(n_initial, rng, sl_sequence)
    logger.info("  -> %d designs validos obtidos", len(initial_designs))

    # Incluir MRL-ASO-001 como referencia obrigatoria
    mrl_params, mrl_design = mrl_aso_001_params(sl_sequence)

    # Coletar todos os designs avaliados
    all_params: list[np.ndarray] = []
    all_designs: list[ASODesign] = []
    all_objectives: list[dict[str, float]] = []
    all_composite: list[float] = []

    # Avaliar designs iniciais
    for params, design in initial_designs:
        obj = evaluate_design(design)
        all_params.append(params)
        all_designs.append(design)
        all_objectives.append(obj)
        all_composite.append(obj["composite_score"])

    # Avaliar MRL-ASO-001
    mrl_obj = evaluate_design(mrl_design)
    mrl_composite = mrl_obj["composite_score"]
    all_params.append(mrl_params)
    all_designs.append(mrl_design)
    all_objectives.append(mrl_obj)
    all_composite.append(mrl_composite)

    logger.info(
        "  MRL-ASO-001: composite = %.6f, dG = %.2f, Tm = %.2f",
        mrl_composite, mrl_obj["dg_binding_kcal"], mrl_obj["tm_celsius"],
    )

    mrl_index = len(all_composite) - 1  # indice do MRL-ASO-001

    # --- Fase 2: Otimizacao bayesiana ---
    logger.info("Fase 2: Iniciando %d iteracoes de otimizacao bayesiana...", n_iterations)

    # Montar matrizes para o GP
    X = np.array(all_params)
    y = np.array(all_composite)

    # Normalizar y para media 0 e std 1 (melhora estabilidade do GP)
    y_mean = np.mean(y)
    y_std = max(np.std(y), 1e-8)
    y_norm = (y - y_mean) / y_std

    # Ajustar GP inicial com otimizacao de hiperparametros
    gp = optimize_hyperparameters(X, y_norm)
    logger.info(
        "  GP inicial: length_scale=%.3f, signal_var=%.3f, log_ml=%.2f",
        gp.length_scale, gp.signal_variance, gp.log_marginal_likelihood(),
    )

    # Rastreamento do progresso
    best_composite_history: list[float] = [float(np.min(y))]
    ei_history: list[float] = []

    for iteration in range(n_iterations):
        # Melhor valor observado (no espaco normalizado)
        f_best_norm = float(np.min(y_norm))

        # Selecionar proximo ponto via Expected Improvement
        next_params = select_next_point(
            gp=gp,
            f_best=f_best_norm,
            bounds=PARAM_BOUNDS,
            n_candidates=N_EI_CANDIDATES,
            rng=rng,
        )

        # Calcular EI do ponto selecionado (para diagnostico)
        ei_val = expected_improvement(
            next_params.reshape(1, -1), gp, f_best_norm,
        )
        ei_history.append(float(ei_val[0]))

        # Converter parametros em design e avaliar
        design = params_to_design(next_params, sl_sequence)
        if design is None:
            # Design invalido — amostrar novamente
            fallback = sample_random_designs(1, rng, sl_sequence)
            if fallback:
                next_params, design = fallback[0]
            else:
                continue

        obj = evaluate_design(design)

        # Adicionar aos dados
        all_params.append(next_params)
        all_designs.append(design)
        all_objectives.append(obj)
        all_composite.append(obj["composite_score"])

        # Atualizar matrizes
        X = np.vstack([X, next_params.reshape(1, -1)])
        y = np.append(y, obj["composite_score"])
        y_mean = np.mean(y)
        y_std = max(np.std(y), 1e-8)
        y_norm = (y - y_mean) / y_std

        # Rastrear progresso (minimo acumulado, monotonico)
        running_min = min(best_composite_history[-1], obj["composite_score"])
        best_composite_history.append(running_min)

        # Refit periodico dos hiperparametros
        if (iteration + 1) % REFIT_EVERY == 0:
            gp = optimize_hyperparameters(X, y_norm)
            logger.info(
                "  Iteracao %d/%d: best=%.6f, EI=%.6e, refit GP (ls=%.3f, sv=%.3f)",
                iteration + 1, n_iterations, running_min,
                ei_val[0], gp.length_scale, gp.signal_variance,
            )
        else:
            # Refit rapido com mesmos hiperparametros
            gp.fit(X, y_norm)

            if (iteration + 1) % 50 == 0:
                logger.info(
                    "  Iteracao %d/%d: best=%.6f, EI=%.6e",
                    iteration + 1, n_iterations, running_min, ei_val[0],
                )

    logger.info("Otimizacao concluida. %d avaliacoes totais.", len(all_composite))

    return {
        "all_params": all_params,
        "all_designs": all_designs,
        "all_objectives": all_objectives,
        "all_composite": all_composite,
        "mrl_index": mrl_index,
        "mrl_objective": mrl_obj,
        "mrl_composite": mrl_composite,
        "best_composite_history": best_composite_history,
        "ei_history": ei_history,
        "gp_final": {
            "length_scale": gp.length_scale,
            "signal_variance": gp.signal_variance,
            "noise_variance": gp.noise_variance,
            "log_marginal_likelihood": gp.log_marginal_likelihood(),
        },
        "n_total_evaluations": len(all_composite),
    }


# ---------------------------------------------------------------------------
# Analise dos resultados
# ---------------------------------------------------------------------------


def analyze_results(
    opt_results: dict[str, Any],
) -> dict[str, Any]:
    """Analisa os resultados da otimizacao: top variantes, Pareto, etc.

    Identifica os 5 melhores designs encontrados, constroi a frente de
    Pareto para os pares de objetivos mais relevantes, e avalia a posicao
    do MRL-ASO-001 em relacao ao otimo encontrado.

    Args:
        opt_results: Resultado de run_bayesian_optimization().

    Returns:
        Dicionario com analise completa.
    """
    all_objectives = opt_results["all_objectives"]
    all_composite = opt_results["all_composite"]
    all_designs = opt_results["all_designs"]
    mrl_index = opt_results["mrl_index"]
    mrl_composite = opt_results["mrl_composite"]
    mrl_obj = opt_results["mrl_objective"]

    n_total = len(all_composite)

    # --- Top 5 variantes (pelo score composto) ---
    sorted_indices = np.argsort(all_composite)
    top_5_indices = sorted_indices[:5].tolist()

    top_5_variants: list[dict[str, Any]] = []
    for rank, idx in enumerate(top_5_indices, 1):
        design = all_designs[idx]
        obj = all_objectives[idx]
        improvement = mrl_composite - all_composite[idx]

        top_5_variants.append({
            "rank": rank,
            "index": idx,
            "sequence": design.sequence,
            "length": design.length,
            "architecture": design.architecture_str,
            "start_pos": design.start_pos,
            "composite_score": round(all_composite[idx], 6),
            "improvement_over_mrl": round(improvement, 6),
            "objectives": obj,
            "is_mrl_aso_001": idx == mrl_index,
        })

    # --- Rank do MRL-ASO-001 ---
    mrl_rank = int(np.searchsorted(np.sort(all_composite), mrl_composite)) + 1

    # --- Analise de Pareto (dG binding vs hairpin penalty) ---
    # Construir matriz de objetivos para Pareto: [f1, f2, f3, f4, f5]
    obj_matrix = np.array([
        [
            o["f1_binding_norm"],
            o["f2_hairpin_penalty"],
            o["f3_dimer_penalty"],
            o["f4_tm_penalty"],
            o["f5_offtarget_penalty"],
        ]
        for o in all_objectives
    ])

    pareto_indices = compute_pareto_front(obj_matrix)
    pareto_points = obj_matrix[pareto_indices]

    # MRL-ASO-001 na frente de Pareto?
    mrl_is_pareto = mrl_index in pareto_indices

    # Distancia do MRL-ASO-001 a frente de Pareto
    mrl_point = obj_matrix[mrl_index]
    if mrl_is_pareto:
        mrl_pareto_distance = 0.0
    else:
        mrl_pareto_distance = distance_to_pareto_front(mrl_point, pareto_points)

    # --- Pareto front para visualizacao (f1 vs f2) ---
    # Subconjunto 2D: dG binding vs hairpin
    obj_2d = obj_matrix[:, :2]
    pareto_2d_indices = compute_pareto_front(obj_2d)

    pareto_2d_data: list[dict[str, Any]] = []
    for idx in pareto_2d_indices:
        design = all_designs[idx]
        pareto_2d_data.append({
            "index": int(idx),
            "sequence": design.sequence,
            "architecture": design.architecture_str,
            "f1_binding": float(obj_2d[idx, 0]),
            "f2_hairpin": float(obj_2d[idx, 1]),
            "composite_score": round(all_composite[idx], 6),
            "is_mrl_aso_001": int(idx) == mrl_index,
        })

    # --- Melhor vs MRL-ASO-001: comparacao detalhada ---
    best_idx = sorted_indices[0]
    best_design = all_designs[best_idx]
    best_obj = all_objectives[best_idx]
    best_composite = all_composite[best_idx]

    comparison = {
        "best_found": {
            "sequence": best_design.sequence,
            "architecture": best_design.architecture_str,
            "start_pos": best_design.start_pos,
            "composite_score": round(best_composite, 6),
            "objectives": best_obj,
        },
        "mrl_aso_001": {
            "sequence": all_designs[mrl_index].sequence,
            "architecture": all_designs[mrl_index].architecture_str,
            "start_pos": all_designs[mrl_index].start_pos,
            "composite_score": round(mrl_composite, 6),
            "objectives": mrl_obj,
        },
        "best_is_mrl": best_idx == mrl_index,
        "composite_difference": round(mrl_composite - best_composite, 6),
        "objective_differences": {
            "dg_binding": round(
                mrl_obj["dg_binding_kcal"] - best_obj["dg_binding_kcal"], 4
            ),
            "hairpin_dg": round(
                mrl_obj["hairpin_dg_kcal"] - best_obj["hairpin_dg_kcal"], 4
            ),
            "self_dimer_dg": round(
                mrl_obj["self_dimer_dg_kcal"] - best_obj["self_dimer_dg_kcal"], 4
            ),
            "tm": round(
                mrl_obj["tm_celsius"] - best_obj["tm_celsius"], 2
            ),
        },
    }

    # --- Estatisticas de convergencia ---
    history = opt_results["best_composite_history"]
    convergence = {
        "initial_best": round(history[0], 6),
        "final_best": round(history[-1], 6),
        "improvement": round(history[0] - history[-1], 6),
        "converged_at_iteration": _find_convergence_point(history),
        "history_length": len(history),
    }

    return {
        "top_5_variants": top_5_variants,
        "mrl_rank": mrl_rank,
        "mrl_rank_total": n_total,
        "mrl_percentile": round((n_total - mrl_rank) / max(n_total - 1, 1) * 100, 2),
        "pareto_analysis": {
            "n_pareto_optimal": len(pareto_indices),
            "pareto_indices": pareto_indices.tolist(),
            "mrl_is_pareto_optimal": mrl_is_pareto,
            "mrl_distance_to_pareto": round(mrl_pareto_distance, 6),
        },
        "pareto_2d_binding_vs_hairpin": pareto_2d_data,
        "comparison_best_vs_mrl": comparison,
        "convergence": convergence,
    }


def _find_convergence_point(history: list[float], threshold: float = 1e-4) -> int:
    """Encontra a iteracao em que a otimizacao convergiu.

    Convergencia: as ultimas 20 iteracoes tem variacao < threshold.
    Retorna o indice da primeira iteracao da janela estavel.
    """
    window = 20
    if len(history) < window:
        return len(history) - 1

    for i in range(len(history) - window):
        window_vals = history[i:i + window]
        if max(window_vals) - min(window_vals) < threshold:
            return i

    return len(history) - 1


# ---------------------------------------------------------------------------
# Gerador de relatorio Markdown
# ---------------------------------------------------------------------------


def _generate_report(
    analysis: dict[str, Any],
    opt_results: dict[str, Any],
    runtime_seconds: float,
) -> str:
    """Gera relatorio BAYESIAN_REPORT.md com os resultados da otimizacao."""
    now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    top_5 = analysis["top_5_variants"]
    comparison = analysis["comparison_best_vs_mrl"]
    pareto = analysis["pareto_analysis"]
    convergence = analysis["convergence"]
    gp_info = opt_results["gp_final"]

    lines: list[str] = []
    lines.append("# Bayesian Optimization of ASO Design Parameters")
    lines.append("")
    lines.append(f"**Generated:** {now}")
    lines.append(f"**Runtime:** {runtime_seconds:.2f} seconds")
    lines.append(f"**Total evaluations:** {opt_results['n_total_evaluations']}")
    lines.append(f"**Protocol:** {N_INITIAL_SAMPLES} initial + {N_OPTIMIZATION_ITERS} BO iterations")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- Resumo executivo ---
    lines.append("## Executive Summary")
    lines.append("")

    best_is_mrl = comparison["best_is_mrl"]
    if best_is_mrl:
        lines.append(
            "**MRL-ASO-001 is the best design found by Bayesian optimization.** "
            "No variant in the search space achieves a better composite score, "
            "confirming its optimality within the explored parameter space."
        )
    else:
        diff = comparison["composite_difference"]
        best_seq = comparison["best_found"]["sequence"]
        best_arch = comparison["best_found"]["architecture"]
        lines.append(
            f"**Bayesian optimization found a potentially better variant.** "
            f"The best design ({best_arch}, composite improvement: {diff:.6f}) "
            f"outperforms MRL-ASO-001 in the composite objective. "
            f"However, this result requires experimental validation — "
            f"the composite score is a simplified mathematical model."
        )
    lines.append("")

    mrl_rank = analysis["mrl_rank"]
    mrl_total = analysis["mrl_rank_total"]
    mrl_pct = analysis["mrl_percentile"]
    lines.append(
        f"MRL-ASO-001 ranks **#{mrl_rank}** out of {mrl_total} designs "
        f"evaluated (percentile: {mrl_pct:.1f}%)."
    )
    lines.append("")

    # --- Espaco de busca ---
    lines.append("## 1. Design Space")
    lines.append("")
    lines.append("| Parameter | Range | Description |")
    lines.append("|-----------|-------|-------------|")
    lines.append("| ASO length | 20-30 nt | Total oligonucleotide length |")
    lines.append("| Gap size | 10-20 nt | Central DNA gap for RNase H |")
    lines.append("| Start position | 0-14 | Position on SL RNA (for 25-mer) |")
    lines.append("| LNA asymmetry | [0, 1] | Distribution of LNA flanks |")
    lines.append("")
    lines.append("**Constraint:** LNA modifications only in flanking regions (gapmer design).")
    lines.append(f"**Minimum LNA flank:** 2 nt per side.")
    lines.append("")

    # --- Objetivos ---
    lines.append("## 2. Multi-Objective Function")
    lines.append("")
    lines.append("| Objective | Weight | Goal | Source |")
    lines.append("|-----------|--------|------|--------|")
    lines.append(f"| dG binding | {WEIGHTS['dg_binding']:.0%} | Minimize (stronger binding) | NN model |")
    lines.append(f"| Hairpin dG | {WEIGHTS['hairpin']:.0%} | Maximize (less self-structure) | NN model |")
    lines.append(f"| Self-dimer dG | {WEIGHTS['self_dimer']:.0%} | Maximize (less dimerization) | NN model |")
    lines.append(f"| Tm penalty | {WEIGHTS['tm_penalty']:.0%} | Zero in 55-75 C range | NN model |")
    lines.append(f"| Off-target | {WEIGHTS['off_target']:.0%} | Minimize (fewer off-targets) | GC/length proxy |")
    lines.append("")

    # --- GP model ---
    lines.append("## 3. Gaussian Process Surrogate Model")
    lines.append("")
    lines.append("- **Kernel:** RBF (Squared Exponential)")
    lines.append(f"- **Final length scale:** {gp_info['length_scale']:.4f}")
    lines.append(f"- **Final signal variance:** {gp_info['signal_variance']:.4f}")
    lines.append(f"- **Noise variance:** {gp_info['noise_variance']:.2e}")
    lines.append(f"- **Log marginal likelihood:** {gp_info['log_marginal_likelihood']:.2f}")
    lines.append("- **Hyperparameter optimization:** Grid search on log marginal likelihood")
    lines.append(f"- **Refit frequency:** Every {REFIT_EVERY} iterations")
    lines.append("")

    # --- Top 5 variantes ---
    lines.append("## 4. Top 5 Variants Found")
    lines.append("")
    lines.append(
        "| Rank | Sequence | Arch | Start | Composite | dG (kcal/mol) | Tm (C) | "
        "Hairpin dG | Self-dimer dG | Is MRL? |"
    )
    lines.append(
        "|------|----------|------|-------|-----------|---------------|--------|"
        "------------|--------------|---------|"
    )
    for v in top_5:
        obj = v["objectives"]
        is_mrl = "YES" if v["is_mrl_aso_001"] else "no"
        lines.append(
            f"| {v['rank']} | `{v['sequence']}` | {v['architecture']} | "
            f"{v['start_pos']} | {v['composite_score']:.6f} | "
            f"{obj['dg_binding_kcal']:.2f} | {obj['tm_celsius']:.1f} | "
            f"{obj['hairpin_dg_kcal']:.2f} | {obj['self_dimer_dg_kcal']:.2f} | "
            f"{is_mrl} |"
        )
    lines.append("")

    # --- Comparacao detalhada ---
    lines.append("## 5. Detailed Comparison: Best Found vs MRL-ASO-001")
    lines.append("")

    best = comparison["best_found"]
    mrl = comparison["mrl_aso_001"]

    lines.append("| Metric | Best Found | MRL-ASO-001 | Difference |")
    lines.append("|--------|-----------|-------------|------------|")
    lines.append(
        f"| Sequence | `{best['sequence']}` | `{mrl['sequence']}` | — |"
    )
    lines.append(
        f"| Architecture | {best['architecture']} | {mrl['architecture']} | — |"
    )
    lines.append(
        f"| Composite score | {best['composite_score']:.6f} | "
        f"{mrl['composite_score']:.6f} | "
        f"{comparison['composite_difference']:.6f} |"
    )
    lines.append(
        f"| dG binding | {best['objectives']['dg_binding_kcal']:.2f} | "
        f"{mrl['objectives']['dg_binding_kcal']:.2f} | "
        f"{comparison['objective_differences']['dg_binding']:.2f} kcal/mol |"
    )
    lines.append(
        f"| Hairpin dG | {best['objectives']['hairpin_dg_kcal']:.2f} | "
        f"{mrl['objectives']['hairpin_dg_kcal']:.2f} | "
        f"{comparison['objective_differences']['hairpin_dg']:.2f} kcal/mol |"
    )
    lines.append(
        f"| Self-dimer dG | {best['objectives']['self_dimer_dg_kcal']:.2f} | "
        f"{mrl['objectives']['self_dimer_dg_kcal']:.2f} | "
        f"{comparison['objective_differences']['self_dimer_dg']:.2f} kcal/mol |"
    )
    lines.append(
        f"| Tm | {best['objectives']['tm_celsius']:.1f} C | "
        f"{mrl['objectives']['tm_celsius']:.1f} C | "
        f"{comparison['objective_differences']['tm']:.1f} C |"
    )
    lines.append("")

    # --- Pareto analysis ---
    lines.append("## 6. Pareto Optimality Assessment")
    lines.append("")
    lines.append(
        f"- **Pareto-optimal designs (5D):** {pareto['n_pareto_optimal']} out of {mrl_total}"
    )
    lines.append(
        f"- **MRL-ASO-001 is Pareto-optimal:** "
        f"{'YES' if pareto['mrl_is_pareto_optimal'] else 'NO'}"
    )
    if not pareto["mrl_is_pareto_optimal"]:
        lines.append(
            f"- **Distance to Pareto front:** {pareto['mrl_distance_to_pareto']:.6f} "
            f"(normalized objective space)"
        )
    lines.append("")

    pareto_2d = analysis["pareto_2d_binding_vs_hairpin"]
    lines.append("### Pareto Front: dG Binding vs Hairpin Penalty (2D projection)")
    lines.append("")
    if pareto_2d:
        lines.append("| Arch | f1 (binding) | f2 (hairpin) | Composite | Is MRL? |")
        lines.append("|------|-------------|-------------|-----------|---------|")
        for p in sorted(pareto_2d, key=lambda x: x["f1_binding"]):
            is_mrl = "YES" if p["is_mrl_aso_001"] else "no"
            lines.append(
                f"| {p['architecture']} | {p['f1_binding']:.4f} | "
                f"{p['f2_hairpin']:.4f} | {p['composite_score']:.6f} | {is_mrl} |"
            )
        lines.append("")

    # --- Convergencia ---
    lines.append("## 7. Convergence Diagnostics")
    lines.append("")
    lines.append(f"- **Initial best composite:** {convergence['initial_best']:.6f}")
    lines.append(f"- **Final best composite:** {convergence['final_best']:.6f}")
    lines.append(f"- **Total improvement:** {convergence['improvement']:.6f}")
    lines.append(f"- **Converged at iteration:** ~{convergence['converged_at_iteration']}")
    lines.append("")

    # --- Caveats ---
    lines.append("## 8. Limitations & Caveats")
    lines.append("")
    lines.append(
        "1. **Simplified model:** The composite objective uses nearest-neighbor "
        "thermodynamics and heuristic off-target scoring. Real ASO efficacy depends "
        "on cellular uptake, RNase H recruitment, metabolic stability, and other "
        "factors not modeled here."
    )
    lines.append(
        "2. **LNA thermodynamics not modeled:** LNA modifications increase Tm by "
        "~3-5 C per modification, but the NN model uses DNA/DNA parameters. This "
        "means absolute Tm values are underestimates for LNA gapmers."
    )
    lines.append(
        "3. **4D parameter space:** The search explores a continuous relaxation of "
        "the discrete design space. Some promising designs in intermediate "
        "parameter regions may have been missed."
    )
    lines.append(
        "4. **No experimental validation:** All results are computational predictions. "
        "Any promising variant must be synthesized and tested in vitro before "
        "drawing conclusions about superiority."
    )
    lines.append("")

    # --- Metodo ---
    lines.append("## 9. Methods")
    lines.append("")
    lines.append("**Surrogate model:** Gaussian Process with RBF kernel, implemented from scratch.")
    lines.append("**Acquisition function:** Expected Improvement (EI) with xi=0.01.")
    lines.append("**Hyperparameter optimization:** Grid search maximizing log marginal likelihood.")
    lines.append(f"**Seed designs:** {N_INITIAL_SAMPLES} random + MRL-ASO-001 reference.")
    lines.append(f"**Optimization iterations:** {N_OPTIMIZATION_ITERS}.")
    lines.append(f"**EI maximization:** Random sampling with {N_EI_CANDIDATES} candidates per iteration.")
    lines.append("**Pareto analysis:** Exhaustive pairwise dominance check (5 objectives).")
    lines.append("")
    lines.append("**References:**")
    lines.append("- Rasmussen & Williams (2006) *Gaussian Processes for Machine Learning*")
    lines.append("- Jones, Schonlau & Welch (1998) *Efficient Global Optimization*")
    lines.append("- SantaLucia (1998) *PNAS* 95(4):1460-1465")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Funcao utilitaria para salvar JSON
# ---------------------------------------------------------------------------


def _save_json(data: Any, name: str) -> Path:
    """Salva dados como JSON no diretorio de resultados."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / f"{name}.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False, default=_json_default)
    return path


def _json_default(obj: Any) -> Any:
    """Serializador customizado para tipos nao-JSON (numpy, etc)."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    raise TypeError(f"Tipo nao serializavel: {type(obj)}")


# ---------------------------------------------------------------------------
# Entry point principal
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a otimizacao bayesiana completa de design ASO.

    Fluxo:
        1. Gerar designs aleatorios iniciais e avaliar MRL-ASO-001
        2. Ajustar GP surrogado e executar otimizacao via EI
        3. Analisar resultados: top variantes, Pareto, convergencia
        4. Gerar relatorio Markdown e salvar JSON
        5. Retornar resultados consolidados

    Args:
        config: Configuracao do organismo alvo. Se None, usa defaults MRL-ASO-001.

    Returns:
        Dicionario com todos os resultados da analise.
    """
    # --- Config padrao para L. infantum (retrocompativel) ---
    if config is None:
        config = TargetConfig(
            aso_sequence=ASO_SEQUENCE,
            aso_name="MRL-ASO-001",
            sl_sequence=SL_SEQUENCE,
            aso_target_start=ASO_TARGET_START,
            aso_target_end=ASO_TARGET_END,
            known_dg=ASO_KNOWN_DG,
            known_tm=ASO_KNOWN_TM,
        )

    sl_seq = config.sl_sequence.upper()
    aso_name = config.aso_name or "ASO"

    logger.info("=" * 70)
    logger.info("MATH 5 — Otimizacao Bayesiana de Design ASO — %s", aso_name)
    logger.info("=" * 70)

    with Timer() as timer:
        # ---------------------------------------------------------------
        # Passo 1: Otimizacao bayesiana
        # ---------------------------------------------------------------
        logger.info("Passo 1/3: Executando otimizacao bayesiana...")
        opt_results = run_bayesian_optimization(sl_sequence=sl_seq)

        # ---------------------------------------------------------------
        # Passo 2: Analise dos resultados
        # ---------------------------------------------------------------
        logger.info("Passo 2/3: Analisando resultados...")
        analysis = analyze_results(opt_results)

        # ---------------------------------------------------------------
        # Passo 3: Salvar resultados
        # ---------------------------------------------------------------
        logger.info("Passo 3/3: Salvando resultados...")

        # JSON com top variantes e analise (sem dados brutos massivos)
        analysis_serializable = _make_serializable(analysis)
        _save_json(analysis_serializable, "bayesian_analysis")
        logger.info("  -> bayesian_analysis.json salvo")

        # JSON com historico de convergencia
        convergence_data = {
            "best_composite_history": opt_results["best_composite_history"],
            "ei_history": opt_results["ei_history"],
            "gp_final": opt_results["gp_final"],
            "n_total_evaluations": opt_results["n_total_evaluations"],
        }
        _save_json(convergence_data, "bayesian_convergence")
        logger.info("  -> bayesian_convergence.json salvo")

    # --- Gerar relatorio Markdown ---
    report = _generate_report(
        analysis=analysis,
        opt_results=opt_results,
        runtime_seconds=timer.elapsed,
    )
    report_path = RESULTS_DIR / "BAYESIAN_REPORT.md"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    logger.info("Relatorio salvo em: %s", report_path)

    # --- Resultado consolidado ---
    comparison = analysis["comparison_best_vs_mrl"]
    pareto = analysis["pareto_analysis"]

    result = {
        "module": "math_5_bayesian",
        "version": "1.0.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": timer.elapsed,
        "status": "success",
        "aso_name": aso_name,
        "organism": config.species_name,
        "optimization": {
            "n_initial_samples": N_INITIAL_SAMPLES,
            "n_iterations": N_OPTIMIZATION_ITERS,
            "n_total_evaluations": opt_results["n_total_evaluations"],
        },
        "mrl_aso_001": {
            "rank": analysis["mrl_rank"],
            "rank_total": analysis["mrl_rank_total"],
            "percentile": analysis["mrl_percentile"],
            "composite_score": round(opt_results["mrl_composite"], 6),
            "is_pareto_optimal": pareto["mrl_is_pareto_optimal"],
            "distance_to_pareto": pareto["mrl_distance_to_pareto"],
        },
        "best_found": {
            "sequence": comparison["best_found"]["sequence"],
            "architecture": comparison["best_found"]["architecture"],
            "composite_score": comparison["best_found"]["composite_score"],
            "is_mrl": comparison["best_is_mrl"],
            "improvement_over_mrl": comparison["composite_difference"],
        },
        "pareto": {
            "n_pareto_optimal": pareto["n_pareto_optimal"],
            "mrl_is_pareto": pareto["mrl_is_pareto_optimal"],
        },
        "convergence": analysis["convergence"],
        "top_5_sequences": [v["sequence"] for v in analysis["top_5_variants"]],
    }

    logger.info("=" * 70)
    logger.info("MATH 5 COMPLETO — Tempo total: %.2f s", timer.elapsed)
    logger.info(
        "  MRL-ASO-001: rank %d/%d (%.1f%%), Pareto-otimo: %s",
        analysis["mrl_rank"], analysis["mrl_rank_total"],
        analysis["mrl_percentile"], pareto["mrl_is_pareto_optimal"],
    )
    if not comparison["best_is_mrl"]:
        logger.info(
            "  Melhor variante encontrada: %s (%s), melhoria = %.6f",
            comparison["best_found"]["sequence"],
            comparison["best_found"]["architecture"],
            comparison["composite_difference"],
        )
    logger.info("=" * 70)

    return result


def _make_serializable(data: Any) -> Any:
    """Converte recursivamente dados para tipos JSON-serializaveis."""
    if isinstance(data, dict):
        return {k: _make_serializable(v) for k, v in data.items()}
    if isinstance(data, list):
        return [_make_serializable(item) for item in data]
    if isinstance(data, np.ndarray):
        return data.tolist()
    if isinstance(data, (np.integer,)):
        return int(data)
    if isinstance(data, (np.floating,)):
        return float(data)
    if isinstance(data, np.bool_):
        return bool(data)
    if hasattr(data, "to_dict"):
        return data.to_dict()
    return data


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
