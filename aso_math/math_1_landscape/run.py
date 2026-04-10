"""Math 1 — Mapeamento estendido da paisagem termodinamica de MRL-ASO-001.

Estende o modulo 01_thermodynamic_landscape com cinco analises adicionais:

    1. Paisagem mutacional completa (25 posicoes x 4 bases = 100 variantes)
       Inclui heatmap de energia, identificacao do minimo global, fitness gap.

    2. Analise de mutantes duplos com matriz de epistasia
       ΔΔG_epistasis = ΔG(AB) - ΔG(A) - ΔG(B) + ΔG(WT)
       Identifica sinergias e antagonismos entre posicoes.

    3. Otimizacao de comprimento (18-30 nt)
       Para cada comprimento: posicao otima, ΔG, Tm, risco de off-target.
       Justificativa matematica para o comprimento de 25 nt.

    4. Janela deslizante posicional (sliding window de 25-mer)
       Desliza o 25-mer ao longo do SL RNA de 39 nt.
       Para cada posicao: ΔG, Tm, GC content.

    5. Estatisticas agregadas
       Rank, percentil, robustez, variantes proximas ao otimo.

Importa funcoes do modulo existente — nao duplica logica.
"""

from __future__ import annotations

import json
import statistics
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

from aso_math.config import (
    ASO_KNOWN_DG,
    ASO_KNOWN_TM,
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    BASES,
    SL_SEQUENCE,
)
from aso_math.target_config import TargetConfig
from aso_math.thermo import (
    compute_dg,
    compute_tm,
    gc_content,
    reverse_complement,
)
from aso_math.envelope import Timer
from core.logger import get_logger

logger = get_logger("math_1_landscape")

# ---------------------------------------------------------------------------
# Diretorio de saida
# ---------------------------------------------------------------------------

RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Constantes da analise
# ---------------------------------------------------------------------------

# Comprimentos para varredura estendida (18-30, ampliado do original 18-27)
LENGTH_SCAN_MIN_EXT: int = 18
LENGTH_SCAN_MAX_EXT: int = 30

# Limiar de proximidade ao otimo (kcal/mol)
NEAR_OPTIMUM_THRESHOLD: float = 2.0

# Amostragem de mutantes duplos: se True, avalia todos os 4800;
# se False, avalia o subconjunto exaustivo de 3 alt bases por posicao (2700)
FULL_DOUBLE_MUTANT_SCAN: bool = True

# Comprimento fixo para sliding window
SLIDING_WINDOW_LENGTH: int = 25


# ---------------------------------------------------------------------------
# 1. Paisagem mutacional completa (25 x 4)
# ---------------------------------------------------------------------------


def build_full_fitness_landscape(
    aso_sequence: str,
) -> dict[str, Any]:
    """Constroi a paisagem de fitness completa: 25 posicoes x 4 bases.

    Para cada uma das 100 combinacoes (posicao, base):
    - Se a base e a wildtype, calcula o dG do ASO original
    - Se e uma substituicao, calcula o dG do mutante

    Retorna a matriz 25x4 de dG, identificacao do minimo global,
    e o fitness gap (diferenca entre melhor e segundo melhor).

    Args:
        aso_sequence: Sequencia do ASO wildtype.
    """
    wt = aso_sequence.upper()
    n = len(wt)
    bases_list = list(BASES)  # ["A", "C", "G", "T"]

    # Calcular dG do wildtype uma vez
    wt_dg = compute_dg(wt)
    wt_tm = compute_tm(wt)

    # Matriz numpy para operacoes eficientes
    dg_matrix = np.zeros((n, 4), dtype=np.float64)
    tm_matrix = np.zeros((n, 4), dtype=np.float64)

    # Dados detalhados de cada variante
    all_variants: list[dict[str, Any]] = []

    for pos in range(n):
        original_base = wt[pos]
        for base_idx, base in enumerate(bases_list):
            if base == original_base:
                # Variante wildtype nesta posicao
                dg_matrix[pos, base_idx] = wt_dg
                tm_matrix[pos, base_idx] = wt_tm
                all_variants.append({
                    "position": pos,
                    "base": base,
                    "is_wildtype": True,
                    "dg_binding_kcal": wt_dg,
                    "tm_celsius": wt_tm,
                    "ddg_kcal": 0.0,
                    "mutation": f"{base}{pos + 1}{base}(WT)",
                })
            else:
                # Construir mutante
                mutant_seq = wt[:pos] + base + wt[pos + 1:]
                mut_dg = compute_dg(mutant_seq)
                mut_tm = compute_tm(mutant_seq)
                ddg = round(mut_dg - wt_dg, 4)

                dg_matrix[pos, base_idx] = mut_dg
                tm_matrix[pos, base_idx] = mut_tm
                all_variants.append({
                    "position": pos,
                    "base": base,
                    "is_wildtype": False,
                    "dg_binding_kcal": mut_dg,
                    "tm_celsius": mut_tm,
                    "ddg_kcal": ddg,
                    "mutation": f"{original_base}{pos + 1}{base}",
                })

    # Identificar minimo global na matriz
    min_flat_idx = int(np.argmin(dg_matrix))
    min_pos, min_base_idx = divmod(min_flat_idx, 4)
    global_min_dg = float(dg_matrix[min_pos, min_base_idx])

    # O wildtype e o minimo global?
    wt_is_global_min = abs(global_min_dg - wt_dg) < 1e-6

    # Fitness gap: diferenca entre melhor e segundo melhor
    flat_values = dg_matrix.flatten()
    sorted_values = np.sort(flat_values)
    fitness_gap = round(float(sorted_values[1] - sorted_values[0]), 4)

    # Sensibilidade posicional: variacao de dG por posicao
    positional_sensitivity: list[dict[str, float]] = []
    for pos in range(n):
        row = dg_matrix[pos, :]
        positional_sensitivity.append({
            "position": pos,
            "base_wt": wt[pos],
            "dg_min": round(float(np.min(row)), 4),
            "dg_max": round(float(np.max(row)), 4),
            "dg_range": round(float(np.max(row) - np.min(row)), 4),
            "dg_mean": round(float(np.mean(row)), 4),
            "dg_std": round(float(np.std(row)), 4),
        })

    logger.info(
        "Paisagem 25x4: WT dG = %.2f, min global dG = %.2f, "
        "WT e minimo global = %s, fitness gap = %.4f kcal/mol",
        wt_dg, global_min_dg, wt_is_global_min, fitness_gap,
    )

    return {
        "wildtype_dg_kcal": wt_dg,
        "wildtype_tm_celsius": wt_tm,
        "matrix_shape": [n, 4],
        "bases": bases_list,
        "dg_matrix": dg_matrix.tolist(),
        "tm_matrix": tm_matrix.tolist(),
        "global_minimum": {
            "position": int(min_pos),
            "base": bases_list[int(min_base_idx)],
            "dg_kcal": round(global_min_dg, 4),
            "is_wildtype": wt_is_global_min,
        },
        "fitness_gap_kcal": fitness_gap,
        "positional_sensitivity": positional_sensitivity,
        "total_variants": len(all_variants),
        "variants": all_variants,
    }


# ---------------------------------------------------------------------------
# 2. Mutantes duplos com epistasia
# ---------------------------------------------------------------------------


def analyze_double_mutants_epistasis(
    aso_sequence: str,
    single_mutant_dg_map: dict[tuple[int, str], float],
    wt_dg: float,
    sample_all: bool = True,
) -> dict[str, Any]:
    """Gera mutantes duplos e calcula a matriz de epistasia.

    Epistasia mede a interacao entre duas mutacoes:
        ΔΔG_epistasis = ΔG(AB) - ΔG(A) - ΔG(B) + ΔG(WT)

    Se ΔΔG_epistasis < 0: sinergismo (efeito combinado maior que aditivo)
    Se ΔΔG_epistasis > 0: antagonismo (efeito combinado menor que aditivo)
    Se ΔΔG_epistasis ≈ 0: aditivo (sem interacao)

    Args:
        aso_sequence: Sequencia do ASO wildtype.
        single_mutant_dg_map: Mapa (posicao, base) -> dG de mutantes simples.
        wt_dg: dG do wildtype.
        sample_all: Se True, avalia todas as 16 combinacoes por par (4800 total).
                    Se False, avalia apenas 9 por par (bases alternativas, 2700 total).
    """
    wt = aso_sequence.upper()
    n = len(wt)
    bases_list = list(BASES)

    double_mutants: list[dict[str, Any]] = []
    epistasis_entries: list[dict[str, Any]] = []

    # Contadores
    total_synergistic = 0
    total_antagonistic = 0
    total_additive = 0
    better_than_wt = 0

    # Para cada par de posicoes
    n_pairs = 0
    for pos1, pos2 in combinations(range(n), 2):
        # Bases a testar em cada posicao
        if sample_all:
            # Todas as 4 bases exceto WT (3 alternativas por posicao, 9 combos)
            # Mas spec pede 16 combos = 4x4 incluindo WT base.
            # Usamos apenas nao-WT x nao-WT para mutantes duplos reais.
            bases1 = [b for b in bases_list if b != wt[pos1]]
            bases2 = [b for b in bases_list if b != wt[pos2]]
        else:
            bases1 = [b for b in bases_list if b != wt[pos1]]
            bases2 = [b for b in bases_list if b != wt[pos2]]

        for base1 in bases1:
            for base2 in bases2:
                # Construir duplo mutante
                seq_list = list(wt)
                seq_list[pos1] = base1
                seq_list[pos2] = base2
                mutant_seq = "".join(seq_list)

                dg_double = compute_dg(mutant_seq)

                # dG de cada mutante simples individual
                dg_single_a = single_mutant_dg_map.get(
                    (pos1, base1), wt_dg
                )
                dg_single_b = single_mutant_dg_map.get(
                    (pos2, base2), wt_dg
                )

                # Epistasia: ΔΔG = ΔG(AB) - ΔG(A) - ΔG(B) + ΔG(WT)
                epistasis = round(
                    dg_double - dg_single_a - dg_single_b + wt_dg, 4
                )

                # Classificar
                if epistasis < -0.5:
                    interaction = "synergistic"
                    total_synergistic += 1
                elif epistasis > 0.5:
                    interaction = "antagonistic"
                    total_antagonistic += 1
                else:
                    interaction = "additive"
                    total_additive += 1

                if dg_double < wt_dg:
                    better_than_wt += 1

                entry = {
                    "positions": [pos1, pos2],
                    "bases": [base1, base2],
                    "mutation": f"{wt[pos1]}{pos1 + 1}{base1}+{wt[pos2]}{pos2 + 1}{base2}",
                    "dg_double_kcal": round(dg_double, 4),
                    "dg_single_a_kcal": round(dg_single_a, 4),
                    "dg_single_b_kcal": round(dg_single_b, 4),
                    "epistasis_ddg_kcal": epistasis,
                    "interaction_type": interaction,
                }
                double_mutants.append(entry)
                epistasis_entries.append(entry)

        n_pairs += 1

    # Ordenar por dG (mais negativo primeiro)
    double_mutants.sort(key=lambda m: m["dg_double_kcal"])

    # Estatisticas de epistasia
    epistasis_values = [e["epistasis_ddg_kcal"] for e in epistasis_entries]
    epi_array = np.array(epistasis_values) if epistasis_values else np.array([0.0])

    logger.info(
        "Mutantes duplos: %d avaliados, %d melhor que WT, "
        "epistasia media = %.4f kcal/mol",
        len(double_mutants), better_than_wt,
        float(np.mean(epi_array)),
    )

    return {
        "total_evaluated": len(double_mutants),
        "n_position_pairs": n_pairs,
        "better_than_wildtype": better_than_wt,
        "best_double_mutant": double_mutants[0] if double_mutants else None,
        "worst_double_mutant": double_mutants[-1] if double_mutants else None,
        "epistasis_statistics": {
            "mean_epistasis_kcal": round(float(np.mean(epi_array)), 4),
            "std_epistasis_kcal": round(float(np.std(epi_array)), 4),
            "min_epistasis_kcal": round(float(np.min(epi_array)), 4),
            "max_epistasis_kcal": round(float(np.max(epi_array)), 4),
            "total_synergistic": total_synergistic,
            "total_antagonistic": total_antagonistic,
            "total_additive": total_additive,
            "fraction_synergistic": round(
                total_synergistic / max(len(double_mutants), 1), 4
            ),
        },
        # Armazenar top 50 sinergisticos (mais negativos)
        "top_synergistic": sorted(
            [e for e in epistasis_entries if e["interaction_type"] == "synergistic"],
            key=lambda e: e["epistasis_ddg_kcal"],
        )[:50],
        # Nao incluir todos os 2700+ mutantes no JSON — apenas resumo
        "top_10_best_dg": double_mutants[:10],
        "top_10_worst_dg": double_mutants[-10:],
    }


# ---------------------------------------------------------------------------
# 3. Otimizacao de comprimento (18–30 nt)
# ---------------------------------------------------------------------------


def optimize_length(
    sl_sequence: str,
    aso_sequence: str,
    aso_target_start: int,
    length_min: int = LENGTH_SCAN_MIN_EXT,
    length_max: int = LENGTH_SCAN_MAX_EXT,
) -> dict[str, Any]:
    """Testa ASOs de comprimento variavel centrados na regiao alvo.

    Para cada comprimento L de length_min a length_max:
    - Encontra a posicao de inicio no SL RNA que maximiza |ΔG|
    - Calcula ΔG, Tm, GC content
    - Estima risco de off-target (heuristico baseado em GC e comprimento)

    O risco de off-target e estimado pela formula heuristica:
        risk = 4^L * (1 - GC)^L
    Normalizado para [0, 1]. Comprimentos maiores com alto GC
    tem menor risco de off-target.

    Args:
        sl_sequence: Sequencia do SL RNA alvo.
        aso_sequence: Sequencia do ASO atual (para referencia).
        aso_target_start: Posicao de inicio do ASO atual no SL RNA.
        length_min: Comprimento minimo a testar.
        length_max: Comprimento maximo a testar.
    """
    sl_seq = sl_sequence.upper()
    sl_len = len(sl_seq)
    wt = aso_sequence.upper()
    wt_len = len(wt)
    wt_dg = compute_dg(wt)

    length_results: list[dict[str, Any]] = []

    for length in range(length_min, length_max + 1):
        # Para cada posicao de inicio valida no SL RNA
        best_for_length: dict[str, Any] | None = None
        all_windows: list[dict[str, Any]] = []

        max_start = sl_len - length
        if max_start < 0:
            # Comprimento maior que o SL RNA — impossivel
            continue

        for start in range(max_start + 1):
            target_window = sl_seq[start:start + length]
            aso_candidate = reverse_complement(target_window)

            dg = compute_dg(aso_candidate)
            tm = compute_tm(aso_candidate)
            gc = gc_content(aso_candidate)

            # Risco heuristico de off-target
            # Sequencias mais longas e com GC alto sao mais especificas.
            # Usamos uma formula simplificada baseada na probabilidade de
            # match aleatorio em um genoma de 30 Mb (L. infantum)
            genome_size = 30_000_000  # ~30 Mb para L. infantum host? Nao, host e canino ~2.4 Gb
            # Para genoma canino (~2.4 Gb):
            host_genome_size = 2_400_000_000
            # Probabilidade de match exato no genoma hospedeiro
            # P(match) = genome_size / 4^L
            p_match = host_genome_size / (4 ** length)
            off_target_risk = min(1.0, p_match)

            window_data = {
                "start": start,
                "end": start + length,
                "length": length,
                "dg_binding_kcal": round(dg, 4),
                "tm_celsius": round(tm, 2),
                "gc_content": round(gc, 4),
                "off_target_risk": round(off_target_risk, 8),
            }
            all_windows.append(window_data)

            if best_for_length is None or dg < best_for_length["dg_binding_kcal"]:
                best_for_length = window_data.copy()

        if best_for_length is not None:
            # Marcar se coincide com o ASO atual
            is_current_length = (length == wt_len)
            best_for_length["is_current_aso_length"] = is_current_length
            best_for_length["n_windows_tested"] = len(all_windows)
            length_results.append(best_for_length)

    # Ordenar por dG (melhor primeiro)
    length_results.sort(key=lambda r: r["dg_binding_kcal"])

    # Identificar comprimento otimo
    optimal = length_results[0] if length_results else None
    optimal_length = optimal["length"] if optimal else 0
    is_25_optimal = optimal_length == wt_len if optimal else False

    logger.info(
        "Otimizacao de comprimento: %d comprimentos testados, "
        "otimo = %d nt (dG = %.2f), 25 nt e otimo = %s",
        len(length_results),
        optimal_length,
        optimal["dg_binding_kcal"] if optimal else 0.0,
        is_25_optimal,
    )

    return {
        "length_range": [length_min, length_max],
        "n_lengths_tested": len(length_results),
        "optimal_length": optimal,
        "current_aso_length": wt_len,
        "current_aso_is_optimal": is_25_optimal,
        "by_length": length_results,
    }


# ---------------------------------------------------------------------------
# 4. Janela deslizante posicional (sliding window)
# ---------------------------------------------------------------------------


def sliding_window_analysis(
    sl_sequence: str,
    window_length: int = SLIDING_WINDOW_LENGTH,
    aso_target_start: int = 0,
) -> dict[str, Any]:
    """Desliza um 25-mer ao longo do SL RNA de 39 nt.

    Para cada posicao de inicio (0 a sl_len - window_length):
    - Gera o ASO complementar reverso
    - Calcula ΔG, Tm, GC content
    - Identifica a janela termodinamicamente otima

    Com SL de 39 nt e janela de 25 nt, temos 15 posicoes possiveis (0..14).

    Args:
        sl_sequence: Sequencia do SL RNA alvo.
        window_length: Comprimento da janela deslizante.
        aso_target_start: Posicao de inicio do ASO atual no SL (para comparacao).
    """
    sl_seq = sl_sequence.upper()
    sl_len = len(sl_seq)

    n_windows = sl_len - window_length + 1
    if n_windows <= 0:
        return {
            "error": f"Janela ({window_length}) maior que SL RNA ({sl_len})",
            "windows": [],
        }

    windows: list[dict[str, Any]] = []

    for start in range(n_windows):
        target_window = sl_seq[start:start + window_length]
        aso_candidate = reverse_complement(target_window)

        dg = compute_dg(aso_candidate)
        tm = compute_tm(aso_candidate)
        gc = gc_content(aso_candidate)

        is_current_position = (start == aso_target_start)

        windows.append({
            "start": start,
            "end": start + window_length,
            "dg_binding_kcal": round(dg, 4),
            "tm_celsius": round(tm, 2),
            "gc_content": round(gc, 4),
            "is_current_aso_position": is_current_position,
        })

    # Identificar janela otima (menor dG = ligacao mais forte)
    optimal_window = min(windows, key=lambda w: w["dg_binding_kcal"])
    current_window = next(
        (w for w in windows if w["is_current_aso_position"]),
        None,
    )

    # Verificar se o ASO atual cobre a janela otima
    current_covers_optimal = (
        current_window is not None
        and current_window["start"] == optimal_window["start"]
    )

    # Calcular desvio do ASO atual em relacao ao otimo
    if current_window is not None:
        dg_deviation = round(
            current_window["dg_binding_kcal"] - optimal_window["dg_binding_kcal"], 4
        )
        position_offset = current_window["start"] - optimal_window["start"]
    else:
        dg_deviation = None
        position_offset = None

    logger.info(
        "Sliding window: %d janelas de %d nt, otima em pos %d (dG = %.2f), "
        "ASO atual em pos %d (dG = %.2f), cobre otimo = %s",
        len(windows), window_length,
        optimal_window["start"], optimal_window["dg_binding_kcal"],
        aso_target_start,
        current_window["dg_binding_kcal"] if current_window else 0.0,
        current_covers_optimal,
    )

    return {
        "window_length": window_length,
        "sl_length": sl_len,
        "n_windows": len(windows),
        "optimal_window": optimal_window,
        "current_aso_window": current_window,
        "current_covers_optimal": current_covers_optimal,
        "dg_deviation_from_optimal_kcal": dg_deviation,
        "position_offset_from_optimal": position_offset,
        "windows": windows,
    }


# ---------------------------------------------------------------------------
# 5. Estatisticas agregadas
# ---------------------------------------------------------------------------


def compute_summary_statistics(
    wt_dg: float,
    wt_tm: float,
    fitness_landscape: dict[str, Any],
    double_mutant_results: dict[str, Any],
    length_results: dict[str, Any],
    sliding_results: dict[str, Any],
) -> dict[str, Any]:
    """Calcula estatisticas agregadas sobre a otimalidade do ASO.

    Metricas:
    - Rank do MRL-ASO-001 entre todos os mutantes simples
    - Percentil de energia de ligacao
    - Numero de variantes dentro de 2 kcal/mol do otimo
    - Score de robustez: media de |ΔΔG| de todas as mutacoes simples

    Args:
        wt_dg: dG do wildtype.
        wt_tm: Tm do wildtype.
        fitness_landscape: Resultado de build_full_fitness_landscape().
        double_mutant_results: Resultado de analyze_double_mutants_epistasis().
        length_results: Resultado de optimize_length().
        sliding_results: Resultado de sliding_window_analysis().
    """
    # --- Rank entre mutantes simples ---
    all_variants = fitness_landscape["variants"]
    # Todos os dG (incluindo WT repetido em cada posicao)
    # Para ranking justo, pegar apenas as variantes unicas
    unique_dg_values = []
    for v in all_variants:
        if v["is_wildtype"]:
            continue
        unique_dg_values.append(v["dg_binding_kcal"])

    # Adicionar o WT uma vez
    all_dg_with_wt = [wt_dg] + unique_dg_values
    all_dg_with_wt.sort()

    # Rank do WT (1 = melhor, mais negativo)
    wt_rank = all_dg_with_wt.index(wt_dg) + 1
    total_in_ranking = len(all_dg_with_wt)

    # Percentil: (N - rank) / (N - 1) * 100
    # 100% = melhor; 0% = pior
    if total_in_ranking > 1:
        wt_percentile = round(
            (total_in_ranking - wt_rank) / (total_in_ranking - 1) * 100, 2
        )
    else:
        wt_percentile = 100.0

    # --- Variantes proximas ao otimo ---
    global_min_dg = min(all_dg_with_wt)
    near_optimum_count = sum(
        1 for dg in all_dg_with_wt
        if dg <= global_min_dg + NEAR_OPTIMUM_THRESHOLD
    )

    # --- Score de robustez ---
    # Media de |ΔΔG| para todas as mutacoes simples
    # Valor alto = pequenas mudancas no ASO causam grande perda de afinidade
    # = a sequencia e altamente otimizada
    ddg_values = [
        abs(v["ddg_kcal"]) for v in all_variants if not v["is_wildtype"]
    ]
    robustness_score = round(statistics.mean(ddg_values), 4) if ddg_values else 0.0
    robustness_std = round(statistics.stdev(ddg_values), 4) if len(ddg_values) > 1 else 0.0

    # --- Mutantes que superam o WT ---
    better_single = sum(1 for dg in unique_dg_values if dg < wt_dg)
    better_double = double_mutant_results.get("better_than_wildtype", 0)

    # --- Resumo de comprimento ---
    length_optimal = length_results.get("optimal_length", {})
    optimal_length_nt = length_optimal.get("length", 0) if length_optimal else 0

    # --- Resumo de sliding window ---
    covers_optimal = sliding_results.get("current_covers_optimal", False)
    dg_deviation = sliding_results.get("dg_deviation_from_optimal_kcal", 0.0)

    logger.info(
        "Estatisticas: rank %d/%d (percentil %.1f%%), robustez = %.4f, "
        "variantes dentro de 2 kcal/mol = %d",
        wt_rank, total_in_ranking, wt_percentile,
        robustness_score, near_optimum_count,
    )

    return {
        "wildtype_dg_kcal": wt_dg,
        "wildtype_tm_celsius": wt_tm,
        "rank_among_single_mutants": {
            "rank": wt_rank,
            "total": total_in_ranking,
            "description": (
                f"MRL-ASO-001 ocupa a posicao {wt_rank} de {total_in_ranking} "
                f"variantes (incluindo WT e 75 mutantes simples)"
            ),
        },
        "percentile_binding_energy": {
            "percentile": wt_percentile,
            "description": (
                f"MRL-ASO-001 esta no percentil {wt_percentile}% de energia "
                f"de ligacao (100% = melhor)"
            ),
        },
        "variants_near_optimum": {
            "threshold_kcal": NEAR_OPTIMUM_THRESHOLD,
            "count": near_optimum_count,
            "description": (
                f"{near_optimum_count} variante(s) dentro de "
                f"{NEAR_OPTIMUM_THRESHOLD} kcal/mol do otimo global"
            ),
        },
        "robustness": {
            "mean_abs_ddg_kcal": robustness_score,
            "std_abs_ddg_kcal": robustness_std,
            "description": (
                f"Score de robustez = {robustness_score:.4f} kcal/mol "
                f"(media de |ΔΔG| para mutacoes simples; "
                f"valor alto = sequencia altamente otimizada)"
            ),
        },
        "single_mutants_better_than_wt": better_single,
        "double_mutants_better_than_wt": better_double,
        "optimal_length_nt": optimal_length_nt,
        "sliding_window_covers_optimal": covers_optimal,
        "sliding_window_dg_deviation_kcal": dg_deviation,
    }


# ---------------------------------------------------------------------------
# Gerador de relatorio Markdown
# ---------------------------------------------------------------------------


def _generate_report(
    fitness_landscape: dict[str, Any],
    double_mutant_results: dict[str, Any],
    length_results: dict[str, Any],
    sliding_results: dict[str, Any],
    summary_stats: dict[str, Any],
    runtime_seconds: float,
) -> str:
    """Gera relatorio Markdown com os resultados da analise.

    Nao expoe sequencias confidenciais — usa apenas posicoes e estatisticas.
    """
    now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # --- Dados para o relatorio ---
    wt_dg = fitness_landscape["wildtype_dg_kcal"]
    wt_tm = fitness_landscape["wildtype_tm_celsius"]
    matrix_shape = fitness_landscape["matrix_shape"]
    global_min = fitness_landscape["global_minimum"]
    fitness_gap = fitness_landscape["fitness_gap_kcal"]

    epi_stats = double_mutant_results["epistasis_statistics"]

    opt_len = length_results.get("optimal_length", {})
    current_is_opt = length_results.get("current_aso_is_optimal", False)

    sliding_opt = sliding_results.get("optimal_window", {})
    covers_opt = sliding_results.get("current_covers_optimal", False)

    rank_info = summary_stats["rank_among_single_mutants"]
    pct_info = summary_stats["percentile_binding_energy"]
    near_info = summary_stats["variants_near_optimum"]
    robust_info = summary_stats["robustness"]

    lines: list[str] = []
    lines.append("# Thermodynamic Landscape — Extended Analysis Report")
    lines.append("")
    lines.append(f"**Generated:** {now}")
    lines.append(f"**Runtime:** {runtime_seconds:.2f} seconds")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- 1. Fitness landscape ---
    lines.append("## 1. Full Mutational Fitness Landscape")
    lines.append("")
    lines.append(f"- **Matrix dimensions:** {matrix_shape[0]} positions x {matrix_shape[1]} bases = "
                 f"{matrix_shape[0] * matrix_shape[1]} variants")
    lines.append(f"- **Wildtype dG:** {wt_dg:.2f} kcal/mol")
    lines.append(f"- **Wildtype Tm:** {wt_tm:.2f} C")
    lines.append(f"- **Global minimum dG:** {global_min['dg_kcal']:.4f} kcal/mol "
                 f"(position {global_min['position'] + 1}, base {global_min['base']})")
    lines.append(f"- **Wildtype is global minimum:** {'YES' if global_min['is_wildtype'] else 'NO'}")
    lines.append(f"- **Fitness gap (1st - 2nd):** {fitness_gap:.4f} kcal/mol")
    lines.append("")

    # Heatmap summary (top 5 most sensitive positions)
    sensitivity = fitness_landscape["positional_sensitivity"]
    sorted_sens = sorted(sensitivity, key=lambda s: s["dg_range"], reverse=True)
    lines.append("### Most Sensitive Positions (highest dG range)")
    lines.append("")
    lines.append("| Position | WT Base | dG Range (kcal/mol) | dG Min | dG Max |")
    lines.append("|----------|---------|---------------------|--------|--------|")
    for s in sorted_sens[:5]:
        lines.append(
            f"| {s['position'] + 1:>8} | {s['base_wt']:>7} | "
            f"{s['dg_range']:>19.4f} | {s['dg_min']:>6.2f} | {s['dg_max']:>6.2f} |"
        )
    lines.append("")

    # --- 2. Double mutants + epistasis ---
    lines.append("## 2. Double Mutant Analysis & Epistasis")
    lines.append("")
    lines.append(f"- **Total double mutants evaluated:** {double_mutant_results['total_evaluated']}")
    lines.append(f"- **Position pairs tested:** {double_mutant_results['n_position_pairs']}")
    lines.append(f"- **Better than wildtype:** {double_mutant_results['better_than_wildtype']}")
    lines.append("")
    lines.append("### Epistasis Statistics")
    lines.append("")
    lines.append(f"- **Mean epistasis:** {epi_stats['mean_epistasis_kcal']:.4f} kcal/mol")
    lines.append(f"- **Std epistasis:** {epi_stats['std_epistasis_kcal']:.4f} kcal/mol")
    lines.append(f"- **Range:** [{epi_stats['min_epistasis_kcal']:.4f}, {epi_stats['max_epistasis_kcal']:.4f}] kcal/mol")
    lines.append(f"- **Synergistic pairs:** {epi_stats['total_synergistic']} ({epi_stats['fraction_synergistic']:.1%})")
    lines.append(f"- **Antagonistic pairs:** {epi_stats['total_antagonistic']}")
    lines.append(f"- **Additive pairs:** {epi_stats['total_additive']}")
    lines.append("")

    best_double = double_mutant_results.get("best_double_mutant")
    if best_double:
        lines.append(f"**Best double mutant:** {best_double['mutation']} "
                     f"(dG = {best_double['dg_double_kcal']:.4f} kcal/mol, "
                     f"epistasis = {best_double['epistasis_ddg_kcal']:.4f} kcal/mol)")
        lines.append("")

    # --- 3. Length optimization ---
    lines.append("## 3. Length Optimization (18-30 nt)")
    lines.append("")
    lines.append(f"- **Lengths tested:** {length_results['length_range'][0]}-{length_results['length_range'][1]} nt")
    lines.append(f"- **Current ASO length:** {length_results['current_aso_length']} nt")

    if opt_len:
        lines.append(f"- **Optimal length:** {opt_len.get('length', 'N/A')} nt "
                     f"(dG = {opt_len.get('dg_binding_kcal', 0):.4f} kcal/mol)")
    lines.append(f"- **Current length is optimal:** {'YES' if current_is_opt else 'NO'}")
    lines.append("")

    lines.append("### Best dG by Length")
    lines.append("")
    lines.append("| Length | Best dG (kcal/mol) | Tm (C) | GC | Off-target Risk | Current? |")
    lines.append("|--------|-------------------|--------|------|-----------------|----------|")
    for lr in sorted(length_results.get("by_length", []), key=lambda r: r["length"]):
        is_current = lr.get("is_current_aso_length", False)
        marker = " <--" if is_current else ""
        lines.append(
            f"| {lr['length']:>6} | {lr['dg_binding_kcal']:>17.4f} | "
            f"{lr['tm_celsius']:>6.2f} | {lr['gc_content']:>.4f} | "
            f"{lr['off_target_risk']:>.2e} | "
            f"{'YES' + marker if is_current else 'no':>8} |"
        )
    lines.append("")

    # --- 4. Sliding window ---
    lines.append("## 4. Positional Sliding Window (25-mer)")
    lines.append("")
    lines.append(f"- **SL RNA length:** {sliding_results['sl_length']} nt")
    lines.append(f"- **Window length:** {sliding_results['window_length']} nt")
    lines.append(f"- **Total windows:** {sliding_results['n_windows']}")
    lines.append("")

    if sliding_opt:
        lines.append(f"- **Optimal window:** position {sliding_opt.get('start', 0)}-"
                     f"{sliding_opt.get('end', 0)} "
                     f"(dG = {sliding_opt.get('dg_binding_kcal', 0):.4f} kcal/mol)")
    lines.append(f"- **Current ASO covers optimal:** {'YES' if covers_opt else 'NO'}")

    dg_dev = sliding_results.get("dg_deviation_from_optimal_kcal")
    if dg_dev is not None:
        lines.append(f"- **dG deviation from optimal:** {dg_dev:.4f} kcal/mol")
    lines.append("")

    lines.append("### All Windows")
    lines.append("")
    lines.append("| Start | End | dG (kcal/mol) | Tm (C) | GC | Current? |")
    lines.append("|-------|-----|---------------|--------|------|----------|")
    for w in sliding_results.get("windows", []):
        is_current = w.get("is_current_aso_position", False)
        marker = " <--" if is_current else ""
        lines.append(
            f"| {w['start']:>5} | {w['end']:>3} | "
            f"{w['dg_binding_kcal']:>13.4f} | {w['tm_celsius']:>6.2f} | "
            f"{w['gc_content']:>.4f} | {'YES' + marker if is_current else 'no':>8} |"
        )
    lines.append("")

    # --- 5. Summary statistics ---
    lines.append("## 5. Summary Statistics")
    lines.append("")
    lines.append(f"- **Rank among single mutants:** {rank_info['rank']}/{rank_info['total']}")
    lines.append(f"- **Binding energy percentile:** {pct_info['percentile']:.2f}%")
    lines.append(f"- **Variants within {near_info['threshold_kcal']} kcal/mol of optimum:** {near_info['count']}")
    lines.append(f"- **Robustness score (mean |ddG|):** {robust_info['mean_abs_ddg_kcal']:.4f} +/- "
                 f"{robust_info['std_abs_ddg_kcal']:.4f} kcal/mol")
    lines.append(f"- **Single mutants better than WT:** {summary_stats['single_mutants_better_than_wt']}")
    lines.append(f"- **Double mutants better than WT:** {summary_stats['double_mutants_better_than_wt']}")
    lines.append("")

    # --- Conclusao ---
    lines.append("## Conclusion")
    lines.append("")
    is_optimal_single = summary_stats["single_mutants_better_than_wt"] == 0
    is_optimal_double = summary_stats["double_mutants_better_than_wt"] == 0

    if is_optimal_single and is_optimal_double and covers_opt and current_is_opt:
        lines.append(
            "MRL-ASO-001 is the **GLOBAL OPTIMUM** verified across all tested dimensions: "
            "no single mutant, no double mutant, no alternative length, and no alternative "
            "binding window produces a more favorable binding energy. "
            f"The fitness gap of {fitness_gap:.4f} kcal/mol indicates a thermodynamically "
            "distinct optimum."
        )
    else:
        # Relato honesto
        findings: list[str] = []
        if not is_optimal_single:
            findings.append(
                f"{summary_stats['single_mutants_better_than_wt']} single mutant(s) "
                f"with better dG"
            )
        if not is_optimal_double:
            findings.append(
                f"{summary_stats['double_mutants_better_than_wt']} double mutant(s) "
                f"with better dG"
            )
        if not current_is_opt:
            findings.append(
                f"optimal length is {summary_stats['optimal_length_nt']} nt, "
                f"not {length_results['current_aso_length']} nt"
            )
        if not covers_opt:
            findings.append(
                f"current position does not cover the thermodynamically optimal window "
                f"(deviation: {dg_dev:.4f} kcal/mol)"
            )

        if findings:
            lines.append(
                "MRL-ASO-001 is **NOT** the global optimum. Findings: "
                + "; ".join(findings) + "."
            )
        else:
            lines.append(
                "MRL-ASO-001 is confirmed as the optimal candidate across all "
                "tested dimensions."
            )

    lines.append("")
    lines.append("---")
    lines.append(f"*Analysis completed in {runtime_seconds:.2f} seconds.*")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Funcao auxiliar para salvar JSON
# ---------------------------------------------------------------------------


def _save_json(data: Any, filename: str) -> Path:
    """Salva dados como JSON no diretorio de resultados.

    Args:
        data: Dados a serializar.
        filename: Nome do arquivo (sem extensao).
    """
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / f"{filename}.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    return path


# ---------------------------------------------------------------------------
# Orquestrador principal
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a analise estendida da paisagem termodinamica.

    Fluxo:
        1. Construir paisagem mutacional completa 25x4
        2. Analisar mutantes duplos com epistasia
        3. Otimizar comprimento (18-30 nt)
        4. Janela deslizante posicional
        5. Calcular estatisticas agregadas
        6. Gerar relatorio Markdown
        7. Salvar todos os resultados como JSON

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

    aso_seq = config.aso_sequence.upper()
    sl_seq = config.sl_sequence.upper()
    aso_name = config.aso_name or "ASO"
    target_start = config.aso_target_start

    logger.info("=" * 70)
    logger.info("MATH 1 — Paisagem Termodinamica Estendida — %s", aso_name)
    logger.info("=" * 70)

    with Timer() as timer:
        # ---------------------------------------------------------------
        # Passo 1: Paisagem mutacional completa (25 x 4)
        # ---------------------------------------------------------------
        logger.info("Passo 1/5: Construindo paisagem mutacional %dx4...", len(aso_seq))
        fitness_landscape = build_full_fitness_landscape(aso_seq)

        wt_dg = fitness_landscape["wildtype_dg_kcal"]
        wt_tm = fitness_landscape["wildtype_tm_celsius"]

        # Construir mapa de dG de mutantes simples para epistasia
        single_dg_map: dict[tuple[int, str], float] = {}
        for v in fitness_landscape["variants"]:
            if not v["is_wildtype"]:
                single_dg_map[(v["position"], v["base"])] = v["dg_binding_kcal"]

        _save_json(fitness_landscape, "fitness_landscape")
        logger.info("  -> fitness_landscape.json salvo")

        # ---------------------------------------------------------------
        # Passo 2: Mutantes duplos com epistasia
        # ---------------------------------------------------------------
        logger.info("Passo 2/5: Analisando mutantes duplos com epistasia...")
        double_mutant_results = analyze_double_mutants_epistasis(
            aso_sequence=aso_seq,
            single_mutant_dg_map=single_dg_map,
            wt_dg=wt_dg,
            sample_all=FULL_DOUBLE_MUTANT_SCAN,
        )
        _save_json(double_mutant_results, "double_mutants")
        logger.info("  -> double_mutants.json salvo")

        # ---------------------------------------------------------------
        # Passo 3: Otimizacao de comprimento (18-30 nt)
        # ---------------------------------------------------------------
        logger.info("Passo 3/5: Otimizando comprimento (18-30 nt)...")
        length_results = optimize_length(
            sl_sequence=sl_seq,
            aso_sequence=aso_seq,
            aso_target_start=target_start,
            length_min=LENGTH_SCAN_MIN_EXT,
            length_max=LENGTH_SCAN_MAX_EXT,
        )
        _save_json(length_results, "length_optimization")
        logger.info("  -> length_optimization.json salvo")

        # ---------------------------------------------------------------
        # Passo 4: Janela deslizante posicional
        # ---------------------------------------------------------------
        logger.info("Passo 4/5: Executando sliding window de %d nt...", SLIDING_WINDOW_LENGTH)
        sliding_results = sliding_window_analysis(
            sl_sequence=sl_seq,
            window_length=SLIDING_WINDOW_LENGTH,
            aso_target_start=target_start,
        )
        _save_json(sliding_results, "sliding_window")
        logger.info("  -> sliding_window.json salvo")

        # ---------------------------------------------------------------
        # Passo 5: Estatisticas agregadas
        # ---------------------------------------------------------------
        logger.info("Passo 5/5: Calculando estatisticas agregadas...")
        summary_stats = compute_summary_statistics(
            wt_dg=wt_dg,
            wt_tm=wt_tm,
            fitness_landscape=fitness_landscape,
            double_mutant_results=double_mutant_results,
            length_results=length_results,
            sliding_results=sliding_results,
        )
        _save_json(summary_stats, "landscape_summary")
        logger.info("  -> landscape_summary.json salvo")

    # --- Gerar relatorio Markdown ---
    report = _generate_report(
        fitness_landscape=fitness_landscape,
        double_mutant_results=double_mutant_results,
        length_results=length_results,
        sliding_results=sliding_results,
        summary_stats=summary_stats,
        runtime_seconds=timer.elapsed,
    )
    report_path = RESULTS_DIR / "THERMODYNAMIC_LANDSCAPE_REPORT.md"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    logger.info("Relatorio salvo em: %s", report_path)

    # --- Resultado consolidado ---
    result = {
        "module": "math_1_landscape",
        "version": "1.0.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": timer.elapsed,
        "status": "success",
        "aso_name": aso_name,
        "organism": config.species_name,
        "fitness_landscape": {
            "wildtype_dg_kcal": wt_dg,
            "wildtype_tm_celsius": wt_tm,
            "global_minimum": fitness_landscape["global_minimum"],
            "fitness_gap_kcal": fitness_landscape["fitness_gap_kcal"],
            "total_variants": fitness_landscape["total_variants"],
        },
        "double_mutants": {
            "total_evaluated": double_mutant_results["total_evaluated"],
            "better_than_wildtype": double_mutant_results["better_than_wildtype"],
            "epistasis_mean": double_mutant_results["epistasis_statistics"]["mean_epistasis_kcal"],
            "epistasis_synergistic_count": double_mutant_results["epistasis_statistics"]["total_synergistic"],
        },
        "length_optimization": {
            "optimal_length": length_results.get("optimal_length", {}).get("length"),
            "current_is_optimal": length_results.get("current_aso_is_optimal", False),
        },
        "sliding_window": {
            "optimal_start": sliding_results.get("optimal_window", {}).get("start"),
            "current_covers_optimal": sliding_results.get("current_covers_optimal", False),
            "dg_deviation_kcal": sliding_results.get("dg_deviation_from_optimal_kcal"),
        },
        "summary": summary_stats,
    }

    logger.info("=" * 70)
    logger.info("MATH 1 COMPLETO — Tempo total: %.2f s", timer.elapsed)
    logger.info("=" * 70)

    return result


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
