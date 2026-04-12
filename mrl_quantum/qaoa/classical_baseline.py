"""Baseline classico para comparacao com QAOA.

Implementa:
    1. Busca exaustiva para subconjuntos pequenos (N <= 10)
    2. Heuristica gulosa (greedy) para N=25 completo
    3. Comparacao QAOA vs classico

A busca exaustiva garante o minimo global exato para N=10
(2^10 = 1024 configuracoes), servindo como ground truth para
validar a qualidade da solucao QAOA.
"""

from __future__ import annotations

import logging
import time
from itertools import product
from typing import Any

import numpy as np

from mrl_quantum.config import (
    DG_BASE,
    N_POSITIONS,
    SEED,
)
from mrl_quantum.qaoa.qubo_formulation import evaluate_bitstring

logger = logging.getLogger("mrl_quantum.qaoa.classical")


def exhaustive_search(
    n_positions: int = 10,
    dg_base: float | None = None,
    top_k: int = 10,
) -> dict[str, Any]:
    """Busca exaustiva: enumera todas as 2^n configuracoes.

    Limitado a n <= 20 por viabilidade computacional.

    Args:
        n_positions: Numero de posicoes a otimizar (default: 10).
        dg_base: dG base do duplex (default: DG_BASE).
        top_k: Numero de melhores configuracoes a retornar.

    Returns:
        Dicionario com minimo global, top-k, e estatisticas.

    Raises:
        ValueError: Se n_positions > 20 (inviavel computacionalmente).
    """
    if n_positions > 20:
        raise ValueError(
            f"Busca exaustiva inviavel para n={n_positions} "
            f"(2^{n_positions} = {2**n_positions} configuracoes). "
            f"Use n <= 20."
        )

    dg = dg_base if dg_base is not None else DG_BASE
    total_configs = 2 ** n_positions

    logger.info(
        "Busca exaustiva: n=%d, total=%d configuracoes", n_positions, total_configs
    )

    t0 = time.time()

    results: list[dict[str, Any]] = []
    for bits in product([0, 1], repeat=n_positions):
        bitstring = list(bits)
        eval_result = evaluate_bitstring(
            bitstring, n_positions=n_positions, dg_base=dg
        )
        results.append(eval_result)

    # Ordenar por dG total (menor = mais estavel)
    results.sort(key=lambda r: r["dg_total"])

    elapsed = round(time.time() - t0, 4)

    # Estatisticas da distribuicao
    all_dg = [r["dg_total"] for r in results]

    report = {
        "method": "exhaustive_search",
        "n_positions": n_positions,
        "total_configurations": total_configs,
        "runtime_seconds": elapsed,
        "dg_base": dg,
        "global_minimum": results[0],
        "global_maximum": results[-1],
        "top_k": results[:top_k],
        "statistics": {
            "mean_dg": round(float(np.mean(all_dg)), 4),
            "std_dg": round(float(np.std(all_dg)), 4),
            "min_dg": round(float(np.min(all_dg)), 4),
            "max_dg": round(float(np.max(all_dg)), 4),
            "median_dg": round(float(np.median(all_dg)), 4),
        },
    }

    logger.info(
        "Busca exaustiva concluida: dG_min=%.4f, dG_max=%.4f, tempo=%.4fs",
        report["global_minimum"]["dg_total"],
        report["global_maximum"]["dg_total"],
        elapsed,
    )

    return report


def greedy_search(
    n_positions: int | None = None,
    dg_base: float | None = None,
    seed: int = SEED,
) -> dict[str, Any]:
    """Heuristica gulosa: atribui Sp ou Rp posicao a posicao.

    Para cada posicao (em ordem), escolhe a configuracao que
    minimiza o dG parcial ate aquele ponto. Nao garante otimo global,
    mas e rapido e serve como baseline.

    Args:
        n_positions: Numero de posicoes (default: N_POSITIONS = 25).
        dg_base: dG base do duplex (default: DG_BASE).
        seed: Semente (para desempate deterministico).

    Returns:
        Dicionario com solucao greedy e metricas.
    """
    n = n_positions if n_positions is not None else N_POSITIONS
    dg = dg_base if dg_base is not None else DG_BASE

    np.random.seed(seed)

    logger.info("Heuristica gulosa: n=%d, seed=%d", n, seed)

    t0 = time.time()

    bitstring = [0] * n  # Comeca com tudo Rp

    for i in range(n):
        # Testar Rp (0) vs Sp (1) na posicao i
        bitstring[i] = 0
        eval_rp = evaluate_bitstring(bitstring, n_positions=n, dg_base=dg)

        bitstring[i] = 1
        eval_sp = evaluate_bitstring(bitstring, n_positions=n, dg_base=dg)

        # Escolher o que da menor dG (mais estavel)
        if eval_rp["dg_total"] <= eval_sp["dg_total"]:
            bitstring[i] = 0
        else:
            bitstring[i] = 1

    elapsed = round(time.time() - t0, 4)

    final_eval = evaluate_bitstring(bitstring, n_positions=n, dg_base=dg)

    report = {
        "method": "greedy_search",
        "n_positions": n,
        "runtime_seconds": elapsed,
        "dg_base": dg,
        "solution": final_eval,
        "seed": seed,
    }

    logger.info(
        "Heuristica gulosa concluida: dG=%.4f, Sp=%d/Rp=%d, tempo=%.4fs",
        final_eval["dg_total"],
        final_eval["sp_count"],
        final_eval["rp_count"],
        elapsed,
    )

    return report


def compare_methods(
    qaoa_result: dict[str, Any] | None = None,
    n_exhaustive: int = 10,
    n_positions_full: int | None = None,
    dg_base: float | None = None,
    seed: int = SEED,
) -> dict[str, Any]:
    """Compara QAOA com baselines classicos.

    Args:
        qaoa_result: Resultado do QAOA (da run_qaoa ou run_qaoa_single).
        n_exhaustive: Numero de posicoes para busca exaustiva.
        n_positions_full: Numero total de posicoes para greedy.
        dg_base: dG base.
        seed: Semente.

    Returns:
        Dicionario com comparacao entre metodos.
    """
    n_full = n_positions_full if n_positions_full is not None else N_POSITIONS
    dg = dg_base if dg_base is not None else DG_BASE

    logger.info(
        "Comparando metodos: exaustivo(n=%d), greedy(n=%d)",
        n_exhaustive, n_full,
    )

    # 1. Busca exaustiva (N=10)
    exhaustive = exhaustive_search(
        n_positions=n_exhaustive, dg_base=dg, top_k=10
    )

    # 2. Greedy para N completo
    greedy_full = greedy_search(
        n_positions=n_full, dg_base=dg, seed=seed
    )

    # 3. Greedy para N=10 (comparacao justa com exaustivo)
    greedy_small = greedy_search(
        n_positions=n_exhaustive, dg_base=dg, seed=seed
    )

    comparison = {
        "exhaustive_n10": {
            "dg_optimal": exhaustive["global_minimum"]["dg_total"],
            "configuration": exhaustive["global_minimum"]["configuration"],
            "sp_count": exhaustive["global_minimum"]["sp_count"],
        },
        "greedy_n10": {
            "dg_found": greedy_small["solution"]["dg_total"],
            "configuration": greedy_small["solution"]["configuration"],
            "sp_count": greedy_small["solution"]["sp_count"],
            "gap_vs_exhaustive": round(
                greedy_small["solution"]["dg_total"]
                - exhaustive["global_minimum"]["dg_total"],
                4,
            ),
        },
        "greedy_n25": {
            "dg_found": greedy_full["solution"]["dg_total"],
            "configuration": greedy_full["solution"]["configuration"],
            "sp_count": greedy_full["solution"]["sp_count"],
        },
    }

    # Se temos resultado QAOA, adicionar comparacao
    if qaoa_result and "best" in qaoa_result:
        qaoa_best = qaoa_result["best"]
        n_qaoa = qaoa_result.get("n_positions", n_full)

        comparison["qaoa"] = {
            "dg_found": qaoa_best["dg_total"],
            "configuration": qaoa_best["configuration"],
            "sp_count": qaoa_best["sp_count"],
            "n_positions": n_qaoa,
            "n_layers": qaoa_result.get("n_layers", "?"),
            "optimizer": qaoa_result.get("optimizer", "?"),
        }

        # Comparacao QAOA vs exaustivo (se mesmo N)
        if n_qaoa == n_exhaustive:
            gap = round(
                qaoa_best["dg_total"]
                - exhaustive["global_minimum"]["dg_total"],
                4,
            )
            pct = round(
                abs(gap) / abs(exhaustive["global_minimum"]["dg_total"]) * 100,
                2,
            ) if exhaustive["global_minimum"]["dg_total"] != 0 else 0.0

            comparison["qaoa_vs_exhaustive"] = {
                "gap_kcal": gap,
                "gap_percent": pct,
                "qaoa_found_optimal": gap == 0.0,
            }

        # Comparacao QAOA vs greedy (mesmo N)
        if n_qaoa == n_full:
            gap_greedy = round(
                qaoa_best["dg_total"]
                - greedy_full["solution"]["dg_total"],
                4,
            )
            comparison["qaoa_vs_greedy"] = {
                "gap_kcal": gap_greedy,
                "qaoa_better_than_greedy": gap_greedy < 0,
            }

    report = {
        "method": "comparison",
        "comparison": comparison,
        "exhaustive_full": exhaustive,
        "greedy_small": greedy_small,
        "greedy_full": greedy_full,
    }

    # Log do sumario
    logger.info(
        "Comparacao concluida: exaustivo_n10=%.4f, greedy_n10=%.4f, greedy_n25=%.4f",
        exhaustive["global_minimum"]["dg_total"],
        greedy_small["solution"]["dg_total"],
        greedy_full["solution"]["dg_total"],
    )

    return report
