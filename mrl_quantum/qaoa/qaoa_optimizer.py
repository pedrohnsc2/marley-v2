"""QAOA para otimizacao de estereoisomeros fosforotioato do MRL-ASO-001.

Executa QAOA com p=1,2,3 camadas usando COBYLA e SPSA como otimizadores
classicos. Extrai top-10 configuracoes e compara com baseline classico.

Fluxo:
    1. Formular QUBO (qubo_formulation.py)
    2. Obter sampler (backends.py)
    3. Executar QAOA via MinimumEigenOptimizer
    4. Extrair amostras e converter para configuracoes Rp/Sp
    5. Gerar relatorio JSON com metadados

Referencia: Farhi E et al. (2014) arXiv:1411.4028 — QAOA original
"""

from __future__ import annotations

import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA, SPSA
from qiskit_algorithms.utils import algorithm_globals
from qiskit_optimization.algorithms import MinimumEigenOptimizer

from mrl_quantum.config import (
    ASO_SEQUENCE,
    DG_BASE,
    N_POSITIONS,
    QAOA_MAX_ITER,
    QAOA_P_LAYERS,
    QAOA_RESULTS_DIR,
    SEED,
)
from mrl_quantum.qaoa.backends import get_sampler
from mrl_quantum.qaoa.qubo_formulation import (
    bitstring_to_config,
    build_qubo,
    evaluate_bitstring,
    get_aso_position_label,
)

logger = logging.getLogger("mrl_quantum.qaoa.optimizer")

# Limite pratico de qubits para simulacao statevector
# 25 qubits = 2^25 amplitudes -> matriz expm estoura memoria no scipy
# 16 qubits e o maximo seguro para StatevectorSampler em hardware tipico
# 10 qubits no modo light (~30s vs ~10min com 16)
MAX_QAOA_QUBITS: int = 16


def _get_optimizer(name: str, max_iter: int, seed: int) -> Any:
    """Retorna instancia do otimizador classico pelo nome.

    Args:
        name: Nome do otimizador ('COBYLA' ou 'SPSA').
        max_iter: Numero maximo de iteracoes.
        seed: Semente para SPSA (COBYLA e deterministic por natureza).

    Returns:
        Instancia do otimizador configurado.
    """
    if name.upper() == "COBYLA":
        return COBYLA(maxiter=max_iter)
    if name.upper() == "SPSA":
        return SPSA(maxiter=max_iter)
    raise ValueError(f"Otimizador desconhecido: {name}. Use 'COBYLA' ou 'SPSA'.")


def run_qaoa_single(
    n_layers: int = 2,
    optimizer_name: str = "COBYLA",
    backend_name: str = "aer_statevector",
    max_iter: int = QAOA_MAX_ITER,
    n_positions: int | None = None,
    seed: int = SEED,
) -> dict[str, Any]:
    """Executa uma rodada de QAOA com parametros especificos.

    Args:
        n_layers: Numero de camadas QAOA (p).
        optimizer_name: Nome do otimizador classico.
        backend_name: Nome do backend quantico.
        max_iter: Iteracoes maximas do otimizador.
        n_positions: Numero de posicoes (default: N_POSITIONS).
        seed: Semente para reproducibilidade.

    Returns:
        Dicionario com resultados da otimizacao.
    """
    n = n_positions if n_positions is not None else N_POSITIONS

    # Configurar semente global
    algorithm_globals.random_seed = seed
    np.random.seed(seed)

    logger.info(
        "Iniciando QAOA: p=%d, otimizador=%s, backend=%s, n=%d, seed=%d",
        n_layers, optimizer_name, backend_name, n, seed,
    )

    t0 = time.time()

    # 1. Formular QUBO
    qp = build_qubo(n_positions=n)

    # 2. Obter sampler
    sampler, backend_used = get_sampler(backend_name=backend_name, seed=seed)

    # 3. Configurar QAOA
    optimizer = _get_optimizer(optimizer_name, max_iter, seed)
    qaoa = QAOA(
        sampler=sampler,
        optimizer=optimizer,
        reps=n_layers,
    )

    # 4. Executar otimizacao
    min_eigen_optimizer = MinimumEigenOptimizer(qaoa)
    result = min_eigen_optimizer.solve(qp)

    elapsed = round(time.time() - t0, 2)

    # 5. Extrair amostras (top-10 por fval)
    samples = []
    if hasattr(result, "samples") and result.samples:
        # Ordenar por fval (menor = melhor, pois minimizamos dG)
        sorted_samples = sorted(result.samples, key=lambda s: s.fval)
        for s in sorted_samples[:10]:
            bitstring = [int(xi) for xi in s.x]
            eval_result = evaluate_bitstring(bitstring, n_positions=n)
            samples.append({
                "bitstring": bitstring,
                "fval_qubo": round(float(s.fval), 4),
                "dg_total": eval_result["dg_total"],
                "probability": round(float(s.probability), 6),
                "configuration": eval_result["configuration"],
                "sp_count": eval_result["sp_count"],
                "rp_count": eval_result["rp_count"],
                "rp_triples": eval_result["rp_triples"],
            })

    # 6. Resultado principal
    best_bitstring = [int(xi) for xi in result.x]
    best_eval = evaluate_bitstring(best_bitstring, n_positions=n)

    run_result = {
        "n_layers": n_layers,
        "optimizer": optimizer_name,
        "backend_used": backend_used,
        "n_positions": n,
        "seed": seed,
        "max_iter": max_iter,
        "runtime_seconds": elapsed,
        "status": str(result.status),
        "best": {
            "bitstring": best_bitstring,
            "fval_qubo": round(float(result.fval), 4),
            **best_eval,
        },
        "top_10_samples": samples,
    }

    logger.info(
        "QAOA concluido: p=%d, dG_best=%.4f kcal/mol, Sp=%d/Rp=%d, tempo=%.2fs",
        n_layers,
        best_eval["dg_total"],
        best_eval["sp_count"],
        best_eval["rp_count"],
        elapsed,
    )

    return run_result


def run_qaoa(
    n_layers: int | None = None,
    optimizer_name: str | None = None,
    backend_name: str = "aer_statevector",
    max_iter: int = QAOA_MAX_ITER,
    n_positions: int | None = None,
    seed: int = SEED,
    save_results: bool = True,
    max_qubits: int | None = None,
) -> dict[str, Any]:
    """Executa QAOA com varredura de parametros e gera relatorio completo.

    Se n_positions > MAX_QAOA_QUBITS, usa abordagem hibrida:
    QAOA nas primeiras MAX_QAOA_QUBITS posicoes + greedy no restante.
    Documenta honestamente a limitacao.

    Se n_layers e optimizer_name nao forem especificados, executa todas
    as combinacoes de p=[1,2,3] x optimizer=[COBYLA,SPSA].

    Args:
        n_layers: Camadas QAOA (None = varrer QAOA_P_LAYERS).
        optimizer_name: Otimizador (None = varrer QAOA_OPTIMIZERS).
        backend_name: Backend quantico.
        max_iter: Iteracoes maximas.
        n_positions: Numero de posicoes.
        seed: Semente.
        save_results: Se True, salva JSON em QAOA_RESULTS_DIR.

    Returns:
        Dicionario com todos os resultados e sumario.
    """
    from mrl_quantum.qaoa.classical_baseline import (
        compare_methods,
        exhaustive_search,
        greedy_search,
    )

    t_total = time.time()

    n = n_positions if n_positions is not None else N_POSITIONS
    effective_max = max_qubits if max_qubits is not None else MAX_QAOA_QUBITS

    # Para N > effective_max, usamos abordagem hibrida
    # (StatevectorSampler com 25 qubits gera matrizes 2^25 x 2^25
    # que excedem a capacidade do scipy.sparse.linalg.expm)
    n_qaoa = min(n, effective_max)
    is_hybrid = n > effective_max

    if is_hybrid:
        logger.info(
            "N=%d excede limite de qubits (%d). Usando abordagem hibrida: "
            "QAOA(n=%d) + greedy para posicoes restantes.",
            n, MAX_QAOA_QUBITS, n_qaoa,
        )

    # Determinar combinacoes
    layers_list = [n_layers] if n_layers is not None else QAOA_P_LAYERS
    opt_list = [optimizer_name] if optimizer_name is not None else ["COBYLA"]

    logger.info(
        "Iniciando varredura QAOA: camadas=%s, otimizadores=%s, n_qaoa=%d",
        layers_list, opt_list, n_qaoa,
    )

    runs: list[dict[str, Any]] = []

    for p in layers_list:
        for opt in opt_list:
            try:
                run_result = run_qaoa_single(
                    n_layers=p,
                    optimizer_name=opt,
                    backend_name=backend_name,
                    max_iter=max_iter,
                    n_positions=n_qaoa,
                    seed=seed,
                )
                runs.append(run_result)
            except Exception as exc:
                logger.error(
                    "Falha no QAOA p=%d, opt=%s: %s", p, opt, exc
                )
                runs.append({
                    "n_layers": p,
                    "optimizer": opt,
                    "status": "FAILED",
                    "error": str(exc),
                })

    # Encontrar melhor resultado QAOA
    successful_runs = [r for r in runs if r.get("status") != "FAILED"]
    best_run = None
    if successful_runs:
        best_run = min(
            successful_runs,
            key=lambda r: r["best"]["dg_total"],
        )

    # Para abordagem hibrida: estender solucao QAOA com greedy
    hybrid_result = None
    if is_hybrid and best_run:
        qaoa_bits = best_run["best"]["bitstring"]
        # Estender com greedy para as posicoes restantes
        full_bits = list(qaoa_bits) + [0] * (n - n_qaoa)
        # Greedy nas posicoes n_qaoa..n-1
        for i in range(n_qaoa, n):
            full_bits[i] = 0
            eval_rp = evaluate_bitstring(full_bits, n_positions=n)
            full_bits[i] = 1
            eval_sp = evaluate_bitstring(full_bits, n_positions=n)
            full_bits[i] = 0 if eval_rp["dg_total"] <= eval_sp["dg_total"] else 1

        hybrid_result = evaluate_bitstring(full_bits, n_positions=n)
        hybrid_result["method"] = "hybrid_qaoa_greedy"
        hybrid_result["qaoa_positions"] = n_qaoa
        hybrid_result["greedy_positions"] = n - n_qaoa

    # Comparacao com baselines classicos
    n_exhaustive = min(n_qaoa, 10)
    comparison = compare_methods(
        qaoa_result=best_run,
        n_exhaustive=n_exhaustive,
        n_positions_full=n,
        seed=seed,
    )

    # Greedy puro para N completo (baseline)
    greedy_full = greedy_search(n_positions=n, seed=seed)

    total_elapsed = round(time.time() - t_total, 2)

    # Determinar o melhor resultado final
    final_best = None
    final_method = None
    if hybrid_result:
        final_best = hybrid_result
        final_method = "hybrid_qaoa_greedy"
    elif best_run:
        final_best = best_run["best"]
        final_method = f"qaoa_p{best_run['n_layers']}_{best_run['optimizer']}"

    # Verificar se greedy puro e melhor (documentar honestamente)
    greedy_beats_qaoa = False
    if final_best and greedy_full:
        if greedy_full["solution"]["dg_total"] < final_best["dg_total"]:
            greedy_beats_qaoa = True
            logger.info(
                "NOTA: greedy puro (dG=%.4f) superou QAOA (dG=%.4f). "
                "Isso e esperado para problemas onde a estrutura do QUBO "
                "favorece solucoes gulosas.",
                greedy_full["solution"]["dg_total"],
                final_best["dg_total"],
            )

    # Sumario textual
    summary_lines = [
        f"QAOA PS Stereoisomer Optimization — MRL-ASO-001",
        f"Posicoes: {n} (QAOA em {n_qaoa}), Seed: {seed}, Backend: {backend_name}",
        f"Combinacoes testadas: {len(runs)}",
        f"Sucesso: {len(successful_runs)}/{len(runs)}",
        f"Tempo total: {total_elapsed}s",
    ]

    if is_hybrid:
        summary_lines.append(
            f"Abordagem hibrida: QAOA({n_qaoa} qubits) + greedy({n - n_qaoa} posicoes)"
        )

    if final_best:
        summary_lines.extend([
            f"",
            f"--- Melhor resultado ({final_method}) ---",
            f"dG total: {final_best['dg_total']:.4f} kcal/mol",
            f"Configuracao: Sp={final_best['sp_count']}, Rp={final_best['rp_count']}",
            f"Rp triples: {final_best.get('rp_triples', '?')}",
            f"dG base (sem PS): {DG_BASE:.2f} kcal/mol",
            f"Melhoria PS: {final_best['dg_total'] - DG_BASE:.4f} kcal/mol",
        ])

    if greedy_full:
        gsol = greedy_full["solution"]
        summary_lines.extend([
            f"",
            f"--- Baseline greedy (N={n}) ---",
            f"dG greedy: {gsol['dg_total']:.4f} kcal/mol",
            f"Sp={gsol['sp_count']}, Rp={gsol['rp_count']}",
        ])

    if greedy_beats_qaoa:
        summary_lines.append(
            f"NOTA: Greedy superou QAOA por "
            f"{final_best['dg_total'] - greedy_full['solution']['dg_total']:.4f} kcal/mol"
        )

    # Mapa posicao-a-posicao do melhor resultado
    if final_best and n == N_POSITIONS and "bitstring" in final_best:
        summary_lines.append("")
        summary_lines.append("Mapa de estereoisomeros (melhor configuracao):")
        for i, xi in enumerate(final_best["bitstring"]):
            label = get_aso_position_label(i)
            stereo = "Sp" if xi == 1 else "Rp"
            summary_lines.append(f"  {label}: {stereo}")

    summary_text = "\n".join(summary_lines)

    report = {
        "module": "mrl_quantum.qaoa",
        "version": "0.1.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": total_elapsed,
        "seed": seed,
        "backend_requested": backend_name,
        "aso_sequence": ASO_SEQUENCE,
        "n_positions": n,
        "n_qaoa_qubits": n_qaoa,
        "is_hybrid": is_hybrid,
        "dg_base": DG_BASE,
        "parameters": {
            "layers_tested": layers_list,
            "optimizers_tested": opt_list,
            "max_iter": max_iter,
            "max_qaoa_qubits": MAX_QAOA_QUBITS,
        },
        "runs": runs,
        "best_qaoa": best_run["best"] if best_run else None,
        "best_qaoa_config": {
            "n_layers": best_run["n_layers"],
            "optimizer": best_run["optimizer"],
        } if best_run else None,
        "hybrid_result": hybrid_result,
        "greedy_full": greedy_full["solution"] if greedy_full else None,
        "greedy_beats_qaoa": greedy_beats_qaoa,
        "comparison": comparison["comparison"],
        "final_best": final_best,
        "final_method": final_method,
        "summary": summary_text,
    }

    # Salvar resultados
    if save_results:
        QAOA_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        output_path = QAOA_RESULTS_DIR / "qaoa_ps_optimization.json"
        with open(output_path, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2, ensure_ascii=False)
        logger.info("Resultados salvos em: %s", output_path)

    return report
