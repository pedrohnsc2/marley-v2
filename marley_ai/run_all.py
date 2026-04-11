"""Orquestrador do track marley_ai — executa todos os modulos em sequencia.

Uso:
    python -m marley_ai.run_all                 # todos os modulos
    python -m marley_ai.run_all --module 03     # modulo especifico
    python -m marley_ai.run_all --list          # listar modulos disponiveis

Cada modulo e independente e pode ser executado separadamente:
    python -m marley_ai.01_rag.run
    python -m marley_ai.02_leish_kg.run
    ...

A ordem de execucao respeita dependencias entre modulos:
    01_rag, 02_leish_kg (sem dependencias — podem rodar em paralelo)
    03_leish_esm (sem dependencias)
    04_rna_fm (sem dependencias)
    05_rosettafold (depende de 04)
    06_evodiff (depende de 03)
    07_contrastive (depende de 03)
    08_rl_ppo (depende de 07)
    09_sae (depende de 03)
    10_digital_twin (depende de 07, 08)
    11_scientist (depende de todos)
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import sys
import time
from typing import Any

# Modulos em ordem de execucao (respeitando dependencias)
MODULES: list[tuple[str, str]] = [
    ("01", "marley_ai.01_rag.run"),
    ("02", "marley_ai.02_leish_kg.run"),
    ("03", "marley_ai.03_leish_esm.run"),
    ("04", "marley_ai.04_rna_fm.run"),
    ("05", "marley_ai.05_rosettafold.run"),
    ("06", "marley_ai.06_evodiff.run"),
    ("07", "marley_ai.07_contrastive.run"),
    ("08", "marley_ai.08_rl_ppo.run"),
    ("09", "marley_ai.09_sae.run"),
    ("10", "marley_ai.10_digital_twin.run"),
    ("11", "marley_ai.11_scientist.run"),
]

MODULE_NAMES: dict[str, str] = {
    "01": "RAG — Retrieval-Augmented Generation",
    "02": "Leish-KG — Knowledge Graph",
    "03": "Leish-ESM — Protein Embeddings",
    "04": "RNA-FM — RNA Foundation Models",
    "05": "RoseTTAFold — 3D Structure Prediction",
    "06": "EvoDiff — Evolutionary Diffusion",
    "07": "Contrastive — Epitope-MHC Learning",
    "08": "RL-PPO — Construct Optimization",
    "09": "SAE — Sparse Autoencoders",
    "10": "Digital Twin — Immune Simulation",
    "11": "AI Scientist — Autonomous Agent",
}


def _run_module(
    module_id: str,
    module_path: str,
    config: Any = None,
) -> dict[str, Any]:
    """Importa e executa um modulo.

    Segue o mesmo padrao do aso_math.run_all: importa dinamicamente,
    inspeciona assinatura do main(), e passa config se aceito.

    Returns:
        Dict com status, duracao e erros (se houver).
    """
    result: dict[str, Any] = {
        "module": module_id,
        "name": MODULE_NAMES.get(module_id, ""),
        "path": module_path,
        "status": "pending",
        "duration_seconds": 0.0,
        "error": "",
    }

    print("=" * 60)
    print(f"Executando modulo {module_id}: {MODULE_NAMES.get(module_id, module_path)}")
    print("=" * 60)

    start = time.time()
    try:
        mod = importlib.import_module(module_path)
        # Passa config se o modulo aceitar (inspeciona assinatura)
        sig = inspect.signature(mod.main)
        if config is not None and "config" in sig.parameters:
            mod.main(config=config)
        else:
            mod.main()
        result["status"] = "success"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
        print(f"  ERRO: Modulo {module_id} falhou: {exc}")
    finally:
        result["duration_seconds"] = round(time.time() - start, 2)

    print(f"  -> {result['status']} ({result['duration_seconds']:.2f} s)")
    print()

    return result


def _print_summary(results: list[dict[str, Any]]) -> None:
    """Imprime resumo da execucao."""
    print()
    print("=" * 60)
    print("RESUMO DA EXECUCAO — marley_ai")
    print("=" * 60)

    total_time = 0.0
    for r in results:
        total_time += r["duration_seconds"]
        status_icon = "OK" if r["status"] == "success" else "FALHOU"
        print(
            f"  {status_icon:8s}  {r['name']:45s}  {r['duration_seconds']:6.2f} s"
        )

    print("-" * 60)
    print(f"  Tempo total: {total_time:.2f} s")
    print(f"  Modulos: {len(results)} executados, "
          f"{sum(1 for r in results if r['status'] == 'success')} ok, "
          f"{sum(1 for r in results if r['status'] == 'failed')} falhas")
    print("=" * 60)


def main() -> None:
    """Entry point do orquestrador."""
    parser = argparse.ArgumentParser(
        description="Executa os modulos de IA/ML do projeto Marley.",
    )
    parser.add_argument(
        "--module",
        type=str,
        default=None,
        help="Executar apenas um modulo especifico (ex: 01, 03, 11).",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Listar modulos disponiveis e sair.",
    )
    args = parser.parse_args()

    if args.list:
        print("Modulos disponiveis no marley_ai:")
        print("-" * 60)
        for mid, mpath in MODULES:
            print(f"  {mid}  {MODULE_NAMES.get(mid, mpath)}")
        print("-" * 60)
        print(f"Total: {len(MODULES)} modulos")
        return

    print("=" * 60)
    print("marley_ai — Track 3: IA/ML para Leishmaniose Canina")
    print("=" * 60)
    print()

    results: list[dict[str, Any]] = []

    if args.module:
        # Modulo especifico
        target = args.module.zfill(2)
        found = False
        for mid, mpath in MODULES:
            if mid == target:
                results.append(_run_module(mid, mpath))
                found = True
                break
        if not found:
            print(
                f"Modulo '{args.module}' nao encontrado. "
                f"Disponiveis: {[m[0] for m in MODULES]}"
            )
            sys.exit(1)
    else:
        # Todos os modulos
        for mid, mpath in MODULES:
            result = _run_module(mid, mpath)
            results.append(result)
            if result["status"] == "failed":
                print(f"  Aviso: modulo {mid} falhou. Continuando...")

    _print_summary(results)

    # Exit code
    failures = [r for r in results if r["status"] == "failed"]
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
