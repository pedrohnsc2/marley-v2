"""Orquestrador do suite aso_math — executa todos os modulos em sequencia.

Uso:
    python -m aso_math.run_all                 # todos os modulos + certificado
    python -m aso_math.run_all --module 01     # modulo especifico
    python -m aso_math.run_all --certificate   # apenas certificado (requer resultados)

Cada modulo e independente e pode ser executado separadamente:
    python -m aso_math.01_thermodynamic_landscape.run
    python -m aso_math.02_selectivity_proof.run
    ...
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from typing import Any

from core.logger import get_logger

logger = get_logger("aso_math_runner")

# Modulos em ordem de execucao
MODULES: list[tuple[str, str]] = [
    ("01", "aso_math.01_thermodynamic_landscape.run"),
    ("02", "aso_math.02_selectivity_proof.run"),
    ("03", "aso_math.03_evolutionary_conservation.run"),
    ("04", "aso_math.04_exhaustive_optimization.run"),
    ("05", "aso_math.05_resistance_model.run"),
]


def _run_module(
    module_id: str,
    module_path: str,
    config: Any = None,
) -> dict[str, Any]:
    """Importa e executa um modulo.

    Se config for fornecido e o modulo aceitar, passa como parametro.

    Returns:
        Dict com status, duracao e erros (se houver).
    """
    result: dict[str, Any] = {
        "module": module_id,
        "path": module_path,
        "status": "pending",
        "duration_seconds": 0.0,
        "error": "",
    }

    logger.info("=" * 60)
    logger.info("Executando modulo %s (%s)", module_id, module_path)
    logger.info("=" * 60)

    start = time.time()
    try:
        mod = importlib.import_module(module_path)
        # Passa config se o modulo aceitar (inspeciona assinatura)
        import inspect
        sig = inspect.signature(mod.main)
        if config is not None and "config" in sig.parameters:
            mod.main(config=config)
        else:
            mod.main()
        result["status"] = "success"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
        logger.error("Modulo %s falhou: %s", module_id, exc)
    finally:
        result["duration_seconds"] = round(time.time() - start, 2)

    logger.info(
        "Modulo %s: %s (%.2f s)",
        module_id, result["status"], result["duration_seconds"],
    )
    return result


def _run_certificate() -> dict[str, Any]:
    """Executa o certificado matematico final."""
    logger.info("=" * 60)
    logger.info("Gerando Certificado Matematico")
    logger.info("=" * 60)

    result: dict[str, Any] = {
        "module": "certificate",
        "status": "pending",
        "duration_seconds": 0.0,
        "error": "",
    }

    start = time.time()
    try:
        from aso_math.reports.math_certificate import generate_certificate
        cert = generate_certificate()
        result["status"] = "success"
        result["verdict"] = cert.get("verdict", "unknown")
        result["score"] = cert.get("overall_score", 0)
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
        logger.error("Certificado falhou: %s", exc)
    finally:
        result["duration_seconds"] = round(time.time() - start, 2)

    return result


def _print_summary(results: list[dict[str, Any]]) -> None:
    """Imprime resumo da execucao."""
    logger.info("")
    logger.info("=" * 60)
    logger.info("RESUMO DA EXECUCAO")
    logger.info("=" * 60)

    total_time = 0.0
    for r in results:
        total_time += r["duration_seconds"]
        status_icon = "OK" if r["status"] == "success" else "FALHOU"
        extra = ""
        if "verdict" in r:
            extra = f" [{r['verdict']}]"
        logger.info(
            "  %s  %-45s  %6.2f s%s",
            status_icon, r.get("path", r["module"]), r["duration_seconds"], extra,
        )

    logger.info("-" * 60)
    logger.info("  Tempo total: %.2f s", total_time)
    logger.info("=" * 60)


def main() -> None:
    """Entry point do orquestrador."""
    parser = argparse.ArgumentParser(
        description="Executa os modulos de validacao matematica do MRL-ASO-001.",
    )
    parser.add_argument(
        "--module",
        type=str,
        default=None,
        help="Executar apenas um modulo especifico (ex: 01, 02, 03, 04, 05).",
    )
    parser.add_argument(
        "--certificate",
        action="store_true",
        help="Executar apenas o certificado (requer resultados de modulos).",
    )
    parser.add_argument(
        "--organism",
        type=str,
        default=None,
        help="Organismo do registro (ex: trypanosoma_cruzi, brugia_malayi).",
    )
    args = parser.parse_args()

    # Construir TargetConfig se --organism fornecido
    config = None
    if args.organism:
        from aso_math.organisms import get_target_config
        config = get_target_config(args.organism)
        logger.info("=" * 60)
        logger.info("aso_math — Validacao para %s", config.species_name)
        logger.info("  SL RNA: %s (%d nt)", config.sl_sequence, config.sl_length)
        logger.info("  Doenca: %s", config.disease_name)
        logger.info("=" * 60)
    else:
        logger.info("=" * 60)
        logger.info("aso_math — Validacao Matematica do MRL-ASO-001")
        logger.info("=" * 60)

    results: list[dict[str, Any]] = []

    if args.certificate:
        # Apenas certificado
        cert_result = _run_certificate()
        results.append(cert_result)
    elif args.module:
        # Modulo especifico
        target = args.module.zfill(2)
        found = False
        for mid, mpath in MODULES:
            if mid == target:
                results.append(_run_module(mid, mpath))
                found = True
                break
        if not found:
            logger.error("Modulo '%s' nao encontrado. Disponiveis: %s",
                         args.module, [m[0] for m in MODULES])
            sys.exit(1)
    else:
        # Todos os modulos + certificado
        for mid, mpath in MODULES:
            result = _run_module(mid, mpath)
            results.append(result)
            if result["status"] == "failed":
                logger.warning("Modulo %s falhou. Continuando...", mid)

        # Certificado final
        cert_result = _run_certificate()
        results.append(cert_result)

    _print_summary(results)

    # Exit code
    failures = [r for r in results if r["status"] == "failed"]
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
