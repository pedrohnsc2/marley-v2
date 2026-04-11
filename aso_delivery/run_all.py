"""Orquestrador do suite aso_delivery — executa todos os modulos em sequencia.

Uso:
    python -m aso_delivery.run_all                  # todos os modulos + certificado
    python -m aso_delivery.run_all --module a       # modulo especifico
    python -m aso_delivery.run_all --certificate    # apenas certificado (requer resultados)

Cada modulo e independente e pode ser executado separadamente:
    python -m aso_delivery.module_a_stability.run
    python -m aso_delivery.module_b_membrane.run
    ...
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from core.logger import get_logger

logger = get_logger("aso_delivery_runner")

# Diretorio de saida consolidado
_RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"

# Modulos em ordem de execucao — pipeline de delivery
MODULES: list[tuple[str, str]] = [
    ("A", "aso_delivery.module_a_stability.run"),
    ("B", "aso_delivery.module_b_membrane.run"),
    ("C", "aso_delivery.module_c_conjugate.run"),
    ("D", "aso_delivery.module_d_lnp.run"),
    ("E", "aso_delivery.module_e_admet.run"),
    ("F", "aso_delivery.module_f_immune_sde.run"),
]

# Mapeamento de IDs para facilitar busca por --module
_MODULE_MAP: dict[str, str] = {mid.lower(): mpath for mid, mpath in MODULES}


def _run_module(
    module_id: str,
    module_path: str,
) -> dict[str, Any]:
    """Importa e executa um modulo.

    Returns:
        Dict com status, duracao, dados retornados e erros (se houver).
    """
    result: dict[str, Any] = {
        "module": module_id,
        "path": module_path,
        "status": "pending",
        "duration_seconds": 0.0,
        "error": "",
        "data": None,
    }

    logger.info("=" * 60)
    logger.info("Executando modulo %s (%s)", module_id, module_path)
    logger.info("=" * 60)

    start = time.time()
    try:
        mod = importlib.import_module(module_path)
        data = mod.main()
        result["status"] = "success"
        result["data"] = data
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
    """Executa o certificado de delivery final."""
    logger.info("=" * 60)
    logger.info("Gerando Certificado de Delivery")
    logger.info("=" * 60)

    result: dict[str, Any] = {
        "module": "certificate",
        "status": "pending",
        "duration_seconds": 0.0,
        "error": "",
    }

    start = time.time()
    try:
        from aso_delivery.reports.delivery_certificate import generate_certificate
        cert = generate_certificate()
        result["status"] = "success"
        result["verdict"] = cert.get("overall_verdict", "unknown")
        result["score"] = cert.get("composite_score", 0)
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)
        logger.error("Certificado falhou: %s", exc)
    finally:
        result["duration_seconds"] = round(time.time() - start, 2)

    return result


def _save_consolidated_report(
    module_results: list[dict[str, Any]],
    total_time: float,
) -> Path:
    """Grava relatorio consolidado em JSON."""
    n_success = sum(1 for r in module_results if r["status"] == "success")
    n_total = len(module_results)

    report: dict[str, Any] = {
        "title": "MRL-ASO-001 Delivery Pipeline Report",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "total_runtime_seconds": total_time,
        "modules_executed": n_total,
        "modules_succeeded": n_success,
        "modules_failed": n_total - n_success,
        "all_passed": n_success == n_total,
        "module_results": [],
    }

    for r in module_results:
        entry: dict[str, Any] = {
            "module": r["module"],
            "path": r.get("path", ""),
            "status": r["status"],
            "duration_seconds": r["duration_seconds"],
        }
        if r.get("error"):
            entry["error"] = r["error"]
        # Incluir status e executive_summary de cada modulo (nao o data completo)
        data = r.get("data")
        if data and isinstance(data, dict):
            entry["module_status"] = data.get("status", "unknown")
            entry["executive_summary"] = data.get("executive_summary", "")
        report["module_results"].append(entry)

    _RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = _RESULTS_DIR / "delivery_report.json"
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)

    logger.info("Relatorio consolidado gravado em: %s", output_path)
    return output_path


def _print_summary(results: list[dict[str, Any]]) -> None:
    """Imprime resumo da execucao."""
    logger.info("")
    logger.info("=" * 60)
    logger.info("RESUMO DA EXECUCAO — aso_delivery")
    logger.info("=" * 60)

    total_time = 0.0
    for r in results:
        total_time += r["duration_seconds"]
        status_icon = "OK" if r["status"] == "success" else "FALHOU"
        extra = ""
        if "verdict" in r:
            extra = f" [{r['verdict']}]"
        if "score" in r:
            extra += f" (score: {r['score']})"
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
        description="Executa os modulos de delivery do MRL-ASO-001.",
    )
    parser.add_argument(
        "--module",
        type=str,
        default=None,
        help="Executar apenas um modulo especifico (ex: a, b, c, d, e, f).",
    )
    parser.add_argument(
        "--certificate",
        action="store_true",
        help="Executar apenas o certificado (requer resultados de modulos).",
    )
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("aso_delivery — Pipeline de Delivery do MRL-ASO-001")
    logger.info("=" * 60)

    results: list[dict[str, Any]] = []

    if args.certificate:
        # Apenas certificado
        cert_result = _run_certificate()
        results.append(cert_result)
    elif args.module:
        # Modulo especifico
        target = args.module.lower()
        if target in _MODULE_MAP:
            mid = target.upper()
            mpath = _MODULE_MAP[target]
            results.append(_run_module(mid, mpath))
        else:
            logger.error(
                "Modulo '%s' nao encontrado. Disponiveis: %s",
                args.module, [m[0] for m in MODULES],
            )
            sys.exit(1)
    else:
        # Todos os modulos + certificado
        for mid, mpath in MODULES:
            result = _run_module(mid, mpath)
            results.append(result)
            if result["status"] == "failed":
                logger.warning("Modulo %s falhou. Continuando...", mid)

        # Gravar relatorio consolidado
        total_time = sum(r["duration_seconds"] for r in results)
        _save_consolidated_report(results, total_time)

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
