"""Utilidades para o envelope JSON padronizado de resultados.

Todo modulo do aso_math grava resultados usando o mesmo formato envelope,
permitindo que o certificado final leia e agregue de forma consistente.
"""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from aso_math.config import RESULTS_DIR


def create_envelope(module_name: str) -> dict[str, Any]:
    """Cria um envelope JSON vazio para um modulo.

    O envelope contem metadados padronizados que todo modulo deve preencher.
    O campo 'data' deve ser preenchido pelo modulo com seus resultados.

    Args:
        module_name: Identificador do modulo (ex: "01_thermodynamic_landscape").

    Returns:
        Dicionario com a estrutura envelope pronta para preencher.
    """
    return {
        "module": module_name,
        "version": "1.0.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": 0.0,
        "tier_used": 0,
        "status": "pending",
        "warnings": [],
        "summary": {
            "conclusion": "",
            "key_metrics": {},
        },
        "data": {},
    }


def write_result(envelope: dict[str, Any], module_name: str | None = None) -> Path:
    """Grava o envelope de resultados como JSON.

    Sempre grava, mesmo em caso de falha — um JSON com status 'failed'
    e infinitamente mais util que um traceback sem output.

    Args:
        envelope: Dicionario envelope preenchido.
        module_name: Nome do arquivo (sem extensao). Se None, usa envelope['module'].

    Returns:
        Caminho do arquivo JSON gravado.
    """
    name = module_name or envelope.get("module", "unknown")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = RESULTS_DIR / f"{name}.json"

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(envelope, fh, indent=2, ensure_ascii=False)

    return output_path


class Timer:
    """Context manager para medir tempo de execucao de um modulo."""

    def __init__(self) -> None:
        self.start: float = 0.0
        self.elapsed: float = 0.0

    def __enter__(self) -> "Timer":
        self.start = time.time()
        return self

    def __exit__(self, *args: Any) -> None:
        self.elapsed = round(time.time() - self.start, 2)
