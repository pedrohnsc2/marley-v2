"""Utilidades para o envelope JSON padronizado de resultados AI.

Estende o padrao do aso_math.envelope para incluir campos especificos
de modulos de ML: metricas de treinamento, device utilizado, e
referencias a artefatos (embeddings, checkpoints, grafos).

Todo modulo do marley_ai grava resultados usando este envelope,
garantindo consistencia e rastreabilidade.
"""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from marley_ai.config import RESULTS_DIR


def create_envelope(module_name: str, version: str = "0.1.0") -> dict[str, Any]:
    """Cria um envelope JSON vazio para um modulo AI.

    O envelope contem metadados padronizados que todo modulo deve preencher.
    Campos adicionais (vs aso_math): device, artifacts, metrics.

    Args:
        module_name: Identificador do modulo (ex: "03_leish_esm").
        version: Versao do modulo.

    Returns:
        Dicionario com a estrutura envelope pronta para preencher.
    """
    return {
        "module": module_name,
        "version": version,
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": 0.0,
        "device": "",
        "status": "pending",
        "warnings": [],
        "summary": {
            "conclusion": "",
            "key_metrics": {},
        },
        "artifacts": [],    # Lista de caminhos para artefatos (embeddings, modelos)
        "metrics": {},      # Metricas de treinamento/inferencia (loss, accuracy, etc.)
        "dependencies": [], # Modulos que forneceram inputs para este
        "data": {},
    }


def write_result(
    envelope: dict[str, Any],
    module_name: str | None = None,
    output_dir: Path | None = None,
) -> Path:
    """Grava o envelope de resultados como JSON.

    Sempre grava, mesmo em caso de falha — um JSON com status 'failed'
    e infinitamente mais util que um traceback sem output.

    Args:
        envelope: Dicionario envelope preenchido.
        module_name: Nome do arquivo (sem extensao). Se None, usa envelope['module'].
        output_dir: Diretorio de saida. Se None, usa RESULTS_DIR.

    Returns:
        Caminho do arquivo JSON gravado.
    """
    name = module_name or envelope.get("module", "unknown")
    target_dir = output_dir or RESULTS_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    output_path = target_dir / f"{name}.json"

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(envelope, fh, indent=2, ensure_ascii=False)

    return output_path


def load_result(module_name: str, results_dir: Path | None = None) -> dict[str, Any]:
    """Carrega um envelope de resultados de um modulo.

    Util para modulos que dependem de resultados de modulos anteriores.

    Args:
        module_name: Nome do modulo (ex: "03_leish_esm").
        results_dir: Diretorio onde buscar. Se None, usa RESULTS_DIR.

    Returns:
        Dicionario com o envelope de resultados.

    Raises:
        FileNotFoundError: Se o arquivo nao existir.
    """
    target_dir = results_dir or RESULTS_DIR
    path = target_dir / f"{module_name}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"Resultados de '{module_name}' nao encontrados em {path}. "
            f"Execute o modulo antes de usar seus resultados."
        )
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class Timer:
    """Context manager para medir tempo de execucao de um modulo."""

    def __init__(self) -> None:
        self.start: float = 0.0
        self.elapsed: float = 0.0

    def __enter__(self) -> Timer:
        self.start = time.time()
        return self

    def __exit__(self, *args: Any) -> None:
        self.elapsed = round(time.time() - self.start, 2)
