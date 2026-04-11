"""01_rag — Retrieval-Augmented Generation sobre literatura de leishmaniose.

Indexa artigos cientificos do PubMed em um indice TF-IDF vetorial,
permitindo consultas semanticas que alimentam o agente cientifico
(modulo 11) e validam hipoteses contra a literatura existente.

Pipeline:
    1. Busca de abstracts no PubMed via Entrez (Biopython)
    2. Indexacao com TF-IDF (numpy puro, sem dependencias pesadas)
    3. Retrieval por similaridade coseno com citacoes PMID
    4. Geracao de respostas (template ou Claude API se disponivel)

Abordagem em camadas (graceful degradation):
    Tier 1: TF-IDF + coseno (sempre funciona — numpy)
    Tier 2: sentence-transformers (se instalado — melhor qualidade)
    Tier 3: Claude API para geracao natural (se anthropic instalado)
"""

from __future__ import annotations

from typing import Any

from marley_ai.registry import register

__version__ = "0.2.0"


@register("01_rag")
class RAGModule:
    """Modulo RAG para literatura de leishmaniose.

    Implementa o protocolo AIModule do marley_ai.registry,
    permitindo execucao pelo orquestrador central.
    """

    def __init__(self) -> None:
        self._config: Any = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao RAG."""
        self._config = config

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica dependencias: numpy e Biopython."""
        missing: list[str] = []

        try:
            import numpy  # noqa: F401
        except ImportError:
            missing.append("numpy")

        try:
            from Bio import Entrez  # noqa: F401
        except ImportError:
            missing.append("biopython")

        return {
            "valid": len(missing) == 0,
            "missing": missing,
        }

    def run(self) -> dict[str, Any]:
        """Executa pipeline RAG completo."""
        from .run import main
        return main(config=self._config)

    def get_dependencies(self) -> list[str]:
        """RAG nao depende de outros modulos — e o primeiro da cadeia."""
        return []
