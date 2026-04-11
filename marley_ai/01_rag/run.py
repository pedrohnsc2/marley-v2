"""Execucao do modulo 01_rag — RAG sobre literatura de leishmaniose.

Uso:
    python -m marley_ai.01_rag.run
"""

from __future__ import annotations

from typing import Any

from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result


def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Executa o pipeline de RAG.

    TODO: Implementar:
        1. Carregar corpus de documentos (PDFs de leishmaniose, patentes ASO)
        2. Aplicar chunking com overlap para preservar contexto biologico
        3. Gerar embeddings com modelo de linguagem (e5-large-v2)
        4. Construir indice FAISS para busca vetorial
        5. Validar com queries de teste (ex: "mecanismo de acao do SL RNA")
        6. Salvar indice e metadados para uso pelo modulo 11_scientist

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com resultados do RAG.
    """
    envelope = create_envelope("01_rag")

    with Timer() as timer:
        # TODO: Ingestao de documentos
        # Fontes: PubMed (leishmaniasis + ASO), TriTrypDB, WHO NTD reports
        # Formato: PDF -> texto via PyMuPDF ou pdfplumber

        # TODO: Chunking semantico
        # Preservar limites de paragrafo/secao quando possivel
        # chunk_size=512 tokens, overlap=64 tokens

        # TODO: Embedding
        # Modelo: intfloat/e5-large-v2 (1024-dim, bom para retrieval)
        # Prefixo de query: "query: " vs "passage: " (formato e5)

        # TODO: Indexacao
        # FAISS IndexFlatIP para prototipo (forca bruta, exato)
        # Migrar para IndexIVFFlat quando corpus > 10k chunks

        # TODO: Validacao
        # Queries de teste com respostas conhecidas
        # Metricas: recall@k, MRR, tempo de busca

        envelope["status"] = "stub"
        envelope["summary"]["conclusion"] = (
            "Modulo RAG em construcao. "
            "Indexara literatura de leishmaniose para consulta semantica."
        )
        envelope["summary"]["key_metrics"] = {
            "corpus_size": 0,
            "index_size": 0,
            "embedding_dim": 1024,
        }

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"[01_rag] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
