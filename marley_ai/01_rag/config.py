"""Configuracao do modulo 01_rag — RAG sobre literatura de leishmaniose.

Define queries PubMed, parametros de chunking/embedding e caminhos
de armazenamento. Segue o padrao frozen dataclass do projeto.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Caminhos do modulo
# ---------------------------------------------------------------------------

RAG_ROOT: Final[Path] = AI_ROOT / "01_rag"
RAG_RESULTS_DIR: Final[Path] = RAG_ROOT / "results"
PAPERS_INDEX_PATH: Final[Path] = RAG_RESULTS_DIR / "papers_index.json"
TFIDF_INDEX_PATH: Final[Path] = RAG_RESULTS_DIR / "tfidf_index.json"
QUERY_RESULTS_PATH: Final[Path] = RAG_RESULTS_DIR / "query_results.json"

# ---------------------------------------------------------------------------
# Configuracao de acesso ao PubMed (NCBI Entrez)
# ---------------------------------------------------------------------------

# Email obrigatorio para usar a API Entrez — substituir pelo email real
ENTREZ_EMAIL: Final[str] = "user@example.com"

# Queries tematicas para busca no PubMed
# Cobrem os eixos principais do projeto Marley:
#   - ASO / antisense contra Leishmania
#   - Spliced leader RNA (alvo terapeutico)
#   - Vacina canina contra leishmaniose
#   - Biologia molecular de L. infantum
PUBMED_QUERIES: Final[list[str]] = [
    "Leishmania infantum antisense oligonucleotide",
    "spliced leader RNA Leishmania",
    "canine leishmaniasis vaccine",
    "Leishmania infantum drug target",
    "antisense oligonucleotide trypanosomatid",
    "Leishmania RNA interference",
    "canine visceral leishmaniasis immunotherapy",
    "Leishmania infantum genomics proteomics",
]

# Numero maximo de artigos por query (total = queries * max_per_query)
MAX_PAPERS_PER_QUERY: Final[int] = 60
MAX_PAPERS_TOTAL: Final[int] = 500

# ---------------------------------------------------------------------------
# Parametros de indexacao TF-IDF
# ---------------------------------------------------------------------------

# Tamanho minimo de abstract para indexar (caracteres)
MIN_ABSTRACT_LENGTH: Final[int] = 100

# Termos de parada (stop words) minimos para TF-IDF em ingles cientifico
# Nao incluimos termos biologicos aqui — queremos que "Leishmania" tenha peso
STOP_WORDS: Final[frozenset[str]] = frozenset({
    "a", "an", "the", "and", "or", "but", "in", "on", "at", "to", "for",
    "of", "with", "by", "from", "is", "are", "was", "were", "be", "been",
    "being", "have", "has", "had", "do", "does", "did", "will", "would",
    "shall", "should", "may", "might", "can", "could", "not", "no", "nor",
    "so", "if", "then", "than", "that", "this", "these", "those", "it",
    "its", "as", "we", "our", "their", "they", "he", "she", "his", "her",
    "which", "who", "whom", "what", "when", "where", "how", "all", "each",
    "both", "few", "more", "most", "other", "some", "such", "only", "own",
    "same", "very", "also", "into", "over", "after", "before", "between",
    "under", "about", "up", "out", "through", "during", "here", "there",
    "while", "because", "although", "however", "therefore", "thus",
})

# ---------------------------------------------------------------------------
# Parametros de retrieval
# ---------------------------------------------------------------------------

# Numero padrao de resultados retornados por query
DEFAULT_TOP_K: Final[int] = 10

# Limiar minimo de similaridade coseno para considerar relevante
SIMILARITY_THRESHOLD: Final[float] = 0.05


# ---------------------------------------------------------------------------
# Dataclass de configuracao (estende AIModuleConfig)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RAGConfig(AIModuleConfig):
    """Configuracao especifica para o modulo de RAG.

    Campos adicionais controlam chunking, embedding e retrieval.
    """

    # --- Chunking ---
    chunk_size: int = 512          # Tokens por chunk
    chunk_overlap: int = 64        # Overlap entre chunks (preserva contexto)

    # --- Embedding (TF-IDF — nao requer modelo externo) ---
    embedding_model: str = "tfidf"  # Tier 1: TF-IDF puro com numpy
    embedding_dim: int = 0          # Definido dinamicamente pelo vocabulario

    # --- Retrieval ---
    top_k: int = DEFAULT_TOP_K
    rerank: bool = False            # Sem reranking no Tier 1
    similarity_threshold: float = SIMILARITY_THRESHOLD

    # --- Caminhos ---
    corpus_dir: Path = field(
        default_factory=lambda: RAG_ROOT / "corpus"
    )
    index_dir: Path = field(
        default_factory=lambda: RAG_ROOT / "index"
    )
    results_dir: Path = field(
        default_factory=lambda: RAG_RESULTS_DIR
    )
