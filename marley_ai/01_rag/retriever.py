"""Retriever vetorial para busca semantica em artigos de leishmaniose.

Recebe uma query em linguagem natural, computa seu vetor TF-IDF
usando o mesmo vocabulario do indice, e retorna os artigos mais
similares via similaridade coseno.

Resultados incluem PMID, score de relevancia e trecho do abstract,
formatados para uso direto pelo gerador de respostas.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from .config import (
    DEFAULT_TOP_K,
    PAPERS_INDEX_PATH,
    SIMILARITY_THRESHOLD,
    TFIDF_INDEX_PATH,
)
from .indexer import (
    compute_tf,
    normalize_l2,
    tokenize,
)


# ---------------------------------------------------------------------------
# Resultado de busca
# ---------------------------------------------------------------------------

def _truncate_abstract(abstract: str, max_chars: int = 300) -> str:
    """Trunca abstract preservando palavras completas.

    Corta na ultima palavra completa antes do limite e adiciona '...'
    se o texto foi truncado.
    """
    if len(abstract) <= max_chars:
        return abstract
    truncated = abstract[:max_chars].rsplit(" ", 1)[0]
    return truncated + "..."


def format_citation(paper: dict[str, Any]) -> str:
    """Formata citacao curta no estilo PMID para um artigo.

    Formato: "Autor et al. (Ano) Journal. PMID: 12345678"
    """
    authors = paper.get("authors", [])
    if len(authors) >= 2:
        author_str = f"{authors[0]} et al."
    elif authors:
        author_str = authors[0]
    else:
        author_str = "Unknown"

    year = paper.get("year", "n.d.")
    journal = paper.get("journal", "")
    pmid = paper.get("pmid", "")

    parts = [author_str, f"({year})"]
    if journal:
        parts.append(f"{journal}.")
    parts.append(f"PMID: {pmid}")

    return " ".join(parts)


# ---------------------------------------------------------------------------
# Classe Retriever
# ---------------------------------------------------------------------------

class PaperRetriever:
    """Busca vetorial em artigos indexados usando TF-IDF + coseno.

    Carrega o indice do disco e mantem em memoria para queries rapidas.
    Suporta multiplas queries sem recarregar o indice.
    """

    def __init__(
        self,
        papers: list[dict[str, Any]] | None = None,
        tfidf_index: dict[str, Any] | None = None,
    ) -> None:
        """Inicializa retriever com artigos e indice TF-IDF.

        Se nao fornecidos, carrega do disco (gerados pelo indexer).

        Args:
            papers: Lista de artigos com metadados.
            tfidf_index: Indice TF-IDF com vocabulario, IDF e matriz.

        Raises:
            FileNotFoundError: Se indice nao existe no disco.
        """
        if papers is not None and tfidf_index is not None:
            self._papers = papers
            self._vocab = tfidf_index["vocabulary"]
            self._idf = np.array(tfidf_index["idf"], dtype=np.float64)
            self._matrix = np.array(tfidf_index["tfidf_matrix"], dtype=np.float64)
        else:
            self._load_from_disk()

        self._n_docs = len(self._papers)

    def _load_from_disk(self) -> None:
        """Carrega artigos e indice TF-IDF do disco."""
        if not PAPERS_INDEX_PATH.exists():
            raise FileNotFoundError(
                f"Indice de artigos nao encontrado em {PAPERS_INDEX_PATH}. "
                f"Execute o indexer primeiro: python -m marley_ai.01_rag.run"
            )
        if not TFIDF_INDEX_PATH.exists():
            raise FileNotFoundError(
                f"Indice TF-IDF nao encontrado em {TFIDF_INDEX_PATH}. "
                f"Execute o indexer primeiro: python -m marley_ai.01_rag.run"
            )

        with open(PAPERS_INDEX_PATH, encoding="utf-8") as f:
            self._papers = json.load(f)

        with open(TFIDF_INDEX_PATH, encoding="utf-8") as f:
            index = json.load(f)

        self._vocab = index["vocabulary"]
        self._idf = np.array(index["idf"], dtype=np.float64)
        self._matrix = np.array(index["tfidf_matrix"], dtype=np.float64)

    def _query_vector(self, query: str) -> np.ndarray:
        """Computa vetor TF-IDF normalizado para uma query.

        Usa o mesmo vocabulario e IDF do indice de documentos.
        Termos da query fora do vocabulario sao ignorados.
        """
        tokens = tokenize(query)
        tf_vec = compute_tf(tokens, self._vocab)
        tfidf_vec = tf_vec * self._idf
        return normalize_l2(tfidf_vec)

    def search(
        self,
        query: str,
        top_k: int = DEFAULT_TOP_K,
        threshold: float = SIMILARITY_THRESHOLD,
    ) -> list[dict[str, Any]]:
        """Busca artigos mais relevantes para uma query.

        Calcula similaridade coseno entre o vetor da query e todos os
        documentos indexados, retorna os top_k com score acima do limiar.

        Args:
            query: Pergunta em linguagem natural.
            top_k: Numero maximo de resultados.
            threshold: Score minimo para incluir no resultado.

        Returns:
            Lista de dicts com: pmid, title, abstract_snippet, score,
            citation, authors, year, journal.
        """
        if self._n_docs == 0:
            return []

        # Vetor da query
        q_vec = self._query_vector(query)

        # Similaridade coseno: dot product (vetores ja normalizados)
        scores = self._matrix @ q_vec

        # Ordenar por score decrescente
        ranked_indices = np.argsort(scores)[::-1]

        results: list[dict[str, Any]] = []
        for idx in ranked_indices[:top_k]:
            score = float(scores[idx])
            if score < threshold:
                break

            paper = self._papers[idx]
            results.append({
                "rank": len(results) + 1,
                "pmid": paper["pmid"],
                "title": paper["title"],
                "abstract_snippet": _truncate_abstract(paper["abstract"]),
                "score": round(score, 4),
                "citation": format_citation(paper),
                "authors": paper.get("authors", []),
                "year": paper.get("year", ""),
                "journal": paper.get("journal", ""),
            })

        return results

    def get_context_for_generation(
        self,
        query: str,
        top_k: int = 5,
    ) -> str:
        """Retorna contexto formatado para alimentar o gerador de respostas.

        Combina titulo e abstract dos artigos mais relevantes em um
        bloco de texto com citacoes PMID para uso pelo LLM ou template.

        Args:
            query: Pergunta do usuario.
            top_k: Numero de artigos para incluir no contexto.

        Returns:
            String formatada com contexto e citacoes.
        """
        results = self.search(query, top_k=top_k)

        if not results:
            return "Nenhum artigo relevante encontrado no indice."

        context_parts: list[str] = []
        for r in results:
            part = (
                f"[PMID:{r['pmid']}] {r['title']}\n"
                f"  {r['abstract_snippet']}\n"
                f"  Relevancia: {r['score']:.4f}"
            )
            context_parts.append(part)

        header = f"Contexto recuperado para: \"{query}\"\n"
        header += f"Artigos encontrados: {len(results)}\n"
        header += "-" * 60

        return header + "\n" + "\n\n".join(context_parts)

    @property
    def n_papers(self) -> int:
        """Numero total de artigos no indice."""
        return self._n_docs

    @property
    def vocab_size(self) -> int:
        """Tamanho do vocabulario TF-IDF."""
        return len(self._vocab)
