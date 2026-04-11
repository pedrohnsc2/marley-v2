"""Gerador de respostas para o RAG de leishmaniose.

Combina resultados do retriever com um modelo de linguagem para
gerar respostas citando fontes. Implementa abordagem em camadas:

    Tier 1 (sempre disponivel): Template com trechos dos artigos
    Tier 2 (se anthropic instalado): Claude API para geracao natural

O gerador nunca inventa informacao — todas as afirmacoes devem
ser rastreadas a um PMID especifico do corpus indexado.
"""

from __future__ import annotations

from typing import Any

from .retriever import PaperRetriever, format_citation


# ---------------------------------------------------------------------------
# Deteccao de tier disponivel
# ---------------------------------------------------------------------------

def _detect_tier() -> str:
    """Detecta o melhor tier de geracao disponivel.

    Returns:
        "claude" se anthropic SDK disponivel, "template" caso contrario.
    """
    try:
        import anthropic  # noqa: F401
        return "claude"
    except ImportError:
        return "template"


# ---------------------------------------------------------------------------
# Tier 1: Geracao por template (sempre disponivel)
# ---------------------------------------------------------------------------

def _generate_template_answer(
    query: str,
    results: list[dict[str, Any]],
) -> dict[str, Any]:
    """Gera resposta estruturada usando template com citacoes.

    Organiza os artigos recuperados em um formato legivel com
    destaque para relevancia e citacoes PMID. Nao requer LLM.

    Args:
        query: Pergunta original do usuario.
        results: Lista de artigos recuperados pelo retriever.

    Returns:
        Dict com answer, sources, tier e metadata.
    """
    if not results:
        return {
            "answer": (
                f"Nenhum artigo relevante encontrado para: \"{query}\"\n"
                "Tente reformular a pergunta ou indexar mais artigos."
            ),
            "sources": [],
            "tier": "template",
            "n_sources": 0,
        }

    # Construir resposta estruturada
    lines: list[str] = []
    lines.append(f"== Resultados para: \"{query}\" ==\n")
    lines.append(f"Artigos relevantes encontrados: {len(results)}\n")

    for r in results:
        lines.append(f"--- [{r['rank']}] Score: {r['score']:.4f} ---")
        lines.append(f"Titulo: {r['title']}")
        lines.append(f"Citacao: {r['citation']}")
        lines.append(f"Resumo: {r['abstract_snippet']}")
        lines.append("")

    # Sumario das fontes
    lines.append("== Fontes ==")
    for r in results:
        lines.append(f"  [{r['rank']}] PMID:{r['pmid']} — {r['title'][:80]}")

    answer = "\n".join(lines)
    sources = [
        {"pmid": r["pmid"], "title": r["title"], "score": r["score"]}
        for r in results
    ]

    return {
        "answer": answer,
        "sources": sources,
        "tier": "template",
        "n_sources": len(sources),
    }


# ---------------------------------------------------------------------------
# Tier 2: Geracao com Claude API (se disponivel)
# ---------------------------------------------------------------------------

def _generate_claude_answer(
    query: str,
    results: list[dict[str, Any]],
) -> dict[str, Any]:
    """Gera resposta natural usando Claude API com citacoes.

    Envia o contexto dos artigos recuperados como prompt e pede
    uma resposta fundamentada com citacoes PMID. Se a API falhar,
    faz fallback para o template.

    Args:
        query: Pergunta original do usuario.
        results: Lista de artigos recuperados pelo retriever.

    Returns:
        Dict com answer, sources, tier e metadata.
    """
    try:
        import anthropic
    except ImportError:
        return _generate_template_answer(query, results)

    if not results:
        return _generate_template_answer(query, results)

    # Montar contexto dos artigos
    context_parts: list[str] = []
    for r in results:
        context_parts.append(
            f"[PMID:{r['pmid']}] {r['title']}\n{r['abstract_snippet']}"
        )
    context = "\n\n".join(context_parts)

    # Prompt de sistema para geracao com citacoes
    system_prompt = (
        "Voce e um assistente cientifico especializado em leishmaniose e "
        "terapias antisense (ASO). Responda perguntas usando APENAS as "
        "informacoes fornecidas nos artigos abaixo. Cite cada afirmacao "
        "com o PMID correspondente no formato [PMID:XXXXXXXX]. Se a "
        "informacao nao estiver nos artigos, diga explicitamente. "
        "Responda em ingles cientifico."
    )

    user_prompt = (
        f"Artigos recuperados:\n\n{context}\n\n"
        f"Pergunta: {query}\n\n"
        f"Responda citando os PMIDs relevantes."
    )

    try:
        client = anthropic.Anthropic()
        response = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1024,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}],
        )
        answer_text = response.content[0].text

        sources = [
            {"pmid": r["pmid"], "title": r["title"], "score": r["score"]}
            for r in results
        ]

        return {
            "answer": answer_text,
            "sources": sources,
            "tier": "claude",
            "n_sources": len(sources),
            "model": "claude-sonnet-4-20250514",
            "usage": {
                "input_tokens": response.usage.input_tokens,
                "output_tokens": response.usage.output_tokens,
            },
        }

    except Exception as e:
        print(f"[generator] AVISO: Claude API falhou ({e}). Usando template.")
        result = _generate_template_answer(query, results)
        result["tier"] = "template (fallback de claude)"
        result["error"] = str(e)
        return result


# ---------------------------------------------------------------------------
# Classe Generator — interface publica
# ---------------------------------------------------------------------------

class RAGGenerator:
    """Gerador de respostas RAG com selecao automatica de tier.

    Detecta automaticamente o melhor metodo de geracao disponivel
    e fornece interface uniforme independente do tier.
    """

    def __init__(self, retriever: PaperRetriever) -> None:
        """Inicializa com um retriever ja configurado.

        Args:
            retriever: Instancia de PaperRetriever com indice carregado.
        """
        self._retriever = retriever
        self._tier = _detect_tier()
        print(f"[generator] Tier de geracao: {self._tier}")

    def answer(
        self,
        query: str,
        top_k: int = 5,
    ) -> dict[str, Any]:
        """Gera resposta para uma pergunta sobre leishmaniose.

        Busca artigos relevantes via retriever e gera resposta
        usando o melhor tier disponivel.

        Args:
            query: Pergunta em linguagem natural.
            top_k: Numero de artigos para contexto.

        Returns:
            Dict com: answer (str), sources (list), tier (str),
            n_sources (int), query (str).
        """
        # Recuperar artigos relevantes
        results = self._retriever.search(query, top_k=top_k)

        # Gerar resposta no tier adequado
        if self._tier == "claude":
            response = _generate_claude_answer(query, results)
        else:
            response = _generate_template_answer(query, results)

        response["query"] = query
        return response

    def batch_answer(
        self,
        queries: list[str],
        top_k: int = 5,
    ) -> list[dict[str, Any]]:
        """Gera respostas para multiplas perguntas.

        Util para avaliar o sistema com queries de teste.

        Args:
            queries: Lista de perguntas.
            top_k: Numero de artigos por query.

        Returns:
            Lista de respostas (mesma estrutura de answer()).
        """
        return [self.answer(q, top_k=top_k) for q in queries]

    @property
    def tier(self) -> str:
        """Tier de geracao atualmente em uso."""
        return self._tier
