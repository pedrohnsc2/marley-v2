"""Execucao do modulo 01_rag — RAG sobre dados do projeto Marley.

Indexa arquivos JSON, MD e CSV do projeto em um indice vetorial FAISS,
usando embeddings de sentence-transformers (all-MiniLM-L6-v2, 384-dim).
Valida com queries de teste cobrindo os eixos do projeto: ASO, vacina,
drug targets, quantum e resistance.

Uso:
    python -m marley_ai.01_rag.run
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from marley_ai.config import AIModuleConfig, PROJECT_ROOT
from marley_ai.envelope import Timer, create_envelope, write_result


# ---------------------------------------------------------------------------
# Constantes — diretorios do corpus e parametros
# ---------------------------------------------------------------------------

CORPUS_DIRS: list[Path] = [
    PROJECT_ROOT / "results",
    PROJECT_ROOT / "aso_math" / "results",
    PROJECT_ROOT / "aso_delivery" / "results",
    PROJECT_ROOT / "marley_ai" / "results",
    PROJECT_ROOT / "vaccine_platforms" / "reports" / "results",
    PROJECT_ROOT / "docs",
    PROJECT_ROOT / "mrl_quantum" / "results",
]

SUPPORTED_EXTENSIONS: set[str] = {".json", ".md"}

CHUNK_SIZE: int = 500        # caracteres por chunk
CHUNK_OVERLAP: int = 50      # overlap entre chunks consecutivos
EMBEDDING_MODEL: str = "all-MiniLM-L6-v2"
EMBEDDING_DIM: int = 384
TOP_K: int = 3               # resultados por query de validacao

VALIDATION_QUERIES: list[str] = [
    "What is the binding energy of MRL-ASO-001?",
    "Which vaccine platform has the lowest cost?",
    "How many drug targets were identified?",
    "What is the QAOA optimal configuration?",
    "What is the resistance barrier for the ASO therapy?",
]


# ---------------------------------------------------------------------------
# Corpus — coleta e chunking
# ---------------------------------------------------------------------------

def _collect_documents() -> list[dict[str, Any]]:
    """Walk corpus directories and read JSON/MD files.

    Returns list of dicts with keys: source_file, content, extension.
    Skips files that cannot be decoded or are empty.
    """
    documents: list[dict[str, Any]] = []

    for corpus_dir in CORPUS_DIRS:
        if not corpus_dir.exists():
            continue
        for root, _dirs, files in os.walk(corpus_dir):
            for fname in sorted(files):
                fpath = Path(root) / fname
                ext = fpath.suffix.lower()
                if ext not in SUPPORTED_EXTENSIONS:
                    continue
                try:
                    raw = fpath.read_text(encoding="utf-8")
                except (OSError, UnicodeDecodeError):
                    continue
                if not raw.strip():
                    continue

                # Skip our own output to avoid self-referential indexing
                if fpath.name == "01_rag.json":
                    continue

                # Flatten JSON to readable text
                if ext == ".json":
                    content = _flatten_json(raw, str(fpath.relative_to(PROJECT_ROOT)))
                else:
                    content = raw

                documents.append({
                    "source_file": str(fpath.relative_to(PROJECT_ROOT)),
                    "content": content,
                    "extension": ext,
                })

    return documents


def _flatten_json(raw: str, source_label: str) -> str:
    """Convert JSON string to readable key-value text.

    For nested objects, flattens with dot notation. Arrays are summarized
    to avoid extremely long outputs.
    """
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        return raw

    lines: list[str] = [f"[Source: {source_label}]"]
    _flatten_recursive(data, "", lines, depth=0, max_depth=4)
    return "\n".join(lines)


def _flatten_recursive(
    obj: Any,
    prefix: str,
    lines: list[str],
    depth: int,
    max_depth: int,
) -> None:
    """Recursively flatten a JSON object to key: value lines."""
    if depth > max_depth:
        lines.append(f"{prefix}: [nested data truncated]")
        return

    if isinstance(obj, dict):
        for key, value in obj.items():
            new_prefix = f"{prefix}.{key}" if prefix else key
            _flatten_recursive(value, new_prefix, lines, depth + 1, max_depth)
    elif isinstance(obj, list):
        if len(obj) == 0:
            lines.append(f"{prefix}: []")
        elif len(obj) <= 5:
            for i, item in enumerate(obj):
                _flatten_recursive(item, f"{prefix}[{i}]", lines, depth + 1, max_depth)
        else:
            # Summarize long arrays: first 3 + count
            for i in range(3):
                _flatten_recursive(obj[i], f"{prefix}[{i}]", lines, depth + 1, max_depth)
            lines.append(f"{prefix}: ... ({len(obj)} items total)")
    else:
        # Scalar value — truncate very long strings
        s = str(obj)
        if len(s) > 300:
            s = s[:300] + "..."
        lines.append(f"{prefix}: {s}")


def _chunk_documents(
    documents: list[dict[str, Any]],
    chunk_size: int = CHUNK_SIZE,
    overlap: int = CHUNK_OVERLAP,
) -> tuple[list[str], list[dict[str, Any]]]:
    """Split documents into overlapping text chunks.

    Returns:
        Tuple of (chunks_text, chunks_metadata).
        Each chunk_metadata has: source_file, chunk_index, char_start, char_end.
    """
    chunks: list[str] = []
    metadata: list[dict[str, Any]] = []

    for doc in documents:
        text = doc["content"]
        source = doc["source_file"]
        idx = 0
        chunk_i = 0
        while idx < len(text):
            end = min(idx + chunk_size, len(text))
            chunk_text = text[idx:end]
            # Only add non-trivial chunks
            if len(chunk_text.strip()) > 20:
                chunks.append(chunk_text)
                metadata.append({
                    "source_file": source,
                    "chunk_index": chunk_i,
                    "char_start": idx,
                    "char_end": end,
                })
                chunk_i += 1
            idx += chunk_size - overlap

    return chunks, metadata


# ---------------------------------------------------------------------------
# Embedding + FAISS indexing
# ---------------------------------------------------------------------------

def _build_index(
    chunks: list[str],
) -> tuple[Any, Any, Any]:
    """Embed chunks and build FAISS index.

    Returns:
        Tuple of (faiss_index, embeddings_array, model).
    """
    # Fix OpenMP conflict on macOS (torch + faiss both link libomp)
    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

    from sentence_transformers import SentenceTransformer
    import faiss

    print(f"[01_rag] Carregando modelo {EMBEDDING_MODEL}...")
    model = SentenceTransformer(EMBEDDING_MODEL)

    print(f"[01_rag] Gerando embeddings para {len(chunks)} chunks...")
    embeddings = model.encode(
        chunks,
        show_progress_bar=False,
        normalize_embeddings=False,
        convert_to_numpy=True,
    )
    embeddings = embeddings.astype(np.float32)

    # L2-normalize for cosine similarity via inner product
    faiss.normalize_L2(embeddings)

    # Build flat inner-product index (exact search, cosine after normalization)
    dim = embeddings.shape[1]
    index = faiss.IndexFlatIP(dim)
    index.add(embeddings)

    print(f"[01_rag] Indice FAISS criado: {index.ntotal} vetores, {dim}-dim")
    return index, embeddings, model


# ---------------------------------------------------------------------------
# Retrieval + validation
# ---------------------------------------------------------------------------

def _search(
    query: str,
    model: Any,
    index: Any,
    chunks: list[str],
    metadata: list[dict[str, Any]],
    top_k: int = TOP_K,
) -> list[dict[str, Any]]:
    """Search the FAISS index for a query and return top-k results."""
    import faiss as _faiss

    q_emb = model.encode([query], convert_to_numpy=True).astype(np.float32)
    _faiss.normalize_L2(q_emb)

    scores, indices = index.search(q_emb, top_k)

    results: list[dict[str, Any]] = []
    for rank, (score, idx) in enumerate(zip(scores[0], indices[0]), start=1):
        if idx < 0:
            continue
        results.append({
            "rank": rank,
            "score": round(float(score), 4),
            "chunk_preview": chunks[idx][:200],
            "source_file": metadata[idx]["source_file"],
            "chunk_index": metadata[idx]["chunk_index"],
        })

    return results


def _run_validation(
    queries: list[str],
    model: Any,
    index: Any,
    chunks: list[str],
    metadata: list[dict[str, Any]],
    top_k: int = TOP_K,
) -> tuple[list[dict[str, Any]], float]:
    """Run validation queries and compute mean top-1 score.

    Returns:
        Tuple of (validation_results, mean_top1_score).
    """
    validation_results: list[dict[str, Any]] = []
    top1_scores: list[float] = []

    for query in queries:
        results = _search(query, model, index, chunks, metadata, top_k=top_k)
        top1_score = results[0]["score"] if results else 0.0
        top1_scores.append(top1_score)

        validation_results.append({
            "query": query,
            "top_k_results": results,
            "top1_score": round(top1_score, 4),
        })

    mean_score = float(np.mean(top1_scores)) if top1_scores else 0.0
    return validation_results, mean_score


# ---------------------------------------------------------------------------
# Corpus statistics
# ---------------------------------------------------------------------------

def _corpus_stats(documents: list[dict[str, Any]]) -> dict[str, Any]:
    """Compute summary statistics over the corpus."""
    n_json = sum(1 for d in documents if d["extension"] == ".json")
    n_md = sum(1 for d in documents if d["extension"] == ".md")
    total_chars = sum(len(d["content"]) for d in documents)
    return {
        "n_json": n_json,
        "n_md": n_md,
        "n_documents": len(documents),
        "total_chars": total_chars,
    }


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Executa o pipeline de RAG sobre dados do projeto Marley.

    Etapas:
        1. Coleta documentos (JSON + MD) de todos os diretorios do projeto
        2. Divide em chunks de ~500 caracteres com 50 char overlap
        3. Gera embeddings com all-MiniLM-L6-v2 (384-dim)
        4. Constroi indice FAISS (IndexFlatIP para cosine similarity)
        5. Valida com 5 queries de teste cobrindo cada eixo do projeto
        6. Salva envelope com metricas e resultados de validacao

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com resultados do RAG.
    """
    envelope = create_envelope("01_rag")
    envelope["device"] = "cpu"

    with Timer() as timer:
        try:
            # 1. Coleta de documentos
            print("[01_rag] Coletando documentos do corpus...")
            documents = _collect_documents()
            n_docs = len(documents)
            stats = _corpus_stats(documents)
            print(f"[01_rag]   {n_docs} documentos ({stats['n_json']} JSON, {stats['n_md']} MD)")
            print(f"[01_rag]   {stats['total_chars']:,} caracteres totais")

            if n_docs == 0:
                envelope["status"] = "failed"
                envelope["summary"]["conclusion"] = "Nenhum documento encontrado no corpus."
                envelope["runtime_seconds"] = 0.0
                output_path = write_result(envelope)
                print(f"[01_rag] Resultado salvo em {output_path}")
                return envelope

            # 2. Chunking
            print("[01_rag] Dividindo em chunks...")
            chunks, metadata = _chunk_documents(documents)
            n_chunks = len(chunks)
            print(f"[01_rag]   {n_chunks} chunks gerados (size={CHUNK_SIZE}, overlap={CHUNK_OVERLAP})")

            # 3 + 4. Embedding + FAISS indexing
            print("[01_rag] Construindo indice vetorial...")
            faiss_index, embeddings, model = _build_index(chunks)
            dim = embeddings.shape[1]

            # 5. Validacao
            print(f"[01_rag] Executando {len(VALIDATION_QUERIES)} queries de validacao...")
            validation_results, mean_score = _run_validation(
                VALIDATION_QUERIES, model, faiss_index, chunks, metadata,
            )

            for vr in validation_results:
                q = vr["query"][:60]
                s = vr["top1_score"]
                src = vr["top_k_results"][0]["source_file"] if vr["top_k_results"] else "N/A"
                print(f"  Q: {q}...")
                print(f"    Top-1 score={s:.4f}  source={src}")

            # 6. Montagem do envelope
            envelope["status"] = "complete"
            envelope["dependencies"] = []

            envelope["summary"]["conclusion"] = (
                f"RAG indexou {n_docs} documentos em {n_chunks} chunks. "
                f"{dim}-dim embeddings via {EMBEDDING_MODEL}. "
                f"Top query recall demonstrado para {len(VALIDATION_QUERIES)} "
                f"queries de validacao."
            )
            envelope["summary"]["key_metrics"] = {
                "corpus_size": n_docs,
                "n_chunks": n_chunks,
                "index_size": n_chunks,
                "embedding_dim": dim,
                "model": EMBEDDING_MODEL,
                "n_validation_queries": len(VALIDATION_QUERIES),
                "mean_top1_score": round(mean_score, 4),
            }

            envelope["data"] = {
                "validation_results": validation_results,
                "corpus_stats": stats,
                "corpus_dirs": [str(d.relative_to(PROJECT_ROOT)) for d in CORPUS_DIRS],
                "parameters": {
                    "chunk_size": CHUNK_SIZE,
                    "chunk_overlap": CHUNK_OVERLAP,
                    "embedding_model": EMBEDDING_MODEL,
                    "embedding_dim": dim,
                    "top_k": TOP_K,
                },
            }

            envelope["metrics"] = {
                "mean_top1_score": round(mean_score, 4),
                "min_top1_score": round(
                    min(vr["top1_score"] for vr in validation_results), 4
                ),
                "max_top1_score": round(
                    max(vr["top1_score"] for vr in validation_results), 4
                ),
                "corpus_documents": n_docs,
                "corpus_chunks": n_chunks,
                "corpus_chars": stats["total_chars"],
            }

            envelope["artifacts"] = [
                "FAISS IndexFlatIP (in-memory, not persisted)",
                f"Embeddings: {n_chunks}x{dim} float32 matrix",
            ]

        except Exception as exc:
            envelope["status"] = "failed"
            envelope["summary"]["conclusion"] = f"RAG falhou: {exc}"
            envelope["warnings"].append(str(exc))

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"[01_rag] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
