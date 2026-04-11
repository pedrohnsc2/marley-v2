"""Execucao do modulo 03_leish_esm — ESM-2 Protein Embeddings.

Pipeline completo:
    1. Carregar dataset de proteinas de L. infantum
    2. Carregar modelo ESM-2 (650M params) em MPS
    3. Extrair embeddings per-sequence (mean pooling) para todas
    4. Extrair embeddings per-residue para proteinas-fonte (analise de epitopos)
    5. Calcular matriz de similaridade coseno
    6. Clustering via k-means
    7. Analise de regioes epitopicas vs background
    8. Salvar embeddings (.npy) e resultados (JSON envelope)

Os embeddings gerados sao consumidos por:
    - Modulo 06 (EvoDiff): design generativo de peptideos
    - Modulo 07 (Contrastivo): aprendizado contrastivo ASO-proteina
    - Modulo 09 (SAE): interpretabilidade mecanistica

Uso:
    python -m marley_ai.03_leish_esm.run
"""

from __future__ import annotations

import logging
import sys
from typing import Any

from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

from .config import ESMConfig, get_default_config
from .dataset import (
    build_dataset,
    get_category_counts,
)
from .model import ESMModel
from .embeddings import (
    save_embeddings,
    cosine_similarity_matrix,
    find_similar_pairs,
    kmeans_cluster,
    cluster_summary,
    analyze_epitope_regions,
)

# Configurar logging para visibilidade no terminal
logging.basicConfig(
    level=logging.INFO,
    format="[%(name)s] %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger("03_leish_esm")


# ---------------------------------------------------------------------------
# Registro no orquestrador
# ---------------------------------------------------------------------------

@register("03_leish_esm")
class LeishESM:
    """Modulo ESM-2: embeddings proteicos para L. infantum.

    Gera representacoes vetoriais densas (1280-dim) usando o modelo de
    linguagem proteica ESM-2, treinado em ~250M sequencias do UniRef50.
    """

    def __init__(self) -> None:
        self.config: ESMConfig | None = None

    def configure(self, config: Any) -> None:
        if isinstance(config, ESMConfig):
            self.config = config
        else:
            self.config = get_default_config()

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica se torch e esm estao disponiveis."""
        missing: list[str] = []
        try:
            import torch
        except ImportError:
            missing.append("torch")
        try:
            import esm
        except ImportError:
            try:
                from transformers import AutoModel
            except ImportError:
                missing.append("esm ou transformers")
        return {"valid": len(missing) == 0, "missing": missing}

    def run(self) -> dict[str, Any]:
        return main(self.config)

    def get_dependencies(self) -> list[str]:
        """Sem dependencias — modulo de base (feature extraction)."""
        return []


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Gera embeddings ESM-2 para proteinas de L. infantum.

    Pipeline:
        1. Construir dataset (epitopos + fontes + drug targets + controles)
        2. Carregar ESM-2 no device MPS/CPU
        3. Extrair embeddings per-sequence (mean pooling)
        4. Extrair embeddings per-residue (proteinas-fonte apenas)
        5. Analise: similaridade, clustering, regioes epitopicas
        6. Salvar resultados

    Args:
        config: configuracao ESMConfig. Se None, usa defaults detectados.

    Returns:
        Envelope JSON com estatisticas e metadados dos embeddings.
    """
    if config is None or not isinstance(config, ESMConfig):
        config = get_default_config()

    envelope = create_envelope("03_leish_esm", version=config.version)
    envelope["device"] = config.device

    with Timer() as timer:
        try:
            result = _run_pipeline(config, envelope)
            envelope.update(result)
            envelope["status"] = "success"
        except Exception as e:
            logger.error("Pipeline falhou: %s", e, exc_info=True)
            envelope["status"] = "failed"
            envelope["warnings"].append(f"Erro na execucao: {e}")

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope, output_dir=config.results_dir)
    logger.info("Resultado salvo em %s (%.1fs)", output_path, timer.elapsed)

    return envelope


def _run_pipeline(config: ESMConfig, envelope: dict) -> dict[str, Any]:
    """Executa o pipeline interno (separado para tratamento de erros limpo).

    Returns:
        Dict com campos para merge no envelope.
    """
    # --- 1. Dataset ---
    logger.info("=== Fase 1: Construindo dataset ===")
    entries = build_dataset(n_controls=5, seed=config.seed)
    counts = get_category_counts(entries)
    logger.info(
        "Dataset: %d sequencias (%s)",
        len(entries),
        ", ".join(f"{cat}={n}" for cat, n in sorted(counts.items())),
    )

    # Preparar sequencias no formato (nome, sequencia)
    sequences = [(e.name, e.sequence) for e in entries]

    # --- 2. Carregar modelo ---
    logger.info("=== Fase 2: Carregando ESM-2 (%s) ===", config.esm_model)
    model = ESMModel(config)
    model.load()
    logger.info("Modelo carregado via %s no device %s", model.backend, model.device)

    # --- 3. Embeddings per-sequence (mean pooling) ---
    logger.info("=== Fase 3: Extraindo embeddings per-sequence ===")
    seq_embeddings = model.embed_sequences(sequences, pool_strategy="mean")
    logger.info(
        "Embeddings extraidos: %d sequencias, dim=%d",
        len(seq_embeddings), config.embedding_dim,
    )

    # Salvar embeddings per-sequence
    emb_dir = config.embeddings_dir
    saved_paths = save_embeddings(seq_embeddings, emb_dir)
    logger.info("Embeddings salvos em %s (%d arquivos)", emb_dir, len(saved_paths))

    # --- 4. Embeddings per-residue (proteinas-fonte apenas) ---
    logger.info("=== Fase 4: Extraindo embeddings per-residue (source proteins) ===")
    source_sequences = [
        (e.name, e.sequence) for e in entries if e.category == "source_protein"
    ]
    per_residue = model.embed_sequences(source_sequences, pool_strategy="per_residue")
    logger.info("Per-residue extraidos para %d proteinas-fonte", len(per_residue))

    # Salvar per-residue em subdiretorio separado
    per_residue_dir = emb_dir / "per_residue"
    save_embeddings(per_residue, per_residue_dir)

    # --- 5. Analise de similaridade ---
    logger.info("=== Fase 5: Analise de similaridade coseno ===")
    sim_matrix, sim_names = cosine_similarity_matrix(seq_embeddings)
    similar_pairs = find_similar_pairs(sim_matrix, sim_names, threshold=0.80)

    if similar_pairs:
        logger.info("Top 5 pares mais similares:")
        for a, b, sim in similar_pairs[:5]:
            logger.info("  %.3f  %s <-> %s", sim, a, b)
    else:
        logger.info("Nenhum par com similaridade > 0.80")

    # --- 6. Clustering ---
    logger.info("=== Fase 6: Clustering k-means (k=4) ===")
    clusters = kmeans_cluster(seq_embeddings, k=4, seed=config.seed)
    cluster_groups = cluster_summary(clusters, entries)

    for cid, members in sorted(cluster_groups.items()):
        logger.info("  Cluster %d: %d membros", cid, len(members))
        for m in members[:3]:  # mostrar ate 3 membros por cluster
            logger.info("    - %s", m)
        if len(members) > 3:
            logger.info("    ... e mais %d", len(members) - 3)

    # --- 7. Analise de regioes epitopicas ---
    logger.info("=== Fase 7: Analise de regioes epitopicas ===")
    epitope_analysis = analyze_epitope_regions(per_residue, entries)

    for name, analysis in epitope_analysis.items():
        logger.info(
            "  %s: divergencia=%.4f, norm_ratio=%.4f",
            name, analysis["divergence_score"], analysis["norm_ratio"],
        )

    # --- 8. Montar resultados ---
    logger.info("=== Fase 8: Montando envelope de resultados ===")

    # Converter similarity matrix para formato serializavel
    sim_data = {
        "names": sim_names,
        "matrix_shape": list(sim_matrix.shape),
        "mean_similarity": float(sim_matrix.mean()),
        "median_similarity": float(_median(sim_matrix)),
    }

    # Estatisticas dos embeddings
    all_emb = list(seq_embeddings.values())
    emb_norms = [float(e.dot(e) ** 0.5) for e in all_emb]

    return {
        "summary": {
            "conclusion": (
                f"Embeddings ESM-2 gerados para {len(entries)} sequencias de "
                f"L. infantum ({config.esm_model}, dim={config.embedding_dim}). "
                f"Similaridade media={sim_data['mean_similarity']:.3f}. "
                f"{len(similar_pairs)} pares com sim > 0.80. "
                f"Analise epitopica em {len(epitope_analysis)} proteinas-fonte."
            ),
            "key_metrics": {
                "n_sequences": len(entries),
                "embedding_dim": config.embedding_dim,
                "model": config.esm_model,
                "backend": model.backend,
                "device": str(model.device),
                "n_epitopes": counts.get("epitope", 0),
                "n_source_proteins": counts.get("source_protein", 0),
                "n_drug_targets": counts.get("drug_target", 0),
                "n_controls": counts.get("control", 0),
                "mean_similarity": round(sim_data["mean_similarity"], 4),
                "n_similar_pairs": len(similar_pairs),
                "n_clusters": len(cluster_groups),
                "n_epitope_analyses": len(epitope_analysis),
            },
        },
        "artifacts": saved_paths,
        "metrics": {
            "embedding_norms": {
                "mean": round(float(sum(emb_norms) / len(emb_norms)), 4),
                "min": round(float(min(emb_norms)), 4),
                "max": round(float(max(emb_norms)), 4),
            },
            "similarity": sim_data,
        },
        "data": {
            "dataset_counts": counts,
            "similar_pairs": [
                {"a": a, "b": b, "similarity": round(s, 4)}
                for a, b, s in similar_pairs[:20]  # top 20 para o JSON
            ],
            "clusters": {
                str(cid): members
                for cid, members in cluster_groups.items()
            },
            "epitope_analysis": epitope_analysis,
            "sequences": {
                e.name: {
                    "category": e.category,
                    "length": len(e.sequence),
                    "gene_id": e.gene_id,
                    "description": e.description,
                }
                for e in entries
            },
        },
        "dependencies": [],
    }


def _median(arr):
    """Calcula mediana de um array numpy (sem importar scipy)."""
    flat = arr.flatten()
    return float(sorted(flat)[len(flat) // 2])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    result = main()
    status = result.get("status", "unknown")
    if status == "success":
        metrics = result.get("summary", {}).get("key_metrics", {})
        print(f"\n{'='*60}")
        print(f"  03_leish_esm CONCLUIDO — {metrics.get('n_sequences', 0)} sequencias")
        print(f"  Modelo: {metrics.get('model', '?')} via {metrics.get('backend', '?')}")
        print(f"  Device: {metrics.get('device', '?')}")
        print(f"  Similaridade media: {metrics.get('mean_similarity', 0):.4f}")
        print(f"  Pares similares (>0.80): {metrics.get('n_similar_pairs', 0)}")
        print(f"{'='*60}")
    else:
        print(f"\n[FALHA] Status: {status}")
        for w in result.get("warnings", []):
            print(f"  AVISO: {w}")
        sys.exit(1)
