"""Extracao, analise e persistencia de embeddings ESM-2.

Responsavel por:
    1. Extrair embeddings para todas as sequencias do dataset
    2. Salvar como .npy para consumo por outros modulos (06, 07, 09)
    3. Calcular similaridade coseno entre pares de sequencias
    4. Identificar regioes atipicas (potenciais sitios funcionais)
    5. Clustering basico via k-means (implementacao numpy pura)

Os embeddings sao vetores de 1280 dimensoes que capturam propriedades
estruturais e funcionais de cada proteina. Sequencias com funcao similar
tendem a ficar proximas no espaco de embedding, mesmo sem homologia
de sequencia detectavel por BLAST.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

from .dataset import ProteinEntry

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Persistencia de embeddings
# ---------------------------------------------------------------------------

def save_embeddings(
    embeddings: dict[str, np.ndarray],
    output_dir: Path,
) -> list[str]:
    """Salva embeddings como arquivos .npy individuais.

    Cada sequencia gera um arquivo separado para facilitar carregamento
    seletivo por outros modulos. Tambem salva um indice JSON mapeando
    nomes para caminhos.

    Args:
        embeddings: dict nome -> numpy array
        output_dir: diretorio de saida (criado se nao existir)

    Returns:
        Lista de caminhos dos arquivos salvos.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    saved_paths: list[str] = []
    index: dict[str, dict[str, str | list[int]]] = {}

    for name, emb in embeddings.items():
        # Sanitizar nome para uso como filename
        safe_name = name.replace("/", "_").replace(" ", "_")
        path = output_dir / f"{safe_name}.npy"
        np.save(str(path), emb)
        saved_paths.append(str(path))

        index[name] = {
            "path": str(path),
            "shape": list(emb.shape),
            "dtype": str(emb.dtype),
        }

    # Salvar indice para descoberta por outros modulos
    index_path = output_dir / "embedding_index.json"
    with open(index_path, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)
    saved_paths.append(str(index_path))

    logger.info(
        "Salvos %d embeddings em %s (indice: %s)",
        len(embeddings), output_dir, index_path,
    )
    return saved_paths


def load_embedding(name: str, embeddings_dir: Path) -> np.ndarray:
    """Carrega um embedding salvo pelo nome.

    Funcao utilitaria para outros modulos (06, 07, 09) que consomem
    os embeddings gerados por este modulo.

    Args:
        name: nome da sequencia (como usado em save_embeddings)
        embeddings_dir: diretorio onde os .npy estao salvos

    Returns:
        numpy array com o embedding

    Raises:
        FileNotFoundError: se o embedding nao existir
    """
    safe_name = name.replace("/", "_").replace(" ", "_")
    path = embeddings_dir / f"{safe_name}.npy"
    if not path.exists():
        raise FileNotFoundError(
            f"Embedding '{name}' nao encontrado em {path}. "
            f"Execute o modulo 03_leish_esm primeiro."
        )
    return np.load(str(path))


# ---------------------------------------------------------------------------
# Analise de similaridade
# ---------------------------------------------------------------------------

def cosine_similarity_matrix(
    embeddings: dict[str, np.ndarray],
) -> tuple[np.ndarray, list[str]]:
    """Calcula matriz de similaridade coseno entre todas as sequencias.

    Similaridade coseno mede o angulo entre vetores no espaco de 1280 dims.
    Valor 1.0 = identicos, 0.0 = ortogonais, -1.0 = opostos.

    Para proteinas com funcao similar, esperamos sim > 0.8.
    Controles aleatorios devem ter sim < 0.5 contra proteinas reais.

    Args:
        embeddings: dict nome -> embedding (1D array)

    Returns:
        Tupla (matriz_similaridade, lista_de_nomes) onde a matriz e NxN
        e os nomes correspondem aos indices.
    """
    names = sorted(embeddings.keys())
    n = len(names)
    matrix = np.zeros((n, n), dtype=np.float32)

    # Empilhar embeddings e normalizar
    vectors = np.stack([embeddings[name] for name in names])
    # Normalizar para unit vectors (evita divisao por zero)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    norms = np.where(norms == 0, 1.0, norms)
    normalized = vectors / norms

    # Similaridade coseno = produto escalar entre vetores normalizados
    matrix = normalized @ normalized.T

    return matrix, names


def find_similar_pairs(
    sim_matrix: np.ndarray,
    names: list[str],
    threshold: float = 0.85,
) -> list[tuple[str, str, float]]:
    """Encontra pares de sequencias com alta similaridade.

    Util para identificar proteinas com funcao potencialmente similar
    mesmo sem homologia de sequencia detectavel.

    Args:
        sim_matrix: matriz NxN de similaridade coseno
        names: nomes correspondentes aos indices
        threshold: similaridade minima para reportar

    Returns:
        Lista de (nome_a, nome_b, similaridade) ordenada por similaridade desc.
    """
    pairs: list[tuple[str, str, float]] = []
    n = len(names)

    for i in range(n):
        for j in range(i + 1, n):
            sim = float(sim_matrix[i, j])
            if sim >= threshold:
                pairs.append((names[i], names[j], sim))

    pairs.sort(key=lambda x: x[2], reverse=True)
    return pairs


# ---------------------------------------------------------------------------
# Clustering via k-means (implementacao numpy pura)
# ---------------------------------------------------------------------------

def kmeans_cluster(
    embeddings: dict[str, np.ndarray],
    k: int = 4,
    max_iter: int = 100,
    seed: int = 42,
) -> dict[str, int]:
    """Agrupa sequencias em k clusters via k-means.

    Implementacao simples com numpy — nao depende de sklearn.
    Inicializacao k-means++ para convergencia mais rapida.

    Clusters esperados para o dataset Marley:
        - Cluster 0: epitopos curtos (9-mers)
        - Cluster 1: proteinas-fonte (fragmentos longos)
        - Cluster 2: drug targets
        - Cluster 3: controles aleatorios

    Args:
        embeddings: dict nome -> embedding 1D
        k: numero de clusters
        max_iter: maximo de iteracoes
        seed: semente para reproducibilidade

    Returns:
        Dict mapeando nome -> cluster_id (0-indexed)
    """
    rng = np.random.RandomState(seed)
    names = sorted(embeddings.keys())
    X = np.stack([embeddings[name] for name in names])
    n, d = X.shape

    if n <= k:
        # Se temos menos pontos que clusters, cada ponto e seu cluster
        return {name: i for i, name in enumerate(names)}

    # --- Inicializacao k-means++ ---
    centroids = np.zeros((k, d), dtype=np.float32)
    # Primeiro centroide: aleatorio
    idx = rng.randint(0, n)
    centroids[0] = X[idx]

    for c in range(1, k):
        # Distancia de cada ponto ao centroide mais proximo
        dists = np.min(
            np.sum((X[:, None, :] - centroids[None, :c, :]) ** 2, axis=2),
            axis=1,
        )
        # Probabilidade proporcional a distancia^2
        probs = dists / dists.sum()
        idx = rng.choice(n, p=probs)
        centroids[c] = X[idx]

    # --- Iteracoes de Lloyd ---
    labels = np.zeros(n, dtype=np.int32)

    for iteration in range(max_iter):
        # Atribuir cada ponto ao centroide mais proximo
        dists = np.sum((X[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(dists, axis=1).astype(np.int32)

        # Verificar convergencia
        if np.array_equal(labels, new_labels):
            logger.info("K-means convergiu em %d iteracoes", iteration + 1)
            break
        labels = new_labels

        # Recalcular centroides
        for c in range(k):
            mask = labels == c
            if mask.any():
                centroids[c] = X[mask].mean(axis=0)

    return {name: int(labels[i]) for i, name in enumerate(names)}


def cluster_summary(
    clusters: dict[str, int],
    entries: list[ProteinEntry],
) -> dict[int, list[str]]:
    """Resume quais sequencias estao em cada cluster.

    Util para avaliar se o clustering separa biologicamente:
    proteinas com funcao similar devem ficar no mesmo cluster.

    Args:
        clusters: dict nome -> cluster_id
        entries: lista original de ProteinEntry

    Returns:
        Dict cluster_id -> lista de nomes naquele cluster.
    """
    summary: dict[int, list[str]] = {}
    for name, cid in sorted(clusters.items(), key=lambda x: x[1]):
        summary.setdefault(cid, []).append(name)
    return summary


# ---------------------------------------------------------------------------
# Analise de regioes atipicas (per-residue)
# ---------------------------------------------------------------------------

def analyze_epitope_regions(
    per_residue_embeddings: dict[str, np.ndarray],
    entries: list[ProteinEntry],
) -> dict[str, dict]:
    """Analisa como regioes epitopicas diferem do restante da proteina.

    Para cada proteina-fonte, compara o embedding medio dos residuos
    na regiao do epitopo com o embedding medio do restante. Regioes
    com alta divergencia podem ser sitios funcionalmente importantes.

    Esta analise e possivel apenas para proteinas-fonte que contem
    os epitopos como subsequencias.

    Args:
        per_residue_embeddings: dict nome -> array (seq_len, embed_dim)
        entries: lista de ProteinEntry com metadados

    Returns:
        Dict com analise por proteina-fonte:
            - epitope_positions: posicoes encontradas
            - divergence_score: divergencia coseno entre epitopo e background
            - epitope_mean_norm: norma L2 media na regiao do epitopo
            - background_mean_norm: norma L2 media fora do epitopo
    """
    from vaccine_platforms.shared.epitopes import EPITOPES

    # Mapear gene_id -> epitopos
    gene_epitopes: dict[str, list] = {}
    for ep in EPITOPES:
        gene_epitopes.setdefault(ep.gene_id, []).append(ep)

    results: dict[str, dict] = {}

    # Analisar apenas proteinas-fonte com embeddings per-residue
    source_entries = {
        e.name: e for e in entries if e.category == "source_protein"
    }

    for name, entry in source_entries.items():
        if name not in per_residue_embeddings:
            continue

        emb = per_residue_embeddings[name]  # (seq_len, embed_dim)
        seq = entry.sequence
        seq_len = emb.shape[0]

        # Encontrar posicoes dos epitopos na sequencia
        epitope_positions: list[tuple[int, int, str]] = []
        for ep in gene_epitopes.get(entry.gene_id, []):
            pos = seq.find(ep.peptide)
            if pos >= 0 and pos + len(ep.peptide) <= seq_len:
                epitope_positions.append((pos, pos + len(ep.peptide), ep.peptide))

        if not epitope_positions:
            continue

        # Criar mascara para residuos epitopicos
        epitope_mask = np.zeros(seq_len, dtype=bool)
        for start, end, _ in epitope_positions:
            epitope_mask[start:end] = True

        background_mask = ~epitope_mask

        if not epitope_mask.any() or not background_mask.any():
            continue

        # Embeddings medios das regioes
        epitope_mean = emb[epitope_mask].mean(axis=0)
        background_mean = emb[background_mask].mean(axis=0)

        # Divergencia coseno: 1 - cos_sim (0 = identicos, 2 = opostos)
        cos_sim = float(
            np.dot(epitope_mean, background_mean)
            / (np.linalg.norm(epitope_mean) * np.linalg.norm(background_mean) + 1e-8)
        )
        divergence = 1.0 - cos_sim

        # Normas medias — regioes funcionais tendem a ter normas maiores
        epitope_norms = float(np.linalg.norm(emb[epitope_mask], axis=1).mean())
        background_norms = float(np.linalg.norm(emb[background_mask], axis=1).mean())

        results[name] = {
            "epitope_positions": [
                {"start": s, "end": e, "peptide": p}
                for s, e, p in epitope_positions
            ],
            "divergence_score": round(divergence, 4),
            "cosine_similarity": round(cos_sim, 4),
            "epitope_mean_norm": round(epitope_norms, 4),
            "background_mean_norm": round(background_norms, 4),
            "norm_ratio": round(epitope_norms / (background_norms + 1e-8), 4),
        }

    return results
