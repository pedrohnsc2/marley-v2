"""InfoNCE contrastive loss para pares epitopo-DLA.

Implementacao em numpy puro da perda contrastiva InfoNCE (Noise
Contrastive Estimation), inspirada no CLIP (Radford et al., 2021).

A ideia central: dado um par positivo (epitopo, alelo_ligante), o
modelo deve atribuir similaridade alta a esse par e similaridade
baixa a pares negativos (epitopo, alelo_nao_ligante).

A temperatura tau controla a "nitidez" da distribuicao softmax —
valores menores tornam o modelo mais confiante nas decisoes.
"""

from __future__ import annotations

import numpy as np


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Similaridade cosseno entre dois vetores.

    Args:
        a: vetor de embedding do epitopo
        b: vetor de embedding do alelo

    Returns:
        Similaridade no intervalo [-1, 1]
    """
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a < 1e-10 or norm_b < 1e-10:
        return 0.0
    return float(np.dot(a, b) / (norm_a * norm_b))


def cosine_similarity_matrix(embeddings_a: np.ndarray, embeddings_b: np.ndarray) -> np.ndarray:
    """Matriz de similaridade cosseno entre dois conjuntos de embeddings.

    Calcula sim(a_i, b_j) para todo par (i, j), resultando em uma
    matriz N x M onde N = |a| e M = |b|.

    Args:
        embeddings_a: matriz (N, D) de embeddings de epitopos
        embeddings_b: matriz (M, D) de embeddings de alelos

    Returns:
        Matriz (N, M) de similaridades cosseno
    """
    # Normalizar L2 por linha
    norms_a = np.linalg.norm(embeddings_a, axis=1, keepdims=True)
    norms_b = np.linalg.norm(embeddings_b, axis=1, keepdims=True)

    # Evitar divisao por zero
    norms_a = np.maximum(norms_a, 1e-10)
    norms_b = np.maximum(norms_b, 1e-10)

    a_normalized = embeddings_a / norms_a
    b_normalized = embeddings_b / norms_b

    return a_normalized @ b_normalized.T


def info_nce_loss(
    anchor_embeddings: np.ndarray,
    positive_embeddings: np.ndarray,
    negative_embeddings: np.ndarray,
    temperature: float = 0.07,
) -> tuple[float, np.ndarray]:
    """InfoNCE loss para um batch de pares positivos + negativos.

    Para cada ancora (epitopo), calcula:
        L = -log( exp(sim(a, p+) / tau) / sum_j(exp(sim(a, n_j) / tau)) )

    onde p+ e o positivo e n_j inclui o positivo + todos os negativos.

    Args:
        anchor_embeddings: (B, D) embeddings dos epitopos (ancoras)
        positive_embeddings: (B, D) embeddings dos alelos positivos
        negative_embeddings: (B*neg_ratio, D) embeddings dos alelos negativos
        temperature: parametro tau que escala as similaridades

    Returns:
        Tupla (loss_media, gradientes_por_ancora)
    """
    batch_size = anchor_embeddings.shape[0]
    total_loss = 0.0
    # Gradientes da loss em relacao aos embeddings da ancora
    grad_anchors = np.zeros_like(anchor_embeddings)

    for i in range(batch_size):
        anchor = anchor_embeddings[i]  # (D,)
        positive = positive_embeddings[i]  # (D,)

        # Selecionar negativos para esta ancora
        # Negativos sao distribuidos uniformemente entre as ancoras
        neg_per_anchor = negative_embeddings.shape[0] // batch_size
        neg_start = i * neg_per_anchor
        neg_end = neg_start + neg_per_anchor
        negatives = negative_embeddings[neg_start:neg_end]  # (neg_ratio, D)

        # Concatenar positivo + negativos para o denominador
        all_candidates = np.vstack([positive.reshape(1, -1), negatives])  # (1+neg, D)

        # Similaridade cosseno escalada pela temperatura
        sims = np.array([
            cosine_similarity(anchor, cand) for cand in all_candidates
        ]) / temperature

        # Log-sum-exp trick para estabilidade numerica
        max_sim = np.max(sims)
        exp_sims = np.exp(sims - max_sim)
        log_sum_exp = max_sim + np.log(np.sum(exp_sims) + 1e-10)

        # Loss para esta ancora: -log P(positivo)
        positive_logit = sims[0]  # similaridade com o positivo
        loss_i = -positive_logit + log_sum_exp
        total_loss += loss_i

        # Gradiente simplificado da loss em relacao ao anchor embedding:
        # d(loss)/d(anchor) = -d(sim_pos)/d(anchor) + sum_j softmax_j * d(sim_j)/d(anchor)
        # Para similaridade cosseno: d(sim)/d(a) ~= (b/|b| - sim*a/|a|) / |a|
        softmax_probs = exp_sims / (np.sum(exp_sims) + 1e-10)
        anchor_norm = np.linalg.norm(anchor) + 1e-10

        grad_i = np.zeros_like(anchor)
        for j, cand in enumerate(all_candidates):
            cand_norm = np.linalg.norm(cand) + 1e-10
            sim_val = cosine_similarity(anchor, cand)
            # Gradiente da similaridade cosseno em relacao a 'anchor'
            dsim_da = (cand / (anchor_norm * cand_norm)) - (sim_val * anchor / (anchor_norm ** 2))
            weight = softmax_probs[j] - (1.0 if j == 0 else 0.0)
            grad_i += weight * dsim_da / temperature

        grad_anchors[i] = grad_i

    avg_loss = total_loss / batch_size
    grad_anchors /= batch_size

    return avg_loss, grad_anchors


def hard_negative_mining(
    anchor_embedding: np.ndarray,
    negative_embeddings: np.ndarray,
    top_k: int = 3,
) -> np.ndarray:
    """Seleciona os negativos mais dificeis (mais similares a ancora).

    Hard negatives sao fundamentais para o aprendizado contrastivo:
    forcam o modelo a aprender distincoes sutis entre binders e
    non-binders, em vez de apenas separar sequencias muito diferentes.

    Args:
        anchor_embedding: (D,) embedding do epitopo ancora
        negative_embeddings: (N, D) embeddings de todos os candidatos negativos
        top_k: numero de hard negatives a selecionar

    Returns:
        (top_k, D) embeddings dos negativos mais dificeis
    """
    if negative_embeddings.shape[0] <= top_k:
        return negative_embeddings

    sims = np.array([
        cosine_similarity(anchor_embedding, neg)
        for neg in negative_embeddings
    ])

    # Indices dos top_k negativos com maior similaridade (mais dificeis)
    hard_indices = np.argsort(sims)[-top_k:]
    return negative_embeddings[hard_indices]
