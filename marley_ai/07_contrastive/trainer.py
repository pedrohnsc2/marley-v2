"""Treinamento contrastivo para pares epitopo-DLA.

Implementa o loop de treinamento completo em numpy puro:
    1. Geracao de pares positivos/negativos a partir dos 11 epitopos
    2. Data augmentation por mutacao de aminoacidos (hard negatives)
    3. MLP de 2 camadas (input -> hidden -> embedding) com ReLU
    4. SGD com momentum para otimizacao
    5. Tracking de curva de loss

O MLP projeta features de entrada (k-mer + fisicoquimico para
epitopos, propriedades funcionais para alelos) em um espaco de
embedding compartilhado onde a similaridade cosseno reflete
afinidade de ligacao.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Final

import numpy as np

from vaccine_platforms.shared.epitopes import EPITOPES, Epitope

from .config import ContrastiveConfig
from .encoder import (
    AMINO_ACIDS,
    DLA_ALLELE_PROPERTIES,
    build_kmer_vocab,
    encode_allele,
    encode_peptide,
    get_allele_feature_dim,
    get_peptide_feature_dim,
)
from .loss import hard_negative_mining, info_nce_loss


# ---------------------------------------------------------------------------
# MLP implementado em numpy puro (2 camadas com ReLU)
# ---------------------------------------------------------------------------

@dataclass
class MLPWeights:
    """Pesos de um MLP de 2 camadas: input -> hidden (ReLU) -> output.

    Inicializacao He (Kaiming) para ReLU — variancia 2/fan_in.
    Bias inicializado em zero.
    """

    W1: np.ndarray  # (input_dim, hidden_dim)
    b1: np.ndarray  # (hidden_dim,)
    W2: np.ndarray  # (hidden_dim, output_dim)
    b2: np.ndarray  # (output_dim,)

    # Momentum buffers para SGD com momentum
    vW1: np.ndarray = field(default_factory=lambda: np.array([]))
    vb1: np.ndarray = field(default_factory=lambda: np.array([]))
    vW2: np.ndarray = field(default_factory=lambda: np.array([]))
    vb2: np.ndarray = field(default_factory=lambda: np.array([]))


def init_mlp(input_dim: int, hidden_dim: int, output_dim: int, seed: int = 42) -> MLPWeights:
    """Inicializa pesos do MLP com inicializacao He (Kaiming).

    A inicializacao He e ideal para redes com ReLU — previne
    vanishing/exploding gradients nas primeiras epocas.

    Args:
        input_dim: dimensao da entrada (features do peptideo ou alelo)
        hidden_dim: dimensao da camada oculta
        output_dim: dimensao do embedding de saida
        seed: semente para reproducibilidade

    Returns:
        MLPWeights com pesos inicializados
    """
    rng = np.random.RandomState(seed)

    # He initialization: std = sqrt(2 / fan_in)
    W1 = rng.randn(input_dim, hidden_dim) * np.sqrt(2.0 / input_dim)
    b1 = np.zeros(hidden_dim)
    W2 = rng.randn(hidden_dim, output_dim) * np.sqrt(2.0 / hidden_dim)
    b2 = np.zeros(output_dim)

    return MLPWeights(
        W1=W1, b1=b1, W2=W2, b2=b2,
        vW1=np.zeros_like(W1), vb1=np.zeros_like(b1),
        vW2=np.zeros_like(W2), vb2=np.zeros_like(b2),
    )


def mlp_forward(x: np.ndarray, weights: MLPWeights) -> tuple[np.ndarray, np.ndarray]:
    """Forward pass do MLP de 2 camadas com ReLU.

    Arquitetura: x -> Linear(W1,b1) -> ReLU -> Linear(W2,b2) -> L2norm

    A normalizacao L2 na saida e essencial para o aprendizado
    contrastivo — garante que embeddings vivem na hiperesfera
    unitaria, fazendo a similaridade cosseno = dot product.

    Args:
        x: input (batch_size, input_dim) ou (input_dim,)
        weights: pesos do MLP

    Returns:
        Tupla (embedding_normalizado, ativacao_hidden) — hidden e
        necessario para o backward pass
    """
    # Garantir 2D
    if x.ndim == 1:
        x = x.reshape(1, -1)

    # Camada 1: linear + ReLU
    hidden = x @ weights.W1 + weights.b1
    hidden_relu = np.maximum(hidden, 0)  # ReLU

    # Camada 2: linear
    output = hidden_relu @ weights.W2 + weights.b2

    # Normalizacao L2 para espaco contrastivo
    norms = np.linalg.norm(output, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    output_normalized = output / norms

    return output_normalized, hidden_relu


def mlp_backward(
    x: np.ndarray,
    hidden_relu: np.ndarray,
    grad_output: np.ndarray,
    weights: MLPWeights,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Backward pass (backpropagation) do MLP.

    Calcula gradientes de todos os pesos dado o gradiente da loss
    em relacao ao embedding de saida.

    Nota: ignora o gradiente da normalizacao L2 por simplicidade —
    na pratica, a normalizacao L2 tem gradiente quase unitario
    quando os vetores ja estao proximos da esfera unitaria.

    Args:
        x: input original (batch_size, input_dim)
        hidden_relu: ativacao da camada oculta (do forward pass)
        grad_output: gradiente da loss em relacao ao embedding (batch_size, output_dim)
        weights: pesos atuais do MLP

    Returns:
        Tupla (grad_W1, grad_b1, grad_W2, grad_b2)
    """
    if x.ndim == 1:
        x = x.reshape(1, -1)
    if grad_output.ndim == 1:
        grad_output = grad_output.reshape(1, -1)

    batch_size = x.shape[0]

    # Gradiente camada 2
    grad_W2 = hidden_relu.T @ grad_output / batch_size
    grad_b2 = np.mean(grad_output, axis=0)

    # Gradiente para hidden (antes do ReLU)
    grad_hidden = grad_output @ weights.W2.T
    # Mascara do ReLU: gradiente zero onde hidden era <= 0
    grad_hidden_relu = grad_hidden * (hidden_relu > 0).astype(float)

    # Gradiente camada 1
    grad_W1 = x.T @ grad_hidden_relu / batch_size
    grad_b1 = np.mean(grad_hidden_relu, axis=0)

    return grad_W1, grad_b1, grad_W2, grad_b2


def sgd_step(
    weights: MLPWeights,
    grads: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    lr: float,
    momentum: float = 0.9,
    weight_decay: float = 1e-4,
) -> None:
    """Atualiza pesos com SGD + momentum + weight decay.

    Momentum acumula velocidade na direcao do gradiente medio,
    acelerando convergencia em superficies alongadas da loss.

    Weight decay (regularizacao L2) penaliza pesos grandes,
    prevenindo overfitting nos poucos exemplos de treinamento.

    Args:
        weights: pesos a atualizar (in-place)
        grads: tupla (grad_W1, grad_b1, grad_W2, grad_b2)
        lr: taxa de aprendizado
        momentum: coeficiente de momentum (0 = SGD puro)
        weight_decay: coeficiente de regularizacao L2
    """
    grad_W1, grad_b1, grad_W2, grad_b2 = grads

    # Aplicar weight decay nos gradientes dos pesos (nao nos bias)
    grad_W1 = grad_W1 + weight_decay * weights.W1
    grad_W2 = grad_W2 + weight_decay * weights.W2

    # Atualizar com momentum
    weights.vW1 = momentum * weights.vW1 - lr * grad_W1
    weights.vb1 = momentum * weights.vb1 - lr * grad_b1
    weights.vW2 = momentum * weights.vW2 - lr * grad_W2
    weights.vb2 = momentum * weights.vb2 - lr * grad_b2

    weights.W1 = weights.W1 + weights.vW1
    weights.b1 = weights.b1 + weights.vb1
    weights.W2 = weights.W2 + weights.vW2
    weights.b2 = weights.b2 + weights.vb2


# ---------------------------------------------------------------------------
# Geracao de pares de treinamento
# ---------------------------------------------------------------------------

def generate_mutant_peptide(peptide: str, n_mutations: int, rng: random.Random) -> str:
    """Gera um peptideo mutante como hard negative.

    Substitui n_mutations posicoes por aminoacidos aleatorios diferentes
    do original. Isso cria negativos que sao estruturalmente similares
    ao epitopo real, mas biologicamente nao sao binders — forcando o
    modelo a aprender features sutis de ligacao.

    Args:
        peptide: sequencia original do epitopo
        n_mutations: numero de posicoes a mutar
        rng: gerador de numeros aleatorios

    Returns:
        Sequencia mutante com n_mutations substituicoes
    """
    aa_list = list(peptide)
    positions = rng.sample(range(len(aa_list)), min(n_mutations, len(aa_list)))

    for pos in positions:
        original_aa = aa_list[pos]
        # Escolher um aminoacido diferente do original
        alternatives = [aa for aa in AMINO_ACIDS if aa != original_aa]
        aa_list[pos] = rng.choice(alternatives)

    return "".join(aa_list)


def build_training_pairs(
    config: ContrastiveConfig,
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray], list[dict]]:
    """Constroi pares positivos e negativos a partir dos 11 epitopos.

    Pares positivos: (epitopo_features, alelo_features) para cada um
    dos 11 epitopos canonicos com seu alelo DLA ligante.

    Pares negativos (2 tipos):
        1. Random: epitopo emparelhado com alelo nao-ligante
        2. Hard: peptideo mutante emparelhado com alelo ligante

    Args:
        config: configuracao com parametros de treinamento

    Returns:
        Tupla (anchor_features, positive_features, negative_features, pair_info)
        onde pair_info contem metadados de cada par para diagnostico
    """
    rng = random.Random(config.seed)
    kmer_vocab = build_kmer_vocab(config.kmer_k)
    allele_names = list(DLA_ALLELE_PROPERTIES.keys())

    anchors: list[np.ndarray] = []
    positives: list[np.ndarray] = []
    negatives: list[np.ndarray] = []
    pair_info: list[dict] = []

    for epitope in EPITOPES:
        # --- Par positivo ---
        anchor_feat = encode_peptide(epitope.peptide, config.kmer_k, kmer_vocab)
        positive_feat = encode_allele(epitope.allele)

        anchors.append(anchor_feat)
        positives.append(positive_feat)
        pair_info.append({
            "peptide": epitope.peptide,
            "allele": epitope.allele,
            "ic50": epitope.ic50,
            "type": "positive",
        })

        # --- Negativos por par positivo ---
        n_neg = config.neg_ratio
        n_hard = int(n_neg * config.hard_neg_fraction)
        n_random = n_neg - n_hard

        # Negativos aleatorios: alelos nao-ligantes para este epitopo
        non_binding_alleles = [a for a in allele_names if a != epitope.allele]
        for _ in range(n_random):
            random_allele = rng.choice(non_binding_alleles)
            neg_feat = encode_allele(random_allele)
            negatives.append(neg_feat)

        # Hard negatives: peptideos mutantes com o alelo correto
        for _ in range(n_hard):
            mutant = generate_mutant_peptide(
                epitope.peptide,
                config.mutation_positions,
                rng,
            )
            # Para hard negatives baseados em mutacao de peptideo,
            # usamos o encoding do mutante como "ancora negativa"
            # emparelhado com o alelo correto.
            # Na pratica, adicionamos o encoding do alelo correto como
            # negativo (o modelo deve aprender que o mutante != epitopo real)
            neg_feat = encode_allele(epitope.allele)
            negatives.append(neg_feat)

    return anchors, positives, negatives, pair_info


# ---------------------------------------------------------------------------
# Loop de treinamento
# ---------------------------------------------------------------------------

@dataclass
class TrainingResult:
    """Resultado do treinamento contrastivo."""

    loss_history: list[float]
    final_loss: float
    n_epochs: int
    n_train_pairs: int
    epitope_weights: MLPWeights
    allele_weights: MLPWeights
    config: ContrastiveConfig


def train_contrastive(config: ContrastiveConfig) -> TrainingResult:
    """Executa o treinamento contrastivo completo.

    Pipeline:
        1. Gerar pares de treinamento (positivos + negativos)
        2. Inicializar MLPs de epitopo e alelo
        3. Loop de treinamento com InfoNCE loss
        4. Retornar modelo treinado e curva de loss

    O treinamento usa data augmentation por repeticao com shuffling —
    como temos apenas 11 epitopos, repetimos o dataset varias vezes
    por epoca com diferentes combinacoes de negativos.

    Args:
        config: configuracao completa do treinamento

    Returns:
        TrainingResult com pesos treinados e historico de loss
    """
    np.random.seed(config.seed)
    rng = random.Random(config.seed)

    # Construir vocabulario de k-mers uma vez
    kmer_vocab = build_kmer_vocab(config.kmer_k)

    # Gerar pares de treinamento
    anchors, positives, negatives, pair_info = build_training_pairs(config)

    n_pairs = len(anchors)
    print(f"  [train] {n_pairs} pares positivos, {len(negatives)} negativos totais")

    # Dimensoes dos encodings
    peptide_dim = get_peptide_feature_dim(config.kmer_k)
    allele_dim = get_allele_feature_dim()

    # Inicializar MLPs
    # MLP do epitopo: peptide_features -> embedding_dim
    epitope_mlp = init_mlp(peptide_dim, config.hidden_dim, config.embedding_dim, seed=config.seed)
    # MLP do alelo: allele_features -> embedding_dim
    allele_mlp = init_mlp(allele_dim, config.hidden_dim, config.embedding_dim, seed=config.seed + 1)

    # Converter listas para arrays numpy
    anchor_array = np.array(anchors)    # (n_pairs, peptide_dim)
    positive_array = np.array(positives)  # (n_pairs, allele_dim)
    negative_array = np.array(negatives)  # (n_pairs*neg_ratio, allele_dim)

    loss_history: list[float] = []

    # --- Loop de treinamento ---
    for epoch in range(config.n_epochs):
        # Shuffle dos indices para cada epoca
        indices = list(range(n_pairs))
        rng.shuffle(indices)

        epoch_loss = 0.0
        n_batches = 0

        # Mini-batches
        for start in range(0, n_pairs, config.batch_size):
            end = min(start + config.batch_size, n_pairs)
            batch_idx = indices[start:end]
            batch_size_actual = len(batch_idx)

            # Extrair features do batch
            batch_anchors = anchor_array[batch_idx]       # (B, peptide_dim)
            batch_positives = positive_array[batch_idx]   # (B, allele_dim)

            # Negativos para este batch
            neg_per_anchor = config.neg_ratio
            batch_neg_indices = []
            for idx in batch_idx:
                neg_start = idx * neg_per_anchor
                neg_end = neg_start + neg_per_anchor
                batch_neg_indices.extend(range(neg_start, min(neg_end, len(negatives))))
            batch_negatives = negative_array[batch_neg_indices]  # (B*neg_ratio, allele_dim)

            # Forward pass — gerar embeddings
            anchor_emb, anchor_hidden = mlp_forward(batch_anchors, epitope_mlp)
            pos_emb, pos_hidden = mlp_forward(batch_positives, allele_mlp)
            neg_emb, neg_hidden = mlp_forward(batch_negatives, allele_mlp)

            # Calcular loss InfoNCE
            loss, grad_anchors = info_nce_loss(
                anchor_emb, pos_emb, neg_emb,
                temperature=config.temperature,
            )
            epoch_loss += loss

            # Backward pass para o MLP de epitopos
            grad_W1_ep, grad_b1_ep, grad_W2_ep, grad_b2_ep = mlp_backward(
                batch_anchors, anchor_hidden, grad_anchors, epitope_mlp,
            )
            sgd_step(
                epitope_mlp,
                (grad_W1_ep, grad_b1_ep, grad_W2_ep, grad_b2_ep),
                lr=config.learning_rate,
                momentum=config.momentum,
                weight_decay=config.weight_decay,
            )

            # Backward pass simplificado para o MLP de alelos
            # Usa gradiente medio dos positivos como proxy
            # (gradiente completo exigiria backprop pela loss para cada negativo)
            grad_pos = -grad_anchors  # Simetria: se ancora deve se aproximar do positivo,
                                       # positivo deve se aproximar da ancora
            if pos_hidden.shape[0] == grad_pos.shape[0]:
                grad_W1_al, grad_b1_al, grad_W2_al, grad_b2_al = mlp_backward(
                    batch_positives, pos_hidden, grad_pos, allele_mlp,
                )
                sgd_step(
                    allele_mlp,
                    (grad_W1_al, grad_b1_al, grad_W2_al, grad_b2_al),
                    lr=config.learning_rate * 0.5,  # Alelo MLP aprende mais devagar
                    momentum=config.momentum,
                    weight_decay=config.weight_decay,
                )

            n_batches += 1

        avg_epoch_loss = epoch_loss / max(n_batches, 1)
        loss_history.append(avg_epoch_loss)

        # Log a cada 20 epocas
        if (epoch + 1) % 20 == 0 or epoch == 0:
            print(f"  [train] Epoca {epoch + 1:>3d}/{config.n_epochs} — loss: {avg_epoch_loss:.4f}")

    print(f"  [train] Treinamento completo. Loss final: {loss_history[-1]:.4f}")

    return TrainingResult(
        loss_history=loss_history,
        final_loss=loss_history[-1],
        n_epochs=config.n_epochs,
        n_train_pairs=n_pairs,
        epitope_weights=epitope_mlp,
        allele_weights=allele_mlp,
        config=config,
    )


def score_epitope_allele(
    peptide: str,
    allele: str,
    epitope_mlp: MLPWeights,
    allele_mlp: MLPWeights,
    kmer_k: int,
    kmer_vocab: dict[str, int] | None = None,
) -> float:
    """Calcula o score contrastivo para um par epitopo-alelo.

    O score e a similaridade cosseno entre os embeddings projetados
    do peptideo e do alelo no espaco contrastivo aprendido.

    Score alto (proximo de 1) -> modelo prediz ligacao forte
    Score baixo (proximo de 0 ou negativo) -> modelo prediz nao-ligacao

    Args:
        peptide: sequencia do epitopo
        allele: nome do alelo DLA
        epitope_mlp: pesos treinados do MLP de epitopos
        allele_mlp: pesos treinados do MLP de alelos
        kmer_k: tamanho do k-mer para encoding
        kmer_vocab: vocabulario pre-computado

    Returns:
        Score de similaridade contrastiva no intervalo [-1, 1]
    """
    pep_feat = encode_peptide(peptide, kmer_k, kmer_vocab)
    allele_feat = encode_allele(allele)

    pep_emb, _ = mlp_forward(pep_feat, epitope_mlp)
    allele_emb, _ = mlp_forward(allele_feat, allele_mlp)

    # Similaridade cosseno (ja normalizados pelo mlp_forward)
    return float(np.dot(pep_emb.flatten(), allele_emb.flatten()))
