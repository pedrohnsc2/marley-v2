"""Encoder customizado de RNA baseado em Transformer (PyTorch).

Implementa um modelo de linguagem para RNA com treinamento auto-supervisionado
via predicao de nucleotideo mascarado (Masked Nucleotide Modeling — MNM),
inspirado na abordagem BERT/RNA-FM mas treinavel from scratch em dataset
pequeno usando Apple MPS.

Arquitetura:
    - Embedding: one-hot (A=0, U=1, G=2, C=3, MASK=4) + positional encoding
    - Encoder: Transformer encoder (3 camadas, 4 cabecas, dim=128)
    - Saida: embeddings por nucleotideo + embedding de sequencia (mean pooling)
    - Treinamento: mascarar 15% das posicoes, predizer nucleotideo original

O modelo gera representacoes vetoriais densas de sequencias de RNA que
capturam padroes estruturais e composicionais. Usado para comparar SL RNA
de Leishmania contra RNAs humanos/caninos no espaco de embeddings.

Ref: Devlin J et al. (2019) BERT — pre-training of deep bidirectional transformers
Ref: Chen J et al. (2022) arXiv:2204.00300 — RNA-FM
"""

from __future__ import annotations

import math
import random
from typing import Any

import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import ENCODER_PARAMS


# ---------------------------------------------------------------------------
# Mapeamento de nucleotideos para indices
# ---------------------------------------------------------------------------

NUC_TO_IDX: dict[str, int] = {"A": 0, "U": 1, "G": 2, "C": 3}
MASK_IDX: int = 4
PAD_IDX: int = 5
IDX_TO_NUC: dict[int, str] = {v: k for k, v in NUC_TO_IDX.items()}
VOCAB_SIZE: int = 6  # A, U, G, C, MASK, PAD


# ---------------------------------------------------------------------------
# Codificacao posicional sinusoidal
# Ref: Vaswani A et al. (2017) Attention is All You Need
# ---------------------------------------------------------------------------

class SinusoidalPositionalEncoding(nn.Module):
    """Encoding posicional sinusoidal para sequencias de RNA.

    Injeta informacao sobre a posicao absoluta de cada nucleotideo
    na sequencia, permitindo que o transformer distinga posicoes
    sem recorrencia ou convolucao.
    """

    def __init__(self, embed_dim: int, max_len: int = 512) -> None:
        super().__init__()
        # Calcula a matrix de encoding posicional uma unica vez
        pe = torch.zeros(max_len, embed_dim)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, embed_dim, 2, dtype=torch.float)
            * (-math.log(10000.0) / embed_dim)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)  # (1, max_len, embed_dim)
        # Registra como buffer (nao e parametro treinavel)
        self.register_buffer("pe", pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Soma encoding posicional ao embedding de entrada.

        Args:
            x: Tensor (batch, seq_len, embed_dim).

        Returns:
            Tensor com mesma shape, acrescido de informacao posicional.
        """
        return x + self.pe[:, : x.size(1), :]


# ---------------------------------------------------------------------------
# Encoder de RNA baseado em Transformer
# ---------------------------------------------------------------------------

class RNAEncoder(nn.Module):
    """Encoder Transformer customizado para sequencias de RNA.

    Gera embeddings densos por nucleotideo e por sequencia a partir de
    sequencias de RNA codificadas como indices inteiros.

    Atributos:
        embedding: Camada de embedding (vocab_size -> embed_dim)
        pos_encoder: Encoding posicional sinusoidal
        transformer_encoder: Stack de camadas TransformerEncoder
        classifier: Projecao para predicao de nucleotideo mascarado
    """

    def __init__(
        self,
        vocab_size: int = VOCAB_SIZE,
        embed_dim: int = ENCODER_PARAMS["embed_dim"],
        num_heads: int = ENCODER_PARAMS["num_heads"],
        num_layers: int = ENCODER_PARAMS["num_layers"],
        feedforward_dim: int = ENCODER_PARAMS["feedforward_dim"],
        dropout: float = ENCODER_PARAMS["dropout"],
        max_seq_len: int = ENCODER_PARAMS["max_seq_len"],
    ) -> None:
        super().__init__()
        self.embed_dim = embed_dim

        # Embedding de nucleotideos (indices -> vetores densos)
        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=PAD_IDX)

        # Encoding posicional sinusoidal
        self.pos_encoder = SinusoidalPositionalEncoding(embed_dim, max_seq_len)

        # Dropout pos-embedding (regularizacao)
        self.dropout = nn.Dropout(dropout)

        # Stack de camadas Transformer Encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embed_dim,
            nhead=num_heads,
            dim_feedforward=feedforward_dim,
            dropout=dropout,
            batch_first=True,
            norm_first=True,  # Pre-norm (mais estavel para modelos pequenos)
        )
        self.transformer_encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_layers,
            enable_nested_tensor=False,  # Evita warning com norm_first=True
        )

        # Camada de normalizacao final
        self.layer_norm = nn.LayerNorm(embed_dim)

        # Classificador para predicao de nucleotideo mascarado
        # Projeta embedding -> 4 classes (A, U, G, C)
        self.classifier = nn.Sequential(
            nn.Linear(embed_dim, embed_dim),
            nn.GELU(),
            nn.LayerNorm(embed_dim),
            nn.Linear(embed_dim, 4),  # 4 nucleotideos canonicos
        )

    def forward(
        self,
        input_ids: torch.Tensor,
        padding_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        """Forward pass: gera embeddings e logits de classificacao.

        Args:
            input_ids: Tensor (batch, seq_len) com indices de nucleotideos.
            padding_mask: Tensor (batch, seq_len) booleano, True = posicao mascarada
                          (padding). Usado para ignorar posicoes de pad no attention.

        Returns:
            Dicionario com:
                - "token_embeddings": (batch, seq_len, embed_dim) — por nucleotideo
                - "sequence_embedding": (batch, embed_dim) — media da sequencia
                - "logits": (batch, seq_len, 4) — predicao de nucleotideo
        """
        # Embedding + posicional
        x = self.embedding(input_ids)  # (batch, seq_len, embed_dim)
        x = self.pos_encoder(x)
        x = self.dropout(x)

        # Transformer encoder (auto-atencao bidirecional)
        # src_key_padding_mask: True indica posicoes a ignorar
        x = self.transformer_encoder(x, src_key_padding_mask=padding_mask)
        x = self.layer_norm(x)

        # Embedding de sequencia = media sobre posicoes nao-padding
        if padding_mask is not None:
            # Inverte mascara: True = posicao valida
            valid_mask = ~padding_mask
            mask_expanded = valid_mask.unsqueeze(-1).float()  # (batch, seq_len, 1)
            seq_emb = (x * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1).clamp(min=1)
        else:
            seq_emb = x.mean(dim=1)

        # Logits para predicao de nucleotideo mascarado
        logits = self.classifier(x)  # (batch, seq_len, 4)

        return {
            "token_embeddings": x,
            "sequence_embedding": seq_emb,
            "logits": logits,
        }


# ---------------------------------------------------------------------------
# Funcoes de preparacao de dados
# ---------------------------------------------------------------------------

def encode_sequence(sequence: str, max_len: int | None = None) -> list[int]:
    """Converte sequencia de RNA para lista de indices inteiros.

    Converte T -> U (DNA -> RNA) para unificar representacao.

    Args:
        sequence: Sequencia de RNA/DNA (ex: "AACTAACG...").
        max_len: Comprimento maximo. Se None, usa tamanho da sequencia.

    Returns:
        Lista de indices inteiros (A=0, U=1, G=2, C=3).
    """
    # Normaliza: converte T para U (DNA -> RNA)
    seq = sequence.upper().replace("T", "U")
    indices = []
    for nuc in seq:
        if nuc in NUC_TO_IDX:
            indices.append(NUC_TO_IDX[nuc])
        # Ignora caracteres desconhecidos (N, etc.)

    if max_len is not None and len(indices) < max_len:
        # Padding com PAD_IDX
        indices.extend([PAD_IDX] * (max_len - len(indices)))
    elif max_len is not None:
        indices = indices[:max_len]

    return indices


def create_masked_input(
    input_ids: torch.Tensor,
    mask_fraction: float = ENCODER_PARAMS["mask_fraction"],
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Cria input mascarado para treinamento auto-supervisionado (MNM).

    Mascara uma fracao aleatoria de posicoes nao-padding e retorna:
    - Input com posicoes substituidas por MASK_IDX
    - Labels com nucleotideo original nas posicoes mascaradas (-100 em outras)
    - Mascara booleana indicando posicoes mascaradas

    Args:
        input_ids: Tensor (batch, seq_len) com indices originais.
        mask_fraction: Fracao de posicoes a mascarar (default: 0.15).

    Returns:
        Tupla (masked_input, labels, mask_positions).
    """
    masked_input = input_ids.clone()
    labels = torch.full_like(input_ids, -100)  # -100 = ignorar no loss

    # Identifica posicoes validas (nao-padding)
    valid_mask = (input_ids != PAD_IDX) & (input_ids != MASK_IDX)

    for i in range(input_ids.size(0)):
        # Posicoes validas nesta sequencia
        valid_positions = valid_mask[i].nonzero(as_tuple=True)[0].tolist()
        if not valid_positions:
            continue

        # Seleciona posicoes para mascarar
        n_mask = max(1, int(len(valid_positions) * mask_fraction))
        mask_positions = random.sample(valid_positions, min(n_mask, len(valid_positions)))

        for pos in mask_positions:
            labels[i, pos] = input_ids[i, pos]
            masked_input[i, pos] = MASK_IDX

    mask_bool = labels != -100
    return masked_input, labels, mask_bool


# ---------------------------------------------------------------------------
# Funcao de treinamento
# ---------------------------------------------------------------------------

def train_rna_encoder(
    sequences: dict[str, str],
    device: str = "mps",
    num_epochs: int = ENCODER_PARAMS["num_epochs"],
    learning_rate: float = ENCODER_PARAMS["learning_rate"],
    verbose: bool = True,
) -> tuple[RNAEncoder, list[float]]:
    """Treina o encoder de RNA com Masked Nucleotide Modeling.

    Gera dados de treinamento a partir das sequencias fornecidas,
    incluindo variantes com mutacoes pontuais para aumentar o dataset.
    Treina usando cross-entropy loss nas posicoes mascaradas.

    Args:
        sequences: Dicionario {nome: sequencia_rna}.
        device: Dispositivo de computo ("mps", "cuda", "cpu").
        num_epochs: Numero de epocas de treinamento.
        learning_rate: Taxa de aprendizado para Adam.
        verbose: Se True, imprime progresso a cada 10 epocas.

    Returns:
        Tupla (modelo_treinado, lista_de_losses_por_epoca).
    """
    # --- Preparacao do dataset com data augmentation ---
    # Gera variantes com mutacoes pontuais para aumentar diversidade
    all_sequences = list(sequences.values())
    augmented = _augment_sequences(all_sequences, n_variants=5)
    all_sequences.extend(augmented)

    # Calcula comprimento maximo para padding uniforme
    max_len = max(len(s.replace("T", "U")) for s in all_sequences)

    # Codifica todas as sequencias
    encoded = [encode_sequence(seq, max_len=max_len) for seq in all_sequences]
    data_tensor = torch.tensor(encoded, dtype=torch.long)

    # --- Inicializa modelo ---
    model = RNAEncoder(max_seq_len=max(max_len, 256))
    model = model.to(device)
    model.train()

    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=num_epochs, eta_min=learning_rate * 0.01
    )

    losses: list[float] = []
    batch_size = min(ENCODER_PARAMS["batch_size"], len(data_tensor))

    if verbose:
        print(f"  [model] Treinando encoder de RNA:")
        print(f"    Dataset: {len(data_tensor)} sequencias, max_len={max_len}")
        print(f"    Device: {device}, Epochs: {num_epochs}")

    for epoch in range(num_epochs):
        epoch_loss = 0.0
        n_batches = 0

        # Shuffle dos dados a cada epoca
        indices = torch.randperm(len(data_tensor))

        for start in range(0, len(data_tensor), batch_size):
            batch_idx = indices[start : start + batch_size]
            batch = data_tensor[batch_idx].to(device)

            # Cria mascara de padding
            padding_mask = batch == PAD_IDX

            # Cria input mascarado para MNM
            masked_input, labels, _ = create_masked_input(batch)
            masked_input = masked_input.to(device)
            labels = labels.to(device)

            # Forward pass
            outputs = model(masked_input, padding_mask=padding_mask)
            logits = outputs["logits"]  # (batch, seq_len, 4)

            # Cross-entropy loss apenas nas posicoes mascaradas
            # Reshape para (batch*seq_len, 4) e (batch*seq_len,)
            loss = F.cross_entropy(
                logits.reshape(-1, 4),
                labels.reshape(-1),
                ignore_index=-100,
            )

            # Backward + update
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        scheduler.step()
        avg_loss = epoch_loss / max(n_batches, 1)
        losses.append(avg_loss)

        if verbose and (epoch + 1) % 20 == 0:
            print(f"    Epoch {epoch + 1:3d}/{num_epochs}: loss = {avg_loss:.4f}")

    model.eval()
    if verbose:
        print(f"    Loss final: {losses[-1]:.4f}")

    return model, losses


def _augment_sequences(
    sequences: list[str],
    n_variants: int = 5,
) -> list[str]:
    """Gera variantes com mutacoes pontuais para data augmentation.

    Para cada sequencia original, gera n_variants com 1-3 mutacoes
    pontuais aleatorias. Isso aumenta a diversidade do dataset de
    treinamento sem introduzir vieses sistematicos.

    Args:
        sequences: Lista de sequencias originais.
        n_variants: Numero de variantes por sequencia.

    Returns:
        Lista de sequencias mutantes.
    """
    nucleotides = ["A", "U", "G", "C"]
    variants = []

    for seq in sequences:
        seq_rna = seq.upper().replace("T", "U")
        for _ in range(n_variants):
            mutant = list(seq_rna)
            # 1-3 mutacoes pontuais aleatorias
            n_mutations = random.randint(1, min(3, len(mutant)))
            positions = random.sample(range(len(mutant)), n_mutations)
            for pos in positions:
                original = mutant[pos]
                # Substitui por nucleotideo diferente
                alternatives = [n for n in nucleotides if n != original]
                mutant[pos] = random.choice(alternatives)
            variants.append("".join(mutant))

    return variants


# ---------------------------------------------------------------------------
# Geracao de embeddings
# ---------------------------------------------------------------------------

@torch.no_grad()
def get_sequence_embeddings(
    model: RNAEncoder,
    sequences: dict[str, str],
    device: str = "mps",
) -> dict[str, torch.Tensor]:
    """Gera embeddings de sequencia para um conjunto de RNAs.

    Executa forward pass no modelo treinado e retorna o embedding
    de sequencia (mean pooling) para cada RNA.

    Args:
        model: Encoder de RNA treinado.
        sequences: Dicionario {nome: sequencia_rna}.
        device: Dispositivo de computo.

    Returns:
        Dicionario {nome: tensor_embedding (embed_dim,)}.
    """
    model.eval()
    embeddings = {}

    # Calcula max_len para padding uniforme
    max_len = max(len(s.replace("T", "U")) for s in sequences.values())

    for name, seq in sequences.items():
        encoded = encode_sequence(seq, max_len=max_len)
        input_ids = torch.tensor([encoded], dtype=torch.long, device=device)
        padding_mask = input_ids == PAD_IDX

        outputs = model(input_ids, padding_mask=padding_mask)
        embeddings[name] = outputs["sequence_embedding"].squeeze(0).cpu()

    return embeddings
