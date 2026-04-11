"""Configuracao do modulo 07_contrastive — Aprendizado contrastivo.

Parametros para treinamento contrastivo InfoNCE de pares epitopo-DLA
usando numpy puro (sem PyTorch). O modelo aprende um espaco de
embeddings onde epitopos e seus alelos DLA ligantes ficam proximos.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from marley_ai.config import AIModuleConfig, AI_ROOT


@dataclass(frozen=True)
class ContrastiveConfig(AIModuleConfig):
    """Configuracao especifica para aprendizado contrastivo epitopo-MHC.

    O modelo alinha dois espacos de embeddings:
        - Encoder de epitopos (k-mer + propriedades fisicoquimicas -> vetor)
        - Encoder de alelos DLA (property-based -> vetor)

    Treinamento por InfoNCE loss com negativos hard-mined.
    Implementacao inteiramente em numpy (sem PyTorch).
    """

    # --- Arquitetura ---
    embedding_dim: int = 64         # Dimensao do espaco compartilhado
    hidden_dim: int = 128           # Dimensao da camada oculta do MLP
    temperature: float = 0.07       # Temperatura da InfoNCE loss
    kmer_k: int = 3                 # Tamanho do k-mer para encoding de sequencia
    n_physicochemical: int = 5      # Propriedades fisicoquimicas por aminoacido

    # --- Treinamento ---
    learning_rate: float = 0.01     # Taxa de aprendizado (SGD com momentum)
    momentum: float = 0.9           # Momentum do SGD
    batch_size: int = 16            # Tamanho do mini-batch
    n_epochs: int = 200             # Epocas de treinamento
    neg_ratio: int = 5              # Negativos por par positivo
    weight_decay: float = 1e-4      # Regularizacao L2

    # --- Negative mining ---
    hard_neg_fraction: float = 0.5  # Fracao de negativos que sao "hard"
    mutation_positions: int = 2     # Posicoes mutadas para hard negatives

    # --- Dados e saida ---
    training_data_path: Path = field(
        default_factory=lambda: AI_ROOT / "07_contrastive" / "data"
    )
    checkpoints_dir: Path = field(
        default_factory=lambda: AI_ROOT / "07_contrastive" / "checkpoints"
    )
