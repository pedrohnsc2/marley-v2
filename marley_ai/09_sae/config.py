"""Configuracao do modulo 09_sae — Sparse Autoencoders.

Parametros ajustados para representacoes hand-crafted de aminoacidos.
Sem ESM-2 por enquanto — usamos vetores de propriedades fisico-quimicas
com estatisticas por janela deslizante.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Dimensoes da codificacao de sequencias
# ---------------------------------------------------------------------------
# 6 propriedades por aminoacido x 4 estatisticas (mean, std, max, min)
# x 3 janelas = 72 features de entrada
NUM_AA_PROPERTIES: int = 6
NUM_STATISTICS: int = 4        # mean, std, max, min
NUM_WINDOWS: int = 3           # janelas: toda a sequencia, metade N-ter, metade C-ter
INPUT_DIM: int = NUM_AA_PROPERTIES * NUM_STATISTICS * NUM_WINDOWS  # 72


@dataclass(frozen=True)
class SAEConfig(AIModuleConfig):
    """Configuracao especifica para Sparse Autoencoders.

    O SAE expande a dimensionalidade dos vetores de propriedades
    fisico-quimicas para um espaco esparso overcomplete (4-8x),
    onde cada neuronio ativo corresponde a um conceito biologico
    interpretavel (ex: "regiao hidrofobica", "cluster aromatico").
    """

    # --- Arquitetura ---
    input_dim: int = INPUT_DIM          # 72 features de propriedades AA
    hidden_dim: int = 256               # Espaco esparso (~ 3.5x overcomplete)
    sparsity_penalty: float = 0.05      # Coeficiente L1 — moderado, esparsidade via bias adaptativo
    tied_weights: bool = True           # Decoder W_dec = W_enc^T (reduz params)
    target_sparsity: float = 0.15       # Fracao de neuronios ativos alvo (~15%)

    # --- Treinamento ---
    learning_rate: float = 5e-3
    n_epochs: int = 2000
    seed: int = 42

    # --- Interpretabilidade ---
    top_k_features: int = 10            # Features mais relevantes no relatorio
    n_control_peptides: int = 50        # Peptideos aleatorios para contraste

    # --- Caminhos ---
    features_dir: Path = field(
        default_factory=lambda: AI_ROOT / "09_sae" / "features"
    )
    checkpoints_dir: Path = field(
        default_factory=lambda: AI_ROOT / "09_sae" / "checkpoints"
    )
