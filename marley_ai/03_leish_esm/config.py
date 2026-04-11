"""Configuracao do modulo 03_leish_esm — ESM-2 Protein Embeddings.

Define parametros para inferencia com o modelo ESM-2 (650M params)
em Apple Silicon M3 Pro com 36GB RAM. O modelo cabe facilmente em
MPS com batch_size=4.

O ESM-2 t33 (33 camadas transformer, dim=1280) e o maior modelo
que roda confortavelmente em 36GB sem quantizacao. Para maquinas
com menos RAM, usar esm2_t12_35M_UR50D (12 camadas, dim=480).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from marley_ai.config import AIModuleConfig, AI_ROOT, MODELS_DIR


# ---------------------------------------------------------------------------
# Modelo ESM-2 — variantes disponiveis
# ---------------------------------------------------------------------------
# esm2_t6_8M_UR50D    —   8M params, dim= 320, 6  camadas (tiny)
# esm2_t12_35M_UR50D  —  35M params, dim= 480, 12 camadas (small)
# esm2_t30_150M_UR50D — 150M params, dim= 640, 30 camadas (medium)
# esm2_t33_650M_UR50D — 650M params, dim=1280, 33 camadas (large)  <-- usado
# esm2_t36_3B_UR50D   —   3B params, dim=2560, 36 camadas (xl)
# esm2_t48_15B_UR50D  —  15B params, dim=5120, 48 camadas (xxl)


@dataclass(frozen=True)
class ESMConfig(AIModuleConfig):
    """Configuracao especifica para embeddings ESM-2.

    Parametros otimizados para Apple M3 Pro (36GB RAM).
    O modelo de 650M parametros ocupa ~2.5GB em float32,
    deixando ampla margem para batches e embeddings em memoria.
    """

    # --- Modelo ---
    model_name: str = "facebook/esm2_t33_650M_UR50D"
    esm_model: str = "esm2_t33_650M_UR50D"  # nome para fair-esm
    embedding_dim: int = 1280                 # hidden size do t33
    num_layers: int = 33                      # camadas transformer
    layer: int = 33                           # camada para extrair embeddings

    # --- Fallback para maquinas com pouca RAM ---
    fallback_model_name: str = "facebook/esm2_t12_35M_UR50D"
    fallback_esm_model: str = "esm2_t12_35M_UR50D"
    fallback_embedding_dim: int = 480
    fallback_num_layers: int = 12

    # --- Inferencia ---
    batch_size: int = 4              # conservador para 36GB RAM
    max_seq_length: int = 1024       # truncar proteinas longas
    pool_strategy: str = "mean"      # "mean", "cls", ou "per_residue"

    # --- Caminhos ---
    models_dir: Path = field(
        default_factory=lambda: MODELS_DIR
    )
    embeddings_dir: Path = field(
        default_factory=lambda: AI_ROOT / "03_leish_esm" / "results" / "embeddings"
    )
    results_dir: Path = field(
        default_factory=lambda: AI_ROOT / "03_leish_esm" / "results"
    )


def get_default_config() -> ESMConfig:
    """Retorna configuracao padrao para o modulo 03_leish_esm.

    Detecta o melhor device disponivel (MPS > CUDA > CPU)
    e retorna config com parametros otimizados para o hardware.
    """
    from marley_ai.config import detect_device

    device = detect_device()
    return ESMConfig(
        module_slug="03_leish_esm",
        module_name="ESM-2 Protein Embeddings",
        version="0.1.0",
        device=device,
        max_memory_gb=36.0,
    )
