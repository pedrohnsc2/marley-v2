"""Configuracao do modulo 06_evodiff — Geracao por difusao discreta.

Modelo de difusao treinado do zero em PyTorch para gerar variantes
de ASOs (nucleotideos) e epitopos (aminoacidos) com propriedades
otimizadas para o projeto Marley.

Ref: Alamdari S et al. (2023) EvoDiff — NeurIPS 2023
     Austin J et al. (2021) Structured Denoising Diffusion — NeurIPS 2021
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Alfabetos biologicos
# ---------------------------------------------------------------------------

# Nucleotideos para variantes de ASO
NUCLEOTIDE_ALPHABET: Final[str] = "ACGT"

# 20 aminoacidos canonicos para variantes de epitopos
AMINO_ACID_ALPHABET: Final[str] = "ACDEFGHIKLMNPQRSTVWY"

# Token especial de mascara (indice = len(alphabet))
MASK_TOKEN: Final[str] = "<MASK>"


# ---------------------------------------------------------------------------
# Configuracao principal
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EvoDiffConfig(AIModuleConfig):
    """Configuracao para o modelo de difusao discreta.

    O modelo e treinado do zero em sequencias biologicas curtas, usando
    um transformer pequeno que cabe na memoria do M3 Pro (36GB).

    Dois modos de operacao:
        - "nucleotide": gera variantes de ASO (~24 nt, alfabeto ACGT)
        - "amino_acid": gera variantes de epitopos (~9 aa, 20 AAs)
    """

    # --- Difusao ---
    n_diffusion_steps: int = 100          # T: numero de passos de difusao
    beta_start: float = 0.001             # Beta inicial (schedule linear)
    beta_end: float = 0.1                 # Beta final

    # --- Modelo (transformer denoiser) ---
    d_model: int = 128                    # Dimensao do embedding
    n_heads: int = 4                      # Cabecas de atencao (128/4 = 32 dim/head)
    n_layers: int = 2                     # Camadas do transformer encoder
    dropout: float = 0.1                  # Dropout para regularizacao

    # --- Treinamento ---
    n_epochs: int = 200                   # Epocas de treinamento
    batch_size: int = 32                  # Tamanho do batch
    learning_rate: float = 3e-4           # Adam learning rate
    seed: int = 42                        # Semente para reproducibilidade

    # --- Geracao ---
    n_samples: int = 50                   # Sequencias geradas por execucao
    n_denoising_steps: int = 100          # Passos de denoising na geracao

    # --- Sequencia ---
    sequence_type: str = "nucleotide"     # "nucleotide" ou "amino_acid"
    min_length: int = 18                  # Comprimento minimo (ASO)
    max_length: int = 27                  # Comprimento maximo (ASO)

    # --- Filtragem de ASO ---
    gc_min: float = 0.25                  # GC content minimo
    gc_max: float = 0.65                  # GC content maximo
    homopolymer_max: int = 4             # Maximo de bases repetidas consecutivas

    # --- Filtragem de epitopos ---
    min_hydrophobicity: float = -2.0      # Kyte-Doolittle medio minimo
    max_hydrophobicity: float = 2.0       # Kyte-Doolittle medio maximo

    # --- Ranking ---
    top_k: int = 10                       # Top-K candidatos no relatorio final

    # --- Caminhos ---
    candidates_dir: Path = field(
        default_factory=lambda: AI_ROOT / "06_evodiff" / "candidates"
    )
    results_dir: Path = field(
        default_factory=lambda: AI_ROOT / "06_evodiff" / "results"
    )

    def get_alphabet(self) -> str:
        """Retorna o alfabeto correto para o tipo de sequencia."""
        if self.sequence_type == "nucleotide":
            return NUCLEOTIDE_ALPHABET
        elif self.sequence_type == "amino_acid":
            return AMINO_ACID_ALPHABET
        else:
            raise ValueError(
                f"sequence_type invalido: '{self.sequence_type}'. "
                f"Validos: 'nucleotide', 'amino_acid'"
            )

    def get_vocab_size(self) -> int:
        """Tamanho do vocabulario: alfabeto + token de mascara."""
        return len(self.get_alphabet()) + 1  # +1 para MASK_TOKEN
