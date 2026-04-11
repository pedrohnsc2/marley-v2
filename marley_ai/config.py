"""Configuracao centralizada do track de IA/ML (marley_ai).

Define a dataclass imutavel AIModuleConfig que todo modulo recebe.
Segue o padrao frozen dataclass do aso_math para consistencia.

Cada modulo pode estender AIModuleConfig com campos especificos,
mas deve manter os campos base para integracao com o orquestrador.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final


# ---------------------------------------------------------------------------
# Caminhos do projeto
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
AI_ROOT: Final[Path] = Path(__file__).resolve().parent
RESULTS_DIR: Final[Path] = AI_ROOT / "results"

# Caminhos para resultados dos tracks anteriores (inputs do AI track)
ASO_MATH_RESULTS: Final[Path] = PROJECT_ROOT / "aso_math" / "results"
ASO_DELIVERY_RESULTS: Final[Path] = PROJECT_ROOT / "aso_delivery" / "results"
VACCINE_RESULTS: Final[Path] = PROJECT_ROOT / "vaccine_platforms" / "results"

# Cache de modelos e embeddings
MODELS_DIR: Final[Path] = AI_ROOT / ".models"
CACHE_DIR: Final[Path] = AI_ROOT / ".cache"

# ---------------------------------------------------------------------------
# Registro de modulos AI (slug -> nome do pacote)
# ---------------------------------------------------------------------------

MODULES: Final[dict[str, str]] = {
    "01": "01_rag",
    "02": "02_leish_kg",
    "03": "03_leish_esm",
    "04": "04_rna_fm",
    "05": "05_rosettafold",
    "06": "06_evodiff",
    "07": "07_contrastive",
    "08": "08_rl_ppo",
    "09": "09_sae",
    "10": "10_digital_twin",
    "11": "11_scientist",
}


# ---------------------------------------------------------------------------
# Configuracao base para todos os modulos
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AIModuleConfig:
    """Configuracao imutavel base para qualquer modulo do track AI.

    Campos:
        module_slug: Identificador curto do modulo (ex: "03_leish_esm").
        module_name: Nome legivel (ex: "ESM-2 Protein Embeddings").
        version: Versao do modulo.
        device: Dispositivo de computacao — "mps" para Apple Silicon, "cpu" fallback.
        dtype: Precisao numerica — "float32" para estabilidade, "float16" para velocidade.
        max_memory_gb: Limite de memoria para o modulo (evita OOM em M1/M2).
        seed: Semente global para reproducibilidade.
        output_dir: Diretorio para resultados do modulo.
        cache_dir: Diretorio para cache de modelos e embeddings.
        aso_math_results: Caminho para resultados do track aso_math.
        aso_delivery_results: Caminho para resultados do track aso_delivery.
        vaccine_results: Caminho para resultados do track vaccine_platforms.
    """

    # --- Identidade do modulo ---
    module_slug: str
    module_name: str
    version: str = "0.1.0"

    # --- Hardware ---
    device: str = "mps"        # Apple Silicon M1/M2/M3 — fallback para "cpu"
    dtype: str = "float32"     # float32 para estabilidade numerica
    max_memory_gb: float = 8.0  # Limite conservador para MacBook
    seed: int = 42

    # --- Caminhos de saida ---
    output_dir: Path = field(default_factory=lambda: RESULTS_DIR)
    cache_dir: Path = field(default_factory=lambda: CACHE_DIR)

    # --- Caminhos para resultados de tracks anteriores ---
    aso_math_results: Path = field(default_factory=lambda: ASO_MATH_RESULTS)
    aso_delivery_results: Path = field(default_factory=lambda: ASO_DELIVERY_RESULTS)
    vaccine_results: Path = field(default_factory=lambda: VACCINE_RESULTS)

    def __post_init__(self) -> None:
        """Valida device e cria diretorios se necessario."""
        valid_devices = {"mps", "cpu", "cuda"}
        if self.device not in valid_devices:
            raise ValueError(
                f"Device invalido: '{self.device}'. Validos: {valid_devices}"
            )
        valid_dtypes = {"float16", "float32", "bfloat16"}
        if self.dtype not in valid_dtypes:
            raise ValueError(
                f"dtype invalido: '{self.dtype}'. Validos: {valid_dtypes}"
            )


def detect_device() -> str:
    """Detecta o melhor dispositivo disponivel.

    Ordem de preferencia: MPS (Apple Silicon) > CUDA > CPU.
    Nao importa torch diretamente para evitar dependencia no scaffold.
    """
    try:
        import torch
        if torch.backends.mps.is_available():
            return "mps"
        if torch.cuda.is_available():
            return "cuda"
    except ImportError:
        pass
    return "cpu"
