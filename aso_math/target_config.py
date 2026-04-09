"""Configuracao parametrizavel para execucao do pipeline em qualquer organismo.

Substitui os valores hardcoded de config.py. Um TargetConfig e construido
uma vez no inicio da execucao e passado a cada modulo.

Retrocompativel: se nenhum TargetConfig for fornecido, os modulos usam
os defaults de L. infantum (comportamento identico ao anterior).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class TargetConfig:
    """Configuracao imutavel para uma execucao de validacao ASO.

    Cada campo tem um default correspondente ao L. infantum (MRL-ASO-001)
    para manter retrocompatibilidade com o pipeline original.
    """

    # --- Identidade do organismo ---
    organism_slug: str = "leishmania_infantum"
    species_name: str = "Leishmania infantum"
    taxonomic_group: str = "kinetoplastida"  # kinetoplastida, nematoda, platyhelminthes
    disease_name: str = "leishmaniose visceral"

    # --- Sequencia alvo (SL RNA) ---
    sl_sequence: str = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"
    sl_length: int = 0  # auto-calculado no __post_init__

    # --- ASO candidato ---
    aso_sequence: str = ""  # vazio = auto-design pelo pipeline
    aso_name: str = ""
    aso_target_start: int = 0
    aso_target_end: int = 0

    # --- Biologia especifica do organismo ---
    mutation_rate: float = 2.0e-9        # por nt por replicacao
    generation_time_hours: float = 12.0  # tempo de divisao
    sl_copy_number: int = 150            # copias no tandem array
    dg_functional_threshold: float = -15.0  # kcal/mol

    # --- Parametros de design ---
    length_scan_min: int = 18
    length_scan_max: int = 27
    aso_concentration: float = 250e-9  # 250 nM

    # --- Configuracao LNA ---
    lna_5prime: int = 4
    lna_3prime: int = 4

    # --- Conservacao cross-species ---
    related_sl_sequences: dict[str, str] = field(default_factory=dict)
    divergence_time_mya: float = 350.0
    neutral_rate_per_site_per_year: float = 3.0e-9

    # --- Caminhos de arquivos ---
    host_transcriptome_path: Path | None = None
    output_dir: Path = field(default_factory=lambda: Path("aso_math/results"))

    # --- Valores conhecidos para validacao (opcional) ---
    known_tm: float | None = None
    known_dg: float | None = None
    known_gc: float | None = None

    def __post_init__(self) -> None:
        """Auto-calcula campos derivados."""
        if self.sl_length == 0:
            object.__setattr__(self, "sl_length", len(self.sl_sequence))
        if self.aso_sequence and self.aso_target_end == 0:
            object.__setattr__(self, "aso_target_end", self.sl_length)
