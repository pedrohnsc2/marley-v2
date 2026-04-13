"""Configuracoes parametrizaveis para cada pipeline do Marley.

Cada pipeline tem um dataclass frozen com defaults correspondentes
ao L. infantum, mantendo retrocompatibilidade com o pipeline original.

O padrao segue o TargetConfig de aso_math/target_config.py.

Uso:
    from core.pipeline_config import VaccineConfig, load_config

    # Defaults (L. infantum)
    config = VaccineConfig()

    # Custom
    config = VaccineConfig(organism="trypanosoma_cruzi", conservation_threshold=0.7)

    # De arquivo JSON
    config = load_config("vaccine", "configs/presets/vaccine/leishmania_infantum.json")
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, fields, asdict
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Base mixin
# ---------------------------------------------------------------------------

class ConfigMixin:
    """Metodos comuns para todas as configs de pipeline."""

    def to_dict(self) -> dict[str, Any]:
        """Serializa a config para dict (compativel com JSON)."""
        result = {}
        for f in fields(self):  # type: ignore[arg-type]
            val = getattr(self, f.name)
            if isinstance(val, Path):
                val = str(val) if val else None
            result[f.name] = val
        return result

    def to_json(self, path: Path) -> None:
        """Salva a config como JSON."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2, ensure_ascii=False)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Any:
        """Constroi a config de um dict, ignorando chaves desconhecidas."""
        valid_fields = {f.name for f in fields(cls)}  # type: ignore[arg-type]
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered)

    @classmethod
    def from_json(cls, path: Path) -> Any:
        """Carrega config de um arquivo JSON."""
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))


# ---------------------------------------------------------------------------
# Vaccine Pipeline Config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VaccineConfig(ConfigMixin):
    """Configuracao para o pipeline de vacina reversa.

    Defaults correspondem ao pipeline original de L. infantum.
    """

    # --- Organismo alvo ---
    organism: str = "leishmania_infantum"
    species_name: str = "Leishmania infantum"
    proteome_source: str = "tritrypdb"
    proteome_id: str = "JPCM5"  # strain

    # --- Pesos de scoring ---
    conservation_weight: float = 0.4
    immunogenicity_weight: float = 0.6
    conservation_threshold: float = 0.6

    # --- Alelos MHC (canino por default) ---
    mhc_alleles: tuple[str, ...] = (
        "HLA-A*02:01",
        "HLA-A*24:02",
        "HLA-B*07:02",
        "HLA-B*08:01",
        "HLA-B*35:01",
        "HLA-B*44:02",
    )

    # --- Controle de execucao ---
    force_rerun: bool = False
    skip_fetch: bool = False
    dry_run: bool = False

    # --- Output ---
    output_dir: Path = field(default_factory=lambda: Path("results"))


# ---------------------------------------------------------------------------
# Drug Target Pipeline Config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DrugConfig(ConfigMixin):
    """Configuracao para o pipeline de drug target discovery.

    Pesos de druggability correspondem aos definidos em core/models.py.
    """

    # --- Organismo ---
    organism: str = "leishmania_infantum"
    species_name: str = "Leishmania infantum"

    # --- Pesos de druggability ---
    divergence_weight: float = 0.40
    active_site_weight: float = 0.35
    essentiality_weight: float = 0.25

    # --- Organismo hospedeiro para comparacao ---
    host_organism: str = "homo_sapiens"

    # --- Filtros ---
    top_n_targets: int = 20

    # --- Controle de execucao ---
    priority_only: bool = False
    force_rerun: bool = False
    dry_run: bool = False

    # --- Output ---
    output_dir: Path = field(default_factory=lambda: Path("results"))


# ---------------------------------------------------------------------------
# Docking Pipeline Config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DockingConfig(ConfigMixin):
    """Configuracao para o pipeline de molecular docking.

    Pesos de composite score correspondem aos de core/models.py.
    """

    # --- Organismo ---
    organism: str = "leishmania_infantum"
    species_name: str = "Leishmania infantum"

    # --- Alvos ---
    top_n_targets: int = 5
    target_gene_ids: tuple[str, ...] = ()  # vazio = usar top N do drug pipeline

    # --- Compostos ---
    compound_library: str = "default"  # "default" | "chembl" | "repurposing"
    max_compounds_per_target: int = 100

    # --- Parametros de docking ---
    exhaustiveness: int = 16
    blind_docking: bool = False

    # --- Pesos de scoring ---
    affinity_weight: float = 0.50
    admet_weight: float = 0.25
    repurposing_weight: float = 0.15
    selectivity_weight: float = 0.10

    # --- ADMET ---
    predict_admet: bool = False  # requer pkCSM API

    # --- Controle de execucao ---
    force_rerun: bool = False
    dry_run: bool = False

    # --- Output ---
    output_dir: Path = field(default_factory=lambda: Path("results"))


# ---------------------------------------------------------------------------
# RNA Entropy Pipeline Config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RnaConfig(ConfigMixin):
    """Configuracao para o pipeline de analise de entropia RNA.

    Pesos de information score correspondem aos de core/models.py.
    """

    # --- Organismo ---
    organism: str = "leishmania_infantum"
    species_name: str = "Leishmania infantum"

    # --- Organismo hospedeiro ---
    host_organism: str = "homo_sapiens"

    # --- SL RNA ---
    sl_sequence: str = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"

    # --- Pesos de information score ---
    entropy_delta_weight: float = 0.35
    conservation_weight: float = 0.25
    codon_bias_weight: float = 0.20
    sl_rna_weight: float = 0.10
    structure_weight: float = 0.10

    # --- Filtros ---
    top_n_targets: int | None = None  # None = all

    # --- Controle de execucao ---
    skip_fetch: bool = False
    priority_only: bool = False
    force_rerun: bool = False
    dry_run: bool = False

    # --- Output ---
    output_dir: Path = field(default_factory=lambda: Path("results"))


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_CONFIG_CLASSES: dict[str, type] = {
    "vaccine": VaccineConfig,
    "drug": DrugConfig,
    "docking": DockingConfig,
    "rna": RnaConfig,
}


def load_config(pipeline: str, path: str | Path) -> ConfigMixin:
    """Carrega config de um arquivo JSON para o pipeline especificado.

    Args:
        pipeline: ID do pipeline ("vaccine", "drug", "docking", "rna").
        path: Caminho para o arquivo JSON.

    Returns:
        Instancia da config correspondente.

    Raises:
        ValueError: Se o pipeline nao tiver config registrada.
    """
    if pipeline not in _CONFIG_CLASSES:
        available = ", ".join(sorted(_CONFIG_CLASSES.keys()))
        raise ValueError(
            f"Pipeline sem config registrada: '{pipeline}'. Disponiveis: {available}"
        )
    cls = _CONFIG_CLASSES[pipeline]
    return cls.from_json(Path(path))


def get_default_config(pipeline: str) -> ConfigMixin:
    """Retorna a config default (L. infantum) para um pipeline."""
    if pipeline not in _CONFIG_CLASSES:
        # aso_math usa TargetConfig diretamente
        if pipeline == "aso_math":
            from aso_math.target_config import TargetConfig
            return TargetConfig()  # type: ignore[return-value]
        available = ", ".join(sorted(_CONFIG_CLASSES.keys()))
        raise ValueError(
            f"Pipeline sem config registrada: '{pipeline}'. Disponiveis: {available}"
        )
    return _CONFIG_CLASSES[pipeline]()
