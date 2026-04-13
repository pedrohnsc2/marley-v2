"""Registro central de pipelines disponiveis no Marley.

Cada pipeline e definido com seus stages, modulo de configuracao
e preset default. O registro e usado pelo CLI unificado para
descobrir e executar pipelines dinamicamente.

Uso:
    from core.registry import PIPELINES, get_pipeline

    pipe = get_pipeline("vaccine")
    for stage in pipe.stages:
        print(stage.stage_id, stage.module_path)
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class StageDefinition:
    """Definicao de um stage de pipeline (nao uma execucao)."""

    stage_id: str
    name: str
    module_path: str
    entry_function: str
    description: str = ""


@dataclass(frozen=True)
class PipelineDefinition:
    """Definicao de um pipeline completo."""

    pipeline_id: str
    display_name: str
    description: str
    stages: tuple[StageDefinition, ...] = ()
    config_module: str = ""  # dotted path to config dataclass module
    config_class: str = ""   # class name within that module
    default_preset: str = "leishmania_infantum"


# ---------------------------------------------------------------------------
# Pipeline: Vaccine (reverse vaccinology, 5 stages)
# ---------------------------------------------------------------------------

VACCINE_STAGES = (
    StageDefinition(
        stage_id="01_fetch_genome",
        name="Fetch Genome",
        module_path="pipeline.01_fetch_genome",
        entry_function="fetch_genome",
        description="Download proteome from TriTrypDB",
    ),
    StageDefinition(
        stage_id="02_filter_surface",
        name="Filter Surface Proteins",
        module_path="pipeline.02_filter_surface",
        entry_function="filter_surface_proteins",
        description="SignalP 6.0 signal peptide screen",
    ),
    StageDefinition(
        stage_id="03_conservation",
        name="Conservation Analysis",
        module_path="pipeline.03_conservation",
        entry_function="analyze_conservation",
        description="BLAST against Leishmania strains",
    ),
    StageDefinition(
        stage_id="04_immunogenicity",
        name="Immunogenicity Scoring",
        module_path="pipeline.04_immunogenicity",
        entry_function="score_immunogenicity",
        description="IEDB MHC-I binding predictions",
    ),
    StageDefinition(
        stage_id="05_report",
        name="Report Generation",
        module_path="pipeline.05_report",
        entry_function="generate_report",
        description="CSV + Markdown summary of candidates",
    ),
)

# ---------------------------------------------------------------------------
# Pipeline: Drug Target Discovery (v2, 5 stages)
# ---------------------------------------------------------------------------

DRUG_STAGES = (
    StageDefinition(
        stage_id="01_fetch_enzymes",
        name="Fetch Enzymes",
        module_path="drug_targets.01_fetch_enzymes",
        entry_function="fetch_enzymes",
        description="UniProt download by metabolic pathway",
    ),
    StageDefinition(
        stage_id="02_human_comparison",
        name="Human Comparison",
        module_path="drug_targets.02_human_comparison",
        entry_function="compare_with_human",
        description="BLAST against human proteome",
    ),
    StageDefinition(
        stage_id="03_essentiality",
        name="Essentiality Check",
        module_path="drug_targets.03_essentiality",
        entry_function="check_essentiality",
        description="DEG + TriTrypDB knockout data",
    ),
    StageDefinition(
        stage_id="04_druggability",
        name="Druggability Scoring",
        module_path="drug_targets.04_druggability",
        entry_function="score_druggability",
        description="Composite druggability score + AlphaFold links",
    ),
    StageDefinition(
        stage_id="05_report",
        name="Report Generation",
        module_path="drug_targets.05_report",
        entry_function="generate_drug_targets_report",
        description="CSV + Markdown summary of drug targets",
    ),
)

# ---------------------------------------------------------------------------
# Pipeline: Molecular Docking (v3, 5 stages)
# ---------------------------------------------------------------------------

DOCKING_STAGES = (
    StageDefinition(
        stage_id="06_fetch_structures",
        name="Fetch Structures",
        module_path="drug_targets.06_fetch_structures",
        entry_function="fetch_and_prepare_structures",
        description="AlphaFold PDB + PDBQT preparation",
    ),
    StageDefinition(
        stage_id="07_compound_library",
        name="Compound Library",
        module_path="drug_targets.07_compound_library",
        entry_function="build_compound_library",
        description="ChEMBL + repurposing library assembly",
    ),
    StageDefinition(
        stage_id="08_docking",
        name="Molecular Docking",
        module_path="drug_targets.08_docking",
        entry_function="run_docking_campaign",
        description="AutoDock Vina docking simulations",
    ),
    StageDefinition(
        stage_id="09_admet_filter",
        name="ADMET Filter",
        module_path="drug_targets.09_admet_filter",
        entry_function="filter_admet",
        description="Lipinski + pkCSM ADMET predictions",
    ),
    StageDefinition(
        stage_id="10_docking_report",
        name="Docking Report",
        module_path="drug_targets.10_docking_report",
        entry_function="generate_docking_report",
        description="Top hits + PyMOL scripts",
    ),
)

# ---------------------------------------------------------------------------
# Pipeline: RNA Entropy (7 stages)
# ---------------------------------------------------------------------------

RNA_STAGES = (
    StageDefinition(
        stage_id="01_fetch_transcriptome",
        name="Fetch Transcriptome",
        module_path="rna_entropy.01_fetch_transcriptome",
        entry_function="fetch_transcriptome",
        description="UniProt mRNA download",
    ),
    StageDefinition(
        stage_id="02_codon_usage",
        name="Codon Usage Analysis",
        module_path="rna_entropy.02_codon_usage",
        entry_function="analyze_codon_usage",
        description="RSCU computation for codon bias",
    ),
    StageDefinition(
        stage_id="03_shannon_entropy",
        name="Shannon Entropy",
        module_path="rna_entropy.03_shannon_entropy",
        entry_function="calculate_shannon_entropy",
        description="Per-transcript information content",
    ),
    StageDefinition(
        stage_id="04_sl_rna_analysis",
        name="SL RNA Analysis",
        module_path="rna_entropy.04_sl_rna_analysis",
        entry_function="analyze_sl_rna",
        description="Spliced leader detection and analysis",
    ),
    StageDefinition(
        stage_id="05_human_comparison",
        name="Human Comparison",
        module_path="rna_entropy.05_human_comparison",
        entry_function="compare_with_human",
        description="Entropy delta calculation vs human",
    ),
    StageDefinition(
        stage_id="06_structure_prediction",
        name="Structure Prediction",
        module_path="rna_entropy.06_structure_prediction",
        entry_function="predict_rna_structures",
        description="MFE via ViennaRNA or GC fallback",
    ),
    StageDefinition(
        stage_id="07_report",
        name="Report Generation",
        module_path="rna_entropy.07_report",
        entry_function="generate_rna_report",
        description="CSV + Markdown summary of RNA targets",
    ),
)

# ---------------------------------------------------------------------------
# Pipeline: ASO Math Validation (5 modules + certificate)
# ---------------------------------------------------------------------------

ASO_MATH_STAGES = (
    StageDefinition(
        stage_id="01_thermodynamic_landscape",
        name="Thermodynamic Landscape",
        module_path="aso_math.01_thermodynamic_landscape.run",
        entry_function="main",
        description="Tm, dG, dH, dS calculations for ASO-target duplex",
    ),
    StageDefinition(
        stage_id="02_selectivity_proof",
        name="Selectivity Proof",
        module_path="aso_math.02_selectivity_proof.run",
        entry_function="main",
        description="Off-target analysis against host transcriptome",
    ),
    StageDefinition(
        stage_id="03_evolutionary_conservation",
        name="Evolutionary Conservation",
        module_path="aso_math.03_evolutionary_conservation.run",
        entry_function="main",
        description="Cross-species SL RNA conservation analysis",
    ),
    StageDefinition(
        stage_id="04_exhaustive_optimization",
        name="Exhaustive Optimization",
        module_path="aso_math.04_exhaustive_optimization.run",
        entry_function="main",
        description="Full ASO length/position scan for optimal design",
    ),
    StageDefinition(
        stage_id="05_resistance_model",
        name="Resistance Model",
        module_path="aso_math.05_resistance_model.run",
        entry_function="main",
        description="Evolutionary resistance escape probability",
    ),
)

# ---------------------------------------------------------------------------
# Global registry
# ---------------------------------------------------------------------------

PIPELINES: dict[str, PipelineDefinition] = {
    "vaccine": PipelineDefinition(
        pipeline_id="vaccine",
        display_name="Reverse Vaccinology",
        description="Vaccine candidate discovery: genome fetch, surface filter, conservation, immunogenicity, report",
        stages=VACCINE_STAGES,
        config_module="core.pipeline_config",
        config_class="VaccineConfig",
        default_preset="leishmania_infantum",
    ),
    "drug": PipelineDefinition(
        pipeline_id="drug",
        display_name="Drug Target Discovery",
        description="Enzymatic drug target identification: enzyme fetch, human comparison, essentiality, druggability, report",
        stages=DRUG_STAGES,
        config_module="core.pipeline_config",
        config_class="DrugConfig",
        default_preset="leishmania_infantum",
    ),
    "docking": PipelineDefinition(
        pipeline_id="docking",
        display_name="Molecular Docking",
        description="Structure-based drug design: fetch structures, compound library, docking, ADMET, report",
        stages=DOCKING_STAGES,
        config_module="core.pipeline_config",
        config_class="DockingConfig",
        default_preset="leishmania_infantum",
    ),
    "rna": PipelineDefinition(
        pipeline_id="rna",
        display_name="RNA Entropy Analysis",
        description="Information theory RNA target discovery: transcriptome, codon usage, entropy, SL RNA, human comparison, structure, report",
        stages=RNA_STAGES,
        config_module="core.pipeline_config",
        config_class="RnaConfig",
        default_preset="leishmania_infantum",
    ),
    "aso_math": PipelineDefinition(
        pipeline_id="aso_math",
        display_name="ASO Mathematical Validation",
        description="5-module mathematical proof suite for antisense oligonucleotide candidates",
        stages=ASO_MATH_STAGES,
        config_module="aso_math.target_config",
        config_class="TargetConfig",
        default_preset="leishmania_infantum",
    ),
}


def get_pipeline(pipeline_id: str) -> PipelineDefinition:
    """Retorna a definicao de um pipeline pelo ID.

    Raises:
        ValueError: Se o pipeline nao existir no registro.
    """
    if pipeline_id not in PIPELINES:
        available = ", ".join(sorted(PIPELINES.keys()))
        raise ValueError(
            f"Pipeline desconhecido: '{pipeline_id}'. Disponiveis: {available}"
        )
    return PIPELINES[pipeline_id]


def list_pipelines() -> list[dict[str, str]]:
    """Retorna lista resumida de pipelines disponiveis."""
    return [
        {
            "id": p.pipeline_id,
            "name": p.display_name,
            "stages": str(len(p.stages)),
            "description": p.description,
        }
        for p in PIPELINES.values()
    ]
