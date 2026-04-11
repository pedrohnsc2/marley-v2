"""Configuracao do modulo 11_scientist — AI Scientist.

Define caminhos para resultados de todos os modulos, papeis dos agentes,
e parametros do loop de descoberta. Sem dependencia de LLM — opera
como motor de sintese computacional deterministico.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from marley_ai.config import AIModuleConfig, AI_ROOT, RESULTS_DIR


# ---------------------------------------------------------------------------
# Mapeamento de modulos -> caminhos de resultados
# Cada modulo grava seu JSON em local diferente; aqui centralizamos.
# ---------------------------------------------------------------------------

MODULE_RESULT_PATHS: Final[dict[str, Path]] = {
    "01_rag": RESULTS_DIR / "01_rag.json",
    "02_leish_kg": RESULTS_DIR / "02_leish_kg.json",
    "03_leish_esm": AI_ROOT / "03_leish_esm" / "results" / "03_leish_esm.json",
    "04_rna_fm": RESULTS_DIR / "04_rna_fm.json",
    "05_rosettafold": AI_ROOT / "05_rosettafold" / "results" / "05_rosettafold.json",
    "06_evodiff": AI_ROOT / "06_evodiff" / "results" / "06_evodiff.json",
    "07_contrastive": AI_ROOT / "07_contrastive" / "results" / "07_contrastive.json",
    "08_rl_ppo": AI_ROOT / "08_rl_ppo" / "results" / "08_rl_ppo.json",
    "09_sae": RESULTS_DIR / "09_sae.json",
    "10_digital_twin": RESULTS_DIR / "10_digital_twin.json",
}

# ---------------------------------------------------------------------------
# Definicao dos papeis dos agentes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AgentRole:
    """Define papel, fontes de dados e capacidades de um agente."""

    name: str
    description: str
    # Slugs dos modulos cujos resultados este agente consome
    source_modules: tuple[str, ...]
    # Tipo de analise que produz
    output_type: str


AGENT_ROLES: Final[dict[str, AgentRole]] = {
    "literature": AgentRole(
        name="LiteratureAgent",
        description="Consulta resultados do RAG para identificar gaps no conhecimento",
        source_modules=("01_rag", "02_leish_kg"),
        output_type="knowledge_gaps",
    ),
    "design": AgentRole(
        name="DesignAgent",
        description="Propoe variantes de ASO e epitopos usando EvoDiff + RL",
        source_modules=("06_evodiff", "08_rl_ppo", "04_rna_fm"),
        output_type="variant_proposals",
    ),
    "validation": AgentRole(
        name="ValidationAgent",
        description="Avalia candidatos usando termodinamica, estrutura e entrega",
        source_modules=("04_rna_fm", "05_rosettafold", "08_rl_ppo", "06_evodiff"),
        output_type="validation_report",
    ),
    "knowledge": AgentRole(
        name="KnowledgeAgent",
        description="Usa KG e ESM para encontrar conexoes inesperadas",
        source_modules=("02_leish_kg", "03_leish_esm", "09_sae"),
        output_type="unexpected_connections",
    ),
    "report": AgentRole(
        name="ReportAgent",
        description="Sintetiza outputs de todos os agentes em narrativa cientifica",
        source_modules=(
            "01_rag", "02_leish_kg", "03_leish_esm", "04_rna_fm",
            "05_rosettafold", "06_evodiff", "08_rl_ppo", "09_sae",
            "10_digital_twin",
        ),
        output_type="scientific_narrative",
    ),
}


# ---------------------------------------------------------------------------
# Parametros do loop de descoberta
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DiscoveryConfig:
    """Parametros para o loop de descoberta cientifica."""

    # Numero maximo de iteracoes do loop (1 por padrao — expansivel)
    max_iterations: int = 1

    # Limiar minimo de confianca para aceitar uma hipotese (0-1)
    hypothesis_confidence_threshold: float = 0.6

    # Numero maximo de hipoteses a gerar por iteracao
    max_hypotheses_per_iteration: int = 10

    # Limiar de concordancia entre agentes para "consenso" (fracao)
    consensus_threshold: float = 0.6

    # Numero maximo de experimentos propostos por iteracao
    max_proposed_experiments: int = 5

    # Pesos para calculo de score de hipotese
    # (evidencia estrutural, termodinamica, dados de ML, consistencia KG)
    scoring_weights: dict[str, float] = field(default_factory=lambda: {
        "structural_evidence": 0.25,
        "thermodynamic_evidence": 0.25,
        "ml_evidence": 0.25,
        "knowledge_graph_support": 0.25,
    })


# ---------------------------------------------------------------------------
# Configuracao principal do modulo
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ScientistConfig(AIModuleConfig):
    """Configuracao do modulo 11_scientist.

    Herda de AIModuleConfig mas nao depende de LLM — opera como
    motor de sintese computacional que le resultados dos modulos 01-10
    e gera hipoteses, propostas experimentais e relatorio consolidado.
    """

    # Parametros do loop de descoberta
    discovery: DiscoveryConfig = field(default_factory=DiscoveryConfig)

    # Diretorio de saida especifico do modulo
    scientist_output_dir: Path = field(
        default_factory=lambda: AI_ROOT / "11_scientist" / "results"
    )
