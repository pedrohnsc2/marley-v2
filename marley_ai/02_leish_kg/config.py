"""Configuracao do modulo 02_leish_kg — Knowledge Graph."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Caminhos especificos do modulo
# ---------------------------------------------------------------------------

KG_ROOT: Final[Path] = AI_ROOT / "02_leish_kg"
KG_RESULTS_DIR: Final[Path] = KG_ROOT / "results"


# ---------------------------------------------------------------------------
# Tipos de nos e arestas do grafo (ontologia)
# ---------------------------------------------------------------------------

# Tipos de nos aceitos no KG
NODE_TYPES: Final[tuple[str, ...]] = (
    "Gene",
    "Protein",
    "Epitope",
    "Drug",
    "ASO",
    "RNA",
    "Pathway",
    "Paper",
    "Organism",
)

# Tipos de arestas aceitos no KG
EDGE_TYPES: Final[tuple[str, ...]] = (
    "ENCODES",         # Gene -> Protein
    "INHIBITS",        # ASO/Drug -> Protein/RNA
    "PART_OF",         # Protein -> Pathway, Gene -> Organism
    "VALIDATED_BY",    # Protein/Epitope -> Paper
    "BINDS",           # Epitope -> Protein, ASO -> RNA
    "ACTIVATES",       # Protein -> Pathway
    "TARGETS",         # ASO -> RNA, Drug -> Protein
    "DERIVED_FROM",    # Epitope -> Gene
    "EXPRESSED_IN",    # Gene -> Organism
    "INFECTS",         # Organism -> Organism
)


@dataclass(frozen=True)
class KGConfig(AIModuleConfig):
    """Configuracao especifica para o Knowledge Graph de Leishmania.

    Campos adicionais controlam formato do grafo, limites de tamanho,
    e caminhos para resultados.
    """

    # --- Formato do grafo ---
    graph_format: str = "networkx"  # "networkx" (puro Python, sem Neo4j)
    max_nodes: int = 50_000         # Limite de nos (proteoma completo ~8k)
    max_edges: int = 200_000        # Limite de arestas

    # --- Caminhos ---
    results_dir: Path = field(
        default_factory=lambda: KG_RESULTS_DIR
    )

    # --- Exportacao ---
    export_json: bool = True        # Exportar grafo como JSON
    export_graphml: bool = False    # Exportar como GraphML (requer networkx)
