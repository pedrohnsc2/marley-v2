"""Execucao do modulo 02_leish_kg — Knowledge Graph de Leishmania.

Constroi um Knowledge Graph in-memory usando networkx a partir dos
dados CSV/JSON gerados pelo pipeline Marley. Integra genes, pathways,
compostos, epitopos e ASO candidates num grafo unificado que permite
raciocinio sobre relacoes parasita-hospedeiro-terapia.

Uso:
    python -m marley_ai.02_leish_kg.run
"""

from __future__ import annotations

import csv
import importlib as _il
import json
from collections import Counter
from pathlib import Path
from typing import Any

import networkx as nx

from marley_ai.config import AIModuleConfig, PROJECT_ROOT
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

# Importacao via importlib — Python nao permite "from marley_ai.02_leish_kg..."
# porque "02" e interpretado como literal numerico invalido no parser.
_kg_config = _il.import_module("marley_ai.02_leish_kg.config")
_kg_schema = _il.import_module("marley_ai.02_leish_kg.schema")

KGConfig = _kg_config.KGConfig
EdgeType = _kg_schema.EdgeType
GeneNode = _kg_schema.GeneNode
PathwayNode = _kg_schema.PathwayNode
OrganismNode = _kg_schema.OrganismNode
EpitopeNode = _kg_schema.EpitopeNode
ASONode = _kg_schema.ASONode
KGNode = _kg_schema.KGNode
KGEdge = _kg_schema.KGEdge

# Compound node is not in the schema — we use the base DrugNode
# (compounds from docking are drug candidates)
DrugNode = _kg_schema.DrugNode


# ---------------------------------------------------------------------------
# Data file paths (relative to PROJECT_ROOT)
# ---------------------------------------------------------------------------

DRUG_TARGETS_CSV: Path = PROJECT_ROOT / "results" / "drug_targets_top20.csv"
CONSTRUCT_JSON: Path = PROJECT_ROOT / "results" / "construct" / "construct_card.json"
DOCKING_CSV: Path = PROJECT_ROOT / "results" / "docking_scores.csv"
ASO_CSV: Path = PROJECT_ROOT / "results" / "aso" / "aso_candidates.csv"
RNA_TARGETS_JSON: Path = PROJECT_ROOT / "results" / "rna" / "rna_targets_scored.json"


# ---------------------------------------------------------------------------
# Helper: add node and edge to graph
# ---------------------------------------------------------------------------

def _add_node(graph: nx.DiGraph, node: KGNode) -> None:
    """Add a typed node to the graph with all its attributes."""
    graph.add_node(
        node.node_id,
        node_type=node.node_type,
        label=node.label,
        **node.properties,
    )


def _add_edge(graph: nx.DiGraph, edge: KGEdge) -> None:
    """Add a typed edge to the graph with all its properties."""
    graph.add_edge(
        edge.source_id,
        edge.target_id,
        edge_type=edge.edge_type.value,
        **edge.properties,
    )


# ---------------------------------------------------------------------------
# Data loaders — each reads one data source and populates the graph
# ---------------------------------------------------------------------------

def _load_organism(graph: nx.DiGraph) -> None:
    """Create the single L. infantum organism node."""
    org = OrganismNode(
        taxon_id="5671",
        scientific_name="Leishmania infantum",
        common_name="L. infantum",
        role="parasita",
    )
    _add_node(graph, org)


def _load_drug_targets(
    graph: nx.DiGraph, warnings: list[str],
) -> tuple[int, int]:
    """Load gene and pathway nodes from drug_targets_top20.csv.

    Returns (n_genes, n_pathways) loaded.
    """
    if not DRUG_TARGETS_CSV.exists():
        warnings.append(f"Arquivo nao encontrado: {DRUG_TARGETS_CSV}")
        return 0, 0

    genes_created: set[str] = set()
    pathways_created: set[str] = set()

    with open(DRUG_TARGETS_CSV, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gene_id = row["gene_id"]
            gene_name = row["gene_name"]
            pathway = row["pathway"]

            # Gene node (deduplicate by gene_id)
            if gene_id not in genes_created:
                gene_node = GeneNode(
                    gene_id=gene_id,
                    gene_name=gene_name,
                    organism="L. infantum",
                )
                _add_node(graph, gene_node)
                genes_created.add(gene_id)

                # Gene -> Organism: part_of
                _add_edge(graph, KGEdge(
                    source_id=f"gene:{gene_id}",
                    target_id="organism:5671",
                    edge_type=EdgeType.PART_OF,
                ))

            # Pathway node (deduplicate)
            if pathway and pathway not in pathways_created:
                pw_node = PathwayNode(
                    pathway_id=pathway,
                    pathway_name=pathway.replace("_", " ").title(),
                )
                _add_node(graph, pw_node)
                pathways_created.add(pathway)

            # Gene -> Pathway: belongs_to_pathway (PART_OF)
            if pathway:
                gene_node_id = f"gene:{gene_id}"
                pathway_node_id = f"pathway:{pathway}"
                if not graph.has_edge(gene_node_id, pathway_node_id):
                    _add_edge(graph, KGEdge(
                        source_id=gene_node_id,
                        target_id=pathway_node_id,
                        edge_type=EdgeType.PART_OF,
                        properties={"relation": "belongs_to_pathway"},
                    ))

    return len(genes_created), len(pathways_created)


def _load_docking_compounds(
    graph: nx.DiGraph, warnings: list[str],
) -> tuple[int, int]:
    """Load compound nodes and docking edges from docking_scores.csv.

    Creates Compound nodes (using DrugNode) and edges to Gene nodes
    with binding_affinity as weight.

    Returns (n_compounds, n_docking_edges) loaded.
    """
    if not DOCKING_CSV.exists():
        warnings.append(f"Arquivo nao encontrado: {DOCKING_CSV}")
        return 0, 0

    compounds_created: set[str] = set()
    n_edges = 0

    with open(DOCKING_CSV, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            compound_id = row["compound_id"]
            compound_name = row["compound_name"]
            target_gene_id = row["target_gene_id"]
            target_gene_name = row["target_gene_name"]
            binding_affinity = float(row["binding_affinity"])

            # Compound node (deduplicate)
            if compound_id not in compounds_created:
                compound_node = DrugNode(
                    drug_id=compound_id,
                    drug_name=compound_name,
                    mechanism="docking_candidate",
                    status="virtual_screening",
                )
                # Override node_type to "Compound" in the graph
                graph.add_node(
                    compound_node.node_id,
                    node_type="Compound",
                    label=compound_name,
                    **compound_node.properties,
                )
                compounds_created.add(compound_id)

            # Ensure target gene exists in graph; create if absent
            gene_node_id = f"gene:{target_gene_id}"
            if not graph.has_node(gene_node_id):
                gene_node = GeneNode(
                    gene_id=target_gene_id,
                    gene_name=target_gene_name,
                    organism="L. infantum",
                )
                _add_node(graph, gene_node)
                _add_edge(graph, KGEdge(
                    source_id=gene_node_id,
                    target_id="organism:5671",
                    edge_type=EdgeType.PART_OF,
                ))

            # Compound -> Gene: docks_with (TARGETS with binding affinity)
            _add_edge(graph, KGEdge(
                source_id=compound_node.node_id,
                target_id=gene_node_id,
                edge_type=EdgeType.TARGETS,
                properties={
                    "relation": "docks_with",
                    "binding_affinity": binding_affinity,
                },
            ))
            n_edges += 1

    return len(compounds_created), n_edges


def _extract_gene_id_from_construct(raw_gene_id: str) -> str:
    """Extract the core LINF gene identifier from construct_card gene_id.

    Examples:
        "LINF_240013900-T1-p1" -> "LINF_240013900"
        "LINF_080015500-T1-p1" -> "LINF_080015500"
    """
    if raw_gene_id.startswith("LINF_"):
        parts = raw_gene_id.split("-")
        return parts[0]
    return raw_gene_id


def _extract_gene_name_from_construct(raw_gene_name: str) -> str:
    """Extract a readable gene product name from construct_card gene_name.

    The construct_card gene_name field is a long annotation string.
    We extract the gene_product field from it.
    """
    if "gene_product=" in raw_gene_name:
        start = raw_gene_name.index("gene_product=") + len("gene_product=")
        rest = raw_gene_name[start:]
        end = rest.find("|")
        if end != -1:
            return rest[:end].strip()
        return rest.strip()
    return raw_gene_name[:50] if len(raw_gene_name) > 50 else raw_gene_name


def _load_epitopes(
    graph: nx.DiGraph, warnings: list[str],
) -> int:
    """Load epitope nodes from construct_card.json.

    Creates Epitope -> Gene edges (derived_from).

    Returns n_epitopes loaded.
    """
    if not CONSTRUCT_JSON.exists():
        warnings.append(f"Arquivo nao encontrado: {CONSTRUCT_JSON}")
        return 0

    with open(CONSTRUCT_JSON, encoding="utf-8") as fh:
        data = json.load(fh)

    epitopes = data.get("epitopes", [])
    if not epitopes:
        warnings.append("Nenhum epitopo encontrado em construct_card.json")
        return 0

    epitopes_created: set[str] = set()
    genes_created: set[str] = set()

    for ep in epitopes:
        peptide = ep["peptide"]

        # Deduplicate epitopes by peptide sequence (e.g., LLTANVCYK appears
        # multiple times for different gene copies)
        if peptide in epitopes_created:
            continue

        allele = ep.get("allele", "")
        ic50 = float(ep.get("ic50", 0.0))
        raw_gene_id = ep.get("gene_id", "")
        raw_gene_name = ep.get("gene_name", "")

        gene_id = _extract_gene_id_from_construct(raw_gene_id)
        gene_name = _extract_gene_name_from_construct(raw_gene_name)

        # Epitope node
        epitope_node = EpitopeNode(
            peptide=peptide,
            allele=allele,
            ic50=ic50,
            gene_id=gene_id,
        )
        _add_node(graph, epitope_node)
        epitopes_created.add(peptide)

        # Ensure gene exists
        gene_node_id = f"gene:{gene_id}"
        if gene_id not in genes_created and not graph.has_node(gene_node_id):
            gene_node = GeneNode(
                gene_id=gene_id,
                gene_name=gene_name,
                organism="L. infantum",
            )
            _add_node(graph, gene_node)
            genes_created.add(gene_id)

            # Gene -> Organism
            _add_edge(graph, KGEdge(
                source_id=gene_node_id,
                target_id="organism:5671",
                edge_type=EdgeType.PART_OF,
            ))

        # Epitope -> Gene: derived_from
        _add_edge(graph, KGEdge(
            source_id=epitope_node.node_id,
            target_id=gene_node_id,
            edge_type=EdgeType.DERIVED_FROM,
            properties={"relation": "derived_from"},
        ))

    return len(epitopes_created)


def _load_aso_candidates(
    graph: nx.DiGraph, warnings: list[str],
) -> int:
    """Load top ASO candidate as a single node.

    We load only the top-ranked MRL-ASO-001, since all ASOs target the
    same SL RNA region. The ASO is connected to a virtual SL RNA gene node.

    Returns 1 if loaded, 0 otherwise.
    """
    if not ASO_CSV.exists():
        warnings.append(f"Arquivo nao encontrado: {ASO_CSV}")
        return 0

    with open(ASO_CSV, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        top_aso = next(reader, None)

    if top_aso is None:
        warnings.append("Nenhum ASO candidate encontrado em aso_candidates.csv")
        return 0

    aso_id = top_aso["aso_id"]
    aso_sequence = top_aso["aso_sequence"]
    tm = float(top_aso.get("tm_celsius", 0.0))
    dg = float(top_aso.get("delta_g_kcal", 0.0))
    gc_content = float(top_aso.get("gc_content", 0.0))
    target_start = int(top_aso.get("target_start", 0))
    target_end = int(top_aso.get("target_end", 0))

    # ASO node
    aso_node = ASONode(
        aso_id=aso_id,
        sequence=aso_sequence,
        target_start=target_start,
        target_end=target_end,
        tm=tm,
        dg=dg,
        gc_content=gc_content,
    )
    _add_node(graph, aso_node)

    # Create a virtual gene node for SL RNA (the ASO target)
    sl_gene_id = "SL_RNA"
    sl_gene_node_id = f"gene:{sl_gene_id}"
    if not graph.has_node(sl_gene_node_id):
        sl_gene = GeneNode(
            gene_id=sl_gene_id,
            gene_name="Spliced Leader RNA",
            organism="L. infantum",
        )
        _add_node(graph, sl_gene)

        # SL RNA Gene -> Organism
        _add_edge(graph, KGEdge(
            source_id=sl_gene_node_id,
            target_id="organism:5671",
            edge_type=EdgeType.PART_OF,
        ))

    # ASO -> Gene (SL RNA): targets
    _add_edge(graph, KGEdge(
        source_id=aso_node.node_id,
        target_id=sl_gene_node_id,
        edge_type=EdgeType.TARGETS,
        properties={
            "relation": "targets",
            "mechanism": "Antisense hybridization + RNase H cleavage",
            "tm_celsius": tm,
            "delta_g_kcal": dg,
        },
    ))

    return 1


# ---------------------------------------------------------------------------
# Graph analysis functions
# ---------------------------------------------------------------------------

def _get_node_type_counts(graph: nx.DiGraph) -> dict[str, int]:
    """Count nodes by type."""
    types: list[str] = [
        data.get("node_type", "Unknown")
        for _, data in graph.nodes(data=True)
    ]
    return dict(Counter(types).most_common())


def _get_edge_type_counts(graph: nx.DiGraph) -> dict[str, int]:
    """Count edges by type."""
    types: list[str] = [
        data.get("edge_type", "Unknown")
        for _, _, data in graph.edges(data=True)
    ]
    return dict(Counter(types).most_common())


def _get_top_hubs(graph: nx.DiGraph, n: int = 10) -> list[dict[str, Any]]:
    """Return top N nodes by degree (most connected)."""
    degree_list = sorted(
        graph.degree(), key=lambda x: x[1], reverse=True,
    )
    result = []
    for node_id, degree in degree_list[:n]:
        node_data = graph.nodes[node_id]
        result.append({
            "node_id": node_id,
            "label": node_data.get("label", node_id),
            "node_type": node_data.get("node_type", "Unknown"),
            "degree": degree,
        })
    return result


def _get_betweenness_top(
    graph: nx.DiGraph, n: int = 5,
) -> list[dict[str, Any]]:
    """Return top N nodes by betweenness centrality.

    Uses the undirected view for betweenness since the graph is directed
    but we want to measure bridging importance regardless of direction.
    """
    undirected = graph.to_undirected()
    centrality = nx.betweenness_centrality(undirected)
    sorted_centrality = sorted(
        centrality.items(), key=lambda x: x[1], reverse=True,
    )
    result = []
    for node_id, score in sorted_centrality[:n]:
        node_data = graph.nodes[node_id]
        result.append({
            "node_id": node_id,
            "label": node_data.get("label", node_id),
            "node_type": node_data.get("node_type", "Unknown"),
            "betweenness_centrality": round(score, 6),
        })
    return result


def _get_pathway_distribution(graph: nx.DiGraph) -> dict[str, int]:
    """Count how many genes belong to each pathway via PART_OF edges."""
    distribution: dict[str, int] = {}

    for source, target, data in graph.edges(data=True):
        edge_type = data.get("edge_type", "")
        relation = data.get("relation", "")

        # Only count Gene -> Pathway edges
        if edge_type == "PART_OF" and relation == "belongs_to_pathway":
            source_data = graph.nodes.get(source, {})
            target_data = graph.nodes.get(target, {})
            if (source_data.get("node_type") == "Gene"
                    and target_data.get("node_type") == "Pathway"):
                pathway_label = target_data.get("label", target)
                distribution[pathway_label] = distribution.get(pathway_label, 0) + 1

    # Sort by count descending
    return dict(sorted(
        distribution.items(), key=lambda x: x[1], reverse=True,
    ))


def _get_degree_distribution(graph: nx.DiGraph) -> dict[str, Any]:
    """Compute degree distribution summary statistics."""
    degrees = [d for _, d in graph.degree()]
    if not degrees:
        return {"min": 0, "max": 0, "mean": 0.0, "median": 0.0}

    degrees_sorted = sorted(degrees)
    n = len(degrees_sorted)
    mean_deg = sum(degrees) / n
    median_deg = (
        degrees_sorted[n // 2]
        if n % 2 == 1
        else (degrees_sorted[n // 2 - 1] + degrees_sorted[n // 2]) / 2
    )

    return {
        "min": degrees_sorted[0],
        "max": degrees_sorted[-1],
        "mean": round(mean_deg, 2),
        "median": round(float(median_deg), 2),
    }


# ---------------------------------------------------------------------------
# Module registration
# ---------------------------------------------------------------------------

@register("02_leish_kg")
class LeishKGModule:
    """Modulo Knowledge Graph de Leishmania para o sistema Marley AI."""

    def __init__(self) -> None:
        self._config: KGConfig | None = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao do orquestrador."""
        if isinstance(config, KGConfig):
            self._config = config
        else:
            self._config = KGConfig(
                module_slug="02_leish_kg",
                module_name="Leishmania Knowledge Graph",
            )

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica se os arquivos de dados estao disponiveis."""
        missing = []
        for path in [DRUG_TARGETS_CSV, CONSTRUCT_JSON, DOCKING_CSV, ASO_CSV]:
            if not path.exists():
                missing.append(str(path))
        return {"valid": len(missing) == 0, "missing": missing}

    def run(self) -> dict[str, Any]:
        """Executa pipeline de construcao do Knowledge Graph."""
        return main(self._config)

    def get_dependencies(self) -> list[str]:
        """Sem dependencias de outros modulos AI — usa dados do pipeline."""
        return []


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Constroi o Knowledge Graph de Leishmania a partir dos dados do pipeline.

    Etapas:
        1. Cria o no Organism (L. infantum)
        2. Carrega genes e pathways de drug_targets_top20.csv
        3. Carrega compostos e docking edges de docking_scores.csv
        4. Carrega epitopos de construct_card.json
        5. Carrega ASO candidate de aso_candidates.csv
        6. Executa analise do grafo (hubs, centralidade, distribuicao)
        7. Monta envelope de resultados

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com estatisticas e analise do grafo.
    """
    if config is None:
        config = KGConfig(
            module_slug="02_leish_kg",
            module_name="Leishmania Knowledge Graph",
        )

    envelope = create_envelope("02_leish_kg")
    envelope["device"] = "cpu"

    with Timer() as timer:
        graph = nx.DiGraph()
        warnings: list[str] = []

        # -------------------------------------------------------------------
        # 1. Organism node — base of the graph
        # -------------------------------------------------------------------
        print("[02_leish_kg] Criando no base do organismo...")
        _load_organism(graph)

        # -------------------------------------------------------------------
        # 2. Genes and pathways from drug targets
        # -------------------------------------------------------------------
        print("[02_leish_kg] Carregando drug targets...")
        n_genes, n_pathways = _load_drug_targets(graph, warnings)
        print(f"  Genes: {n_genes}, Pathways: {n_pathways}")

        # -------------------------------------------------------------------
        # 3. Compounds from docking scores
        # -------------------------------------------------------------------
        print("[02_leish_kg] Carregando docking compounds...")
        n_compounds, n_docking_edges = _load_docking_compounds(graph, warnings)
        print(f"  Compounds: {n_compounds}, Docking edges: {n_docking_edges}")

        # -------------------------------------------------------------------
        # 4. Epitopes from vaccine construct
        # -------------------------------------------------------------------
        print("[02_leish_kg] Carregando epitopos...")
        n_epitopes = _load_epitopes(graph, warnings)
        print(f"  Epitopos: {n_epitopes}")

        # -------------------------------------------------------------------
        # 5. ASO candidate
        # -------------------------------------------------------------------
        print("[02_leish_kg] Carregando ASO candidate...")
        n_aso = _load_aso_candidates(graph, warnings)
        print(f"  ASO nodes: {n_aso}")

        # -------------------------------------------------------------------
        # 6. Graph analysis
        # -------------------------------------------------------------------
        print("\n[02_leish_kg] Analisando grafo...")

        n_nodes = graph.number_of_nodes()
        n_edges = graph.number_of_edges()
        density = nx.density(graph)
        n_connected = nx.number_connected_components(graph.to_undirected())

        node_type_counts = _get_node_type_counts(graph)
        edge_type_counts = _get_edge_type_counts(graph)
        top_hubs = _get_top_hubs(graph, n=10)
        betweenness_top5 = _get_betweenness_top(graph, n=5)
        pathway_dist = _get_pathway_distribution(graph)
        degree_dist = _get_degree_distribution(graph)

        n_node_types = len(node_type_counts)
        n_edge_types = len(edge_type_counts)

        top_hub_label = top_hubs[0]["label"] if top_hubs else "N/A"

        print(f"  Nodes: {n_nodes}")
        print(f"  Edges: {n_edges}")
        print(f"  Density: {density:.4f}")
        print(f"  Connected components: {n_connected}")
        print(f"  Node types: {n_node_types}")
        print(f"  Edge types: {n_edge_types}")
        print(f"  Top hub: {top_hub_label}")

        if top_hubs:
            print("\n  Top 5 hubs by degree:")
            for hub in top_hubs[:5]:
                print(
                    f"    {hub['label']} ({hub['node_type']}): "
                    f"degree={hub['degree']}"
                )

        if betweenness_top5:
            print("\n  Top 5 by betweenness centrality:")
            for node in betweenness_top5:
                print(
                    f"    {node['label']} ({node['node_type']}): "
                    f"centrality={node['betweenness_centrality']:.4f}"
                )

        # -------------------------------------------------------------------
        # 7. Assemble result envelope
        # -------------------------------------------------------------------
        envelope["status"] = "complete"
        envelope["warnings"] = warnings
        envelope["dependencies"] = [
            "results/drug_targets_top20.csv",
            "results/construct/construct_card.json",
            "results/docking_scores.csv",
            "results/aso/aso_candidates.csv",
        ]

        envelope["summary"]["conclusion"] = (
            f"Knowledge Graph construido com {n_nodes} nodes e {n_edges} edges. "
            f"{n_node_types} tipos de nodes, {n_edge_types} tipos de relacoes. "
            f"Hub gene: {top_hub_label}."
        )

        envelope["summary"]["key_metrics"] = {
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "n_relation_types": n_edge_types,
            "n_node_types": n_node_types,
            "n_genes": node_type_counts.get("Gene", 0),
            "n_pathways": node_type_counts.get("Pathway", 0),
            "n_compounds": node_type_counts.get("Compound", 0),
            "n_epitopes": node_type_counts.get("Epitope", 0),
            "density": round(density, 6),
            "n_connected_components": n_connected,
        }

        envelope["data"] = {
            "node_type_counts": node_type_counts,
            "edge_type_counts": edge_type_counts,
            "top_hubs": top_hubs,
            "betweenness_top5": betweenness_top5,
            "pathway_distribution": pathway_dist,
            "degree_distribution": degree_dist,
        }

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"\n[02_leish_kg] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
