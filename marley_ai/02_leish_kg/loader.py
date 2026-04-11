"""Construcao do Knowledge Graph a partir dos dados do projeto Marley.

Carrega entidades biologicas de tres fontes internas:
    1. vaccine_platforms/shared/epitopes.py — 11 epitopos + genes de origem
    2. drug_targets/04_druggability.py — 7 alvos farmacologicos validados
    3. aso_math/config.py — ASO MRL-ASO-001 + SL RNA alvo

Conecta as entidades com arestas tipadas formando um grafo integrado
que permite raciocinio sobre as relacoes parasita-hospedeiro-terapia.
"""

from __future__ import annotations

from typing import Any

import networkx as nx

from marley_ai.02_leish_kg.schema import (
    ASONode,
    DrugNode,
    EdgeType,
    EpitopeNode,
    GeneNode,
    KGEdge,
    KGNode,
    OrganismNode,
    PathwayNode,
    ProteinNode,
    RNANode,
)


# ---------------------------------------------------------------------------
# Dados dos 7 alvos farmacologicos validados
# (extraidos de drug_targets/04_druggability.py para evitar dependencia
#  de core.db e core.models que requerem Supabase)
# ---------------------------------------------------------------------------

_DRUG_TARGETS: list[dict[str, Any]] = [
    {
        "gene_id": "PRIORITY_TryS",
        "gene_name": "TryS",
        "enzyme_class": "trypanothione_synthetase",
        "pathway": "trypanothione",
        "ec_number": "6.3.1.9",
        "druggability_score": 0.98,
        "description": "Completamente ausente em humanos. Alvo prioritario global.",
    },
    {
        "gene_id": "PRIORITY_TryR",
        "gene_name": "TryR",
        "enzyme_class": "trypanothione_reductase",
        "pathway": "trypanothione",
        "ec_number": "1.8.1.12",
        "druggability_score": 0.96,
        "description": "Ausente em humanos. Inibidores ja identificados in vitro.",
    },
    {
        "gene_id": "PRIORITY_ADL",
        "gene_name": "ADL",
        "enzyme_class": "adenylosuccinate_lyase",
        "pathway": "purine_salvage",
        "ec_number": "4.3.2.2",
        "druggability_score": 0.95,
        "description": "Unica enzima de purinas essencial sozinha.",
    },
    {
        "gene_id": "PRIORITY_SMT",
        "gene_name": "SMT",
        "enzyme_class": "sterol_24C_methyltransferase",
        "pathway": "sterol_biosynthesis",
        "ec_number": "2.1.1.41",
        "druggability_score": 0.94,
        "description": "Ausente em mamiferos. Base dos azolicos antileishmania.",
    },
    {
        "gene_id": "PRIORITY_GMPS",
        "gene_name": "GMPS",
        "enzyme_class": "GMP_synthetase",
        "pathway": "purine_salvage",
        "ec_number": "6.3.5.2",
        "druggability_score": 0.92,
        "description": "Calcanhar de Aquiles da via de purinas.",
    },
    {
        "gene_id": "PRIORITY_6PGDH",
        "gene_name": "6PGDH",
        "enzyme_class": "6_phosphogluconate_dehydrogenase",
        "pathway": "pentose_phosphate",
        "ec_number": "1.1.1.44",
        "druggability_score": 0.90,
        "description": "Menos de 35% identidade com humano. Modelo 3D disponivel.",
    },
    {
        "gene_id": "PRIORITY_XPRT",
        "gene_name": "XPRT",
        "enzyme_class": "xanthine_phosphoribosyltransferase",
        "pathway": "purine_salvage",
        "ec_number": "2.4.2.22",
        "druggability_score": 0.88,
        "description": "Ausente em humanos. Incorpora analogos de purinas.",
    },
]

# Descricoes das vias metabolicas do parasita
_PATHWAY_DESCRIPTIONS: dict[str, str] = {
    "trypanothione": (
        "Via de tripanotiona — sistema redox exclusivo de tripanosomatideos, "
        "substitui glutationa/tioredoxina dos mamiferos. Essencial para "
        "defesa contra estresse oxidativo do macrofago."
    ),
    "purine_salvage": (
        "Via de salvamento de purinas — Leishmania nao sintetiza purinas "
        "de novo, dependendo inteiramente desta via para obter nucleotideos "
        "do hospedeiro. Alvo classico para quimioterapia."
    ),
    "sterol_biosynthesis": (
        "Biossintese de esterois — Leishmania produz ergosterol ao inves "
        "de colesterol, usando enzimas ausentes em mamiferos. Base para "
        "tratamento com azolicos (fluconazol, itraconazol)."
    ),
    "pentose_phosphate": (
        "Via das pentoses-fosfato — produz NADPH e ribose-5-fosfato. "
        "A 6PGDH de Leishmania tem divergencia estrutural suficiente "
        "do homologo humano para ser alvo seletivo."
    ),
    "sl_rna_trans_splicing": (
        "Trans-splicing de SL RNA — mecanismo unico de processamento "
        "de mRNA em tripanosomatideos. Ausente em mamiferos, alvo "
        "do ASO MRL-ASO-001 desenvolvido pelo pipeline Marley."
    ),
}

# Farmacos conhecidos contra leishmaniose (referencia para o grafo)
_KNOWN_DRUGS: list[dict[str, str]] = [
    {
        "drug_id": "miltefosine",
        "drug_name": "Miltefosine",
        "mechanism": "Disruptura de membrana lipidica e apoptose do parasita",
        "status": "approved",
        "targets_pathway": "sterol_biosynthesis",
    },
    {
        "drug_id": "amphotericin_b",
        "drug_name": "Amphotericin B",
        "mechanism": "Liga-se a ergosterol, formando poros na membrana",
        "status": "approved",
        "targets_pathway": "sterol_biosynthesis",
    },
    {
        "drug_id": "antimonial_pentavalente",
        "drug_name": "Antimoniato de meglumina",
        "mechanism": "Inibicao de tripanotiona redutase e topoisomerase",
        "status": "approved",
        "targets_pathway": "trypanothione",
    },
    {
        "drug_id": "allopurinol",
        "drug_name": "Allopurinol",
        "mechanism": "Analogo de purina incorporado via XPRT, inibe sintese de RNA",
        "status": "approved",
        "targets_pathway": "purine_salvage",
    },
]

# Publicacoes-chave do projeto (referencia para validacao)
_KEY_PAPERS: list[dict[str, Any]] = [
    {
        "paper_id": "crooke2017",
        "title": "Antisense technology: an overview and prospectus",
        "doi": "10.1038/nrd.2016.199",
        "year": 2017,
        "authors": "Crooke ST et al.",
    },
    {
        "paper_id": "rogers2011",
        "title": "Chromosome and gene copy number variation in L. infantum",
        "doi": "10.1371/journal.pgen.1002237",
        "year": 2011,
        "authors": "Rogers MB et al.",
    },
    {
        "paper_id": "liang2003",
        "title": "Trans and spliced leader RNA genes in L. infantum",
        "doi": "10.1016/j.ijpara.2003.09.003",
        "year": 2003,
        "authors": "Liang XH et al.",
    },
    {
        "paper_id": "santalucia1998",
        "title": "Unified view of polymer, dumbbell and oligonucleotide DNA NN thermodynamics",
        "doi": "10.1073/pnas.95.4.1460",
        "year": 1998,
        "authors": "SantaLucia J Jr.",
    },
]


# ---------------------------------------------------------------------------
# Funcoes de construcao do grafo
# ---------------------------------------------------------------------------

def _add_node(graph: nx.DiGraph, node: KGNode) -> None:
    """Adiciona um no ao grafo com seus atributos."""
    graph.add_node(
        node.node_id,
        node_type=node.node_type,
        label=node.label,
        **node.properties,
    )


def _add_edge(graph: nx.DiGraph, edge: KGEdge) -> None:
    """Adiciona uma aresta ao grafo com tipo e propriedades."""
    graph.add_edge(
        edge.source_id,
        edge.target_id,
        edge_type=edge.edge_type.value,
        **edge.properties,
    )


def _load_organisms(graph: nx.DiGraph) -> list[OrganismNode]:
    """Cria nos para os organismos relevantes ao projeto.

    Retorna a lista de nos criados para referencia.
    """
    organisms = [
        OrganismNode(
            taxon_id="5671",
            scientific_name="Leishmania infantum",
            common_name="L. infantum",
            role="parasita",
        ),
        OrganismNode(
            taxon_id="9606",
            scientific_name="Homo sapiens",
            common_name="Humano",
            role="hospedeiro",
        ),
        OrganismNode(
            taxon_id="9615",
            scientific_name="Canis lupus familiaris",
            common_name="Cao domestico",
            role="hospedeiro_reservatorio",
        ),
    ]

    for org in organisms:
        _add_node(graph, org)

    # Relacao de infeccao: L. infantum infecta cao e humano
    _add_edge(graph, KGEdge(
        source_id="organism:5671",
        target_id="organism:9615",
        edge_type=EdgeType.INFECTS,
        properties={"disease": "Leishmaniose visceral canina"},
    ))
    _add_edge(graph, KGEdge(
        source_id="organism:5671",
        target_id="organism:9606",
        edge_type=EdgeType.INFECTS,
        properties={"disease": "Leishmaniose visceral humana (calazar)"},
    ))

    return organisms


def _load_epitopes(graph: nx.DiGraph) -> int:
    """Carrega os 11 epitopos do pipeline Marley e seus genes de origem.

    Cria nos EpitopeNode e GeneNode, conecta com arestas DERIVED_FROM
    (epitopo derivado do gene) e ENCODES (gene codifica proteina).

    Retorna o numero de epitopos carregados.
    """
    # Importacao local para evitar dependencia circular
    from vaccine_platforms.shared.epitopes import EPITOPES

    # Rastreio de genes ja criados (evitar duplicatas)
    genes_created: set[str] = set()

    for ep in EPITOPES:
        # Criar no do epitopo
        epitope_node = EpitopeNode(
            peptide=ep.peptide,
            allele=ep.allele,
            ic50=ep.ic50,
            rank=ep.rank,
            position=ep.position,
            gene_id=ep.gene_id,
        )
        _add_node(graph, epitope_node)

        # Criar no do gene de origem (se ainda nao existe)
        gene_node_id = f"gene:{ep.gene_id}"
        if ep.gene_id not in genes_created:
            # Extrair cromossomo do gene_id (ex: LINF_240013900 -> chr24)
            chromosome = ""
            if ep.gene_id.startswith("LINF_"):
                chromosome = f"chr{ep.gene_id[5:7]}"

            gene_node = GeneNode(
                gene_id=ep.gene_id,
                gene_name=ep.gene_name,
                chromosome=chromosome,
                organism="L. infantum",
            )
            _add_node(graph, gene_node)
            genes_created.add(ep.gene_id)

            # Gene pertence ao organismo L. infantum
            _add_edge(graph, KGEdge(
                source_id=gene_node_id,
                target_id="organism:5671",
                edge_type=EdgeType.EXPRESSED_IN,
            ))

        # Epitopo derivado do gene
        _add_edge(graph, KGEdge(
            source_id=epitope_node.node_id,
            target_id=gene_node_id,
            edge_type=EdgeType.DERIVED_FROM,
            properties={"position": ep.position},
        ))

    return len(EPITOPES)


def _load_drug_targets(graph: nx.DiGraph) -> int:
    """Carrega os 7 alvos farmacologicos validados e suas vias metabolicas.

    Cria nos ProteinNode e PathwayNode, conecta com arestas PART_OF
    (proteina participa da via) e EXPRESSED_IN (no organismo).

    Retorna o numero de alvos carregados.
    """
    # Rastreio de vias ja criadas
    pathways_created: set[str] = set()

    for target in _DRUG_TARGETS:
        # Criar no da proteina/enzima
        protein_node = ProteinNode(
            protein_id=target["gene_id"],
            protein_name=target["gene_name"],
            enzyme_class=target["enzyme_class"],
            ec_number=target["ec_number"],
            druggability_score=target["druggability_score"],
        )
        _add_node(graph, protein_node)

        # Criar no da via metabolica (se ainda nao existe)
        pathway_id = target["pathway"]
        if pathway_id not in pathways_created:
            pathway_node = PathwayNode(
                pathway_id=pathway_id,
                pathway_name=pathway_id.replace("_", " ").title(),
                description=_PATHWAY_DESCRIPTIONS.get(pathway_id, ""),
            )
            _add_node(graph, pathway_node)
            pathways_created.add(pathway_id)

        # Proteina participa da via
        _add_edge(graph, KGEdge(
            source_id=protein_node.node_id,
            target_id=f"pathway:{pathway_id}",
            edge_type=EdgeType.PART_OF,
            properties={"role": target["enzyme_class"]},
        ))

        # Proteina expressa em L. infantum
        _add_edge(graph, KGEdge(
            source_id=protein_node.node_id,
            target_id="organism:5671",
            edge_type=EdgeType.EXPRESSED_IN,
        ))

    return len(_DRUG_TARGETS)


def _load_aso_data(graph: nx.DiGraph) -> None:
    """Carrega o ASO MRL-ASO-001 e o SL RNA alvo do aso_math.

    Cria nos ASONode e RNANode, conecta com arestas TARGETS e BINDS.
    Tambem cria a via de trans-splicing e conecta o SL RNA a ela.
    """
    from aso_math.config import (
        ASO_KNOWN_DG,
        ASO_KNOWN_GC,
        ASO_KNOWN_TM,
        ASO_SEQUENCE,
        ASO_TARGET_END,
        ASO_TARGET_START,
        SL_RNA_COPY_NUMBER,
        SL_SEQUENCE,
    )

    # No do SL RNA
    sl_rna = RNANode(
        rna_id="SL_RNA_LINF",
        rna_name="Spliced Leader RNA (L. infantum)",
        sequence=SL_SEQUENCE,
        rna_type="SL",
        copy_number=SL_RNA_COPY_NUMBER,
    )
    _add_node(graph, sl_rna)

    # No do ASO candidato
    aso = ASONode(
        aso_id="MRL-ASO-001",
        sequence=ASO_SEQUENCE,
        target_start=ASO_TARGET_START,
        target_end=ASO_TARGET_END,
        tm=ASO_KNOWN_TM,
        dg=ASO_KNOWN_DG,
        gc_content=ASO_KNOWN_GC,
    )
    _add_node(graph, aso)

    # ASO alvo o SL RNA
    _add_edge(graph, KGEdge(
        source_id=aso.node_id,
        target_id=sl_rna.node_id,
        edge_type=EdgeType.TARGETS,
        properties={
            "mechanism": "Antisense hybridization + RNase H cleavage",
            "target_region": f"{ASO_TARGET_START}-{ASO_TARGET_END}",
        },
    ))

    # ASO se liga ao SL RNA
    _add_edge(graph, KGEdge(
        source_id=aso.node_id,
        target_id=sl_rna.node_id,
        edge_type=EdgeType.BINDS,
        properties={"dg_kcal_mol": ASO_KNOWN_DG, "tm_celsius": ASO_KNOWN_TM},
    ))

    # ASO inibe o SL RNA (consequencia funcional)
    _add_edge(graph, KGEdge(
        source_id=aso.node_id,
        target_id=sl_rna.node_id,
        edge_type=EdgeType.INHIBITS,
        properties={"mechanism": "Blocks trans-splicing of all mRNAs"},
    ))

    # Via de trans-splicing
    ts_pathway = PathwayNode(
        pathway_id="sl_rna_trans_splicing",
        pathway_name="SL RNA Trans-Splicing",
        description=_PATHWAY_DESCRIPTIONS.get("sl_rna_trans_splicing", ""),
    )
    _add_node(graph, ts_pathway)

    # SL RNA participa da via de trans-splicing
    _add_edge(graph, KGEdge(
        source_id=sl_rna.node_id,
        target_id=ts_pathway.node_id,
        edge_type=EdgeType.PART_OF,
    ))

    # SL RNA expresso em L. infantum
    _add_edge(graph, KGEdge(
        source_id=sl_rna.node_id,
        target_id="organism:5671",
        edge_type=EdgeType.EXPRESSED_IN,
    ))


def _load_drugs(graph: nx.DiGraph) -> int:
    """Carrega farmacos conhecidos contra leishmaniose.

    Cria nos DrugNode e conecta com as vias metabolicas que alvejam.

    Retorna o numero de farmacos carregados.
    """
    for drug in _KNOWN_DRUGS:
        drug_node = DrugNode(
            drug_id=drug["drug_id"],
            drug_name=drug["drug_name"],
            mechanism=drug["mechanism"],
            status=drug["status"],
        )
        _add_node(graph, drug_node)

        # Farmaco alvo a via metabolica
        pathway_id = drug["targets_pathway"]
        pathway_node_id = f"pathway:{pathway_id}"
        if graph.has_node(pathway_node_id):
            _add_edge(graph, KGEdge(
                source_id=drug_node.node_id,
                target_id=pathway_node_id,
                edge_type=EdgeType.TARGETS,
                properties={"mechanism": drug["mechanism"]},
            ))

    return len(_KNOWN_DRUGS)


def _load_papers(graph: nx.DiGraph) -> int:
    """Carrega publicacoes-chave e conecta com entidades relevantes.

    Retorna o numero de papers carregados.
    """
    for paper in _KEY_PAPERS:
        from marley_ai.02_leish_kg.schema import PaperNode
        paper_node = PaperNode(
            paper_id=paper["paper_id"],
            title=paper["title"],
            doi=paper["doi"],
            year=paper["year"],
            authors=paper["authors"],
        )
        _add_node(graph, paper_node)

    # Conectar papers com entidades que validam

    # Crooke 2017 — valida a abordagem ASO
    if graph.has_node("aso:MRL-ASO-001"):
        _add_edge(graph, KGEdge(
            source_id="aso:MRL-ASO-001",
            target_id="paper:crooke2017",
            edge_type=EdgeType.VALIDATED_BY,
            properties={"context": "Referencia para design de ASOs terapeuticos"},
        ))

    # Liang 2003 — valida o SL RNA como alvo
    if graph.has_node("rna:SL_RNA_LINF"):
        _add_edge(graph, KGEdge(
            source_id="rna:SL_RNA_LINF",
            target_id="paper:liang2003",
            edge_type=EdgeType.VALIDATED_BY,
            properties={"context": "Numero de copias do SL RNA (~150)"},
        ))

    # SantaLucia 1998 — valida os calculos termodinamicos
    if graph.has_node("aso:MRL-ASO-001"):
        _add_edge(graph, KGEdge(
            source_id="aso:MRL-ASO-001",
            target_id="paper:santalucia1998",
            edge_type=EdgeType.VALIDATED_BY,
            properties={"context": "Parametros NN para calculo de Tm e dG"},
        ))

    return len(_KEY_PAPERS)


def _wire_cross_references(graph: nx.DiGraph) -> int:
    """Conecta entidades entre os diferentes dominios do grafo.

    Cria arestas adicionais que ligam:
        - Genes de epitopos com proteinas (ENCODES)
        - Vias metabolicas com o organismo

    Retorna o numero de arestas criadas.
    """
    edges_added = 0

    # Genes que codificam proteinas reconhecidas como alvos de epitopos
    # (os genes dos epitopos codificam proteinas que sao apresentadas via MHC)
    gene_nodes = [
        nid for nid in graph.nodes()
        if graph.nodes[nid].get("node_type") == "Gene"
    ]

    for gene_nid in gene_nodes:
        gene_name = graph.nodes[gene_nid].get("gene_name", "")
        # Criar no de proteina para cada gene de epitopo
        # (a proteina e o produto do gene, e o epitopo e um fragmento dela)
        protein_id = f"product_{gene_nid.replace('gene:', '')}"
        protein_node_id = f"protein:{protein_id}"

        if not graph.has_node(protein_node_id):
            protein_node = ProteinNode(
                protein_id=protein_id,
                protein_name=f"{gene_name} (proteina)",
            )
            _add_node(graph, protein_node)

            # Gene codifica proteina
            _add_edge(graph, KGEdge(
                source_id=gene_nid,
                target_id=protein_node_id,
                edge_type=EdgeType.ENCODES,
            ))
            edges_added += 1

            # Conectar epitopos a esta proteina (epitopo se liga a proteina)
            for ep_nid in graph.predecessors(gene_nid):
                ep_data = graph.nodes.get(ep_nid, {})
                if ep_data.get("node_type") == "Epitope":
                    _add_edge(graph, KGEdge(
                        source_id=ep_nid,
                        target_id=protein_node_id,
                        edge_type=EdgeType.BINDS,
                        properties={"context": "Epitopo derivado desta proteina"},
                    ))
                    edges_added += 1

    # Vias metabolicas pertencem ao parasita
    pathway_nodes = [
        nid for nid in graph.nodes()
        if graph.nodes[nid].get("node_type") == "Pathway"
    ]
    for pw_nid in pathway_nodes:
        if not graph.has_edge(pw_nid, "organism:5671"):
            _add_edge(graph, KGEdge(
                source_id=pw_nid,
                target_id="organism:5671",
                edge_type=EdgeType.PART_OF,
                properties={"context": "Via metabolica do parasita"},
            ))
            edges_added += 1

    return edges_added


# ---------------------------------------------------------------------------
# API publica
# ---------------------------------------------------------------------------

def build_knowledge_graph() -> nx.DiGraph:
    """Constroi o Knowledge Graph completo do projeto Marley.

    Carrega todas as entidades biologicas do projeto e conecta
    com arestas tipadas. O grafo resultante e um DiGraph do NetworkX
    com atributos em cada no e aresta.

    Returns:
        Grafo dirigido NetworkX com o KG completo.
    """
    graph = nx.DiGraph()

    # 1. Organismos (base do grafo)
    _load_organisms(graph)

    # 2. Epitopos e genes de origem
    n_epitopes = _load_epitopes(graph)

    # 3. Alvos farmacologicos e vias metabolicas
    n_targets = _load_drug_targets(graph)

    # 4. ASO e SL RNA
    _load_aso_data(graph)

    # 5. Farmacos conhecidos
    n_drugs = _load_drugs(graph)

    # 6. Publicacoes de referencia
    n_papers = _load_papers(graph)

    # 7. Conexoes cruzadas entre dominios
    n_cross = _wire_cross_references(graph)

    print(
        f"[02_leish_kg] Grafo construido: "
        f"{graph.number_of_nodes()} nos, {graph.number_of_edges()} arestas"
    )
    print(
        f"  Fontes: {n_epitopes} epitopos, {n_targets} alvos, "
        f"{n_drugs} farmacos, {n_papers} papers, {n_cross} conexoes cruzadas"
    )

    return graph
