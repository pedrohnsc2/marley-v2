"""Ontologia do Knowledge Graph — tipos de nos e arestas.

Define dataclasses imutaveis para cada tipo de no biologico no grafo,
e um Enum para os tipos de relacao (arestas). Cada no possui um ID unico,
um tipo, e um dicionario de propriedades especificas.

A validacao garante que todo no inserido no grafo segue o esquema
definido aqui, prevenindo dados malformados e facilitando queries
tipadas sobre o grafo.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, unique
from typing import Any


# ---------------------------------------------------------------------------
# Tipos de aresta (relacoes biologicas)
# ---------------------------------------------------------------------------

@unique
class EdgeType(Enum):
    """Tipos de relacao entre entidades biologicas no KG.

    Cada aresta do grafo deve ter exatamente um destes tipos.
    A semantica segue convencoes de ontologias biomedicas (GO, SIO).
    """

    ENCODES = "ENCODES"              # Gene -> Protein
    INHIBITS = "INHIBITS"            # ASO/Drug -> Protein/RNA
    PART_OF = "PART_OF"              # Protein -> Pathway, Gene -> Organism
    VALIDATED_BY = "VALIDATED_BY"    # Protein/Epitope -> Paper
    BINDS = "BINDS"                  # Epitope -> Protein, ASO -> RNA
    ACTIVATES = "ACTIVATES"          # Protein -> Pathway
    TARGETS = "TARGETS"              # ASO -> RNA, Drug -> Protein
    DERIVED_FROM = "DERIVED_FROM"    # Epitope -> Gene
    EXPRESSED_IN = "EXPRESSED_IN"    # Gene -> Organism
    INFECTS = "INFECTS"              # Organism -> Organism


# ---------------------------------------------------------------------------
# Classe base para nos
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class KGNode:
    """No base do Knowledge Graph.

    Todo no deve ter um ID unico (node_id), um tipo (node_type),
    e um rotulo legivel (label). Propriedades adicionais ficam
    no dicionario `properties`.
    """

    node_id: str
    node_type: str
    label: str
    properties: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Valida campos obrigatorios."""
        if not self.node_id:
            raise ValueError("node_id nao pode ser vazio")
        if not self.node_type:
            raise ValueError("node_type nao pode ser vazio")
        if not self.label:
            raise ValueError("label nao pode ser vazio")

    def to_dict(self) -> dict[str, Any]:
        """Serializa o no para dicionario (compativel com JSON)."""
        return {
            "node_id": self.node_id,
            "node_type": self.node_type,
            "label": self.label,
            "properties": dict(self.properties),
        }


# ---------------------------------------------------------------------------
# Tipos especificos de nos — cada um com propriedades relevantes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GeneNode(KGNode):
    """No representando um gene de L. infantum ou outro organismo.

    Campos:
        gene_id: Identificador no genoma (ex: LINF_240013900)
        gene_name: Nome do produto genico anotado
        chromosome: Cromossomo onde o gene esta localizado
        organism: Organismo de origem
    """

    def __init__(
        self,
        gene_id: str,
        gene_name: str,
        *,
        chromosome: str = "",
        organism: str = "L. infantum",
    ) -> None:
        props = {
            "gene_id": gene_id,
            "gene_name": gene_name,
            "chromosome": chromosome,
            "organism": organism,
        }
        # Frozen dataclass — usamos object.__setattr__ via super().__init__
        super().__init__(
            node_id=f"gene:{gene_id}",
            node_type="Gene",
            label=gene_name,
            properties=props,
        )


@dataclass(frozen=True)
class ProteinNode(KGNode):
    """No representando uma proteina (enzima, receptor, etc).

    Campos:
        protein_id: UniProt accession ou ID interno
        protein_name: Nome da proteina
        enzyme_class: Classe enzimatica (se aplicavel)
        ec_number: Numero EC (se aplicavel)
        druggability_score: Score de druggabilidade (0-1)
    """

    def __init__(
        self,
        protein_id: str,
        protein_name: str,
        *,
        enzyme_class: str = "",
        ec_number: str = "",
        druggability_score: float = 0.0,
    ) -> None:
        props = {
            "protein_id": protein_id,
            "protein_name": protein_name,
            "enzyme_class": enzyme_class,
            "ec_number": ec_number,
            "druggability_score": druggability_score,
        }
        super().__init__(
            node_id=f"protein:{protein_id}",
            node_type="Protein",
            label=protein_name,
            properties=props,
        )


@dataclass(frozen=True)
class EpitopeNode(KGNode):
    """No representando um epitopo predito pelo pipeline Marley.

    Campos:
        peptide: Sequencia aminoacidica do epitopo (9-mer)
        allele: Alelo DLA que apresenta o epitopo
        ic50: Afinidade de ligacao predita (nM)
        rank: Rank percentual no NetMHCpan
        position: Posicao no gene de origem
    """

    def __init__(
        self,
        peptide: str,
        *,
        allele: str = "",
        ic50: float = 0.0,
        rank: float = 0.0,
        position: int = 0,
        gene_id: str = "",
    ) -> None:
        props = {
            "peptide": peptide,
            "allele": allele,
            "ic50": ic50,
            "rank": rank,
            "position": position,
            "gene_id": gene_id,
        }
        super().__init__(
            node_id=f"epitope:{peptide}",
            node_type="Epitope",
            label=f"Epitope {peptide}",
            properties=props,
        )


@dataclass(frozen=True)
class ASONode(KGNode):
    """No representando um oligonucleotideo antisenso (ASO).

    Campos:
        sequence: Sequencia do ASO (DNA)
        target_start: Posicao inicial no alvo
        target_end: Posicao final no alvo
        tm: Temperatura de melting predita (Celsius)
        dg: Energia livre de Gibbs (kcal/mol)
        gc_content: Conteudo GC fracional
    """

    def __init__(
        self,
        aso_id: str,
        sequence: str,
        *,
        target_start: int = 0,
        target_end: int = 0,
        tm: float = 0.0,
        dg: float = 0.0,
        gc_content: float = 0.0,
    ) -> None:
        props = {
            "sequence": sequence,
            "length": len(sequence),
            "target_start": target_start,
            "target_end": target_end,
            "tm": tm,
            "dg": dg,
            "gc_content": gc_content,
        }
        super().__init__(
            node_id=f"aso:{aso_id}",
            node_type="ASO",
            label=f"ASO {aso_id}",
            properties=props,
        )


@dataclass(frozen=True)
class RNANode(KGNode):
    """No representando um RNA (SL RNA, mRNA, etc).

    Campos:
        sequence: Sequencia nucleotidica
        rna_type: Tipo de RNA (SL, mRNA, rRNA, etc.)
        copy_number: Numero de copias no genoma
    """

    def __init__(
        self,
        rna_id: str,
        rna_name: str,
        *,
        sequence: str = "",
        rna_type: str = "SL",
        copy_number: int = 0,
    ) -> None:
        props = {
            "sequence": sequence,
            "rna_type": rna_type,
            "length": len(sequence),
            "copy_number": copy_number,
        }
        super().__init__(
            node_id=f"rna:{rna_id}",
            node_type="RNA",
            label=rna_name,
            properties=props,
        )


@dataclass(frozen=True)
class PathwayNode(KGNode):
    """No representando uma via metabolica.

    Campos:
        pathway_id: Identificador da via (slug)
        pathway_name: Nome legivel da via
        description: Descricao da via
    """

    def __init__(
        self,
        pathway_id: str,
        pathway_name: str,
        *,
        description: str = "",
    ) -> None:
        props = {
            "pathway_id": pathway_id,
            "pathway_name": pathway_name,
            "description": description,
        }
        super().__init__(
            node_id=f"pathway:{pathway_id}",
            node_type="Pathway",
            label=pathway_name,
            properties=props,
        )


@dataclass(frozen=True)
class OrganismNode(KGNode):
    """No representando um organismo (parasita, hospedeiro, vetor).

    Campos:
        taxon_id: ID taxonomico (NCBI taxonomy)
        scientific_name: Nome cientifico
        common_name: Nome comum
        role: Papel no contexto da doenca (parasita, hospedeiro, vetor)
    """

    def __init__(
        self,
        taxon_id: str,
        scientific_name: str,
        *,
        common_name: str = "",
        role: str = "",
    ) -> None:
        props = {
            "taxon_id": taxon_id,
            "scientific_name": scientific_name,
            "common_name": common_name,
            "role": role,
        }
        super().__init__(
            node_id=f"organism:{taxon_id}",
            node_type="Organism",
            label=scientific_name,
            properties=props,
        )


@dataclass(frozen=True)
class DrugNode(KGNode):
    """No representando um farmaco ou composto quimico.

    Campos:
        drug_name: Nome do farmaco
        mechanism: Mecanismo de acao
        status: Status de desenvolvimento (aprovado, experimental, etc.)
    """

    def __init__(
        self,
        drug_id: str,
        drug_name: str,
        *,
        mechanism: str = "",
        status: str = "experimental",
    ) -> None:
        props = {
            "drug_name": drug_name,
            "mechanism": mechanism,
            "status": status,
        }
        super().__init__(
            node_id=f"drug:{drug_id}",
            node_type="Drug",
            label=drug_name,
            properties=props,
        )


@dataclass(frozen=True)
class PaperNode(KGNode):
    """No representando uma publicacao cientifica.

    Campos:
        doi: Digital Object Identifier
        title: Titulo do artigo
        year: Ano de publicacao
        authors: Lista de autores (primeiro + et al.)
    """

    def __init__(
        self,
        paper_id: str,
        title: str,
        *,
        doi: str = "",
        year: int = 0,
        authors: str = "",
    ) -> None:
        props = {
            "doi": doi,
            "title": title,
            "year": year,
            "authors": authors,
        }
        super().__init__(
            node_id=f"paper:{paper_id}",
            node_type="Paper",
            label=title,
            properties=props,
        )


# ---------------------------------------------------------------------------
# Aresta tipada
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class KGEdge:
    """Aresta tipada do Knowledge Graph.

    Conecta dois nos com uma relacao semantica (EdgeType).
    Propriedades adicionais (peso, evidencia, etc.) ficam em `properties`.
    """

    source_id: str
    target_id: str
    edge_type: EdgeType
    properties: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Valida campos obrigatorios."""
        if not self.source_id:
            raise ValueError("source_id nao pode ser vazio")
        if not self.target_id:
            raise ValueError("target_id nao pode ser vazio")
        if not isinstance(self.edge_type, EdgeType):
            raise TypeError(
                f"edge_type deve ser EdgeType, recebeu {type(self.edge_type)}"
            )

    def to_dict(self) -> dict[str, Any]:
        """Serializa a aresta para dicionario (compativel com JSON)."""
        return {
            "source_id": self.source_id,
            "target_id": self.target_id,
            "edge_type": self.edge_type.value,
            "properties": dict(self.properties),
        }
