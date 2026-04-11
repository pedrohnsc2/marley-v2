"""Schemas de dados para comunicacao entre modulos AI.

Define contratos imutaveis (frozen dataclasses) para os tipos de dados
que fluem entre modulos. Garante tipagem forte e validacao na fronteira
entre modulos, evitando erros silenciosos.

Cada schema corresponde a um tipo de artefato biologico que os modulos
produzem e consomem.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class ProteinEmbedding:
    """Embedding de uma proteina gerado pelo ESM-2 ou modelo similar.

    Usado pelos modulos 03_leish_esm, 07_contrastive e 10_digital_twin.
    O embedding captura propriedades estruturais e funcionais da proteina
    em um vetor denso de alta dimensionalidade.
    """
    protein_id: str              # Identificador no proteoma (ex: LINF_240013900)
    gene_name: str               # Nome do gene anotado
    sequence: str                # Sequencia aminoacidica completa
    embedding_dim: int           # Dimensionalidade do embedding (ex: 1280 para ESM-2)
    model_name: str              # Modelo que gerou (ex: "esm2_t33_650M_UR50D")
    embedding_path: str          # Caminho para o .npy com o vetor


@dataclass(frozen=True)
class RNAStructure:
    """Estrutura predita de RNA (2D ou 3D).

    Usado pelos modulos 04_rna_fm e 05_rosettafold.
    O SL RNA de Leishmania forma hairpins criticos para o trans-splicing;
    a estrutura 3D informa o design de ASOs que desestabilizam o alvo.
    """
    rna_id: str                  # Identificador (ex: "sl_rna_linf")
    sequence: str                # Sequencia nucleotidica
    dot_bracket: str             # Estrutura secundaria em notacao dot-bracket
    free_energy: float           # Energia livre de Gibbs da estrutura (kcal/mol)
    structure_path: str          # Caminho para o .pdb (se 3D disponivel)
    model_name: str              # Modelo que gerou (ex: "RoseTTAFold2NA")
    confidence: float            # pLDDT medio ou score de confianca


@dataclass(frozen=True)
class KnowledgeTriple:
    """Tripla do Knowledge Graph de Leishmania.

    Usado pelo modulo 02_leish_kg. Cada tripla conecta duas entidades
    biologicas por uma relacao, formando o grafo de conhecimento que
    alimenta o RAG e o agente cientifico.

    Exemplo: (GP63, "cliva", "Complemento C3b")
    """
    subject: str                 # Entidade de origem
    predicate: str               # Tipo de relacao
    obj: str                     # Entidade de destino (object)
    source: str                  # Referencia bibliografica ou database
    confidence: float = 1.0      # Score de confianca (1.0 = manual, <1 = predito)


@dataclass(frozen=True)
class EpitopeScore:
    """Score de um epitopo gerado por modelo contrastivo ou RL.

    Usado pelos modulos 07_contrastive e 08_rl_ppo.
    Combina afinidade de ligacao predita com imunogenicidade estimada
    para ranking de candidatos vacinais.
    """
    peptide: str                 # Sequencia do epitopo
    allele: str                  # Alelo MHC/DLA
    binding_score: float         # Score de afinidade (quanto maior, melhor)
    immunogenicity: float        # Score de imunogenicidade estimada
    source_gene: str             # Gene de origem no parasita
    model_name: str              # Modelo que gerou o score


@dataclass(frozen=True)
class DesignCandidate:
    """Candidato gerado por design computacional (difusao ou RL).

    Usado pelos modulos 06_evodiff e 08_rl_ppo.
    Representa uma sequencia (proteina ou ASO) gerada de novo ou
    otimizada por modelos generativos.
    """
    candidate_id: str            # Identificador unico
    sequence: str                # Sequencia gerada
    sequence_type: str           # "protein" ou "nucleic_acid"
    score: float                 # Score de qualidade do modelo gerador
    properties: dict[str, Any] = field(default_factory=dict)  # Propriedades preditas
    generation_method: str = ""  # Metodo de geracao (ex: "evodiff_oa_dm_38M")
