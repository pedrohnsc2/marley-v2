"""Encoders de sequencia para epitopos e alelos DLA.

Implementa encodings numericos para representar peptideos e alelos
como vetores de features, que serao inputs do MLP contrastivo.

Tres estrategias de encoding:
    1. k-mer: frequencia de subsequencias de tamanho k
    2. AAindex: propriedades fisicoquimicas dos 20 aminoacidos
    3. Combinado: concatenacao de k-mer + fisicoquimico

Para alelos DLA, usa-se encoding baseado em propriedades funcionais
conhecidas (especificidade de ancoragem, tamanho do bolso, etc).
"""

from __future__ import annotations

from itertools import product
from typing import Final

import numpy as np


# ---------------------------------------------------------------------------
# Tabela de propriedades fisicoquimicas dos 20 aminoacidos padrao.
#
# Propriedades (normalizadas para ~[-1, 1]):
#   0: Hidrofobicidade (escala Kyte-Doolittle, normalizada)
#   1: Carga a pH 7.0 (+1, 0, -1)
#   2: Volume (Angstrom^3, normalizado)
#   3: Polaridade (0 = apolar, 1 = polar neutro, 2 = polar carregado, norm.)
#   4: Flexibilidade (B-factor relativo, normalizado)
#
# Fontes: Kyte & Doolittle (1982), Zimmerman et al. (1968),
#          Bhaskaran & Ponnuswamy (1988)
# ---------------------------------------------------------------------------

AA_PROPERTIES: Final[dict[str, list[float]]] = {
    # AA:  [hidrofob,  carga,  volume,  polaridade, flexib.]
    "A": [ 0.700,  0.0,  -0.733,  -1.0,  -0.238],
    "R": [-1.000,  1.0,   0.596,   1.0,   0.379],
    "N": [-0.556,  0.0,  -0.244,   0.5,   0.517],
    "D": [-0.556,  -1.0, -0.311,   1.0,   0.586],
    "C": [ 0.389,  0.0,  -0.356,   0.5,  -0.172],
    "E": [-0.556,  -1.0,  0.022,   1.0,   0.345],
    "Q": [-0.556,  0.0,   0.022,   0.5,   0.483],
    "G": [-0.111,  0.0,  -1.000,  -1.0,   0.690],
    "H": [-0.500,  0.5,   0.089,   0.5,   0.103],
    "I": [ 1.000,  0.0,   0.156,  -1.0,  -0.690],
    "L": [ 0.667,  0.0,   0.156,  -1.0,  -0.517],
    "K": [-0.611,  1.0,   0.289,   1.0,   0.414],
    "M": [ 0.333,  0.0,   0.156,  -1.0,  -0.103],
    "F": [ 0.444,  0.0,   0.556,  -1.0,  -0.586],
    "P": [-0.278,  0.0,  -0.378,  -1.0,   0.103],
    "S": [-0.167,  0.0,  -0.622,   0.5,   0.552],
    "T": [-0.056,  0.0,  -0.356,   0.5,   0.172],
    "W": [-0.167,  0.0,   1.000,  -1.0,  -1.000],
    "Y": [-0.222,  0.0,   0.644,   0.5,  -0.345],
    "V": [ 0.778,  0.0,  -0.178,  -1.0,  -0.414],
}

# Aminoacidos padrao em ordem alfabetica (para vocabulario de k-mers)
AMINO_ACIDS: Final[str] = "ACDEFGHIKLMNPQRSTVWY"

# ---------------------------------------------------------------------------
# Propriedades funcionais dos 3 alelos DLA caninos.
#
# Encodings baseados em caracteristicas conhecidas de apresentacao:
#   0: Preferencia por ancoragem hidrofobica na posicao 2 (P2)
#   1: Preferencia por ancoragem C-terminal (P9/Pomega)
#   2: Tamanho relativo do bolso de ligacao
#   3: Especificidade ampla vs restrita (breadth)
#   4: Frequencia populacional relativa em caes
#
# Valores derivados de Venkataraman et al. (2007) e Ross et al. (2012)
# sobre polimorfismo do MHC-I canino.
# ---------------------------------------------------------------------------

DLA_ALLELE_PROPERTIES: Final[dict[str, list[float]]] = {
    # Alelo:          [P2_hid, Pomega, bolso, breadth, freq_pop]
    "DLA-8803401": [ 0.8,   0.7,   0.5,   0.6,     0.3],
    "DLA-8850101": [ 0.6,   0.9,   0.7,   0.8,     0.5],
    "DLA-8850801": [ 0.7,   0.6,   0.6,   0.7,     0.2],
}


def build_kmer_vocab(k: int) -> dict[str, int]:
    """Constroi vocabulario de k-mers para os 20 aminoacidos.

    Args:
        k: tamanho do k-mer (tipicamente 2 ou 3)

    Returns:
        Mapeamento k-mer -> indice no vetor de frequencias
    """
    vocab = {}
    for i, kmer in enumerate(product(AMINO_ACIDS, repeat=k)):
        vocab["".join(kmer)] = i
    return vocab


def encode_kmer(peptide: str, k: int, vocab: dict[str, int] | None = None) -> np.ndarray:
    """Encoding de frequencia de k-mers para um peptideo.

    Conta a frequencia relativa de cada k-mer na sequencia.
    Para peptideos curtos (9-mers), k=3 gera vetores esparsos
    mas informativos sobre padroes locais de sequencia.

    Args:
        peptide: sequencia aminoacidica (ex: "RMMRSLTPF")
        k: tamanho do k-mer
        vocab: vocabulario pre-computado (se None, computa internamente)

    Returns:
        Vetor de frequencias relativas (tamanho = 20^k)
    """
    if vocab is None:
        vocab = build_kmer_vocab(k)

    vec = np.zeros(len(vocab), dtype=np.float64)
    n_kmers = len(peptide) - k + 1

    if n_kmers <= 0:
        return vec

    for i in range(n_kmers):
        kmer = peptide[i:i + k]
        if kmer in vocab:
            vec[vocab[kmer]] += 1.0

    # Normalizar por contagem total de k-mers
    vec /= n_kmers
    return vec


def encode_physicochemical(peptide: str) -> np.ndarray:
    """Encoding fisicoquimico baseado em AAindex.

    Para cada posicao do peptideo, extrai 5 propriedades fisicoquimicas.
    O vetor final e a concatenacao das propriedades de todas as posicoes.

    Para peptideos de tamanho variavel, faz padding com zeros ate o
    tamanho maximo (15 residuos). Os 11 epitopos canonicos sao 9-mers.

    Args:
        peptide: sequencia aminoacidica

    Returns:
        Vetor de propriedades (tamanho = max_len * n_properties)
    """
    max_len = 15  # tamanho maximo de peptideo (com margem)
    n_props = len(next(iter(AA_PROPERTIES.values())))

    vec = np.zeros(max_len * n_props, dtype=np.float64)

    for i, aa in enumerate(peptide[:max_len]):
        props = AA_PROPERTIES.get(aa)
        if props is not None:
            start = i * n_props
            vec[start:start + n_props] = props

    return vec


def encode_peptide(peptide: str, k: int, kmer_vocab: dict[str, int] | None = None) -> np.ndarray:
    """Encoding combinado de um peptideo: k-mer + fisicoquimico.

    Concatena o vetor de frequencia de k-mers com o vetor de
    propriedades fisicoquimicas, gerando uma representacao rica
    que captura tanto padroes de sequencia quanto propriedades
    biofisicas.

    Args:
        peptide: sequencia aminoacidica
        k: tamanho do k-mer
        kmer_vocab: vocabulario pre-computado

    Returns:
        Vetor concatenado [kmer_features | physicochemical_features]
    """
    kmer_feat = encode_kmer(peptide, k, kmer_vocab)
    physchem_feat = encode_physicochemical(peptide)
    return np.concatenate([kmer_feat, physchem_feat])


def encode_allele(allele: str) -> np.ndarray:
    """Encoding de um alelo DLA baseado em propriedades funcionais.

    Cada alelo DLA tem um vetor de 5 propriedades que descrevem
    suas caracteristicas de ligacao e apresentacao antigenica.

    Para alelos desconhecidos, retorna vetor zero com warning.

    Args:
        allele: nome do alelo DLA (ex: "DLA-8850101")

    Returns:
        Vetor de propriedades do alelo (tamanho = 5)
    """
    props = DLA_ALLELE_PROPERTIES.get(allele)
    if props is None:
        # Alelo desconhecido — retorna vetor zero
        return np.zeros(5, dtype=np.float64)
    return np.array(props, dtype=np.float64)


def get_peptide_feature_dim(k: int) -> int:
    """Calcula a dimensao total do encoding combinado de peptideo.

    Util para inicializar os pesos do MLP antes de ver dados.

    Args:
        k: tamanho do k-mer

    Returns:
        Dimensao total = 20^k (k-mer) + 15*5 (fisicoquimico)
    """
    kmer_dim = len(AMINO_ACIDS) ** k
    physchem_dim = 15 * 5  # max_len * n_properties
    return kmer_dim + physchem_dim


def get_allele_feature_dim() -> int:
    """Retorna a dimensao do encoding de alelo.

    Returns:
        Dimensao = 5 (propriedades funcionais)
    """
    return 5
