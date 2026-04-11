"""Validadores biofisicos e de sequencia para construtos vacinais.

Implementacao em Python puro (sem dependencias externas) dos principais
parametros usados na avaliacao de candidatos vacinais recombinantes:

- Peso molecular (soma dos residuos menos agua de condensacao)
- Ponto isoeletrico (metodo de Henderson-Hasselbalch iterativo)
- Indice de instabilidade (metodo de Guruprasad, DIWV weights)
- GRAVY (media da hidrofobicidade de Kyte-Doolittle)
- Conteudo GC (para sequencias de DNA/mRNA)
- CAI (Codon Adaptation Index)
- Busca por sitios de restricao indesejados

Todos os valores de referencia (pesos moleculares, pKa, hidrofobicidade)
sao extraidos de tabelas classicas da bioquimica -- Lehninger, Creighton,
e artigos originais citados inline.
"""

from __future__ import annotations

import math
from typing import Final

# ===========================================================================
# Tabelas de propriedades de aminoacidos
# ===========================================================================

# Peso molecular de cada aminoacido em Daltons (massa do residuo + H2O).
# Fonte: NIST / Lehninger Principles of Biochemistry, 8th ed.
# Para calcular peso do peptideo: soma dos residuos - (n-1) * 18.015 (agua)
AMINO_ACID_MW: Final[dict[str, float]] = {
    "A": 89.094,   # Alanina
    "R": 174.203,  # Arginina
    "N": 132.119,  # Asparagina
    "D": 133.104,  # Acido aspartico
    "C": 121.154,  # Cisteina
    "E": 147.130,  # Acido glutamico
    "Q": 146.146,  # Glutamina
    "G": 75.032,   # Glicina
    "H": 155.156,  # Histidina
    "I": 131.175,  # Isoleucina
    "L": 131.175,  # Leucina
    "K": 146.189,  # Lisina
    "M": 149.208,  # Metionina
    "F": 165.192,  # Fenilalanina
    "P": 115.132,  # Prolina
    "S": 105.093,  # Serina
    "T": 119.119,  # Treonina
    "W": 204.228,  # Triptofano
    "Y": 181.191,  # Tirosina
    "V": 117.148,  # Valina
}

# Massa da molecula de agua removida na ligacao peptidica
_WATER_MW: Final[float] = 18.015

# Escala de hidrofobicidade de Kyte-Doolittle (1982).
# Valores positivos = hidrofobico, negativos = hidrofilico.
# Usada para calcular GRAVY (Grand Average of Hydropathy).
KYTE_DOOLITTLE: Final[dict[str, float]] = {
    "A":  1.800,
    "R": -4.500,
    "N": -3.500,
    "D": -3.500,
    "C":  2.500,
    "E": -3.500,
    "Q": -3.500,
    "G": -0.400,
    "H": -3.200,
    "I":  4.500,
    "L":  3.800,
    "K": -3.900,
    "M":  1.900,
    "F":  2.800,
    "P": -1.600,
    "S": -0.800,
    "T": -0.700,
    "W": -0.900,
    "Y": -1.300,
    "V":  4.200,
}

# Valores de pKa para calculo do ponto isoeletrico.
# Fonte: Lehninger / EMBOSS / DTASelect.
# Estrutura: {grupo: pKa}
#   - N-terminal (alfa-amino) e C-terminal (alfa-carboxila)
#   - Cadeias laterais ionizaveis: D, E, C, Y, H, K, R
_PKA_VALUES: Final[dict[str, float]] = {
    "N_TERM": 9.69,    # alfa-amino do N-terminal
    "C_TERM": 2.34,    # alfa-carboxila do C-terminal
    "D": 3.90,         # Aspartato (cadeia lateral)
    "E": 4.07,         # Glutamato (cadeia lateral)
    "C": 8.18,         # Cisteina (tiol)
    "Y": 10.46,        # Tirosina (fenol)
    "H": 6.04,         # Histidina (imidazol)
    "K": 10.54,        # Lisina (epsilon-amino)
    "R": 12.48,        # Arginina (guanidinio)
}

# Tabela DIWV (Dipeptide Instability Weight Values) de Guruprasad et al. (1990).
# Proteinas com indice de instabilidade > 40 sao classificadas como instaveis.
# A tabela contem pesos para todos os 400 dipeptideos possiveis.
# Referencia: Guruprasad K, Reddy BV, Pandit MW. Protein Eng. 1990;4(2):155-61.
_DIWV: Final[dict[str, dict[str, float]]] = {
    "A": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 1.0, "C": 44.94, "E": 1.0, "Q": 1.0, "G": 1.0, "H": -7.49, "I": 1.0, "L": 1.0, "K": 1.0, "M": 1.0, "F": 1.0, "P": 20.26, "S": 1.0, "T": 1.0, "W": 1.0, "Y": 1.0, "V": 1.0},
    "R": {"A":  1.0, "R": 58.28, "N": 1.0, "D": 1.0, "C": 1.0, "E": 1.0, "Q": 1.0, "G": -7.49, "H": 1.0, "I": 1.0, "L": 1.0, "K": 1.0, "M": 1.0, "F": 1.0, "P": 20.26, "S": 1.0, "T": 1.0, "W": 58.28, "Y": 1.0, "V": 1.0},
    "N": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 1.0, "C": -1.88, "E": 1.0, "Q": -6.54, "G": -7.49, "H": 1.0, "I": 44.94, "L": 1.0, "K": 24.68, "M": 1.0, "F": -14.03, "P": -1.88, "S": 1.0, "T": -7.49, "W": -9.37, "Y": 1.0, "V": 1.0},
    "D": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 1.0, "C": 1.0, "E": 1.0, "Q": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "L": 1.0, "K": -7.49, "M": 1.0, "F": -6.54, "P": 1.0, "S": 20.26, "T": -14.03, "W": 1.0, "Y": 1.0, "V": 1.0},
    "C": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 20.26, "C": 1.0, "E": 1.0, "Q": -6.54, "G": 1.0, "H": 1.0, "I": 1.0, "L": 20.26, "K": 1.0, "M": 33.60, "F": 1.0, "P": 20.26, "S": 1.0, "T": 33.60, "W": 24.68, "Y": 1.0, "V": -6.54},
    "E": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 20.26, "C": 44.94, "E": 33.60, "Q": 20.26, "G": 1.0, "H": -6.54, "I": 20.26, "L": 1.0, "K": 1.0, "M": 1.0, "F": 1.0, "P": 20.26, "S": 20.26, "T": 1.0, "W": -14.03, "Y": 1.0, "V": 1.0},
    "Q": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 20.26, "C": -6.54, "E": 20.26, "Q": 20.26, "G": 1.0, "H": 1.0, "I": 1.0, "L": 1.0, "K": 1.0, "M": 1.0, "F": -6.54, "P": 20.26, "S": 44.94, "T": 1.0, "W": 1.0, "Y": -6.54, "V": -6.54},
    "G": {"A":  -7.49, "R": 1.0, "N": -7.49, "D": 1.0, "C": 1.0, "E": -6.54, "Q": 1.0, "G": 13.34, "H": 1.0, "I": -7.49, "L": 1.0, "K": -7.49, "M": 1.0, "F": 1.0, "P": 1.0, "S": 1.0, "T": -7.49, "W": 13.34, "Y": -7.49, "V": 1.0},
    "H": {"A":  1.0, "R": 1.0, "N": 24.68, "D": 1.0, "C": 1.0, "E": 1.0, "Q": 1.0, "G": -9.37, "H": 1.0, "I": 44.94, "L": 1.0, "K": 24.68, "M": 1.0, "F": -9.37, "P": -1.88, "S": 1.0, "T": -1.88, "W": -1.88, "Y": 44.94, "V": 1.0},
    "I": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 1.0, "C": 1.0, "E": 44.94, "Q": 1.0, "G": 1.0, "H": 13.34, "I": 1.0, "L": 20.26, "K": 1.0, "M": 1.0, "F": 1.0, "P": -1.88, "S": 1.0, "T": 1.0, "W": 1.0, "Y": 1.0, "V": -7.49},
    "L": {"A":  1.0, "R": 20.26, "N": 1.0, "D": 1.0, "C": 1.0, "E": 1.0, "Q": 33.60, "G": 20.26, "H": 1.0, "I": 1.0, "L": 1.0, "K": -7.49, "M": 1.0, "F": 1.0, "P": 20.26, "S": 1.0, "T": 1.0, "W": 24.68, "Y": 1.0, "V": 1.0},
    "K": {"A":  1.0, "R": 33.60, "N": 1.0, "D": 1.0, "C": 1.0, "E": 1.0, "Q": 24.64, "G": -7.49, "H": 1.0, "I": -7.49, "L": -7.49, "K": 1.0, "M": 33.60, "F": 1.0, "P": -6.54, "S": 1.0, "T": 1.0, "W": 1.0, "Y": 1.0, "V": -7.49},
    "M": {"A":  13.34, "R": -1.88, "N": 1.0, "D": 1.0, "C": 1.0, "E": 1.0, "Q": -6.54, "G": 1.0, "H": 58.28, "I": 1.0, "L": 1.0, "K": 1.0, "M": -1.88, "F": 1.0, "P": 44.94, "S": 44.94, "T": -1.88, "W": 1.0, "Y": 24.68, "V": 1.0},
    "F": {"A":  1.0, "R": 1.0, "N": 1.0, "D": 13.34, "C": 1.0, "E": 1.0, "Q": 1.0, "G": 1.0, "H": 1.0, "I": 1.0, "L": 1.0, "K": -14.03, "M": 1.0, "F": 1.0, "P": 20.26, "S": 1.0, "T": 1.0, "W": 1.0, "Y": 33.60, "V": 1.0},
    "P": {"A":  20.26, "R": -6.54, "N": 1.0, "D": -6.54, "C": -6.54, "E": 18.38, "Q": 20.26, "G": 1.0, "H": 1.0, "I": 1.0, "L": 1.0, "K": 1.0, "M": -6.54, "F": 20.26, "P": 20.26, "S": 20.26, "T": 1.0, "W": -1.88, "Y": 1.0, "V": 20.26},
    "S": {"A":  1.0, "R": 20.26, "N": 1.0, "D": 1.0, "C": 33.60, "E": 20.26, "Q": 20.26, "G": 1.0, "H": 1.0, "I": 1.0, "L": 1.0, "K": 1.0, "M": 1.0, "F": 1.0, "P": 44.94, "S": 20.26, "T": 1.0, "W": 1.0, "Y": 1.0, "V": 1.0},
    "T": {"A":  1.0, "R": 1.0, "N": -14.03, "D": 1.0, "C": 1.0, "E": 20.26, "Q": -6.54, "G": -7.49, "H": 1.0, "I": 1.0, "L": 1.0, "K": 1.0, "M": 1.0, "F": 13.34, "P": 1.0, "S": 1.0, "T": 1.0, "W": -14.03, "Y": 1.0, "V": 1.0},
    "W": {"A":  -14.03, "R": 1.0, "N": 13.34, "D": 1.0, "C": 1.0, "E": 1.0, "Q": 1.0, "G": -9.37, "H": 24.68, "I": 1.0, "L": 13.34, "K": 1.0, "M": 24.68, "F": 1.0, "P": 1.0, "S": 1.0, "T": -14.03, "W": 1.0, "Y": 1.0, "V": -7.49},
    "Y": {"A":  24.68, "R": -15.91, "N": 1.0, "D": 24.68, "C": 1.0, "E": -6.54, "Q": 1.0, "G": -7.49, "H": 13.34, "I": 1.0, "L": 1.0, "K": 1.0, "M": 44.94, "F": 1.0, "P": 13.34, "S": 1.0, "T": -7.49, "W": -9.37, "Y": 13.34, "V": 1.0},
    "V": {"A":  1.0, "R": 1.0, "N": 1.0, "D": -14.03, "C": 1.0, "E": 1.0, "Q": 1.0, "G": -7.49, "H": 1.0, "I": 1.0, "L": 1.0, "K": -1.88, "M": 1.0, "F": 1.0, "P": 20.26, "S": 1.0, "T": -7.49, "W": 1.0, "Y": -6.54, "V": 1.0},
}

# Tabela de uso de codons para organismos comuns.
# Valores = frequencia relativa do codon mais usado para cada aminoacido.
# Fonte: Kazusa Codon Usage Database (https://www.kazusa.or.jp/codon/)
_CODON_TABLE: Final[dict[str, str]] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Frequencia relativa de sinonimos (RSCU) para Canis lupus familiaris.
# Codons mais frequentes tem valor proximo de 1.0; raros proximo de 0.
# Fonte: Kazusa Codon Usage Database / GenBank RefSeq canino.
_RSCU_CANIS: Final[dict[str, float]] = {
    "TTT": 0.45, "TTC": 0.55, "TTA": 0.07, "TTG": 0.13,
    "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40,
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16, "ATG": 1.00,
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    "TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06,
    "CCT": 0.29, "CCC": 0.32, "CCA": 0.28, "CCG": 0.11,
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "TAT": 0.43, "TAC": 0.57, "TAA": 0.28, "TAG": 0.20,
    "CAT": 0.41, "CAC": 0.59, "CAA": 0.27, "CAG": 0.73,
    "AAT": 0.46, "AAC": 0.54, "AAA": 0.42, "AAG": 0.58,
    "GAT": 0.46, "GAC": 0.54, "GAA": 0.42, "GAG": 0.58,
    "TGT": 0.45, "TGC": 0.55, "TGA": 0.52, "TGG": 1.00,
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21,
    "AGT": 0.15, "AGC": 0.24, "AGA": 0.21, "AGG": 0.20,
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}

# Frequencia relativa de sinonimos para E. coli K-12 (expressao procariotica).
_RSCU_ECOLI: Final[dict[str, float]] = {
    "TTT": 0.58, "TTC": 0.42, "TTA": 0.14, "TTG": 0.13,
    "CTT": 0.12, "CTC": 0.10, "CTA": 0.04, "CTG": 0.47,
    "ATT": 0.49, "ATC": 0.39, "ATA": 0.11, "ATG": 1.00,
    "GTT": 0.28, "GTC": 0.20, "GTA": 0.17, "GTG": 0.35,
    "TCT": 0.17, "TCC": 0.15, "TCA": 0.14, "TCG": 0.14,
    "CCT": 0.18, "CCC": 0.13, "CCA": 0.20, "CCG": 0.49,
    "ACT": 0.19, "ACC": 0.40, "ACA": 0.17, "ACG": 0.25,
    "GCT": 0.18, "GCC": 0.27, "GCA": 0.21, "GCG": 0.34,
    "TAT": 0.59, "TAC": 0.41, "TAA": 0.61, "TAG": 0.09,
    "CAT": 0.57, "CAC": 0.43, "CAA": 0.34, "CAG": 0.66,
    "AAT": 0.49, "AAC": 0.51, "AAA": 0.74, "AAG": 0.26,
    "GAT": 0.63, "GAC": 0.37, "GAA": 0.68, "GAG": 0.32,
    "TGT": 0.46, "TGC": 0.54, "TGA": 0.30, "TGG": 1.00,
    "CGT": 0.36, "CGC": 0.36, "CGA": 0.07, "CGG": 0.11,
    "AGT": 0.16, "AGC": 0.25, "AGA": 0.07, "AGG": 0.04,
    "GGT": 0.35, "GGC": 0.37, "GGA": 0.13, "GGG": 0.15,
}

# Mapa de organismo para tabela RSCU correspondente
_RSCU_TABLES: Final[dict[str, dict[str, float]]] = {
    "canis": _RSCU_CANIS,
    "canis_lupus_familiaris": _RSCU_CANIS,
    "dog": _RSCU_CANIS,
    "ecoli": _RSCU_ECOLI,
    "e_coli": _RSCU_ECOLI,
    "escherichia_coli": _RSCU_ECOLI,
}


# ===========================================================================
# Funcoes de validacao -- proteinas
# ===========================================================================

def calculate_molecular_weight(sequence: str) -> float:
    """Calcula o peso molecular de uma sequencia proteica em Daltons.

    O peso molecular de um peptideo e a soma dos pesos dos aminoacidos
    individuais menos (n-1) moleculas de agua perdidas na formacao das
    ligacoes peptidicas.

    MW = sum(MW_residuos) - (n-1) * 18.015

    Args:
        sequence: sequencia de aminoacidos em letra unica (ex: 'MKLLVV')

    Returns:
        Peso molecular em Daltons (Da)

    Raises:
        ValueError: se a sequencia contem caracteres nao reconhecidos
    """
    sequence = sequence.upper().strip()
    if not sequence:
        return 0.0

    unknown = set(sequence) - set(AMINO_ACID_MW)
    if unknown:
        raise ValueError(
            f"Aminoacido(s) nao reconhecido(s): {', '.join(sorted(unknown))}"
        )

    total = sum(AMINO_ACID_MW[aa] for aa in sequence)
    # Subtrair agua de condensacao para cada ligacao peptidica
    water_loss = (len(sequence) - 1) * _WATER_MW
    return round(total - water_loss, 2)


def calculate_isoelectric_point(sequence: str) -> float:
    """Calcula o ponto isoeletrico (pI) pelo metodo de biseccao.

    O pI e o pH onde a carga liquida da proteina e zero. Usamos o
    metodo iterativo de biseccao sobre a funcao de carga liquida,
    calculada pela equacao de Henderson-Hasselbalch:

    Para grupos acidos (C-term, D, E, C, Y):
        carga = -1 / (1 + 10^(pKa - pH))

    Para grupos basicos (N-term, K, R, H):
        carga = +1 / (1 + 10^(pH - pKa))

    A precisao e de 0.01 unidades de pH.

    Args:
        sequence: sequencia de aminoacidos em letra unica

    Returns:
        Ponto isoeletrico (pH onde carga liquida = 0)
    """
    sequence = sequence.upper().strip()
    if not sequence:
        return 0.0

    # Contar residuos ionizaveis
    counts: dict[str, int] = {}
    for aa in sequence:
        if aa in _PKA_VALUES:
            counts[aa] = counts.get(aa, 0) + 1

    def _net_charge(ph: float) -> float:
        """Calcula carga liquida a um dado pH usando Henderson-Hasselbalch."""
        charge = 0.0

        # N-terminal (basico): +1 quando protonado
        charge += 1.0 / (1.0 + 10.0 ** (ph - _PKA_VALUES["N_TERM"]))

        # C-terminal (acido): -1 quando desprotonado
        charge -= 1.0 / (1.0 + 10.0 ** (_PKA_VALUES["C_TERM"] - ph))

        # Cadeias laterais basicas: H, K, R (contribuem carga positiva)
        for aa in ("H", "K", "R"):
            n = counts.get(aa, 0)
            if n > 0:
                charge += n / (1.0 + 10.0 ** (ph - _PKA_VALUES[aa]))

        # Cadeias laterais acidas: D, E, C, Y (contribuem carga negativa)
        for aa in ("D", "E", "C", "Y"):
            n = counts.get(aa, 0)
            if n > 0:
                charge -= n / (1.0 + 10.0 ** (_PKA_VALUES[aa] - ph))

        return charge

    # Biseccao entre pH 0 e pH 14
    ph_low = 0.0
    ph_high = 14.0
    while (ph_high - ph_low) > 0.001:
        ph_mid = (ph_low + ph_high) / 2.0
        if _net_charge(ph_mid) > 0.0:
            ph_low = ph_mid
        else:
            ph_high = ph_mid

    return round((ph_low + ph_high) / 2.0, 2)


def calculate_instability_index(sequence: str) -> float:
    """Calcula o indice de instabilidade pelo metodo de Guruprasad (1990).

    O indice quantifica a estabilidade in vivo de uma proteina com base
    na composicao de dipeptideos. Proteinas com II > 40 sao classificadas
    como instaveis.

    Formula:
        II = (10 / L) * sum(DIWV[xi][xi+1]) para i = 0..L-2

    onde L e o comprimento da sequencia e DIWV e a tabela de pesos de
    instabilidade de dipeptideos.

    Args:
        sequence: sequencia de aminoacidos em letra unica

    Returns:
        Indice de instabilidade (< 40 = estavel, >= 40 = instavel)
    """
    sequence = sequence.upper().strip()
    n = len(sequence)
    if n < 2:
        return 0.0

    total = 0.0
    for i in range(n - 1):
        aa1 = sequence[i]
        aa2 = sequence[i + 1]
        if aa1 in _DIWV and aa2 in _DIWV[aa1]:
            total += _DIWV[aa1][aa2]

    return round((10.0 / n) * total, 2)


def calculate_gravy(sequence: str) -> float:
    """Calcula o GRAVY (Grand Average of Hydropathy) pela escala de Kyte-Doolittle.

    O GRAVY e a media aritmetica da hidrofobicidade de todos os residuos.
    Valores positivos indicam proteina globular hidrofobica; negativos indicam
    proteina soluvel/hidrofilica.

    Para vacinas recombinantes, GRAVY levemente negativo e desejavel
    (boa solubilidade para formulacao).

    Args:
        sequence: sequencia de aminoacidos em letra unica

    Returns:
        Valor GRAVY (media de hidrofobicidade Kyte-Doolittle)
    """
    sequence = sequence.upper().strip()
    if not sequence:
        return 0.0

    values = []
    for aa in sequence:
        if aa in KYTE_DOOLITTLE:
            values.append(KYTE_DOOLITTLE[aa])

    if not values:
        return 0.0

    return round(sum(values) / len(values), 4)


def calculate_aromaticity(sequence: str) -> float:
    """Calcula a aromaticidade relativa da proteina.

    Frequencia relativa de aminoacidos aromaticos (F, W, Y) na sequencia.
    Importante para estabilidade estrutural e interacoes pi-stacking.

    Args:
        sequence: sequencia de aminoacidos em letra unica

    Returns:
        Fracao de aminoacidos aromaticos (0.0 a 1.0)
    """
    sequence = sequence.upper().strip()
    if not sequence:
        return 0.0

    aromatic_count = sum(1 for aa in sequence if aa in ("F", "W", "Y"))
    return round(aromatic_count / len(sequence), 4)


# ===========================================================================
# Funcoes de validacao -- DNA/mRNA
# ===========================================================================

def calculate_gc_content(dna: str) -> float:
    """Calcula o conteudo GC de uma sequencia de DNA ou mRNA.

    O conteudo GC ideal para expressao em mamiferos e entre 40-60%.
    Valores muito altos ou muito baixos prejudicam a estabilidade do
    mRNA e a eficiencia de traducao.

    Para mRNA terapeutico, GC ~ 55-65% e considerado otimo apos
    otimizacao de codons (substituicao de U por pseudouridina nao
    altera a contagem GC).

    Args:
        dna: sequencia de nucleotideos (aceita T ou U)

    Returns:
        Fracao GC (0.0 a 1.0)
    """
    dna = dna.upper().strip()
    if not dna:
        return 0.0

    gc_count = sum(1 for nt in dna if nt in ("G", "C"))
    return round(gc_count / len(dna), 4)


def calculate_cai(dna: str, organism: str = "canis") -> float:
    """Calcula o Codon Adaptation Index (CAI) para um dado organismo.

    O CAI mede quao bem os codons de uma sequencia correspondem ao
    uso preferencial de codons do organismo hospedeiro. Varia de 0 a 1;
    valores > 0.8 indicam boa adaptacao.

    O calculo usa a media geometrica das frequencias relativas (RSCU)
    de cada codon, normalizada pelo codon mais frequente de cada
    familia de sinonimos.

    Formula:
        CAI = exp((1/L) * sum(ln(w_i)))

    onde w_i e o peso relativo do codon i (RSCU_i / max_RSCU para
    aquele aminoacido) e L e o numero de codons.

    Args:
        dna: sequencia codificante (CDS), comprimento multiplo de 3
        organism: organismo alvo ('canis', 'ecoli', etc.)

    Returns:
        CAI entre 0.0 e 1.0

    Raises:
        ValueError: se o organismo nao e suportado ou sequencia invalida
    """
    dna = dna.upper().strip().replace("U", "T")
    organism_key = organism.lower().strip()

    if organism_key not in _RSCU_TABLES:
        supported = ", ".join(sorted(_RSCU_TABLES.keys()))
        raise ValueError(
            f"Organismo '{organism}' nao suportado. Opcoes: {supported}"
        )

    if len(dna) < 3:
        return 0.0

    rscu = _RSCU_TABLES[organism_key]

    # Calcular peso relativo (w) para cada codon:
    # w_i = RSCU_i / max(RSCU para codons sinonimos do mesmo aminoacido)
    # Primeiro, agrupar codons por aminoacido e encontrar o max RSCU
    aa_max_rscu: dict[str, float] = {}
    for codon, amino in _CODON_TABLE.items():
        if amino == "*":
            continue
        current_max = aa_max_rscu.get(amino, 0.0)
        codon_rscu = rscu.get(codon, 0.0)
        if codon_rscu > current_max:
            aa_max_rscu[amino] = codon_rscu

    # Calcular w para cada codon
    codon_weights: dict[str, float] = {}
    for codon, amino in _CODON_TABLE.items():
        if amino == "*":
            continue
        max_rscu = aa_max_rscu.get(amino, 1.0)
        if max_rscu > 0:
            codon_weights[codon] = rscu.get(codon, 0.0) / max_rscu
        else:
            codon_weights[codon] = 0.0

    # Calcular CAI como media geometrica dos pesos
    log_sum = 0.0
    codon_count = 0

    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        if len(codon) < 3:
            break
        if codon not in _CODON_TABLE:
            continue
        amino = _CODON_TABLE[codon]
        if amino == "*":
            continue
        w = codon_weights.get(codon, 0.0)
        if w > 0:
            log_sum += math.log(w)
            codon_count += 1

    if codon_count == 0:
        return 0.0

    cai = math.exp(log_sum / codon_count)
    return round(cai, 4)


def check_restriction_sites(
    dna: str,
    sites: dict[str, str] | None = None,
) -> list[str]:
    """Busca sitios de restricao indesejados na sequencia de DNA.

    Sitios de restricao dentro da CDS podem causar problemas durante
    a clonagem. Esta funcao verifica a presenca de enzimas comuns
    usadas em vetores de expressao.

    Se nenhum dicionario de sitios for fornecido, usa um conjunto
    padrao de enzimas comuns em clonagem molecular.

    Args:
        dna: sequencia de DNA a ser verificada
        sites: dicionario {nome_enzima: sequencia_de_reconhecimento}

    Returns:
        Lista de nomes de enzimas cujos sitios foram encontrados
    """
    # Sitios de restricao comuns em vetores de expressao (clonagem)
    default_sites: dict[str, str] = {
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "NdeI": "CATATG",
        "XhoI": "CTCGAG",
        "NcoI": "CCATGG",
        "BglII": "AGATCT",
        "SalI": "GTCGAC",
        "NotI": "GCGGCCGC",
        "XbaI": "TCTAGA",
        "SpeI": "ACTAGT",
        "KpnI": "GGTACC",
        "SacI": "GAGCTC",
        "PstI": "CTGCAG",
        "ApaI": "GGGCCC",
        "BsaI": "GGTCTC",
        "BbsI": "GAAGAC",
    }

    if sites is None:
        sites = default_sites

    dna = dna.upper().strip()
    found: list[str] = []

    for enzyme, recognition_seq in sites.items():
        recognition_seq = recognition_seq.upper()
        if recognition_seq in dna:
            found.append(enzyme)

    return sorted(found)


# ===========================================================================
# Validacao integrada
# ===========================================================================

def validate_construct(sequence: str) -> dict[str, float | str | bool]:
    """Executa todos os validadores proteicos e retorna um resumo.

    Esta funcao e o ponto de entrada principal para validacao rapida
    de qualquer construto vacinal. Retorna um dicionario com todos
    os parametros biofisicos e classificacoes de estabilidade.

    Criterios de classificacao (baseados na literatura):
    - MW: informativo (sem limite estrito)
    - pI: ideal entre 6-8 para solubilidade em pH fisiologico
    - II < 40: proteina estavel in vivo (Guruprasad, 1990)
    - GRAVY negativo: boa solubilidade aquosa
    - Aromaticidade: informativo

    Args:
        sequence: sequencia aminoacidica do construto

    Returns:
        Dicionario com todos os parametros calculados e classificacoes
    """
    sequence = sequence.upper().strip()

    mw = calculate_molecular_weight(sequence)
    pi = calculate_isoelectric_point(sequence)
    ii = calculate_instability_index(sequence)
    gravy = calculate_gravy(sequence)
    arom = calculate_aromaticity(sequence)

    # Classificacoes baseadas em criterios da literatura
    is_stable = ii < 40.0
    is_soluble = gravy < 0.0
    pi_physiological = 6.0 <= pi <= 8.0

    return {
        "length_aa": len(sequence),
        "molecular_weight_da": mw,
        "isoelectric_point": pi,
        "instability_index": ii,
        "is_stable": is_stable,
        "gravy": gravy,
        "is_soluble": is_soluble,
        "aromaticity": arom,
        "pi_physiological_range": pi_physiological,
        "classification": (
            "ESTAVEL e SOLUVEL"
            if is_stable and is_soluble
            else "ESTAVEL mas HIDROFOBICO"
            if is_stable
            else "INSTAVEL e SOLUVEL"
            if is_soluble
            else "INSTAVEL e HIDROFOBICO"
        ),
    }
