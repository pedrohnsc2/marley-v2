"""Avaliacao de sequencias geradas pelo modelo de difusao.

Para ASOs: usa funcoes termodinamicas do aso_math (dG, Tm, GC%)
Para epitopos: calcula propriedades fisicoquimicas (MW, hidrofobicidade, carga)

Filtra candidatos por criterios de qualidade e rankeia por score
multi-objetivo. Calcula novidade via distancia de edicao (Levenshtein)
em relacao as sequencias de treinamento.
"""

from __future__ import annotations

import re
from typing import Any, Final

from aso_math.thermo import compute_dg, compute_tm, gc_content


# ---------------------------------------------------------------------------
# Tabelas de propriedades de aminoacidos
# ---------------------------------------------------------------------------

# Pesos moleculares medios dos aminoacidos (Da) — residuo (sem agua)
# Fonte: NIST Standard Reference Database
AA_MW: Final[dict[str, float]] = {
    "A": 71.04, "R": 156.10, "N": 114.04, "D": 115.03, "C": 103.01,
    "E": 129.04, "Q": 128.06, "G": 57.02, "H": 137.06, "I": 113.08,
    "L": 113.08, "K": 128.09, "M": 131.04, "F": 147.07, "P": 97.05,
    "S": 87.03, "T": 101.05, "W": 186.08, "Y": 163.06, "V": 99.07,
}

# Escala de hidrofobicidade Kyte-Doolittle (1982)
# Valores positivos = hidrofobico, negativos = hidrofilico
KYTE_DOOLITTLE: Final[dict[str, float]] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Carga a pH 7.0 (simplificada)
# K, R, H = +1; D, E = -1; outros = 0
AA_CHARGE: Final[dict[str, float]] = {
    "A": 0, "R": 1, "N": 0, "D": -1, "C": 0,
    "E": -1, "Q": 0, "G": 0, "H": 0.1, "I": 0,
    "L": 0, "K": 1, "M": 0, "F": 0, "P": 0,
    "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
}


# ---------------------------------------------------------------------------
# Funcoes auxiliares
# ---------------------------------------------------------------------------

def levenshtein_distance(s1: str, s2: str) -> int:
    """Calcula distancia de edicao de Levenshtein entre duas strings.

    Implementacao classica com programacao dinamica O(n*m).
    Usada para medir novidade das sequencias geradas.

    Args:
        s1: Primeira string.
        s2: Segunda string.

    Returns:
        Numero minimo de edicoes (insercao, delecao, substituicao).
    """
    n, m = len(s1), len(s2)
    # Otimizacao: usar apenas duas linhas da matriz
    prev = list(range(m + 1))
    curr = [0] * (m + 1)

    for i in range(1, n + 1):
        curr[0] = i
        for j in range(1, m + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            curr[j] = min(
                prev[j] + 1,       # Delecao
                curr[j - 1] + 1,   # Insercao
                prev[j - 1] + cost, # Substituicao
            )
        prev, curr = curr, prev

    return prev[m]


def min_edit_distance(seq: str, reference_seqs: list[str]) -> int:
    """Calcula a menor distancia de edicao entre seq e qualquer referencia.

    Args:
        seq: Sequencia candidata.
        reference_seqs: Lista de sequencias de referencia (treinamento).

    Returns:
        Menor distancia de edicao encontrada.
    """
    if not reference_seqs:
        return len(seq)
    return min(levenshtein_distance(seq, ref) for ref in reference_seqs)


def has_homopolymer(seq: str, max_repeat: int = 4) -> bool:
    """Verifica se a sequencia contem homopolymers acima do limite.

    Homopolymers (ex: AAAA, GGGG) causam problemas na sintese de
    oligonucleotideos e podem gerar estruturas secundarias indesejaveis.

    Args:
        seq: Sequencia de nucleotideos.
        max_repeat: Numero maximo de repeticoes consecutivas permitidas.

    Returns:
        True se contem homopolymer acima do limite.
    """
    pattern = r"(.)\1{" + str(max_repeat) + r",}"
    return bool(re.search(pattern, seq))


# ---------------------------------------------------------------------------
# Avaliacao de ASOs
# ---------------------------------------------------------------------------

def evaluate_aso(
    sequence: str,
    training_seqs: list[str],
    gc_min: float = 0.25,
    gc_max: float = 0.65,
    homopolymer_max: int = 4,
) -> dict[str, Any]:
    """Avalia uma variante de ASO com metricas termodinamicas.

    Calcula dG (energia livre de hibridizacao), Tm (temperatura de melting),
    GC content, e verifica filtros de qualidade.

    Args:
        sequence: Sequencia do ASO candidato.
        training_seqs: Sequencias de treinamento (para novidade).
        gc_min: GC content minimo aceitavel.
        gc_max: GC content maximo aceitavel.
        homopolymer_max: Maximo de bases repetidas consecutivas.

    Returns:
        Dicionario com metricas e flags de filtragem.
    """
    seq = sequence.upper()
    gc = gc_content(seq)
    dg = compute_dg(seq)
    tm = compute_tm(seq)
    novelty = min_edit_distance(seq, training_seqs)

    # Filtros de qualidade
    gc_ok = gc_min <= gc <= gc_max
    homopolymer_ok = not has_homopolymer(seq, homopolymer_max)
    valid_bases = all(c in "ACGT" for c in seq)

    # Score multi-objetivo para ASO:
    # - dG mais negativo = melhor (ligacao mais forte)
    # - Tm entre 60-75C = ideal para ASO
    # - GC entre 0.25-0.65 = aceitavel
    # - Novidade > 0 = nao e copia direta
    dg_score = min(abs(dg) / 30.0, 1.0)  # Normaliza: -30 kcal/mol = 1.0
    tm_score = max(0.0, 1.0 - abs(tm - 67.5) / 20.0)  # Otimo em ~67.5C
    gc_score = 1.0 if gc_ok else 0.3
    novelty_score = min(novelty / 5.0, 1.0)  # 5+ edicoes = novidade maxima
    homopolymer_score = 1.0 if homopolymer_ok else 0.3

    # Score final ponderado
    composite_score = (
        0.35 * dg_score
        + 0.25 * tm_score
        + 0.15 * gc_score
        + 0.15 * novelty_score
        + 0.10 * homopolymer_score
    )

    passes_filter = gc_ok and homopolymer_ok and valid_bases and novelty > 0

    return {
        "sequence": seq,
        "length": len(seq),
        "dg": dg,
        "tm": tm,
        "gc": round(gc, 4),
        "edit_distance": novelty,
        "gc_ok": gc_ok,
        "homopolymer_ok": homopolymer_ok,
        "valid_bases": valid_bases,
        "passes_filter": passes_filter,
        "score": round(composite_score, 4),
        "dg_score": round(dg_score, 4),
        "tm_score": round(tm_score, 4),
        "gc_score": round(gc_score, 4),
        "novelty_score": round(novelty_score, 4),
    }


# ---------------------------------------------------------------------------
# Avaliacao de epitopos
# ---------------------------------------------------------------------------

def compute_mw(sequence: str) -> float:
    """Calcula peso molecular aproximado de um peptideo (Da).

    MW = soma dos residuos + 18.02 (agua adicionada na hidrolise).

    Args:
        sequence: Sequencia de aminoacidos.

    Returns:
        Peso molecular em Daltons.
    """
    mw = sum(AA_MW.get(aa, 110.0) for aa in sequence.upper())
    mw += 18.02  # Agua
    return round(mw, 2)


def compute_hydrophobicity(sequence: str) -> float:
    """Calcula hidrofobicidade media (Kyte-Doolittle) do peptideo.

    Args:
        sequence: Sequencia de aminoacidos.

    Returns:
        Score medio de hidrofobicidade.
    """
    if not sequence:
        return 0.0
    scores = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in sequence.upper()]
    return round(sum(scores) / len(scores), 4)


def compute_charge(sequence: str) -> float:
    """Calcula carga liquida aproximada a pH 7.0.

    Args:
        sequence: Sequencia de aminoacidos.

    Returns:
        Carga liquida do peptideo.
    """
    return round(sum(AA_CHARGE.get(aa, 0.0) for aa in sequence.upper()), 1)


def evaluate_epitope(
    sequence: str,
    training_seqs: list[str],
    min_hydro: float = -2.0,
    max_hydro: float = 2.0,
) -> dict[str, Any]:
    """Avalia uma variante de epitopo com propriedades fisicoquimicas.

    Calcula MW, hidrofobicidade, carga e novidade em relacao aos epitopos
    originais. Epitopos ideais para vacinas tendem a ser anfipáticos
    (hidrofobicidade moderada) com carga levemente positiva.

    Args:
        sequence: Sequencia do epitopo candidato.
        training_seqs: Sequencias de treinamento (para novidade).
        min_hydro: Hidrofobicidade minima aceitavel.
        max_hydro: Hidrofobicidade maxima aceitavel.

    Returns:
        Dicionario com metricas e flags de filtragem.
    """
    seq = sequence.upper()
    mw = compute_mw(seq)
    hydro = compute_hydrophobicity(seq)
    charge = compute_charge(seq)
    novelty = min_edit_distance(seq, training_seqs)

    # Filtros de qualidade
    valid_aas = all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in seq)
    hydro_ok = min_hydro <= hydro <= max_hydro

    # Score multi-objetivo para epitopos vacinais:
    # - MW moderado (900-1200 Da para 9-mers) = boa apresentacao no MHC
    # - Hidrofobicidade moderada = anfipático (bom para binding MHC)
    # - Carga levemente positiva = favorece interacao com celula
    # - Novidade > 0 = nao e copia de epitopo existente
    target_mw = 1050.0  # MW tipico de 9-mer
    mw_score = max(0.0, 1.0 - abs(mw - target_mw) / 500.0)
    hydro_score = max(0.0, 1.0 - abs(hydro) / 3.0)  # Otimo perto de 0
    charge_score = max(0.0, 1.0 - abs(charge - 0.5) / 4.0)  # Levemente positivo
    novelty_score = min(novelty / 4.0, 1.0)

    composite_score = (
        0.25 * mw_score
        + 0.25 * hydro_score
        + 0.20 * charge_score
        + 0.30 * novelty_score
    )

    passes_filter = valid_aas and hydro_ok and novelty > 0

    return {
        "sequence": seq,
        "length": len(seq),
        "mw": mw,
        "hydrophobicity": hydro,
        "charge": charge,
        "edit_distance": novelty,
        "valid_aas": valid_aas,
        "hydro_ok": hydro_ok,
        "passes_filter": passes_filter,
        "score": round(composite_score, 4),
        "mw_score": round(mw_score, 4),
        "hydro_score": round(hydro_score, 4),
        "charge_score": round(charge_score, 4),
        "novelty_score": round(novelty_score, 4),
    }


# ---------------------------------------------------------------------------
# Avaliacao em batch
# ---------------------------------------------------------------------------

def evaluate_batch(
    sequences: list[str],
    training_seqs: list[str],
    mode: str = "aso",
    **kwargs: Any,
) -> list[dict[str, Any]]:
    """Avalia um batch de sequencias geradas e retorna ranking.

    Args:
        sequences: Lista de sequencias geradas.
        training_seqs: Lista de sequencias de treinamento (para novidade).
        mode: "aso" para nucleotideos, "epitope" para aminoacidos.
        **kwargs: Parametros extras passados ao avaliador especifico.

    Returns:
        Lista de avaliacoes, ordenada por score decrescente.
    """
    results = []

    for seq in sequences:
        if mode == "aso":
            result = evaluate_aso(seq, training_seqs, **kwargs)
        elif mode == "epitope":
            result = evaluate_epitope(seq, training_seqs, **kwargs)
        else:
            raise ValueError(f"Modo invalido: '{mode}'. Use 'aso' ou 'epitope'.")
        results.append(result)

    # Ordena por score decrescente
    results.sort(key=lambda r: r["score"], reverse=True)

    return results


def filter_and_rank(
    evaluations: list[dict[str, Any]],
    top_k: int = 10,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Filtra candidatos validos e retorna top-K.

    Args:
        evaluations: Lista de avaliacoes (ja ordenadas por score).
        top_k: Numero de candidatos a retornar.

    Returns:
        Tupla (top_k_validos, todos_validos).
    """
    valid = [e for e in evaluations if e["passes_filter"]]
    return valid[:top_k], valid
