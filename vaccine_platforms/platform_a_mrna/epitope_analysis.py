"""A1 -- Analise de redundancia e cobertura de epitopos.

Compara os 11 epitopos canonicos par-a-par usando:
    - Distancia de Hamming normalizada (similaridade de sequencia)
    - Correlacao de IC50 (ambos ligantes fortes?)
    - Mesmo gene de origem? Mesmo alelo DLA?

Identifica pares redundantes (>70% similaridade + mesmo alelo) e testa
variantes do construto com menos epitopos, removendo os mais redundantes
primeiro. Para cada variante, recalcula parametros biofisicos via
vaccine_platforms.shared.validators.

Os 11 epitopos sao IMUTAVEIS -- esta analise recomenda, nunca altera.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Final

from vaccine_platforms.shared.epitopes import (
    ADJUVANT_SEQUENCE,
    ADJUVANT_EPITOPE_LINKER,
    EPITOPES,
    LINKERS,
    SIGNAL_PEPTIDE_TPA,
    Epitope,
)
from vaccine_platforms.shared.validators import (
    calculate_instability_index,
    calculate_molecular_weight,
)

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

# Limiar para classificar pares como redundantes
SIMILARITY_THRESHOLD: Final[float] = 0.70
IC50_STRONG_BINDER: Final[float] = 100.0  # nM -- abaixo e considerado forte

# Variantes de construto a testar (numero de epitopos)
VARIANT_SIZES: Final[tuple[int, ...]] = (11, 9, 7)

# Limiar minimo de cobertura teorica para aceitar reducao
COVERAGE_THRESHOLD: Final[float] = 0.95


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PairwiseComparison:
    """Resultado da comparacao entre dois epitopos."""

    epitope_a: str
    epitope_b: str
    hamming_similarity: float
    same_gene: bool
    same_allele: bool
    both_strong_binders: bool
    is_redundant: bool


@dataclass(frozen=True, slots=True)
class ConstructVariant:
    """Variante do construto com N epitopos selecionados."""

    epitope_count: int
    epitopes_included: tuple[str, ...]
    epitopes_removed: tuple[str, ...]
    protein_length_aa: int
    molecular_weight_da: float
    instability_index: float
    is_stable: bool
    antigenicity_proxy: float
    allele_coverage: dict[str, int]
    gene_coverage: dict[str, int]
    theoretical_coverage_pct: float


@dataclass
class RedundancyAnalysis:
    """Resultado completo da analise de redundancia."""

    total_epitopes: int
    unique_genes: int
    unique_alleles: int
    pairwise_comparisons: list[PairwiseComparison]
    redundant_pairs: list[PairwiseComparison]
    construct_variants: list[ConstructVariant]
    minimum_epitopes_95pct: int
    recommendation: str


# ---------------------------------------------------------------------------
# Funcoes de calculo
# ---------------------------------------------------------------------------


def _hamming_similarity(seq_a: str, seq_b: str) -> float:
    """Calcula similaridade por distancia de Hamming normalizada.

    Para sequencias de comprimentos diferentes, preenche a menor
    com gaps e conta como mismatch. A similaridade e 1.0 - (hamming / max_len).

    Args:
        seq_a: primeira sequencia peptidica
        seq_b: segunda sequencia peptidica

    Returns:
        Similaridade entre 0.0 (totalmente diferente) e 1.0 (identico)
    """
    max_len = max(len(seq_a), len(seq_b))
    if max_len == 0:
        return 1.0

    # Preencher sequencia menor com gaps para alinhar
    a = seq_a.ljust(max_len, "-")
    b = seq_b.ljust(max_len, "-")

    matches = sum(1 for ca, cb in zip(a, b) if ca == cb)
    return round(matches / max_len, 4)


def _antigenicity_proxy(epitopes: list[Epitope]) -> float:
    """Proxy de antigenicidade baseado na composicao de aminoacidos.

    Modelo simplificado inspirado no VaxiJen: pontua a frequencia de
    aminoacidos associados a antigenicidade em proteinas parasitarias
    (residuos aromaticos, carregados e polares especificos).

    A pontuacao e normalizada para o intervalo ~0.4-0.7 (onde > 0.4
    seria considerado antigenico pelo VaxiJen com threshold default).

    NOTA: Este e apenas um proxy computacional. A validacao real requer
    o servidor VaxiJen v2.0 (http://www.ddg-pharmfac.net/vaxijen/).

    Args:
        epitopes: lista de epitopos incluidos na variante

    Returns:
        Pontuacao proxy de antigenicidade (escala ~0-1)
    """
    # Pesos derivados da analise de auto/cross-covariance (ACC) em
    # proteinas antigenicas de parasitas vs nao-antigenicas.
    # Aminoacidos com peso positivo aparecem mais em antigenos validados.
    antigenicity_weights: dict[str, float] = {
        "A": 0.02, "R": 0.15, "N": 0.08, "D": -0.05,
        "C": 0.10, "E": -0.03, "Q": 0.04, "G": -0.01,
        "H": 0.06, "I": 0.12, "L": 0.10, "K": 0.13,
        "M": 0.14, "F": 0.16, "P": -0.02, "S": 0.03,
        "T": 0.05, "W": 0.18, "Y": 0.15, "V": 0.11,
    }

    # Concatenar todos os peptideos da variante
    all_residues = "".join(ep.peptide for ep in epitopes)
    if not all_residues:
        return 0.0

    # Calcular pontuacao media ponderada
    score = sum(
        antigenicity_weights.get(aa, 0.0) for aa in all_residues
    ) / len(all_residues)

    # Normalizar para faixa VaxiJen-like (~0.3-0.8)
    # Os pesos brutos dao ~0.05-0.15; escalar para faixa biologica
    normalized = 0.4 + score * 2.5
    return round(max(0.0, min(1.0, normalized)), 4)


def _build_variant_sequence(epitopes: list[Epitope]) -> str:
    """Monta sequencia proteica de uma variante do construto.

    Usa a mesma arquitetura do construto original:
        [tPA]-[L7/L12]-EAAAK-[Epi1]-AAY-...-AAY-[EpiN]

    Args:
        epitopes: epitopos a incluir na variante

    Returns:
        Sequencia aminoacidica completa
    """
    linker = LINKERS["CTL"]
    parts = [
        SIGNAL_PEPTIDE_TPA,
        ADJUVANT_SEQUENCE,
        ADJUVANT_EPITOPE_LINKER,
        linker.join(ep.peptide for ep in epitopes),
    ]
    return "".join(parts)


def _allele_coverage(epitopes: list[Epitope]) -> dict[str, int]:
    """Conta epitopos por alelo DLA.

    Args:
        epitopes: lista de epitopos

    Returns:
        Dicionario {alelo: contagem}
    """
    coverage: dict[str, int] = {}
    for ep in epitopes:
        coverage[ep.allele] = coverage.get(ep.allele, 0) + 1
    return coverage


def _gene_coverage(epitopes: list[Epitope]) -> dict[str, int]:
    """Conta epitopos por gene de origem.

    Args:
        epitopes: lista de epitopos

    Returns:
        Dicionario {gene_name: contagem}
    """
    coverage: dict[str, int] = {}
    for ep in epitopes:
        coverage[ep.gene_name] = coverage.get(ep.gene_name, 0) + 1
    return coverage


def _theoretical_coverage(
    variant_epitopes: list[Epitope],
    all_epitopes: tuple[Epitope, ...],
) -> float:
    """Calcula cobertura teorica da variante vs construto completo.

    Considera tres dimensoes:
    1. Cobertura alelica -- fracao de alelos DLA representados
    2. Cobertura genica -- fracao de genes de L. infantum representados
    3. Cobertura de afinidade -- soma de 1/IC50 normalizada

    A cobertura final e a media ponderada das tres dimensoes.

    Args:
        variant_epitopes: epitopos na variante reduzida
        all_epitopes: todos os 11 epitopos canonicos

    Returns:
        Fracao de cobertura teorica (0.0 a 1.0)
    """
    # Dimensao 1: alelos DLA
    all_alleles = set(ep.allele for ep in all_epitopes)
    var_alleles = set(ep.allele for ep in variant_epitopes)
    allele_cov = len(var_alleles) / len(all_alleles) if all_alleles else 1.0

    # Dimensao 2: genes-alvo
    all_genes = set(ep.gene_name for ep in all_epitopes)
    var_genes = set(ep.gene_name for ep in variant_epitopes)
    gene_cov = len(var_genes) / len(all_genes) if all_genes else 1.0

    # Dimensao 3: soma de afinidade (1/IC50)
    # Quanto menor o IC50, maior a afinidade -- usamos 1/IC50 como peso
    all_affinity = sum(1.0 / ep.ic50 for ep in all_epitopes)
    var_affinity = sum(1.0 / ep.ic50 for ep in variant_epitopes)
    affinity_cov = var_affinity / all_affinity if all_affinity > 0 else 1.0

    # Media ponderada: afinidade tem peso 40%, alelo 35%, gene 25%
    coverage = 0.40 * affinity_cov + 0.35 * allele_cov + 0.25 * gene_cov
    return round(min(1.0, coverage), 4)


def _rank_epitopes_by_redundancy(
    epitopes: tuple[Epitope, ...],
    redundant_pairs: list[PairwiseComparison],
) -> list[Epitope]:
    """Ordena epitopos do mais redundante ao menos redundante.

    O "score de redundancia" de cada epitopo e o numero de vezes que
    aparece em pares redundantes, ponderado pela similaridade do par
    e penalizado se o epitopo for um binder forte (baixo IC50).

    Epitopos com maior score sao candidatos a remocao.

    Args:
        epitopes: todos os 11 epitopos
        redundant_pairs: pares classificados como redundantes

    Returns:
        Lista de epitopos ordenada do mais ao menos redundante
    """
    # Acumular score de redundancia por epitopo
    scores: dict[str, float] = {ep.peptide: 0.0 for ep in epitopes}

    for pair in redundant_pairs:
        # Contribuicao proporcional a similaridade do par
        scores[pair.epitope_a] += pair.hamming_similarity
        scores[pair.epitope_b] += pair.hamming_similarity

    # Penalizar epitopos com alta afinidade (nao queremos remove-los)
    # IC50 baixo = binder forte = penalidade maior (subtrair do score)
    ic50_map = {ep.peptide: ep.ic50 for ep in epitopes}
    max_ic50 = max(ep.ic50 for ep in epitopes)

    for peptide in scores:
        # Normalizar IC50: binder forte (IC50=10) -> penalidade ~0.9
        affinity_penalty = 1.0 - (ic50_map[peptide] / max_ic50)
        scores[peptide] -= affinity_penalty * 0.5

    # Ordenar do mais redundante (maior score) ao menos
    sorted_peptides = sorted(scores, key=lambda p: scores[p], reverse=True)

    # Mapear de volta para objetos Epitope
    ep_map = {ep.peptide: ep for ep in epitopes}
    return [ep_map[p] for p in sorted_peptides]


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def run_redundancy_analysis() -> RedundancyAnalysis:
    """Executa a analise completa de redundancia dos 11 epitopos.

    Etapas:
        1. Compara todos os pares (55 combinacoes de 11 choose 2)
        2. Identifica pares redundantes (>70% similaridade + mesmo alelo)
        3. Testa variantes com 11, 9, 7 epitopos
        4. Para cada variante: recalcula propriedades biofisicas
        5. Determina numero minimo mantendo >=95% cobertura teorica

    Returns:
        RedundancyAnalysis com resultados completos
    """
    epitopes = EPITOPES
    all_comparisons: list[PairwiseComparison] = []

    # --- Passo 1: Comparacoes par-a-par ---
    for ep_a, ep_b in combinations(epitopes, 2):
        similarity = _hamming_similarity(ep_a.peptide, ep_b.peptide)
        same_gene = ep_a.gene_id == ep_b.gene_id
        same_allele = ep_a.allele == ep_b.allele
        both_strong = (
            ep_a.ic50 < IC50_STRONG_BINDER and ep_b.ic50 < IC50_STRONG_BINDER
        )

        # Criterio de redundancia: alta similaridade + mesmo alelo
        is_redundant = similarity >= SIMILARITY_THRESHOLD and same_allele

        comp = PairwiseComparison(
            epitope_a=ep_a.peptide,
            epitope_b=ep_b.peptide,
            hamming_similarity=similarity,
            same_gene=same_gene,
            same_allele=same_allele,
            both_strong_binders=both_strong,
            is_redundant=is_redundant,
        )
        all_comparisons.append(comp)

    redundant_pairs = [c for c in all_comparisons if c.is_redundant]

    # --- Passo 2: Ranquear epitopos por redundancia ---
    ranked = _rank_epitopes_by_redundancy(epitopes, redundant_pairs)

    # --- Passo 3: Testar variantes ---
    variants: list[ConstructVariant] = []

    for size in VARIANT_SIZES:
        if size >= len(epitopes):
            # Construto completo -- incluir todos os 11
            included = list(epitopes)
            removed: list[Epitope] = []
        else:
            # Remover os epitopos mais redundantes ate atingir o tamanho alvo
            to_remove = len(epitopes) - size
            removed = ranked[:to_remove]
            removed_peptides = set(ep.peptide for ep in removed)
            included = [ep for ep in epitopes if ep.peptide not in removed_peptides]

        sequence = _build_variant_sequence(included)
        mw = calculate_molecular_weight(sequence)
        ii = calculate_instability_index(sequence)
        ag_proxy = _antigenicity_proxy(included)
        alleles = _allele_coverage(included)
        genes = _gene_coverage(included)
        coverage = _theoretical_coverage(included, epitopes)

        variant = ConstructVariant(
            epitope_count=len(included),
            epitopes_included=tuple(ep.peptide for ep in included),
            epitopes_removed=tuple(ep.peptide for ep in removed),
            protein_length_aa=len(sequence),
            molecular_weight_da=mw,
            instability_index=ii,
            is_stable=ii < 40.0,
            antigenicity_proxy=ag_proxy,
            allele_coverage=alleles,
            gene_coverage=genes,
            theoretical_coverage_pct=round(coverage * 100, 2),
        )
        variants.append(variant)

    # --- Passo 4: Determinar minimo com >=95% cobertura ---
    # Testar tamanhos decrescentes ate a cobertura cair abaixo de 95%
    min_epitopes = len(epitopes)
    for variant in sorted(variants, key=lambda v: v.epitope_count):
        if variant.theoretical_coverage_pct >= COVERAGE_THRESHOLD * 100:
            min_epitopes = variant.epitope_count
            break

    # --- Passo 5: Gerar recomendacao ---
    unique_genes = len(set(ep.gene_name for ep in epitopes))
    unique_alleles = len(set(ep.allele for ep in epitopes))

    if len(redundant_pairs) == 0:
        recommendation = (
            "Nenhum par redundante encontrado. Todos os 11 epitopos sao "
            "suficientemente distintos em sequencia e especificidade alelica. "
            "Recomenda-se manter o construto completo com 11 epitopos."
        )
    elif min_epitopes == len(epitopes):
        recommendation = (
            f"Encontrados {len(redundant_pairs)} par(es) com alta similaridade, "
            f"mas a cobertura teorica cai abaixo de 95% ao remover qualquer epitopo. "
            f"Recomenda-se manter todos os 11 epitopos no construto."
        )
    else:
        recommendation = (
            f"Encontrados {len(redundant_pairs)} par(es) redundantes. "
            f"O construto pode ser reduzido para {min_epitopes} epitopos "
            f"mantendo >= 95% da cobertura teorica. Contudo, os 11 epitopos "
            f"originais oferecem redundancia intencional contra diversidade "
            f"genetica do parasita e da populacao canina."
        )

    return RedundancyAnalysis(
        total_epitopes=len(epitopes),
        unique_genes=unique_genes,
        unique_alleles=unique_alleles,
        pairwise_comparisons=all_comparisons,
        redundant_pairs=redundant_pairs,
        construct_variants=variants,
        minimum_epitopes_95pct=min_epitopes,
        recommendation=recommendation,
    )


def to_dict(analysis: RedundancyAnalysis) -> dict:
    """Serializa a analise de redundancia para JSON.

    Args:
        analysis: resultado da analise de redundancia

    Returns:
        Dicionario serializavel para JSON
    """
    return {
        "total_epitopes": analysis.total_epitopes,
        "unique_source_genes": analysis.unique_genes,
        "unique_alleles": analysis.unique_alleles,
        "pairwise_comparisons_count": len(analysis.pairwise_comparisons),
        "redundant_pairs_count": len(analysis.redundant_pairs),
        "redundant_pairs": [
            {
                "epitope_a": p.epitope_a,
                "epitope_b": p.epitope_b,
                "hamming_similarity": p.hamming_similarity,
                "same_gene": p.same_gene,
                "same_allele": p.same_allele,
                "both_strong_binders": p.both_strong_binders,
            }
            for p in analysis.redundant_pairs
        ],
        "top_similar_pairs": sorted(
            [
                {
                    "epitope_a": p.epitope_a,
                    "epitope_b": p.epitope_b,
                    "hamming_similarity": p.hamming_similarity,
                    "same_gene": p.same_gene,
                    "same_allele": p.same_allele,
                }
                for p in analysis.pairwise_comparisons
            ],
            key=lambda x: x["hamming_similarity"],
            reverse=True,
        )[:10],
        "construct_variants": [
            {
                "epitope_count": v.epitope_count,
                "epitopes_included": list(v.epitopes_included),
                "epitopes_removed": list(v.epitopes_removed),
                "protein_length_aa": v.protein_length_aa,
                "molecular_weight_da": v.molecular_weight_da,
                "instability_index": v.instability_index,
                "is_stable": v.is_stable,
                "antigenicity_proxy": v.antigenicity_proxy,
                "allele_coverage": v.allele_coverage,
                "gene_coverage": v.gene_coverage,
                "theoretical_coverage_pct": v.theoretical_coverage_pct,
            }
            for v in analysis.construct_variants
        ],
        "minimum_epitopes_95pct_coverage": analysis.minimum_epitopes_95pct,
        "recommendation": analysis.recommendation,
    }
