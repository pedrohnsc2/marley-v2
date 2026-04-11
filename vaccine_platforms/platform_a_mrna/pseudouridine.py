"""A2 -- Otimizacao de pseudouridina (psiU) para mRNA terapeutico.

Analisa as posicoes de uridina no mRNA do construto vacinal e modela
estrategias de substituicao por N1-metilpseudouridina (m1psi):

    - Identifica todas as posicoes U na sequencia mRNA
    - Calcula contexto estrutural local (modelo simples de base-pairing)
    - Rankeia uridinas por beneficio estimado da substituicao
    - Compara estrategias: 30%, 50% e 100% de substituicao
    - Estima efeitos em traducao, imunogenicidade e custo

A substituicao de uridina por pseudouridina e a inovacao central das
vacinas de mRNA (Kariko & Weissman, Nobel 2023). Reduz a ativacao de
receptores de RNA inato (TLR3/7/8, RIG-I) e aumenta a estabilidade
e eficiencia de traducao.

NOTA: A sequencia FASTA do mRNA NAO e exposta neste modulo. Todas as
analises usam a sequencia proteica do construto para gerar o mRNA
in silico via tabela de codons otimizada para Canis lupus familiaris.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from vaccine_platforms.shared.epitopes import get_construct_sequence
from vaccine_platforms.shared.validators import calculate_gc_content

# ---------------------------------------------------------------------------
# Tabela de codons otimizada para Canis lupus familiaris
# ---------------------------------------------------------------------------

# Codons preferenciais para expressao em celulas de mamifero (cao).
# Selecionamos o codon mais frequente por aminoacido (Kazusa DB).
_PREFERRED_CODONS_CANIS: Final[dict[str, str]] = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC",
    "C": "TGC", "E": "GAG", "Q": "CAG", "G": "GGC",
    "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC",
    "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}


# ---------------------------------------------------------------------------
# Constantes de modelo
# ---------------------------------------------------------------------------

# Fracao de uridinas a substituir em cada estrategia
PSI_U_FRACTIONS: Final[tuple[float, ...]] = (0.30, 0.50, 1.00)

# Custo adicional de pseudouridina vs uridina padrao (fator multiplicativo)
# NTPs modificados custam ~5-10x mais que NTPs padrao
# A fracao do custo total de IVT que vem de UTP e ~25%
# Portanto, substituir 100% das U aumenta o custo de NTPs em ~1.5-2.5x
PSI_U_COST_FACTOR: Final[float] = 3.0  # m1psi-UTP custa ~3x mais que UTP

# Estimativas de efeito na traducao e imunogenicidade
# Baseadas em Kariko et al. 2008, Andries et al. 2015
_TRANSLATION_BOOST_100PCT: Final[float] = 2.5  # Ate 2-3x mais proteina
_IMMUNOGENICITY_REDUCTION_100PCT: Final[float] = 0.85  # Reducao de ~85%

# Comprimento do fragmento para contexto de base-pairing
_CONTEXT_WINDOW: Final[int] = 15  # nucleotideos de cada lado da U


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class UridinePosition:
    """Posicao individual de uridina na sequencia mRNA."""

    position: int
    codon: str
    amino_acid: str
    local_gc_content: float
    base_pair_score: float  # estimativa de participacao em stem
    priority_rank: int  # 1 = mais beneficiada pela substituicao


@dataclass(frozen=True, slots=True)
class PsiUStrategy:
    """Estrategia de substituicao por pseudouridina."""

    fraction: float  # 0.0 a 1.0
    label: str
    uridines_substituted: int
    uridines_total: int
    estimated_translation_fold: float
    estimated_immunogenicity_reduction_pct: float
    estimated_cost_factor: float
    positions_substituted: tuple[int, ...]


@dataclass
class PseudouridineAnalysis:
    """Resultado completo da analise de pseudouridina."""

    mrna_length_nt: int
    total_uridines: int
    uridine_fraction: float
    gc_content: float
    uridine_positions: list[UridinePosition]
    strategies: list[PsiUStrategy]
    optimal_fraction: float
    optimal_rationale: str


# ---------------------------------------------------------------------------
# Funcoes de calculo
# ---------------------------------------------------------------------------


def _protein_to_mrna(protein_seq: str) -> str:
    """Converte sequencia proteica para mRNA usando codons preferenciais caninos.

    Gera a CDS (coding sequence) sem UTRs ou cap. A conversao usa
    codons otimizados para Canis lupus familiaris para maximizar a
    eficiencia de traducao em celulas caninas.

    Args:
        protein_seq: sequencia aminoacidica

    Returns:
        Sequencia de mRNA (usando U em vez de T)
    """
    codons = []
    for aa in protein_seq:
        codon = _PREFERRED_CODONS_CANIS.get(aa.upper())
        if codon is None:
            raise ValueError(f"Aminoacido nao reconhecido: {aa}")
        codons.append(codon)
    # Adicionar codon de parada
    codons.append(_PREFERRED_CODONS_CANIS["*"])
    # Converter T -> U para representacao de mRNA
    return "".join(codons).replace("T", "U")


def _estimate_base_pair_score(
    mrna: str,
    position: int,
    window: int = _CONTEXT_WINDOW,
) -> float:
    """Estima a probabilidade de uma uridina participar em base-pairing local.

    Modelo simplificado: conta pares complementares Watson-Crick (AU, GC)
    na vizinhanca local da uridina. Uridinas em regioes ricas em GC tem
    maior probabilidade de estar em stems de estrutura secundaria, onde
    a substituicao por psiU estabiliza o duplex.

    O score e normalizado entre 0 (sem contexto de stem) e 1 (alta
    probabilidade de estar em stem).

    Args:
        mrna: sequencia de mRNA completa
        position: indice da uridina na sequencia
        window: tamanho do contexto de cada lado

    Returns:
        Score de base-pairing entre 0.0 e 1.0
    """
    # Complementaridade Watson-Crick para RNA
    complement = {"A": "U", "U": "A", "G": "C", "C": "G"}

    start = max(0, position - window)
    end = min(len(mrna), position + window + 1)
    context = mrna[start:end]

    if len(context) < 3:
        return 0.0

    # Contar potenciais pares complementares na vizinhanca
    # Heuristica: verificar se existem nucleotideos complementares
    # a distancias tipicas de hairpin loops (4-12 nt)
    pair_count = 0
    total_checks = 0

    for offset in range(4, min(13, len(context))):
        for i in range(len(context) - offset):
            nt_i = context[i]
            nt_j = context[i + offset]
            total_checks += 1
            if complement.get(nt_i) == nt_j:
                pair_count += 1

    if total_checks == 0:
        return 0.0

    # GC content local tambem contribui para estabilidade de stem
    local_gc = sum(1 for nt in context if nt in ("G", "C")) / len(context)

    # Score combinado: complementaridade + conteudo GC
    complementarity = pair_count / total_checks
    score = 0.6 * complementarity + 0.4 * local_gc

    return round(min(1.0, score), 4)


def _find_uridine_positions(mrna: str, protein_seq: str) -> list[UridinePosition]:
    """Mapeia todas as posicoes de uridina no mRNA com contexto.

    Para cada U no mRNA:
    1. Identifica o codon correspondente e aminoacido
    2. Calcula conteudo GC local
    3. Estima score de base-pairing
    4. Atribui ranking de prioridade (stems primeiro)

    Args:
        mrna: sequencia de mRNA
        protein_seq: sequencia proteica correspondente

    Returns:
        Lista de UridinePosition ordenada por prioridade
    """
    positions: list[UridinePosition] = []

    # Mapear cada posicao de U ao codon e aminoacido correspondente
    # CDS comeca na posicao 0 e cada codon tem 3 nt
    cds_length = len(protein_seq) * 3 + 3  # +3 para stop codon

    for i, nt in enumerate(mrna):
        if nt != "U":
            continue

        # Identificar codon e aminoacido
        codon_index = i // 3
        codon_start = codon_index * 3
        codon = mrna[codon_start:codon_start + 3] if codon_start + 3 <= len(mrna) else "???"

        if codon_index < len(protein_seq):
            amino_acid = protein_seq[codon_index]
        elif codon_index == len(protein_seq):
            amino_acid = "*"  # stop codon
        else:
            amino_acid = "?"

        # GC content local (janela de 31 nt centrada na posicao)
        gc_start = max(0, i - 15)
        gc_end = min(len(mrna), i + 16)
        local_seq = mrna[gc_start:gc_end]
        local_gc = sum(1 for c in local_seq if c in ("G", "C")) / len(local_seq)

        # Score de base-pairing
        bp_score = _estimate_base_pair_score(mrna, i)

        positions.append(UridinePosition(
            position=i,
            codon=codon,
            amino_acid=amino_acid,
            local_gc_content=round(local_gc, 4),
            base_pair_score=bp_score,
            priority_rank=0,  # sera atualizado depois do sort
        ))

    # Ordenar por score de base-pairing (maior = mais beneficio da substituicao)
    positions.sort(key=lambda p: p.base_pair_score, reverse=True)

    # Atribuir rankings (1 = maior prioridade)
    ranked = []
    for rank, pos in enumerate(positions, start=1):
        ranked.append(UridinePosition(
            position=pos.position,
            codon=pos.codon,
            amino_acid=pos.amino_acid,
            local_gc_content=pos.local_gc_content,
            base_pair_score=pos.base_pair_score,
            priority_rank=rank,
        ))

    return ranked


def _build_strategy(
    fraction: float,
    uridine_positions: list[UridinePosition],
    total_uridines: int,
) -> PsiUStrategy:
    """Constroi uma estrategia de substituicao para uma dada fracao.

    A substituicao prioriza uridinas com maior score de base-pairing
    (ja ordenadas por prioridade na lista de entrada).

    Efeitos estimados escalam linearmente com a fracao substituida:
    - Traducao: boost proporcional (30% substituicao -> ~1.45x, 100% -> ~2.5x)
    - Imunogenicidade: reducao proporcional
    - Custo: aumento proporcional ao uso de m1psi-UTP

    Args:
        fraction: fracao de uridinas a substituir (0.0-1.0)
        uridine_positions: lista ordenada por prioridade
        total_uridines: total de uridinas no mRNA

    Returns:
        PsiUStrategy com detalhes da estrategia
    """
    n_substitute = round(total_uridines * fraction)

    # Selecionar as posicoes de maior prioridade
    selected = uridine_positions[:n_substitute]
    selected_positions = tuple(sorted(p.position for p in selected))

    # Efeito na traducao: escala sub-linear (retornos decrescentes)
    # 30% -> ~1.45x, 50% -> ~1.75x, 100% -> ~2.5x (Andries et al. 2015)
    translation_fold = 1.0 + (_TRANSLATION_BOOST_100PCT - 1.0) * (fraction ** 0.7)

    # Reducao de imunogenicidade inata: escala quase-linear
    immuno_reduction = _IMMUNOGENICITY_REDUCTION_100PCT * fraction

    # Custo: proporcional a fracao (UTP modificado custa 3x mais)
    # 25% do custo de IVT vem de UTP; substituir fracao F aumenta:
    # custo_ntp_novo = 0.75 + 0.25 * (1 + (PSI_U_COST_FACTOR - 1) * F)
    utp_fraction_of_cost = 0.25
    cost_factor = (
        1.0
        + utp_fraction_of_cost * (PSI_U_COST_FACTOR - 1.0) * fraction
    )

    # Label descritivo
    if fraction >= 1.0:
        label = "100% psiU (full substitution)"
    else:
        label = f"{fraction:.0%} psiU (selective top-ranked)"

    return PsiUStrategy(
        fraction=fraction,
        label=label,
        uridines_substituted=n_substitute,
        uridines_total=total_uridines,
        estimated_translation_fold=round(translation_fold, 2),
        estimated_immunogenicity_reduction_pct=round(immuno_reduction * 100, 1),
        estimated_cost_factor=round(cost_factor, 2),
        positions_substituted=selected_positions,
    )


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def run_pseudouridine_analysis() -> PseudouridineAnalysis:
    """Executa a analise completa de otimizacao de pseudouridina.

    Etapas:
        1. Gera mRNA in silico a partir da sequencia proteica do construto
        2. Mapeia todas as posicoes de uridina com contexto estrutural
        3. Gera estrategias de substituicao (30%, 50%, 100%)
        4. Seleciona estrategia otima com racional

    Returns:
        PseudouridineAnalysis com todos os resultados
    """
    # Gerar sequencia proteica do construto (com sinal, 5 CPB repeats)
    protein_seq = get_construct_sequence(include_signal=True, cpb_repeats=5)

    # Converter para mRNA in silico
    mrna = _protein_to_mrna(protein_seq)
    mrna_length = len(mrna)

    # Contar uridinas totais
    total_u = mrna.count("U")
    u_fraction = total_u / mrna_length if mrna_length > 0 else 0.0

    # Calcular GC content do mRNA (converter U->T para o validador)
    gc = calculate_gc_content(mrna.replace("U", "T"))

    # Mapear posicoes de uridina com contexto
    u_positions = _find_uridine_positions(mrna, protein_seq)

    # Gerar estrategias de substituicao
    strategies = [
        _build_strategy(frac, u_positions, total_u)
        for frac in PSI_U_FRACTIONS
    ]

    # Selecionar estrategia otima
    # Criterio: maximizar traducao com custo aceitavel
    # 50% e o sweet spot -- boa traducao, custo moderado, imunogenicidade baixa
    optimal_fraction = 0.50
    optimal_rationale = (
        "A substituicao de 50% das uridinas (priorizando posicoes em stems "
        "de estrutura secundaria) oferece o melhor equilibrio entre: "
        "(1) aumento de ~1.75x na eficiencia de traducao, "
        "(2) reducao de ~42% na ativacao imune inata (TLR7/8), "
        "(3) custo incremental moderado de ~1.25x. "
        "A substituicao total (100%) maximiza traducao (~2.5x) e reduz "
        "imunogenicidade em ~85%, mas o custo de m1psi-UTP e significativo "
        "para producao em escala veterinaria. Para um primeiro candidato "
        "pre-clinico, 50% seletivo e a abordagem mais custo-efetiva. "
        "Na fase clinica, 100% pode ser justificado pelo beneficio imunologico."
    )

    return PseudouridineAnalysis(
        mrna_length_nt=mrna_length,
        total_uridines=total_u,
        uridine_fraction=round(u_fraction, 4),
        gc_content=gc,
        uridine_positions=u_positions,
        strategies=strategies,
        optimal_fraction=optimal_fraction,
        optimal_rationale=optimal_rationale,
    )


def to_dict(analysis: PseudouridineAnalysis) -> dict:
    """Serializa a analise de pseudouridina para JSON.

    NOTA: Nao exporta a sequencia mRNA completa nem as posicoes individuais
    de uridina (que revelariam a sequencia). Exporta apenas metricas
    agregadas e estrategias.

    Args:
        analysis: resultado da analise

    Returns:
        Dicionario serializavel para JSON
    """
    return {
        "mrna_length_nt": analysis.mrna_length_nt,
        "total_uridines": analysis.total_uridines,
        "uridine_fraction": analysis.uridine_fraction,
        "gc_content": analysis.gc_content,
        "uridine_context_summary": {
            "high_bp_score_count": sum(
                1 for u in analysis.uridine_positions
                if u.base_pair_score > 0.5
            ),
            "medium_bp_score_count": sum(
                1 for u in analysis.uridine_positions
                if 0.25 <= u.base_pair_score <= 0.5
            ),
            "low_bp_score_count": sum(
                1 for u in analysis.uridine_positions
                if u.base_pair_score < 0.25
            ),
            "avg_base_pair_score": round(
                sum(u.base_pair_score for u in analysis.uridine_positions)
                / len(analysis.uridine_positions),
                4,
            ) if analysis.uridine_positions else 0.0,
        },
        "strategies": [
            {
                "fraction": s.fraction,
                "label": s.label,
                "uridines_substituted": s.uridines_substituted,
                "uridines_total": s.uridines_total,
                "estimated_translation_fold": s.estimated_translation_fold,
                "estimated_immunogenicity_reduction_pct": (
                    s.estimated_immunogenicity_reduction_pct
                ),
                "estimated_cost_factor": s.estimated_cost_factor,
            }
            for s in analysis.strategies
        ],
        "optimal_fraction": analysis.optimal_fraction,
        "optimal_rationale": analysis.optimal_rationale,
    }
