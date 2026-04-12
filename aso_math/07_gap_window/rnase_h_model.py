"""Modelo de eficiencia de clivagem por RNase H e scoring composto de gapmers.

Modela a eficiencia de RNase H como funcao do tamanho do gap de DNA central,
combinando com propriedades termodinamicas para um score composto que rankeia
configuracoes de gapmer.

Eficiencia de RNase H por tamanho de gap (Crooke et al. 2017):
    - Gap  8-9 nt:  0.60 (minimo funcional, clivagem reduzida)
    - Gap 10-12 nt: 0.90 (otimo para maioria dos ASOs)
    - Gap 13-15 nt: 1.00 (eficiencia maxima)
    - Gap 16+ nt:   0.95 (leve reducao por perda de protecao dos flancos)

Scoring composto:
    score = w_dg * norm_dg + w_tm * norm_tm + w_rnase_h * rnase_h
          + w_nuclease * nuclease_resistance

    Pesos: dG (0.3), Tm (0.2), RNase H (0.3), resistencia a nucleases (0.2)
"""

from __future__ import annotations

from dataclasses import dataclass

from aso_math.thermo import compute_dg, compute_tm

from .gap_enumerator import GapWindowConfig

# ---------------------------------------------------------------------------
# Constantes do modelo
# ---------------------------------------------------------------------------

# Boost de Tm por posicao LNA (McTigue et al. 2004)
TM_BOOST_PER_LNA: float = 3.0

# Pesos do score composto
W_DG: float = 0.30
W_TM: float = 0.20
W_RNASE_H: float = 0.30
W_NUCLEASE: float = 0.20


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ScoredConfig:
    """Resultado do scoring de uma configuracao gapmer.

    Attributes:
        config: Configuracao LNA-DNA-LNA avaliada.
        dg_binding: Energia livre de ligacao (kcal/mol).
        tm_base: Tm base sem ajuste LNA (Celsius).
        tm_adjusted: Tm ajustada com boost LNA (Celsius).
        rnase_h_efficiency: Eficiencia de clivagem RNase H [0, 1].
        nuclease_resistance: Resistencia a nucleases [0, 1].
        composite_score: Score composto final [0, 1].
    """

    config: GapWindowConfig
    dg_binding: float
    tm_base: float
    tm_adjusted: float
    rnase_h_efficiency: float
    nuclease_resistance: float
    composite_score: float


# ---------------------------------------------------------------------------
# Modelo de eficiencia RNase H
# ---------------------------------------------------------------------------


def rnase_h_efficiency(gap_size: int) -> float:
    """Calcula eficiencia de clivagem por RNase H como funcao do gap.

    Modelo baseado em Crooke et al. (2017) e Kurreck (2003):
        - Gap  8-9:  0.60 (funcional mas subotimo)
        - Gap 10-12: 0.90 (faixa otima para gapmers curtos)
        - Gap 13-15: 1.00 (eficiencia maxima)
        - Gap 16+:   0.95 (leve reducao)
        - Gap < 8:   0.00 (nao funcional)

    Args:
        gap_size: Numero de nucleotideos no gap de DNA central.

    Returns:
        Eficiencia em [0.0, 1.0].
    """
    if gap_size < 8:
        return 0.0
    if gap_size <= 9:
        return 0.6
    if gap_size <= 12:
        return 0.9
    if gap_size <= 15:
        return 1.0
    return 0.95


# ---------------------------------------------------------------------------
# Calculo de resistencia a nucleases
# ---------------------------------------------------------------------------


def nuclease_resistance(total_lna: int, max_lna: int = 10) -> float:
    """Calcula resistencia a nucleases proporcional ao total de LNA.

    LNA nos flancos protege contra exonucleases 3' e 5'. A resistencia
    e proporcional ao total de posicoes LNA, normalizada pelo maximo.

    Args:
        total_lna: Total de posicoes LNA (5' + 3').
        max_lna: Maximo esperado para normalizacao. Default: 10.

    Returns:
        Resistencia em [0.0, 1.0].
    """
    if max_lna <= 0:
        return 0.0
    return min(1.0, total_lna / max_lna)


# ---------------------------------------------------------------------------
# Normalizacao de componentes
# ---------------------------------------------------------------------------


def _normalize_dg(dg: float, dg_max_abs: float = 50.0) -> float:
    """Normaliza dG em [0, 1]. Mais negativo = melhor = mais proximo de 1."""
    return min(1.0, abs(dg) / dg_max_abs)


def _normalize_tm(tm_adjusted: float, tm_optimal: float = 75.0, half_width: float = 20.0) -> float:
    """Normaliza Tm ajustada via funcao triangular centrada no otimo.

    Pico em tm_optimal, zero em (tm_optimal - half_width) e (tm_optimal + half_width).

    Args:
        tm_adjusted: Tm com boost LNA (Celsius).
        tm_optimal: Tm ideal para ASO. Default: 75 C.
        half_width: Meia-largura da funcao triangular. Default: 20 C.

    Returns:
        Score em [0.0, 1.0].
    """
    return max(0.0, 1.0 - abs(tm_adjusted - tm_optimal) / half_width)


# ---------------------------------------------------------------------------
# Scoring composto
# ---------------------------------------------------------------------------


def score_config(
    config: GapWindowConfig,
    aso_seq: str,
    sl_seq: str,
) -> ScoredConfig:
    """Calcula score composto para uma configuracao gapmer.

    Integra 4 dimensoes com pesos pre-definidos:
        score = 0.3 * norm_dg + 0.2 * norm_tm + 0.3 * rnase_h + 0.2 * nuclease

    Args:
        config: Configuracao LNA-DNA-LNA a avaliar.
        aso_seq: Sequencia do ASO (5' -> 3').
        sl_seq: Sequencia do SL RNA alvo.

    Returns:
        ScoredConfig com todos os componentes e score final.
    """
    # Propriedades termodinamicas base
    dg = compute_dg(aso_seq)
    tm_base = compute_tm(aso_seq)

    # Ajuste de Tm por LNA
    tm_boost = config.total_lna * TM_BOOST_PER_LNA
    tm_adj = tm_base + tm_boost

    # Componentes do score
    rnase_h = rnase_h_efficiency(config.dna_gap)
    nuc_resist = nuclease_resistance(config.total_lna)

    # Normalizacao
    norm_dg = _normalize_dg(dg)
    norm_tm = _normalize_tm(tm_adj)

    # Score composto
    composite = round(
        W_DG * norm_dg
        + W_TM * norm_tm
        + W_RNASE_H * rnase_h
        + W_NUCLEASE * nuc_resist,
        4,
    )

    return ScoredConfig(
        config=config,
        dg_binding=dg,
        tm_base=tm_base,
        tm_adjusted=round(tm_adj, 2),
        rnase_h_efficiency=rnase_h,
        nuclease_resistance=round(nuc_resist, 4),
        composite_score=composite,
    )
