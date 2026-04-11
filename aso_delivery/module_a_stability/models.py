"""Modelos de correcao por pH para estabilidade de ASOs em ambiente acido.

Implementa tres modelos independentes:
1. Correcao termodinamica por protonacao de bases (Henderson-Hasselbalb)
2. Cinetica de degradacao por nucleases (primeira ordem)
3. Estabilidade conformacional de nucleotideos LNA

A motivacao biologica: amastigotas de L. infantum residem DENTRO do
fagolisossomo de macrofagos caninos, onde o pH varia de 4.5 a 5.0 e
ha abundancia de nucleases acidas. O ASO precisa sobreviver nesse
ambiente hostil por tempo suficiente para atingir o SL RNA alvo.

Referencias:
- SantaLucia J Jr (1998) PNAS 95(4):1460-1465 — parametros nearest-neighbor
- Siegfried NA et al. (2010) Biochemistry 49(15):3225-3236 — protonacao de bases
- Eckstein F (2014) Nucleic Acids Res 42(6):3777-3788 — backbone fosforotioato
- Vester B, Wengel J (2004) Biochemistry 43(42):13233-13241 — estabilidade LNA
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes fisico-quimicas
# ---------------------------------------------------------------------------

# Constante dos gases em kcal/(mol*K)
R_KCAL: Final[float] = 1.987e-3

# Temperatura fisiologica (37 C)
T_PHYSIOLOGICAL_K: Final[float] = 310.15

# ---------------------------------------------------------------------------
# pKa de bases nucleotidicas relevantes em ambiente acido
# Ref: Siegfried NA et al. (2010) Biochemistry 49(15):3225-3236
#      Legault P, Pardi A (1997) JACS 119(28):6621-6628
# ---------------------------------------------------------------------------

# Citosina N3: protonacao forma par C+*G Hoogsteen, desestabiliza Watson-Crick
PKA_CYTOSINE_N3: Final[float] = 4.2

# Adenina N1: protonacao mais rara em pH biologico, mas relevante em pH 4.5
PKA_ADENINE_N1: Final[float] = 3.5

# Penalidade energetica por protonacao completa de cada base
# Protonacao de C desestabiliza o par C-G em ~1.5-2.5 kcal/mol
# Protonacao de A desestabiliza o par A-T em ~1.0-2.0 kcal/mol
# Ref: Moody EM et al. (2004) RNA 10(8):1174-1183
DDG_PROTONATION_C: Final[float] = 2.0   # kcal/mol por C protonada
DDG_PROTONATION_A: Final[float] = 1.5   # kcal/mol por A protonada


# ---------------------------------------------------------------------------
# 1. Modelo de protonacao por pH (Henderson-Hasselbalch)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ProtonationCorrection:
    """Correcao de energia livre por protonacao de bases em pH acido.

    Para cada base C e A na sequencia do ASO, calcula a fracao protonada
    via Henderson-Hasselbalch e aplica a penalidade energetica proporcional.

    Biologicamente: no fagolisossomo (pH 4.5), uma fracao significativa
    das citosinas tera N3 protonado, desestabilizando os pares C-G do
    duplex ASO:SL RNA. Adeninas sao menos afetadas porque pKa(A) = 3.5
    esta abaixo do pH do fagolisossomo.
    """

    ph: float
    n_cytosines: int
    n_adenines: int
    fraction_protonated_c: float
    fraction_protonated_a: float
    ddg_cytosines: float       # penalidade total por citosinas (kcal/mol)
    ddg_adenines: float        # penalidade total por adeninas (kcal/mol)
    ddg_total: float           # penalidade total (positivo = desestabiliza)


def fraction_protonated(ph: float, pka: float) -> float:
    """Calcula fracao de moleculas protonadas via Henderson-Hasselbalch.

    f = 1 / (1 + 10^(pH - pKa))

    Em pH << pKa: f -> 1.0 (totalmente protonada)
    Em pH >> pKa: f -> 0.0 (totalmente desprotonada)
    Em pH == pKa: f = 0.5

    A intuicao: quanto mais acido o ambiente, mais protons disponiveis
    para protonate a base, desestabilizando as pontes de hidrogenio
    Watson-Crick que sustentam o duplex ASO:alvo.
    """
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def compute_protonation_correction(
    ph: float,
    n_cytosines: int,
    n_adenines: int,
) -> ProtonationCorrection:
    """Calcula a correcao de delta_G por protonacao de bases em dado pH.

    Modelo:
        ddG(pH) = sum_C[ f_C(pH) * DDG_C ] + sum_A[ f_A(pH) * DDG_A ]

    Onde f_C e f_A sao fracoes protonadas de citosina e adenina,
    e DDG_C/DDG_A sao as penalidades por protonacao completa.

    Para o ASO MRL-ASO-001 (5 LNA + 15 DNA + 5 LNA, backbone PS):
    as bases LNA sao mais resistentes a protonacao porque o acucar travado
    restringe a geometria, mas a protonacao da base nitrogenada ocorre
    independentemente do acucar. Portanto, aplicamos a correcao a todas
    as C e A da sequencia.

    Args:
        ph: pH do ambiente.
        n_cytosines: Numero de citosinas na sequencia do ASO.
        n_adenines: Numero de adeninas na sequencia do ASO.

    Returns:
        ProtonationCorrection com todas as contribuicoes detalhadas.
    """
    f_c = fraction_protonated(ph, PKA_CYTOSINE_N3)
    f_a = fraction_protonated(ph, PKA_ADENINE_N1)

    ddg_c = n_cytosines * f_c * DDG_PROTONATION_C
    ddg_a = n_adenines * f_a * DDG_PROTONATION_A

    return ProtonationCorrection(
        ph=ph,
        n_cytosines=n_cytosines,
        n_adenines=n_adenines,
        fraction_protonated_c=round(f_c, 6),
        fraction_protonated_a=round(f_a, 6),
        ddg_cytosines=round(ddg_c, 4),
        ddg_adenines=round(ddg_a, 4),
        ddg_total=round(ddg_c + ddg_a, 4),
    )


def compute_dg_at_ph(
    dg_neutral: float,
    ph: float,
    n_cytosines: int,
    n_adenines: int,
) -> float:
    """Calcula delta_G corrigido para um dado pH.

    dG(pH) = dG(7.4) + ddG_protonacao(pH)

    Nota: ddG_protonacao e sempre positivo (desestabiliza),
    portanto dG(acido) > dG(neutro) (menos negativo = mais fraco).

    Args:
        dg_neutral: delta_G a pH 7.4 (kcal/mol).
        ph: pH do ambiente alvo.
        n_cytosines: Numero de citosinas no ASO.
        n_adenines: Numero de adeninas no ASO.

    Returns:
        delta_G corrigido em kcal/mol.
    """
    correction = compute_protonation_correction(ph, n_cytosines, n_adenines)
    return round(dg_neutral + correction.ddg_total, 4)


def compute_tm_at_ph(
    tm_neutral: float,
    dg_neutral: float,
    dg_at_ph: float,
    n_bp: int,
) -> float:
    """Estima Tm corrigido para um dado pH.

    Aproximacao linear: a variacao de Tm e proporcional a variacao de dG
    por par de bases.

    dTm = (dG_pH - dG_7.4) / n_bp * fator_escala

    O fator de escala (~5-8 C por kcal/mol*bp) vem da relacao empirica
    entre estabilidade termodinamica e Tm para duplexes de 20-30 nt.
    Ref: Owczarzy R et al. (2004) Biochemistry 43(12):3537-3554

    Args:
        tm_neutral: Tm a pH 7.4 (Celsius).
        dg_neutral: delta_G a pH 7.4 (kcal/mol).
        dg_at_ph: delta_G no pH alvo (kcal/mol).
        n_bp: Numero de pares de bases no duplex.

    Returns:
        Tm estimado no pH alvo (Celsius).
    """
    # Fator empirico: ~2 C por kcal/mol de desestabilizacao total do duplex
    # Ref: Owczarzy R et al. (2004) Biochemistry 43(12):3537-3554
    # Para um duplex de 25 bp, dG ~ -28 kcal/mol, Tm ~ 108 C
    # Proporcionalmente: dTm/dG_total ~ Tm/|dG| ~ 108/28 ~ 3.9 C/(kcal/mol)
    # Usamos valor conservador de 3.0
    SCALE_FACTOR = 3.0
    delta_dg = dg_at_ph - dg_neutral  # positivo quando desestabiliza
    # delta_dg positivo -> desestabilizacao -> Tm DIMINUI
    delta_tm = -delta_dg * SCALE_FACTOR
    return round(tm_neutral + delta_tm, 2)


def compute_fraction_bound(
    dg: float,
    temperature_k: float = T_PHYSIOLOGICAL_K,
) -> float:
    """Calcula fracao de ASO ligado ao alvo no equilibrio termodinamico.

    f_bound = 1 / (1 + e^(dG / RT))

    Para dG muito negativo: f -> 1.0 (quase todo ASO ligado)
    Para dG = 0: f = 0.5
    Para dG positivo: f -> 0.0 (quase todo ASO livre)

    A intuicao: mesmo com desestabilizacao por pH, se dG permanece
    suficientemente negativo, a fracao ligada e alta o bastante para
    funcao terapeutica. O limiar pratico e ~90% bound.
    """
    exponent = dg / (R_KCAL * temperature_k)
    # Protecao contra overflow: se exponent > 500, f_bound ~ 0
    if exponent > 500:
        return 0.0
    if exponent < -500:
        return 1.0
    return round(1.0 / (1.0 + math.exp(exponent)), 6)


# ---------------------------------------------------------------------------
# 2. Cinetica de degradacao por nucleases (primeira ordem)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NucleaseResistanceProfile:
    """Perfil de resistencia a nucleases para diferentes quimicas de backbone.

    Modela a degradacao como cinetica de primeira ordem:
        [ASO](t) = [ASO]_0 * exp(-k * t)

    Onde k depende da quimica do backbone e do pH:
    - PO (fosfodiester nativo): k ~ 0.5/h a pH 4.5 — degradacao rapida
    - PS (fosforotioato):       k ~ 0.001/h a pH 4.5 — 500x mais resistente
    - LNA flanking:             k ~ 0.0001/h — 10x mais resistente que PS

    Ref: Eckstein F (2014) Nucleic Acids Res 42(6):3777-3788
         Kurreck J (2003) Eur J Biochem 270(8):1628-1644
    """

    backbone_type: str
    k_degradation: float          # constante de degradacao (1/hora)
    half_life_hours: float        # meia-vida em horas
    fraction_at_24h: float        # fracao intacta apos 24 horas
    fraction_at_48h: float        # fracao intacta apos 48 horas
    therapeutic_window_met: bool  # True se meia-vida >= 24 horas


# Constantes de degradacao em pH acido (4.5)
# Valores derivados de literatura consolidada
# Ref: Kurreck J (2003) Eur J Biochem 270(8):1628-1644
#      Stein CA et al. (1988) Nucleic Acids Res 16(8):3209-3221
K_PO_PH45: Final[float] = 0.5       # backbone fosfodiester nativo
K_PS_PH45: Final[float] = 0.001     # backbone fosforotioato
K_LNA_FACTOR: Final[float] = 0.1    # fator de reducao para regioes LNA


@dataclass(frozen=True)
class GapmerDegradationModel:
    """Modelo de degradacao segmentada para gapmer LNA-DNA-LNA.

    O gapmer tem tres regioes com resistencias diferentes:
    - Flancos LNA (5' e 3'): k_LNA = k_PS * 0.1
    - Gap DNA central: k_DNA = k_PS

    A degradacao do gapmer inteiro e dominada pelo segmento mais vulneravel
    (o gap DNA), mas o ASO so perde funcao quando o gap inteiro e degradado.
    Modelamos com a constante efetiva do gap central.

    Para MRL-ASO-001 (5+15+5 gapmer):
    - Os flancos LNA protegem as extremidades contra exonucleases 3'->5' e 5'->3'
    - O gap central e vulneravel a endonucleases, mas PS confere resistencia
    - A constante efetiva e aproximada pela do gap PS
    """

    lna_5prime_length: int
    dna_gap_length: int
    lna_3prime_length: int
    k_lna_flank: float
    k_dna_gap: float
    k_effective: float     # constante efetiva para o gapmer inteiro
    half_life_hours: float


def compute_half_life(k: float) -> float:
    """Calcula meia-vida a partir da constante de degradacao.

    t_1/2 = ln(2) / k

    Para k = 0 (infinitamente estavel), retorna infinito.
    """
    if k <= 0:
        return float("inf")
    return round(math.log(2) / k, 2)


def compute_fraction_remaining(k: float, t_hours: float) -> float:
    """Calcula fracao de ASO intacto apos t horas.

    f(t) = exp(-k * t)

    Cinetica de primeira ordem: a taxa de degradacao e proporcional
    a concentracao atual de ASO intacto.
    """
    return round(math.exp(-k * t_hours), 6)


def compute_nuclease_resistance(
    backbone_type: str,
    k_degradation: float,
) -> NucleaseResistanceProfile:
    """Calcula perfil completo de resistencia a nucleases.

    Args:
        backbone_type: Tipo de backbone ("PO", "PS", "LNA-DNA-LNA gapmer").
        k_degradation: Constante de degradacao em 1/hora.

    Returns:
        Perfil de resistencia com meia-vida e janela terapeutica.
    """
    half_life = compute_half_life(k_degradation)
    f_24h = compute_fraction_remaining(k_degradation, 24.0)
    f_48h = compute_fraction_remaining(k_degradation, 48.0)

    # Janela terapeutica: ASO precisa estar ativo por >= 24 horas
    # Consideramos ativo se meia-vida >= 24 horas
    therapeutic_ok = half_life >= 24.0

    return NucleaseResistanceProfile(
        backbone_type=backbone_type,
        k_degradation=k_degradation,
        half_life_hours=half_life,
        fraction_at_24h=f_24h,
        fraction_at_48h=f_48h,
        therapeutic_window_met=therapeutic_ok,
    )


def compute_gapmer_degradation(
    lna_5prime: int = 5,
    dna_gap: int = 15,
    lna_3prime: int = 5,
    k_ps: float = K_PS_PH45,
) -> GapmerDegradationModel:
    """Modela degradacao segmentada de um gapmer LNA-DNA-LNA.

    O modelo considera que:
    1. Flancos LNA sao 10x mais resistentes que PS puro
       (bridge metileno impede acesso de nucleases ao backbone)
    2. Gap DNA central tem resistencia PS padrao
    3. A constante efetiva do gapmer e dominada pelo gap central,
       pois e o segmento rate-limiting para perda de funcao

    No contexto do fagolisossomo: as nucleases acidas (catepsinas,
    DNase II) atacam preferencialmente terminais 3' e ligacoes PO.
    O backbone PS resiste porque o enxofre impede a coordenacao do
    ion metalico catalitico no sitio ativo da nuclease.

    Args:
        lna_5prime: Numero de nucleotideos LNA no flanco 5'.
        dna_gap: Numero de nucleotideos DNA no gap central.
        lna_3prime: Numero de nucleotideos LNA no flanco 3'.
        k_ps: Constante de degradacao para backbone PS puro.

    Returns:
        Modelo de degradacao detalhado do gapmer.
    """
    k_lna = k_ps * K_LNA_FACTOR
    k_dna = k_ps

    # A constante efetiva: media ponderada pelo comprimento de cada segmento,
    # mas dominada pelo gap central (que e o elo mais fraco da cadeia)
    total_length = lna_5prime + dna_gap + lna_3prime
    k_eff = (
        lna_5prime * k_lna + dna_gap * k_dna + lna_3prime * k_lna
    ) / total_length

    half_life = compute_half_life(k_eff)

    return GapmerDegradationModel(
        lna_5prime_length=lna_5prime,
        dna_gap_length=dna_gap,
        lna_3prime_length=lna_3prime,
        k_lna_flank=round(k_lna, 8),
        k_dna_gap=round(k_dna, 6),
        k_effective=round(k_eff, 8),
        half_life_hours=half_life,
    )


# ---------------------------------------------------------------------------
# 3. Estabilidade conformacional LNA
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LNAConformationProfile:
    """Perfil de estabilidade conformacional de nucleotideos LNA.

    LNA (Locked Nucleic Acid) tem uma ponte metileno C2'-C4' que trava
    o acucar na conformacao C3'-endo (North). Em pH acido:

    - A protonacao de bases pode enfraquecer H-bonds Watson-Crick
    - MAS a restricao covalente da ponte metileno MANTEM a geometria C3'-endo
    - Isso significa que LNA e estruturalmente mais resistente a pH acido
      do que DNA convencional (cujo acucar pode repuckering para C2'-endo)

    O parametro chave e a fracao de nucleotideos LNA que mantem C3'-endo
    em funcao do pH. Para LNA, essa fracao e essencialmente 100% porque
    a ponte covalente nao e afetada por protonacao.
    """

    ph: float
    n_lna_residues: int
    n_dna_residues: int
    c3_endo_fraction_lna: float    # fracao LNA em C3'-endo
    c3_endo_fraction_dna: float    # fracao DNA em C3'-endo (para comparacao)
    geometry_maintained: bool       # True se >= 95% LNA mantendo C3'-endo
    free_energy_penalty: float     # penalidade de repuckering (kcal/mol)


# Energia livre de restricao conformacional
# Ref: Vester B, Wengel J (2004) Biochemistry 43(42):13233-13241
#      Bondensgaard K et al. (2000) Chemistry 6(15):2687-2695
# A ponte metileno adiciona ~2-4 kcal/mol de estabilidade por residuo LNA
DG_LNA_LOCK: Final[float] = -3.0   # kcal/mol por residuo LNA

# Fracao de DNA livre que adota C3'-endo em funcao do pH
# Em pH neutro: ~30-40% (equilibrio N/S)
# Em pH acido: protonacao pode favorecer levemente C2'-endo
DNA_C3_ENDO_NEUTRAL: Final[float] = 0.36
DNA_C3_ENDO_PH_SHIFT_PER_UNIT: Final[float] = 0.03  # reducao por unidade de pH


def compute_lna_conformation(
    ph: float,
    n_lna: int,
    n_dna: int,
) -> LNAConformationProfile:
    """Calcula perfil conformacional LNA em funcao do pH.

    O modelo:
    - LNA: C3'-endo e fixado pela ponte metileno (fracao ~ 0.99-1.00)
      independente do pH. Apenas desnaturacao termica extrema pode
      quebrar a ponte covalente.
    - DNA: C3'-endo depende do equilibrio conformacional, que e
      levemente perturbado por protonacao em pH acido.

    A vantagem biologica e clara: no fagolisossomo acido, os flancos
    LNA do gapmer mantem geometria ideal para hibridizacao com o alvo,
    enquanto os flancos DNA equivalentes perderiam conformacao.

    Args:
        ph: pH do ambiente.
        n_lna: Numero de residuos LNA (flancos do gapmer).
        n_dna: Numero de residuos DNA (gap central do gapmer).

    Returns:
        Perfil conformacional com fracoes C3'-endo e penalidade energetica.
    """
    # LNA: ponte metileno garante C3'-endo independente do pH
    # Pequena correcao estatistica: ~1% pode ter distorcao por crowding
    lna_c3_endo = 0.99

    # DNA: equilibrio N/S perturbado por pH acido
    # Em pH neutro: ~36% C3'-endo (Altona & Sundaralingam 1972)
    # Cada unidade de pH abaixo de 7.4 reduz ~3%
    ph_shift = max(0.0, 7.4 - ph)
    dna_c3_endo = max(0.10, DNA_C3_ENDO_NEUTRAL - ph_shift * DNA_C3_ENDO_PH_SHIFT_PER_UNIT)

    # Penalidade energetica: residuos DNA que perdem C3'-endo pagam
    # uma penalidade de ~0.5 kcal/mol por residuo
    # (perda de pre-organizacao para hibridizacao tipo A)
    dna_penalty_per_residue = 0.5
    fraction_dna_lost = max(0.0, DNA_C3_ENDO_NEUTRAL - dna_c3_endo)
    penalty = n_dna * fraction_dna_lost * dna_penalty_per_residue

    # Geometria considerada mantida se >= 95% dos LNA em C3'-endo
    geometry_ok = lna_c3_endo >= 0.95

    return LNAConformationProfile(
        ph=ph,
        n_lna_residues=n_lna,
        n_dna_residues=n_dna,
        c3_endo_fraction_lna=round(lna_c3_endo, 4),
        c3_endo_fraction_dna=round(dna_c3_endo, 4),
        geometry_maintained=geometry_ok,
        free_energy_penalty=round(penalty, 4),
    )
