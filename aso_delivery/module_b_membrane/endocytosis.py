"""Modelos de endocitose para captacao de ASOs por macrofagos.

Implementa tres vias de uptake:
1. Pinocitose de fase fluida (nao-especifica, baixa eficiencia)
2. Endocitose mediada por receptor (receptores scavenger ligam PS)
3. Gymnosis (captacao livre sem transfeccao, demonstrada para PS-ASOs)

A motivacao biologica: macrofagos sao fagocitos profissionais com
atividade endocitica 5-10x maior que a maioria das celulas. Isso
torna macrofagos o tipo celular IDEAL para delivery de PS-ASOs
sem necessidade de formulacao complexa.

Alem disso, receptores scavenger classe A (SR-A) e estabilina-1/2
em macrofagos reconhecem especificamente o backbone fosforotioato,
proporcionando um mecanismo de captacao ativa e direcionada.

Para Leishmania, isso e perfeito: as amastigotas residem DENTRO
dos macrofagos infectados — exatamente as celulas que captam
preferencialmente PS-ASOs.

Referencias:
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
- Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
- Stein CA et al. (2010) Nucleic Acids Res 38(5):e3 — gymnosis
- Butler M et al. (2000) J Pharmacol Exp Ther 292(2):547-555 — SR-A binding
- Koller E et al. (2011) Nucleic Acids Res 39(11):4795-4807
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Constantes de endocitose
# ---------------------------------------------------------------------------

# Constante de Boltzmann em kcal/(mol*K)
KB_KCAL: Final[float] = 1.987e-3

# Temperatura fisiologica
T_PHYSIOLOGICAL_K: Final[float] = 310.15

# Volume tipico de um macrofago (microlitros convertido para litros)
# Ref: Macrofagos caninos: diametro ~15-20 um -> volume ~2-4 pL
MACROPHAGE_VOLUME_L: Final[float] = 3.0e-12  # 3 pL

# Volume do fagolisossomo (fração do volume celular: ~5-10%)
PHAGOLYSOSOME_VOLUME_FRACTION: Final[float] = 0.07

# Taxa de pinocitose em macrofagos (volume internalizado por hora)
# Ref: Steinman RM et al. (1976) J Cell Biol 68(3):665-687
# Macrofagos internalizam ~100% do volume celular por hora em pinocitose
# (isso e 5-10x maior que fibroblastos)
PINOCYTOSIS_VOLUME_FRACTION_PER_HOUR: Final[float] = 1.0

# Concentracao extracelular padrao de ASO em ensaios in vitro
# Ref: Crooke et al. (2017) — concentracao tipica: 1-10 uM para gymnosis
ASO_EXTRACELLULAR_CONC_M: Final[float] = 1.0e-6  # 1 uM

# ---------------------------------------------------------------------------
# Constantes de receptores scavenger
# ---------------------------------------------------------------------------

# Kd de ligacao PS-ASO a receptores scavenger classe A (SR-A)
# Ref: Butler M et al. (2000) — Kd ~10-50 nM para PS-ASOs em macrofagos
KD_SRA_BINDING_M: Final[float] = 30.0e-9  # 30 nM

# Numero de receptores scavenger por macrofago
# Ref: Platt N, Gordon S (2001) J Clin Invest 108(5):649-654
# SR-A: ~50,000-200,000 por macrofago
N_SRA_PER_MACROPHAGE: Final[int] = 100_000

# Taxa de internalizacao de receptor-ligando (reciclagem de receptores)
# Ref: Brown MS, Goldstein JL (1986) — taxa tipica: 3-5 min por ciclo
# ~12-20 ciclos/hora
RECEPTOR_INTERNALIZATION_CYCLES_PER_HOUR: Final[float] = 15.0

# Fator de concentracao endossomal (volume extracelular >> volume endossomal)
# Cada vesicula endocitica concentra o conteudo ~100-1000x
ENDOSOMAL_CONCENTRATION_FACTOR: Final[float] = 500.0

# ---------------------------------------------------------------------------
# Constantes de gymnosis
# ---------------------------------------------------------------------------

# Taxa de uptake por gymnosis (captacao livre sem transfeccao)
# Ref: Stein CA et al. (2010) — gymnosis efetiva em 1-10 uM apos 24-72h
# Mecanismo proposto: combinacao de pinocitose + adsorptive endocytosis
# A taxa e mais lenta que receptor-mediada mas constante ao longo do tempo
GYMNOSIS_RATE_FRACTION_PER_HOUR: Final[float] = 0.005  # 0.5% por hora


# ---------------------------------------------------------------------------
# Dataclasses de resultado
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EndocytosisPathway:
    """Modelo de uma via de endocitose para captacao de ASO.

    Cada via e caracterizada por:
    - Nome e mecanismo
    - Constante de taxa (k, 1/hora)
    - Eficiencia (fracao internalizada em 24 horas)
    - Concentracao intracelular resultante
    - Relevancia para macrofagos especificamente
    """

    name: str
    mechanism: str
    rate_constant_per_hour: float
    efficiency_24h: float            # fracao internalizada em 24h (0-1)
    intracellular_conc_nm: float     # concentracao no fagolisossomo (nM)
    is_dominant_in_macrophages: bool  # True se essa via domina em macrofagos
    key_reference: str
    notes: str


@dataclass(frozen=True)
class EndocytosisComparison:
    """Comparacao das tres vias de endocitose para um ASO especifico.

    Inclui as tres vias individuais mais a estimativa combinada
    de captacao total (soma das vias independentes, cap em 1.0).
    """

    aso_name: str
    extracellular_conc_um: float
    pinocytosis: EndocytosisPathway
    receptor_mediated: EndocytosisPathway
    gymnosis: EndocytosisPathway
    total_efficiency_24h: float
    total_intracellular_conc_nm: float
    dominant_pathway: str
    macrophage_advantage_factor: float  # vantagem vs celula nao-macrofago


# ---------------------------------------------------------------------------
# 1. Pinocitose de fase fluida
# ---------------------------------------------------------------------------


def compute_pinocytosis(
    extracellular_conc_m: float = ASO_EXTRACELLULAR_CONC_M,
    volume_fraction_per_hour: float = PINOCYTOSIS_VOLUME_FRACTION_PER_HOUR,
    time_hours: float = 24.0,
) -> EndocytosisPathway:
    """Modela captacao por pinocitose de fase fluida.

    Pinocitose e a internalizacao nao-seletiva de fluido extracelular
    em vesiculas. O ASO e captado simplesmente por estar presente
    no volume internalizado.

    Modelo:
    A cada hora, o macrofago internaliza ~100% do seu volume em vesiculas.
    A concentracao dentro da vesicula e igual a concentracao extracelular.
    As vesiculas sao entregues ao fagolisossomo, onde o conteudo e concentrado.

    k_pinocytosis = volume_fraction_per_hour * (V_cell / V_phagolysosome)

    Eficiencia apos t horas:
    f(t) = 1 - exp(-k * t)

    Nota: pinocitose e MUITO mais ativa em macrofagos (~5-10x) do que
    em celulas epiteliais ou fibroblastos. Isso e uma vantagem natural
    para delivery de ASOs anti-Leishmania.

    Ref: Steinman RM et al. (1976) — pinocitose em macrofagos

    Args:
        extracellular_conc_m: Concentracao extracelular em mol/L.
        volume_fraction_per_hour: Fracao do volume celular internalizado/hora.
        time_hours: Tempo de exposicao em horas.

    Returns:
        Modelo de pinocitose com taxa e eficiencia.
    """
    # Taxa efetiva: volume internalizado contribui para concentracao endossomal
    # Cada hora, V_cell de fluido e internalizado e entregue ao phagolysosome
    # que tem volume V_cell * 0.07. Fator de concentracao: 1/0.07 ~ 14x
    concentration_factor = 1.0 / PHAGOLYSOSOME_VOLUME_FRACTION

    # Constante de taxa de acumulo (considerando saturacao)
    # k = taxa_de_internalizacao * eficiencia_de_entrega
    # A eficiencia de entrega e ~50% (parte do conteudo e reciclado)
    delivery_efficiency = 0.5
    k_pino = volume_fraction_per_hour * delivery_efficiency * 0.01  # normalizado

    # Eficiencia em t horas (modelo de primeiro ordem com saturacao)
    efficiency = 1.0 - math.exp(-k_pino * time_hours)

    # Concentracao no fagolisossomo
    # [intra] = [extra] * eficiencia * fator_concentracao
    conc_intra_m = extracellular_conc_m * efficiency * concentration_factor
    conc_intra_nm = conc_intra_m * 1.0e9

    return EndocytosisPathway(
        name="Fluid-phase pinocytosis",
        mechanism=(
            "Non-selective internalization of extracellular fluid in vesicles. "
            "ASO enters passively with bulk fluid. Vesicles delivered to "
            "phagolysosome. Macrophages internalize ~100% cell volume per hour."
        ),
        rate_constant_per_hour=round(k_pino, 6),
        efficiency_24h=round(efficiency, 4),
        intracellular_conc_nm=round(conc_intra_nm, 2),
        is_dominant_in_macrophages=False,
        key_reference="Steinman RM et al. (1976) J Cell Biol 68(3):665-687",
        notes=(
            "Low efficiency per unit time but constitutive and non-saturable. "
            "5-10x more active in macrophages than in fibroblasts. "
            "Contributes ~10-15% of total ASO uptake in macrophages."
        ),
    )


# ---------------------------------------------------------------------------
# 2. Endocitose mediada por receptor
# ---------------------------------------------------------------------------


def compute_receptor_mediated_endocytosis(
    extracellular_conc_m: float = ASO_EXTRACELLULAR_CONC_M,
    kd_m: float = KD_SRA_BINDING_M,
    n_receptors: int = N_SRA_PER_MACROPHAGE,
    cycles_per_hour: float = RECEPTOR_INTERNALIZATION_CYCLES_PER_HOUR,
    time_hours: float = 24.0,
) -> EndocytosisPathway:
    """Modela captacao mediada por receptores scavenger.

    O backbone fosforotioato (PS) e reconhecido por receptores scavenger
    de classe A (SR-A) e estabilina-1/2 em macrofagos. Essa interacao
    e ESPECIFICA e de alta afinidade (Kd ~30 nM).

    Modelo Michaelis-Menten:
    Fracao_ocupada = [ASO] / ([ASO] + Kd)

    Receptores ocupados por ciclo = n_receptores * fracao_ocupada
    Moleculas internalizadas por hora = receptores_ocupados * ciclos/hora

    Concentracao intracelular:
    [intra] = moleculas_internalizadas_total / (N_A * V_phagolysosome)

    Este e o mecanismo DOMINANTE de uptake de PS-ASOs em macrofagos.
    A alta densidade de SR-A em macrofagos (~100,000 por celula) e
    a rapida reciclagem (~15 ciclos/hora) permitem captacao eficiente
    mesmo em concentracoes sub-micromolares.

    Refs:
    - Butler M et al. (2000) — SR-A liga PS com Kd ~30 nM
    - Koller E et al. (2011) — uptake via receptores em macrofagos

    Args:
        extracellular_conc_m: Concentracao extracelular (M).
        kd_m: Constante de dissociacao do receptor (M).
        n_receptors: Receptores por celula.
        cycles_per_hour: Ciclos de internalizacao por hora.
        time_hours: Tempo de exposicao.

    Returns:
        Modelo de endocitose mediada por receptor.
    """
    # Fracao de receptores ocupados (Michaelis-Menten)
    fraction_occupied = extracellular_conc_m / (extracellular_conc_m + kd_m)

    # Moleculas internalizadas por hora
    molecules_per_hour = n_receptors * fraction_occupied * cycles_per_hour

    # Modelo de estado estacionario com turnover intracelular
    # O ASO se acumula no phagolysosome mas tambem e degradado/exocitado.
    # Influxo: J_in = moleculas_por_hora
    # Efluxo: J_out = k_turnover * N_intracelular
    # Estado estacionario: N_ss = J_in / k_turnover
    #
    # k_turnover para PS-ASOs em lisossomos: ~0.03/h (meia-vida ~23h)
    # Ref: Geary RS et al. (2015) — PS-ASOs sao resistentes a nucleases
    # mas sofrem degradacao lenta e exocitose parcial
    k_turnover = 0.03  # 1/hora (combina degradacao + exocitose)

    avogadro = 6.022e23
    v_phago = MACROPHAGE_VOLUME_L * PHAGOLYSOSOME_VOLUME_FRACTION

    # Estado estacionario: N_ss moleculas no phagolysosome
    n_ss_molecules = molecules_per_hour / k_turnover

    # Concentracao no estado estacionario
    conc_ss_m = n_ss_molecules / (avogadro * v_phago)

    # Abordagem temporal: C(t) = C_ss * (1 - exp(-k_turnover * t))
    # O sistema atinge ~63% do SS em 1/k = 33h, ~95% em 3/k = 100h
    conc_at_time_m = conc_ss_m * (1.0 - math.exp(-k_turnover * time_hours))
    conc_intra_nm = conc_at_time_m * 1.0e9

    # Taxa efetiva de acumulo (para comparacao entre vias)
    k_receptor = k_turnover  # a cinetica e dominada pelo turnover

    # Eficiencia: fracao do estado estacionario alcancada em t horas
    # C(t)/C_ss = 1 - exp(-k_turnover * t)
    # Em 24h: 1 - exp(-0.03*24) = 1 - exp(-0.72) = 0.513 (51.3% do SS)
    efficiency_eff = 1.0 - math.exp(-k_turnover * time_hours)

    return EndocytosisPathway(
        name="Receptor-mediated endocytosis (scavenger receptors)",
        mechanism=(
            "PS backbone is recognized by scavenger receptor class A (SR-A) "
            "and stabilin-1/2 on macrophages. High-affinity binding (Kd ~30 nM) "
            "leads to clathrin-coated pit formation and internalization. "
            "Receptors recycle every ~4 minutes for ~15 cycles/hour."
        ),
        rate_constant_per_hour=round(k_receptor, 6),
        efficiency_24h=round(efficiency_eff, 4),
        intracellular_conc_nm=round(conc_intra_nm, 2),
        is_dominant_in_macrophages=True,
        key_reference="Butler M et al. (2000) J Pharmacol Exp Ther 292(2):547-555",
        notes=(
            f"Dominant uptake mechanism in macrophages. "
            f"Receptor occupancy at {extracellular_conc_m*1e6:.1f} uM: "
            f"{fraction_occupied*100:.1f}%. "
            f"~{n_receptors:,} SR-A receptors per macrophage with "
            f"{cycles_per_hour:.0f} internalization cycles/hour. "
            f"PS backbone provides NATURAL tropism for macrophages — "
            f"no targeting ligand needed."
        ),
    )


# ---------------------------------------------------------------------------
# 3. Gymnosis (captacao livre)
# ---------------------------------------------------------------------------


def compute_gymnosis(
    extracellular_conc_m: float = ASO_EXTRACELLULAR_CONC_M,
    rate_per_hour: float = GYMNOSIS_RATE_FRACTION_PER_HOUR,
    time_hours: float = 24.0,
) -> EndocytosisPathway:
    """Modela captacao por gymnosis (free uptake sem transfeccao).

    Gymnosis e a captacao celular de ASOs sem uso de agentes de
    transfeccao, demonstrada experimentalmente por Stein CA et al.
    (2010) para PS-ASOs em concentracoes de 1-10 uM.

    O mecanismo exato nao e completamente elucidado, mas envolve:
    - Adsorptive endocytosis (ligacao nao-especifica a superficie celular)
    - Contribuicao de pinocitose estimulada por ligacao de PS a membrana
    - Possivel internalizacao via caveolae em alguns tipos celulares

    Modelo simplificado:
    f(t) = 1 - exp(-k_gym * t)

    Para gymnosis, k_gym ~= 0.005/hora (0.5% por hora), resultando em
    ~11% de eficiencia em 24h e ~30% em 72h.

    A gymnosis e particularmente relevante porque:
    1. Nao requer formulacao complexa (LNP, conjugados)
    2. Funciona em concentracoes clinicamente alcancaveis
    3. E mais eficiente em macrofagos (maior atividade endocitica)

    Ref: Stein CA et al. (2010) Nucleic Acids Res 38(5):e3

    Args:
        extracellular_conc_m: Concentracao extracelular (M).
        rate_per_hour: Fracao internalizada por hora.
        time_hours: Tempo de exposicao.

    Returns:
        Modelo de gymnosis.
    """
    # Eficiencia (modelo exponencial)
    efficiency = 1.0 - math.exp(-rate_per_hour * time_hours)

    # Concentracao intracelular
    # Volume phagolysosome = V_cell * 7%
    v_phago = MACROPHAGE_VOLUME_L * PHAGOLYSOSOME_VOLUME_FRACTION
    concentration_factor = 1.0 / PHAGOLYSOSOME_VOLUME_FRACTION

    conc_intra_m = extracellular_conc_m * efficiency * concentration_factor
    conc_intra_nm = conc_intra_m * 1.0e9

    return EndocytosisPathway(
        name="Gymnosis (free uptake)",
        mechanism=(
            "Free cellular uptake without transfection reagents. Demonstrated "
            "for PS-ASOs at 1-10 uM concentrations. Mechanism involves "
            "adsorptive endocytosis and membrane-stimulated pinocytosis. "
            "Effective after 24-72 hours of exposure."
        ),
        rate_constant_per_hour=round(rate_per_hour, 6),
        efficiency_24h=round(efficiency, 4),
        intracellular_conc_nm=round(conc_intra_nm, 2),
        is_dominant_in_macrophages=False,
        key_reference="Stein CA et al. (2010) Nucleic Acids Res 38(5):e3",
        notes=(
            "Clinically relevant: no formulation needed. "
            "Works at achievable plasma concentrations. "
            "Particularly effective for PS-ASOs due to protein binding "
            "and membrane adsorption properties. More active in phagocytic cells."
        ),
    )


# ---------------------------------------------------------------------------
# Comparacao integrada das tres vias
# ---------------------------------------------------------------------------


def compute_endocytosis_comparison(
    aso_name: str,
    extracellular_conc_um: float = 1.0,
) -> EndocytosisComparison:
    """Compara as tres vias de endocitose para captacao de ASO.

    Calcula eficiencia individual e combinada de pinocitose,
    endocitose mediada por receptor e gymnosis. Identifica a via
    dominante e estima a vantagem de macrofagos vs outras celulas.

    A vantagem de macrofagos e calculada como:
    - SR-A e 10-50x mais expresso em macrofagos vs epiteliais
    - Pinocitose e 5-10x mais ativa em macrofagos
    - Gymnosis e 2-3x mais eficiente em fagocitos profissionais
    - Fator total estimado: ~5-10x mais uptake em macrofagos

    Ref: Crooke ST et al. (2017) — tropismo natural de PS para macrofagos

    Args:
        aso_name: Identificador do ASO.
        extracellular_conc_um: Concentracao extracelular em micromolar.

    Returns:
        Comparacao integrada das tres vias.
    """
    conc_m = extracellular_conc_um * 1.0e-6

    # Calcular cada via
    pino = compute_pinocytosis(conc_m)
    receptor = compute_receptor_mediated_endocytosis(conc_m)
    gym = compute_gymnosis(conc_m)

    # Eficiencia total (as vias sao parcialmente independentes)
    # P(total) = 1 - P(nenhuma via captura) = 1 - (1-p1)(1-p2)(1-p3)
    total_eff = 1.0 - (
        (1.0 - pino.efficiency_24h)
        * (1.0 - receptor.efficiency_24h)
        * (1.0 - gym.efficiency_24h)
    )
    total_eff = min(1.0, total_eff)

    # Concentracao total = soma das contribuicoes
    total_conc = (
        pino.intracellular_conc_nm
        + receptor.intracellular_conc_nm
        + gym.intracellular_conc_nm
    )

    # Via dominante
    pathways = {
        "pinocytosis": pino.intracellular_conc_nm,
        "receptor_mediated": receptor.intracellular_conc_nm,
        "gymnosis": gym.intracellular_conc_nm,
    }
    dominant = max(pathways, key=pathways.get)  # type: ignore[arg-type]

    # Vantagem de macrofagos
    # Em celulas nao-macrofago: SR-A e 20x menor, pinocitose 7x menor, gymnosis 2.5x menor
    non_macrophage_factor = (
        pino.intracellular_conc_nm / 7.0
        + receptor.intracellular_conc_nm / 20.0
        + gym.intracellular_conc_nm / 2.5
    )
    if non_macrophage_factor > 0:
        advantage = total_conc / non_macrophage_factor
    else:
        advantage = 1.0

    return EndocytosisComparison(
        aso_name=aso_name,
        extracellular_conc_um=extracellular_conc_um,
        pinocytosis=pino,
        receptor_mediated=receptor,
        gymnosis=gym,
        total_efficiency_24h=round(total_eff, 4),
        total_intracellular_conc_nm=round(total_conc, 2),
        dominant_pathway=dominant,
        macrophage_advantage_factor=round(advantage, 1),
    )
