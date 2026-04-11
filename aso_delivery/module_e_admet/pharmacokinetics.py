"""Modelos farmacocineticos para PS-ASOs em Canis lupus familiaris.

Implementa modelagem PK de dois compartimentos para MRL-ASO-001
(gapmer LNA-DNA-LNA, backbone PS completo) apos injecao subcutanea
em caes de 15 kg com leishmaniose visceral.

Os PS-ASOs tem farmacocinetica unica comparada a small molecules:
- Absorcao SC rapida e quase completa (~80-95% biodisponibilidade)
- Ligacao extensa a proteinas plasmaticas (albumina >85%)
- Distribuicao preferencial para figado > rim > baco (tropismo PS)
- NAO metabolizados por CYP450 (vantagem sobre small molecules)
- Degradados por endo/exonucleases (flancos LNA protegem)
- Excrecao renal de metabolitos encurtados

O tropismo hepatico/esplenico e uma VANTAGEM para leishmaniose:
L. infantum reside primariamente em macrofagos do figado e baco,
exatamente os tecidos onde PS-ASOs se acumulam naturalmente.

Modelo PK de dois compartimentos:
- Compartimento central: plasma + tecidos bem perfundidos
- Compartimento periferico: tecidos-alvo (figado, baco, rim)
- Absorcao SC: primeira ordem com ka ~ 0.5-1.0 h^-1

Referencias:
- Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
- Yu RZ et al. (2007) J Pharmacol Exp Ther 320(1):108-116
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
- Bennett CF (2019) Annu Rev Pharmacol Toxicol 59:447-464
- Mipomersen label (FDA, 2013)
- Inotersen label (FDA, 2018)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np

# ---------------------------------------------------------------------------
# Constantes do animal modelo
# ---------------------------------------------------------------------------

# Cao de referencia para leishmaniose visceral canina
BODY_WEIGHT_KG: Final[float] = 15.0

# Volume sanguineo: ~80 mL/kg em caes (vs ~70 mL/kg em humanos)
# Ref: Doxey DL (1983) Vet Clin Pathol
BLOOD_VOLUME_ML_PER_KG: Final[float] = 80.0
BLOOD_VOLUME_ML: Final[float] = BODY_WEIGHT_KG * BLOOD_VOLUME_ML_PER_KG

# Albumina serica canina: ~2.5-3.5 g/dL (similar a humana)
PLASMA_ALBUMIN_G_DL: Final[float] = 3.0

# ---------------------------------------------------------------------------
# Constantes do ASO
# ---------------------------------------------------------------------------

# Peso molecular: ~8500 Da para 25-mer com LNA + PS
# Calculo: 25 nt * ~330 Da/nt (media) + contribuicoes LNA e PS
ASO_MW_DA: Final[float] = 8500.0

# Concentracao de formulacao para injecao SC
# Standard veterinario: 200 mg/mL para ASOs concentrados
FORMULATION_CONC_MG_ML: Final[float] = 200.0

# Numero de ligacoes fosforotioato (24 para 25-mer)
N_PS_LINKAGES: Final[int] = 24

# ---------------------------------------------------------------------------
# Parametros de absorcao (SC)
# ---------------------------------------------------------------------------

# Biodisponibilidade SC de PS-ASOs: 80-95%
# Ref: Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
# Mipomersen SC bioavailability: ~84% (FDA label)
# Inotersen SC bioavailability: ~90% (FDA label)
BIOAVAILABILITY_SC_LOW: Final[float] = 0.80
BIOAVAILABILITY_SC_HIGH: Final[float] = 0.95
BIOAVAILABILITY_SC_MEAN: Final[float] = 0.87

# Constante de absorcao SC: primeira ordem, ka ~ 0.5-1.0 h^-1
# Ref: Yu RZ et al. (2007) J Pharmacol Exp Ther 320(1):108-116
# Tmax SC = 2-4 horas -> ka = ln(2)/Tmax ~ 0.23-0.35 (Tmax-based)
# Modelo direto: ka ~ 0.5-1.0 h^-1 (dados PK populacionais)
KA_SC_LOW: Final[float] = 0.5    # h^-1
KA_SC_HIGH: Final[float] = 1.0   # h^-1
KA_SC_MEAN: Final[float] = 0.7   # h^-1

# Tmax pos-SC em caes e primatas
TMAX_SC_HOURS_LOW: Final[float] = 2.0
TMAX_SC_HOURS_HIGH: Final[float] = 4.0

# ---------------------------------------------------------------------------
# Parametros de distribuicao
# ---------------------------------------------------------------------------

# Ligacao a proteinas plasmaticas: PS-ASOs ligam extensivamente a albumina
# Ref: Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
# PS backbone confere carga negativa -> interacao com albumina
# Fracao ligada: >85% em concentracoes terapeuticas
PROTEIN_BINDING_FRACTION: Final[float] = 0.90
FRACTION_UNBOUND: Final[float] = 1.0 - PROTEIN_BINDING_FRACTION

# Volume de distribuicao no estado estacionario (Vdss)
# Para PS-ASOs: tipicamente 0.1-0.5 L/kg (ligacao a proteinas limita Vd)
# Mipomersen: Vdss = ~0.2 L/kg
# Inotersen: Vdss = ~0.26 L/kg
# Ref: Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
VD_SS_L_PER_KG: Final[float] = 0.25

# Coeficientes de particionamento tecido:plasma para PS-ASOs
# Ref: Geary RS et al. (2015) — dados de autoradiografia em macacos/caes
# PS-ASOs acumulam preferencialmente em orgaos com sinusoides fenestrados
TISSUE_PARTITION_COEFFICIENTS: Final[dict[str, float]] = {
    "liver": 30.0,         # maior concentracao — sinusoides hepáticos fenestrados
    "kidney_cortex": 25.0, # proximal tubule reabsorption de oligos filtrados
    "spleen": 15.0,        # sinusoides esplénicos — RELEVANTE para leishmaniose
    "lymph_nodes": 10.0,   # drenagem linfatica — RELEVANTE para leishmaniose
    "bone_marrow": 5.0,    # celulas reticuloendoteliais
    "heart": 1.0,          # minimal accumulation
    "muscle": 0.5,         # baixa perfusao relativa
    "brain": 0.01,         # barreira hematoencefalica impede passagem
}

# ---------------------------------------------------------------------------
# Parametros de metabolismo
# ---------------------------------------------------------------------------

# PS-ASOs NAO sao metabolizados por enzimas CYP450
# Degradados por nucleases (3' exonuclease primariamente, endonucleases)
# LNA flanking regions conferem resistencia a exonucleases
# Ref: Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563

# Meia-vida em tecido (nao plasma) — domina a PK terminal
# PS-ASOs em caes: meia-vida tecidual de 2-4 semanas
# Ref: Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
TISSUE_HALF_LIFE_DAYS_LOW: Final[float] = 14.0   # 2 semanas
TISSUE_HALF_LIFE_DAYS_HIGH: Final[float] = 30.0  # ~4 semanas
# MRL-ASO-001 com LNA: espera-se extremo superior pela protecao LNA
TISSUE_HALF_LIFE_DAYS_MEAN: Final[float] = 21.0  # 3 semanas

# Meia-vida plasmatica (fase alpha de distribuicao rapida)
# PS-ASOs: t1/2 plasma = 1-4 horas (distribuicao rapida para tecidos)
PLASMA_HALF_LIFE_HOURS: Final[float] = 2.5

# ---------------------------------------------------------------------------
# Parametros de excrecao
# ---------------------------------------------------------------------------

# Clearance renal de metabolitos encurtados
# PS-ASOs intactos sao pouco filtrados (ligacao a proteinas)
# Metabolitos menores (<10 nt) sao filtrados e excretados
# Ref: Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
# CL total em caes: ~1-3 mL/min/kg
CL_TOTAL_ML_MIN_KG: Final[float] = 2.0
CL_RENAL_FRACTION: Final[float] = 0.70   # ~70% excrecao renal (metabolitos)
CL_HEPATIC_FRACTION: Final[float] = 0.30  # ~30% via biliar


# ===========================================================================
# Dataclasses de resultado
# ===========================================================================


@dataclass(frozen=True)
class AbsorptionProfile:
    """Perfil de absorcao SC do ASO.

    Modela a absorcao pos-injecao subcutanea como cinetica de
    primeira ordem com biodisponibilidade conhecida de PS-ASOs.
    """

    route: str
    dose_mg_kg: float
    dose_absolute_mg: float
    injection_volume_ml: float
    formulation_concentration_mg_ml: float
    bioavailability_percent: float
    dose_absorbed_mg: float
    ka_per_hour: float
    tmax_hours_range: tuple[float, float]
    cmax_estimated_ng_ml: float
    auc_estimated_ng_h_ml: float


@dataclass(frozen=True)
class DistributionProfile:
    """Perfil de distribuicao do ASO em tecidos caninos.

    Modelo de dois compartimentos com coeficientes de particionamento
    tecido:plasma derivados de dados de autoradiografia de PS-ASOs.
    """

    vd_ss_l_per_kg: float
    vd_ss_l_total: float
    protein_binding_percent: float
    fraction_unbound: float
    tissue_concentrations: dict[str, float]
    tissue_partition_coefficients: dict[str, float]
    two_compartment_params: dict[str, float]


@dataclass(frozen=True)
class MetabolismProfile:
    """Perfil metabolico do ASO.

    PS-ASOs sao degradados por nucleases, NAO por CYP450.
    Sem interacoes farmacologicas com small molecules.
    """

    cyp450_metabolism: bool
    primary_degradation_pathway: str
    nuclease_protection: str
    tissue_half_life_days: float
    plasma_half_life_hours: float
    metabolites: list[str]
    drug_drug_interactions: str


@dataclass(frozen=True)
class ExcretionProfile:
    """Perfil de excrecao do ASO.

    Eliminacao primaria via renal de metabolitos encurtados.
    ASO intacto e minimamente excretado (retido por proteinas).
    """

    cl_total_ml_min_kg: float
    cl_total_ml_min: float
    cl_renal_ml_min: float
    cl_hepatic_ml_min: float
    renal_fraction: float
    terminal_half_life_days: float
    terminal_half_life_hours: float
    time_to_steady_state_days: float


@dataclass(frozen=True)
class PKTimePoint:
    """Ponto temporal da simulacao PK de dois compartimentos."""

    time_hours: float
    plasma_conc_ng_ml: float
    tissue_conc_ng_ml: float
    amount_absorbed_mg: float
    amount_eliminated_mg: float


@dataclass(frozen=True)
class DosingRegimen:
    """Regime posologico recomendado para MRL-ASO-001 em caes.

    Baseado em escalonamento alometrico de mipomersen/inotersen
    e consideracoes especificas de leishmaniose.
    """

    loading_dose_mg_kg: float
    loading_dose_mg: float
    loading_frequency: str
    loading_duration_weeks: int
    maintenance_dose_mg_kg: float
    maintenance_dose_mg: float
    maintenance_frequency: str
    treatment_duration_weeks: int
    total_treatment_weeks: int
    injection_volume_loading_ml: float
    injection_volume_maintenance_ml: float
    estimated_trough_ng_ml: float
    therapeutic_index: float
    comparison_miltefosine: dict[str, str]


# ===========================================================================
# Funcoes de calculo
# ===========================================================================


def compute_absorption(
    dose_mg_kg: float = 5.0,
    body_weight_kg: float = BODY_WEIGHT_KG,
) -> AbsorptionProfile:
    """Calcula perfil de absorcao SC para dose especificada.

    A absorcao SC de PS-ASOs e rapida e previsivel:
    - Deposito SC -> drenagem linfatica + capilar -> plasma
    - Ka ~ 0.5-1.0 h^-1 (primeira ordem)
    - Biodisponibilidade ~87% (media de mipomersen e inotersen)
    - Tmax = 2-4 horas

    Cmax estimado via modelo de 1 compartimento simplificado:
        Cmax = (F * Dose * ka) / (Vd * (ka - ke)) * (exp(-ke*Tmax) - exp(-ka*Tmax))

    Simplificacao para estimativa:
        Cmax ~ (F * Dose) / Vd (pico de distribuicao inicial)

    Args:
        dose_mg_kg: Dose em mg/kg.
        body_weight_kg: Peso do animal em kg.

    Returns:
        Perfil de absorcao completo.
    """
    dose_abs = dose_mg_kg * body_weight_kg
    inj_vol = dose_abs / FORMULATION_CONC_MG_ML
    dose_absorbed = dose_abs * BIOAVAILABILITY_SC_MEAN

    # Cmax estimado (modelo simplificado)
    # Vd em mL para conversao de unidades
    vd_ml = VD_SS_L_PER_KG * body_weight_kg * 1000.0
    # Cmax ~ F * Dose / Vd (em ng/mL: dose em ng = dose_mg * 1e6)
    cmax_ng_ml = (dose_absorbed * 1e6) / vd_ml

    # AUC estimado (modelo de 1 compartimento):
    # AUC = F * Dose / CL
    cl_ml_min = CL_TOTAL_ML_MIN_KG * body_weight_kg
    cl_ml_h = cl_ml_min * 60.0  # converter para mL/h
    auc_ng_h_ml = (dose_absorbed * 1e6) / cl_ml_h

    return AbsorptionProfile(
        route="subcutaneous injection",
        dose_mg_kg=dose_mg_kg,
        dose_absolute_mg=round(dose_abs, 1),
        injection_volume_ml=round(inj_vol, 2),
        formulation_concentration_mg_ml=FORMULATION_CONC_MG_ML,
        bioavailability_percent=round(BIOAVAILABILITY_SC_MEAN * 100, 1),
        dose_absorbed_mg=round(dose_absorbed, 1),
        ka_per_hour=KA_SC_MEAN,
        tmax_hours_range=(TMAX_SC_HOURS_LOW, TMAX_SC_HOURS_HIGH),
        cmax_estimated_ng_ml=round(cmax_ng_ml, 1),
        auc_estimated_ng_h_ml=round(auc_ng_h_ml, 1),
    )


def compute_distribution(
    dose_absorbed_mg: float,
    body_weight_kg: float = BODY_WEIGHT_KG,
) -> DistributionProfile:
    """Calcula perfil de distribuicao em tecidos caninos.

    PS-ASOs distribuem-se preferencialmente para orgaos com
    endotelio fenestrado (figado, baco, rim). Esta e uma
    vantagem natural para leishmaniose visceral:

    L. infantum coloniza macrofagos no figado (celulas de Kupffer),
    baco (zona marginal) e medula ossea — os mesmos tecidos onde
    PS-ASOs acumulam em concentracoes 15-30x plasmaticas.

    Modelo de dois compartimentos:
    - k12 (plasma -> tecido): ~0.5 h^-1 (distribuicao rapida)
    - k21 (tecido -> plasma): ~0.01 h^-1 (retorno lento)
    - ke (eliminacao do central): CL/Vc

    Args:
        dose_absorbed_mg: Dose absorvida em mg.
        body_weight_kg: Peso do animal em kg.

    Returns:
        Perfil de distribuicao com concentracoes teciduais.
    """
    vd_total = VD_SS_L_PER_KG * body_weight_kg

    # Concentracao plasmatica de referencia (Cmax)
    vd_ml = vd_total * 1000.0
    c_plasma_ng_ml = (dose_absorbed_mg * 1e6) / vd_ml

    # Concentracoes teciduais (Kp * C_plasma)
    tissue_concs: dict[str, float] = {}
    for tissue, kp in TISSUE_PARTITION_COEFFICIENTS.items():
        tissue_concs[tissue] = round(c_plasma_ng_ml * kp, 1)

    # Parametros do modelo de dois compartimentos
    # Vc (volume central) ~ 0.1 L/kg (plasma + tecidos bem perfundidos)
    vc_l_kg = 0.10
    # Vt (volume periferico) = Vdss - Vc
    vt_l_kg = VD_SS_L_PER_KG - vc_l_kg

    # Constantes de transferencia intercompartimentais
    # k12: transferencia plasma -> tecido (rapida para PS-ASOs)
    # k21: transferencia tecido -> plasma (lenta — acumulacao tecidual)
    # ke: eliminacao do compartimento central
    cl_l_h_kg = (CL_TOTAL_ML_MIN_KG * 60.0) / 1000.0
    ke = cl_l_h_kg / vc_l_kg
    k12 = 0.50  # h^-1 — distribuicao rapida
    k21 = 0.01  # h^-1 — retorno muito lento

    two_comp = {
        "vc_l_per_kg": vc_l_kg,
        "vt_l_per_kg": round(vt_l_kg, 3),
        "k12_per_hour": k12,
        "k21_per_hour": k21,
        "ke_per_hour": round(ke, 4),
        "alpha_per_hour": round(
            0.5 * (k12 + k21 + ke + math.sqrt((k12 + k21 + ke) ** 2 - 4 * k21 * ke)),
            4,
        ),
        "beta_per_hour": round(
            0.5 * (k12 + k21 + ke - math.sqrt((k12 + k21 + ke) ** 2 - 4 * k21 * ke)),
            4,
        ),
    }

    return DistributionProfile(
        vd_ss_l_per_kg=VD_SS_L_PER_KG,
        vd_ss_l_total=round(vd_total, 2),
        protein_binding_percent=round(PROTEIN_BINDING_FRACTION * 100, 1),
        fraction_unbound=FRACTION_UNBOUND,
        tissue_concentrations=tissue_concs,
        tissue_partition_coefficients=dict(TISSUE_PARTITION_COEFFICIENTS),
        two_compartment_params=two_comp,
    )


def compute_metabolism() -> MetabolismProfile:
    """Calcula perfil metabolico do ASO.

    Diferenca fundamental entre ASOs e small molecules:
    - Small molecules: metabolismo hepatico via CYP450, meia-vida curta
    - PS-ASOs: degradacao por nucleases, sem CYP450, sem DDI

    Para MRL-ASO-001 (LNA-DNA-LNA gapmer, PS):
    - Flancos LNA protegem contra exonucleases 3'->5' e 5'->3'
    - Gap DNA central e vulneravel a endonucleases
    - Metabolitos: oligonucleotideos progressivamente mais curtos
    - Eventualmente: nucleotideos individuais (reciclados)

    A ausencia de metabolismo CYP450 e uma vantagem clinica:
    - Caes com leishmaniose frequentemente recebem outros farmacos
    - Alopurinol (xantina oxidase), antimoniais, miltefosina
    - MRL-ASO-001 pode ser combinado sem risco de DDI farmacocinetica

    Returns:
        Perfil metabolico completo.
    """
    return MetabolismProfile(
        cyp450_metabolism=False,
        primary_degradation_pathway=(
            "Endonuclease-mediated cleavage of DNA gap region, "
            "followed by 3'-exonuclease trimming of fragments. "
            "LNA flanking regions (5+5 nt) resist exonuclease degradation."
        ),
        nuclease_protection=(
            "LNA-DNA-LNA gapmer design: methylene-bridged LNA at 5' and 3' "
            "flanks blocks exonuclease access. PS backbone throughout "
            "resists endonuclease cleavage (sulfur displaces catalytic metal ion). "
            "Estimated tissue half-life ~21 days (Module A: t1/2 = 1083 h in "
            "phagolysosomal conditions)."
        ),
        tissue_half_life_days=TISSUE_HALF_LIFE_DAYS_MEAN,
        plasma_half_life_hours=PLASMA_HALF_LIFE_HOURS,
        metabolites=[
            "Shortened oligonucleotides (15-20 nt) — initial endonuclease products",
            "Truncated fragments (8-14 nt) — secondary exonuclease products",
            "Short oligomers (<8 nt) — renally filtered and excreted",
            "Individual nucleotides/nucleosides — recycled via salvage pathway",
        ],
        drug_drug_interactions=(
            "No CYP450 metabolism — no expected pharmacokinetic DDI with "
            "co-administered drugs (allopurinol, miltefosine, antimonials). "
            "This is a significant clinical advantage for combination therapy "
            "in canine leishmaniasis, where multimodal treatment is standard."
        ),
    )


def compute_excretion(
    body_weight_kg: float = BODY_WEIGHT_KG,
) -> ExcretionProfile:
    """Calcula perfil de excrecao do ASO.

    PS-ASOs intactos sao pouco excretados (retidos por ligacao proteica).
    A excrecao primaria e de metabolitos encurtados (<10 nt) via renal.

    Sequencia de eliminacao:
    1. ASO intacto liga a proteinas -> nao filtrado
    2. Nucleases degradam ASO em fragmentos menores
    3. Fragmentos <10 nt perdem afinidade por proteinas
    4. Fragmentos livres sao filtrados glomerularmente
    5. Excrecao urinaria de oligonucleotideos curtos e nucleotideos

    Tempo para estado estacionario: ~5 meias-vidas terminais
    Com t1/2 terminal = 21 dias: steady state em ~105 dias (~15 semanas)

    Args:
        body_weight_kg: Peso do animal em kg.

    Returns:
        Perfil de excrecao completo.
    """
    cl_total = CL_TOTAL_ML_MIN_KG * body_weight_kg
    cl_renal = cl_total * CL_RENAL_FRACTION
    cl_hepatic = cl_total * CL_HEPATIC_FRACTION

    t_half_terminal_days = TISSUE_HALF_LIFE_DAYS_MEAN
    t_half_terminal_hours = t_half_terminal_days * 24.0

    # Tempo para estado estacionario: ~5 meias-vidas
    t_ss_days = 5.0 * t_half_terminal_days

    return ExcretionProfile(
        cl_total_ml_min_kg=CL_TOTAL_ML_MIN_KG,
        cl_total_ml_min=round(cl_total, 1),
        cl_renal_ml_min=round(cl_renal, 1),
        cl_hepatic_ml_min=round(cl_hepatic, 1),
        renal_fraction=CL_RENAL_FRACTION,
        terminal_half_life_days=t_half_terminal_days,
        terminal_half_life_hours=t_half_terminal_hours,
        time_to_steady_state_days=t_ss_days,
    )


def compute_pk_simulation(
    dose_mg: float,
    bioavailability: float = BIOAVAILABILITY_SC_MEAN,
    body_weight_kg: float = BODY_WEIGHT_KG,
    duration_hours: float = 168.0,
    dt_hours: float = 0.5,
) -> list[PKTimePoint]:
    """Simula concentracao plasmatica e tecidual ao longo do tempo.

    Modelo de dois compartimentos com absorcao SC de primeira ordem:

        dA_sc/dt  = -ka * A_sc
        dA_c/dt   = ka * A_sc - (k12 + ke) * A_c + k21 * A_t
        dA_t/dt   = k12 * A_c - k21 * A_t

    Onde:
        A_sc = quantidade no deposito SC
        A_c  = quantidade no compartimento central
        A_t  = quantidade no compartimento periferico (tecidos)
        ka   = constante de absorcao SC
        k12  = constante de distribuicao central -> periferico
        k21  = constante de redistribuicao periferico -> central
        ke   = constante de eliminacao do central

    Integracao por Euler explicito (suficiente para dt = 0.5h).

    Args:
        dose_mg: Dose absoluta em mg.
        bioavailability: Fracao absorvida (0-1).
        body_weight_kg: Peso do animal.
        duration_hours: Duracao da simulacao em horas.
        dt_hours: Passo temporal em horas.

    Returns:
        Lista de pontos temporais com concentracoes.
    """
    # Parametros PK
    ka = KA_SC_MEAN
    k12 = 0.50
    k21 = 0.01
    vc_l = 0.10 * body_weight_kg  # volume central em L
    cl_l_h = (CL_TOTAL_ML_MIN_KG * body_weight_kg * 60.0) / 1000.0
    ke = cl_l_h / vc_l
    vt_l = (VD_SS_L_PER_KG - 0.10) * body_weight_kg  # volume periferico em L

    # Condicoes iniciais (quantidades em mg)
    a_sc = dose_mg * bioavailability  # deposito SC
    a_c = 0.0   # compartimento central
    a_t = 0.0   # compartimento periferico
    a_eliminated = 0.0

    n_steps = int(duration_hours / dt_hours)
    results: list[PKTimePoint] = []

    for i in range(n_steps + 1):
        t = i * dt_hours

        # Concentracoes (ng/mL = ug/L)
        c_plasma = (a_c / vc_l) * 1e3 if vc_l > 0 else 0.0  # mg/L -> ug/L = ng/mL
        c_tissue = (a_t / vt_l) * 1e3 if vt_l > 0 else 0.0

        # Registrar ponto a cada hora inteira (ou primeiro e ultimo)
        if i == 0 or abs(t - round(t)) < dt_hours / 2 and t == int(t):
            results.append(PKTimePoint(
                time_hours=round(t, 1),
                plasma_conc_ng_ml=round(max(0.0, c_plasma), 2),
                tissue_conc_ng_ml=round(max(0.0, c_tissue), 2),
                amount_absorbed_mg=round(dose_mg * bioavailability - a_sc, 3),
                amount_eliminated_mg=round(a_eliminated, 3),
            ))

        # Euler step
        da_sc = -ka * a_sc
        da_c = ka * a_sc - (k12 + ke) * a_c + k21 * a_t
        da_t = k12 * a_c - k21 * a_t
        da_elim = ke * a_c

        a_sc += da_sc * dt_hours
        a_c += da_c * dt_hours
        a_t += da_t * dt_hours
        a_eliminated += da_elim * dt_hours

        # Protecao contra valores negativos (artefato numerico)
        a_sc = max(0.0, a_sc)
        a_c = max(0.0, a_c)
        a_t = max(0.0, a_t)

    return results


def compute_dosing_regimen(
    body_weight_kg: float = BODY_WEIGHT_KG,
) -> DosingRegimen:
    """Calcula regime posologico recomendado para caes com leishmaniose.

    Baseado em:
    1. Escalonamento alometrico de mipomersen (200 mg/semana em humanos de 70 kg)
       Dose canina = dose_humana * (BW_cao / BW_humano)^0.75 / BW_cao
    2. Dados de PK de PS-ASOs em caes (Geary RS et al. 2015)
    3. Consideracoes especificas de leishmaniose:
       - Carga parasitaria alta no inicio -> loading dose necessaria
       - Tratamento longo (semanas a meses) para clearance parasitario
       - Monitoramento de parasitemia guia duracao

    Comparacao com miltefosine (tratamento padrao canino):
    - Miltefosine: 2 mg/kg/dia VO por 28 dias
    - Problemas: teratogenicidade, GI toxicity, resistencia emergente
    - Vantagem do ASO: mecanismo diferente, sem resistencia cruzada

    Args:
        body_weight_kg: Peso do animal em kg.

    Returns:
        Regime posologico completo com comparacao.
    """
    # Escalonamento alometrico de mipomersen
    # Mipomersen humano: 200 mg/semana SC em adulto ~70 kg = ~2.86 mg/kg/semana
    # Escalonamento para cao de 15 kg:
    # Dose_cao (mg/kg) = Dose_humana (mg/kg) * (BW_humano / BW_cao)^0.25
    # Fator: (70/15)^0.25 = 1.47
    # Dose_cao = 2.86 * 1.47 = 4.2 mg/kg/semana
    # Arredondamento clinico: 5 mg/kg/semana
    human_dose_mg_kg = 200.0 / 70.0
    allometric_factor = (70.0 / body_weight_kg) ** 0.25
    scaled_dose_mg_kg = human_dose_mg_kg * allometric_factor

    # Arredondamento clinico para dose praticavel
    maintenance_dose_mg_kg = 5.0
    maintenance_dose_mg = maintenance_dose_mg_kg * body_weight_kg

    # Loading dose: 2x manutencao por 2 semanas (saturar tecidos rapidamente)
    loading_dose_mg_kg = 10.0
    loading_dose_mg = loading_dose_mg_kg * body_weight_kg

    # Volumes de injecao
    inj_vol_loading = loading_dose_mg / FORMULATION_CONC_MG_ML
    inj_vol_maint = maintenance_dose_mg / FORMULATION_CONC_MG_ML

    # Duracao do tratamento
    # Leishmaniose visceral canina: clearance parasitario lento
    # Tipicamente 3-6 meses de tratamento para reduzir carga parasitaria
    loading_weeks = 2
    maintenance_weeks = 10
    total_weeks = loading_weeks + maintenance_weeks

    # Concentracao trough estimada no steady state
    # Ctrough ~ (F * Dose) / (CL * tau) onde tau = intervalo entre doses
    tau_hours = 7.0 * 24.0  # 1 semana = 168 horas
    cl_ml_h = CL_TOTAL_ML_MIN_KG * body_weight_kg * 60.0
    # Fator de acumulacao: 1 / (1 - exp(-ke * tau))
    ke_terminal = math.log(2) / (TISSUE_HALF_LIFE_DAYS_MEAN * 24.0)
    accumulation = 1.0 / (1.0 - math.exp(-ke_terminal * tau_hours))
    ctrough_ng_ml = (
        BIOAVAILABILITY_SC_MEAN * maintenance_dose_mg * 1e6
    ) / (cl_ml_h * tau_hours) * accumulation

    # Indice terapeutico estimado
    # NOAEL em caes para PS-ASOs: ~40-50 mg/kg/semana (dados de Ionis Pharmaceuticals)
    # Dose terapeutica: 5 mg/kg/semana
    # TI = NOAEL / dose_terapeutica
    noael_mg_kg = 40.0
    therapeutic_index = noael_mg_kg / maintenance_dose_mg_kg

    # Comparacao com miltefosine
    comparison: dict[str, str] = {
        "miltefosine_dose": "2 mg/kg/day PO for 28 days",
        "miltefosine_cost": "Moderate (oral formulation)",
        "miltefosine_efficacy": "~70-80% clinical improvement; relapse common",
        "miltefosine_toxicity": "GI toxicity, teratogenicity, nephrotoxicity",
        "miltefosine_resistance": "Emerging resistance documented in endemic areas",
        "mrl_aso_001_dose": f"{maintenance_dose_mg_kg} mg/kg/week SC for {total_weeks} weeks",
        "mrl_aso_001_cost": "Higher (oligonucleotide synthesis)",
        "mrl_aso_001_efficacy": "Predicted: dual action (SL RNA knockdown + TLR9 immune activation)",
        "mrl_aso_001_toxicity": "Expected class effects: injection site reactions, mild thrombocytopenia",
        "mrl_aso_001_resistance": "Extremely unlikely — SL RNA is essential and 100% conserved",
        "key_advantage": (
            "MRL-ASO-001 targets a pan-trypanosomatid essential RNA (SL RNA) that is "
            "100% conserved and cannot mutate without loss of viability. Combined with "
            "TLR9-mediated immune activation, this provides a dual mechanism of action "
            "unavailable to any current leishmaniasis treatment."
        ),
    }

    return DosingRegimen(
        loading_dose_mg_kg=loading_dose_mg_kg,
        loading_dose_mg=loading_dose_mg,
        loading_frequency="twice weekly SC",
        loading_duration_weeks=loading_weeks,
        maintenance_dose_mg_kg=maintenance_dose_mg_kg,
        maintenance_dose_mg=maintenance_dose_mg,
        maintenance_frequency="once weekly SC",
        treatment_duration_weeks=maintenance_weeks,
        total_treatment_weeks=total_weeks,
        injection_volume_loading_ml=round(inj_vol_loading, 2),
        injection_volume_maintenance_ml=round(inj_vol_maint, 2),
        estimated_trough_ng_ml=round(ctrough_ng_ml, 1),
        therapeutic_index=round(therapeutic_index, 1),
        comparison_miltefosine=comparison,
    )
