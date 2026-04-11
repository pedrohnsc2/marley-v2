"""Modelo PBPK (Physiologically-Based Pharmacokinetic) de 4 compartimentos.

Simula a distribuicao de MRL-ASO-001 em um cao de 15 kg com leishmaniose
visceral apos injecao subcutanea. O modelo captura:

Compartimentos:
    1. Plasma (central) — absorve a dose SC, distribui para orgaos
    2. Figado — hepatomegalia 1.5x, sinusoides fenestrados, Kp=30
    3. Baco — esplenomegalia 3x, principal reservatorio do parasita, Kp=15
    4. Rim — cortex renal, Kp=25, excrecao de metabolitos

O modelo e derivado dos parametros de module_e_admet.pharmacokinetics
(constantes de absorcao, distribuicao, eliminacao) sem duplicacao.

Leishmaniose modifica a PK:
- Esplenomegalia aumenta volume do baco 3x -> mais farmaco acumula
- Hepatomegalia aumenta volume do figado 1.5x -> sequestro hepatico
- Macrofagos infectados captam mais ASO (receptores scavenger upregulated)

Integracao: RK4 (Runge-Kutta 4a ordem) implementado com numpy.

Refs:
- Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51 — PK de PS-ASOs
- Davies B, Morris T (1993) Pharm Res — volumes de orgaos caninos
- Koutinas AF et al. (1999) Vet Immunol Immunopathol — patologia canine VL
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

# Importar parametros PK do modulo E (sem duplicar)
from aso_delivery.module_e_admet.pharmacokinetics import (
    BIOAVAILABILITY_SC_MEAN,
    KA_SC_MEAN,
    CL_TOTAL_ML_MIN_KG,
    TISSUE_PARTITION_COEFFICIENTS,
)

from marley_ai.config import AI_ROOT as _AI_ROOT  # noqa: F401 — caminho base

# Importar config do digital twin via importacao relativa
# (nome do pacote comeca com digito, nao pode ser importado diretamente)
from . config import (
    CARDIAC_OUTPUT_L_H,
    DOG_WEIGHT_KG,
    KIDNEY_BLOOD_FLOW_FRACTION,
    KIDNEY_VOLUME_L,
    LIVER_BLOOD_FLOW_FRACTION,
    LIVER_VOLUME_L,
    LOADING_DOSE_MG_KG,
    LOADING_DURATION_WEEKS,
    LOADING_FREQUENCY_HOURS,
    MAINTENANCE_DOSE_MG_KG,
    MAINTENANCE_FREQUENCY_HOURS,
    PLASMA_VOLUME_L,
    SIMULATION_HOURS,
    SPLEEN_BLOOD_FLOW_FRACTION,
    SPLEEN_VOLUME_L,
)


# ---------------------------------------------------------------------------
# Indices das variaveis de estado PBPK
# ---------------------------------------------------------------------------

IDX_SC: Final[int] = 0       # deposito subcutaneo (mg)
IDX_PLASMA: Final[int] = 1   # plasma (mg)
IDX_LIVER: Final[int] = 2    # figado (mg)
IDX_SPLEEN: Final[int] = 3   # baco (mg)
IDX_KIDNEY: Final[int] = 4   # rim (mg)
N_PK_VARS: Final[int] = 5


# ---------------------------------------------------------------------------
# Parametros PBPK
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PBPKParams:
    """Parametros do modelo PBPK de 4 compartimentos para PS-ASOs em caes.

    Volumes em litros, fluxos em L/h, constantes em 1/h.
    Coeficientes de particionamento (Kp) do modulo E.
    """

    # Peso do animal
    body_weight_kg: float = DOG_WEIGHT_KG

    # Absorcao SC
    ka: float = KA_SC_MEAN                        # constante de absorcao (1/h)
    bioavailability: float = BIOAVAILABILITY_SC_MEAN  # fracao absorvida

    # Volumes de compartimento (L)
    v_plasma: float = PLASMA_VOLUME_L
    v_liver: float = LIVER_VOLUME_L
    v_spleen: float = SPLEEN_VOLUME_L
    v_kidney: float = KIDNEY_VOLUME_L

    # Fluxo sanguineo para cada orgao (L/h)
    q_liver: float = CARDIAC_OUTPUT_L_H * LIVER_BLOOD_FLOW_FRACTION
    q_spleen: float = CARDIAC_OUTPUT_L_H * SPLEEN_BLOOD_FLOW_FRACTION
    q_kidney: float = CARDIAC_OUTPUT_L_H * KIDNEY_BLOOD_FLOW_FRACTION

    # Coeficientes de particionamento tecido:plasma (do modulo E)
    kp_liver: float = TISSUE_PARTITION_COEFFICIENTS["liver"]
    kp_spleen: float = TISSUE_PARTITION_COEFFICIENTS["spleen"]
    kp_kidney: float = TISSUE_PARTITION_COEFFICIENTS["kidney_cortex"]

    # Eliminacao
    # CL total em L/h = CL_TOTAL_ML_MIN_KG * BW * 60 / 1000
    cl_total_l_h: float = CL_TOTAL_ML_MIN_KG * DOG_WEIGHT_KG * 60.0 / 1000.0

    # Fracao de clearance renal vs hepatico
    cl_renal_fraction: float = 0.70
    cl_hepatic_fraction: float = 0.30


# ---------------------------------------------------------------------------
# Protocolo de dosagem
# ---------------------------------------------------------------------------

def generate_dose_schedule(
    duration_hours: float = SIMULATION_HOURS,
    body_weight_kg: float = DOG_WEIGHT_KG,
) -> list[tuple[float, float]]:
    """Gera cronograma de doses (hora, dose_mg) para o protocolo completo.

    Protocolo MRL-ASO-001 (derivado de module_e):
    - Semanas 1-2: loading dose 10 mg/kg SC, 2x/semana (dias 0, 3.5, 7, 10.5)
    - Semanas 3+: manutencao 5 mg/kg SC, 1x/semana (dias 14, 21, 28)

    Args:
        duration_hours: Duracao total da simulacao em horas.
        body_weight_kg: Peso do animal em kg.

    Returns:
        Lista de tuplas (hora_da_dose, dose_em_mg).
    """
    schedule: list[tuple[float, float]] = []
    loading_dose_mg = LOADING_DOSE_MG_KG * body_weight_kg
    maintenance_dose_mg = MAINTENANCE_DOSE_MG_KG * body_weight_kg
    loading_end_hours = LOADING_DURATION_WEEKS * 7 * 24  # 336 horas

    # Fase de loading: 2x por semana por 2 semanas
    t = 0.0
    while t < loading_end_hours and t < duration_hours:
        schedule.append((t, loading_dose_mg))
        t += LOADING_FREQUENCY_HOURS

    # Fase de manutencao: 1x por semana
    t = loading_end_hours
    while t < duration_hours:
        schedule.append((t, maintenance_dose_mg))
        t += MAINTENANCE_FREQUENCY_HOURS

    return schedule


# ---------------------------------------------------------------------------
# Sistema de ODEs PBPK
# ---------------------------------------------------------------------------

def pbpk_rhs(
    y: NDArray[np.float64],
    params: PBPKParams,
) -> NDArray[np.float64]:
    """Calcula o lado direito do sistema PBPK.

    Modelo de perfusao limitada (flow-limited) para cada orgao:
        dA_organ/dt = Q_organ * (C_plasma - C_organ/Kp) - CL_organ * C_organ/Kp

    Onde:
        A = quantidade de farmaco no orgao (mg)
        Q = fluxo sanguineo (L/h)
        C_plasma = A_plasma / V_plasma (mg/L)
        C_organ = A_organ / V_organ (mg/L)
        Kp = coeficiente de particionamento tecido:plasma
        CL_organ = clearance no orgao (L/h)

    Args:
        y: Vetor de estado [A_sc, A_plasma, A_liver, A_spleen, A_kidney] em mg.
        params: Parametros PBPK.

    Returns:
        Vetor de derivadas dy/dt (5 componentes).
    """
    # Extrair quantidades (com floor de zero)
    a_sc = max(y[IDX_SC], 0.0)
    a_plasma = max(y[IDX_PLASMA], 0.0)
    a_liver = max(y[IDX_LIVER], 0.0)
    a_spleen = max(y[IDX_SPLEEN], 0.0)
    a_kidney = max(y[IDX_KIDNEY], 0.0)

    # Concentracoes (mg/L)
    c_plasma = a_plasma / params.v_plasma if params.v_plasma > 0 else 0.0
    c_liver = a_liver / params.v_liver if params.v_liver > 0 else 0.0
    c_spleen = a_spleen / params.v_spleen if params.v_spleen > 0 else 0.0
    c_kidney = a_kidney / params.v_kidney if params.v_kidney > 0 else 0.0

    # Concentracoes livres no tecido (para retorno ao plasma)
    # C_livre = C_organ / Kp
    c_liver_free = c_liver / params.kp_liver
    c_spleen_free = c_spleen / params.kp_spleen
    c_kidney_free = c_kidney / params.kp_kidney

    # Clearance por orgao (L/h)
    cl_hepatic = params.cl_total_l_h * params.cl_hepatic_fraction
    cl_renal = params.cl_total_l_h * params.cl_renal_fraction

    dydt = np.zeros(N_PK_VARS, dtype=np.float64)

    # dA_sc/dt = -ka * A_sc (absorcao SC de primeira ordem)
    dydt[IDX_SC] = -params.ka * a_sc

    # dA_plasma/dt = ka*A_sc - Q_total_saida*C_plasma + retornos dos orgaos
    # Entrada: absorcao SC + retorno dos orgaos
    # Saida: distribuicao para orgaos
    influx_from_organs = (
        params.q_liver * c_liver_free
        + params.q_spleen * c_spleen_free
        + params.q_kidney * c_kidney_free
    )
    efflux_to_organs = (
        params.q_liver + params.q_spleen + params.q_kidney
    ) * c_plasma

    dydt[IDX_PLASMA] = (
        params.ka * a_sc
        + influx_from_organs
        - efflux_to_organs
    )

    # dA_liver/dt = Q_liver * C_plasma - Q_liver * C_liver/Kp - CL_hepatic * C_liver/Kp
    dydt[IDX_LIVER] = (
        params.q_liver * c_plasma
        - params.q_liver * c_liver_free
        - cl_hepatic * c_liver_free
    )

    # dA_spleen/dt = Q_spleen * C_plasma - Q_spleen * C_spleen/Kp
    # Baco nao tem clearance proprio — farmaco retorna ao plasma
    dydt[IDX_SPLEEN] = (
        params.q_spleen * c_plasma
        - params.q_spleen * c_spleen_free
    )

    # dA_kidney/dt = Q_kidney * C_plasma - Q_kidney * C_kidney/Kp - CL_renal * C_kidney/Kp
    dydt[IDX_KIDNEY] = (
        params.q_kidney * c_plasma
        - params.q_kidney * c_kidney_free
        - cl_renal * c_kidney_free
    )

    return dydt


# ---------------------------------------------------------------------------
# Integrador RK4 com dosagem multipla
# ---------------------------------------------------------------------------

def _rk4_step(
    y: NDArray[np.float64],
    params: PBPKParams,
    dt: float,
) -> NDArray[np.float64]:
    """Executa um passo RK4 (Runge-Kutta 4a ordem).

    Args:
        y: Estado atual.
        params: Parametros PBPK.
        dt: Passo temporal (horas).

    Returns:
        Novo estado apos o passo.
    """
    k1 = pbpk_rhs(y, params)
    k2 = pbpk_rhs(y + 0.5 * dt * k1, params)
    k3 = pbpk_rhs(y + 0.5 * dt * k2, params)
    k4 = pbpk_rhs(y + dt * k3, params)

    y_new = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    # Nao-negatividade (concentracoes biologicas)
    np.maximum(y_new, 0.0, out=y_new)
    return y_new


def integrate_pbpk(
    params: PBPKParams | None = None,
    duration_hours: float = SIMULATION_HOURS,
    dt: float = 0.1,
    dose_schedule: list[tuple[float, float]] | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integra o modelo PBPK com dosagem multipla.

    A cada dose no cronograma, adiciona a quantidade de farmaco ao
    deposito SC (multiplicada pela biodisponibilidade).

    Args:
        params: Parametros PBPK. Se None, usa defaults.
        duration_hours: Duracao total em horas.
        dt: Passo temporal (horas).
        dose_schedule: Lista de (hora, dose_mg). Se None, gera automaticamente.

    Returns:
        Tupla (time_array, solution) onde solution tem shape (n_steps+1, 5).
        Colunas: [SC, plasma, liver, spleen, kidney] em mg.
    """
    if params is None:
        params = PBPKParams()
    if dose_schedule is None:
        dose_schedule = generate_dose_schedule(
            duration_hours=duration_hours,
            body_weight_kg=params.body_weight_kg,
        )

    n_steps = int(duration_hours / dt)
    t = np.linspace(0.0, duration_hours, n_steps + 1)
    sol = np.zeros((n_steps + 1, N_PK_VARS), dtype=np.float64)

    # Condicoes iniciais: tudo zero (nenhum farmaco no corpo)
    # A primeira dose sera adicionada no loop

    # Converter dose_schedule para dict de indices temporais
    # Para cada dose, encontrar o indice mais proximo
    dose_indices: dict[int, float] = {}
    for dose_time, dose_mg in dose_schedule:
        idx = int(round(dose_time / dt))
        idx = min(idx, n_steps)
        # Acumular se multiplas doses no mesmo timestep
        dose_indices[idx] = dose_indices.get(idx, 0.0) + dose_mg

    for i in range(n_steps):
        # Aplicar dose no deposito SC (se houver)
        if i in dose_indices:
            sol[i, IDX_SC] += dose_indices[i] * params.bioavailability

        # Passo RK4
        sol[i + 1] = _rk4_step(sol[i], params, dt)

    # Aplicar dose final se existir
    if n_steps in dose_indices:
        sol[n_steps, IDX_SC] += dose_indices[n_steps] * params.bioavailability

    return t, sol


# ---------------------------------------------------------------------------
# Funcoes auxiliares para extrair concentracoes
# ---------------------------------------------------------------------------

def amounts_to_concentrations(
    sol: NDArray[np.float64],
    params: PBPKParams | None = None,
) -> NDArray[np.float64]:
    """Converte quantidades (mg) para concentracoes (ng/mL = ug/L).

    Args:
        sol: Array de quantidades, shape (n_steps+1, 5).
        params: Parametros PBPK (para volumes).

    Returns:
        Array de concentracoes (ng/mL), mesmo shape.
    """
    if params is None:
        params = PBPKParams()

    conc = np.zeros_like(sol)
    # SC nao tem concentracao significativa (deposito)
    conc[:, IDX_SC] = sol[:, IDX_SC]  # mantemos em mg para referencia

    # Concentracao = (quantidade em mg) / (volume em L) * 1000 (-> ng/mL = ug/L)
    # Na verdade: mg/L = ug/mL = 1000 ng/mL
    # Entao: C(ng/mL) = A(mg) / V(L) * 1e3
    conc[:, IDX_PLASMA] = sol[:, IDX_PLASMA] / params.v_plasma * 1e3
    conc[:, IDX_LIVER] = sol[:, IDX_LIVER] / params.v_liver * 1e3
    conc[:, IDX_SPLEEN] = sol[:, IDX_SPLEEN] / params.v_spleen * 1e3
    conc[:, IDX_KIDNEY] = sol[:, IDX_KIDNEY] / params.v_kidney * 1e3

    return conc


def get_spleen_concentration_uM(
    sol: NDArray[np.float64],
    params: PBPKParams | None = None,
) -> NDArray[np.float64]:
    """Extrai concentracao esplenica em microM (para coupling com imune).

    Converte ng/mL -> uM usando peso molecular do ASO.
    1 ng/mL = 1e-9 g/mL = 1e-6 g/L
    C(uM) = C(g/L) / MW(g/mol) * 1e6
    C(uM) = C(ng/mL) * 1e-6 / 8500 * 1e6 = C(ng/mL) / 8500

    Args:
        sol: Array de quantidades, shape (n_steps+1, 5).
        params: Parametros PBPK.

    Returns:
        Array 1D de concentracoes esplenicas em uM.
    """
    if params is None:
        params = PBPKParams()

    from marley_ai.ten_digital_twin.config import ASO_MW_DA

    c_spleen_ng_ml = sol[:, IDX_SPLEEN] / params.v_spleen * 1e3
    c_spleen_uM = c_spleen_ng_ml / ASO_MW_DA
    return c_spleen_uM
