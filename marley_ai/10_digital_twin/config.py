"""Configuracao do modulo 10_digital_twin — Digital twin de cao infectado sob tratamento.

Parametros para simulacao integrada PBPK + imune + estocastica de um cao
de 15 kg com leishmaniose visceral recebendo MRL-ASO-001 por 28 dias.

Os parametros de PK sao importados de aso_delivery.module_e_admet.pharmacokinetics
e os parametros imunes de aso_delivery.module_f_immune_sde.ode_model,
evitando duplicacao. Aqui definimos apenas parametros especificos da
integracao (coupling) e da simulacao unificada.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Parametros do animal modelo infectado
# ---------------------------------------------------------------------------

# Cao medio com leishmaniose visceral
DOG_WEIGHT_KG: float = 15.0

# Modificacoes de orgaos por leishmaniose visceral
# Ref: Koutinas AF et al. (1999) Vet Immunol Immunopathol — canine VL pathology
# Esplenomegalia: volume esplenico aumenta 2-4x em caes com LVC
SPLEEN_VOLUME_FACTOR: float = 3.0

# Hepatomegalia: volume hepatico aumenta 1.3-1.7x
LIVER_VOLUME_FACTOR: float = 1.5

# Carga parasitaria inicial (amastigotas por macrofago, escala do orgao)
INITIAL_PARASITE_LOAD: float = 1e6


# ---------------------------------------------------------------------------
# Parametros da simulacao
# ---------------------------------------------------------------------------

SIMULATION_DAYS: int = 28
SIMULATION_HOURS: int = SIMULATION_DAYS * 24  # 672 horas
DT_HOURS: float = 0.1  # passo de integracao (6 minutos)

# Numero de simulacoes Monte Carlo para camada SDE
N_MONTE_CARLO: int = 100

# Semente para reproducibilidade
SEED: int = 42


# ---------------------------------------------------------------------------
# Protocolo de dosagem (derivado de module_e)
# ---------------------------------------------------------------------------

# Loading dose: 10 mg/kg SC, 2x por semana, por 2 semanas
LOADING_DOSE_MG_KG: float = 10.0
LOADING_FREQUENCY_HOURS: float = 84.0  # ~3.5 dias (2x/semana)
LOADING_DURATION_WEEKS: int = 2

# Manutencao: 5 mg/kg SC, 1x por semana
MAINTENANCE_DOSE_MG_KG: float = 5.0
MAINTENANCE_FREQUENCY_HOURS: float = 168.0  # 7 dias (1x/semana)


# ---------------------------------------------------------------------------
# Parametros PBPK (4 compartimentos)
# ---------------------------------------------------------------------------

# Volumes de compartimento para cao de 15 kg (L)
# Ref: Davies B, Morris T (1993) Pharm Res — volumes de orgaos caninos
# Volume central (plasma): ~80 mL/kg * 0.55 (hematocrito) -> ~0.66 L plasma
PLASMA_VOLUME_L: float = 0.10 * DOG_WEIGHT_KG  # Vc = 0.1 L/kg (modulo E)

# Figado: ~32 g/kg em caes -> ~480 g -> ~0.48 L (densidade ~1.0)
# Com hepatomegalia 1.5x -> ~0.72 L
LIVER_VOLUME_L: float = 0.032 * DOG_WEIGHT_KG * LIVER_VOLUME_FACTOR

# Baco: ~3.5 g/kg em caes saudaveis -> ~52.5 g -> ~0.053 L
# Com esplenomegalia 3x -> ~0.16 L
SPLEEN_VOLUME_L: float = 0.0035 * DOG_WEIGHT_KG * SPLEEN_VOLUME_FACTOR

# Rim: ~8 g/kg -> ~120 g -> ~0.12 L (ambos rins)
KIDNEY_VOLUME_L: float = 0.008 * DOG_WEIGHT_KG

# Fluxo sanguineo para cada orgao (L/h)
# Ref: Davies B, Morris T (1993) — cardiac output em caes ~5.5 L/min/m^2
# Convertido para fracao do debito cardiaco
CARDIAC_OUTPUT_L_H: float = 5.5 * 60.0 * (DOG_WEIGHT_KG / 70.0) ** 0.75

# Fracoes do debito cardiaco para cada orgao (cao)
LIVER_BLOOD_FLOW_FRACTION: float = 0.25   # 25% (portal + arterial hepatico)
SPLEEN_BLOOD_FLOW_FRACTION: float = 0.04  # 4% (arteria esplenica)
KIDNEY_BLOOD_FLOW_FRACTION: float = 0.20  # 20% (arterias renais)


# ---------------------------------------------------------------------------
# Parametros de coupling PK -> imune
# ---------------------------------------------------------------------------

# Fator de conversao: concentracao tecidual (ng/mL) -> concentracao efetiva (uM)
# MW do ASO = 8500 Da = 8500 g/mol
# 1 ng/mL = 1 ug/L = (1e-6 g/L) / (8500 g/mol) = 1.176e-10 mol/L = 0.0001176 uM
ASO_MW_DA: float = 8500.0
NG_ML_TO_UM: float = 1e-3 / ASO_MW_DA * 1e6  # ng/mL -> uM (~0.0001176)


# ---------------------------------------------------------------------------
# Intensidade do ruido estocastico
# ---------------------------------------------------------------------------

# Sigma para Euler-Maruyama (variabilidade entre caes individuais)
# 0.15 = 15% variabilidade — maior que intra-celular (0.10) do modulo F
# porque inclui variabilidade farmacocinetica entre individuos
SDE_SIGMA: float = 0.15


# ---------------------------------------------------------------------------
# Dataclass de configuracao
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DigitalTwinConfig(AIModuleConfig):
    """Configuracao especifica para o digital twin de tratamento ASO.

    Integra tres camadas de simulacao:
    1. PBPK: distribuicao do farmaco em 4 compartimentos
    2. ODE imune: resposta de macrofagos ao ASO (deterministico)
    3. SDE estocastico: variabilidade entre caes individuais (Monte Carlo)

    A simulacao responde: se injetarmos MRL-ASO-001 em um cao infectado,
    o que acontece hora-a-hora por 28 dias?
    """

    # --- Simulacao ---
    simulation_days: int = SIMULATION_DAYS
    simulation_hours: int = SIMULATION_HOURS
    dt_hours: float = DT_HOURS
    n_monte_carlo: int = N_MONTE_CARLO
    seed: int = SEED

    # --- Animal ---
    dog_weight_kg: float = DOG_WEIGHT_KG
    initial_parasite_load: float = INITIAL_PARASITE_LOAD

    # --- Dosagem ---
    loading_dose_mg_kg: float = LOADING_DOSE_MG_KG
    maintenance_dose_mg_kg: float = MAINTENANCE_DOSE_MG_KG

    # --- Ruido SDE ---
    sde_sigma: float = SDE_SIGMA

    # --- Caminhos ---
    simulations_dir: Path = field(
        default_factory=lambda: AI_ROOT / "10_digital_twin" / "simulations"
    )
