"""Modelo ODE deterministico de macrofago infectado sob tratamento com MRL-ASO-001.

Simula a dinamica intracelular de um macrofago canino infectado por
amastigotas de L. infantum, exposto ao ASO MRL-ASO-001. O modelo
captura dois mecanismos de acao simultaneos:

1. Antisense: bloqueio do trans-splicing do SL RNA, impedindo
   maturacao de mRNAs do parasita (k_aso)
2. Imunoestimulatorio: ativacao de TLR9 por motivos CpG-like
   no backbone PS, induzindo IFN-gamma -> iNOS -> NO (k_tlr9)

Estado do sistema (5 variaveis):
    P(t) — carga parasitaria intracelular (amastigotas)
    A(t) — concentracao de ASO ativo no fagolisossomo
    I(t) — nivel de IFN-gamma (citocina Th1)
    T(t) — nivel de TNF-alpha (citocina pro-inflamatoria)
    N(t) — nivel de oxido nitrico (NO, leishmanicida)

Referencias:
- Liew FY et al. (1990) J Exp Med 172(5):1557-1559 — NO mata amastigotas
- Klinman DM (2004) Nat Rev Immunol 4(4):249-258 — CpG e TLR9
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO PS e TLR9
- Gantt KR et al. (2001) J Immunol 167(2):893-901 — IFN-gamma e iNOS
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Indices das variaveis de estado (para legibilidade)
# ---------------------------------------------------------------------------

IDX_P: Final[int] = 0   # parasitas
IDX_A: Final[int] = 1   # ASO ativo
IDX_I: Final[int] = 2   # IFN-gamma
IDX_T: Final[int] = 3   # TNF-alpha
IDX_N: Final[int] = 4   # NO
N_VARS: Final[int] = 5


# ---------------------------------------------------------------------------
# Parametros do modelo
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MacrophageParams:
    """Parametros cineticos do modelo de macrofago infectado.

    Todos os parametros sao calibrados a partir da literatura e dos
    resultados de modulos anteriores do pipeline (Module A para k_deg).

    Unidades: concentracoes em unidades arbitrarias normalizadas,
    taxas em 1/hora.
    """

    # --- Crescimento do parasita ---
    # Taxa intrinseca de replicacao de amastigotas
    # Ref: tempo de duplicacao ~35h -> r = ln(2)/35 ~ 0.02/h
    r_P: float = 0.02

    # Capacidade de carga por macrofago
    # Ref: macrofagos caninos suportam ~100-200 amastigotas antes de lise
    K_P: float = 200.0

    # --- Matanca do parasita ---
    # Taxa de matanca mediada por NO
    # Ref: Liew et al. (1990) — NO e o principal mecanismo leishmanicida
    k_kill: float = 0.1

    # Taxa de matanca direta por bloqueio de splicing (antisense)
    # Ref: sem trans-splicing, mRNAs do parasita nao sao processados
    k_aso: float = 0.05

    # --- Farmacocinetica do ASO ---
    # Taxa de captacao do ASO extracelular para o fagolisossomo
    # Ref: macrofagos fagocitam PS-ASOs via receptores scavenger
    k_uptake: float = 0.05

    # Concentracao extracelular de ASO (microM)
    A_ext: float = 10.0

    # Taxa de degradacao do ASO (Module A: t1/2 = 1083h)
    # k_deg = ln(2) / 1083 ~ 0.00064/h
    k_deg: float = math.log(2) / 1083.0

    # --- Ativacao imune via TLR9 ---
    # Taxa de inducao de IFN-gamma por TLR9 ativado pelo ASO
    # Ref: Klinman (2004) — CpG DNA ativa TLR9 em fagolisossomos
    k_tlr9: float = 0.03

    # Taxa de inducao de TNF-alpha por TLR9
    k_tnf: float = 0.02

    # Taxa de decaimento de IFN-gamma
    # Ref: meia-vida de citocinas intracelulares ~4-8h
    d_I: float = 0.1

    # Taxa de decaimento de TNF-alpha
    d_T: float = 0.12

    # --- Cascata NO ---
    # Taxa de producao de NO via IFN-gamma -> iNOS
    # Ref: Gantt et al. (2001) — IFN-gamma induz iNOS em macrofagos
    k_no: float = 0.08

    # Taxa de decaimento de NO (meia-vida ~5-10 segundos in vivo,
    # mas usamos meia-vida efetiva no contexto intracelular ~2-4h)
    d_N: float = 0.2


@dataclass
class SimulationConfig:
    """Configuracao temporal da simulacao."""

    # Duracao total em horas (7 dias)
    t_end: float = 168.0

    # Passo de tempo para integracao (horas)
    dt: float = 0.1

    # Condicoes iniciais
    P0: float = 30.0    # amastigotas por macrofago (10-50, media 30)
    A0: float = 0.0     # ASO comeca em zero (acumula por captacao)
    I0: float = 0.0     # IFN-gamma basal ~ 0
    T0: float = 0.0     # TNF-alpha basal ~ 0
    N0: float = 0.0     # NO basal ~ 0

    @property
    def y0(self) -> NDArray[np.float64]:
        """Vetor de condicoes iniciais."""
        return np.array([self.P0, self.A0, self.I0, self.T0, self.N0],
                        dtype=np.float64)

    @property
    def n_steps(self) -> int:
        """Numero de passos temporais."""
        return int(self.t_end / self.dt)

    @property
    def time_array(self) -> NDArray[np.float64]:
        """Vetor de tempos."""
        return np.linspace(0.0, self.t_end, self.n_steps + 1)


# ---------------------------------------------------------------------------
# Sistema de ODEs
# ---------------------------------------------------------------------------


def ode_rhs(
    y: NDArray[np.float64],
    params: MacrophageParams,
) -> NDArray[np.float64]:
    """Calcula o lado direito do sistema de ODEs.

    Recebe o vetor de estado y = [P, A, I, T, N] e retorna dy/dt.
    Aplica floor de zero para evitar valores negativos numericos.

    Args:
        y: Vetor de estado (5 componentes).
        params: Parametros cineticos do modelo.

    Returns:
        Vetor de derivadas dy/dt (5 componentes).
    """
    # Extrair variaveis (com floor de zero)
    P = max(y[IDX_P], 0.0)
    A = max(y[IDX_A], 0.0)
    I = max(y[IDX_I], 0.0)
    T = max(y[IDX_T], 0.0)
    N = max(y[IDX_N], 0.0)

    dydt = np.zeros(N_VARS, dtype=np.float64)

    # dP/dt: crescimento logistico - matanca por NO - matanca por ASO
    dydt[IDX_P] = (
        params.r_P * P * (1.0 - P / params.K_P)
        - params.k_kill * N * P
        - params.k_aso * A * P
    )

    # dA/dt: captacao do meio extracelular - degradacao
    dydt[IDX_A] = params.k_uptake * params.A_ext - params.k_deg * A

    # dI/dt: inducao por TLR9 - decaimento
    dydt[IDX_I] = params.k_tlr9 * A - params.d_I * I

    # dT/dt: inducao por TLR9 - decaimento
    dydt[IDX_T] = params.k_tnf * A - params.d_T * T

    # dN/dt: producao via IFN-gamma/iNOS - decaimento
    dydt[IDX_N] = params.k_no * I - params.d_N * N

    return dydt


# ---------------------------------------------------------------------------
# Integrador deterministico (Euler explicito para consistencia com SDE)
# ---------------------------------------------------------------------------


def integrate_ode(
    params: MacrophageParams,
    config: SimulationConfig,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integra o sistema de ODEs pelo metodo de Euler explicito.

    Usa Euler explicito (nao RK4) para manter consistencia com o
    integrador SDE (Euler-Maruyama), permitindo comparacao direta
    entre trajetorias deterministicas e estocasticas.

    Args:
        params: Parametros cineticos do modelo.
        config: Configuracao temporal da simulacao.

    Returns:
        Tupla (time_array, solution) onde solution tem shape (n_steps+1, 5).
    """
    t = config.time_array
    n_steps = config.n_steps
    dt = config.dt

    sol = np.zeros((n_steps + 1, N_VARS), dtype=np.float64)
    sol[0] = config.y0

    for i in range(n_steps):
        dydt = ode_rhs(sol[i], params)
        sol[i + 1] = sol[i] + dydt * dt
        # Impor nao-negatividade (concentracoes biologicas)
        np.maximum(sol[i + 1], 0.0, out=sol[i + 1])

    return t, sol


# ---------------------------------------------------------------------------
# Cenarios de simulacao
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ScenarioResult:
    """Resultado de um cenario de simulacao deterministica."""

    name: str
    description: str
    time: NDArray[np.float64]
    solution: NDArray[np.float64]      # shape (n_steps+1, 5)
    final_parasite_load: float
    initial_parasite_load: float
    parasite_reduction_pct: float
    clearance_achieved: bool           # P < 1 ao final
    time_to_90pct_clearance_hours: float  # tempo para 90% de reducao (-1 se nao atingido)

    def to_dict(self) -> dict:
        """Converte para dicionario serializavel (sem arrays numpy)."""
        return {
            "name": self.name,
            "description": self.description,
            "final_parasite_load": round(self.final_parasite_load, 4),
            "initial_parasite_load": round(self.initial_parasite_load, 4),
            "parasite_reduction_pct": round(self.parasite_reduction_pct, 2),
            "clearance_achieved": self.clearance_achieved,
            "time_to_90pct_clearance_hours": round(self.time_to_90pct_clearance_hours, 2),
        }


def compute_clearance_time(
    time: NDArray[np.float64],
    parasite_trajectory: NDArray[np.float64],
    threshold_fraction: float = 0.1,
) -> float:
    """Calcula o tempo para atingir um limiar de reducao parasitaria.

    Args:
        time: Vetor de tempos.
        parasite_trajectory: Vetor P(t).
        threshold_fraction: Fracao residual (0.1 = 90% de clearance).

    Returns:
        Tempo em horas para atingir o limiar, ou -1.0 se nao atingido.
    """
    P0 = parasite_trajectory[0]
    if P0 <= 0:
        return 0.0

    threshold = P0 * threshold_fraction
    below = np.where(parasite_trajectory <= threshold)[0]

    if len(below) == 0:
        return -1.0

    return float(time[below[0]])


def run_scenario(
    name: str,
    description: str,
    params: MacrophageParams,
    config: SimulationConfig,
) -> ScenarioResult:
    """Executa um cenario de simulacao deterministica.

    Args:
        name: Identificador curto do cenario.
        description: Descricao legivel do cenario.
        params: Parametros (podem ser modificados para o cenario).
        config: Configuracao temporal.

    Returns:
        Resultado completo do cenario.
    """
    t, sol = integrate_ode(params, config)

    P = sol[:, IDX_P]
    P0 = P[0]
    P_final = P[-1]
    reduction = (1.0 - P_final / P0) * 100.0 if P0 > 0 else 0.0
    clearance_time = compute_clearance_time(t, P)

    return ScenarioResult(
        name=name,
        description=description,
        time=t,
        solution=sol,
        final_parasite_load=float(P_final),
        initial_parasite_load=float(P0),
        parasite_reduction_pct=float(max(reduction, 0.0)),
        clearance_achieved=bool(P_final < 1.0),
        time_to_90pct_clearance_hours=float(clearance_time),
    )
