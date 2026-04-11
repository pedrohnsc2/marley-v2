"""Extensao estocastica do modelo de macrofago infectado (SDE via Euler-Maruyama).

Adiciona ruido Wiener multiplicativo a cada equacao do sistema ODE,
modelando a variabilidade biologica inerente entre macrofagos individuais:
- Variacao na carga parasitaria inicial
- Heterogeneidade na captacao de ASO
- Diferencias na resposta imune inata entre celulas

O metodo de Euler-Maruyama e o analogo estocastico do metodo de Euler:
    X_{n+1} = X_n + f(X_n) * dt + sigma * X_n * dW_n
onde dW_n ~ N(0, dt) sao incrementos Wiener independentes.

Monte Carlo com N=1000 simulacoes gera:
- Trajetoria media +/- 95% CI para cada variavel
- Probabilidade de clearance parasitaria (P < 1)
- Distribuicao de tempos de clearance

Referencias:
- Kloeden PE, Platen E (1992) Numerical Solution of Stochastic
  Differential Equations. Springer. — Metodo Euler-Maruyama
- Higham DJ (2001) SIAM Review 43(3):525-546 — Tutorial SDE numerico
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

from aso_delivery.module_f_immune_sde.ode_model import (
    IDX_A,
    IDX_I,
    IDX_N,
    IDX_P,
    IDX_T,
    N_VARS,
    MacrophageParams,
    SimulationConfig,
    compute_clearance_time,
    ode_rhs,
)

# ---------------------------------------------------------------------------
# Constantes do modelo estocastico
# ---------------------------------------------------------------------------

# Intensidade do ruido multiplicativo (sigma)
# Valor de 0.1 = 10% de variabilidade relativa por passo,
# representando heterogeneidade celular tipica
DEFAULT_SIGMA: Final[float] = 0.1

# Numero de simulacoes Monte Carlo
DEFAULT_N_SIMS: Final[int] = 1000

# Semente para reprodutibilidade
DEFAULT_SEED: Final[int] = 42


# ---------------------------------------------------------------------------
# Resultado do ensemble estocastico
# ---------------------------------------------------------------------------


@dataclass
class SDEEnsembleResult:
    """Resultado agregado de multiplas simulacoes SDE.

    Armazena estatisticas do ensemble (media, percentis, probabilidades)
    sem manter todas as trajetorias individuais em memoria.
    """

    name: str
    n_simulations: int
    sigma: float
    time: NDArray[np.float64]                # (n_steps+1,)

    # Estatisticas por variavel: shape (n_steps+1, 5)
    mean: NDArray[np.float64]
    std: NDArray[np.float64]
    percentile_2_5: NDArray[np.float64]      # limite inferior 95% CI
    percentile_97_5: NDArray[np.float64]     # limite superior 95% CI
    median: NDArray[np.float64]

    # Metricas de clearance
    clearance_probability: float             # fracao de sims com P_final < 1
    mean_clearance_time_hours: float         # media do tempo de clearance (apenas sims com clearance)
    std_clearance_time_hours: float
    mean_final_parasite_load: float
    std_final_parasite_load: float
    mean_parasite_reduction_pct: float

    def to_dict(self) -> dict:
        """Converte para dicionario serializavel (sem arrays numpy)."""
        return {
            "name": self.name,
            "n_simulations": self.n_simulations,
            "sigma": self.sigma,
            "clearance_probability": round(self.clearance_probability, 4),
            "mean_clearance_time_hours": round(self.mean_clearance_time_hours, 2),
            "std_clearance_time_hours": round(self.std_clearance_time_hours, 2),
            "mean_final_parasite_load": round(self.mean_final_parasite_load, 4),
            "std_final_parasite_load": round(self.std_final_parasite_load, 4),
            "mean_parasite_reduction_pct": round(self.mean_parasite_reduction_pct, 2),
        }

    def time_series_dict(self, var_names: list[str] | None = None) -> dict:
        """Exporta series temporais medias +/- SD para JSON.

        Amostra a cada ~1 hora para manter o JSON manejavel
        (168 pontos em vez de 1681).

        Args:
            var_names: Nomes das variaveis (default: P, A, I, T, N).

        Returns:
            Dicionario com series temporais para cada variavel.
        """
        if var_names is None:
            var_names = ["P_parasites", "A_aso", "I_ifn_gamma",
                         "T_tnf_alpha", "N_nitric_oxide"]

        # Amostrar a cada ~1h (indice a cada 10 passos se dt=0.1)
        dt = self.time[1] - self.time[0] if len(self.time) > 1 else 0.1
        step = max(1, int(1.0 / dt))
        indices = list(range(0, len(self.time), step))
        # Garantir que o ultimo ponto esta incluido
        if indices[-1] != len(self.time) - 1:
            indices.append(len(self.time) - 1)

        series: dict[str, list[dict[str, float]]] = {}
        for var_idx, var_name in enumerate(var_names):
            points: list[dict[str, float]] = []
            for i in indices:
                points.append({
                    "time_hours": round(float(self.time[i]), 1),
                    "mean": round(float(self.mean[i, var_idx]), 6),
                    "std": round(float(self.std[i, var_idx]), 6),
                    "ci_lower": round(float(self.percentile_2_5[i, var_idx]), 6),
                    "ci_upper": round(float(self.percentile_97_5[i, var_idx]), 6),
                    "median": round(float(self.median[i, var_idx]), 6),
                })
            series[var_name] = points

        return series


# ---------------------------------------------------------------------------
# Integrador Euler-Maruyama
# ---------------------------------------------------------------------------


def euler_maruyama_step(
    y: NDArray[np.float64],
    params: MacrophageParams,
    dt: float,
    sigma: float,
    dW: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Executa um passo do metodo Euler-Maruyama.

    X_{n+1} = X_n + f(X_n)*dt + sigma * X_n * dW_n

    O ruido multiplicativo (proporcional ao estado) modela variabilidade
    biologica: celulas com mais parasitas tem maior incerteza absoluta.

    Args:
        y: Estado atual (5 componentes).
        params: Parametros cineticos.
        dt: Passo temporal.
        sigma: Intensidade do ruido.
        dW: Incremento Wiener (5 componentes), pré-gerado.

    Returns:
        Novo estado apos o passo.
    """
    drift = ode_rhs(y, params)
    diffusion = sigma * y * dW

    y_new = y + drift * dt + diffusion
    # Impor nao-negatividade
    np.maximum(y_new, 0.0, out=y_new)

    return y_new


def run_single_sde(
    params: MacrophageParams,
    config: SimulationConfig,
    sigma: float,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Executa uma unica trajetoria SDE.

    Args:
        params: Parametros cineticos.
        config: Configuracao temporal.
        sigma: Intensidade do ruido.
        rng: Gerador de numeros aleatorios.

    Returns:
        Array de shape (n_steps+1, 5) com a trajetoria completa.
    """
    n_steps = config.n_steps
    dt = config.dt
    sqrt_dt = np.sqrt(dt)

    sol = np.zeros((n_steps + 1, N_VARS), dtype=np.float64)
    sol[0] = config.y0

    for i in range(n_steps):
        # Incremento Wiener: dW ~ N(0, dt) = sqrt(dt) * N(0, 1)
        dW = sqrt_dt * rng.standard_normal(N_VARS)
        sol[i + 1] = euler_maruyama_step(sol[i], params, dt, sigma, dW)

    return sol


# ---------------------------------------------------------------------------
# Ensemble Monte Carlo
# ---------------------------------------------------------------------------


def run_sde_ensemble(
    name: str,
    params: MacrophageParams,
    config: SimulationConfig,
    n_sims: int = DEFAULT_N_SIMS,
    sigma: float = DEFAULT_SIGMA,
    seed: int = DEFAULT_SEED,
) -> SDEEnsembleResult:
    """Executa ensemble Monte Carlo de simulacoes SDE.

    Gera n_sims trajetorias independentes e calcula estatisticas
    agregadas (media, desvio padrao, percentis, probabilidade de clearance).

    Estrategia de memoria: acumula estatisticas incrementalmente
    para n_sims grande, mas para 1000 sims x 1681 passos x 5 vars
    o uso total e ~67 MB, aceitavel.

    Args:
        name: Identificador do cenario.
        params: Parametros cineticos.
        config: Configuracao temporal.
        n_sims: Numero de simulacoes Monte Carlo.
        sigma: Intensidade do ruido.
        seed: Semente para reprodutibilidade.

    Returns:
        Resultado agregado do ensemble.
    """
    rng = np.random.default_rng(seed)
    t = config.time_array
    n_steps = config.n_steps

    # Armazenar todas as trajetorias para calcular percentis
    # shape: (n_sims, n_steps+1, N_VARS)
    all_trajectories = np.zeros((n_sims, n_steps + 1, N_VARS), dtype=np.float64)

    clearance_times: list[float] = []
    final_loads: list[float] = []

    for sim_idx in range(n_sims):
        sol = run_single_sde(params, config, sigma, rng)
        all_trajectories[sim_idx] = sol

        # Metricas de clearance para esta simulacao
        P_traj = sol[:, IDX_P]
        P_final = P_traj[-1]
        final_loads.append(P_final)

        ct = compute_clearance_time(t, P_traj)
        if ct >= 0:
            clearance_times.append(ct)

    # Estatisticas do ensemble
    mean = np.mean(all_trajectories, axis=0)
    std = np.std(all_trajectories, axis=0)
    pct_2_5 = np.percentile(all_trajectories, 2.5, axis=0)
    pct_97_5 = np.percentile(all_trajectories, 97.5, axis=0)
    median = np.median(all_trajectories, axis=0)

    # Clearance
    final_loads_arr = np.array(final_loads)
    clearance_prob = float(np.mean(final_loads_arr < 1.0))

    if len(clearance_times) > 0:
        mean_ct = float(np.mean(clearance_times))
        std_ct = float(np.std(clearance_times))
    else:
        mean_ct = -1.0
        std_ct = 0.0

    P0 = config.P0
    mean_reduction = (1.0 - float(np.mean(final_loads_arr)) / P0) * 100.0 if P0 > 0 else 0.0

    return SDEEnsembleResult(
        name=name,
        n_simulations=n_sims,
        sigma=sigma,
        time=t,
        mean=mean,
        std=std,
        percentile_2_5=pct_2_5,
        percentile_97_5=pct_97_5,
        median=median,
        clearance_probability=clearance_prob,
        mean_clearance_time_hours=mean_ct,
        std_clearance_time_hours=std_ct,
        mean_final_parasite_load=float(np.mean(final_loads_arr)),
        std_final_parasite_load=float(np.std(final_loads_arr)),
        mean_parasite_reduction_pct=max(mean_reduction, 0.0),
    )


# ---------------------------------------------------------------------------
# Dose-resposta
# ---------------------------------------------------------------------------


@dataclass
class DoseResponsePoint:
    """Resultado de uma dose na curva dose-resposta."""

    dose_uM: float
    clearance_probability_72h: float
    mean_time_to_90pct_clearance_hours: float
    mean_final_parasite_load: float
    mean_parasite_reduction_pct: float

    def to_dict(self) -> dict:
        return {
            "dose_uM": round(self.dose_uM, 2),
            "clearance_probability_72h": round(self.clearance_probability_72h, 4),
            "mean_time_to_90pct_clearance_hours": round(self.mean_time_to_90pct_clearance_hours, 2),
            "mean_final_parasite_load": round(self.mean_final_parasite_load, 4),
            "mean_parasite_reduction_pct": round(self.mean_parasite_reduction_pct, 2),
        }


def run_dose_response(
    doses_uM: list[float],
    base_params: MacrophageParams,
    config: SimulationConfig,
    n_sims: int = DEFAULT_N_SIMS,
    sigma: float = DEFAULT_SIGMA,
    seed: int = DEFAULT_SEED,
) -> list[DoseResponsePoint]:
    """Executa curva dose-resposta completa.

    Para cada dose, executa ensemble SDE e calcula:
    - Probabilidade de clearance em 72h
    - Tempo medio para 90% de clearance
    - Carga parasitaria final media

    A avaliacao em 72h usa config com t_end=72h para clearance_probability,
    mas o ensemble padrao usa t_end completo para tempo de clearance.

    Args:
        doses_uM: Lista de doses extracelulares em microM.
        base_params: Parametros base (A_ext sera substituido).
        config: Configuracao temporal.
        n_sims: Numero de simulacoes por dose.
        sigma: Intensidade do ruido.
        seed: Semente base (incrementada por dose).

    Returns:
        Lista de resultados por dose.
    """
    results: list[DoseResponsePoint] = []

    # Config de 72h para avaliar clearance nesse horizonte
    config_72h = SimulationConfig(
        t_end=72.0,
        dt=config.dt,
        P0=config.P0,
        A0=config.A0,
        I0=config.I0,
        T0=config.T0,
        N0=config.N0,
    )

    for dose_idx, dose in enumerate(doses_uM):
        # Modificar A_ext para esta dose
        params = MacrophageParams(
            r_P=base_params.r_P,
            K_P=base_params.K_P,
            k_kill=base_params.k_kill,
            k_aso=base_params.k_aso,
            k_uptake=base_params.k_uptake,
            A_ext=dose,
            k_deg=base_params.k_deg,
            k_tlr9=base_params.k_tlr9,
            k_tnf=base_params.k_tnf,
            d_I=base_params.d_I,
            d_T=base_params.d_T,
            k_no=base_params.k_no,
            d_N=base_params.d_N,
        )

        # Ensemble para tempo completo (metricas gerais)
        ensemble_full = run_sde_ensemble(
            name=f"dose_{dose:.1f}uM",
            params=params,
            config=config,
            n_sims=n_sims,
            sigma=sigma,
            seed=seed + dose_idx,
        )

        # Ensemble para 72h (probabilidade de clearance nesse horizonte)
        ensemble_72h = run_sde_ensemble(
            name=f"dose_{dose:.1f}uM_72h",
            params=params,
            config=config_72h,
            n_sims=n_sims,
            sigma=sigma,
            seed=seed + dose_idx + 1000,
        )

        results.append(DoseResponsePoint(
            dose_uM=dose,
            clearance_probability_72h=ensemble_72h.clearance_probability,
            mean_time_to_90pct_clearance_hours=ensemble_full.mean_clearance_time_hours,
            mean_final_parasite_load=ensemble_full.mean_final_parasite_load,
            mean_parasite_reduction_pct=ensemble_full.mean_parasite_reduction_pct,
        ))

    return results


def estimate_ec_values(
    dose_response: list[DoseResponsePoint],
) -> dict[str, float]:
    """Estima EC50 e EC90 por interpolacao linear da curva dose-resposta.

    EC50/EC90: dose que produz 50%/90% de reducao parasitaria media.

    Args:
        dose_response: Resultados da curva dose-resposta.

    Returns:
        Dicionario com EC50 e EC90 (em microM), ou -1 se nao estimavel.
    """
    doses = [p.dose_uM for p in dose_response]
    reductions = [p.mean_parasite_reduction_pct for p in dose_response]

    ec50 = _interpolate_ec(doses, reductions, 50.0)
    ec90 = _interpolate_ec(doses, reductions, 90.0)

    return {
        "EC50_uM": round(ec50, 3) if ec50 >= 0 else -1.0,
        "EC90_uM": round(ec90, 3) if ec90 >= 0 else -1.0,
    }


def _interpolate_ec(
    doses: list[float],
    reductions: list[float],
    target_pct: float,
) -> float:
    """Interpolacao linear para encontrar dose que atinge target_pct de reducao.

    Args:
        doses: Doses em ordem crescente.
        reductions: Reducao percentual correspondente.
        target_pct: Percentual alvo (ex: 50 para EC50).

    Returns:
        Dose interpolada, ou -1.0 se nao encontrada.
    """
    for i in range(len(doses) - 1):
        r1, r2 = reductions[i], reductions[i + 1]
        d1, d2 = doses[i], doses[i + 1]

        if r1 <= target_pct <= r2:
            # Interpolacao linear
            if abs(r2 - r1) < 1e-10:
                return d1
            frac = (target_pct - r1) / (r2 - r1)
            return d1 + frac * (d2 - d1)

    # Se a primeira dose ja supera o alvo
    if len(reductions) > 0 and reductions[0] >= target_pct:
        return doses[0]

    # Se nenhuma dose atinge o alvo
    return -1.0
