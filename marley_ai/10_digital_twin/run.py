"""Execucao do modulo 10_digital_twin --- Digital twin imunologico canino.

Simula resposta imune canina ao construto vacinal + terapia ASO usando
modelo de 7 compartimentos (ODEs) integrado com scipy.integrate.solve_ivp.

Compartimentos:
    N  - Naive T cells
    E  - Effector T cells (Th1)
    M  - Memory T cells
    B  - B cells (antibody-producing)
    A  - Antibody titer
    P  - Parasite load (L. infantum amastigotes)
    ASO_c - ASO concentration (MRL-ASO-001)

Cenarios:
    1. Vaccine only (prime + boost, sem ASO)
    2. Vaccine + ASO dual therapy

Uso:
    python -m marley_ai.10_digital_twin.run
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register


# ---------------------------------------------------------------------------
# ODE parameters --- biologically plausible for canine VL
# ---------------------------------------------------------------------------

PARAMS: dict[str, float] = {
    # --- Naive T cells ---
    "lambda_N": 50.0,       # naive T cell production rate (cells/day)
    "d_N": 0.01,            # naive T cell death rate (1/day)
    # --- T cell activation ---
    "beta_vax": 0.0002,     # T cell activation rate by vaccine antigen
    "beta_par": 5e-8,       # T cell activation rate by parasite (cross-presentation)
    "d_E": 0.05,            # effector T cell death rate
    "alpha_M": 0.02,        # memory cell formation rate
    # --- Memory ---
    "d_M": 0.001,           # memory cell death rate (long-lived)
    "sigma": 0.05,          # memory reactivation rate (by antigen)
    "sigma_par": 1e-7,      # memory reactivation by parasite
    # --- B cells / antibodies ---
    "beta_B": 0.0001,       # B cell activation rate
    "d_B": 0.05,            # B cell death rate
    "gamma_A": 5.0,         # antibody production rate per B cell
    "d_A": 0.02,            # antibody decay rate (IgG half-life ~21d in dogs)
    # --- Parasite ---
    "r_P": 0.15,            # parasite growth rate (L. infantum intracellular)
    "K_P": 1e7,             # carrying capacity (visceral organs)
    # --- Killing rates ---
    "kill_E": 1e-5,         # effector killing rate (Th1/macrophage mediated)
    "kill_A": 5e-6,         # antibody-mediated killing (opsonization + ADCC)
    "kill_ASO": 0.002,      # ASO-mediated killing (SL RNA disruption)
    # --- ASO pharmacokinetics ---
    "d_ASO": 0.04,          # ASO clearance rate (tissue half-life ~17h)
    # --- Vaccine antigen kinetics ---
    "Ag_decay": 0.15,       # antigen decay from vaccine depot (1/day)
}

# ---------------------------------------------------------------------------
# Dosing schedules
# ---------------------------------------------------------------------------

# Vaccine: prime day 0, boost day 21 (antigen pulse = 1000 units each)
VACCINE_SCHEDULE: list[tuple[float, float]] = [
    (0.0, 1000.0),    # prime
    (21.0, 1000.0),   # boost
]

# ASO therapy: starts day 28, SC injection dose=500, every 7 days for 4 weeks
ASO_SCHEDULE: list[tuple[float, float]] = [
    (28.0, 500.0),
    (35.0, 500.0),
    (42.0, 500.0),
    (49.0, 500.0),
]

# Initial conditions: N=1000, E=0, M=0, B=0, A=0, P=1e5, ASO_c=0
# P=1e5: moderate visceral Leishmania burden at treatment start
Y0: list[float] = [1000.0, 0.0, 0.0, 0.0, 0.0, 1e5, 0.0]

# Simulation span
T_SPAN: tuple[float, float] = (0.0, 365.0)
N_EVAL_POINTS: int = 366


# ---------------------------------------------------------------------------
# ODE system
# ---------------------------------------------------------------------------

def _antigen_signal(
    t: float,
    vaccine_schedule: list[tuple[float, float]],
    ag_decay: float = 0.3,
) -> float:
    """Compute antigen signal at time t from decaying vaccine pulses.

    Args:
        t: Current time in days.
        vaccine_schedule: List of (day, dose) tuples.
        ag_decay: Antigen decay rate (1/day). Controls depot persistence.

    Returns:
        Total antigen level at time t.
    """
    ag = 0.0
    for t_dose, dose in vaccine_schedule:
        if t >= t_dose:
            ag += dose * np.exp(-ag_decay * (t - t_dose))
    return ag


def _aso_input(
    t: float,
    aso_schedule: list[tuple[float, float]],
) -> float:
    """Compute ASO input rate at time t (bolus over 0.1 days)."""
    for t_dose, dose in aso_schedule:
        if t_dose <= t < t_dose + 0.1:
            return dose / 0.1
    return 0.0


def deriv(
    t: float,
    y: list[float],
    params: dict[str, float],
    vaccine_schedule: list[tuple[float, float]],
    aso_schedule: list[tuple[float, float]],
) -> list[float]:
    """Right-hand side of the 7-compartment canine immune response ODE.

    The model includes dual immune activation: vaccine antigen drives the
    initial T cell priming, while parasite load provides ongoing immune
    stimulation via cross-presentation. This coupling ensures the immune
    response tracks parasite burden and memory cells can maintain control
    after ASO-driven reduction.

    Args:
        t: Current time (days).
        y: State vector [N, E, M, B, A, P, ASO_c].
        params: Dictionary of kinetic parameters.
        vaccine_schedule: List of (day, dose) for vaccine antigen pulses.
        aso_schedule: List of (day, dose) for ASO injections.

    Returns:
        Derivatives [dN, dE, dM, dB, dA, dP, dASO_c].
    """
    N, E, M, B, A, P, ASO_c = y

    # Enforce non-negativity for biological quantities
    N = max(N, 0.0)
    E = max(E, 0.0)
    M = max(M, 0.0)
    B = max(B, 0.0)
    A = max(A, 0.0)
    P = max(P, 0.0)
    ASO_c = max(ASO_c, 0.0)

    # Antigen signal from vaccine (decaying exponential pulse)
    Ag = _antigen_signal(t, vaccine_schedule, ag_decay=params["Ag_decay"])

    # ASO dosing input
    aso_in = _aso_input(t, aso_schedule)

    # Combined activation signal: vaccine antigen + parasite cross-presentation
    # Vaccine antigen drives initial priming; parasite provides ongoing stimulus
    activation_vax = params["beta_vax"] * Ag
    activation_par = params["beta_par"] * P
    total_activation = activation_vax + activation_par

    # Memory reactivation: by vaccine antigen or parasite encounter
    reactivation = params["sigma"] * Ag + params["sigma_par"] * P

    # ODEs
    dN = params["lambda_N"] - params["d_N"] * N - total_activation * N
    dE = (
        total_activation * N
        + reactivation * M
        - params["d_E"] * E
        - params["alpha_M"] * E
    )
    dM = params["alpha_M"] * E - params["d_M"] * M - reactivation * M
    dB = params["beta_B"] * (Ag + params["beta_par"] * P) * E - params["d_B"] * B
    dA = params["gamma_A"] * B - params["d_A"] * A
    dP = (
        params["r_P"] * P * (1.0 - P / params["K_P"])
        - params["kill_E"] * E * P
        - params["kill_A"] * A * P
        - params["kill_ASO"] * ASO_c * P
    )
    dASO = aso_in - params["d_ASO"] * ASO_c

    return [dN, dE, dM, dB, dA, dP, dASO]


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------

def run_simulation(
    params: dict[str, float],
    vaccine_schedule: list[tuple[float, float]],
    aso_schedule: list[tuple[float, float]],
    y0: list[float],
    t_span: tuple[float, float],
    t_eval: np.ndarray,
) -> Any:
    """Integrate the ODE system using scipy solve_ivp (RK45).

    Args:
        params: Kinetic parameters.
        vaccine_schedule: Vaccine dosing schedule.
        aso_schedule: ASO dosing schedule.
        y0: Initial conditions.
        t_span: (t_start, t_end).
        t_eval: Array of time points for output.

    Returns:
        scipy OdeResult object.
    """
    sol = solve_ivp(
        fun=deriv,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        args=(params, vaccine_schedule, aso_schedule),
        method="RK45",
        max_step=0.5,
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")
    return sol


# ---------------------------------------------------------------------------
# Metric extraction
# ---------------------------------------------------------------------------

def extract_metrics(
    sol: Any,
    label: str,
) -> dict[str, Any]:
    """Extract key immunological metrics from an ODE solution.

    Args:
        sol: scipy OdeResult with .t and .y arrays.
        label: Human-readable label for this scenario.

    Returns:
        Dictionary of key metrics.
    """
    t = sol.t
    E_arr = sol.y[1]
    M_arr = sol.y[2]
    A_arr = sol.y[4]
    P_arr = sol.y[5]

    # Peak effector T cells
    peak_E_idx = int(np.argmax(E_arr))
    peak_E_value = float(E_arr[peak_E_idx])
    peak_E_day = float(t[peak_E_idx])

    # Parasite clearance time (day when P < 1)
    clearance_indices = np.where(P_arr < 1.0)[0]
    if len(clearance_indices) > 0:
        clearance_day = float(t[clearance_indices[0]])
    else:
        clearance_day = float("inf")

    # Time to 50% parasite reduction
    p0 = P_arr[0]
    half_indices = np.where(P_arr < p0 * 0.5)[0]
    if len(half_indices) > 0:
        half_reduction_day = float(t[half_indices[0]])
    else:
        half_reduction_day = float("inf")

    # Values at day 365 (last time point)
    memory_365 = float(M_arr[-1])
    antibody_365 = float(A_arr[-1])
    parasite_365 = float(P_arr[-1])

    return {
        "label": label,
        "peak_effector_cells": round(peak_E_value, 2),
        "peak_effector_day": round(peak_E_day, 1),
        "clearance_day": round(clearance_day, 1) if clearance_day != float("inf") else None,
        "half_reduction_day": round(half_reduction_day, 1) if half_reduction_day != float("inf") else None,
        "memory_cells_day365": round(memory_365, 2),
        "antibody_titer_day365": round(antibody_365, 2),
        "parasite_load_day365": round(parasite_365, 4),
    }


def _sol_to_dict(sol: Any) -> dict[str, list[float]]:
    """Convert ODE solution arrays to JSON-serializable dicts.

    Args:
        sol: scipy OdeResult.

    Returns:
        Dictionary mapping compartment names to lists of floats.
    """
    return {
        "naive": [round(float(v), 4) for v in sol.y[0]],
        "effector": [round(float(v), 4) for v in sol.y[1]],
        "memory": [round(float(v), 4) for v in sol.y[2]],
        "b_cells": [round(float(v), 4) for v in sol.y[3]],
        "antibody": [round(float(v), 4) for v in sol.y[4]],
        "parasite": [round(float(v), 4) for v in sol.y[5]],
        "aso_concentration": [round(float(v), 4) for v in sol.y[6]],
    }


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

@register("10_digital_twin")
class DigitalTwinModule:
    """Modulo de digital twin imunologico canino com ODE de 7 compartimentos."""

    def __init__(self) -> None:
        self._config: AIModuleConfig | None = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao do orquestrador."""
        self._config = config

    def validate_inputs(self) -> dict[str, Any]:
        """Sem dependencias externas --- usa parametros internos."""
        return {"valid": True, "missing": []}

    def run(self) -> dict[str, Any]:
        """Executa pipeline do digital twin."""
        return main(self._config)

    def get_dependencies(self) -> list[str]:
        """Sem dependencias de outros modulos."""
        return []


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Simula resposta imune canina ao construto vacinal + ASO.

    Executa dois cenarios:
        1. Vaccine only (prime + boost, sem ASO)
        2. Vaccine + ASO dual therapy (prime + boost + 4 doses ASO)

    Compara tempos de clearance parasitario e extrai metricas
    imunologicas chave para cada cenario.

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com resultados da simulacao.
    """
    envelope = create_envelope("10_digital_twin")
    envelope["device"] = "cpu"  # numpy + scipy --- sem GPU

    with Timer() as timer:

        params = PARAMS.copy()
        t_eval = np.linspace(T_SPAN[0], T_SPAN[1], N_EVAL_POINTS)

        # -----------------------------------------------------------
        # Scenario 1: Vaccine only (no ASO)
        # -----------------------------------------------------------
        print("[10_digital_twin] Cenario 1: Vaccine only...")
        empty_aso: list[tuple[float, float]] = []
        sol_vax = run_simulation(
            params=params,
            vaccine_schedule=VACCINE_SCHEDULE,
            aso_schedule=empty_aso,
            y0=Y0,
            t_span=T_SPAN,
            t_eval=t_eval,
        )
        metrics_vax = extract_metrics(sol_vax, "vaccine_only")
        print(f"  Peak effector: {metrics_vax['peak_effector_cells']:.1f} "
              f"cells at day {metrics_vax['peak_effector_day']:.0f}")
        print(f"  Clearance day: {metrics_vax['clearance_day']}")

        # -----------------------------------------------------------
        # Scenario 2: Vaccine + ASO dual therapy
        # -----------------------------------------------------------
        print("\n[10_digital_twin] Cenario 2: Vaccine + ASO dual therapy...")
        sol_dual = run_simulation(
            params=params,
            vaccine_schedule=VACCINE_SCHEDULE,
            aso_schedule=ASO_SCHEDULE,
            y0=Y0,
            t_span=T_SPAN,
            t_eval=t_eval,
        )
        metrics_dual = extract_metrics(sol_dual, "dual_therapy")
        print(f"  Peak effector: {metrics_dual['peak_effector_cells']:.1f} "
              f"cells at day {metrics_dual['peak_effector_day']:.0f}")
        print(f"  Clearance day: {metrics_dual['clearance_day']}")

        # -----------------------------------------------------------
        # Comparison
        # -----------------------------------------------------------
        clearance_vax = metrics_vax["clearance_day"]
        clearance_dual = metrics_dual["clearance_day"]

        if clearance_vax is not None and clearance_dual is not None and clearance_vax > 0:
            improvement_pct = round(
                (clearance_vax - clearance_dual) / clearance_vax * 100.0, 1,
            )
        elif clearance_dual is not None and clearance_vax is None:
            # Vaccine alone never cleared, dual did
            improvement_pct = 100.0
        else:
            improvement_pct = 0.0

        print(f"\n[10_digital_twin] ASO improvement: {improvement_pct}%")

        # -----------------------------------------------------------
        # Build envelope
        # -----------------------------------------------------------
        envelope["status"] = "complete"
        envelope["dependencies"] = []

        # Conclusion
        clearance_vax_str = (
            f"{clearance_vax} dias"
            if clearance_vax is not None
            else "nao atingido em 365 dias"
        )
        clearance_dual_str = (
            f"{clearance_dual} dias"
            if clearance_dual is not None
            else "nao atingido em 365 dias"
        )

        envelope["summary"]["conclusion"] = (
            f"Digital twin simulou resposta imune canina por 365 dias. "
            f"Vacina + ASO: clearance em {clearance_dual_str} vs "
            f"vacina sozinha: {clearance_vax_str}. "
            f"Melhoria de {improvement_pct}%."
        )

        envelope["summary"]["key_metrics"] = {
            "simulation_days": 365,
            "n_compartments": 7,
            "n_doses": len(VACCINE_SCHEDULE) + len(ASO_SCHEDULE),
            "clearance_day_vaccine_only": clearance_vax,
            "clearance_day_dual_therapy": clearance_dual,
            "peak_effector_cells": metrics_dual["peak_effector_cells"],
            "peak_effector_day": metrics_dual["peak_effector_day"],
            "memory_cells_day365": metrics_dual["memory_cells_day365"],
            "antibody_titer_day365": metrics_dual["antibody_titer_day365"],
            "aso_improvement_pct": improvement_pct,
        }

        # Metrics for both scenarios
        envelope["metrics"] = {
            "vaccine_only": metrics_vax,
            "dual_therapy": metrics_dual,
        }

        # Full time-series data
        time_points = [round(float(v), 2) for v in sol_dual.t]

        envelope["data"] = {
            "time_points": time_points,
            "vaccine_only": _sol_to_dict(sol_vax),
            "dual_therapy": _sol_to_dict(sol_dual),
            "vaccine_schedule": [
                {"day": d, "dose": dose} for d, dose in VACCINE_SCHEDULE
            ],
            "aso_schedule": [
                {"day": d, "dose": dose} for d, dose in ASO_SCHEDULE
            ],
            "parameters": params,
            "initial_conditions": {
                "naive": Y0[0],
                "effector": Y0[1],
                "memory": Y0[2],
                "b_cells": Y0[3],
                "antibody": Y0[4],
                "parasite": Y0[5],
                "aso_concentration": Y0[6],
            },
        }

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"\n[10_digital_twin] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
