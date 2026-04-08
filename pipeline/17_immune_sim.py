"""Simplified ODE-based immune response simulation for the vaccine construct.

Since C-ImmSim may not compile on Apple Silicon, this module implements
a compartmental ODE model of the immune cascade following a three-dose
vaccination protocol.  The model tracks antigen-presenting cells, T helper
cells (CD4+), cytotoxic T cells (CD8+), B cells, antibodies, and memory
cells over 365 days.

Parameters are calibrated from literature values and adjusted dynamically
based on the vaccine construct properties (epitope count, average IC50,
adjuvant type).

Usage:
    python -m pipeline.17_immune_sim
    python -m pipeline.17_immune_sim --force
    python -m pipeline.17_immune_sim --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Final

import numpy as np
from scipy.integrate import solve_ivp

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("immune_sim")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONSTRUCT_CARD_PATH: Final[str] = "results/construct/construct_card.json"
OUTPUT_DIR: Final[str] = "results/immune_sim"

# Vaccination protocol: (day, relative dose amount)
DOSE_SCHEDULE: Final[list[tuple[float, float]]] = [
    (0.0, 1.0),    # Day 0  -- priming
    (28.0, 1.0),   # Day 28 -- first boost
    (56.0, 1.0),   # Day 56 -- second boost
]

SIMULATION_DAYS: Final[int] = 365

# ODE parameters (calibrated from literature)
DEFAULT_PARAMS: Final[dict[str, float]] = {
    "k_uptake": 0.1,           # DC antigen uptake rate
    "k_th_activation": 0.05,   # Th cell activation by DCs
    "k_tc_activation": 0.08,   # Tc activation (CTL epitopes)
    "k_b_activation": 0.03,    # B cell activation by Th
    "k_ab_production": 0.1,    # Antibody production by B cells
    "k_memory": 0.01,          # Memory cell formation rate
    "k_decay_ag": 0.5,         # Antigen decay
    "k_decay_tc": 0.02,        # Tc cell decay
    "k_decay_ab": 0.01,        # Antibody decay
    "th1_bias": 0.7,           # Th1/Th2 ratio (L7/L12 adjuvant pushes Th1)
    "ifn_gamma_factor": 1.5,   # IFN-gamma boost from Th1 response
}

# State variable names and their indices
STATE_NAMES: Final[list[str]] = ["Ag", "A", "Th", "Tc", "B", "Ab", "M"]
IDX_AG: Final[int] = 0
IDX_A: Final[int] = 1
IDX_TH: Final[int] = 2
IDX_TC: Final[int] = 3
IDX_B: Final[int] = 4
IDX_AB: Final[int] = 5
IDX_M: Final[int] = 6

# Initial conditions: all populations at baseline
INITIAL_CONDITIONS: Final[list[float]] = [
    0.0,    # Ag  -- no antigen before vaccination
    0.01,   # A   -- baseline DCs
    0.01,   # Th  -- baseline CD4+
    0.01,   # Tc  -- baseline CD8+
    0.01,   # B   -- baseline B cells
    0.01,   # Ab  -- baseline antibodies
    0.01,   # M   -- baseline memory
]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> dict | list:
    """Load and return a JSON file."""
    with open(path) as fh:
        return json.load(fh)


def _write_json(data: object, path: Path) -> None:
    """Write *data* as pretty-printed JSON to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
    logger.info("Wrote %s", path)


def _write_text(text: str, path: Path) -> None:
    """Write plain text to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    logger.info("Wrote %s", path)


# ---------------------------------------------------------------------------
# Vaccine parameter loading
# ---------------------------------------------------------------------------


def _load_vaccine_params(card_path: Path) -> dict:
    """Extract simulation-relevant parameters from the construct card.

    Returns a dict with:
        epitope_count   -- total number of epitopes
        avg_ic50        -- average IC50 across all epitopes (nM)
        ctl_count       -- number of CTL (9-mer) epitopes
        htl_count       -- number of HTL (15-mer) epitopes
        adjuvant        -- adjuvant name string
    """
    card = _load_json(card_path)
    epitopes: list[dict] = card.get("epitopes", [])

    ic50_values = [ep.get("ic50", 500.0) for ep in epitopes]
    avg_ic50 = sum(ic50_values) / len(ic50_values) if ic50_values else 500.0

    # Classify by peptide length: 9-mers are typically CTL, 15-mers HTL.
    # All epitopes in this construct are 9-mers (MHC-I / CTL).
    ctl_count = sum(1 for ep in epitopes if len(ep.get("peptide", "")) <= 11)
    htl_count = sum(1 for ep in epitopes if len(ep.get("peptide", "")) > 11)

    return {
        "epitope_count": card.get("epitope_count", len(epitopes)),
        "avg_ic50": round(avg_ic50, 2),
        "ctl_count": ctl_count,
        "htl_count": htl_count,
        "adjuvant": card.get("adjuvant", "L7L12"),
    }


def _adjust_params(
    base_params: dict[str, float],
    vaccine: dict,
) -> dict[str, float]:
    """Adjust ODE parameters based on vaccine construct properties.

    Rules:
        - More CTL epitopes increase k_tc_activation.
        - Lower average IC50 (stronger binding) increases k_th_activation.
        - L7/L12 adjuvant sets th1_bias = 0.7 (strong Th1).
        - CpG co-adjuvant boosts ifn_gamma_factor to 2.0.
    """
    params = dict(base_params)

    # CTL epitope boost: scale from baseline (10 epitopes as reference)
    ctl_ratio = vaccine["ctl_count"] / 10.0
    params["k_tc_activation"] = base_params["k_tc_activation"] * max(0.5, min(ctl_ratio, 2.0))

    # IC50 binding strength: lower IC50 = stronger binding = higher activation
    # Reference IC50 = 100 nM; half that doubles the activation rate
    ic50_factor = 100.0 / max(vaccine["avg_ic50"], 1.0)
    ic50_factor = max(0.5, min(ic50_factor, 3.0))
    params["k_th_activation"] = base_params["k_th_activation"] * ic50_factor

    # Adjuvant effects
    adjuvant = vaccine.get("adjuvant", "").upper()
    if "L7" in adjuvant or "L12" in adjuvant:
        params["th1_bias"] = 0.7
    else:
        params["th1_bias"] = 0.5

    # CpG co-adjuvant is assumed for all Leishmania constructs
    params["ifn_gamma_factor"] = 2.0

    return params


# ---------------------------------------------------------------------------
# ODE system
# ---------------------------------------------------------------------------


def immune_odes(
    t: float,
    y: list[float],
    params: dict[str, float],
    doses: list[tuple[float, float]],
) -> list[float]:
    """System of ODEs modelling the immune response after vaccination.

    State variables (normalised to [0, ~max]):
        Ag -- Antigen concentration
        A  -- Activated antigen-presenting cells (DCs)
        Th -- T helper cells (CD4+)
        Tc -- Cytotoxic T cells (CD8+)
        B  -- B cells
        Ab -- Antibodies (titre proxy)
        M  -- Memory cells

    Args:
        t: Current time (days).
        y: State vector [Ag, A, Th, Tc, B, Ab, M].
        params: Kinetic parameters dict.
        doses: List of (day, amount) tuples for each injection.

    Returns:
        Derivatives [dAg, dA, dTh, dTc, dB, dAb, dM].
    """
    Ag, A, Th, Tc, B, Ab, M = y

    # Dose injection: modelled as a smoothed pulse (half-day window)
    dose_input = sum(
        dose_amount
        for dose_time, dose_amount in doses
        if abs(t - dose_time) < 0.5
    )

    # Antigen dynamics: injected, taken up by DCs, and decayed
    dAg = dose_input - params["k_decay_ag"] * Ag - params["k_uptake"] * A * Ag

    # DC activation: stimulated by antigen, self-limited, decays
    dA = params["k_uptake"] * Ag * (1.0 - A) - 0.1 * A

    # T helper (CD4+): activated by DCs, boosted by memory recall
    dTh = (
        params["k_th_activation"] * A * (1.0 - Th)
        - 0.01 * Th
        + 0.05 * M
    )

    # Cytotoxic T (CD8+): activated by DCs + Th help + IFN-gamma boost
    dTc = (
        params["k_tc_activation"] * A * Th * params["ifn_gamma_factor"]
        - params["k_decay_tc"] * Tc
        + 0.05 * M
    )

    # B cells: activated by Th help, self-limited
    dB = params["k_b_activation"] * Th * (1.0 - B) - 0.01 * B

    # Antibodies: produced by B cells, decay slowly
    dAb = params["k_ab_production"] * B - params["k_decay_ab"] * Ab

    # Memory cells: formed from all adaptive lineages, very slow decay
    dM = params["k_memory"] * (Th + Tc + B) * (1.0 - M) - 0.001 * M

    return [dAg, dA, dTh, dTc, dB, dAb, dM]


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------


def _run_simulation(
    params: dict[str, float],
    doses: list[tuple[float, float]],
    t_end: int = SIMULATION_DAYS,
) -> dict:
    """Integrate the ODE system and return structured results.

    Returns a dict with:
        t        -- time array (days)
        states   -- dict mapping state name to its time series
        params   -- the parameters used
        doses    -- the dose schedule used
    """
    t_span = (0.0, float(t_end))
    t_eval = np.linspace(0.0, float(t_end), t_end * 10 + 1)

    sol = solve_ivp(
        fun=lambda t, y: immune_odes(t, y, params, doses),
        t_span=t_span,
        y0=INITIAL_CONDITIONS,
        method="RK45",
        t_eval=t_eval,
        max_step=0.5,
        rtol=1e-8,
        atol=1e-10,
    )

    if not sol.success:
        logger.error("ODE integration failed: %s", sol.message)
        sys.exit(1)

    states: dict[str, np.ndarray] = {}
    for idx, name in enumerate(STATE_NAMES):
        states[name] = sol.y[idx]

    return {
        "t": sol.t,
        "states": states,
        "params": params,
        "doses": doses,
    }


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------


def _peak_info(t: np.ndarray, values: np.ndarray) -> dict:
    """Find the peak value and its corresponding day."""
    idx = int(np.argmax(values))
    return {
        "peak_level": round(float(values[idx]), 4),
        "day_of_peak": round(float(t[idx]), 1),
        "final_level": round(float(values[-1]), 4),
    }


def _compute_th1_th2_timeseries(
    sim: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute Th1 and Th2 relative levels over time.

    The Th1/Th2 balance is modelled as a fixed split governed by the
    th1_bias parameter, modulated by the IFN-gamma factor for Th1.
    """
    th = sim["states"]["Th"]
    bias = sim["params"]["th1_bias"]
    ifn = sim["params"]["ifn_gamma_factor"]

    th1 = th * bias * ifn / (1.0 + ifn * bias)
    th2 = th * (1.0 - bias) / (1.0 + ifn * bias)
    return th1, th2


def _analyse_results(sim: dict) -> dict:
    """Extract summary statistics from the simulation."""
    t = sim["t"]
    states = sim["states"]

    # Peak info for each adaptive cell type
    peaks: dict[str, dict] = {}
    for name in ["Th", "Tc", "B", "Ab", "M"]:
        peaks[name] = _peak_info(t, states[name])

    # Th1/Th2 balance
    th1, th2 = _compute_th1_th2_timeseries(sim)
    th1_pct = float(th1[-1] / (th1[-1] + th2[-1] + 1e-12)) * 100.0

    th1_peak_idx = int(np.argmax(th1))
    th1_peak_day = round(float(t[th1_peak_idx]), 1)

    # Memory after each dose
    memory = states["M"]
    dose_days = [d[0] for d in sim["doses"]]
    memory_after_doses: list[float] = []
    for dose_day in dose_days:
        # Measure memory 14 days after each dose
        measure_day = dose_day + 14.0
        idx = int(np.argmin(np.abs(t - measure_day)))
        memory_after_doses.append(round(float(memory[idx]), 4))

    # Memory stability (final / peak)
    m_peak = peaks["M"]["peak_level"]
    m_final = peaks["M"]["final_level"]
    memory_stability = round(m_final / max(m_peak, 1e-12) * 100.0, 1)

    # Estimated protection duration based on memory half-life
    # Memory decay rate is 0.001 per day; half-life = ln(2)/0.001
    memory_half_life_days = round(np.log(2) / 0.001, 0)
    # Protection threshold: memory > 50% of peak
    if m_final > 0.5 * m_peak:
        estimated_protection = f">{int(memory_half_life_days)} days (~{int(memory_half_life_days / 365)} years)"
    elif m_final > 0.25 * m_peak:
        estimated_protection = f"~{int(memory_half_life_days * 0.5)} days"
    else:
        estimated_protection = "Short-lived; booster recommended within 6 months"

    return {
        "peaks": peaks,
        "th1_dominance_pct": round(th1_pct, 1),
        "th1_peak_day": th1_peak_day,
        "memory_after_doses": memory_after_doses,
        "memory_stability_pct": memory_stability,
        "estimated_protection": estimated_protection,
        "memory_half_life_days": int(memory_half_life_days),
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _generate_plots(sim: dict, analysis: dict, output_dir: Path) -> list[str]:
    """Generate immune response plots and save to output_dir.

    Returns list of generated file paths.  If matplotlib is unavailable
    the function logs a warning and returns an empty list.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning(
            "matplotlib not available; skipping plot generation."
        )
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[str] = []
    t = sim["t"]
    dose_days = [d[0] for d in sim["doses"]]

    # -- Plot 1: Immune kinetics --
    fig, ax = plt.subplots(figsize=(12, 7))

    plot_map = {
        "CD4+ T helper (Th)": sim["states"]["Th"],
        "CD8+ Cytotoxic (Tc)": sim["states"]["Tc"],
        "B cells": sim["states"]["B"],
        "Antibodies": sim["states"]["Ab"],
        "Memory cells": sim["states"]["M"],
    }
    colors = ["#2196F3", "#F44336", "#4CAF50", "#FF9800", "#9C27B0"]

    for (label, values), color in zip(plot_map.items(), colors):
        ax.plot(t, values, label=label, color=color, linewidth=2)

    for day in dose_days:
        ax.axvline(x=day, color="gray", linestyle="--", alpha=0.5, linewidth=1)
    ax.axvline(x=dose_days[0], color="gray", linestyle="--", alpha=0.5,
               linewidth=1, label="Dose days")

    ax.set_xlabel("Days post-vaccination", fontsize=12)
    ax.set_ylabel("Relative cell population / titre", fontsize=12)
    ax.set_title("Immune Response Kinetics After Vaccination", fontsize=14)
    ax.legend(loc="upper right", fontsize=10)
    ax.set_xlim(0, SIMULATION_DAYS)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    path = output_dir / "immune_kinetics.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    generated.append(str(path))
    logger.info("Saved %s", path)

    # -- Plot 2: Th1/Th2 balance --
    th1, th2 = _compute_th1_th2_timeseries(sim)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    ax1.plot(t, th1, label="Th1 (IFN-gamma, TNF-alpha)", color="#D32F2F",
             linewidth=2)
    ax1.plot(t, th2, label="Th2 (IL-4, IL-10)", color="#1976D2", linewidth=2)
    for day in dose_days:
        ax1.axvline(x=day, color="gray", linestyle="--", alpha=0.5)
    ax1.set_ylabel("Relative level", fontsize=12)
    ax1.set_title("Th1/Th2 Balance Over Time", fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Ratio subplot
    ratio = th1 / (th2 + 1e-12)
    ax2.plot(t, ratio, color="#388E3C", linewidth=2)
    ax2.axhline(y=1.0, color="gray", linestyle=":", alpha=0.7,
                label="Th1 = Th2")
    for day in dose_days:
        ax2.axvline(x=day, color="gray", linestyle="--", alpha=0.5)
    ax2.set_xlabel("Days post-vaccination", fontsize=12)
    ax2.set_ylabel("Th1 / Th2 ratio", fontsize=12)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()

    path = output_dir / "th1_th2_balance.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    generated.append(str(path))
    logger.info("Saved %s", path)

    # -- Plot 3: Memory formation --
    memory = sim["states"]["M"]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(t, memory, color="#9C27B0", linewidth=2.5, label="Memory cells")
    ax.fill_between(t, 0, memory, color="#CE93D8", alpha=0.3)

    for i, day in enumerate(dose_days):
        ax.axvline(x=day, color="gray", linestyle="--", alpha=0.5)
        label = ["Dose 1 (prime)", "Dose 2 (boost)", "Dose 3 (boost)"][i]
        ax.annotate(
            label,
            xy=(day, float(memory[int(np.argmin(np.abs(t - day)))])),
            xytext=(day + 15, float(np.max(memory)) * (0.3 + 0.2 * i)),
            fontsize=9,
            arrowprops={"arrowstyle": "->", "color": "gray"},
            color="gray",
        )

    # Mark memory after each dose (+14 days)
    for i, level in enumerate(analysis["memory_after_doses"]):
        measure_day = dose_days[i] + 14.0
        ax.plot(measure_day, level, "o", color="#6A1B9A", markersize=8)
        ax.annotate(
            f"{level:.3f}",
            xy=(measure_day, level),
            xytext=(measure_day + 5, level + 0.02),
            fontsize=9,
            color="#6A1B9A",
        )

    ax.set_xlabel("Days post-vaccination", fontsize=12)
    ax.set_ylabel("Memory cell level (normalised)", fontsize=12)
    ax.set_title("Memory Cell Formation After Each Dose", fontsize=14)
    ax.legend(loc="lower right", fontsize=10)
    ax.set_xlim(0, SIMULATION_DAYS)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    path = output_dir / "memory_formation.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    generated.append(str(path))
    logger.info("Saved %s", path)

    return generated


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _build_report(
    vaccine: dict,
    params: dict[str, float],
    analysis: dict,
    plot_paths: list[str],
) -> str:
    """Generate the immune simulation Markdown report."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    peaks = analysis["peaks"]

    # Determine checkmarks for ideal profile
    th1_pct = analysis["th1_dominance_pct"]
    strong_th1 = th1_pct > 60.0
    moderate_ctl = peaks["Tc"]["peak_level"] > 0.05
    low_th2 = th1_pct > 60.0  # Th2 is low if Th1 dominates
    lasting_memory = analysis["memory_stability_pct"] > 50.0

    check = lambda ok: "Yes" if ok else "No"

    # Adjuvant effect description
    adj = vaccine.get("adjuvant", "").upper()
    if "L7" in adj or "L12" in adj:
        adj_effect = "strong Th1 bias"
    else:
        adj_effect = "moderate Th1 bias"

    # Memory after dose 2 as percentage of final
    m_after_d2 = analysis["memory_after_doses"][1] if len(analysis["memory_after_doses"]) > 1 else 0.0
    m_final = peaks["M"]["final_level"]
    m_d2_pct = round(m_after_d2 / max(m_final, 1e-12) * 100.0, 1)

    lines: list[str] = []
    lines.append("# Immune Response Simulation")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")

    # -- Vaccination Protocol --
    lines.append("## Vaccination Protocol")
    lines.append("- Dose 1: Day 0 (priming)")
    lines.append("- Dose 2: Day 28 (boosting)")
    lines.append("- Dose 3: Day 56 (boosting)")
    lines.append(f"- Simulation: {SIMULATION_DAYS} days")
    lines.append(f"- Epitopes in construct: {vaccine['epitope_count']} "
                 f"({vaccine['ctl_count']} CTL, {vaccine['htl_count']} HTL)")
    lines.append(f"- Average IC50: {vaccine['avg_ic50']} nM")
    lines.append(f"- Adjuvant: {vaccine.get('adjuvant', 'L7L12')}")
    lines.append("")

    # -- Model Parameters --
    lines.append("## Adjusted Model Parameters")
    lines.append("")
    lines.append("| Parameter | Value | Description |")
    lines.append("|-----------|:-----:|-------------|")
    param_descriptions = {
        "k_uptake": "DC antigen uptake rate",
        "k_th_activation": "Th cell activation by DCs",
        "k_tc_activation": "Tc activation (CTL epitopes)",
        "k_b_activation": "B cell activation by Th",
        "k_ab_production": "Antibody production by B cells",
        "k_memory": "Memory cell formation rate",
        "k_decay_ag": "Antigen decay",
        "k_decay_tc": "Tc cell decay",
        "k_decay_ab": "Antibody decay",
        "th1_bias": "Th1/Th2 ratio bias",
        "ifn_gamma_factor": "IFN-gamma boost factor",
    }
    for key, desc in param_descriptions.items():
        val = params.get(key, 0.0)
        lines.append(f"| {key} | {val:.4f} | {desc} |")
    lines.append("")

    # -- Peak Immune Responses --
    lines.append("## Results")
    lines.append("")
    lines.append("### Peak Immune Responses")
    lines.append("| Cell Type | Peak Level | Day of Peak | After 365 days |")
    lines.append("|-----------|:----------:|:-----------:|:--------------:|")

    display_names = {
        "Th": "CD4+ T (Th1)",
        "Tc": "CD8+ T (CTL)",
        "B": "B cells",
        "Ab": "Antibodies",
        "M": "Memory cells",
    }
    for key in ["Th", "Tc", "B", "Ab", "M"]:
        p = peaks[key]
        lines.append(
            f"| {display_names[key]} | {p['peak_level']:.4f} "
            f"| Day {p['day_of_peak']:.0f} | {p['final_level']:.4f} |"
        )
    lines.append("")

    # -- Th1/Th2 Balance --
    lines.append("### Th1/Th2 Balance")
    lines.append(f"- Th1 dominance: {th1_pct:.1f}% (target: >60% for Leishmania)")
    lines.append(f"- Peak IFN-gamma response: Day {analysis['th1_peak_day']:.0f}")
    lines.append(f"- L7/L12 adjuvant effect: {adj_effect}")
    lines.append("")

    # -- Memory Response --
    lines.append("### Memory Response")
    lines.append(f"- Memory cells established after dose 2: {m_d2_pct:.1f}%")
    lines.append(f"- Stable memory at Day 365: {analysis['memory_stability_pct']:.1f}%")
    lines.append(f"- Predicted duration of protection: {analysis['estimated_protection']}")
    lines.append("")

    # -- Comparison with Ideal Profile --
    lines.append("### Comparison with Ideal Profile")
    lines.append("For anti-Leishmania vaccination, the ideal immune profile is:")
    lines.append(f"- Strong Th1 (IFN-gamma, TNF-alpha): {check(strong_th1)}")
    lines.append(f"- Moderate CD8+ cytotoxic response: {check(moderate_ctl)}")
    lines.append(f"- Low Th2 (no IL-4/IL-10 dominance): {check(low_th2)}")
    lines.append(f"- Long-lasting memory: {check(lasting_memory)}")
    lines.append("")

    # -- Plot references --
    if plot_paths:
        lines.append("## Generated Plots")
        for p in plot_paths:
            fname = Path(p).name
            lines.append(f"- ![{fname}]({fname})")
        lines.append("")

    # -- Methods note --
    lines.append("## Methods")
    lines.append(
        "This simulation uses a simplified compartmental ODE model with "
        "7 state variables (antigen, DCs, CD4+ T, CD8+ T, B cells, "
        "antibodies, memory cells). The system is integrated using the "
        "RK45 method via scipy.integrate.solve_ivp. Parameters are "
        "calibrated from published immunological literature and adjusted "
        "based on the vaccine construct properties (epitope count, average "
        "IC50, adjuvant type). This model provides a qualitative overview "
        "of the expected immune dynamics; quantitative predictions require "
        "validation with experimental data."
    )
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def simulate_immune_response(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the immune response simulation pipeline.

    Steps:
        1. Load vaccine parameters from the construct card.
        2. Adjust ODE parameters based on vaccine properties.
        3. Integrate the ODE system over 365 days with 3 doses.
        4. Analyse peak responses, Th1/Th2 balance, and memory.
        5. Generate plots (unless matplotlib is unavailable).
        6. Write the simulation report.

    Args:
        force: If True, overwrite existing output even if present.
        dry_run: If True, run the simulation (pure math, no API calls)
            but skip plot generation if matplotlib is unavailable.

    Returns:
        Path to the generated immune simulation report.
    """
    output_base = Path(OUTPUT_DIR)
    report_path = output_base / "immune_report.md"

    if report_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            report_path,
        )
        return str(report_path)

    # ------------------------------------------------------------------
    # 1. Load vaccine parameters
    # ------------------------------------------------------------------
    card_path = Path(CONSTRUCT_CARD_PATH)
    if not card_path.exists():
        logger.error("Construct card not found at %s.", card_path)
        sys.exit(1)

    vaccine = _load_vaccine_params(card_path)
    logger.info(
        "Loaded vaccine: %d epitopes (%d CTL, %d HTL), avg IC50=%.1f nM, "
        "adjuvant=%s",
        vaccine["epitope_count"],
        vaccine["ctl_count"],
        vaccine["htl_count"],
        vaccine["avg_ic50"],
        vaccine["adjuvant"],
    )

    # ------------------------------------------------------------------
    # 2. Adjust parameters
    # ------------------------------------------------------------------
    params = _adjust_params(DEFAULT_PARAMS, vaccine)
    logger.info(
        "Adjusted params: k_th_activation=%.4f, k_tc_activation=%.4f, "
        "th1_bias=%.2f, ifn_gamma=%.1f",
        params["k_th_activation"],
        params["k_tc_activation"],
        params["th1_bias"],
        params["ifn_gamma_factor"],
    )

    # ------------------------------------------------------------------
    # 3. Run simulation
    # ------------------------------------------------------------------
    logger.info(
        "Running ODE simulation: %d days, %d doses ...",
        SIMULATION_DAYS,
        len(DOSE_SCHEDULE),
    )
    sim = _run_simulation(params, list(DOSE_SCHEDULE))
    logger.info("Simulation complete: %d time points.", len(sim["t"]))

    # ------------------------------------------------------------------
    # 4. Analyse results
    # ------------------------------------------------------------------
    analysis = _analyse_results(sim)
    logger.info(
        "Peak Th=%.4f (day %.0f), Tc=%.4f (day %.0f), Ab=%.4f (day %.0f)",
        analysis["peaks"]["Th"]["peak_level"],
        analysis["peaks"]["Th"]["day_of_peak"],
        analysis["peaks"]["Tc"]["peak_level"],
        analysis["peaks"]["Tc"]["day_of_peak"],
        analysis["peaks"]["Ab"]["peak_level"],
        analysis["peaks"]["Ab"]["day_of_peak"],
    )
    logger.info(
        "Th1 dominance: %.1f%%, memory stability: %.1f%%",
        analysis["th1_dominance_pct"],
        analysis["memory_stability_pct"],
    )

    # ------------------------------------------------------------------
    # 5. Generate plots
    # ------------------------------------------------------------------
    plot_paths = _generate_plots(sim, analysis, output_base)

    # ------------------------------------------------------------------
    # 6. Save simulation data and report
    # ------------------------------------------------------------------
    sim_data = {
        "generated": datetime.now().isoformat(),
        "vaccine_params": vaccine,
        "ode_params": params,
        "dose_schedule": [
            {"day": d, "amount": a} for d, a in DOSE_SCHEDULE
        ],
        "simulation_days": SIMULATION_DAYS,
        "peaks": analysis["peaks"],
        "th1_dominance_pct": analysis["th1_dominance_pct"],
        "th1_peak_day": analysis["th1_peak_day"],
        "memory_after_doses": analysis["memory_after_doses"],
        "memory_stability_pct": analysis["memory_stability_pct"],
        "estimated_protection": analysis["estimated_protection"],
        "plots": [str(p) for p in plot_paths],
    }
    _write_json(sim_data, output_base / "immune_sim_data.json")

    report = _build_report(vaccine, params, analysis, plot_paths)
    _write_text(report, report_path)

    logger.info("Immune simulation complete. Report: %s", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley -- Simplified ODE-based immune response simulation "
            "for the multi-epitope vaccine construct."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing outputs even if present.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Run the simulation (pure math) but skip plot generation "
            "if matplotlib is unavailable."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = simulate_immune_response(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
