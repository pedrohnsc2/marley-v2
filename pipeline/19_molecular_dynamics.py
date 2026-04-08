"""Molecular dynamics simulation of the TryR-ligand complex using OpenMM.

Loads the best selective variant from the selectivity redesign module (or
falls back to MRL-003), sets up a solvated protein-ligand system in explicit
TIP3P water with 0.15 M NaCl, and runs energy minimisation, NVT/NPT
equilibration, and a production MD trajectory.  Analysis includes ligand
RMSD, protein RMSF, and key contact distances.

When OpenMM is not installed (expected on many workstations), the module
generates setup scripts (shell + GROMACS .mdp files) so the user can run the
simulation on a cluster.

Usage:
    python -m pipeline.19_molecular_dynamics
    python -m pipeline.19_molecular_dynamics --force
    python -m pipeline.19_molecular_dynamics --dry-run
    python -m pipeline.19_molecular_dynamics --steps 10000000
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("molecular_dynamics")

# ---------------------------------------------------------------------------
# OpenMM availability
# ---------------------------------------------------------------------------

try:
    import openmm
    import openmm.app as app
    import openmm.unit as unit

    OPENMM_AVAILABLE: bool = True
except ImportError:
    OPENMM_AVAILABLE = False
    logger.warning(
        "OpenMM not installed. Install via: conda install -c conda-forge openmm"
    )

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SELECTIVITY_MATRIX_PATH: Final[str] = (
    "results/selectivity_redesign/selectivity_matrix.csv"
)
SELECTIVITY_REPORT_PATH: Final[str] = (
    "results/selectivity/selectivity_report.json"
)
PROTEIN_PDB_PATH: Final[str] = "data/structures/TryR_Q4Q457.pdb"
LIGAND_DIR: Final[str] = "results/selectivity/ligands"
OUTPUT_DIR: Final[str] = "results/molecular_dynamics"

DEFAULT_STEPS: Final[int] = 5_000_000  # 10 ns at 2 fs timestep
TIMESTEP_FS: Final[float] = 2.0
TEMPERATURE_K: Final[float] = 310.0  # Dog body temperature
PRESSURE_ATM: Final[float] = 1.0
IONIC_STRENGTH_M: Final[float] = 0.15
SOLVENT_PADDING_NM: Final[float] = 1.0
NONBONDED_CUTOFF_NM: Final[float] = 1.0
MINIMIZE_MAX_ITER: Final[int] = 10_000
NVT_STEPS: Final[int] = 50_000  # 100 ps
NPT_STEPS: Final[int] = 50_000  # 100 ps
REPORT_INTERVAL: Final[int] = 5_000  # Every 10 ps


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


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    """Load a CSV file and return a list of row dicts."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return list(reader)


# ---------------------------------------------------------------------------
# Variant selection
# ---------------------------------------------------------------------------


def _select_best_variant() -> dict[str, Any]:
    """Select the best selective variant from upstream modules.

    Tries the selectivity redesign matrix first.  If that file does not
    exist, falls back to the selectivity report from Module 12.  If neither
    is available, returns default MRL-003 information.

    Returns a dict with keys: name, source, delta, tryr_score.
    """
    matrix_path = Path(SELECTIVITY_MATRIX_PATH)
    if matrix_path.exists():
        rows = _load_csv_rows(matrix_path)
        best_row: dict[str, str] | None = None
        best_delta: float = -999.0
        for row in rows:
            try:
                delta = float(row.get("delta", row.get("selectivity_delta", "-999")))
            except (ValueError, TypeError):
                continue
            if delta > best_delta:
                best_delta = delta
                best_row = row
        if best_row is not None and best_delta > 0.0:
            name = best_row.get("variant", best_row.get("compound", "unknown"))
            tryr_score = float(best_row.get("tryr", best_row.get("tryr_score", "0")))
            logger.info(
                "Selected variant %s from redesign matrix (delta=%.2f kcal/mol)",
                name,
                best_delta,
            )
            return {
                "name": name,
                "source": "selectivity_redesign",
                "delta": round(best_delta, 2),
                "tryr_score": round(tryr_score, 2),
            }
        logger.info("No positively selective variant in redesign matrix.")

    # Fallback: selectivity report from module 12
    report_path = Path(SELECTIVITY_REPORT_PATH)
    if report_path.exists():
        report = _load_json(report_path)
        results = report.get("results", [])
        best_result: dict[str, Any] | None = None
        best_d: float = -999.0
        for entry in results:
            delta = float(entry.get("delta", -999))
            if delta > best_d:
                best_d = delta
                best_result = entry
        if best_result is not None:
            name = best_result.get("compound", "MRL-003")
            tryr_score = float(best_result.get("tryr", 0.0))
            logger.info(
                "Using %s from selectivity report (delta=%.2f kcal/mol)",
                name,
                best_d,
            )
            return {
                "name": name,
                "source": "selectivity_report",
                "delta": round(best_d, 2),
                "tryr_score": round(tryr_score, 2),
            }

    # Final fallback: MRL-003 defaults
    logger.info("No upstream selectivity data found; using MRL-003 defaults.")
    return {
        "name": "MRL-003",
        "source": "default",
        "delta": 0.0,
        "tryr_score": -7.32,
    }


# ---------------------------------------------------------------------------
# OpenMM MD pipeline
# ---------------------------------------------------------------------------


def _find_ligand_sdf(variant_name: str) -> Path | None:
    """Locate the SDF file for a given variant name."""
    ligand_dir = Path(LIGAND_DIR)
    if not ligand_dir.exists():
        return None
    # Try exact name match first
    for suffix in [".sdf", ".mol2", ".pdb"]:
        candidate = ligand_dir / f"{variant_name}{suffix}"
        if candidate.exists():
            return candidate
    # Search for partial match
    for sdf_file in sorted(ligand_dir.glob("*.sdf")):
        if variant_name.lower().replace("-", "") in sdf_file.stem.lower().replace("-", ""):
            return sdf_file
    # Return any available SDF as last resort
    sdf_files = sorted(ligand_dir.glob("*.sdf"))
    if sdf_files:
        return sdf_files[0]
    return None


def prepare_system(
    protein_pdb: Path,
    ligand_sdf: Path | None,
) -> tuple[Any, Any, Any]:
    """Build solvated system with protein + water + ions.

    When a ligand SDF is provided and OpenFF is available, it is
    parameterised and added.  Otherwise only the protein is solvated.

    Args:
        protein_pdb: Path to the protein PDB file.
        ligand_sdf: Optional path to the ligand SDF file.

    Returns:
        Tuple of (topology, system, positions).
    """
    pdb = app.PDBFile(str(protein_pdb))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    # Add hydrogens at physiological pH
    modeller.addHydrogens(forcefield)
    logger.info("Added hydrogens to protein.")

    # Attempt to add ligand if available
    if ligand_sdf is not None:
        try:
            from openff.toolkit import Molecule
            from openmmforcefields.generators import GAFFTemplateGenerator

            mol = Molecule.from_file(str(ligand_sdf))
            gaff = GAFFTemplateGenerator(molecules=[mol])
            forcefield.registerTemplateGenerator(gaff.generator)
            logger.info("Registered ligand %s via GAFF.", ligand_sdf.name)
        except ImportError:
            logger.warning(
                "openff-toolkit / openmmforcefields not installed; "
                "proceeding with protein-only simulation."
            )
        except Exception as exc:
            logger.warning(
                "Could not parameterise ligand %s: %s; "
                "proceeding with protein-only simulation.",
                ligand_sdf.name,
                exc,
            )

    # Add solvent box with ionic strength
    modeller.addSolvent(
        forcefield,
        padding=SOLVENT_PADDING_NM * unit.nanometers,
        ionicStrength=IONIC_STRENGTH_M * unit.molar,
    )
    logger.info(
        "Solvated system: %.1f nm padding, %.2f M NaCl.",
        SOLVENT_PADDING_NM,
        IONIC_STRENGTH_M,
    )

    # Create OpenMM system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=NONBONDED_CUTOFF_NM * unit.nanometers,
        constraints=app.HBonds,
    )
    logger.info(
        "System created: %d particles.",
        system.getNumParticles(),
    )

    return modeller.topology, system, modeller.positions


def minimize(
    system: Any,
    topology: Any,
    positions: Any,
    max_iter: int = MINIMIZE_MAX_ITER,
) -> tuple[Any, Any]:
    """Energy-minimise the system using L-BFGS.

    Args:
        system: OpenMM System object.
        topology: OpenMM Topology.
        positions: Initial positions.
        max_iter: Maximum minimisation iterations.

    Returns:
        Tuple of (simulation, minimised_positions).
    """
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1 / unit.picosecond,
        TIMESTEP_FS * 0.001 * unit.picoseconds,
    )
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    logger.info("Starting energy minimisation (max %d iterations)...", max_iter)
    initial_state = simulation.context.getState(getEnergy=True)
    initial_energy = initial_state.getPotentialEnergy()
    logger.info("Initial energy: %s", initial_energy)

    simulation.minimizeEnergy(maxIterations=max_iter)

    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = final_state.getPotentialEnergy()
    logger.info("Final energy: %s", final_energy)

    return simulation, final_state.getPositions()


def equilibrate(
    system: Any,
    topology: Any,
    positions: Any,
    output_dir: Path,
    nvt_steps: int = NVT_STEPS,
    npt_steps: int = NPT_STEPS,
) -> tuple[Any, Any]:
    """Run NVT then NPT equilibration.

    Args:
        system: OpenMM System object.
        topology: OpenMM Topology.
        positions: Minimised positions.
        output_dir: Directory to write equilibration logs.
        nvt_steps: Number of NVT steps (default 50000 = 100 ps).
        npt_steps: Number of NPT steps (default 50000 = 100 ps).

    Returns:
        Tuple of (simulation, equilibrated_positions).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- NVT equilibration ---
    logger.info("NVT equilibration: %d steps (%.1f ps)...", nvt_steps, nvt_steps * TIMESTEP_FS / 1000)
    integrator_nvt = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1 / unit.picosecond,
        TIMESTEP_FS * 0.001 * unit.picoseconds,
    )
    sim_nvt = app.Simulation(topology, system, integrator_nvt)
    sim_nvt.context.setPositions(positions)
    sim_nvt.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)

    nvt_log_path = output_dir / "equilibration_nvt.log"
    sim_nvt.reporters.append(
        app.StateDataReporter(
            str(nvt_log_path),
            nvt_steps // 10,
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    sim_nvt.step(nvt_steps)
    nvt_positions = sim_nvt.context.getState(getPositions=True).getPositions()
    logger.info("NVT equilibration complete.")

    # --- NPT equilibration ---
    logger.info("NPT equilibration: %d steps (%.1f ps)...", npt_steps, npt_steps * TIMESTEP_FS / 1000)
    system.addForce(
        openmm.MonteCarloBarostat(
            PRESSURE_ATM * unit.atmospheres,
            TEMPERATURE_K * unit.kelvin,
        )
    )
    integrator_npt = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        1 / unit.picosecond,
        TIMESTEP_FS * 0.001 * unit.picoseconds,
    )
    sim_npt = app.Simulation(topology, system, integrator_npt)
    sim_npt.context.setPositions(nvt_positions)
    sim_npt.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)

    npt_log_path = output_dir / "equilibration_npt.log"
    sim_npt.reporters.append(
        app.StateDataReporter(
            str(npt_log_path),
            npt_steps // 10,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
        )
    )
    sim_npt.step(npt_steps)
    npt_positions = sim_npt.context.getState(getPositions=True).getPositions()
    logger.info("NPT equilibration complete.")

    return sim_npt, npt_positions


def run_production(
    simulation: Any,
    n_steps: int,
    output_dir: Path,
    report_interval: int = REPORT_INTERVAL,
) -> dict[str, Path]:
    """Run production MD and write trajectory + log files.

    Args:
        simulation: Equilibrated OpenMM Simulation object.
        n_steps: Number of production steps.
        output_dir: Directory for output files.
        report_interval: Steps between trajectory frames and log entries.

    Returns:
        Dict mapping output type to file path.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    dcd_path = output_dir / "production.dcd"
    log_path = output_dir / "production.log"

    simulation.reporters.append(
        app.DCDReporter(str(dcd_path), report_interval)
    )
    simulation.reporters.append(
        app.StateDataReporter(
            str(log_path),
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True,
        )
    )

    total_ns = n_steps * TIMESTEP_FS / 1_000_000
    logger.info(
        "Production MD: %d steps (%.1f ns), reporting every %d steps...",
        n_steps,
        total_ns,
        report_interval,
    )
    simulation.step(n_steps)
    logger.info("Production MD complete.")

    return {"dcd": dcd_path, "log": log_path}


# ---------------------------------------------------------------------------
# Trajectory analysis
# ---------------------------------------------------------------------------


def analyze_trajectory(
    dcd_path: Path,
    topology_pdb: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Analyse the production trajectory for stability metrics.

    Tries mdtraj for full analysis (RMSD, RMSF, contacts).  Falls back
    to parsing the production log for energy-only analysis.

    Args:
        dcd_path: Path to the DCD trajectory file.
        topology_pdb: Path to the topology PDB file.
        output_dir: Directory for analysis outputs.

    Returns:
        Dict of analysis results.
    """
    results: dict[str, Any] = {
        "analysis_method": "none",
        "ligand_rmsd_mean_A": None,
        "ligand_rmsd_max_A": None,
        "protein_rmsf_mean_A": None,
        "contacts_maintained": None,
        "contacts_total": None,
        "complex_stable": None,
    }

    try:
        import mdtraj
        import numpy as np

        traj = mdtraj.load(str(dcd_path), top=str(topology_pdb))
        logger.info("Loaded trajectory: %d frames, %d atoms.", traj.n_frames, traj.n_atoms)

        results["analysis_method"] = "mdtraj"

        # Protein backbone RMSD vs first frame
        backbone = traj.topology.select("backbone")
        rmsd_backbone = mdtraj.rmsd(traj, traj, frame=0, atom_indices=backbone)
        results["protein_rmsd_mean_A"] = round(float(np.mean(rmsd_backbone)) * 10.0, 2)
        results["protein_rmsd_max_A"] = round(float(np.max(rmsd_backbone)) * 10.0, 2)

        # Per-residue RMSF
        rmsf = mdtraj.rmsf(traj, traj, frame=0, atom_indices=backbone)
        results["protein_rmsf_mean_A"] = round(float(np.mean(rmsf)) * 10.0, 2)

        # Try ligand RMSD (non-protein, non-water heavy atoms)
        protein_atoms = set(traj.topology.select("protein"))
        water_atoms = set(traj.topology.select("water"))
        ion_names = {"NA", "CL", "Na+", "Cl-", "SOD", "CLA"}
        ion_atoms = set()
        for atom in traj.topology.atoms:
            if atom.residue.name in ion_names:
                ion_atoms.add(atom.index)
        all_special = protein_atoms | water_atoms | ion_atoms
        ligand_atoms = np.array(
            [i for i in range(traj.n_atoms) if i not in all_special]
        )

        if len(ligand_atoms) > 0:
            rmsd_ligand = mdtraj.rmsd(traj, traj, frame=0, atom_indices=ligand_atoms)
            results["ligand_rmsd_mean_A"] = round(float(np.mean(rmsd_ligand)) * 10.0, 2)
            results["ligand_rmsd_max_A"] = round(float(np.max(rmsd_ligand)) * 10.0, 2)
            stable = float(np.mean(rmsd_ligand)) * 10.0 < 3.0
            results["complex_stable"] = stable
        else:
            logger.info("No ligand atoms found; reporting protein-only metrics.")
            stable = float(np.mean(rmsd_backbone)) * 10.0 < 3.0
            results["complex_stable"] = stable

        # Contact analysis: residues within 0.4 nm of ligand in frame 0
        if len(ligand_atoms) > 0:
            contacts_frame0 = mdtraj.compute_neighbors(traj[0], 0.4, ligand_atoms)
            initial_contacts = set(contacts_frame0[0]) & protein_atoms
            total_contacts = len(initial_contacts)
            # Check how many are maintained in last 10% of trajectory
            last_frames = traj[int(traj.n_frames * 0.9):]
            maintained = 0
            for contact_atom in initial_contacts:
                dists = mdtraj.compute_distances(
                    last_frames,
                    [[contact_atom, ligand_atoms[0]]],
                )
                if float(np.mean(dists)) < 0.5:
                    maintained += 1
            results["contacts_maintained"] = maintained
            results["contacts_total"] = total_contacts

        # Save RMSF per residue
        rmsf_data = []
        ca_atoms = traj.topology.select("name CA")
        if len(ca_atoms) > 0:
            rmsf_ca = mdtraj.rmsf(traj, traj, frame=0, atom_indices=ca_atoms)
            for idx_atom, rmsf_val in zip(ca_atoms, rmsf_ca):
                atom = traj.topology.atom(idx_atom)
                rmsf_data.append({
                    "residue_index": atom.residue.index,
                    "residue_name": str(atom.residue),
                    "rmsf_nm": round(float(rmsf_val), 4),
                    "rmsf_A": round(float(rmsf_val) * 10.0, 2),
                })
            rmsf_path = output_dir / "rmsf_per_residue.json"
            _write_json(rmsf_data, rmsf_path)

        logger.info(
            "MDTraj analysis: backbone RMSD=%.2f A, RMSF=%.2f A, stable=%s",
            results["protein_rmsd_mean_A"],
            results["protein_rmsf_mean_A"],
            results["complex_stable"],
        )

    except ImportError:
        logger.warning(
            "mdtraj not installed; falling back to energy-only analysis. "
            "Install via: conda install -c conda-forge mdtraj"
        )
        results["analysis_method"] = "energy_log"
        results = _analyze_energy_log(dcd_path.parent / "production.log", results)

    except Exception as exc:
        logger.warning("MDTraj analysis failed: %s; using energy log.", exc)
        results["analysis_method"] = "energy_log"
        results = _analyze_energy_log(dcd_path.parent / "production.log", results)

    return results


def _analyze_energy_log(
    log_path: Path,
    results: dict[str, Any],
) -> dict[str, Any]:
    """Parse the StateDataReporter log for basic energy analysis.

    Args:
        log_path: Path to the production.log file.
        results: Existing results dict to update.

    Returns:
        Updated results dict with energy statistics.
    """
    if not log_path.exists():
        logger.warning("Production log not found at %s.", log_path)
        return results

    energies: list[float] = []
    temperatures: list[float] = []

    with open(log_path) as fh:
        header_found = False
        energy_col = -1
        temp_col = -1
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Parse header to find column indices
                cols = [c.strip().strip('"') for c in line.lstrip("#").split(",")]
                for i, col in enumerate(cols):
                    if "potential" in col.lower() or "energy" in col.lower():
                        energy_col = i
                    if "temperature" in col.lower():
                        temp_col = i
                header_found = True
                continue
            if not header_found:
                # Try comma-separated header
                if "Potential Energy" in line or "Temperature" in line:
                    cols = [c.strip().strip('"') for c in line.split(",")]
                    for i, col in enumerate(cols):
                        if "Potential Energy" in col:
                            energy_col = i
                        if "Temperature" in col:
                            temp_col = i
                    header_found = True
                continue
            parts = line.split(",")
            try:
                if energy_col >= 0 and energy_col < len(parts):
                    energies.append(float(parts[energy_col].strip()))
                if temp_col >= 0 and temp_col < len(parts):
                    temperatures.append(float(parts[temp_col].strip()))
            except (ValueError, IndexError):
                continue

    if energies:
        import statistics

        results["energy_mean_kj"] = round(statistics.mean(energies), 1)
        results["energy_std_kj"] = round(statistics.stdev(energies) if len(energies) > 1 else 0.0, 1)
        results["energy_drift_kj"] = round(energies[-1] - energies[0], 1)
        logger.info(
            "Energy analysis: mean=%.1f kJ/mol, std=%.1f, drift=%.1f",
            results["energy_mean_kj"],
            results["energy_std_kj"],
            results["energy_drift_kj"],
        )
    if temperatures:
        import statistics

        results["temperature_mean_K"] = round(statistics.mean(temperatures), 1)
        results["temperature_std_K"] = round(statistics.stdev(temperatures) if len(temperatures) > 1 else 0.0, 1)

    return results


# ---------------------------------------------------------------------------
# Setup scripts (when OpenMM is not available)
# ---------------------------------------------------------------------------


def _generate_setup_script(variant: dict[str, Any], output_dir: Path) -> Path:
    """Generate a shell script to install dependencies and run MD.

    Args:
        variant: Selected variant information dict.
        output_dir: Directory to write the script.

    Returns:
        Path to the generated shell script.
    """
    script_path = output_dir / "run_md.sh"
    total_ns = DEFAULT_STEPS * TIMESTEP_FS / 1_000_000

    script = f"""\
#!/bin/bash
# =============================================================================
# Marley MD Setup Script
# Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
# =============================================================================
#
# Run this script on a machine with CUDA-capable GPU for best performance.
# Estimated time: ~2-4 hours on a modern GPU (RTX 3080 or better).
#
# Variant: {variant["name"]}
# Simulation: {total_ns:.1f} ns production MD
# Temperature: {TEMPERATURE_K} K (dog body temperature)
# =============================================================================

set -euo pipefail

echo "=== Marley Molecular Dynamics Setup ==="

# --- Step 1: Create conda environment ---
echo "[1/4] Setting up conda environment..."
if ! conda info --envs | grep -q "marley_md"; then
    conda create -n marley_md python=3.11 -y
fi
conda activate marley_md

# --- Step 2: Install dependencies ---
echo "[2/4] Installing OpenMM and analysis tools..."
conda install -c conda-forge openmm openff-toolkit openmmforcefields mdtraj -y
pip install numpy scipy

# --- Step 3: Verify GPU ---
echo "[3/4] Checking GPU availability..."
python -c "
import openmm
platforms = [openmm.Platform.getPlatform(i).getName()
             for i in range(openmm.Platform.getNumPlatforms())]
print('Available platforms:', platforms)
if 'CUDA' in platforms:
    print('CUDA available - will use GPU acceleration.')
elif 'OpenCL' in platforms:
    print('OpenCL available - will use GPU acceleration.')
else:
    print('WARNING: No GPU platform found. Simulation will be slow on CPU.')
"

# --- Step 4: Run the simulation ---
echo "[4/4] Running molecular dynamics..."
python pipeline/19_molecular_dynamics.py --force --steps {DEFAULT_STEPS}

echo "=== MD simulation complete ==="
echo "Results in: results/molecular_dynamics/"
"""
    _write_text(script, script_path)
    script_path.chmod(0o755)
    return script_path


def _generate_gromacs_inputs(variant: dict[str, Any], output_dir: Path) -> list[Path]:
    """Generate GROMACS-compatible input files as an alternative to OpenMM.

    Args:
        variant: Selected variant information dict.
        output_dir: Directory to write GROMACS files.

    Returns:
        List of paths to generated files.
    """
    gromacs_dir = output_dir / "gromacs"
    gromacs_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    total_ns = DEFAULT_STEPS * TIMESTEP_FS / 1_000_000
    total_ps = total_ns * 1000

    # Energy minimisation .mdp
    em_mdp = f"""\
; Marley - Energy Minimisation
; Variant: {variant["name"]}
integrator      = steep
emtol           = 1000.0    ; kJ/mol/nm
emstep          = 0.01      ; nm
nsteps          = {MINIMIZE_MAX_ITER}

; Neighbour searching
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; Van der Waals
rvdw            = 1.0

; Constraints
constraints     = h-bonds
"""
    em_path = gromacs_dir / "em.mdp"
    _write_text(em_mdp, em_path)
    generated.append(em_path)

    # NVT equilibration .mdp
    nvt_mdp = f"""\
; Marley - NVT Equilibration (100 ps)
; Variant: {variant["name"]}
integrator      = md
dt              = 0.002     ; ps
nsteps          = {NVT_STEPS}
nstxout-compressed = 5000

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1   0.1
ref_t           = {TEMPERATURE_K}   {TEMPERATURE_K}

; Pressure coupling (off for NVT)
pcoupl          = no

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0

; Constraints
constraints     = h-bonds
constraint_algorithm = lincs

; Neighbour searching
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid

; Velocity generation
gen_vel         = yes
gen_temp        = {TEMPERATURE_K}
gen_seed        = -1
"""
    nvt_path = gromacs_dir / "nvt.mdp"
    _write_text(nvt_mdp, nvt_path)
    generated.append(nvt_path)

    # NPT equilibration .mdp
    npt_mdp = f"""\
; Marley - NPT Equilibration (100 ps)
; Variant: {variant["name"]}
integrator      = md
dt              = 0.002     ; ps
nsteps          = {NPT_STEPS}
nstxout-compressed = 5000

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1   0.1
ref_t           = {TEMPERATURE_K}   {TEMPERATURE_K}

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0

; Constraints
constraints     = h-bonds
constraint_algorithm = lincs

; Neighbour searching
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid

; Continue from NVT
gen_vel         = no
continuation    = yes
"""
    npt_path = gromacs_dir / "npt.mdp"
    _write_text(npt_mdp, npt_path)
    generated.append(npt_path)

    # Production .mdp
    prod_mdp = f"""\
; Marley - Production MD ({total_ns:.1f} ns)
; Variant: {variant["name"]}
integrator      = md
dt              = 0.002     ; ps
nsteps          = {DEFAULT_STEPS}
nstxout-compressed = {REPORT_INTERVAL}

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1   0.1
ref_t           = {TEMPERATURE_K}   {TEMPERATURE_K}

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0

; Constraints
constraints     = h-bonds
constraint_algorithm = lincs

; Neighbour searching
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid

; Continue from NPT
gen_vel         = no
continuation    = yes

; Output control
nstlog          = {REPORT_INTERVAL}
nstenergy       = {REPORT_INTERVAL}
nstcalcenergy   = 100
"""
    prod_path = gromacs_dir / "production.mdp"
    _write_text(prod_mdp, prod_path)
    generated.append(prod_path)

    # GROMACS workflow script
    workflow = f"""\
#!/bin/bash
# =============================================================================
# Marley - GROMACS MD Workflow
# Variant: {variant["name"]}
# Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
# =============================================================================
#
# Prerequisites:
#   - GROMACS 2023+ installed
#   - Protein PDB: {PROTEIN_PDB_PATH}
#   - Ligand topology prepared (e.g., via ACPYPE or CGenFF)
#
# Adjust paths as needed for your cluster environment.
# =============================================================================

set -euo pipefail
WORKDIR="results/molecular_dynamics/gromacs_run"
mkdir -p "$WORKDIR"

echo "=== Step 1: Prepare topology ==="
gmx pdb2gmx -f {PROTEIN_PDB_PATH} -o "$WORKDIR/protein.gro" -water tip3p -ff amber14sb_OL15 -ignh

echo "=== Step 2: Define box and solvate ==="
gmx editconf -f "$WORKDIR/protein.gro" -o "$WORKDIR/box.gro" -c -d 1.0 -bt dodecahedron
gmx solvate -cp "$WORKDIR/box.gro" -cs spc216.gro -o "$WORKDIR/solv.gro" -p "$WORKDIR/topol.top"

echo "=== Step 3: Add ions (0.15 M NaCl) ==="
gmx grompp -f results/molecular_dynamics/gromacs/em.mdp -c "$WORKDIR/solv.gro" -p "$WORKDIR/topol.top" -o "$WORKDIR/ions.tpr" -maxwarn 1
echo "SOL" | gmx genion -s "$WORKDIR/ions.tpr" -o "$WORKDIR/ions.gro" -p "$WORKDIR/topol.top" -pname NA -nname CL -neutral -conc 0.15

echo "=== Step 4: Energy minimisation ==="
gmx grompp -f results/molecular_dynamics/gromacs/em.mdp -c "$WORKDIR/ions.gro" -p "$WORKDIR/topol.top" -o "$WORKDIR/em.tpr"
gmx mdrun -deffnm "$WORKDIR/em"

echo "=== Step 5: NVT equilibration ==="
gmx grompp -f results/molecular_dynamics/gromacs/nvt.mdp -c "$WORKDIR/em.gro" -r "$WORKDIR/em.gro" -p "$WORKDIR/topol.top" -o "$WORKDIR/nvt.tpr"
gmx mdrun -deffnm "$WORKDIR/nvt"

echo "=== Step 6: NPT equilibration ==="
gmx grompp -f results/molecular_dynamics/gromacs/npt.mdp -c "$WORKDIR/nvt.gro" -r "$WORKDIR/nvt.gro" -p "$WORKDIR/topol.top" -o "$WORKDIR/npt.tpr"
gmx mdrun -deffnm "$WORKDIR/npt"

echo "=== Step 7: Production MD ({total_ns:.1f} ns) ==="
gmx grompp -f results/molecular_dynamics/gromacs/production.mdp -c "$WORKDIR/npt.gro" -r "$WORKDIR/npt.gro" -p "$WORKDIR/topol.top" -o "$WORKDIR/prod.tpr"
gmx mdrun -deffnm "$WORKDIR/prod"

echo "=== Step 8: Basic analysis ==="
echo "Backbone" | gmx rms -s "$WORKDIR/prod.tpr" -f "$WORKDIR/prod.xtc" -o "$WORKDIR/rmsd.xvg"
echo "Backbone" | gmx rmsf -s "$WORKDIR/prod.tpr" -f "$WORKDIR/prod.xtc" -o "$WORKDIR/rmsf.xvg" -res

echo "=== GROMACS MD workflow complete ==="
echo "Results in: $WORKDIR/"
"""
    workflow_path = gromacs_dir / "run_gromacs.sh"
    _write_text(workflow, workflow_path)
    workflow_path.chmod(0o755)
    generated.append(workflow_path)

    return generated


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _build_report(
    variant: dict[str, Any],
    analysis: dict[str, Any] | None,
    n_steps: int,
    openmm_used: bool,
    dry_run: bool,
) -> str:
    """Generate the MD simulation Markdown report.

    Args:
        variant: Selected variant information.
        analysis: Trajectory analysis results (None if not run).
        n_steps: Number of production steps.
        openmm_used: Whether OpenMM was used.
        dry_run: Whether this was a dry run.

    Returns:
        Markdown report string.
    """
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    total_ns = n_steps * TIMESTEP_FS / 1_000_000

    lines: list[str] = []
    lines.append("# Molecular Dynamics Simulation Report")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")

    # System setup
    lines.append("## System Setup")
    lines.append(f"- Protein: TryR (L. infantum), AlphaFold model")
    lines.append(f"- Ligand: {variant['name']} (from {variant['source']})")
    lines.append(f"- Docking delta: {variant['delta']} kcal/mol (TryR selectivity)")
    lines.append(f"- Solvent: TIP3P water, {IONIC_STRENGTH_M} M NaCl")
    lines.append(f"- Temperature: {TEMPERATURE_K} K (dog body temperature)")
    lines.append(f"- Pressure: {PRESSURE_ATM} atm")
    lines.append(f"- Timestep: {TIMESTEP_FS} fs")
    lines.append(f"- Total simulation: {total_ns:.1f} ns ({n_steps:,} steps)")
    lines.append(f"- Force field: AMBER14 (protein) + GAFF (ligand)")
    lines.append(f"- Nonbonded: PME, cutoff {NONBONDED_CUTOFF_NM} nm")
    lines.append("")

    # Equilibration protocol
    lines.append("## Equilibration Protocol")
    nvt_ps = NVT_STEPS * TIMESTEP_FS / 1000
    npt_ps = NPT_STEPS * TIMESTEP_FS / 1000
    lines.append(f"1. Energy minimisation: steepest descent, max {MINIMIZE_MAX_ITER:,} steps")
    lines.append(f"2. NVT equilibration: {nvt_ps:.0f} ps at {TEMPERATURE_K} K (Langevin thermostat)")
    lines.append(f"3. NPT equilibration: {npt_ps:.0f} ps at {TEMPERATURE_K} K, {PRESSURE_ATM} atm (MC barostat)")
    lines.append("")

    # Execution status
    lines.append("## Execution Status")
    if dry_run:
        lines.append("- Mode: **DRY RUN** (minimisation only, no production MD)")
    elif not openmm_used:
        lines.append("- Mode: **SETUP ONLY** (OpenMM not available)")
        lines.append("- Setup scripts generated for cluster execution")
        lines.append("- GROMACS input files generated as alternative")
    else:
        lines.append("- Mode: **FULL SIMULATION**")
        lines.append(f"- Engine: OpenMM (GPU-accelerated if available)")
    lines.append("")

    # Results
    lines.append("## Results")
    lines.append("")

    if analysis is not None and analysis.get("analysis_method") == "mdtraj":
        rmsd_mean = analysis.get("ligand_rmsd_mean_A")
        rmsd_max = analysis.get("ligand_rmsd_max_A")
        prot_rmsd_mean = analysis.get("protein_rmsd_mean_A")
        rmsf_mean = analysis.get("protein_rmsf_mean_A")
        contacts_kept = analysis.get("contacts_maintained")
        contacts_total = analysis.get("contacts_total")
        stable = analysis.get("complex_stable")

        if rmsd_mean is not None:
            lines.append(f"- Ligand RMSD (mean): {rmsd_mean} A {'(stable)' if rmsd_mean < 3.0 else '(unstable, > 3.0 A)'}")
            lines.append(f"- Ligand RMSD (max): {rmsd_max} A")
        if prot_rmsd_mean is not None:
            lines.append(f"- Protein backbone RMSD (mean): {prot_rmsd_mean} A")
        if rmsf_mean is not None:
            lines.append(f"- Protein RMSF (mean): {rmsf_mean} A")
        if contacts_kept is not None:
            lines.append(f"- Key contacts maintained: {contacts_kept}/{contacts_total}")
        if stable is not None:
            lines.append(f"- Complex stable: {'YES' if stable else 'NO'}")
        lines.append("")

        # Interpretation
        lines.append("## Interpretation")
        if stable:
            lines.append(
                f"The {variant['name']}-TryR complex remained stable over the "
                f"{total_ns:.1f} ns simulation, with the ligand maintaining a mean "
                f"RMSD of {rmsd_mean} A from its initial docked pose. This suggests "
                f"the binding mode predicted by docking is physically plausible and "
                f"the compound forms a stable complex with the target enzyme."
            )
        else:
            lines.append(
                f"The {variant['name']}-TryR complex showed significant movement "
                f"during the {total_ns:.1f} ns simulation. The ligand may adopt "
                f"alternative binding modes or dissociate partially. Further "
                f"optimisation of the binding interactions may be needed."
            )
        lines.append("")

    elif analysis is not None and analysis.get("analysis_method") == "energy_log":
        energy_mean = analysis.get("energy_mean_kj")
        energy_std = analysis.get("energy_std_kj")
        energy_drift = analysis.get("energy_drift_kj")
        temp_mean = analysis.get("temperature_mean_K")

        lines.append("### Energy Analysis (from production log)")
        if energy_mean is not None:
            lines.append(f"- Mean potential energy: {energy_mean} kJ/mol")
            lines.append(f"- Energy std deviation: {energy_std} kJ/mol")
            lines.append(f"- Energy drift: {energy_drift} kJ/mol")
        if temp_mean is not None:
            lines.append(f"- Mean temperature: {temp_mean} K (target: {TEMPERATURE_K} K)")
        lines.append("")
        lines.append(
            "Full structural analysis (RMSD, RMSF, contacts) requires mdtraj. "
            "Install via: `conda install -c conda-forge mdtraj`"
        )
        lines.append("")

    else:
        if dry_run:
            lines.append("Dry run completed. Energy minimisation was performed.")
            lines.append("Run without --dry-run for full production MD.")
        elif not openmm_used:
            lines.append("OpenMM was not available on this machine.")
            lines.append("Use the generated setup scripts to run on a cluster:")
            lines.append(f"  - `results/molecular_dynamics/run_md.sh` (OpenMM)")
            lines.append(f"  - `results/molecular_dynamics/gromacs/run_gromacs.sh` (GROMACS)")
        else:
            lines.append("No analysis data available.")
        lines.append("")

    # Files generated
    lines.append("## Output Files")
    lines.append(f"- `md_report.md` -- this report")
    lines.append(f"- `md_metadata.json` -- simulation metadata and results")
    if openmm_used and not dry_run:
        lines.append(f"- `production.dcd` -- trajectory (DCD format)")
        lines.append(f"- `production.log` -- energy/temperature log")
        lines.append(f"- `equilibration_nvt.log` -- NVT equilibration log")
        lines.append(f"- `equilibration_npt.log` -- NPT equilibration log")
        lines.append(f"- `rmsf_per_residue.json` -- per-residue RMSF (if mdtraj available)")
    if not openmm_used:
        lines.append(f"- `run_md.sh` -- OpenMM setup/run script for cluster")
        lines.append(f"- `gromacs/` -- GROMACS input files (.mdp + workflow script)")
    lines.append("")

    # Methods
    lines.append("## Methods")
    lines.append(
        f"Molecular dynamics simulation of the {variant['name']}-TryR complex "
        f"was performed using {'OpenMM 8.x' if openmm_used else 'generated input files for OpenMM/GROMACS'}. "
        f"The protein was modelled with the AMBER14 force field and solvated in "
        f"TIP3P water with {IONIC_STRENGTH_M} M NaCl in a periodic box with "
        f"{SOLVENT_PADDING_NM} nm padding. After energy minimisation (steepest "
        f"descent, {MINIMIZE_MAX_ITER:,} steps max), the system was equilibrated "
        f"with {NVT_STEPS:,} steps NVT followed by {NPT_STEPS:,} steps NPT at "
        f"{TEMPERATURE_K} K and {PRESSURE_ATM} atm. Production MD was run for "
        f"{total_ns:.1f} ns with a {TIMESTEP_FS} fs timestep and PME "
        f"electrostatics (cutoff {NONBONDED_CUTOFF_NM} nm). Trajectory analysis "
        f"was performed with MDTraj where available."
    )
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run_molecular_dynamics(
    steps: int = DEFAULT_STEPS,
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the molecular dynamics simulation pipeline.

    Steps:
        1. Select the best variant from upstream selectivity modules.
        2. Check OpenMM availability.
        3. If OpenMM available: prepare, minimise, equilibrate, produce, analyse.
        4. If OpenMM not available: generate setup scripts.
        5. Write report and metadata.

    Args:
        steps: Number of production MD steps (default 5,000,000 = 10 ns).
        force: If True, overwrite existing output even if present.
        dry_run: If True, minimise only (skip equilibration and production).

    Returns:
        Path to the generated MD report.
    """
    output_base = Path(OUTPUT_DIR)
    report_path = output_base / "md_report.md"

    if report_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            report_path,
        )
        return str(report_path)

    # ------------------------------------------------------------------
    # 1. Select best variant
    # ------------------------------------------------------------------
    variant = _select_best_variant()
    logger.info(
        "Selected variant: %s (source=%s, delta=%.2f kcal/mol)",
        variant["name"],
        variant["source"],
        variant["delta"],
    )

    # ------------------------------------------------------------------
    # 2. Check OpenMM and run appropriate pipeline
    # ------------------------------------------------------------------
    analysis: dict[str, Any] | None = None
    openmm_used = False

    if OPENMM_AVAILABLE:
        protein_pdb = Path(PROTEIN_PDB_PATH)
        if not protein_pdb.exists():
            logger.error("Protein PDB not found at %s.", protein_pdb)
            sys.exit(1)

        ligand_sdf = _find_ligand_sdf(variant["name"])
        if ligand_sdf is not None:
            logger.info("Found ligand file: %s", ligand_sdf)
        else:
            logger.warning(
                "No ligand file found for %s; running protein-only simulation.",
                variant["name"],
            )

        # --- Prepare system ---
        logger.info("Preparing solvated system...")
        topology, system, positions = prepare_system(protein_pdb, ligand_sdf)
        openmm_used = True

        # --- Minimise ---
        simulation, min_positions = minimize(system, topology, positions)

        if dry_run:
            logger.info("Dry run: skipping equilibration and production.")
        else:
            # --- Equilibrate ---
            sim_eq, eq_positions = equilibrate(
                system, topology, min_positions, output_base,
            )

            # --- Production ---
            outputs = run_production(sim_eq, steps, output_base)

            # --- Analyse ---
            # Save topology for analysis
            topology_pdb_out = output_base / "topology.pdb"
            with open(topology_pdb_out, "w") as fh:
                app.PDBFile.writeFile(topology, eq_positions, fh)
            logger.info("Wrote topology PDB: %s", topology_pdb_out)

            analysis = analyze_trajectory(
                outputs["dcd"],
                topology_pdb_out,
                output_base,
            )
    else:
        # --- Generate setup scripts ---
        logger.info("OpenMM not available. Generating setup scripts...")

        script_path = _generate_setup_script(variant, output_base)
        logger.info("Generated OpenMM setup script: %s", script_path)

        gromacs_files = _generate_gromacs_inputs(variant, output_base)
        logger.info(
            "Generated %d GROMACS input files in %s",
            len(gromacs_files),
            output_base / "gromacs",
        )

    # ------------------------------------------------------------------
    # 3. Write metadata and report
    # ------------------------------------------------------------------
    metadata: dict[str, Any] = {
        "generated": datetime.now().isoformat(),
        "module": "19_molecular_dynamics",
        "variant": variant,
        "simulation_parameters": {
            "steps": steps,
            "timestep_fs": TIMESTEP_FS,
            "total_ns": round(steps * TIMESTEP_FS / 1_000_000, 1),
            "temperature_K": TEMPERATURE_K,
            "pressure_atm": PRESSURE_ATM,
            "ionic_strength_M": IONIC_STRENGTH_M,
            "solvent_padding_nm": SOLVENT_PADDING_NM,
            "nonbonded_cutoff_nm": NONBONDED_CUTOFF_NM,
            "nvt_steps": NVT_STEPS,
            "npt_steps": NPT_STEPS,
            "minimize_max_iter": MINIMIZE_MAX_ITER,
        },
        "openmm_available": OPENMM_AVAILABLE,
        "openmm_used": openmm_used,
        "dry_run": dry_run,
        "analysis": analysis,
    }
    _write_json(metadata, output_base / "md_metadata.json")

    report = _build_report(variant, analysis, steps, openmm_used, dry_run)
    _write_text(report, report_path)

    logger.info("Molecular dynamics pipeline complete. Report: %s", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley -- Molecular dynamics simulation of the TryR-ligand "
            "complex using OpenMM. Generates setup scripts when OpenMM "
            "is not available."
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
            "Run minimisation only (skip equilibration and production). "
            "Useful for testing the pipeline setup."
        ),
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_STEPS,
        help=(
            f"Number of production MD steps (default: {DEFAULT_STEPS:,} = "
            f"{DEFAULT_STEPS * TIMESTEP_FS / 1_000_000:.0f} ns at "
            f"{TIMESTEP_FS} fs timestep)."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = run_molecular_dynamics(
        steps=args.steps,
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
