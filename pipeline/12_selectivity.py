"""Selectivity analysis: compare MRL-003 binding to parasite TryR vs human GR.

Downloads the human glutathione reductase (GR) structure from AlphaFold,
converts it to PDBQT, docks MRL-003 and reference compounds against both
parasite TryR and human GR, and produces a selectivity report showing
differential binding affinities.

A selective drug candidate should bind more tightly to the parasite enzyme
(more negative kcal/mol) than to the human homolog, with a delta of at
least 1.5 kcal/mol.

Usage:
    python -m pipeline.12_selectivity
    python -m pipeline.12_selectivity --force
    python -m pipeline.12_selectivity --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("selectivity")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Paths
STRUCTURES_DIR: Final[Path] = Path("data/structures")
OUTPUT_DIR: Final[Path] = Path("results/selectivity")
HUMAN_GR_PDB: Final[Path] = STRUCTURES_DIR / "human_GR_P00390.pdb"
HUMAN_GR_PDBQT: Final[Path] = STRUCTURES_DIR / "human_GR_P00390_receptor.pdbqt"
TRYR_RECEPTOR: Final[Path] = STRUCTURES_DIR / "TryR_Q4Q457_receptor.pdbqt"

# AlphaFold download URL for human glutathione reductase (UniProt P00390)
ALPHAFOLD_URL: Final[str] = (
    "https://alphafold.ebi.ac.uk/files/AF-P00390-F1-model_v6.pdb"
)

# Tool paths (project virtual environment)
OBABEL_BIN: Final[str] = ".venv/bin/obabel"
VINA_BIN: Final[str] = ".venv/bin/vina"

# Docking parameters
EXHAUSTIVENESS: Final[int] = 8
NUM_MODES: Final[int] = 5
VINA_TIMEOUT: Final[int] = 600

# TryR known binding site (from drug_targets/08_docking.py)
TRYR_GRID: Final[dict[str, float]] = {
    "center_x": 11.0,
    "center_y": 24.0,
    "center_z": 15.0,
    "size_x": 25.0,
    "size_y": 25.0,
    "size_z": 25.0,
}

# Human GR uses blind docking (centroid + 40x40x40 box)
BLIND_BOX_SIZE: Final[float] = 40.0

# Selectivity threshold
SELECTIVITY_THRESHOLD: Final[float] = 1.5  # kcal/mol

# Compound library
COMPOUNDS: Final[dict[str, dict[str, str]]] = {
    "MRL-003": {
        "smiles": (
            "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)"
            "NC(CCC(=O)N)C(=O)N"
        ),
        "description": "Lead compound designed against L. infantum TryR",
    },
    "Pemetrexed": {
        "smiles": (
            "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)"
            "NC(CCC(=O)O)C(=O)O"
        ),
        "description": "Antifolate reference (structural analog of MRL-003)",
    },
    "Menadione": {
        "smiles": "O=C1C=CC(=O)C2=CC=CC=C12",
        "description": "Known GR inhibitor (positive control for human GR binding)",
    },
}

# Dry-run placeholder affinities (kcal/mol)
DRY_RUN_AFFINITIES: Final[dict[str, dict[str, float]]] = {
    "MRL-003": {"TryR": -7.8, "HumanGR": -5.2},
    "Pemetrexed": {"TryR": -7.1, "HumanGR": -6.3},
    "Menadione": {"TryR": -5.5, "HumanGR": -6.8},
}

# HTTP download settings
DOWNLOAD_TIMEOUT: Final[int] = 60
MAX_RETRIES: Final[int] = 3
RETRY_BACKOFF: Final[float] = 2.0


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------


def _download_file(url: str, dest: Path) -> None:
    """Download a file from *url* to *dest* with retry and backoff.

    Args:
        url: Remote URL to fetch.
        dest: Local path to write the downloaded content.

    Raises:
        requests.HTTPError: If the download fails after all retries.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            logger.info(
                "Downloading %s (attempt %d/%d) ...", url, attempt, MAX_RETRIES
            )
            resp = requests.get(url, timeout=DOWNLOAD_TIMEOUT)
            resp.raise_for_status()
            dest.write_bytes(resp.content)
            logger.info("Saved %s (%d bytes).", dest, len(resp.content))
            return
        except requests.RequestException as exc:
            if attempt == MAX_RETRIES:
                logger.error("Download failed after %d attempts: %s", MAX_RETRIES, exc)
                raise
            wait = RETRY_BACKOFF ** attempt
            logger.warning(
                "Attempt %d failed (%s). Retrying in %.1fs ...",
                attempt,
                exc,
                wait,
            )
            time.sleep(wait)


# ---------------------------------------------------------------------------
# Structure preparation
# ---------------------------------------------------------------------------


def _pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> None:
    """Convert a PDB file to PDBQT using Open Babel.

    Uses the ``-xrh`` flag to add hydrogens and generate a receptor-ready
    PDBQT file.

    Args:
        pdb_path: Input PDB file.
        pdbqt_path: Output PDBQT file.

    Raises:
        RuntimeError: If Open Babel conversion fails.
    """
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        OBABEL_BIN,
        "-ipdb", str(pdb_path),
        "-opdbqt",
        "-O", str(pdbqt_path),
        "-xrh",
    ]
    logger.info("Converting %s -> %s", pdb_path.name, pdbqt_path.name)
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=120,
        )
        if result.stderr:
            logger.debug("obabel stderr: %s", result.stderr.strip())
        logger.info("Wrote receptor PDBQT: %s", pdbqt_path)
    except FileNotFoundError:
        msg = f"Open Babel not found at {OBABEL_BIN}. Install with: pip install openbabel-wheel"
        raise RuntimeError(msg)
    except subprocess.CalledProcessError as exc:
        msg = f"Open Babel conversion failed: {exc.stderr}"
        raise RuntimeError(msg)


def _smiles_to_pdbqt(smiles: str, name: str, output_dir: Path) -> Path:
    """Convert a SMILES string to a 3D PDBQT ligand file via Open Babel.

    Generates 3D coordinates with ``--gen3d`` and writes a PDBQT suitable
    for AutoDock Vina docking.

    Args:
        smiles: SMILES string for the compound.
        name: Compound name (used for the output filename).
        output_dir: Directory to write the PDBQT file.

    Returns:
        Path to the generated PDBQT file.

    Raises:
        RuntimeError: If Open Babel conversion fails.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{name}.pdbqt"

    cmd = [
        OBABEL_BIN,
        "-:" + smiles,
        "-opdbqt",
        "-O", str(output_path),
        "--gen3d",
        "-h",
    ]
    logger.info("Generating ligand PDBQT for %s ...", name)
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=120,
        )
        if result.stderr:
            logger.debug("obabel stderr: %s", result.stderr.strip())
        logger.info("Wrote ligand PDBQT: %s", output_path)
    except FileNotFoundError:
        msg = f"Open Babel not found at {OBABEL_BIN}. Install with: pip install openbabel-wheel"
        raise RuntimeError(msg)
    except subprocess.CalledProcessError as exc:
        msg = f"Ligand PDBQT generation failed for {name}: {exc.stderr}"
        raise RuntimeError(msg)

    return output_path


# ---------------------------------------------------------------------------
# Grid box helpers
# ---------------------------------------------------------------------------


def _compute_centroid(pdbqt_path: Path) -> dict[str, float]:
    """Parse PDBQT ATOM records and compute the protein centroid.

    Returns a blind docking grid box centered on the centroid with a
    uniform size of 40 Angstroms per axis.

    Args:
        pdbqt_path: Path to a receptor PDBQT file.

    Returns:
        Dictionary with center_x/y/z and size_x/y/z keys.

    Raises:
        ValueError: If no ATOM records are found.
    """
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []

    with open(pdbqt_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except (ValueError, IndexError):
                    continue

    if not xs:
        msg = f"No ATOM records found in {pdbqt_path}"
        raise ValueError(msg)

    return {
        "center_x": round(sum(xs) / len(xs), 3),
        "center_y": round(sum(ys) / len(ys), 3),
        "center_z": round(sum(zs) / len(zs), 3),
        "size_x": BLIND_BOX_SIZE,
        "size_y": BLIND_BOX_SIZE,
        "size_z": BLIND_BOX_SIZE,
    }


# ---------------------------------------------------------------------------
# Vina docking
# ---------------------------------------------------------------------------


def _parse_vina_output(stdout: str) -> float:
    """Extract the best-pose binding affinity from Vina CLI output.

    Args:
        stdout: Standard output captured from the vina command.

    Returns:
        Binding affinity in kcal/mol for the top-ranked pose.

    Raises:
        RuntimeError: If no scoring lines are found.
    """
    pattern = re.compile(
        r"^\s*1\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)",
    )
    for line in stdout.splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1))

    msg = "Could not parse docking scores from Vina CLI output."
    raise RuntimeError(msg)


def _run_vina(
    receptor_path: Path,
    ligand_path: Path,
    output_path: Path,
    grid_box: dict[str, float],
) -> float:
    """Execute a single AutoDock Vina docking job.

    Tries the Vina Python bindings first, then falls back to the CLI.

    Args:
        receptor_path: Path to the receptor PDBQT file.
        ligand_path: Path to the ligand PDBQT file.
        output_path: Path to write the docked pose PDBQT.
        grid_box: Dictionary with center_x/y/z and size_x/y/z.

    Returns:
        Best binding affinity in kcal/mol.

    Raises:
        RuntimeError: If both Vina Python bindings and CLI fail.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    center = [grid_box["center_x"], grid_box["center_y"], grid_box["center_z"]]
    box_size = [grid_box["size_x"], grid_box["size_y"], grid_box["size_z"]]

    # --- Attempt 1: Vina Python bindings ---
    try:
        from vina import Vina

        v = Vina(sf_name="vina")
        v.set_receptor(str(receptor_path))
        v.set_ligand_from_file(str(ligand_path))
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=NUM_MODES)

        energies = v.energies(n_poses=1)
        v.write_poses(str(output_path), n_poses=1)
        return float(energies[0][0])
    except ImportError:
        logger.debug("Vina Python bindings not available, falling back to CLI.")
    except Exception as exc:
        logger.warning(
            "Vina Python bindings failed (%s). Trying CLI.", exc,
        )

    # --- Attempt 2: Vina CLI ---
    cmd = [
        VINA_BIN,
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--center_x", str(grid_box["center_x"]),
        "--center_y", str(grid_box["center_y"]),
        "--center_z", str(grid_box["center_z"]),
        "--size_x", str(grid_box["size_x"]),
        "--size_y", str(grid_box["size_y"]),
        "--size_z", str(grid_box["size_z"]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--out", str(output_path),
    ]

    logger.info(
        "Running Vina CLI: %s vs %s", receptor_path.name, ligand_path.name,
    )

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=VINA_TIMEOUT,
        )
        return _parse_vina_output(result.stdout)
    except FileNotFoundError:
        msg = (
            f"Vina not found at {VINA_BIN}. "
            "Install with: pip install vina"
        )
        raise RuntimeError(msg)
    except subprocess.CalledProcessError as exc:
        msg = (
            f"Vina CLI failed for {receptor_path.name} vs "
            f"{ligand_path.name}: {exc.stderr}"
        )
        raise RuntimeError(msg)
    except subprocess.TimeoutExpired:
        msg = (
            f"Vina CLI timed out for {receptor_path.name} vs "
            f"{ligand_path.name} ({VINA_TIMEOUT}s limit)."
        )
        raise RuntimeError(msg)


# ---------------------------------------------------------------------------
# Docking orchestration
# ---------------------------------------------------------------------------


def _dock_compound_against_target(
    compound_name: str,
    ligand_pdbqt: Path,
    receptor_pdbqt: Path,
    target_label: str,
    grid_box: dict[str, float],
) -> float:
    """Dock a single compound against a single target receptor.

    Args:
        compound_name: Name of the compound (for filenames).
        ligand_pdbqt: Path to the ligand PDBQT file.
        receptor_pdbqt: Path to the receptor PDBQT file.
        target_label: Short label for the target (e.g. "TryR", "HumanGR").
        grid_box: Docking grid box parameters.

    Returns:
        Best binding affinity in kcal/mol.
    """
    output_path = (
        OUTPUT_DIR / "poses" / f"{compound_name}_{target_label}_out.pdbqt"
    )
    affinity = _run_vina(
        receptor_path=receptor_pdbqt,
        ligand_path=ligand_pdbqt,
        output_path=output_path,
        grid_box=grid_box,
    )
    logger.info(
        "%s vs %s: %.2f kcal/mol", compound_name, target_label, affinity,
    )
    return affinity


def _run_all_dockings(
    ligand_pdbqts: dict[str, Path],
    human_gr_pdbqt: Path,
    human_gr_grid: dict[str, float],
) -> dict[str, dict[str, float]]:
    """Dock all compounds against both TryR and human GR.

    Args:
        ligand_pdbqts: Mapping of compound name to ligand PDBQT path.
        human_gr_pdbqt: Path to the human GR receptor PDBQT.
        human_gr_grid: Grid box for human GR (blind docking centroid).

    Returns:
        Nested dict: {compound_name: {"TryR": affinity, "HumanGR": affinity}}.
    """
    results: dict[str, dict[str, float]] = {}

    for compound_name, ligand_path in ligand_pdbqts.items():
        logger.info("--- Docking %s ---", compound_name)
        affinities: dict[str, float] = {}

        # Dock against parasite TryR
        try:
            affinities["TryR"] = _dock_compound_against_target(
                compound_name=compound_name,
                ligand_pdbqt=ligand_path,
                receptor_pdbqt=TRYR_RECEPTOR,
                target_label="TryR",
                grid_box=TRYR_GRID,
            )
        except RuntimeError as exc:
            logger.error("Failed docking %s vs TryR: %s", compound_name, exc)
            affinities["TryR"] = 0.0

        # Dock against human GR
        try:
            affinities["HumanGR"] = _dock_compound_against_target(
                compound_name=compound_name,
                ligand_pdbqt=ligand_path,
                receptor_pdbqt=human_gr_pdbqt,
                target_label="HumanGR",
                grid_box=human_gr_grid,
            )
        except RuntimeError as exc:
            logger.error("Failed docking %s vs HumanGR: %s", compound_name, exc)
            affinities["HumanGR"] = 0.0

        results[compound_name] = affinities

    return results


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _compute_selectivity_table(
    affinities: dict[str, dict[str, float]],
) -> list[dict[str, Any]]:
    """Build a selectivity table from docking affinities.

    Delta is computed as HumanGR - TryR. A positive delta means the
    compound binds more tightly to TryR (desirable for selectivity).

    Args:
        affinities: Nested dict of compound -> target -> affinity.

    Returns:
        List of row dicts with keys: compound, tryr, human_gr, delta,
        selective.
    """
    rows: list[dict[str, Any]] = []
    for compound, targets in affinities.items():
        tryr = targets.get("TryR", 0.0)
        human_gr = targets.get("HumanGR", 0.0)
        delta = round(human_gr - tryr, 2)
        selective = delta >= SELECTIVITY_THRESHOLD
        rows.append({
            "compound": compound,
            "tryr": round(tryr, 2),
            "human_gr": round(human_gr, 2),
            "delta": delta,
            "selective": selective,
        })
    return rows


def _write_json_report(
    table: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the selectivity report as JSON.

    Args:
        table: List of selectivity row dicts.
        output_path: Path to write the JSON file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mrl003_row = next((r for r in table if r["compound"] == "MRL-003"), None)
    mrl003_delta = mrl003_row["delta"] if mrl003_row else 0.0
    mrl003_selective = mrl003_row["selective"] if mrl003_row else False

    report = {
        "module": "12_selectivity",
        "generated_at": datetime.now().isoformat(),
        "selectivity_threshold_kcal": SELECTIVITY_THRESHOLD,
        "mrl003_selectivity_index": mrl003_delta,
        "mrl003_selective": mrl003_selective,
        "targets": {
            "parasite": {
                "name": "L. infantum trypanothione reductase (TryR)",
                "uniprot": "Q4Q457",
                "receptor_pdbqt": str(TRYR_RECEPTOR),
                "grid_box": TRYR_GRID,
            },
            "human": {
                "name": "Human glutathione reductase (GR)",
                "uniprot": "P00390",
                "receptor_pdbqt": str(HUMAN_GR_PDBQT),
                "grid_box": "blind docking (centroid, 40x40x40)",
            },
        },
        "results": table,
    }

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    logger.info("Wrote JSON report: %s", output_path)


def _write_markdown_report(
    table: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the selectivity report as Markdown.

    Args:
        table: List of selectivity row dicts.
        output_path: Path to write the Markdown file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mrl003_row = next((r for r in table if r["compound"] == "MRL-003"), None)
    mrl003_delta = mrl003_row["delta"] if mrl003_row else 0.0

    lines: list[str] = [
        "# MRL-003 Selectivity Analysis",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Targets",
        "",
        "- **Parasite**: *L. infantum* trypanothione reductase (TryR, UniProt Q4Q457)",
        "- **Human**: Glutathione reductase (GR, UniProt P00390)",
        "",
        "Human GR is the closest human homolog to TryR. A selective inhibitor",
        "should bind TryR significantly more tightly than human GR.",
        "",
        "## Binding Affinities",
        "",
        "| Compound | L. infantum TryR | Human GR | Delta | Selective? |",
        "|----------|:----------------:|:--------:|:-----:|:----------:|",
    ]

    for row in table:
        selective_str = "Yes" if row["selective"] else "No"
        lines.append(
            f"| {row['compound']} "
            f"| {row['tryr']:.2f} kcal/mol "
            f"| {row['human_gr']:.2f} kcal/mol "
            f"| {row['delta']:+.2f} "
            f"| {selective_str} |"
        )

    lines.extend([
        "",
        "## Interpretation",
        "",
        f"- Delta = HumanGR - TryR (positive = selective for parasite enzyme)",
        f"- Selectivity threshold: >= {SELECTIVITY_THRESHOLD:.1f} kcal/mol",
        f"- **MRL-003 selectivity index: {mrl003_delta:+.2f} kcal/mol**",
        "",
    ])

    if mrl003_delta >= SELECTIVITY_THRESHOLD:
        lines.append(
            "MRL-003 shows favorable selectivity for the parasite enzyme, "
            "binding TryR more tightly than human GR by "
            f"{mrl003_delta:.2f} kcal/mol."
        )
    else:
        lines.append(
            "MRL-003 does NOT meet the selectivity threshold. "
            "Further optimization may be needed to improve differential "
            "binding between TryR and human GR."
        )

    lines.extend([
        "",
        "## Reference Compounds",
        "",
        "- **Menadione**: Known human GR inhibitor. Expected to bind human GR",
        "  more tightly than TryR (negative delta).",
        "- **Pemetrexed**: Structural analog of MRL-003. Included for comparison.",
        "",
        "## Methods",
        "",
        "- Docking: AutoDock Vina (exhaustiveness=8)",
        "- TryR grid: Known binding site (center 11.0, 24.0, 15.0; size 25x25x25)",
        "- Human GR grid: Blind docking (centroid, 40x40x40 box)",
        "- Human GR structure: AlphaFold prediction (AF-P00390-F1-model_v6)",
        "",
    ])

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    logger.info("Wrote Markdown report: %s", output_path)


# ---------------------------------------------------------------------------
# PyMOL comparison script
# ---------------------------------------------------------------------------


def _write_pymol_script(output_path: Path) -> None:
    """Generate a PyMOL script to compare TryR and human GR binding poses.

    Loads both receptors and MRL-003 docked poses side by side with
    distinct coloring for easy visual comparison.

    Args:
        output_path: Path to write the .pml file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tryr_pdb = STRUCTURES_DIR / "TryR_Q4Q457.pdb"
    human_gr_pdb = HUMAN_GR_PDB
    mrl003_tryr_pose = OUTPUT_DIR / "poses" / "MRL-003_TryR_out.pdbqt"
    mrl003_gr_pose = OUTPUT_DIR / "poses" / "MRL-003_HumanGR_out.pdbqt"

    script = f"""\
# MRL-003 Selectivity Comparison -- PyMOL Script
# Generated by pipeline.12_selectivity
# Usage: pymol {output_path.name}

# --- Load structures ---
load {tryr_pdb}, TryR
load {human_gr_pdb}, HumanGR
load {mrl003_tryr_pose}, MRL003_TryR_pose
load {mrl003_gr_pose}, MRL003_GR_pose

# --- Style receptors ---
hide everything, all
show cartoon, TryR
show cartoon, HumanGR
color marine, TryR
color salmon, HumanGR
set cartoon_transparency, 0.5, TryR
set cartoon_transparency, 0.5, HumanGR

# --- Style ligands ---
show sticks, MRL003_TryR_pose
show sticks, MRL003_GR_pose
color yellow, MRL003_TryR_pose
color cyan, MRL003_GR_pose
set stick_radius, 0.2

# --- Arrange side by side ---
# Move HumanGR and its ligand to the right
translate [80, 0, 0], HumanGR
translate [80, 0, 0], MRL003_GR_pose

# --- Labels ---
pseudoatom label_tryr, pos=[0, 50, 0]
pseudoatom label_gr, pos=[80, 50, 0]
label label_tryr, "L. infantum TryR"
label label_gr, "Human GR"
set label_size, 20

# --- Binding site surface ---
select tryr_site, TryR within 5 of MRL003_TryR_pose
select gr_site, HumanGR within 5 of MRL003_GR_pose
show surface, tryr_site
show surface, gr_site
set surface_transparency, 0.7, tryr_site
set surface_transparency, 0.7, gr_site

# --- Camera ---
zoom all
set ray_shadows, 0
set antialias, 2
bg_color white

# --- Save session ---
# save {OUTPUT_DIR}/selectivity_comparison.pse
"""

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(script)

    logger.info("Wrote PyMOL script: %s", output_path)


# ---------------------------------------------------------------------------
# Dry run
# ---------------------------------------------------------------------------


def _dry_run() -> str:
    """Execute a dry run that skips downloads and docking.

    Generates reports with placeholder values and logs what would be done.

    Returns:
        Path to the output directory.
    """
    logger.info("[DRY RUN] Starting selectivity analysis (dry run mode).")

    # Log what would happen
    logger.info("[DRY RUN] Would download human GR structure from: %s", ALPHAFOLD_URL)
    logger.info("[DRY RUN] Would save PDB to: %s", HUMAN_GR_PDB)
    logger.info("[DRY RUN] Would convert PDB to PDBQT: %s", HUMAN_GR_PDBQT)

    for name, info in COMPOUNDS.items():
        logger.info(
            "[DRY RUN] Would generate ligand PDBQT for %s (SMILES: %s)",
            name,
            info["smiles"][:40] + "...",
        )
        logger.info("[DRY RUN] Would dock %s against TryR and human GR.", name)

    # Use placeholder affinities
    table = _compute_selectivity_table(DRY_RUN_AFFINITIES)

    # Write reports
    json_path = OUTPUT_DIR / "selectivity_report.json"
    md_path = OUTPUT_DIR / "selectivity_report.md"
    pml_path = OUTPUT_DIR / "compare_binding.pml"

    _write_json_report(table, json_path)
    _write_markdown_report(table, md_path)
    _write_pymol_script(pml_path)

    logger.info("[DRY RUN] Generated reports with placeholder affinities.")
    logger.info("[DRY RUN] Complete. Output directory: %s", OUTPUT_DIR)

    return str(OUTPUT_DIR)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_selectivity(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the full selectivity analysis pipeline.

    Steps:
        1. Download human GR structure from AlphaFold.
        2. Convert human GR PDB to PDBQT.
        3. Generate ligand PDBQTs for all compounds.
        4. Dock all compounds against parasite TryR.
        5. Dock all compounds against human GR.
        6. Compute selectivity metrics.
        7. Generate JSON/Markdown reports and PyMOL script.

    Args:
        force: If True, re-run even if outputs already exist.
        dry_run: If True, skip downloads/docking and use placeholder values.

    Returns:
        Path to the output directory (``results/selectivity``).
    """
    logger.info("=== Module 12: Selectivity Analysis ===")

    # Early exit for dry run
    if dry_run:
        return _dry_run()

    # Check for existing output
    json_path = OUTPUT_DIR / "selectivity_report.json"
    if json_path.exists() and not force:
        logger.info(
            "Selectivity report already exists at %s. Use --force to re-run.",
            json_path,
        )
        return str(OUTPUT_DIR)

    # ------------------------------------------------------------------
    # Step 1: Download human GR structure
    # ------------------------------------------------------------------
    if HUMAN_GR_PDB.exists() and not force:
        logger.info("Human GR PDB already exists: %s", HUMAN_GR_PDB)
    else:
        _download_file(ALPHAFOLD_URL, HUMAN_GR_PDB)

    # ------------------------------------------------------------------
    # Step 2: Convert human GR PDB to PDBQT
    # ------------------------------------------------------------------
    if HUMAN_GR_PDBQT.exists() and not force:
        logger.info("Human GR PDBQT already exists: %s", HUMAN_GR_PDBQT)
    else:
        _pdb_to_pdbqt(HUMAN_GR_PDB, HUMAN_GR_PDBQT)

    # Verify TryR receptor exists
    if not TRYR_RECEPTOR.exists():
        msg = f"TryR receptor PDBQT not found: {TRYR_RECEPTOR}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    # ------------------------------------------------------------------
    # Step 3: Compute human GR grid box (blind docking centroid)
    # ------------------------------------------------------------------
    logger.info("Computing blind docking grid box for human GR ...")
    human_gr_grid = _compute_centroid(HUMAN_GR_PDBQT)
    logger.info(
        "Human GR centroid: (%.1f, %.1f, %.1f), box size: %.0fx%.0fx%.0f",
        human_gr_grid["center_x"],
        human_gr_grid["center_y"],
        human_gr_grid["center_z"],
        human_gr_grid["size_x"],
        human_gr_grid["size_y"],
        human_gr_grid["size_z"],
    )

    # ------------------------------------------------------------------
    # Step 4: Generate ligand PDBQTs
    # ------------------------------------------------------------------
    ligand_dir = OUTPUT_DIR / "ligands"
    ligand_pdbqts: dict[str, Path] = {}

    for compound_name, compound_info in COMPOUNDS.items():
        pdbqt_path = ligand_dir / f"{compound_name}.pdbqt"
        if pdbqt_path.exists() and not force:
            logger.info("Ligand PDBQT already exists: %s", pdbqt_path)
            ligand_pdbqts[compound_name] = pdbqt_path
        else:
            ligand_pdbqts[compound_name] = _smiles_to_pdbqt(
                smiles=compound_info["smiles"],
                name=compound_name,
                output_dir=ligand_dir,
            )

    # ------------------------------------------------------------------
    # Step 5: Dock all compounds against both targets
    # ------------------------------------------------------------------
    logger.info("Starting docking campaign (%d compounds x 2 targets) ...", len(COMPOUNDS))

    affinities = _run_all_dockings(
        ligand_pdbqts=ligand_pdbqts,
        human_gr_pdbqt=HUMAN_GR_PDBQT,
        human_gr_grid=human_gr_grid,
    )

    # ------------------------------------------------------------------
    # Step 6: Compute selectivity and generate reports
    # ------------------------------------------------------------------
    table = _compute_selectivity_table(affinities)

    json_report_path = OUTPUT_DIR / "selectivity_report.json"
    md_report_path = OUTPUT_DIR / "selectivity_report.md"
    pml_path = OUTPUT_DIR / "compare_binding.pml"

    _write_json_report(table, json_report_path)
    _write_markdown_report(table, md_report_path)
    _write_pymol_script(pml_path)

    # ------------------------------------------------------------------
    # Step 7: Log summary
    # ------------------------------------------------------------------
    mrl003_row = next((r for r in table if r["compound"] == "MRL-003"), None)
    if mrl003_row:
        logger.info(
            "MRL-003 selectivity index: %+.2f kcal/mol (selective: %s)",
            mrl003_row["delta"],
            mrl003_row["selective"],
        )

    logger.info("Selectivity analysis complete. Output: %s", OUTPUT_DIR)
    return str(OUTPUT_DIR)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace with ``force`` and ``dry_run`` flags.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Marley -- Selectivity analysis: compare MRL-003 binding to "
            "parasite TryR vs human glutathione reductase (GR)."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run all steps even if outputs already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Skip downloads and docking; generate reports with placeholder "
            "values."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = validate_selectivity(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
