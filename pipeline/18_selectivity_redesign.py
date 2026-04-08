"""Redesign MRL-003 variants for selective TryR binding over human GR.

MRL-003 binds human GR (-8.68 kcal/mol) better than L. infantum TryR
(-7.32 kcal/mol).  This module generates 20 structural variants that
exploit charge, size, and spermidine-mimicry differences between the
two binding pockets, validates them with RDKit, docks against both
targets, and ranks by selectivity index (delta = GR - TryR).

Usage:
    python -m pipeline.18_selectivity_redesign
    python -m pipeline.18_selectivity_redesign --force
    python -m pipeline.18_selectivity_redesign --dry-run
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("selectivity_redesign")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Paths
STRUCTURES_DIR: Final[Path] = Path("data/structures")
OUTPUT_DIR: Final[Path] = Path("results/selectivity_redesign")
TRYR_RECEPTOR: Final[Path] = STRUCTURES_DIR / "TryR_Q4Q457_receptor.pdbqt"
HUMAN_GR_PDBQT: Final[Path] = STRUCTURES_DIR / "human_GR_P00390_receptor.pdbqt"

# Tool paths (project virtual environment)
OBABEL_BIN: Final[str] = ".venv/bin/obabel"
VINA_BIN: Final[str] = ".venv/bin/vina"

# Docking parameters
EXHAUSTIVENESS: Final[int] = 8
NUM_MODES: Final[int] = 5
VINA_TIMEOUT: Final[int] = 600
BLIND_BOX_SIZE: Final[float] = 40.0

# TryR known binding site
TRYR_GRID: Final[dict[str, float]] = {
    "center_x": 11.0,
    "center_y": 24.0,
    "center_z": 15.0,
    "size_x": 25.0,
    "size_y": 25.0,
    "size_z": 25.0,
}

# Selectivity threshold (kcal/mol)
SELECTIVITY_THRESHOLD: Final[float] = 1.5

# Lipinski violation tolerance for antiparasitic compounds
MAX_LIPINSKI_VIOLATIONS: Final[int] = 2

# Original MRL-003 reference values
MRL003_TRYR_AFFINITY: Final[float] = -7.32
MRL003_GR_AFFINITY: Final[float] = -8.68
MRL003_DELTA: Final[float] = MRL003_GR_AFFINITY - MRL003_TRYR_AFFINITY

# MRL-003 base SMILES
MRL003_SMILES: Final[str] = (
    "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)"
    "NC(CCC(=O)N)C(=O)N"
)

# ---------------------------------------------------------------------------
# Variant library
# ---------------------------------------------------------------------------

VARIANTS: Final[list[dict[str, str]]] = [
    # CATEGORY 1: Add positive charges (attract TryR negative pocket,
    #             repel GR positive pocket)
    {
        "name": "MRL-101",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCN",
        "modification": "charge+",
        "description": "Ethylamine tail",
    },
    {
        "name": "MRL-102",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCCN",
        "modification": "charge+",
        "description": "Propylamine tail",
    },
    {
        "name": "MRL-103",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)NCCN)C(=O)N",
        "modification": "charge+",
        "description": "Amine on glutamine chain",
    },
    {
        "name": "MRL-104",
        "smiles": "NC1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)N",
        "modification": "charge+",
        "description": "Amine on aromatic ring",
    },
    {
        "name": "MRL-105",
        "smiles": "C1=CC(=C(N)C=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)N",
        "modification": "charge+",
        "description": "Ortho-amine on ring",
    },
    # CATEGORY 2: Guanidinium groups (strongly cationic)
    {
        "name": "MRL-106",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NC(=N)N",
        "modification": "guanidinium",
        "description": "Guanidinium tail",
    },
    {
        "name": "MRL-107",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)NC(=N)N)C(=O)N",
        "modification": "guanidinium",
        "description": "Guanidinium on side chain",
    },
    {
        "name": "MRL-108",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCNC(=N)N",
        "modification": "guanidinium",
        "description": "Extended guanidinium tail",
    },
    # CATEGORY 3: Increased size (fill TryR wider pocket, clash with
    #             GR narrow pocket)
    {
        "name": "MRL-109",
        "smiles": "C1=CC(=CC=C1CCCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)N",
        "modification": "size",
        "description": "Extended linker (+1 CH2)",
    },
    {
        "name": "MRL-110",
        "smiles": "C1=CC(=CC=C1CCCCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)N",
        "modification": "size",
        "description": "Extended linker (+2 CH2)",
    },
    {
        "name": "MRL-111",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCCCC(=O)N)C(=O)N",
        "modification": "size",
        "description": "Extended glutamine chain (+2 CH2)",
    },
    {
        "name": "MRL-112",
        "smiles": "C1=CC(=CC=C1CC(CC)C2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)N",
        "modification": "size",
        "description": "Branched linker (ethyl)",
    },
    {
        "name": "MRL-113",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)NC1CCCCC1)C(=O)N",
        "modification": "size",
        "description": "Cyclohexyl on tail",
    },
    # CATEGORY 4: Spermidine-like tail (mimics trypanothione unique feature)
    {
        "name": "MRL-114",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCCNCCCCN",
        "modification": "spermidine",
        "description": "Spermidine tail (3+4)",
    },
    {
        "name": "MRL-115",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)NCCCNCCCCN)C(=O)N",
        "modification": "spermidine",
        "description": "Spermidine on side chain",
    },
    {
        "name": "MRL-116",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCNCCN",
        "modification": "spermidine",
        "description": "Short diamine tail (2+2)",
    },
    {
        "name": "MRL-117",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCCN",
        "modification": "spermidine",
        "description": "Aminopropyl tail",
    },
    # CATEGORY 5: Combined (charge + size or spermidine + charge)
    {
        "name": "MRL-118",
        "smiles": "NC1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCCN",
        "modification": "combined",
        "description": "Ring amine + aminopropyl tail",
    },
    {
        "name": "MRL-119",
        "smiles": "C1=CC(=C(N)C=C1CCCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)N)C(=O)NCCN",
        "modification": "combined",
        "description": "Ortho-amine + extended linker + ethylamine",
    },
    {
        "name": "MRL-120",
        "smiles": "C1=CC(=CC=C1CCC2=CNC3=C2C(=O)NC(=N3)N)C(=O)NC(CCC(=O)NC(=N)N)C(=O)NCCCNCCCCN",
        "modification": "combined",
        "description": "Guanidinium side + spermidine tail",
    },
]


# ---------------------------------------------------------------------------
# SMILES validation and Lipinski properties
# ---------------------------------------------------------------------------


def _validate_and_compute_properties(
    variants: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Validate SMILES and compute Lipinski descriptors for each variant.

    Uses RDKit to parse SMILES, compute molecular weight, LogP, hydrogen
    bond donors (HBD) and acceptors (HBA).  Lipinski violations are
    counted and flagged against the antiparasitic tolerance of 2.

    Args:
        variants: List of variant dicts with ``name`` and ``smiles`` keys.

    Returns:
        List of property dicts with keys: name, smiles, modification,
        description, valid, mw, logp, hbd, hba, lipinski_violations,
        lipinski_pass.
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    results: list[dict[str, Any]] = []

    for variant in variants:
        name = variant["name"]
        smiles = variant["smiles"]
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            logger.warning("Invalid SMILES for %s: %s", name, smiles)
            results.append({
                "name": name,
                "smiles": smiles,
                "modification": variant["modification"],
                "description": variant["description"],
                "valid": False,
                "mw": 0.0,
                "logp": 0.0,
                "hbd": 0,
                "hba": 0,
                "lipinski_violations": 99,
                "lipinski_pass": False,
            })
            continue

        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1

        lipinski_pass = violations <= MAX_LIPINSKI_VIOLATIONS

        logger.info(
            "%s: MW=%.1f, LogP=%.2f, HBD=%d, HBA=%d, violations=%d (%s)",
            name, mw, logp, hbd, hba, violations,
            "PASS" if lipinski_pass else "FAIL",
        )

        results.append({
            "name": name,
            "smiles": smiles,
            "modification": variant["modification"],
            "description": variant["description"],
            "valid": True,
            "mw": mw,
            "logp": logp,
            "hbd": hbd,
            "hba": hba,
            "lipinski_violations": violations,
            "lipinski_pass": lipinski_pass,
        })

    valid_count = sum(1 for r in results if r["valid"])
    logger.info(
        "Validated %d / %d variants (SMILES OK).", valid_count, len(results),
    )
    return results


# ---------------------------------------------------------------------------
# Structure preparation
# ---------------------------------------------------------------------------


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
        msg = (
            f"Open Babel not found at {OBABEL_BIN}. "
            "Install with: pip install openbabel-wheel"
        )
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


def _dock_variant(
    variant_name: str,
    ligand_pdbqt: Path,
    receptor_pdbqt: Path,
    target_label: str,
    grid_box: dict[str, float],
) -> float:
    """Dock a single variant against a single target receptor.

    Args:
        variant_name: Name of the variant (for filenames).
        ligand_pdbqt: Path to the ligand PDBQT file.
        receptor_pdbqt: Path to the receptor PDBQT file.
        target_label: Short label for the target (e.g. "TryR", "HumanGR").
        grid_box: Docking grid box parameters.

    Returns:
        Best binding affinity in kcal/mol.
    """
    output_path = (
        OUTPUT_DIR / "poses" / f"{variant_name}_{target_label}_out.pdbqt"
    )
    affinity = _run_vina(
        receptor_path=receptor_pdbqt,
        ligand_path=ligand_pdbqt,
        output_path=output_path,
        grid_box=grid_box,
    )
    logger.info(
        "%s vs %s: %.2f kcal/mol", variant_name, target_label, affinity,
    )
    return affinity


def _dock_all_variants(
    ligand_pdbqts: dict[str, Path],
    human_gr_grid: dict[str, float],
) -> dict[str, dict[str, float]]:
    """Dock all variants against both TryR and human GR.

    Args:
        ligand_pdbqts: Mapping of variant name to ligand PDBQT path.
        human_gr_grid: Grid box for human GR (blind docking centroid).

    Returns:
        Nested dict: {variant_name: {"TryR": affinity, "HumanGR": affinity}}.
    """
    results: dict[str, dict[str, float]] = {}

    for variant_name, ligand_path in ligand_pdbqts.items():
        logger.info("--- Docking %s ---", variant_name)
        affinities: dict[str, float] = {}

        # Dock against parasite TryR
        try:
            affinities["TryR"] = _dock_variant(
                variant_name=variant_name,
                ligand_pdbqt=ligand_path,
                receptor_pdbqt=TRYR_RECEPTOR,
                target_label="TryR",
                grid_box=TRYR_GRID,
            )
        except RuntimeError as exc:
            logger.error("Failed docking %s vs TryR: %s", variant_name, exc)
            affinities["TryR"] = 0.0

        # Dock against human GR
        try:
            affinities["HumanGR"] = _dock_variant(
                variant_name=variant_name,
                ligand_pdbqt=ligand_path,
                receptor_pdbqt=HUMAN_GR_PDBQT,
                target_label="HumanGR",
                grid_box=human_gr_grid,
            )
        except RuntimeError as exc:
            logger.error("Failed docking %s vs HumanGR: %s", variant_name, exc)
            affinities["HumanGR"] = 0.0

        results[variant_name] = affinities

    return results


# ---------------------------------------------------------------------------
# Selectivity computation
# ---------------------------------------------------------------------------


def _compute_selectivity(
    affinities: dict[str, dict[str, float]],
    properties: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Compute selectivity metrics for each variant.

    Delta = GR_affinity - TryR_affinity.  A positive delta means the
    compound binds TryR more tightly than human GR (desirable).

    Args:
        affinities: Nested dict of variant -> target -> affinity.
        properties: List of property dicts from validation step.

    Returns:
        List of row dicts with selectivity and Lipinski data.
    """
    prop_lookup = {p["name"]: p for p in properties}
    rows: list[dict[str, Any]] = []

    for variant_name, targets in affinities.items():
        tryr = targets.get("TryR", 0.0)
        human_gr = targets.get("HumanGR", 0.0)
        delta = round(human_gr - tryr, 2)
        selective = delta >= SELECTIVITY_THRESHOLD

        props = prop_lookup.get(variant_name, {})
        lipinski_pass = props.get("lipinski_pass", False)
        lipinski_violations = props.get("lipinski_violations", 99)

        rows.append({
            "name": variant_name,
            "modification": props.get("modification", ""),
            "description": props.get("description", ""),
            "tryr_affinity": round(tryr, 2),
            "gr_affinity": round(human_gr, 2),
            "delta": delta,
            "selective": selective,
            "lipinski_violations": lipinski_violations,
            "lipinski_pass": lipinski_pass,
        })

    # Sort by delta descending (best selectivity first)
    rows.sort(key=lambda r: r["delta"], reverse=True)
    return rows


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _write_variants_library_csv(
    properties: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the variant library with Lipinski properties to CSV.

    Args:
        properties: List of property dicts from validation step.
        output_path: Path to write the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "name", "smiles", "modification", "description",
        "mw", "logp", "hbd", "hba", "lipinski_violations", "lipinski_pass",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for prop in properties:
            if not prop["valid"]:
                continue
            writer.writerow({k: prop[k] for k in fieldnames})

    logger.info("Wrote variants library CSV: %s", output_path)


def _write_selectivity_matrix_csv(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the selectivity matrix to CSV.

    Args:
        rows: List of selectivity row dicts.
        output_path: Path to write the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "name", "tryr_affinity", "gr_affinity", "delta",
        "selective", "lipinski_pass",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})

    logger.info("Wrote selectivity matrix CSV: %s", output_path)


def _write_markdown_report(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the selectivity redesign report as Markdown.

    Includes full results table, top 5 selective variants, and a
    comparison with the original MRL-003.

    Args:
        rows: List of selectivity row dicts (sorted by delta descending).
        output_path: Path to write the Markdown file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Find the best variant
    best_row = rows[0] if rows else None
    selective_count = sum(1 for r in rows if r["selective"])
    top5 = [r for r in rows if r["selective"]][:5]

    lines: list[str] = [
        "# MRL-003 Selectivity Redesign",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Strategy",
        "",
        "Exploiting charge, size, and spermidine-mimicry differences between "
        "TryR and human GR.",
        "",
        "Key differences between the two binding pockets:",
        "",
        "- **TryR**: net negative charge in binding site (Glu18, Glu466, Asp327)",
        "- **Human GR**: net positive charge (Arg37, Arg347)",
        "- **TryR**: wider pocket (accepts trypanothione MW~723)",
        "- **Human GR**: narrower pocket (accepts glutathione MW~307)",
        "",
        "## Results",
        "",
        "| # | Variant | Modification | TryR (kcal/mol) | Human GR | Delta "
        "| Selective? | Lipinski |",
        "|---|---------|-------------|:----------------:|:--------:|:-----:"
        "|:----------:|:--------:|",
    ]

    for idx, row in enumerate(rows, 1):
        selective_str = "YES" if row["selective"] else "no"
        lipinski_str = (
            f"PASS ({row['lipinski_violations']})"
            if row["lipinski_pass"]
            else f"FAIL ({row['lipinski_violations']})"
        )
        lines.append(
            f"| {idx} "
            f"| {row['name']} "
            f"| {row['modification']} "
            f"| {row['tryr_affinity']:.2f} "
            f"| {row['gr_affinity']:.2f} "
            f"| {row['delta']:+.2f} "
            f"| {selective_str} "
            f"| {lipinski_str} |"
        )

    lines.extend([
        "",
        f"**{selective_count} / {len(rows)} variants meet the selectivity "
        f"threshold (delta >= {SELECTIVITY_THRESHOLD:.1f} kcal/mol).**",
        "",
    ])

    # Top 5 selective variants
    if top5:
        lines.extend([
            "## Best Selective Variants",
            "",
        ])
        for idx, row in enumerate(top5, 1):
            lines.append(
                f"{idx}. **{row['name']}** ({row['description']}): "
                f"delta = {row['delta']:+.2f} kcal/mol, "
                f"TryR = {row['tryr_affinity']:.2f}, "
                f"GR = {row['gr_affinity']:.2f}"
            )
        lines.append("")
    else:
        lines.extend([
            "## Best Selective Variants",
            "",
            "No variants met the selectivity threshold.",
            "",
        ])

    # Comparison with original MRL-003
    lines.extend([
        "## Comparison with Original MRL-003",
        "",
        "| Metric | MRL-003 | Best Variant |",
        "|--------|:-------:|:------------:|",
    ])

    if best_row:
        lines.extend([
            f"| TryR affinity | {MRL003_TRYR_AFFINITY:.2f} kcal/mol "
            f"| {best_row['tryr_affinity']:.2f} kcal/mol |",
            f"| Human GR affinity | {MRL003_GR_AFFINITY:.2f} kcal/mol "
            f"| {best_row['gr_affinity']:.2f} kcal/mol |",
            f"| Delta | {MRL003_DELTA:+.2f} "
            f"| {best_row['delta']:+.2f} |",
            f"| Selective? | NO "
            f"| {'YES' if best_row['selective'] else 'NO'} |",
        ])
    else:
        lines.extend([
            f"| TryR affinity | {MRL003_TRYR_AFFINITY:.2f} kcal/mol | N/A |",
            f"| Human GR affinity | {MRL003_GR_AFFINITY:.2f} kcal/mol | N/A |",
            f"| Delta | {MRL003_DELTA:+.2f} | N/A |",
            "| Selective? | NO | N/A |",
        ])

    lines.extend([
        "",
        "## Methods",
        "",
        "- Variant design: 5 categories (charge+, guanidinium, size, "
        "spermidine, combined)",
        "- SMILES validation: RDKit",
        "- Lipinski filter: up to 2 violations allowed (antiparasitic context)",
        "- 3D generation: Open Babel (--gen3d)",
        f"- Docking: AutoDock Vina (exhaustiveness={EXHAUSTIVENESS})",
        "- TryR grid: Known binding site (center 11.0, 24.0, 15.0; "
        "size 25x25x25)",
        "- Human GR grid: Blind docking (centroid, 40x40x40 box)",
        f"- Selectivity threshold: delta >= {SELECTIVITY_THRESHOLD:.1f} kcal/mol",
        "",
    ])

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    logger.info("Wrote Markdown report: %s", output_path)


def _write_json_report(
    rows: list[dict[str, Any]],
    properties: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the full redesign report as JSON.

    Args:
        rows: List of selectivity row dicts.
        properties: List of property dicts from validation step.
        output_path: Path to write the JSON file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    best_row = rows[0] if rows else None
    selective_variants = [r for r in rows if r["selective"]]

    report = {
        "module": "18_selectivity_redesign",
        "generated_at": datetime.now().isoformat(),
        "selectivity_threshold_kcal": SELECTIVITY_THRESHOLD,
        "total_variants": len(rows),
        "selective_count": len(selective_variants),
        "original_mrl003": {
            "tryr_affinity": MRL003_TRYR_AFFINITY,
            "gr_affinity": MRL003_GR_AFFINITY,
            "delta": MRL003_DELTA,
            "selective": False,
        },
        "best_variant": best_row if best_row else {},
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
        "variant_properties": [
            p for p in properties if p["valid"]
        ],
        "selectivity_results": rows,
    }

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    logger.info("Wrote JSON report: %s", output_path)


# ---------------------------------------------------------------------------
# Dry run
# ---------------------------------------------------------------------------


def _generate_dry_run_affinities(
    properties: list[dict[str, Any]],
) -> dict[str, dict[str, float]]:
    """Generate plausible fake docking affinities for dry-run mode.

    Charged and spermidine variants get better TryR and worse GR scores
    to simulate the intended selectivity shift.  Size variants get modest
    improvements.

    Args:
        properties: List of property dicts from validation step.

    Returns:
        Nested dict: {variant_name: {"TryR": affinity, "HumanGR": affinity}}.
    """
    random.seed(42)  # Reproducible fake data
    affinities: dict[str, dict[str, float]] = {}

    for prop in properties:
        if not prop["valid"]:
            continue

        name = prop["name"]
        mod = prop["modification"]

        # Base: similar to MRL-003
        base_tryr = -7.32
        base_gr = -8.68

        if mod == "charge+":
            tryr_shift = random.uniform(-1.2, -0.3)
            gr_shift = random.uniform(0.5, 1.8)
        elif mod == "guanidinium":
            tryr_shift = random.uniform(-1.5, -0.5)
            gr_shift = random.uniform(0.8, 2.2)
        elif mod == "size":
            tryr_shift = random.uniform(-0.8, 0.2)
            gr_shift = random.uniform(0.3, 1.5)
        elif mod == "spermidine":
            tryr_shift = random.uniform(-1.8, -0.5)
            gr_shift = random.uniform(0.6, 2.0)
        elif mod == "combined":
            tryr_shift = random.uniform(-1.5, -0.3)
            gr_shift = random.uniform(1.0, 2.5)
        else:
            tryr_shift = random.uniform(-0.5, 0.5)
            gr_shift = random.uniform(-0.5, 0.5)

        affinities[name] = {
            "TryR": round(base_tryr + tryr_shift, 2),
            "HumanGR": round(base_gr + gr_shift, 2),
        }

    return affinities


def _dry_run_pipeline(properties: list[dict[str, Any]]) -> str:
    """Execute a dry run that skips docking and uses fake affinities.

    Validates SMILES and computes Lipinski properties (real), then
    generates reports with simulated docking scores.

    Args:
        properties: List of property dicts from validation step.

    Returns:
        Path to the output directory.
    """
    logger.info("[DRY RUN] Skipping 3D generation and docking.")

    for prop in properties:
        if prop["valid"]:
            logger.info(
                "[DRY RUN] Would generate PDBQT for %s", prop["name"],
            )
            logger.info(
                "[DRY RUN] Would dock %s against TryR and human GR.",
                prop["name"],
            )

    # Generate fake affinities
    affinities = _generate_dry_run_affinities(properties)
    logger.info(
        "[DRY RUN] Generated simulated affinities for %d variants.",
        len(affinities),
    )

    # Compute selectivity
    rows = _compute_selectivity(affinities, properties)

    # Write reports
    _write_variants_library_csv(
        properties, OUTPUT_DIR / "variants_library.csv",
    )
    _write_selectivity_matrix_csv(
        rows, OUTPUT_DIR / "selectivity_matrix.csv",
    )
    _write_markdown_report(
        rows, OUTPUT_DIR / "selectivity_redesign_report.md",
    )
    _write_json_report(
        rows, properties, OUTPUT_DIR / "selectivity_redesign_report.json",
    )

    selective_count = sum(1 for r in rows if r["selective"])
    logger.info(
        "[DRY RUN] %d / %d variants selective (simulated data).",
        selective_count, len(rows),
    )
    logger.info("[DRY RUN] Complete. Output directory: %s", OUTPUT_DIR)

    return str(OUTPUT_DIR)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def redesign_for_selectivity(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the selectivity redesign pipeline.

    Steps:
        1. Define 20 MRL-003 structural variants.
        2. Validate SMILES and compute Lipinski properties (RDKit).
        3. Generate 3D PDBQT ligand files (Open Babel).
        4. Dock each variant against TryR and human GR (Vina).
        5. Compute selectivity index (delta = GR - TryR).
        6. Generate CSV, Markdown, and JSON reports.

    Args:
        force: If True, re-run even if outputs already exist.
        dry_run: If True, skip docking and use simulated affinities.

    Returns:
        Path to the output directory (``results/selectivity_redesign``).
    """
    logger.info("=== Module 18: Selectivity Redesign ===")
    logger.info("Evaluating %d MRL-003 variants.", len(VARIANTS))

    # Check for existing output
    json_path = OUTPUT_DIR / "selectivity_redesign_report.json"
    if json_path.exists() and not force and not dry_run:
        logger.info(
            "Redesign report already exists at %s. Use --force to re-run.",
            json_path,
        )
        return str(OUTPUT_DIR)

    # ------------------------------------------------------------------
    # Step 1 & 2: Validate SMILES + compute Lipinski properties
    # ------------------------------------------------------------------
    logger.info("Validating SMILES and computing Lipinski properties ...")
    properties = _validate_and_compute_properties(VARIANTS)

    valid_variants = [p for p in properties if p["valid"]]
    if not valid_variants:
        logger.error("No valid SMILES found. Cannot proceed.")
        return str(OUTPUT_DIR)

    logger.info(
        "%d / %d variants have valid SMILES.", len(valid_variants), len(VARIANTS),
    )

    # Early exit for dry run
    if dry_run:
        return _dry_run_pipeline(properties)

    # ------------------------------------------------------------------
    # Step 3: Verify receptor files exist
    # ------------------------------------------------------------------
    if not TRYR_RECEPTOR.exists():
        msg = f"TryR receptor PDBQT not found: {TRYR_RECEPTOR}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    if not HUMAN_GR_PDBQT.exists():
        msg = f"Human GR receptor PDBQT not found: {HUMAN_GR_PDBQT}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    # ------------------------------------------------------------------
    # Step 4: Compute human GR grid box (blind docking centroid)
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
    # Step 5: Generate ligand PDBQTs
    # ------------------------------------------------------------------
    ligand_dir = OUTPUT_DIR / "ligands"
    ligand_pdbqts: dict[str, Path] = {}

    for prop in valid_variants:
        name = prop["name"]
        pdbqt_path = ligand_dir / f"{name}.pdbqt"

        if pdbqt_path.exists() and not force:
            logger.info("Ligand PDBQT already exists: %s", pdbqt_path)
            ligand_pdbqts[name] = pdbqt_path
        else:
            try:
                ligand_pdbqts[name] = _smiles_to_pdbqt(
                    smiles=prop["smiles"],
                    name=name,
                    output_dir=ligand_dir,
                )
            except RuntimeError as exc:
                logger.error(
                    "Failed to generate PDBQT for %s: %s", name, exc,
                )

    if not ligand_pdbqts:
        logger.error("No ligand PDBQTs generated. Cannot proceed with docking.")
        return str(OUTPUT_DIR)

    logger.info("Generated %d ligand PDBQTs.", len(ligand_pdbqts))

    # ------------------------------------------------------------------
    # Step 6: Dock all variants against both targets
    # ------------------------------------------------------------------
    logger.info(
        "Starting docking campaign (%d variants x 2 targets) ...",
        len(ligand_pdbqts),
    )
    affinities = _dock_all_variants(
        ligand_pdbqts=ligand_pdbqts,
        human_gr_grid=human_gr_grid,
    )

    # ------------------------------------------------------------------
    # Step 7: Compute selectivity and generate reports
    # ------------------------------------------------------------------
    rows = _compute_selectivity(affinities, properties)

    _write_variants_library_csv(
        properties, OUTPUT_DIR / "variants_library.csv",
    )
    _write_selectivity_matrix_csv(
        rows, OUTPUT_DIR / "selectivity_matrix.csv",
    )
    _write_markdown_report(
        rows, OUTPUT_DIR / "selectivity_redesign_report.md",
    )
    _write_json_report(
        rows, properties, OUTPUT_DIR / "selectivity_redesign_report.json",
    )

    # ------------------------------------------------------------------
    # Step 8: Log summary
    # ------------------------------------------------------------------
    selective_count = sum(1 for r in rows if r["selective"])
    logger.info(
        "%d / %d variants meet selectivity threshold (delta >= %.1f).",
        selective_count, len(rows), SELECTIVITY_THRESHOLD,
    )

    if rows:
        best = rows[0]
        logger.info(
            "Best variant: %s (delta = %+.2f, TryR = %.2f, GR = %.2f)",
            best["name"], best["delta"],
            best["tryr_affinity"], best["gr_affinity"],
        )
        logger.info(
            "Original MRL-003: delta = %+.2f (TryR = %.2f, GR = %.2f)",
            MRL003_DELTA, MRL003_TRYR_AFFINITY, MRL003_GR_AFFINITY,
        )

    logger.info("Selectivity redesign complete. Output: %s", OUTPUT_DIR)
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
            "Marley -- Selectivity redesign: generate MRL-003 variants "
            "optimized for TryR selectivity over human GR."
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
            "Skip 3D generation and docking; generate reports with "
            "simulated affinities."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = redesign_for_selectivity(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
