"""Predict resistance risk for MRL-003 via alanine scanning mutagenesis.

Identifies binding site residues from the TryR-MRL-003 docking complex,
performs in-silico alanine scanning by re-docking against each mutant,
classifies mutations by impact on binding affinity, and cross-references
with catalytically essential residues to assess whether the parasite
could realistically evolve resistance.

Usage:
    python -m pipeline.13_resistance
    python -m pipeline.13_resistance --force
    python -m pipeline.13_resistance --dry-run
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("resistance")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent

DOCKING_OUTPUT_PATH: Final[Path] = Path(
    "data/docking/TryR/MRL-003_amide_tail_out.pdbqt"
)
TRYR_PDB_PATH: Final[Path] = Path("data/structures/TryR_Q4Q457.pdb")
TRYR_RECEPTOR_PDBQT_PATH: Final[Path] = Path(
    "data/structures/TryR_Q4Q457_receptor.pdbqt"
)
MRL003_LIGAND_PATH: Final[Path] = Path(
    "data/compounds/TryR/MRL-003_amide_tail.pdbqt"
)
OUTPUT_DIR: Final[Path] = Path("results/resistance")

OBABEL_BIN: Final[str] = ".venv/bin/obabel"
VINA_BIN: Final[str] = ".venv/bin/vina"

# Distance threshold (angstroms) for defining binding site residues.
CONTACT_DISTANCE: Final[float] = 5.0

# Docking grid box (same as original TryR docking campaign).
GRID_BOX: Final[dict[str, float]] = {
    "center_x": 11.0,
    "center_y": 24.0,
    "center_z": 15.0,
    "size_x": 25.0,
    "size_y": 25.0,
    "size_z": 25.0,
}
EXHAUSTIVENESS: Final[int] = 8

# Alanine scanning classification thresholds (kcal/mol).
CRITICAL_THRESHOLD: Final[float] = 2.0
SIGNIFICANT_THRESHOLD: Final[float] = 1.0
BENEFICIAL_THRESHOLD: Final[float] = -1.0

# Backbone atom names to keep when mutating to Ala (plus CB).
BACKBONE_ATOMS: Final[set[str]] = {"N", "CA", "C", "O", "CB", "OXT", "H", "HA"}

# Known catalytically essential residues of TryR across Leishmania species.
# These residues are strictly conserved and indispensable for enzymatic
# activity, meaning the parasite cannot mutate them without losing fitness.
CATALYTIC_RESIDUES: Final[dict[int, str]] = {
    52: "CYS",
    57: "CYS",
    461: "HIS",
    466: "GLU",
}

# Three-letter to one-letter amino acid code mapping.
AA3_TO_1: Final[dict[str, str]] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


class BindingSiteResidue:
    """A residue identified as part of the MRL-003 binding site."""

    __slots__ = ("resname", "resnum", "chain", "min_distance")

    def __init__(
        self,
        resname: str,
        resnum: int,
        chain: str,
        min_distance: float,
    ) -> None:
        self.resname = resname
        self.resnum = resnum
        self.chain = chain
        self.min_distance = min_distance

    def __repr__(self) -> str:
        return (
            f"{self.resname}{self.resnum}:{self.chain} "
            f"(d={self.min_distance:.2f} A)"
        )


class ScanResult:
    """Result of alanine scanning for a single residue."""

    __slots__ = (
        "resname", "resnum", "chain", "wt_affinity", "mut_affinity",
        "delta", "classification", "is_conserved", "is_catalytic",
    )

    def __init__(
        self,
        resname: str,
        resnum: int,
        chain: str,
        wt_affinity: float,
        mut_affinity: float,
    ) -> None:
        self.resname = resname
        self.resnum = resnum
        self.chain = chain
        self.wt_affinity = wt_affinity
        self.mut_affinity = mut_affinity
        self.delta = mut_affinity - wt_affinity
        self.classification = _classify_delta(self.delta)
        self.is_catalytic = resnum in CATALYTIC_RESIDUES
        self.is_conserved = self.is_catalytic


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _classify_delta(delta: float) -> str:
    """Classify the effect of an alanine mutation on binding affinity.

    Args:
        delta: Mutant affinity minus wildtype affinity (kcal/mol).
              Positive values mean weaker binding (worse for drug).

    Returns:
        One of CRITICAL, SIGNIFICANT, NEUTRAL, or BENEFICIAL.
    """
    if delta > CRITICAL_THRESHOLD:
        return "CRITICAL"
    if delta > SIGNIFICANT_THRESHOLD:
        return "SIGNIFICANT"
    if delta < BENEFICIAL_THRESHOLD:
        return "BENEFICIAL"
    return "NEUTRAL"


def _euclidean_distance(
    x1: float, y1: float, z1: float,
    x2: float, y2: float, z2: float,
) -> float:
    """Compute the Euclidean distance between two 3D points."""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


# ---------------------------------------------------------------------------
# Step 1: Parse structures and identify binding site residues
# ---------------------------------------------------------------------------


def _parse_ligand_coords(pdbqt_path: Path) -> list[tuple[float, float, float]]:
    """Extract all atom coordinates from the first MODEL of a ligand PDBQT.

    Only reads the first docking pose (MODEL 1).

    Args:
        pdbqt_path: Path to the docked ligand PDBQT output file.

    Returns:
        List of (x, y, z) tuples for every ATOM/HETATM record in MODEL 1.

    Raises:
        FileNotFoundError: If the PDBQT file does not exist.
        ValueError: If no atom coordinates are found.
    """
    if not pdbqt_path.exists():
        raise FileNotFoundError(f"Ligand PDBQT not found: {pdbqt_path}")

    coords: list[tuple[float, float, float]] = []
    in_model_1 = False

    with open(pdbqt_path, encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()

            if stripped.startswith("MODEL"):
                model_num = stripped.split()[-1]
                if model_num == "1":
                    in_model_1 = True
                else:
                    # Past model 1, stop reading.
                    break

            if stripped.startswith("ENDMDL") and in_model_1:
                break

            if in_model_1 and stripped.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue

    if not coords:
        raise ValueError(f"No atom coordinates found in {pdbqt_path}")

    logger.info(
        "Parsed %d ligand atom coordinates from %s.",
        len(coords),
        pdbqt_path.name,
    )
    return coords


def _parse_pdb_atoms(
    pdb_path: Path,
) -> list[dict[str, object]]:
    """Parse ATOM records from a PDB file.

    Args:
        pdb_path: Path to the PDB file.

    Returns:
        List of dicts with keys: atom_name, resname, chain, resnum,
        x, y, z, line (the original PDB line).

    Raises:
        FileNotFoundError: If the PDB file does not exist.
    """
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    atoms: list[dict[str, object]] = []

    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            try:
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip() or "A"
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue

            atoms.append({
                "atom_name": atom_name,
                "resname": resname,
                "chain": chain,
                "resnum": resnum,
                "x": x,
                "y": y,
                "z": z,
                "line": line,
            })

    logger.info("Parsed %d ATOM records from %s.", len(atoms), pdb_path.name)
    return atoms


def _find_binding_residues(
    ligand_coords: list[tuple[float, float, float]],
    pdb_atoms: list[dict[str, object]],
    distance_cutoff: float = CONTACT_DISTANCE,
) -> list[BindingSiteResidue]:
    """Identify protein residues within a distance cutoff of any ligand atom.

    Args:
        ligand_coords: Coordinates of the docked ligand atoms.
        pdb_atoms: Parsed ATOM records from the protein PDB.
        distance_cutoff: Maximum distance in angstroms to consider a contact.

    Returns:
        Sorted list of unique BindingSiteResidue objects.
    """
    # Track minimum distance per (chain, resnum).
    residue_distances: dict[tuple[str, int], tuple[str, float]] = {}

    for atom in pdb_atoms:
        ax = float(atom["x"])
        ay = float(atom["y"])
        az = float(atom["z"])
        chain = str(atom["chain"])
        resnum = int(atom["resnum"])
        resname = str(atom["resname"])
        key = (chain, resnum)

        for lx, ly, lz in ligand_coords:
            d = _euclidean_distance(ax, ay, az, lx, ly, lz)
            if d <= distance_cutoff:
                if key not in residue_distances or d < residue_distances[key][1]:
                    residue_distances[key] = (resname, d)
                # Found at least one contact atom; no need to check further
                # ligand atoms for a closer distance from this protein atom,
                # but we continue to find the global minimum for this residue.

    residues = [
        BindingSiteResidue(
            resname=resname,
            resnum=resnum,
            chain=chain,
            min_distance=min_dist,
        )
        for (chain, resnum), (resname, min_dist) in residue_distances.items()
    ]
    residues.sort(key=lambda r: (r.chain, r.resnum))

    logger.info(
        "Found %d binding site residues within %.1f A of MRL-003.",
        len(residues),
        distance_cutoff,
    )
    for res in residues:
        logger.info("  %s", res)

    return residues


# ---------------------------------------------------------------------------
# Step 2: Alanine scanning mutagenesis
# ---------------------------------------------------------------------------


def _parse_wildtype_affinity(pdbqt_path: Path) -> float:
    """Extract the best binding affinity from a docked PDBQT output.

    Reads the VINA RESULT remark from MODEL 1.

    Args:
        pdbqt_path: Path to a docked ligand PDBQT file.

    Returns:
        Binding affinity in kcal/mol (negative = favorable).

    Raises:
        ValueError: If no VINA RESULT line is found.
    """
    with open(pdbqt_path, encoding="utf-8") as fh:
        for line in fh:
            if "VINA RESULT" in line:
                parts = line.split()
                # Format: REMARK VINA RESULT:    -7.737      0.000      0.000
                for i, token in enumerate(parts):
                    if token == "RESULT:":
                        return float(parts[i + 1])
                # Fallback: try parsing the 4th whitespace-delimited token.
                return float(parts[3])

    raise ValueError(f"No VINA RESULT found in {pdbqt_path}")


def _mutate_to_alanine(
    pdb_path: Path,
    target_chain: str,
    target_resnum: int,
    output_path: Path,
) -> None:
    """Create a mutant PDB by replacing a residue with alanine.

    For the target residue:
    - Changes the residue name to ALA
    - Removes all side chain atoms beyond CB (keeps N, CA, C, O, CB, OXT, H, HA)

    Glycine residues are handled by simply renaming to ALA (CB is absent but
    this is acceptable for an approximate scanning approach).

    Args:
        pdb_path: Path to the wildtype PDB file.
        target_chain: Chain identifier of the residue to mutate.
        target_resnum: Residue number of the residue to mutate.
        output_path: Path to write the mutant PDB.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mutant_lines: list[str] = []

    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                mutant_lines.append(line)
                continue

            try:
                chain = line[21].strip() or "A"
                resnum = int(line[22:26])
                atom_name = line[12:16].strip()
            except (ValueError, IndexError):
                mutant_lines.append(line)
                continue

            if chain == target_chain and resnum == target_resnum:
                # Skip side chain atoms that are not part of alanine.
                if atom_name not in BACKBONE_ATOMS:
                    continue
                # Rename residue to ALA.
                line = line[:17] + "ALA" + line[20:]

            mutant_lines.append(line)

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.writelines(mutant_lines)

    logger.debug(
        "Wrote mutant PDB (chain %s, residue %d -> ALA) to %s.",
        target_chain,
        target_resnum,
        output_path.name,
    )


def _pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> None:
    """Convert a PDB file to PDBQT format using Open Babel.

    Args:
        pdb_path: Input PDB file path.
        pdbqt_path: Output PDBQT file path.

    Raises:
        RuntimeError: If obabel fails or is not found.
    """
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        OBABEL_BIN,
        str(pdb_path),
        "-O", str(pdbqt_path),
        "-xr",  # receptor mode (no torsion tree)
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=120,
        )
        logger.debug("obabel: %s", result.stderr.strip())
    except FileNotFoundError:
        raise RuntimeError(
            f"Open Babel not found at {OBABEL_BIN}. "
            "Ensure it is installed in the virtual environment."
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"obabel conversion failed for {pdb_path.name}: {exc.stderr}"
        )


def _run_vina_docking(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    output_pdbqt: Path,
) -> float:
    """Run AutoDock Vina and return the best binding affinity.

    Uses the standard TryR grid box parameters and reduced exhaustiveness
    for the scanning campaign.

    Args:
        receptor_pdbqt: Path to the receptor PDBQT file.
        ligand_pdbqt: Path to the ligand PDBQT file.
        output_pdbqt: Path to write the docked pose output.

    Returns:
        Best binding affinity in kcal/mol.

    Raises:
        RuntimeError: If vina fails or is not found.
    """
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        VINA_BIN,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(GRID_BOX["center_x"]),
        "--center_y", str(GRID_BOX["center_y"]),
        "--center_z", str(GRID_BOX["center_z"]),
        "--size_x", str(GRID_BOX["size_x"]),
        "--size_y", str(GRID_BOX["size_y"]),
        "--size_z", str(GRID_BOX["size_z"]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", "1",
        "--out", str(output_pdbqt),
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"Vina not found at {VINA_BIN}. "
            "Ensure it is installed in the virtual environment."
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Vina docking failed for {receptor_pdbqt.name}: {exc.stderr}"
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"Vina docking timed out for {receptor_pdbqt.name} (600s limit)."
        )

    # Parse best affinity from Vina output.
    pattern = re.compile(
        r"^\s*1\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)\s+([-+]?\d+\.?\d*)",
    )
    for line in result.stdout.splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1))

    # Fallback: parse the output PDBQT for VINA RESULT.
    return _parse_wildtype_affinity(output_pdbqt)


def _run_alanine_scan(
    residues: list[BindingSiteResidue],
    wt_affinity: float,
    pdb_path: Path,
    ligand_path: Path,
    work_dir: Path,
) -> list[ScanResult]:
    """Perform alanine scanning for each binding site residue.

    For each residue:
    1. Mutate the residue to alanine in the PDB.
    2. Convert the mutant PDB to PDBQT via obabel.
    3. Re-dock MRL-003 against the mutant receptor.
    4. Record the change in binding affinity.

    Args:
        residues: Binding site residues to scan.
        wt_affinity: Wildtype binding affinity (kcal/mol).
        pdb_path: Path to the wildtype TryR PDB.
        ligand_path: Path to the MRL-003 ligand PDBQT.
        work_dir: Directory for intermediate mutant files.

    Returns:
        List of ScanResult objects, one per residue.
    """
    results: list[ScanResult] = []

    for i, res in enumerate(residues, start=1):
        label = f"{res.resname}{res.resnum}:{res.chain}"
        logger.info(
            "[%d/%d] Scanning %s -> ALA ...",
            i,
            len(residues),
            label,
        )

        # 1. Mutate
        mutant_pdb = work_dir / f"mutant_{res.chain}_{res.resnum}_ALA.pdb"
        _mutate_to_alanine(pdb_path, res.chain, res.resnum, mutant_pdb)

        # 2. Convert to PDBQT
        mutant_pdbqt = work_dir / f"mutant_{res.chain}_{res.resnum}_ALA.pdbqt"
        try:
            _pdb_to_pdbqt(mutant_pdb, mutant_pdbqt)
        except RuntimeError as exc:
            logger.error("  Failed to convert mutant PDB: %s", exc)
            continue

        # 3. Re-dock
        docked_output = work_dir / f"docked_{res.chain}_{res.resnum}_ALA.pdbqt"
        try:
            mut_affinity = _run_vina_docking(
                mutant_pdbqt,
                ligand_path,
                docked_output,
            )
        except RuntimeError as exc:
            logger.error("  Docking failed for %s: %s", label, exc)
            continue

        # 4. Record
        scan = ScanResult(
            resname=res.resname,
            resnum=res.resnum,
            chain=res.chain,
            wt_affinity=wt_affinity,
            mut_affinity=mut_affinity,
        )
        results.append(scan)

        logger.info(
            "  %s: WT=%.3f, MUT=%.3f, delta=%+.3f kcal/mol -> %s%s",
            label,
            scan.wt_affinity,
            scan.mut_affinity,
            scan.delta,
            scan.classification,
            " (catalytic)" if scan.is_catalytic else "",
        )

    return results


# ---------------------------------------------------------------------------
# Step 3-4: Classification and conservation (integrated in ScanResult)
# ---------------------------------------------------------------------------


def _assess_resistance_risk(
    results: list[ScanResult],
) -> tuple[str, dict[str, int]]:
    """Determine overall resistance risk from alanine scanning results.

    Logic:
    - Count CRITICAL residues.
    - Of those, count how many are catalytically essential (cannot mutate).
    - Peripheral critical residues (could mutate) drive resistance risk.

    Args:
        results: Completed alanine scan results.

    Returns:
        Tuple of (risk_level, stats_dict).
        risk_level is one of LOW, MEDIUM, or HIGH.
        stats_dict contains counts for reporting.
    """
    critical = [r for r in results if r.classification == "CRITICAL"]
    significant = [r for r in results if r.classification == "SIGNIFICANT"]
    catalytic_critical = [r for r in critical if r.is_catalytic]
    peripheral_critical = [r for r in critical if not r.is_catalytic]

    stats = {
        "total_residues": len(results),
        "critical": len(critical),
        "significant": len(significant),
        "catalytic_critical": len(catalytic_critical),
        "peripheral_critical": len(peripheral_critical),
    }

    # Risk assessment logic:
    # - If all critical residues are catalytic, risk is LOW.
    # - If some peripheral critical residues exist, risk scales with count.
    if len(peripheral_critical) == 0:
        risk = "LOW"
    elif len(peripheral_critical) <= 2:
        risk = "MEDIUM"
    else:
        risk = "HIGH"

    logger.info(
        "Resistance risk: %s (critical=%d, catalytic=%d, peripheral=%d).",
        risk,
        stats["critical"],
        stats["catalytic_critical"],
        stats["peripheral_critical"],
    )
    return risk, stats


# ---------------------------------------------------------------------------
# Step 5: Generate outputs
# ---------------------------------------------------------------------------


def _write_scan_csv(results: list[ScanResult], output_path: Path) -> None:
    """Write alanine scanning results to CSV.

    Args:
        results: Completed scan results.
        output_path: Path for the output CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "residue",
        "chain",
        "resnum",
        "wildtype_affinity",
        "mutant_affinity",
        "delta",
        "classification",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in sorted(results, key=lambda x: x.delta, reverse=True):
            writer.writerow({
                "residue": r.resname,
                "chain": r.chain,
                "resnum": r.resnum,
                "wildtype_affinity": round(r.wt_affinity, 3),
                "mutant_affinity": round(r.mut_affinity, 3),
                "delta": round(r.delta, 3),
                "classification": r.classification,
            })

    logger.info("Wrote scan CSV to %s.", output_path)


def _write_resistance_report(
    results: list[ScanResult],
    residues: list[BindingSiteResidue],
    risk: str,
    stats: dict[str, int],
    wt_affinity: float,
    output_path: Path,
) -> None:
    """Generate and write the Markdown resistance report.

    Args:
        results: Completed scan results (may be empty in dry-run).
        residues: Binding site residues identified.
        risk: Overall resistance risk level.
        stats: Summary statistics dict.
        wt_affinity: Wildtype binding affinity.
        output_path: Path for the output Markdown file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines: list[str] = []
    lines.append("# MRL-003 Resistance Prediction")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")
    lines.append(f"Wildtype binding affinity: **{wt_affinity:.3f} kcal/mol**")
    lines.append(f"Binding site residues identified: **{len(residues)}**")
    lines.append(f"Contact distance cutoff: **{CONTACT_DISTANCE:.1f} A**")
    lines.append("")

    # -- Alanine Scanning Results --
    lines.append("## Alanine Scanning Results")
    lines.append("")

    if results:
        lines.append(
            "| Residue | Position | WT Affinity | Mutant Affinity "
            "| Delta | Impact | Conserved? |"
        )
        lines.append(
            "|---------|----------|:-----------:|:---------------:"
            "|:-----:|:------:|:----------:|"
        )

        sorted_results = sorted(results, key=lambda r: r.delta, reverse=True)
        for r in sorted_results:
            aa1 = AA3_TO_1.get(r.resname, "?")
            conserved = "Yes (catalytic)" if r.is_catalytic else "No"
            lines.append(
                f"| {r.resname} ({aa1}) | {r.chain}:{r.resnum} "
                f"| {r.wt_affinity:.3f} | {r.mut_affinity:.3f} "
                f"| {r.delta:+.3f} | {r.classification} | {conserved} |"
            )
    else:
        lines.append("*Dry run -- alanine scanning was not performed.*")
        lines.append("")
        lines.append("Binding site residues that would be scanned:")
        lines.append("")
        lines.append("| Residue | Position | Min Distance (A) |")
        lines.append("|---------|----------|:-----------------:|")
        for res in residues:
            aa1 = AA3_TO_1.get(res.resname, "?")
            catalytic_note = " **(catalytic)**" if res.resnum in CATALYTIC_RESIDUES else ""
            lines.append(
                f"| {res.resname} ({aa1}) | {res.chain}:{res.resnum} "
                f"| {res.min_distance:.2f}{catalytic_note} |"
            )

    lines.append("")

    # -- Risk Assessment --
    lines.append("## Risk Assessment")
    lines.append("")

    if results:
        lines.append(f"- Total binding site residues scanned: {stats['total_residues']}")
        lines.append(f"- Critical residues (delta > +{CRITICAL_THRESHOLD:.1f}): {stats['critical']}")
        lines.append(f"- Significant residues (delta > +{SIGNIFICANT_THRESHOLD:.1f}): {stats['significant']}")
        lines.append(
            f"- Of which catalytically essential (cannot mutate): "
            f"{stats['catalytic_critical']}"
        )
        lines.append(
            f"- Peripheral (could mutate): {stats['peripheral_critical']}"
        )
        lines.append(f"- **Overall resistance risk: {risk}**")
    else:
        lines.append("*Risk assessment requires completed alanine scanning.*")

    lines.append("")

    # -- Methodology --
    lines.append("## Methodology")
    lines.append("")
    lines.append(
        "1. Binding site residues were identified as protein atoms within "
        f"{CONTACT_DISTANCE:.1f} A of any MRL-003 atom in the docked pose."
    )
    lines.append(
        "2. Each binding site residue was mutated to alanine by removing "
        "side chain atoms beyond CB and renaming the residue."
    )
    lines.append(
        "3. Mutant structures were re-docked against MRL-003 using AutoDock "
        f"Vina (exhaustiveness={EXHAUSTIVENESS})."
    )
    lines.append(
        "4. Affinity changes were classified as CRITICAL (>{0:+.1f}), "
        "SIGNIFICANT (>{1:+.1f}), NEUTRAL, or BENEFICIAL (<{2:+.1f}).".format(
            CRITICAL_THRESHOLD,
            SIGNIFICANT_THRESHOLD,
            BENEFICIAL_THRESHOLD,
        )
    )
    lines.append(
        "5. Critical residues were cross-referenced with known catalytic "
        "residues of TryR (Cys52, Cys57, His461, Glu466) to determine "
        "whether resistance mutations are biologically viable."
    )
    lines.append("")

    report_text = "\n".join(lines)

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(report_text)

    logger.info("Wrote resistance report to %s.", output_path)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def predict_resistance(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the resistance prediction pipeline.

    1. Parse the docked MRL-003 pose and TryR PDB to identify binding site
       residues within 5 A of the ligand.
    2. Extract the wildtype binding affinity from the docking output.
    3. For each binding site residue, perform alanine scanning mutagenesis
       (mutate -> convert to PDBQT -> re-dock -> measure affinity change).
    4. Classify mutations by impact and cross-reference with conserved
       catalytic residues.
    5. Generate a resistance report and CSV.

    Args:
        force: If True, overwrite existing outputs.
        dry_run: If True, identify binding residues but skip mutagenesis
            and re-docking. Generate a partial report.

    Returns:
        Path to the resistance report Markdown file.
    """
    output_base = OUTPUT_DIR
    report_path = output_base / "resistance_report.md"
    csv_path = output_base / "resistance_scan.csv"

    if report_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            report_path,
        )
        return str(report_path)

    # ------------------------------------------------------------------
    # 1. Parse structures and find binding site residues
    # ------------------------------------------------------------------
    logger.info("Parsing ligand coordinates from %s ...", DOCKING_OUTPUT_PATH)
    ligand_coords = _parse_ligand_coords(DOCKING_OUTPUT_PATH)

    logger.info("Parsing protein structure from %s ...", TRYR_PDB_PATH)
    pdb_atoms = _parse_pdb_atoms(TRYR_PDB_PATH)

    binding_residues = _find_binding_residues(ligand_coords, pdb_atoms)

    if not binding_residues:
        logger.error(
            "No binding site residues found within %.1f A. "
            "Check input files and distance cutoff.",
            CONTACT_DISTANCE,
        )
        sys.exit(1)

    # ------------------------------------------------------------------
    # 2. Extract wildtype binding affinity
    # ------------------------------------------------------------------
    wt_affinity = _parse_wildtype_affinity(DOCKING_OUTPUT_PATH)
    logger.info("Wildtype MRL-003 binding affinity: %.3f kcal/mol.", wt_affinity)

    # ------------------------------------------------------------------
    # 3. Alanine scanning (skip in dry-run mode)
    # ------------------------------------------------------------------
    scan_results: list[ScanResult] = []

    if dry_run:
        logger.info(
            "[DRY RUN] Skipping alanine scanning for %d residues.",
            len(binding_residues),
        )
        for res in binding_residues:
            logger.info(
                "[DRY RUN] Would mutate %s%d:%s -> ALA and re-dock.",
                res.resname,
                res.resnum,
                res.chain,
            )
        risk = "N/A"
        stats: dict[str, int] = {
            "total_residues": len(binding_residues),
            "critical": 0,
            "significant": 0,
            "catalytic_critical": 0,
            "peripheral_critical": 0,
        }
    else:
        # Create a working directory for mutant structures.
        work_dir = output_base / "mutants"
        work_dir.mkdir(parents=True, exist_ok=True)

        scan_results = _run_alanine_scan(
            residues=binding_residues,
            wt_affinity=wt_affinity,
            pdb_path=TRYR_PDB_PATH,
            ligand_path=MRL003_LIGAND_PATH,
            work_dir=work_dir,
        )

        if not scan_results:
            logger.error("All alanine scanning jobs failed.")
            sys.exit(1)

        # 4. Assess resistance risk
        risk, stats = _assess_resistance_risk(scan_results)

        # Write CSV
        _write_scan_csv(scan_results, csv_path)

    # ------------------------------------------------------------------
    # 5. Generate resistance report
    # ------------------------------------------------------------------
    _write_resistance_report(
        results=scan_results,
        residues=binding_residues,
        risk=risk,
        stats=stats,
        wt_affinity=wt_affinity,
        output_path=report_path,
    )

    logger.info("Resistance prediction complete. Report: %s", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Predict resistance risk for MRL-003 against TryR via "
            "in-silico alanine scanning mutagenesis."
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
            "Identify binding residues but skip mutagenesis and re-docking. "
            "Generate a partial report with binding site information only."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = predict_resistance(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
