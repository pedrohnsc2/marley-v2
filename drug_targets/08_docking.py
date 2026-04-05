"""Molecular docking campaign using AutoDock Vina.

Scans prepared receptor and ligand PDBQT files, defines search grid boxes
(known binding sites or blind docking via centroid), executes Vina docking
in parallel, and exports ranked binding affinities to CSV.

Usage:
    python -m drug_targets.08_docking
    python -m drug_targets.08_docking --top-n 10 --exhaustiveness 32
    python -m drug_targets.08_docking --blind-docking
    python -m drug_targets.08_docking --dry-run
"""

from __future__ import annotations

import csv
import random
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STRUCTURES_DIR: Final[Path] = Path("data/structures")
COMPOUNDS_DIR: Final[Path] = Path("data/compounds")
DOCKING_OUTPUT_DIR: Final[Path] = Path("data/docking")
RESULTS_DIR: Final[Path] = Path("results")
OUTPUT_FILE: Final[Path] = RESULTS_DIR / "docking_scores.csv"

CSV_COLUMNS: Final[list[str]] = [
    "target_gene_id",
    "target_gene_name",
    "compound_id",
    "compound_name",
    "binding_affinity",
    "rmsd_lb",
    "rmsd_ub",
    "pdbqt_path",
]

BLIND_BOX_SIZE: Final[float] = 40.0
DEFAULT_N_POSES: Final[int] = 5

KNOWN_BINDING_SITES: Final[dict[str, dict[str, float]]] = {
    "TryR": {
        "center_x": 11.0, "center_y": 24.0, "center_z": 15.0,
        "size_x": 25.0, "size_y": 25.0, "size_z": 25.0,
    },
    "TryS": {
        "center_x": 20.0, "center_y": 30.0, "center_z": 20.0,
        "size_x": 30.0, "size_y": 30.0, "size_z": 30.0,
    },
    "ADL": {
        "center_x": 15.0, "center_y": 20.0, "center_z": 18.0,
        "size_x": 25.0, "size_y": 25.0, "size_z": 25.0,
    },
    "SMT": {
        "center_x": 18.0, "center_y": 22.0, "center_z": 16.0,
        "size_x": 25.0, "size_y": 25.0, "size_z": 25.0,
    },
    "6PGDH": {
        "center_x": 12.0, "center_y": 18.0, "center_z": 14.0,
        "size_x": 25.0, "size_y": 25.0, "size_z": 25.0,
    },
}

logger = get_logger("docking")

# ---------------------------------------------------------------------------
# Grid box helpers
# ---------------------------------------------------------------------------


def _compute_centroid(pdbqt_path: Path) -> dict[str, float]:
    """Parse PDBQT ATOM records and compute the protein centroid.

    Returns a grid box dictionary centered on the centroid with a uniform
    size of 40 Angstroms per axis, suitable for blind docking.

    Args:
        pdbqt_path: Path to a receptor PDBQT file.

    Returns:
        Dictionary with keys center_x, center_y, center_z, size_x,
        size_y, size_z.

    Raises:
        ValueError: If no ATOM records are found in the file.
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
        "center_x": sum(xs) / len(xs),
        "center_y": sum(ys) / len(ys),
        "center_z": sum(zs) / len(zs),
        "size_x": BLIND_BOX_SIZE,
        "size_y": BLIND_BOX_SIZE,
        "size_z": BLIND_BOX_SIZE,
    }


def _get_grid_box(
    gene_name: str,
    receptor_path: Path,
    blind_docking: bool,
) -> dict[str, float]:
    """Determine the docking grid box for a given target.

    Uses pre-defined coordinates for well-studied targets unless
    blind docking is requested. Falls back to centroid computation
    for unknown targets.

    Args:
        gene_name: Short name of the target gene (e.g. "TryR").
        receptor_path: Path to the receptor PDBQT file.
        blind_docking: If True, always compute centroid-based box.

    Returns:
        Dictionary with center and size coordinates.
    """
    if not blind_docking and gene_name in KNOWN_BINDING_SITES:
        logger.debug(
            "Using known binding site for %s.", gene_name,
        )
        return KNOWN_BINDING_SITES[gene_name]

    logger.info(
        "Computing centroid grid box for %s (blind_docking=%s).",
        gene_name,
        blind_docking,
    )
    return _compute_centroid(receptor_path)


# ---------------------------------------------------------------------------
# Target / ligand discovery
# ---------------------------------------------------------------------------


def _extract_gene_name(receptor_filename: str) -> str:
    """Extract the gene name from a receptor PDBQT filename.

    Expects filenames of the form ``<gene_name>_receptor.pdbqt``.

    Args:
        receptor_filename: Stem or full name of the receptor file.

    Returns:
        The gene name portion of the filename.
    """
    return receptor_filename.replace("_receptor.pdbqt", "")


def _discover_pairs(
    top_n: int,
) -> list[dict[str, str | Path]]:
    """Scan data directories for receptor-ligand pairs.

    Args:
        top_n: Maximum number of ligands per target to include.

    Returns:
        List of dicts with keys: target_gene_name, receptor_path,
        ligand_path, compound_id.

    Raises:
        FileNotFoundError: If the structures directory does not exist.
    """
    if not STRUCTURES_DIR.exists():
        msg = f"Structures directory not found: {STRUCTURES_DIR}"
        raise FileNotFoundError(msg)

    receptors = sorted(STRUCTURES_DIR.glob("*_receptor.pdbqt"))
    if not receptors:
        logger.warning("No receptor PDBQT files found in %s.", STRUCTURES_DIR)
        return []

    pairs: list[dict[str, str | Path]] = []

    for receptor_path in receptors:
        full_name = _extract_gene_name(receptor_path.name)
        # Try full name first (e.g. TryS_Q4QJ30), then gene-only (e.g. TryS).
        gene_short = full_name.split("_")[0] if "_" in full_name else full_name
        ligand_dir = COMPOUNDS_DIR / full_name
        gene_name = full_name

        if not ligand_dir.exists():
            ligand_dir = COMPOUNDS_DIR / gene_short
            gene_name = gene_short

        if not ligand_dir.exists():
            logger.warning(
                "No compounds directory for target %s at %s.",
                full_name,
                COMPOUNDS_DIR / full_name,
            )
            continue

        ligands = sorted(ligand_dir.glob("*.pdbqt"))[:top_n]
        if not ligands:
            logger.warning(
                "No ligand PDBQT files in %s.", ligand_dir,
            )
            continue

        for ligand_path in ligands:
            compound_id = ligand_path.stem
            pairs.append({
                "target_gene_name": gene_name,
                "receptor_path": receptor_path,
                "ligand_path": ligand_path,
                "compound_id": compound_id,
            })

    logger.info(
        "Discovered %d receptor-ligand pairs across %d receptors.",
        len(pairs),
        len(receptors),
    )
    return pairs


# ---------------------------------------------------------------------------
# Vina execution
# ---------------------------------------------------------------------------


def _parse_vina_cli_output(stdout: str) -> tuple[float, float, float]:
    """Extract best-pose affinity and RMSD values from Vina CLI output.

    Args:
        stdout: Standard output captured from the vina command.

    Returns:
        Tuple of (binding_affinity, rmsd_lb, rmsd_ub) for the top pose.

    Raises:
        RuntimeError: If no scoring lines are found in the output.
    """
    pattern = re.compile(
        r"^\s*1\s+"
        r"([-+]?\d+\.?\d*)\s+"
        r"([-+]?\d+\.?\d*)\s+"
        r"([-+]?\d+\.?\d*)",
    )
    for line in stdout.splitlines():
        match = pattern.match(line)
        if match:
            return (
                float(match.group(1)),
                float(match.group(2)),
                float(match.group(3)),
            )

    msg = "Could not parse docking scores from Vina CLI output."
    raise RuntimeError(msg)


def _run_single_docking(
    receptor_path: Path,
    ligand_path: Path,
    output_path: Path,
    grid_box: dict[str, float],
    exhaustiveness: int,
) -> dict[str, float]:
    """Execute a single Vina docking job.

    Attempts to use the Vina Python bindings first. If not installed,
    falls back to the ``vina`` command-line tool.

    Args:
        receptor_path: Path to the receptor PDBQT file.
        ligand_path: Path to the ligand PDBQT file.
        output_path: Path to write the docked pose PDBQT.
        grid_box: Dictionary with center_x/y/z and size_x/y/z.
        exhaustiveness: Vina exhaustiveness parameter.

    Returns:
        Dictionary with binding_affinity, rmsd_lb, rmsd_ub.

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
        v.dock(exhaustiveness=exhaustiveness, n_poses=DEFAULT_N_POSES)

        energies = v.energies(n_poses=1)
        best = energies[0]
        v.write_poses(str(output_path), n_poses=1)

        return {
            "binding_affinity": float(best[0]),
            "rmsd_lb": float(best[1]),
            "rmsd_ub": float(best[2]),
        }
    except ImportError:
        logger.debug(
            "Vina Python bindings not available, falling back to CLI.",
        )
    except Exception as exc:
        logger.warning(
            "Vina Python bindings failed for %s vs %s: %s. Trying CLI.",
            receptor_path.name,
            ligand_path.name,
            exc,
        )

    # --- Attempt 2: Vina CLI ---
    cmd = [
        "vina",
        "--receptor", str(receptor_path),
        "--ligand", str(ligand_path),
        "--center_x", str(grid_box["center_x"]),
        "--center_y", str(grid_box["center_y"]),
        "--center_z", str(grid_box["center_z"]),
        "--size_x", str(grid_box["size_x"]),
        "--size_y", str(grid_box["size_y"]),
        "--size_z", str(grid_box["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(DEFAULT_N_POSES),
        "--out", str(output_path),
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=600,
        )
        affinity, rmsd_lb, rmsd_ub = _parse_vina_cli_output(result.stdout)
        return {
            "binding_affinity": affinity,
            "rmsd_lb": rmsd_lb,
            "rmsd_ub": rmsd_ub,
        }
    except FileNotFoundError:
        msg = (
            "Neither Vina Python bindings nor the 'vina' CLI were found. "
            "Install meeko/vina or add vina to PATH."
        )
        raise RuntimeError(msg)
    except subprocess.CalledProcessError as exc:
        msg = (
            f"Vina CLI failed for {receptor_path.name} vs {ligand_path.name}: "
            f"{exc.stderr}"
        )
        raise RuntimeError(msg)
    except subprocess.TimeoutExpired:
        msg = (
            f"Vina CLI timed out for {receptor_path.name} vs "
            f"{ligand_path.name} (600s limit)."
        )
        raise RuntimeError(msg)


def _dock_pair(
    target_gene_name: str,
    receptor_path: Path,
    ligand_path: Path,
    compound_id: str,
    grid_box: dict[str, float],
    exhaustiveness: int,
) -> dict[str, str | float]:
    """Dock a single receptor-ligand pair and return the result dict.

    This function is the unit of work submitted to the process pool.

    Args:
        target_gene_name: Name of the target gene.
        receptor_path: Path to receptor PDBQT.
        ligand_path: Path to ligand PDBQT.
        compound_id: Identifier for the compound.
        grid_box: Docking grid box parameters.
        exhaustiveness: Vina exhaustiveness parameter.

    Returns:
        Dictionary with target_gene_id, target_gene_name, compound_id,
        compound_name, binding_affinity, rmsd_lb, rmsd_ub, pdbqt_path.
    """
    output_dir = DOCKING_OUTPUT_DIR / target_gene_name
    output_path = output_dir / f"{compound_id}_out.pdbqt"

    scores = _run_single_docking(
        receptor_path=receptor_path,
        ligand_path=ligand_path,
        output_path=output_path,
        grid_box=grid_box,
        exhaustiveness=exhaustiveness,
    )

    return {
        "target_gene_id": f"DOCK_{target_gene_name}",
        "target_gene_name": target_gene_name,
        "compound_id": compound_id,
        "compound_name": compound_id,
        "binding_affinity": scores["binding_affinity"],
        "rmsd_lb": scores["rmsd_lb"],
        "rmsd_ub": scores["rmsd_ub"],
        "pdbqt_path": str(output_path),
    }


# ---------------------------------------------------------------------------
# Results handling
# ---------------------------------------------------------------------------


def _write_csv(
    results: list[dict[str, str | float]],
    output_path: Path,
) -> None:
    """Write docking results to CSV sorted by binding affinity (ascending).

    Args:
        results: List of result dictionaries from docking jobs.
        output_path: Path to write the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sorted_results = sorted(
        results,
        key=lambda r: float(r["binding_affinity"]),
    )

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in sorted_results:
            writer.writerow({col: row[col] for col in CSV_COLUMNS})

    logger.info("Wrote %d docking results to %s.", len(sorted_results), output_path)


def _persist_to_database(results: list[dict[str, str | float]]) -> None:
    """Upsert docking results to Supabase.

    Failures are logged but do not halt the pipeline.

    Args:
        results: List of result dictionaries from docking jobs.
    """
    try:
        from core.db import upsert_docking_result
        from core.models import DockingResult
    except ImportError:
        logger.warning("Database modules not available, skipping persistence.")
        return

    persisted = 0
    for row in results:
        try:
            dr = DockingResult(
                target_gene_id=str(row["target_gene_id"]),
                target_gene_name=str(row["target_gene_name"]),
                compound_id=str(row["compound_id"]),
                compound_name=str(row["compound_name"]),
                binding_affinity=float(row["binding_affinity"]),
                rmsd_lb=float(row["rmsd_lb"]),
                rmsd_ub=float(row["rmsd_ub"]),
                pdbqt_path=str(row["pdbqt_path"]),
            )
            upsert_docking_result(dr)
            persisted += 1
        except Exception:
            logger.error(
                "Failed to persist result for %s vs %s — continuing.",
                row["target_gene_name"],
                row["compound_id"],
            )

    logger.info("Persisted %d / %d results to Supabase.", persisted, len(results))


# ---------------------------------------------------------------------------
# Dry run
# ---------------------------------------------------------------------------


def _generate_dry_run_results(
    pairs: list[dict[str, str | Path]],
) -> list[dict[str, str | float]]:
    """Generate fake docking scores for a dry run.

    Args:
        pairs: Discovered receptor-ligand pairs.

    Returns:
        List of result dictionaries with random binding affinities.
    """
    rng = random.Random(42)  # noqa: S311 — reproducible fake scores
    results: list[dict[str, str | float]] = []

    for pair in pairs:
        gene_name = str(pair["target_gene_name"])
        compound_id = str(pair["compound_id"])
        output_path = DOCKING_OUTPUT_DIR / gene_name / f"{compound_id}_out.pdbqt"

        results.append({
            "target_gene_id": f"DOCK_{gene_name}",
            "target_gene_name": gene_name,
            "compound_id": compound_id,
            "compound_name": compound_id,
            "binding_affinity": round(rng.uniform(-8.0, -4.0), 2),
            "rmsd_lb": 0.0,
            "rmsd_ub": 0.0,
            "pdbqt_path": str(output_path),
        })

    return results


def _generate_dry_run_from_csv() -> list[dict[str, str | float]]:
    """Generate fake docking results from the compound library CSV.

    Used when PDBQT files have not been generated yet (full dry-run).

    Returns:
        List of result dictionaries with random binding affinities.
    """
    compound_csv = Path("results/compound_library.csv")
    results: list[dict[str, str | float]] = []
    rng = random.Random(42)  # noqa: S311

    if not compound_csv.exists():
        logger.warning("No compound library CSV found for dry-run fallback.")
        return results

    with open(compound_csv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            target = row.get("target_gene_name", "unknown")
            cpd_id = row.get("compound_id", "unknown")
            cpd_name = row.get("name", cpd_id)
            results.append({
                "target_gene_id": f"DOCK_{target}",
                "target_gene_name": target,
                "compound_id": cpd_id,
                "compound_name": cpd_name,
                "binding_affinity": round(rng.uniform(-8.0, -4.0), 2),
                "rmsd_lb": 0.0,
                "rmsd_ub": 0.0,
                "pdbqt_path": "",
            })

    logger.info("Generated %d fake docking results from compound library.", len(results))
    return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_docking_campaign(
    top_n: int = 5,
    exhaustiveness: int = 16,
    max_workers: int = 4,
    blind_docking: bool = False,
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Execute a full molecular docking campaign against all prepared targets.

    Scans for receptor and ligand PDBQT files, defines docking grid boxes,
    runs AutoDock Vina (via Python bindings or CLI) in parallel, and exports
    ranked results to CSV.

    Args:
        top_n: Maximum number of ligands to dock per target.
        exhaustiveness: Vina exhaustiveness parameter (higher = more thorough).
        max_workers: Number of parallel docking processes.
        blind_docking: If True, use centroid-based grid boxes for all targets
            instead of pre-defined binding site coordinates.
        force: Re-dock and overwrite even if the output CSV already exists.
        dry_run: List pairs and generate fake scores without running Vina.

    Returns:
        Path to the output CSV file (``results/docking_scores.csv``).

    Raises:
        FileNotFoundError: If the structures directory does not exist.
        RuntimeError: If no receptor-ligand pairs are found.
    """
    output_path = OUTPUT_FILE

    if output_path.exists() and not force:
        logger.info(
            "Output file already exists at %s. Use --force to re-dock.",
            output_path,
        )
        return str(output_path)

    # Step 1: Discover receptor-ligand pairs.
    pairs = _discover_pairs(top_n=top_n)
    if not pairs and not dry_run:
        msg = (
            "No receptor-ligand pairs found. Ensure PDBQT files exist in "
            f"{STRUCTURES_DIR} and {COMPOUNDS_DIR}."
        )
        raise RuntimeError(msg)

    logger.info(
        "Docking campaign: %d pairs, exhaustiveness=%d, max_workers=%d, "
        "blind_docking=%s.",
        len(pairs),
        exhaustiveness,
        max_workers,
        blind_docking,
    )

    # Step 2: Dry run — generate fake results and exit early.
    if dry_run:
        if pairs:
            logger.info("[DRY RUN] Listing %d receptor-ligand pairs:", len(pairs))
            for pair in pairs:
                logger.info(
                    "[DRY RUN]   %s <-> %s",
                    Path(str(pair["receptor_path"])).name,
                    Path(str(pair["ligand_path"])).name,
                )
            results = _generate_dry_run_results(pairs)
        else:
            # No PDBQT files yet — generate fake pairs from compound library.
            logger.info("[DRY RUN] No PDBQT files found. Generating fake results from compound library.")
            results = _generate_dry_run_from_csv()

        _write_csv(results, output_path)

        logger.info(
            "[DRY RUN] Wrote %d fake docking scores to %s.",
            len(results),
            output_path,
        )
        return str(output_path)

    # Step 3: Pre-compute grid boxes for each target.
    grid_boxes: dict[str, dict[str, float]] = {}
    for pair in pairs:
        gene_name = str(pair["target_gene_name"])
        if gene_name not in grid_boxes:
            grid_boxes[gene_name] = _get_grid_box(
                gene_name=gene_name,
                receptor_path=Path(str(pair["receptor_path"])),
                blind_docking=blind_docking,
            )

    # Step 4: Execute docking in parallel.
    results: list[dict[str, str | float]] = []
    failed = 0

    try:
        from tqdm import tqdm
        progress_bar = tqdm(total=len(pairs), desc="Docking", unit="pair")
    except ImportError:
        progress_bar = None

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for pair in pairs:
            gene_name = str(pair["target_gene_name"])
            future = executor.submit(
                _dock_pair,
                target_gene_name=gene_name,
                receptor_path=Path(str(pair["receptor_path"])),
                ligand_path=Path(str(pair["ligand_path"])),
                compound_id=str(pair["compound_id"]),
                grid_box=grid_boxes[gene_name],
                exhaustiveness=exhaustiveness,
            )
            futures[future] = pair

        for future in as_completed(futures):
            pair = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception:
                failed += 1
                logger.error(
                    "Docking failed for %s vs %s.",
                    Path(str(pair["receptor_path"])).name,
                    Path(str(pair["ligand_path"])).name,
                    exc_info=True,
                )
            finally:
                if progress_bar is not None:
                    progress_bar.update(1)

    if progress_bar is not None:
        progress_bar.close()

    if not results:
        msg = f"All {len(pairs)} docking jobs failed."
        raise RuntimeError(msg)

    if failed:
        logger.warning("%d / %d docking jobs failed.", failed, len(pairs))

    # Step 5: Save results.
    _write_csv(results, output_path)
    _persist_to_database(results)

    # Step 6: Log summary.
    best = min(results, key=lambda r: float(r["binding_affinity"]))
    logger.info(
        "Completed %d dockings. Best hit: %s against %s at %.2f kcal/mol.",
        len(results),
        best["compound_id"],
        best["target_gene_name"],
        float(best["binding_affinity"]),
    )

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Run molecular docking campaign against L. infantum drug targets "
            "using AutoDock Vina."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Maximum number of ligands to dock per target (default: 5).",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=16,
        help="Vina exhaustiveness parameter (default: 16).",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Number of parallel docking processes (default: 4).",
    )
    parser.add_argument(
        "--blind-docking",
        action="store_true",
        help="Use centroid-based grid boxes instead of known binding sites.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-dock even if docking_scores.csv already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List pairs and generate fake scores without running Vina.",
    )
    args = parser.parse_args()

    result_path = run_docking_campaign(
        top_n=args.top_n,
        exhaustiveness=args.exhaustiveness,
        max_workers=args.max_workers,
        blind_docking=args.blind_docking,
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Complete: %s", result_path)
