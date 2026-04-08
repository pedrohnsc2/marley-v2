"""Predict 3D protein structure for vaccine construct variants.

Runs ESMFold (or another backend configured in ``core.structure``) on each
vaccine construct variant, extracts per-residue pLDDT confidence scores,
maps functional regions onto the predicted structure, and generates
visualization scripts for PyMOL and ChimeraX.

Usage:
    python -m pipeline.07_structure
    python -m pipeline.07_structure --force
    python -m pipeline.07_structure --variant variant_A
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

from tqdm import tqdm

from core.logger import get_logger
from core.structure import (
    extract_plddt_scores,
    generate_chimerax_script,
    generate_pymol_script,
    map_construct_regions,
    predict_structure,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONSTRUCT_DIR = "results/construct"
STRUCTURE_SUBDIR = "structure"

log = get_logger("07_structure")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_fasta_sequence(fasta_path: Path) -> str:
    """Read a single-record FASTA file and return the plain amino-acid sequence."""
    lines: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines)


def _load_construct_card(card_path: Path) -> dict:
    """Load and return the construct_card.json as a dict."""
    with open(card_path) as fh:
        return json.load(fh)


def _write_json(data: object, path: Path) -> None:
    """Write *data* as pretty-printed JSON."""
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)
    log.info("Wrote %s", path)


def _write_confidence_csv(
    residues: list[int], scores: list[float], path: Path
) -> None:
    """Write a two-column CSV (residue, plddt)."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["residue", "plddt"])
        for res, score in zip(residues, scores):
            writer.writerow([res, round(score, 2)])
    log.info("Wrote %s", path)


# ---------------------------------------------------------------------------
# Single-variant prediction
# ---------------------------------------------------------------------------


def predict_variant_structure(
    variant_dir: str, force: bool = False
) -> dict | None:
    """Run structure prediction for one variant.

    Steps:
        1. Read ``vaccine_construct.fasta`` from *variant_dir*.
        2. Read ``construct_card.json`` for metadata.
        3. Create a ``structure/`` subdirectory.
        4. Call ``core.structure.predict_structure`` to obtain a PDB file.
        5. Call ``core.structure.map_construct_regions`` for region annotations.
        6. Call ``core.structure.extract_plddt_scores`` for per-residue confidence.
        7. Call ``core.structure.generate_pymol_script`` for a ``.pml`` file.
        8. Call ``core.structure.generate_chimerax_script`` for a ``.cxc`` file.
        9. Write ``regions.json`` (for the web viewer).
       10. Write ``confidence.csv`` (residue, plddt).
       11. Update ``construct_card.json`` with structure paths and mean pLDDT.
       12. Return a dict of generated file paths, or ``None`` on failure.

    If ``structure/predicted.pdb`` already exists and *force* is ``False``,
    the prediction is skipped and the existing paths are returned.
    """
    vdir = Path(variant_dir)
    variant_name = vdir.name

    fasta_path = vdir / "vaccine_construct.fasta"
    card_path = vdir / "construct_card.json"

    if not fasta_path.exists():
        log.warning("%s: vaccine_construct.fasta not found, skipping.", variant_name)
        return None
    if not card_path.exists():
        log.warning("%s: construct_card.json not found, skipping.", variant_name)
        return None

    struct_dir = vdir / STRUCTURE_SUBDIR
    pdb_path = struct_dir / "predicted.pdb"

    # --- Skip if already done ------------------------------------------------
    if pdb_path.exists() and not force:
        log.info(
            "%s: predicted.pdb already exists; use --force to re-run.", variant_name
        )
        return _collect_result_paths(struct_dir, variant_name)

    # --- Read inputs ---------------------------------------------------------
    sequence = _read_fasta_sequence(fasta_path)
    card = _load_construct_card(card_path)

    log.info(
        "%s: predicting structure for %d-residue construct ...",
        variant_name,
        len(sequence),
    )

    struct_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Predict structure ------------------------------------------------
    try:
        predict_structure(protein_seq=sequence, output_path=str(pdb_path))
    except Exception as exc:
        log.warning(
            "%s: structure prediction failed (%s: %s). Skipping.",
            variant_name,
            type(exc).__name__,
            exc,
        )
        return None

    if not pdb_path.exists():
        log.warning("%s: PDB file was not created. Skipping.", variant_name)
        return None

    log.info("%s: PDB written to %s", variant_name, pdb_path)

    # --- 2. Map construct regions --------------------------------------------
    regions_path = struct_dir / "regions.json"
    try:
        regions = map_construct_regions(
            protein_seq=sequence,
            signal_peptide_name=card.get("signal_peptide", "tPA"),
            adjuvant_name=card.get("adjuvant_name", "L7L12"),
            epitopes=card.get("epitopes", []),
        )
        regions_dicts = [
            r.to_dict() if hasattr(r, "to_dict") else r for r in regions
        ]
        _write_json(regions_dicts, regions_path)
    except Exception as exc:
        log.warning(
            "%s: region mapping failed (%s: %s). Continuing without regions.",
            variant_name,
            type(exc).__name__,
            exc,
        )
        regions_dicts = []

    # --- 3. Extract pLDDT scores ---------------------------------------------
    confidence_path = struct_dir / "confidence.csv"
    mean_plddt: float | None = None
    try:
        pdb_content = pdb_path.read_text()
        plddt_tuples = extract_plddt_scores(pdb_content)
        residues = [t[0] for t in plddt_tuples]
        scores = [t[1] for t in plddt_tuples]
        _write_confidence_csv(residues, scores, confidence_path)
        mean_plddt = round(sum(scores) / len(scores), 2) if scores else None
        log.info("%s: mean pLDDT = %s", variant_name, mean_plddt)
    except Exception as exc:
        log.warning(
            "%s: pLDDT extraction failed (%s: %s).",
            variant_name,
            type(exc).__name__,
            exc,
        )

    # --- 4. Generate PyMOL script --------------------------------------------
    pymol_path = struct_dir / "visualize.pml"
    try:
        generate_pymol_script(
            pdb_filename="predicted.pdb",
            regions=regions if regions else [],
            output_path=str(pymol_path),
        )
        log.info("%s: PyMOL script written to %s", variant_name, pymol_path)
    except Exception as exc:
        log.warning(
            "%s: PyMOL script generation failed (%s: %s).",
            variant_name,
            type(exc).__name__,
            exc,
        )

    # --- 5. Generate ChimeraX script -----------------------------------------
    chimerax_path = struct_dir / "visualize.cxc"
    try:
        generate_chimerax_script(
            pdb_filename="predicted.pdb",
            regions=regions if regions else [],
            output_path=str(chimerax_path),
        )
        log.info("%s: ChimeraX script written to %s", variant_name, chimerax_path)
    except Exception as exc:
        log.warning(
            "%s: ChimeraX script generation failed (%s: %s).",
            variant_name,
            type(exc).__name__,
            exc,
        )

    # --- 6. Update construct_card.json ---------------------------------------
    card["structure"] = {
        "predicted_pdb": str(pdb_path.relative_to(vdir)),
        "regions_json": str(regions_path.relative_to(vdir)),
        "confidence_csv": str(confidence_path.relative_to(vdir)),
        "pymol_script": str(pymol_path.relative_to(vdir)),
        "chimerax_script": str(chimerax_path.relative_to(vdir)),
        "mean_plddt": mean_plddt,
        "predicted_at": datetime.now().isoformat(),
    }
    _write_json(card, card_path)

    return _collect_result_paths(struct_dir, variant_name)


def _collect_result_paths(struct_dir: Path, variant_name: str) -> dict:
    """Return a summary dict of generated file paths for *variant_name*."""
    return {
        "variant": variant_name,
        "pdb": str(struct_dir / "predicted.pdb"),
        "regions": str(struct_dir / "regions.json"),
        "confidence": str(struct_dir / "confidence.csv"),
        "pymol": str(struct_dir / "visualize.pml"),
        "chimerax": str(struct_dir / "visualize.cxc"),
    }


# ---------------------------------------------------------------------------
# All-variants prediction
# ---------------------------------------------------------------------------


def predict_all_structures(force: bool = False) -> list[dict]:
    """Run structure prediction for every variant in *CONSTRUCT_DIR*.

    Discovers directories matching ``variant_*``, runs
    :func:`predict_variant_structure` on each, and returns a list of result
    dicts (one per variant that succeeded).
    """
    construct_root = Path(CONSTRUCT_DIR)
    if not construct_root.exists():
        log.error("Construct directory %s does not exist.", CONSTRUCT_DIR)
        return []

    variant_dirs = sorted(construct_root.glob("variant_*"))
    if not variant_dirs:
        log.warning("No variant_* directories found in %s.", CONSTRUCT_DIR)
        return []

    log.info("Found %d variant(s) to process.", len(variant_dirs))

    results: list[dict] = []
    for vdir in tqdm(variant_dirs, desc="Structure prediction", unit="variant"):
        result = predict_variant_structure(str(vdir), force=force)
        if result is not None:
            results.append(result)

    succeeded = len(results)
    total = len(variant_dirs)
    log.info(
        "Structure prediction complete: %d/%d variant(s) succeeded.",
        succeeded,
        total,
    )
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 7: Predict 3D structure for vaccine construct variants.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run prediction even if predicted.pdb already exists.",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="Run only this variant directory name (e.g. variant_A).",
    )
    args = parser.parse_args()

    if args.variant:
        vpath = Path(CONSTRUCT_DIR) / args.variant
        if not vpath.is_dir():
            log.error("Variant directory %s does not exist.", vpath)
            sys.exit(1)
        result = predict_variant_structure(str(vpath), force=args.force)
        if result is None:
            log.error("Structure prediction failed for %s.", args.variant)
            sys.exit(1)
        log.info("Result: %s", json.dumps(result, indent=2))
    else:
        results = predict_all_structures(force=args.force)
        if not results:
            log.warning("No structures were predicted successfully.")
            sys.exit(1)
        for r in results:
            log.info("  %s: %s", r["variant"], r["pdb"])
