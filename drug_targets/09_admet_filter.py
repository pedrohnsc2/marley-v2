"""ADMET filtering and composite scoring of docked compounds.

Reads docking results from module 08, applies Lipinski Rule of 5 filtering
via RDKit, optionally queries pkCSM for ADMET predictions, computes a
composite ranking score, and exports filtered results to CSV.

Usage:
    python -m drug_targets.09_admet_filter
    python -m drug_targets.09_admet_filter --predict-admet
    python -m drug_targets.09_admet_filter --force
    python -m drug_targets.09_admet_filter --dry-run
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESULTS_DIR: Final[Path] = Path("results")
DOCKING_SCORES_FILE: Final[Path] = RESULTS_DIR / "docking_scores.csv"
COMPOUND_LIBRARY_FILE: Final[Path] = RESULTS_DIR / "compound_library.csv"
OUTPUT_FILE: Final[Path] = RESULTS_DIR / "admet_filtered.csv"

PKCSM_API_URL: Final[str] = (
    "https://biosig.lab.uq.edu.au/pkcsm/prediction"
)
PKCSM_TIMEOUT: Final[int] = 60
PKCSM_TOP_N: Final[int] = 20

LIPINSKI_MAX_MW: Final[float] = 500.0
LIPINSKI_MAX_LOGP: Final[float] = 5.0
LIPINSKI_MAX_HBD: Final[int] = 5
LIPINSKI_MAX_HBA: Final[int] = 10
LIPINSKI_MAX_VIOLATIONS: Final[int] = 1

AFFINITY_WEIGHT: Final[float] = 0.50
ADMET_WEIGHT: Final[float] = 0.25
REPURPOSING_WEIGHT: Final[float] = 0.15
SELECTIVITY_WEIGHT: Final[float] = 0.10
AFFINITY_NORMALIZER: Final[float] = 12.0

CSV_COLUMNS: Final[list[str]] = [
    "target_gene_id",
    "target_gene_name",
    "compound_id",
    "compound_name",
    "smiles",
    "binding_affinity",
    "lipinski_violations",
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "admet_score",
    "is_approved_drug",
    "composite_score",
]

logger = get_logger("admet_filter")

# ---------------------------------------------------------------------------
# RDKit availability
# ---------------------------------------------------------------------------

_RDKIT_AVAILABLE: bool = False

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    _RDKIT_AVAILABLE = True
except ImportError:
    logger.warning(
        "RDKit not available. Lipinski filtering will be skipped; "
        "all compounds pass through with default descriptors."
    )


# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------


def _load_docking_scores(docking_path: Path) -> list[dict[str, str]]:
    """Read the docking scores CSV produced by module 08.

    Args:
        docking_path: Path to ``results/docking_scores.csv``.

    Returns:
        List of row dictionaries from the CSV.

    Raises:
        FileNotFoundError: If the docking scores file does not exist.
    """
    if not docking_path.exists():
        msg = f"Docking scores file not found: {docking_path}"
        raise FileNotFoundError(msg)

    rows: list[dict[str, str]] = []
    with open(docking_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)

    logger.info("Loaded %d docking results from %s.", len(rows), docking_path)
    return rows


def _load_compound_library(library_path: Path) -> dict[str, dict[str, str]]:
    """Load the compound library CSV and index by compound_id.

    The compound library contains SMILES and is_approved_drug fields that
    are not present in the docking scores CSV.

    Args:
        library_path: Path to ``results/compound_library.csv``.

    Returns:
        Dictionary mapping compound_id to row dictionaries.
    """
    index: dict[str, dict[str, str]] = {}

    if not library_path.exists():
        logger.warning(
            "Compound library not found at %s. "
            "SMILES and approval status will be unavailable.",
            library_path,
        )
        return index

    with open(library_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            cid = row.get("compound_id", "")
            if cid:
                index[cid] = row

    logger.info(
        "Loaded %d compounds from library %s.", len(index), library_path,
    )
    return index


# ---------------------------------------------------------------------------
# Lipinski Rule of 5
# ---------------------------------------------------------------------------


def _compute_lipinski(smiles: str) -> dict[str, Any]:
    """Compute Lipinski Rule of 5 properties and additional descriptors.

    If RDKit is not available, returns default values with zero violations
    so the compound passes through the filter.

    Args:
        smiles: SMILES string of the compound.

    Returns:
        Dictionary with keys: mw, logp, hbd, hba, tpsa,
        rotatable_bonds, fsp3, lipinski_violations.
    """
    defaults: dict[str, Any] = {
        "mw": 0.0,
        "logp": 0.0,
        "hbd": 0,
        "hba": 0,
        "tpsa": 0.0,
        "rotatable_bonds": 0,
        "fsp3": 0.0,
        "lipinski_violations": 0,
    }

    if not _RDKIT_AVAILABLE:
        return defaults

    if not smiles or not smiles.strip():
        return defaults

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.debug("RDKit could not parse SMILES: %s", smiles)
        return defaults

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    fsp3 = Lipinski.FractionCSP3(mol)

    violations = 0
    if mw > LIPINSKI_MAX_MW:
        violations += 1
    if logp > LIPINSKI_MAX_LOGP:
        violations += 1
    if hbd > LIPINSKI_MAX_HBD:
        violations += 1
    if hba > LIPINSKI_MAX_HBA:
        violations += 1

    return {
        "mw": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotatable_bonds,
        "fsp3": round(fsp3, 3),
        "lipinski_violations": violations,
    }


def _passes_lipinski(violations: int) -> bool:
    """Check whether a compound passes Lipinski with at most 1 violation.

    Args:
        violations: Number of Lipinski Rule of 5 violations.

    Returns:
        True if the compound has at most LIPINSKI_MAX_VIOLATIONS violations.
    """
    return violations <= LIPINSKI_MAX_VIOLATIONS


# ---------------------------------------------------------------------------
# ADMET prediction via pkCSM
# ---------------------------------------------------------------------------


def _predict_admet_pkcsm(smiles: str) -> dict[str, Any] | None:
    """Query the pkCSM API for ADMET predictions on a single compound.

    This function sends a POST request to the pkCSM prediction endpoint.
    If the service is unavailable or the request fails, returns None.

    Args:
        smiles: SMILES string of the compound.

    Returns:
        Dictionary with prediction results, or None on failure.
    """
    try:
        import requests
    except ImportError:
        logger.warning("requests library not available. Skipping pkCSM query.")
        return None

    try:
        resp = requests.post(
            PKCSM_API_URL,
            data={"smiles": smiles},
            timeout=PKCSM_TIMEOUT,
        )
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        logger.debug("pkCSM request failed for SMILES %s: %s", smiles[:50], exc)
        return None


def _compute_admet_score_from_pkcsm(prediction: dict[str, Any]) -> float:
    """Compute a weighted ADMET score from pkCSM prediction results.

    Evaluates three ADMET checks:
    - Intestinal absorption (>= 30% is a pass)
    - BBB permeability (log BB > -1.0 is a pass, but we treat BBB
      permeability cautiously -- not crossing BBB is acceptable)
    - AMES toxicity (negative / non-mutagenic is a pass)

    The score is the weighted average of passed checks, scaled to [0.0, 1.0].

    Args:
        prediction: Parsed JSON response from pkCSM.

    Returns:
        ADMET score between 0.0 and 1.0.
    """
    checks_passed = 0
    total_checks = 3

    # Intestinal absorption: value >= 30 is acceptable.
    absorption = prediction.get("intestinal_absorption")
    if absorption is not None:
        try:
            if float(absorption) >= 30.0:
                checks_passed += 1
        except (ValueError, TypeError):
            pass
    else:
        # If not available, give benefit of the doubt.
        checks_passed += 1

    # BBB permeability: for non-CNS drugs, NOT crossing BBB is fine.
    # We pass this check by default (leishmaniasis drugs are peripheral).
    bbb = prediction.get("bbb_permeability")
    if bbb is not None:
        # Pass regardless -- peripheral target.
        checks_passed += 1
    else:
        checks_passed += 1

    # AMES toxicity: "No" or 0 means non-mutagenic (pass).
    ames = prediction.get("ames_toxicity")
    if ames is not None:
        ames_str = str(ames).strip().lower()
        if ames_str in ("no", "0", "false", "non-mutagenic"):
            checks_passed += 1
    else:
        # If not available, give partial credit.
        checks_passed += 0

    return round(checks_passed / total_checks, 2)


def _assign_default_admet_score(lipinski_violations: int) -> float:
    """Assign a default ADMET score based on Lipinski violation count.

    Used when --predict-admet is not specified.

    Args:
        lipinski_violations: Number of Lipinski Rule of 5 violations.

    Returns:
        Default ADMET score (0.8, 0.6, or 0.3).
    """
    if lipinski_violations == 0:
        return 0.8
    if lipinski_violations == 1:
        return 0.6
    return 0.3


def _assign_dry_run_admet_score(is_approved_drug: bool) -> float:
    """Assign a default ADMET score for dry-run mode.

    Args:
        is_approved_drug: Whether the compound is an approved drug.

    Returns:
        0.8 for approved drugs, 0.5 otherwise.
    """
    return 0.8 if is_approved_drug else 0.5


# ---------------------------------------------------------------------------
# Composite score
# ---------------------------------------------------------------------------


def _compute_composite_score(
    binding_affinity: float,
    admet_score: float,
    is_approved_drug: bool,
) -> float:
    """Compute the composite ranking score for a compound.

    Formula:
        composite = norm_affinity * 0.50
                  + admet_score * 0.25
                  + repurposing_bonus * 0.15
                  + 0.10

    Where norm_affinity = min(1.0, max(0.0, -binding_affinity / 12.0)).

    Args:
        binding_affinity: Docking binding affinity in kcal/mol (negative).
        admet_score: ADMET score between 0.0 and 1.0.
        is_approved_drug: Whether the compound is an approved drug.

    Returns:
        Composite score between 0.0 and 1.0.
    """
    norm_affinity = min(1.0, max(0.0, -binding_affinity / AFFINITY_NORMALIZER))
    repurposing_bonus = 1.0 if is_approved_drug else 0.0

    composite = (
        norm_affinity * AFFINITY_WEIGHT
        + admet_score * ADMET_WEIGHT
        + repurposing_bonus * REPURPOSING_WEIGHT
        + SELECTIVITY_WEIGHT
    )

    return round(composite, 4)


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def _write_csv(
    rows: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write ADMET-filtered results to CSV, sorted by composite_score DESC.

    Args:
        rows: List of result dictionaries to write.
        output_path: Destination file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sorted_rows = sorted(
        rows,
        key=lambda r: float(r["composite_score"]),
        reverse=True,
    )

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow({col: row.get(col, "") for col in CSV_COLUMNS})

    logger.info("Wrote %d filtered compounds to %s.", len(sorted_rows), output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def filter_admet(
    predict_admet: bool = False,
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Filter docked compounds by ADMET properties and compute composite scores.

    Reads docking results from module 08, enriches with SMILES and approval
    status from the compound library, applies Lipinski Rule of 5 filtering
    via RDKit, optionally queries pkCSM for ADMET predictions, computes a
    composite ranking score, and writes filtered results to CSV.

    Args:
        predict_admet: If True, query pkCSM for ADMET predictions on the
            top 20 Lipinski-passing compounds. Otherwise use default scores.
        force: Overwrite the output CSV even if it already exists.
        dry_run: Skip RDKit and pkCSM; assign default ADMET scores based
            on approval status and write CSV.

    Returns:
        Path to the output CSV file (``results/admet_filtered.csv``).

    Raises:
        FileNotFoundError: If the docking scores file does not exist.
        RuntimeError: If no compounds remain after filtering.
    """
    output_path = OUTPUT_FILE

    if output_path.exists() and not force:
        logger.info(
            "Output file already exists at %s. Use --force to rerun.",
            output_path,
        )
        return str(output_path)

    # ------------------------------------------------------------------
    # Step 1: Load docking results and compound library
    # ------------------------------------------------------------------

    docking_rows = _load_docking_scores(DOCKING_SCORES_FILE)
    compound_index = _load_compound_library(COMPOUND_LIBRARY_FILE)

    if not docking_rows:
        msg = (
            "No docking results found. Run module 08 (docking) first."
        )
        raise RuntimeError(msg)

    # ------------------------------------------------------------------
    # Step 2: Enrich rows with SMILES and approval status
    # ------------------------------------------------------------------

    enriched: list[dict[str, Any]] = []

    for row in docking_rows:
        compound_id = row.get("compound_id", "")
        lib_entry = compound_index.get(compound_id, {})

        smiles = lib_entry.get("smiles", "")
        is_approved_str = lib_entry.get("is_approved_drug", "False")
        is_approved = str(is_approved_str).strip().lower() in (
            "true", "1", "yes",
        )

        enriched.append({
            "target_gene_id": row.get("target_gene_id", ""),
            "target_gene_name": row.get("target_gene_name", ""),
            "compound_id": compound_id,
            "compound_name": row.get("compound_name", compound_id),
            "smiles": smiles,
            "binding_affinity": float(row.get("binding_affinity", 0.0)),
            "is_approved_drug": is_approved,
        })

    # ------------------------------------------------------------------
    # Step 3: Dry run -- skip RDKit/pkCSM, assign defaults, write CSV
    # ------------------------------------------------------------------

    if dry_run:
        logger.info(
            "[DRY RUN] Skipping Lipinski and ADMET checks for %d compounds.",
            len(enriched),
        )

        results: list[dict[str, Any]] = []
        for entry in enriched:
            admet_score = _assign_dry_run_admet_score(entry["is_approved_drug"])
            composite = _compute_composite_score(
                binding_affinity=entry["binding_affinity"],
                admet_score=admet_score,
                is_approved_drug=entry["is_approved_drug"],
            )
            results.append({
                **entry,
                "lipinski_violations": 0,
                "mw": 0.0,
                "logp": 0.0,
                "hbd": 0,
                "hba": 0,
                "tpsa": 0.0,
                "admet_score": admet_score,
                "composite_score": composite,
            })

        _write_csv(results, output_path)

        if results:
            best = max(results, key=lambda r: r["composite_score"])
            logger.info(
                "[DRY RUN] %d compounds written. "
                "Top composite: %s (score=%.4f).",
                len(results),
                best["compound_name"],
                best["composite_score"],
            )

        return str(output_path)

    # ------------------------------------------------------------------
    # Step 4: Lipinski filtering
    # ------------------------------------------------------------------

    lipinski_passed: list[dict[str, Any]] = []
    lipinski_failed: list[dict[str, Any]] = []

    for entry in enriched:
        props = _compute_lipinski(entry["smiles"])
        entry.update({
            "mw": props["mw"],
            "logp": props["logp"],
            "hbd": props["hbd"],
            "hba": props["hba"],
            "tpsa": props["tpsa"],
            "lipinski_violations": props["lipinski_violations"],
        })

        if _passes_lipinski(props["lipinski_violations"]):
            lipinski_passed.append(entry)
        else:
            lipinski_failed.append(entry)

    logger.info(
        "%d of %d compounds pass Lipinski filter.",
        len(lipinski_passed),
        len(enriched),
    )

    if not lipinski_passed:
        logger.warning(
            "No compounds passed Lipinski filter. "
            "Including all compounds with penalty scores."
        )
        lipinski_passed = lipinski_failed + lipinski_passed

    # ------------------------------------------------------------------
    # Step 5: ADMET scoring
    # ------------------------------------------------------------------

    if predict_admet:
        # Sort by binding affinity (best first) and take top N.
        candidates = sorted(
            lipinski_passed,
            key=lambda r: r["binding_affinity"],
        )[:PKCSM_TOP_N]

        pkcsm_success = 0
        pkcsm_failed = 0

        for entry in candidates:
            smiles = entry.get("smiles", "")
            if not smiles:
                entry["admet_score"] = _assign_default_admet_score(
                    entry["lipinski_violations"],
                )
                continue

            prediction = _predict_admet_pkcsm(smiles)
            if prediction is not None:
                entry["admet_score"] = _compute_admet_score_from_pkcsm(
                    prediction,
                )
                pkcsm_success += 1
            else:
                entry["admet_score"] = _assign_default_admet_score(
                    entry["lipinski_violations"],
                )
                pkcsm_failed += 1

        # Assign defaults for compounds not sent to pkCSM.
        pkcsm_ids = {e["compound_id"] for e in candidates}
        for entry in lipinski_passed:
            if entry["compound_id"] not in pkcsm_ids:
                entry["admet_score"] = _assign_default_admet_score(
                    entry["lipinski_violations"],
                )

        logger.info(
            "pkCSM predictions: %d succeeded, %d failed (out of %d).",
            pkcsm_success,
            pkcsm_failed,
            len(candidates),
        )
    else:
        for entry in lipinski_passed:
            entry["admet_score"] = _assign_default_admet_score(
                entry["lipinski_violations"],
            )

    # ------------------------------------------------------------------
    # Step 6: Compute composite scores
    # ------------------------------------------------------------------

    results = []
    for entry in lipinski_passed:
        composite = _compute_composite_score(
            binding_affinity=entry["binding_affinity"],
            admet_score=entry["admet_score"],
            is_approved_drug=entry["is_approved_drug"],
        )
        entry["composite_score"] = composite
        results.append(entry)

    if not results:
        msg = "No compounds remain after ADMET filtering."
        raise RuntimeError(msg)

    # ------------------------------------------------------------------
    # Step 7: Write output CSV
    # ------------------------------------------------------------------

    _write_csv(results, output_path)

    # ------------------------------------------------------------------
    # Step 8: Summary
    # ------------------------------------------------------------------

    best = max(results, key=lambda r: r["composite_score"])
    logger.info(
        "%d of %d compounds pass Lipinski filter. "
        "Top composite: %s (score=%.4f).",
        len([r for r in results if r["lipinski_violations"] <= LIPINSKI_MAX_VIOLATIONS]),
        len(enriched),
        best["compound_name"],
        best["composite_score"],
    )

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Filter docked compounds by Lipinski/ADMET properties and "
            "compute composite ranking scores for L. infantum drug targets."
        ),
    )
    parser.add_argument(
        "--predict-admet",
        action="store_true",
        help=(
            "Query pkCSM for ADMET predictions on top Lipinski-passing "
            "compounds. Without this flag, default scores are assigned."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun even if admet_filtered.csv already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Skip RDKit and pkCSM; assign default ADMET scores and "
            "write CSV for review."
        ),
    )
    args = parser.parse_args()

    result_path = filter_admet(
        predict_admet=args.predict_admet,
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Complete: %s", result_path)
