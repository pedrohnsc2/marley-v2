"""Predict RNA secondary structures and compute minimum free energy.

Uses ViennaRNA Python bindings when available, with a GC-content-based
fallback for environments where ViennaRNA is not installed.  Computes
the final ``information_score`` for each RNA target and saves scored
results for downstream report generation.

Usage:
    python -m rna_entropy.06_structure_prediction
    python -m rna_entropy.06_structure_prediction --force
    python -m rna_entropy.06_structure_prediction --dry-run
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MERGED_TARGETS_PATH: Final[str] = "results/rna/rna_targets_merged.json"
STRUCTURE_CSV_PATH: Final[str] = "results/rna/structure_predictions.csv"
SCORED_TARGETS_PATH: Final[str] = "results/rna/rna_targets_scored.json"

INFO_SCORE_THRESHOLD: Final[float] = 0.70
TOP_N_FALLBACK: Final[int] = 20

# Weights for the information_score formula.
W_ENTROPY_DELTA: Final[float] = 0.30
W_CODON_BIAS: Final[float] = 0.20
W_SL_RNA: Final[float] = 0.15
W_MFE: Final[float] = 0.15
W_CONSERVATION: Final[float] = 0.20

# Empirical coefficient for GC-based MFE estimation.
MFE_GC_COEFFICIENT: Final[float] = -0.3

logger = get_logger("rna_structure")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _estimate_mfe(sequence: str, gc_content: float) -> float:
    """Approximate MFE from GC content and sequence length.

    Higher GC content implies more stable secondary structures and thus
    a more negative minimum free energy.

    Args:
        sequence: The RNA nucleotide sequence.
        gc_content: Fraction of G+C bases (0.0 -- 1.0).

    Returns:
        Estimated MFE in kcal/mol (negative value).
    """
    length = len(sequence) if sequence else 100
    return MFE_GC_COEFFICIENT * length * gc_content


def _fold_sequence(sequence: str, gc_content: float, use_fallback: bool = False) -> tuple[float, str]:
    """Predict secondary structure and MFE for a single RNA sequence.

    Attempts to use ViennaRNA Python bindings first.  Falls back to the
    GC-content-based estimation when the library is unavailable or when
    *use_fallback* is ``True``.

    Args:
        sequence: The RNA nucleotide sequence.
        gc_content: GC fraction for fallback estimation.
        use_fallback: If ``True``, skip ViennaRNA and use estimation.

    Returns:
        Tuple of ``(mfe, dot_bracket)`` where *dot_bracket* is ``""``
        when using the fallback estimator.
    """
    if use_fallback or not sequence:
        return _estimate_mfe(sequence, gc_content), ""

    try:
        import RNA  # type: ignore[import-untyped]

        structure, mfe = RNA.fold(sequence)
        return float(mfe), structure
    except ImportError:
        logger.info("ViennaRNA not available -- using GC-based MFE estimation.")
        return _estimate_mfe(sequence, gc_content), ""


def _compute_information_score(target: dict[str, Any]) -> float:
    """Compute the composite information score for a target.

    The score is a weighted combination of:
        - entropy_delta (normalised to [0, 1] by dividing by 2.0)
        - codon_bias_score
        - SL RNA presence (1.0 if True, 0.0 otherwise)
        - MFE stability (normalised: min(-mfe, 50) / 50)
        - conservation_score

    Args:
        target: Dictionary with target metrics.

    Returns:
        A score between 0.0 and 1.0.
    """
    entropy_delta = min(target.get("entropy_delta", 0.0) / 2.0, 1.0)
    codon_bias = target.get("codon_bias_score", 0.0)
    sl_rna = 1.0 if target.get("has_sl_rna", False) else 0.0
    mfe_raw = target.get("min_free_energy", 0.0)
    mfe_norm = min(abs(mfe_raw), 50.0) / 50.0
    conservation = target.get("conservation_score", 0.0)

    score = (
        W_ENTROPY_DELTA * entropy_delta
        + W_CODON_BIAS * codon_bias
        + W_SL_RNA * sl_rna
        + W_MFE * mfe_norm
        + W_CONSERVATION * conservation
    )
    return round(min(score, 1.0), 4)


def _filter_targets(targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Select targets eligible for structure prediction.

    Prefers targets with ``information_score > 0.70``.  Falls back to the
    top 20 targets sorted by ``entropy_delta`` when no pre-computed scores
    are available.

    Args:
        targets: Full list of merged RNA target dictionaries.

    Returns:
        Filtered subset of targets.
    """
    scored = [t for t in targets if t.get("information_score", 0.0) > INFO_SCORE_THRESHOLD]
    if scored:
        return scored

    logger.info(
        "No targets with information_score > %.2f found; "
        "selecting top %d by entropy_delta.",
        INFO_SCORE_THRESHOLD,
        TOP_N_FALLBACK,
    )
    sorted_targets = sorted(
        targets,
        key=lambda t: t.get("entropy_delta", 0.0),
        reverse=True,
    )
    return sorted_targets[:TOP_N_FALLBACK]


def _save_structure_csv(targets: list[dict[str, Any]], output_path: Path) -> None:
    """Write structure predictions to a CSV file.

    Args:
        targets: List of target dicts with structure prediction fields.
        output_path: Destination CSV path.
    """
    fieldnames: Final[list[str]] = [
        "gene_id",
        "gene_name",
        "length",
        "gc_content",
        "min_free_energy",
        "dot_bracket",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for t in targets:
            seq = t.get("sequence_rna", "")
            writer.writerow(
                {
                    "gene_id": t.get("gene_id", ""),
                    "gene_name": t.get("gene_name", ""),
                    "length": len(seq) if seq else 0,
                    "gc_content": round(t.get("gc_content", 0.0), 4),
                    "min_free_energy": round(t.get("min_free_energy", 0.0), 2),
                    "dot_bracket": t.get("dot_bracket", ""),
                }
            )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def predict_rna_structures(force: bool = False, dry_run: bool = False) -> str:
    """Predict RNA secondary structures and score all targets.

    Loads merged RNA targets, predicts secondary structures via ViennaRNA
    (or a GC-based fallback), computes the final ``information_score``,
    and writes both a CSV of structure predictions and a scored JSON for
    downstream report generation.

    Args:
        force: Re-run predictions even if output files already exist.
        dry_run: Use MFE estimation fallback for all targets and skip
                 ViennaRNA.  No external dependencies required.

    Returns:
        Path to the scored targets JSON file.

    Raises:
        FileNotFoundError: If the merged targets file does not exist.
    """
    merged_path = Path(MERGED_TARGETS_PATH)
    csv_path = Path(STRUCTURE_CSV_PATH)
    scored_path = Path(SCORED_TARGETS_PATH)

    # --- Dry run mode --------------------------------------------------------
    if dry_run:
        logger.info("[DRY RUN] Would load targets from %s", merged_path)
        logger.info("[DRY RUN] Using GC-based MFE estimation (skipping ViennaRNA)")

        if merged_path.exists():
            with open(merged_path, "r", encoding="utf-8") as fh:
                targets: list[dict[str, Any]] = json.load(fh)
            selected = _filter_targets(targets)
        else:
            logger.info("[DRY RUN] Merged targets file not found; using empty set.")
            selected = []

        for target in selected:
            seq = target.get("sequence_rna", "")
            gc = target.get("gc_content", 0.5)
            mfe, dot = _fold_sequence(seq, gc, use_fallback=True)
            target["min_free_energy"] = round(mfe, 2)
            target["dot_bracket"] = dot
            target["information_score"] = _compute_information_score(target)

        if selected:
            most_stable = min(selected, key=lambda t: t.get("min_free_energy", 0.0))
            logger.info(
                "[DRY RUN] Predicted structures for %d targets. "
                "Most stable: %s (MFE = %.1f kcal/mol)",
                len(selected),
                most_stable.get("gene_name", "unknown"),
                most_stable.get("min_free_energy", 0.0),
            )
        else:
            logger.info("[DRY RUN] No targets to predict.")

        return str(scored_path)

    # --- Check for existing files --------------------------------------------
    if csv_path.exists() and scored_path.exists() and not force:
        logger.info(
            "Structure predictions already exist at %s. Use --force to re-run.",
            csv_path,
        )
        return str(scored_path)

    # --- Load merged targets -------------------------------------------------
    if not merged_path.exists():
        raise FileNotFoundError(
            f"Merged targets file not found: {merged_path}. "
            "Run stages 01-05 first."
        )

    with open(merged_path, "r", encoding="utf-8") as fh:
        all_targets: list[dict[str, Any]] = json.load(fh)

    logger.info("Loaded %d merged targets from %s", len(all_targets), merged_path)

    selected = _filter_targets(all_targets)
    logger.info("Selected %d targets for structure prediction.", len(selected))

    # --- Predict structures --------------------------------------------------
    for target in selected:
        seq = target.get("sequence_rna", "")
        gc = target.get("gc_content", 0.5)
        mfe, dot = _fold_sequence(seq, gc)
        target["min_free_energy"] = round(mfe, 2)
        target["dot_bracket"] = dot

    # --- Compute information scores for ALL targets --------------------------
    for target in all_targets:
        if "min_free_energy" not in target:
            seq = target.get("sequence_rna", "")
            gc = target.get("gc_content", 0.5)
            target["min_free_energy"] = round(
                _estimate_mfe(seq, gc), 2
            )
        target["information_score"] = _compute_information_score(target)

    # --- Save structure CSV --------------------------------------------------
    _save_structure_csv(selected, csv_path)
    logger.info("Saved structure predictions to %s", csv_path)

    # --- Save scored targets JSON --------------------------------------------
    scored_path.parent.mkdir(parents=True, exist_ok=True)
    with open(scored_path, "w", encoding="utf-8") as fh:
        json.dump(all_targets, fh, indent=2, ensure_ascii=False)
    logger.info("Saved scored targets to %s", scored_path)

    # --- Summary log ---------------------------------------------------------
    if selected:
        most_stable = min(selected, key=lambda t: t.get("min_free_energy", 0.0))
        logger.info(
            "Predicted structures for %d targets. "
            "Most stable: %s (MFE = %.1f kcal/mol)",
            len(selected),
            most_stable.get("gene_name", "unknown"),
            most_stable.get("min_free_energy", 0.0),
        )

    return str(scored_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Predict RNA secondary structures and compute information scores.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run predictions even if output files already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Use GC-based MFE estimation only (no ViennaRNA).",
    )
    args = parser.parse_args()

    result = predict_rna_structures(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
