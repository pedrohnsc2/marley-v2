"""Optimize vaccine epitopes via Altered Peptide Ligand (APL) design.

Generates single-residue anchor-position variants of each construct epitope,
scores them against IEDB MHC-I binding predictions, and keeps only those
variants that improve predicted IC50 over the original peptide.

Usage:
    python -m pipeline.08_epitope_optimize
    python -m pipeline.08_epitope_optimize --force
    python -m pipeline.08_epitope_optimize --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

logger = get_logger("epitope_optimize")

CONSTRUCT_CARD_PATH: Final[str] = "results/construct/construct_card.json"
CONSTRUCT_CARD_FALLBACKS: Final[list[str]] = [
    "results/construct/variant_A/construct_card.json",
    "results/construct/variant_B/construct_card.json",
    "results/construct/variant_C/construct_card.json",
]

OUTPUT_DIR: Final[str] = "results/v4_optimization"
OUTPUT_FILE: Final[str] = "optimized_epitopes.json"

IEDB_MHCI_URL: Final[str] = (
    "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
)

# Preferred anchor residues for DLA class-I (9-mer CTL epitopes).
# Position indices are 0-based.
CTL_ANCHOR_POS_2_RESIDUES: Final[list[str]] = ["L", "M", "I", "V", "A"]
CTL_ANCHOR_POS_9_RESIDUES: Final[list[str]] = ["L", "V", "I", "F", "M"]
CTL_ANCHORS: Final[dict[int, list[str]]] = {
    1: CTL_ANCHOR_POS_2_RESIDUES,
    8: CTL_ANCHOR_POS_9_RESIDUES,
}

# Preferred anchor residues for DLA class-II (15-mer HTL epitopes).
HTL_HYDROPHOBIC_RESIDUES: Final[list[str]] = ["L", "I", "V", "F", "M", "A"]
HTL_ANCHORS: Final[dict[int, list[str]]] = {
    0: HTL_HYDROPHOBIC_RESIDUES,
    3: HTL_HYDROPHOBIC_RESIDUES,
    5: HTL_HYDROPHOBIC_RESIDUES,
    8: HTL_HYDROPHOBIC_RESIDUES,
}

IEDB_REQUEST_DELAY_S: Final[float] = 1.0
IEDB_MAX_RETRIES: Final[int] = 3
IEDB_RETRY_BACKOFF_S: Final[float] = 5.0


# ---------------------------------------------------------------------------
# Helpers -- I/O
# ---------------------------------------------------------------------------


def _find_construct_card() -> Path:
    """Locate the construct_card.json file, trying fallbacks if necessary."""
    primary = Path(CONSTRUCT_CARD_PATH)
    if primary.exists():
        return primary
    for fallback in CONSTRUCT_CARD_FALLBACKS:
        p = Path(fallback)
        if p.exists():
            logger.info("Primary card not found; using fallback %s", p)
            return p
    raise FileNotFoundError(
        f"construct_card.json not found at {primary} or any fallback location."
    )


def _load_construct_card(path: Path) -> dict:
    """Load and return the construct card JSON."""
    with open(path) as fh:
        return json.load(fh)


def _write_json(data: object, path: Path) -> None:
    """Write *data* as pretty-printed JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)
    logger.info("Wrote %s", path)


# ---------------------------------------------------------------------------
# Helpers -- Epitope classification
# ---------------------------------------------------------------------------


def _classify_epitope(peptide: str) -> str:
    """Return 'CTL' for 9-mers and 'HTL' for 15-mers.

    Raises ValueError for unexpected lengths.
    """
    length = len(peptide)
    if length == 9:
        return "CTL"
    if length == 15:
        return "HTL"
    raise ValueError(
        f"Unexpected peptide length {length} for '{peptide}'. "
        "Expected 9 (CTL) or 15 (HTL)."
    )


# ---------------------------------------------------------------------------
# Variant generation
# ---------------------------------------------------------------------------


def _generate_variants(peptide: str, epitope_type: str) -> list[str]:
    """Generate single-residue anchor-position APL variants of *peptide*.

    Each variant mutates exactly one anchor position to one preferred
    residue.  Variants identical to the original are skipped.

    Returns:
        A list of unique variant peptide strings.
    """
    anchors = CTL_ANCHORS if epitope_type == "CTL" else HTL_ANCHORS
    seen: set[str] = {peptide}
    variants: list[str] = []

    for pos, preferred_residues in anchors.items():
        if pos >= len(peptide):
            continue
        for aa in preferred_residues:
            if peptide[pos] == aa:
                continue
            variant = peptide[:pos] + aa + peptide[pos + 1 :]
            if variant not in seen:
                seen.add(variant)
                variants.append(variant)

    return variants


# ---------------------------------------------------------------------------
# IEDB scoring
# ---------------------------------------------------------------------------


def _query_iedb_mhci(
    peptide: str, allele: str, length: int = 9
) -> float | None:
    """Query the IEDB MHC-I binding prediction API.

    Returns the predicted IC50 (nM) or ``None`` on failure.
    """
    payload = {
        "method": "recommended",
        "sequence_text": peptide,
        "allele": allele,
        "length": str(length),
    }

    last_error: Exception | None = None
    for attempt in range(1, IEDB_MAX_RETRIES + 1):
        try:
            resp = requests.post(IEDB_MHCI_URL, data=payload, timeout=60)
            resp.raise_for_status()
            return _parse_iedb_response(resp.text, peptide)
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt < IEDB_MAX_RETRIES:
                wait = IEDB_RETRY_BACKOFF_S * attempt
                logger.warning(
                    "IEDB request failed (attempt %d/%d): %s. "
                    "Retrying in %.0fs ...",
                    attempt,
                    IEDB_MAX_RETRIES,
                    exc,
                    wait,
                )
                time.sleep(wait)

    logger.error(
        "IEDB request for %s failed after %d attempts: %s",
        peptide,
        IEDB_MAX_RETRIES,
        last_error,
    )
    return None


def _parse_iedb_response(text: str, peptide: str) -> float:
    """Extract the IC50 value from the IEDB tab-separated response.

    The IEDB API returns a header line followed by data lines.  The IC50
    column is typically the last column.  We look for the row matching
    *peptide* and return its IC50.

    Raises ValueError if parsing fails.
    """
    lines = text.strip().splitlines()
    if len(lines) < 2:
        raise ValueError(f"IEDB response too short: {text!r}")

    header = lines[0].split("\t")
    # Find the ic50 column (often labelled 'ic50' or last column)
    ic50_idx: int | None = None
    for idx, col in enumerate(header):
        if "ic50" in col.lower():
            ic50_idx = idx
            break
    if ic50_idx is None:
        # Fall back to last column
        ic50_idx = len(header) - 1

    for line in lines[1:]:
        fields = line.split("\t")
        if peptide in line and len(fields) > ic50_idx:
            try:
                return float(fields[ic50_idx])
            except ValueError:
                continue

    raise ValueError(
        f"Could not find IC50 for peptide '{peptide}' in IEDB response."
    )


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def _optimize_single_epitope(
    epitope: dict,
    dry_run: bool = False,
) -> dict | None:
    """Optimize a single epitope by generating and scoring APL variants.

    Returns a dict describing the best variant, or ``None`` if no
    improvement was found.
    """
    peptide: str = epitope["peptide"]
    allele: str = epitope["allele"]
    original_ic50: float = epitope["ic50"]
    epitope_type = _classify_epitope(peptide)

    variants = _generate_variants(peptide, epitope_type)
    if not variants:
        logger.debug("No variants generated for %s", peptide)
        return None

    logger.info(
        "Epitope %s (%s, IC50=%.2f): %d variants generated.",
        peptide,
        epitope_type,
        original_ic50,
        len(variants),
    )

    best_variant: str | None = None
    best_ic50: float = original_ic50
    all_scored: list[dict] = []

    for variant in variants:
        if dry_run:
            # Assign a synthetic IC50: 50% of original for the first variant,
            # then progressively worse for the rest.
            idx = variants.index(variant)
            fake_ic50 = original_ic50 * (0.5 + 0.1 * idx)
            scored = {
                "variant_peptide": variant,
                "predicted_ic50": round(fake_ic50, 2),
                "improvement_ratio": round(original_ic50 / fake_ic50, 4)
                if fake_ic50 > 0
                else 0.0,
            }
            all_scored.append(scored)
            if fake_ic50 < best_ic50:
                best_ic50 = fake_ic50
                best_variant = variant
        else:
            time.sleep(IEDB_REQUEST_DELAY_S)
            predicted_ic50 = _query_iedb_mhci(
                variant, allele, length=len(variant)
            )
            if predicted_ic50 is None:
                continue
            ratio = (
                round(original_ic50 / predicted_ic50, 4)
                if predicted_ic50 > 0
                else 0.0
            )
            scored = {
                "variant_peptide": variant,
                "predicted_ic50": round(predicted_ic50, 2),
                "improvement_ratio": ratio,
            }
            all_scored.append(scored)
            if predicted_ic50 < best_ic50:
                best_ic50 = predicted_ic50
                best_variant = variant

    # Filter to only improved variants
    improved = [v for v in all_scored if v["improvement_ratio"] > 1.0]
    improved.sort(key=lambda v: v["improvement_ratio"], reverse=True)

    if not best_variant or best_ic50 >= original_ic50:
        logger.info("  No improvement found for %s.", peptide)
        return None

    best_entry = improved[0] if improved else None
    if best_entry is None:
        return None

    return {
        "original_peptide": peptide,
        "original_ic50": original_ic50,
        "gene_id": epitope.get("gene_id", ""),
        "gene_name": epitope.get("gene_name", ""),
        "allele": allele,
        "epitope_type": epitope_type,
        "best_variant": best_entry["variant_peptide"],
        "best_variant_ic50": best_entry["predicted_ic50"],
        "improvement_ratio": best_entry["improvement_ratio"],
        "all_improved_variants": improved,
        "total_variants_tested": len(all_scored),
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def optimize_epitopes(
    force: bool = False, dry_run: bool = False
) -> str:
    """Run the epitope optimization pipeline.

    1. Load epitopes from construct_card.json.
    2. Generate APL variants for each epitope.
    3. Score variants via IEDB API (or fake scores in dry-run mode).
    4. Rank by improvement ratio and save results.

    Args:
        force: If True, overwrite existing output even if present.
        dry_run: If True, skip IEDB API calls and use synthetic IC50 values.

    Returns:
        Path to the output JSON file.
    """
    output_dir = Path(OUTPUT_DIR)
    output_path = output_dir / OUTPUT_FILE

    if output_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            output_path,
        )
        return str(output_path)

    # --- Load epitopes -------------------------------------------------------
    card_path = _find_construct_card()
    card = _load_construct_card(card_path)
    epitopes: list[dict] = card.get("epitopes", [])

    if not epitopes:
        logger.error("No epitopes found in %s.", card_path)
        sys.exit(1)

    # Deduplicate epitopes by peptide sequence to avoid redundant API calls
    seen_peptides: set[str] = set()
    unique_epitopes: list[dict] = []
    for ep in epitopes:
        pep = ep["peptide"]
        if pep not in seen_peptides:
            seen_peptides.add(pep)
            unique_epitopes.append(ep)

    logger.info(
        "Loaded %d epitopes (%d unique) from %s.",
        len(epitopes),
        len(unique_epitopes),
        card_path,
    )

    if dry_run:
        logger.info("DRY RUN mode: no IEDB API calls will be made.")

    # --- Optimize each epitope -----------------------------------------------
    results: list[dict] = []
    for ep in unique_epitopes:
        result = _optimize_single_epitope(ep, dry_run=dry_run)
        if result is not None:
            results.append(result)

    # --- Sort by improvement ratio -------------------------------------------
    results.sort(key=lambda r: r["improvement_ratio"], reverse=True)

    # --- Build output document -----------------------------------------------
    optimized_count = len(results)
    total_count = len(unique_epitopes)
    best_ratio = results[0]["improvement_ratio"] if results else 0.0
    best_sequence = results[0]["best_variant"] if results else "N/A"

    output = {
        "generated": datetime.now().isoformat(),
        "source_construct_card": str(card_path),
        "dry_run": dry_run,
        "total_epitopes": total_count,
        "optimized_epitopes": optimized_count,
        "best_improvement_ratio": best_ratio,
        "best_variant_peptide": best_sequence,
        "results": results,
    }

    _write_json(output, output_path)

    logger.info(
        "Optimized %d of %d epitopes. Best improvement: %s (%.2fx better).",
        optimized_count,
        total_count,
        best_sequence,
        best_ratio,
    )

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 8: Optimize vaccine epitopes via APL design.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate variants without calling IEDB API (synthetic IC50 values).",
    )
    args = parser.parse_args()

    result_path = optimize_epitopes(force=args.force, dry_run=args.dry_run)
    logger.info("Output: %s", result_path)
