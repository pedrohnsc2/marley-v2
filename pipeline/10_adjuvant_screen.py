"""Screen adjuvant candidates for the Marley vaccine construct.

Evaluates a curated library of adjuvants by building test constructs with
safe epitopes and scoring them via VaxiJen for predicted antigenicity.
Protein adjuvants are assembled into a simplified construct; non-protein
adjuvants (CpG ODN, MPLA) are listed as co-administered formulation
components with literature-based recommendations.

Usage:
    python -m pipeline.10_adjuvant_screen
    python -m pipeline.10_adjuvant_screen --force
    python -m pipeline.10_adjuvant_screen --dry-run
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Final

import requests

from core.logger import get_logger

logger = get_logger("adjuvant_screen")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SAFETY_REPORT: Final[str] = "results/v4_optimization/safety_report.json"
OUTPUT_FILE: Final[str] = "results/v4_optimization/adjuvant_ranking.json"

SIGNAL_PEPTIDE_TPA: Final[str] = "MDAMKRGLCCVLLLCGAVFVSAS"

LINKER_ADJUVANT: Final[str] = "EAAAK"
LINKER_CTL: Final[str] = "AAY"
LINKER_HTL: Final[str] = "GPGPG"

VAXIJEN_URL: Final[str] = (
    "http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html"
)
VAXIJEN_THRESHOLD: Final[float] = 0.4
VAXIJEN_TIMEOUT: Final[int] = 30

DRY_RUN_SCORES: Final[dict[str, float | None]] = {
    "L7L12": 0.65,
    "RS09": 0.55,
    "beta_defensin_2": 0.60,
    "flagellin_D0D1": 0.58,
    "CpG_ODN_2006": None,
    "MPLA": None,
}

ADJUVANTS: Final[dict[str, dict[str, Any]]] = {
    "L7L12": {
        "sequence": (
            "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAEE"
            "QSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAK"
            "AKLEAAGATVTVK"
        ),
        "th_bias": "Th1",
        "is_protein": True,
        "notes": (
            "50S ribosomal protein. Strong Th1 inducer via IFN-gamma. "
            "Used in v1."
        ),
    },
    "RS09": {
        "sequence": "APPHALS",
        "th_bias": "Th1",
        "is_protein": True,
        "notes": "Synthetic TLR4 agonist peptide. Short, minimal. Used in v1.",
    },
    "beta_defensin_2": {
        "sequence": (
            "MRINTLILVCLLVQACSLAGIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRK"
            "CCRRK"
        ),
        "th_bias": "Th1",
        "is_protein": True,
        "notes": (
            "Canine beta-defensin-2. Chemoattractant for dendritic cells. "
            "Species-matched."
        ),
    },
    "flagellin_D0D1": {
        "sequence": (
            "MAQVINTNSLSLLTQNNLNKSQSALGTAIERLSSGLRINSAKDDAAGQAIANRFTANI"
            "KGLTQASRNANDGISIAQTTEGALNEINNNLQRVRELAVQSANSTNSQSDLDSIQAEI"
            "TQRLNEID"
        ),
        "th_bias": "mixed",
        "is_protein": True,
        "notes": (
            "Salmonella FliC D0-D1 domains. TLR5 agonist. Strong innate "
            "response."
        ),
    },
    "CpG_ODN_2006": {
        "sequence": "TCGTCGTTTTGTCGTTTTGTCGTT",
        "th_bias": "Th1",
        "is_protein": False,
        "notes": (
            "TLR9 agonist DNA. Strong Th1 in dogs. Co-administered, not "
            "fused."
        ),
    },
    "MPLA": {
        "sequence": "",
        "th_bias": "Th1",
        "is_protein": False,
        "notes": (
            "Monophosphoryl lipid A. TLR4 agonist. Formulation component, "
            "not fused."
        ),
    },
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_safe_epitopes(safety_path: str = SAFETY_REPORT) -> list[dict[str, Any]]:
    """Load epitopes marked as safe from the safety report.

    Args:
        safety_path: Path to the safety_report.json file.

    Returns:
        List of epitope dicts where ``is_safe`` is ``True``.

    Raises:
        FileNotFoundError: If the safety report does not exist.
        ValueError: If no safe epitopes are found.
    """
    path = Path(safety_path)
    if not path.exists():
        msg = f"Safety report not found: {safety_path}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    with open(path, "r") as fh:
        report = json.load(fh)

    # Support both top-level list and {"epitopes": [...]} formats.
    epitopes: list[dict[str, Any]]
    if isinstance(report, list):
        epitopes = report
    elif isinstance(report, dict) and "epitopes" in report:
        epitopes = report["epitopes"]
    else:
        msg = "Unrecognised safety report format."
        logger.error(msg)
        raise ValueError(msg)

    safe = [ep for ep in epitopes if ep.get("is_safe", False)]

    if not safe:
        msg = "No safe epitopes found in the safety report."
        logger.error(msg)
        raise ValueError(msg)

    logger.info(
        "Loaded %d safe epitope(s) from %d total in %s.",
        len(safe),
        len(epitopes),
        safety_path,
    )
    return safe


# ---------------------------------------------------------------------------
# Construct assembly
# ---------------------------------------------------------------------------


def _classify_epitope(epitope: dict[str, Any]) -> str:
    """Classify an epitope as CTL or HTL based on available metadata.

    Heuristic: CTL epitopes are typically 8-11 aa; HTL epitopes are 15+ aa.
    If the epitope dict contains an ``epitope_type`` key, that is used
    directly.

    Args:
        epitope: Single epitope dict from the safety report.

    Returns:
        ``"CTL"`` or ``"HTL"``.
    """
    if "epitope_type" in epitope:
        return epitope["epitope_type"].upper()

    peptide = epitope.get("peptide", epitope.get("sequence", ""))
    return "CTL" if len(peptide) <= 14 else "HTL"


def _get_peptide(epitope: dict[str, Any]) -> str:
    """Extract the peptide sequence from an epitope dict.

    Tries ``peptide`` first, falls back to ``sequence``.

    Args:
        epitope: Single epitope dict.

    Returns:
        Amino-acid sequence string.
    """
    return epitope.get("peptide", epitope.get("sequence", ""))


def _build_construct(
    adjuvant_name: str,
    adjuvant_seq: str,
    safe_epitopes: list[dict[str, Any]],
) -> str:
    """Assemble a simplified test construct for VaxiJen scoring.

    Layout::

        [tPA signal] + [adjuvant] + EAAAK +
        [CTL epitopes joined by AAY] + [HTL epitopes joined by GPGPG]

    Args:
        adjuvant_name: Name of the adjuvant (for logging).
        adjuvant_seq: Amino-acid sequence of the adjuvant.
        safe_epitopes: Safe epitope dicts from the safety report.

    Returns:
        Full amino-acid sequence of the test construct.
    """
    ctl_peptides: list[str] = []
    htl_peptides: list[str] = []

    for ep in safe_epitopes:
        peptide = _get_peptide(ep)
        if not peptide:
            continue
        if _classify_epitope(ep) == "CTL":
            ctl_peptides.append(peptide)
        else:
            htl_peptides.append(peptide)

    parts: list[str] = [SIGNAL_PEPTIDE_TPA, adjuvant_seq, LINKER_ADJUVANT]

    if ctl_peptides:
        parts.append(LINKER_CTL.join(ctl_peptides))

    if htl_peptides:
        if ctl_peptides:
            parts.append(LINKER_ADJUVANT)
        parts.append(LINKER_HTL.join(htl_peptides))

    construct = "".join(parts)
    logger.info(
        "Built construct for %s: %d aa (%d CTL, %d HTL epitopes).",
        adjuvant_name,
        len(construct),
        len(ctl_peptides),
        len(htl_peptides),
    )
    return construct


# ---------------------------------------------------------------------------
# VaxiJen scoring
# ---------------------------------------------------------------------------


def _query_vaxijen(sequence: str) -> float | None:
    """Submit a protein sequence to VaxiJen and parse the antigenicity score.

    The VaxiJen web server returns an HTML page containing the predicted
    score.  This function posts the sequence and scrapes the numeric
    result from the response body.

    Args:
        sequence: Amino-acid sequence of the construct.

    Returns:
        Predicted antigenicity score, or ``None`` if the request fails.
    """
    payload = {
        "sequence": sequence,
        "organism": "parasite",
    }

    try:
        response = requests.post(
            VAXIJEN_URL,
            data=payload,
            timeout=VAXIJEN_TIMEOUT,
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("VaxiJen request failed: %s", exc)
        return None

    # VaxiJen returns the score in a line like:
    #   "Overall Prediction for the Protective Antigen= 0.6543 (Probable ANTIGEN)"
    match = re.search(
        r"Overall Prediction[^=]*=\s*([\d.]+)",
        response.text,
    )
    if match:
        score = float(match.group(1))
        logger.info("VaxiJen score: %.4f", score)
        return score

    logger.warning("Could not parse VaxiJen score from response.")
    return None


# ---------------------------------------------------------------------------
# Ranking and output
# ---------------------------------------------------------------------------


def _build_ranking(
    results: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Sort adjuvant results by VaxiJen score descending.

    Entries with ``None`` scores are placed at the end.

    Args:
        results: List of per-adjuvant result dicts.

    Returns:
        Sorted list (best score first).
    """
    scored = [r for r in results if r.get("vaxijen_score") is not None]
    unscored = [r for r in results if r.get("vaxijen_score") is None]

    scored.sort(key=lambda r: r["vaxijen_score"], reverse=True)
    return scored + unscored


def _save_ranking(
    ranking: list[dict[str, Any]],
    output_path: str = OUTPUT_FILE,
) -> str:
    """Write the adjuvant ranking to JSON.

    Args:
        ranking: Sorted list of adjuvant result dicts.
        output_path: Destination file path.

    Returns:
        Absolute path to the written file.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as fh:
        json.dump(ranking, fh, indent=2, default=str)

    logger.info("Adjuvant ranking saved to %s.", path)
    return str(path.resolve())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def screen_adjuvants(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Screen adjuvant candidates and rank by predicted antigenicity.

    Workflow:
        1. Load safe epitopes from the safety report.
        2. For each protein adjuvant, build a simplified test construct.
        3. Score the construct via VaxiJen (or use estimated scores in
           dry-run mode).
        4. List non-protein adjuvants as co-administered components.
        5. Rank by VaxiJen score and save to JSON.

    Args:
        force: Re-run even if the output file already exists.
        dry_run: Skip VaxiJen API calls and use estimated scores.

    Returns:
        Path to the saved adjuvant ranking JSON.
    """
    output_path = Path(OUTPUT_FILE)

    if output_path.exists() and not force:
        logger.info(
            "Adjuvant ranking already exists at %s. Use --force to re-run.",
            output_path,
        )
        return str(output_path.resolve())

    if dry_run:
        logger.info("DRY RUN: VaxiJen API calls will be skipped.")

    # -- Load safe epitopes -------------------------------------------------
    safe_epitopes = _load_safe_epitopes()

    # -- Screen each adjuvant -----------------------------------------------
    results: list[dict[str, Any]] = []

    for name, info in ADJUVANTS.items():
        logger.info("-" * 50)
        logger.info("Evaluating adjuvant: %s", name)

        entry: dict[str, Any] = {
            "adjuvant": name,
            "th_bias": info["th_bias"],
            "is_protein": info["is_protein"],
            "notes": info["notes"],
        }

        if not info["is_protein"]:
            # Non-protein adjuvant: co-administered, no construct scoring.
            entry["construct_length"] = None
            entry["vaxijen_score"] = None
            entry["is_antigenic"] = None
            entry["recommendation"] = "co-administered"
            logger.info(
                "%s is a non-protein adjuvant. Listed as co-administered.",
                name,
            )
            results.append(entry)
            continue

        # Build test construct for protein adjuvant.
        construct = _build_construct(name, info["sequence"], safe_epitopes)
        entry["construct_length"] = len(construct)

        # Score via VaxiJen.
        if dry_run:
            score = DRY_RUN_SCORES.get(name)
            logger.info(
                "DRY RUN: Assigned estimated VaxiJen score %.2f for %s.",
                score if score is not None else 0.0,
                name,
            )
        else:
            score = _query_vaxijen(construct)

        entry["vaxijen_score"] = score

        if score is not None:
            entry["is_antigenic"] = score > VAXIJEN_THRESHOLD
            entry["recommendation"] = (
                "fused" if score > VAXIJEN_THRESHOLD else "consider alternatives"
            )
        else:
            entry["is_antigenic"] = None
            entry["recommendation"] = "score unavailable"
            logger.warning(
                "VaxiJen score unavailable for %s. Cannot assess antigenicity.",
                name,
            )

        results.append(entry)

    # -- Rank and save ------------------------------------------------------
    ranking = _build_ranking(results)
    saved_path = _save_ranking(ranking)

    # -- Summary log --------------------------------------------------------
    protein_count = sum(1 for r in ranking if r["is_protein"])
    coadmin_count = sum(1 for r in ranking if not r["is_protein"])

    best = next(
        (r for r in ranking if r.get("vaxijen_score") is not None),
        None,
    )
    if best is not None:
        logger.info(
            "Tested %d adjuvants. Best: %s (VaxiJen=%.4f, %s bias).",
            protein_count + coadmin_count,
            best["adjuvant"],
            best["vaxijen_score"],
            best["th_bias"],
        )
    else:
        logger.info(
            "Tested %d adjuvants. No VaxiJen scores available.",
            protein_count + coadmin_count,
        )

    logger.info(
        "Protein adjuvants scored: %d. Co-administered: %d.",
        protein_count,
        coadmin_count,
    )

    return saved_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Step 10: Screen adjuvant candidates for the Marley vaccine "
            "construct and rank by VaxiJen antigenicity score."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if the output file already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Skip VaxiJen API calls and use estimated scores for a fast "
            "local test."
        ),
    )
    args = parser.parse_args()

    result_path = screen_adjuvants(force=args.force, dry_run=args.dry_run)
    print(f"Adjuvant ranking: {result_path}")
