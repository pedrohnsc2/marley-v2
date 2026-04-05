"""Score druggability of *Leishmania infantum* drug targets.

Upserts pre-validated priority targets, computes druggability scores for
computational targets using host-divergence / active-site / essentiality
weights, attaches AlphaFold structure URLs, and exports a ranked CSV.

Usage:
    python -m drug_targets.04_druggability
    python -m drug_targets.04_druggability --force
    python -m drug_targets.04_druggability --dry-run
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Final

from core.logger import get_logger
from core.models import DrugTarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESULTS_DIR: Final[Path] = Path("results")
OUTPUT_FILE: Final[Path] = RESULTS_DIR / "druggability_scores.csv"
ALPHAFOLD_BASE_URL: Final[str] = "https://alphafold.ebi.ac.uk/entry"
DRUGGABILITY_THRESHOLD: Final[float] = 0.60

CSV_COLUMNS: Final[list[str]] = [
    "gene_id",
    "gene_name",
    "pathway",
    "identity_score",
    "active_site_diff",
    "is_essential",
    "druggability_score",
    "alphafold_url",
    "priority",
]

logger = get_logger("druggability")

# ---------------------------------------------------------------------------
# Priority validated targets
# ---------------------------------------------------------------------------

PRIORITY_DRUG_TARGETS: Final[list[dict]] = [
    {
        "gene_id": "PRIORITY_ADL",
        "gene_name": "ADL",
        "sequence": "",
        "enzyme_class": "adenylosuccinate_lyase",
        "pathway": "purine_salvage",
        "evidence": "Unica enzima de metabolismo de purinas essencial sozinha. Homologia limitada com humano.",
        "druggability_score": 0.95,
        "is_essential": True,
        "identity_score": 0.25,
        "active_site_diff": 0.80,
        "priority": True,
        "status": "priority_validated",
    },
    {
        "gene_id": "PRIORITY_GMPS",
        "gene_name": "GMPS",
        "sequence": "",
        "enzyme_class": "GMP_synthetase",
        "pathway": "purine_salvage",
        "evidence": "Calcanhar de Aquiles da via de purinas. Ponto de convergencia obrigatorio.",
        "druggability_score": 0.92,
        "is_essential": True,
        "identity_score": 0.30,
        "active_site_diff": 0.75,
        "priority": True,
        "status": "priority_validated",
    },
    {
        "gene_id": "PRIORITY_TryS",
        "gene_name": "TryS",
        "sequence": "",
        "enzyme_class": "trypanothione_synthetase",
        "pathway": "trypanothione",
        "evidence": "Completamente ausente em humanos. Alvo prioritario global para tripanosomatideos.",
        "druggability_score": 0.98,
        "is_essential": True,
        "identity_score": 0.0,
        "active_site_diff": 1.0,
        "priority": True,
        "status": "priority_validated",
    },
    {
        "gene_id": "PRIORITY_TryR",
        "gene_name": "TryR",
        "sequence": "",
        "enzyme_class": "trypanothione_reductase",
        "pathway": "trypanothione",
        "evidence": "Ausente em humanos. Inibidores ja identificados in vitro.",
        "druggability_score": 0.96,
        "is_essential": True,
        "identity_score": 0.0,
        "active_site_diff": 1.0,
        "priority": True,
        "status": "priority_validated",
    },
    {
        "gene_id": "PRIORITY_SMT",
        "gene_name": "SMT",
        "sequence": "",
        "enzyme_class": "sterol_24C_methyltransferase",
        "pathway": "sterol_biosynthesis",
        "evidence": "Ausente em mamiferos. Base dos compostos azolicos antileishmania.",
        "druggability_score": 0.94,
        "is_essential": True,
        "identity_score": 0.0,
        "active_site_diff": 1.0,
        "priority": True,
        "status": "priority_validated",
    },
    {
        "gene_id": "PRIORITY_6PGDH",
        "gene_name": "6PGDH",
        "sequence": "",
        "enzyme_class": "6_phosphogluconate_dehydrogenase",
        "pathway": "pentose_phosphate",
        "evidence": "Menos de 35% de identidade com humano. Modelo 3D disponivel. Virtual screening realizado.",
        "druggability_score": 0.90,
        "is_essential": True,
        "identity_score": 0.35,
        "active_site_diff": 0.70,
        "priority": True,
        "status": "priority_validated",
    },
    {
        "gene_id": "PRIORITY_XPRT",
        "gene_name": "XPRT",
        "sequence": "",
        "enzyme_class": "xanthine_phosphoribosyltransferase",
        "pathway": "purine_salvage",
        "evidence": "Ausente em humanos. Incorpora analogos de purinas de forma unica.",
        "druggability_score": 0.88,
        "is_essential": False,
        "identity_score": 0.0,
        "active_site_diff": 1.0,
        "priority": True,
        "status": "priority_validated",
    },
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _build_alphafold_url(target: DrugTarget) -> str:
    """Build an AlphaFold structure URL for a drug target.

    Prefers the human homolog UniProt ID when available, then falls back
    to the target's own gene_id if it looks like a UniProt accession.

    Args:
        target: The drug target to generate a URL for.

    Returns:
        An AlphaFold URL string, or empty string if no suitable ID exists.
    """
    if target.human_homolog_id:
        return f"{ALPHAFOLD_BASE_URL}/{target.human_homolog_id}"

    # gene_id may itself be a UniProt-style accession (6-10 alphanumeric).
    gene_id = target.gene_id
    if gene_id and len(gene_id) <= 10 and gene_id.isalnum():
        return f"{ALPHAFOLD_BASE_URL}/{gene_id}"

    return ""


def _upsert_priority_targets() -> list[DrugTarget]:
    """Upsert all pre-validated priority targets into Supabase.

    Returns:
        List of ``DrugTarget`` instances for the priority targets.

    Raises:
        Exception: Logs and continues on individual upsert failures.
    """
    from core.db import upsert_drug_target

    targets: list[DrugTarget] = []

    for entry in PRIORITY_DRUG_TARGETS:
        target = DrugTarget.from_dict(entry)
        targets.append(target)

        try:
            upsert_drug_target(target)
            logger.info(
                "Upserted priority target %s (%s, score=%.2f).",
                target.gene_id,
                target.gene_name,
                target.druggability_score,
            )
        except Exception:
            logger.error(
                "Failed to upsert priority target %s — continuing.",
                target.gene_id,
            )

    logger.info("Upserted %d priority validated targets.", len(targets))
    return targets


def _score_computational_targets() -> list[DrugTarget]:
    """Load non-priority targets from the database, score, and update them.

    For each target whose ``priority`` flag is ``False``:
    1. Compute the druggability score via ``DrugTarget.compute_druggability_score()``.
    2. Attach an AlphaFold URL if the score exceeds the threshold.
    3. Persist updates back to Supabase.

    Returns:
        List of scored ``DrugTarget`` instances (non-priority only).
    """
    from core.db import get_all_drug_targets, update_drug_target

    all_targets = get_all_drug_targets()
    computational = [t for t in all_targets if not t.priority]

    approved_count = 0

    for target in computational:
        target.compute_druggability_score()

        if target.druggability_score > DRUGGABILITY_THRESHOLD:
            target.alphafold_url = _build_alphafold_url(target)
            target.status = "approved"
            approved_count += 1

        try:
            update_drug_target(
                target.gene_id,
                {
                    "druggability_score": target.druggability_score,
                    "alphafold_url": target.alphafold_url,
                    "status": target.status,
                },
            )
        except Exception:
            logger.error(
                "Failed to update target %s — continuing.", target.gene_id
            )

    logger.info(
        "Scored %d computational targets. %d with druggability > %.2f.",
        len(computational),
        approved_count,
        DRUGGABILITY_THRESHOLD,
    )

    return computational


def _write_csv(targets: list[DrugTarget], output_path: Path) -> None:
    """Write scored targets to a CSV file sorted by druggability_score descending.

    Args:
        targets: List of ``DrugTarget`` instances to export.
        output_path: Destination path for the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sorted_targets = sorted(
        targets, key=lambda t: t.druggability_score, reverse=True
    )

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        for target in sorted_targets:
            writer.writerow({
                "gene_id": target.gene_id,
                "gene_name": target.gene_name,
                "pathway": target.pathway,
                "identity_score": f"{target.identity_score:.4f}",
                "active_site_diff": f"{target.active_site_diff:.4f}",
                "is_essential": target.is_essential,
                "druggability_score": f"{target.druggability_score:.4f}",
                "alphafold_url": target.alphafold_url,
                "priority": target.priority,
            })

    logger.info("Wrote %d targets to %s.", len(sorted_targets), output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def score_druggability(force: bool = False, dry_run: bool = False) -> str:
    """Score all drug targets for druggability and export results.

    Upserts pre-validated priority targets, computes druggability scores
    for computational targets retrieved from the database, attaches
    AlphaFold structure URLs where applicable, and writes a ranked CSV.

    Args:
        force: Re-score and overwrite even if the output CSV already exists.
        dry_run: Only process priority targets, skip database queries for
                 computational targets. Writes CSV from priority targets only.

    Returns:
        Path to the output CSV file.

    Raises:
        RuntimeError: If no targets are available to score.
    """
    output_path = OUTPUT_FILE

    if output_path.exists() and not force:
        logger.info(
            "Output file already exists at %s. Use --force to re-score.",
            output_path,
        )
        return str(output_path)

    # Step 1: Always upsert priority validated targets.
    priority_targets = _upsert_priority_targets()

    if dry_run:
        logger.info(
            "[DRY RUN] Would score computational targets from database."
        )
        logger.info(
            "[DRY RUN] Writing CSV with %d priority targets only.",
            len(priority_targets),
        )
        _write_csv(priority_targets, output_path)
        return str(output_path)

    # Step 2: Score computational targets from the database.
    computational_targets = _score_computational_targets()

    # Step 3: Combine and write CSV.
    all_targets = priority_targets + computational_targets
    _write_csv(all_targets, output_path)

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Score druggability of L. infantum drug targets.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-score even if druggability_scores.csv already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Process priority targets only, skip database queries.",
    )
    args = parser.parse_args()

    result = score_druggability(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
