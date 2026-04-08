"""Generate the drug target discovery report for the Marley v2 pipeline.

Fetches all scored drug targets from the database, ranks them by
druggability score, and produces two outputs:

1. ``results/drug_targets_top20.csv`` -- Top 20 targets in tabular form.
2. ``results/drug_targets_report.md`` -- Full Markdown report with
   experimentally validated and computational targets, plus next steps.

Usage:
    python -m drug_targets.05_report
    python -m drug_targets.05_report --dry-run
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Final

from core.db import get_all_drug_targets
from core.logger import get_logger
from core.models import DrugTarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESULTS_DIR: Final[Path] = Path("results")
CSV_OUTPUT: Final[str] = "results/drug_targets_top20.csv"
REPORT_OUTPUT: Final[str] = "results/drug_targets_report.md"
TOP_N: Final[int] = 20
COMPUTATIONAL_THRESHOLD: Final[float] = 0.60

CSV_COLUMNS: Final[list[str]] = [
    "rank",
    "gene_id",
    "gene_name",
    "enzyme_class",
    "pathway",
    "identity_score",
    "active_site_diff",
    "is_essential",
    "druggability_score",
    "alphafold_url",
    "priority",
    "evidence",
]

logger = get_logger("drug_targets_report")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _load_targets_from_csv(csv_path: Path) -> list[DrugTarget]:
    """Load drug targets from a CSV file as fallback when DB is unavailable.

    Args:
        csv_path: Path to the druggability_scores.csv file.

    Returns:
        List of ``DrugTarget`` instances parsed from the CSV.
    """
    targets: list[DrugTarget] = []

    with open(csv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            targets.append(DrugTarget(
                gene_id=row.get("gene_id", ""),
                gene_name=row.get("gene_name", ""),
                sequence="",
                enzyme_class=row.get("enzyme_class", ""),
                pathway=row.get("pathway", ""),
                identity_score=float(row.get("identity_score", 0.0)),
                active_site_diff=float(row.get("active_site_diff", 0.0)),
                is_essential=row.get("is_essential", "False") == "True",
                druggability_score=float(row.get("druggability_score", 0.0)),
                alphafold_url=row.get("alphafold_url", ""),
                evidence=row.get("evidence", ""),
                priority=row.get("priority", "False") == "True",
            ))

    logger.info("Loaded %d targets from CSV fallback.", len(targets))
    return targets


def _sort_targets(targets: list[DrugTarget]) -> list[DrugTarget]:
    """Sort drug targets by druggability_score in descending order.

    Args:
        targets: Unsorted list of drug targets.

    Returns:
        New list sorted by druggability_score (highest first).
    """
    return sorted(targets, key=lambda t: t.druggability_score, reverse=True)


def _write_csv(targets: list[DrugTarget], output_path: Path) -> None:
    """Write the top N targets to a CSV file.

    Args:
        targets: Sorted list of drug targets (already ranked).
        output_path: Path to the output CSV file.
    """
    top_targets = targets[:TOP_N]

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        for rank, target in enumerate(top_targets, start=1):
            writer.writerow({
                "rank": rank,
                "gene_id": target.gene_id,
                "gene_name": target.gene_name,
                "enzyme_class": target.enzyme_class,
                "pathway": target.pathway,
                "identity_score": f"{target.identity_score:.2f}",
                "active_site_diff": f"{target.active_site_diff:.2f}",
                "is_essential": target.is_essential,
                "druggability_score": f"{target.druggability_score:.2f}",
                "alphafold_url": target.alphafold_url,
                "priority": target.priority,
                "evidence": target.evidence,
            })

    logger.info("Wrote %d rows to %s.", len(top_targets), output_path)


def _format_alphafold_link(url: str) -> str:
    """Format an AlphaFold URL as a Markdown link.

    Args:
        url: The AlphaFold URL string (may be empty).

    Returns:
        Markdown link if URL is present, dash otherwise.
    """
    if url:
        return f"[link]({url})"
    return "-"


def _build_validated_section(targets: list[DrugTarget]) -> str:
    """Build the Markdown table for experimentally validated targets.

    Includes targets where ``priority=True``, sorted by druggability_score.

    Args:
        targets: Full sorted list of drug targets.

    Returns:
        Markdown string for Section 1.
    """
    validated = [t for t in targets if t.priority]

    lines: list[str] = []
    lines.append("## Section 1: Experimentally Validated Targets")
    lines.append("")
    lines.append("Targets with prior experimental evidence (priority=True).")
    lines.append("")

    if not validated:
        lines.append("_No experimentally validated targets found._")
        return "\n".join(lines)

    lines.append(
        "| # | Gene | Enzyme Class | Pathway | Score | AlphaFold | Evidence |"
    )
    lines.append(
        "|---|------|-------------|---------|-------|-----------|----------|"
    )

    for rank, target in enumerate(validated, start=1):
        af_link = _format_alphafold_link(target.alphafold_url)
        evidence = target.evidence if target.evidence else "-"
        lines.append(
            f"| {rank} "
            f"| {target.gene_name} "
            f"| {target.enzyme_class} "
            f"| {target.pathway} "
            f"| {target.druggability_score:.2f} "
            f"| {af_link} "
            f"| {evidence} |"
        )

    return "\n".join(lines)


def _build_computational_section(targets: list[DrugTarget]) -> str:
    """Build the Markdown table for computational targets.

    Includes targets where ``priority=False`` and
    ``druggability_score > COMPUTATIONAL_THRESHOLD``.

    Args:
        targets: Full sorted list of drug targets.

    Returns:
        Markdown string for Section 2.
    """
    computational = [
        t for t in targets
        if not t.priority and t.druggability_score > COMPUTATIONAL_THRESHOLD
    ]

    lines: list[str] = []
    lines.append("## Section 2: Computational Targets")
    lines.append("")
    lines.append(
        "New targets identified by the pipeline "
        f"(priority=False, score > {COMPUTATIONAL_THRESHOLD:.2f})."
    )
    lines.append("")

    if not computational:
        lines.append("_No computational targets above threshold._")
        return "\n".join(lines)

    lines.append(
        "| # | Gene | Enzyme Class | Pathway "
        "| Identity | Active Site Diff | Score | AlphaFold |"
    )
    lines.append(
        "|---|------|-------------|---------|"
        "----------|-----------------|-------|-----------|"
    )

    for rank, target in enumerate(computational, start=1):
        af_link = _format_alphafold_link(target.alphafold_url)
        lines.append(
            f"| {rank} "
            f"| {target.gene_name} "
            f"| {target.enzyme_class} "
            f"| {target.pathway} "
            f"| {target.identity_score:.2f} "
            f"| {target.active_site_diff:.2f} "
            f"| {target.druggability_score:.2f} "
            f"| {af_link} |"
        )

    return "\n".join(lines)


def _build_next_steps_section() -> str:
    """Build the static next-steps section of the report.

    Returns:
        Markdown string for Section 3.
    """
    return """## Section 3: Next Steps

### Molecular Docking
- Perform virtual screening against top 5 targets using AutoDock Vina
- Focus on TryS and TryR (absent in humans, highest selectivity potential)
- Use AlphaFold structures as receptor models

### In Vitro Validation
- Express recombinant enzymes for top 3 targets
- Enzyme inhibition assays (IC50 determination)
- Selectivity index: compare IC50 against human homolog

### Lead Optimization
- ADMET profiling of hit compounds
- Structure-activity relationship (SAR) studies
- Combination therapy assessment with existing antileishmanials

### Pathway-Based Strategy
- Purine salvage: dual targeting ADL + GMPS for synergistic effect
- Trypanothione: TryS + TryR combination for complete pathway blockade
- Sterol: SMT inhibitors as backbone, combine with miltefosine"""


def _build_report(targets: list[DrugTarget], generated_date: str) -> str:
    """Assemble the full Markdown report from all sections.

    Args:
        targets: Sorted list of all drug targets.
        generated_date: ISO-formatted date string for the report header.

    Returns:
        Complete Markdown report as a single string.
    """
    total = len(targets)
    validated_count = sum(1 for t in targets if t.priority)
    computational_count = sum(
        1 for t in targets
        if not t.priority and t.druggability_score > COMPUTATIONAL_THRESHOLD
    )
    pathways = sorted({t.pathway for t in targets if t.pathway})

    sections: list[str] = []

    # Header
    sections.append("# Marley v2 — Drug Target Discovery Report")
    sections.append("")
    sections.append(f"Generated: {generated_date}")
    sections.append("")

    # Summary
    sections.append("## Summary")
    sections.append(f"- Total targets analyzed: {total}")
    sections.append(f"- Priority validated targets: {validated_count}")
    sections.append(
        f"- Computational targets (score > {COMPUTATIONAL_THRESHOLD:.2f}): "
        f"{computational_count}"
    )
    sections.append(f"- Pathways covered: {', '.join(pathways)}")
    sections.append("")
    sections.append("---")
    sections.append("")

    # Section 1
    sections.append(_build_validated_section(targets))
    sections.append("")
    sections.append("---")
    sections.append("")

    # Section 2
    sections.append(_build_computational_section(targets))
    sections.append("")
    sections.append("---")
    sections.append("")

    # Section 3
    sections.append(_build_next_steps_section())
    sections.append("")

    return "\n".join(sections)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_drug_targets_report(dry_run: bool = False) -> str:
    """Generate the drug target discovery report (CSV + Markdown).

    Fetches all drug targets from the database, sorts by druggability
    score, and writes:

    - ``results/drug_targets_top20.csv`` with the top 20 targets.
    - ``results/drug_targets_report.md`` with the full analysis report.

    Args:
        dry_run: If ``True``, log what would be generated without writing
                 any files.

    Returns:
        Path to the generated Markdown report file.

    Raises:
        RuntimeError: If no drug targets are found in the database.
    """
    csv_path = Path(CSV_OUTPUT)
    report_path = Path(REPORT_OUTPUT)

    # Fetch and sort targets.
    logger.info("Fetching drug targets from database...")
    all_targets: list[DrugTarget] = []

    try:
        all_targets = get_all_drug_targets()
    except Exception:
        logger.warning("Could not fetch from database — falling back to CSV.")

    # Fallback: read from the druggability CSV if DB is unavailable.
    if not all_targets:
        druggability_csv = Path("results/druggability_scores.csv")
        if druggability_csv.exists():
            logger.info("Loading targets from %s.", druggability_csv)
            all_targets = _load_targets_from_csv(druggability_csv)

    if not all_targets:
        msg = "No drug targets found. Run earlier pipeline steps first."
        logger.error(msg)
        raise RuntimeError(msg)

    targets = _sort_targets(all_targets)
    logger.info(
        "Loaded %d drug target(s). Top score: %.2f (%s).",
        len(targets),
        targets[0].druggability_score,
        targets[0].gene_name,
    )

    if dry_run:
        logger.info("[DRY RUN] Would write CSV with top %d targets to %s.", TOP_N, csv_path)
        logger.info("[DRY RUN] Would write Markdown report to %s.", report_path)
        validated_count = sum(1 for t in targets if t.priority)
        computational_count = sum(
            1 for t in targets
            if not t.priority and t.druggability_score > COMPUTATIONAL_THRESHOLD
        )
        logger.info(
            "[DRY RUN] Summary: %d total, %d validated, %d computational.",
            len(targets),
            validated_count,
            computational_count,
        )
        return str(report_path)

    # Ensure output directory exists.
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Write CSV.
    _write_csv(targets, csv_path)

    # Write Markdown report.
    generated_date = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")
    report_content = _build_report(targets, generated_date)

    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_content)

    logger.info("Wrote report to %s.", report_path)

    return str(report_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate the drug target discovery report (CSV + Markdown).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would be generated without writing files.",
    )
    args = parser.parse_args()

    result = generate_drug_targets_report(dry_run=args.dry_run)
    logger.info("Complete: %s", result)
