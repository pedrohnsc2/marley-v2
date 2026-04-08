"""Generate a comprehensive docking report from ADMET-filtered results.

Reads the ADMET-filtered docking results, produces a ranked top-hits CSV,
a full Markdown analysis report, and PyMOL visualization scripts for each
top hit whose docked pose PDBQT file exists on disk.

Usage:
    python -m drug_targets.10_docking_report
    python -m drug_targets.10_docking_report --dry-run
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESULTS_DIR: Final[Path] = Path("results")
ADMET_INPUT: Final[Path] = RESULTS_DIR / "admet_filtered.csv"
TOP_HITS_OUTPUT: Final[Path] = RESULTS_DIR / "docking_top_hits.csv"
REPORT_OUTPUT: Final[Path] = RESULTS_DIR / "docking_report.md"
PYMOL_DIR: Final[Path] = RESULTS_DIR / "docking_pymol"
STRUCTURES_DIR: Final[Path] = Path("data/structures")
DOCKING_DIR: Final[Path] = Path("data/docking")

TOP_N: Final[int] = 10
AFFINITY_THRESHOLD: Final[float] = -6.0

TOP_HITS_COLUMNS: Final[list[str]] = [
    "rank",
    "target_gene_name",
    "compound_id",
    "compound_name",
    "smiles",
    "binding_affinity",
    "lipinski_violations",
    "admet_score",
    "is_approved_drug",
    "composite_score",
]

# Map of target gene short names to descriptive names for the report.
TARGET_DESCRIPTIONS: Final[dict[str, str]] = {
    "TryS": "trypanothione synthetase",
    "TryR": "trypanothione reductase",
    "ADL": "adenylosuccinate lyase",
    "SMT": "sterol 24-C-methyltransferase",
    "6PGDH": "6-phosphogluconate dehydrogenase",
    "GMPS": "GMP synthase",
    "HGPRT": "hypoxanthine-guanine phosphoribosyltransferase",
    "XPRT": "xanthine phosphoribosyltransferase",
    "DHFR-TS": "dihydrofolate reductase-thymidylate synthase",
    "PTR1": "pteridine reductase 1",
}

logger = get_logger("docking_report")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_admet_results(input_path: Path) -> list[dict[str, str]]:
    """Load ADMET-filtered docking results from CSV.

    Args:
        input_path: Path to the ``admet_filtered.csv`` file.

    Returns:
        List of row dictionaries as read from the CSV.

    Raises:
        FileNotFoundError: If the input CSV does not exist.
        RuntimeError: If the CSV is empty or has no data rows.
    """
    if not input_path.exists():
        msg = (
            f"ADMET-filtered results not found at {input_path}. "
            "Run module 09 (ADMET filtering) first."
        )
        raise FileNotFoundError(msg)

    rows: list[dict[str, str]] = []
    with open(input_path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)

    if not rows:
        msg = f"ADMET-filtered CSV at {input_path} contains no data rows."
        raise RuntimeError(msg)

    logger.info("Loaded %d ADMET-filtered results from %s.", len(rows), input_path)
    return rows


def _parse_float(value: str, default: float = 0.0) -> float:
    """Safely parse a string to float with a fallback default.

    Args:
        value: String representation of a number.
        default: Value to return if parsing fails.

    Returns:
        Parsed float or the default.
    """
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def _parse_bool(value: str) -> bool:
    """Parse a string boolean value.

    Args:
        value: String like "True", "true", "1", "yes".

    Returns:
        Boolean interpretation of the string.
    """
    return value.strip().lower() in {"true", "1", "yes"}


# ---------------------------------------------------------------------------
# Top hits CSV
# ---------------------------------------------------------------------------


def _rank_top_hits(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Sort rows by composite_score descending and return the top N.

    Args:
        rows: All ADMET-filtered result rows.

    Returns:
        Top N rows sorted by composite_score (highest first), each
        augmented with a ``rank`` field.
    """
    sorted_rows = sorted(
        rows,
        key=lambda r: _parse_float(r.get("composite_score", "0")),
        reverse=True,
    )

    top = sorted_rows[:TOP_N]
    for idx, row in enumerate(top, start=1):
        row["rank"] = str(idx)

    return top


def _write_top_hits_csv(top_hits: list[dict[str, str]], output_path: Path) -> None:
    """Write the ranked top hits to a CSV file.

    Args:
        top_hits: Ranked list of top hit rows (must include ``rank``).
        output_path: Destination CSV file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=TOP_HITS_COLUMNS,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in top_hits:
            writer.writerow(row)

    logger.info("Wrote %d top hits to %s.", len(top_hits), output_path)


# ---------------------------------------------------------------------------
# Markdown report builder
# ---------------------------------------------------------------------------


def _build_summary_section(
    rows: list[dict[str, str]],
    top_hits: list[dict[str, str]],
    generated_date: str,
) -> str:
    """Build the report header and summary section.

    Args:
        rows: All ADMET-filtered result rows.
        top_hits: Ranked top N hits.
        generated_date: ISO-formatted date for the header.

    Returns:
        Markdown string for the header and summary.
    """
    total_pairs = len(rows)
    unique_compounds = len({r.get("compound_id", "") for r in rows})
    approved_in_top = sum(
        1 for r in top_hits if _parse_bool(r.get("is_approved_drug", "False"))
    )
    targets_with_hits = sorted({r.get("target_gene_name", "") for r in top_hits})

    lines: list[str] = [
        "# Marley v3 -- Molecular Docking Report",
        "",
        f"Generated: {generated_date}",
        "",
        "## Summary",
        f"- Total docking pairs evaluated: {total_pairs}",
        f"- Compounds passing ADMET filter: {unique_compounds}",
        f"- Approved drugs in top 10: {approved_in_top}",
        f"- Targets with hits: {', '.join(targets_with_hits)}",
        "",
        "---",
        "",
    ]
    return "\n".join(lines)


def _build_top_hits_table(top_hits: list[dict[str, str]]) -> str:
    """Build the top 10 docking hits Markdown table.

    Args:
        top_hits: Ranked top N hit rows.

    Returns:
        Markdown string for the top hits section.
    """
    lines: list[str] = [
        "## Top 10 Docking Hits",
        "",
        "| # | Target | Compound | Name | Affinity (kcal/mol) "
        "| Lipinski | ADMET | Approved | Score |",
        "|---|--------|----------|------|---------------------"
        "|----------|-------|----------|-------|",
    ]

    for row in top_hits:
        rank = row.get("rank", "")
        target = row.get("target_gene_name", "")
        compound_id = row.get("compound_id", "")
        name = row.get("compound_name", "")
        affinity = _parse_float(row.get("binding_affinity", "0"))
        lipinski = row.get("lipinski_violations", "0")
        admet = _parse_float(row.get("admet_score", "0"))
        approved = "Yes" if _parse_bool(row.get("is_approved_drug", "False")) else "No"
        score = _parse_float(row.get("composite_score", "0"))

        lines.append(
            f"| {rank} | {target} | {compound_id} | {name} "
            f"| {affinity:.1f} | {lipinski} | {admet:.2f} "
            f"| {approved} | {score:.2f} |"
        )

    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _build_repurposing_section(rows: list[dict[str, str]]) -> str:
    """Build the drug repurposing candidates section.

    Includes approved drugs with binding affinity below the threshold
    against priority targets.

    Args:
        rows: All ADMET-filtered result rows.

    Returns:
        Markdown string for the repurposing section.
    """
    candidates = [
        r for r in rows
        if _parse_bool(r.get("is_approved_drug", "False"))
        and _parse_float(r.get("binding_affinity", "0")) < AFFINITY_THRESHOLD
    ]

    # Sort by binding affinity (most negative first).
    candidates.sort(key=lambda r: _parse_float(r.get("binding_affinity", "0")))

    lines: list[str] = [
        "## Drug Repurposing Candidates",
        "",
        f"Approved drugs that showed binding affinity < {AFFINITY_THRESHOLD} "
        "kcal/mol against priority targets.",
        "",
    ]

    if not candidates:
        lines.append(
            "_No approved drugs met the affinity threshold._"
        )
        lines.extend(["", "---", ""])
        return "\n".join(lines)

    lines.extend([
        "| Drug | Target | Affinity | Current Indication | Potential |",
        "|------|--------|----------|--------------------|-----------|",
    ])

    for row in candidates:
        name = row.get("compound_name", row.get("compound_id", ""))
        target = row.get("target_gene_name", "")
        affinity = _parse_float(row.get("binding_affinity", "0"))
        desc = TARGET_DESCRIPTIONS.get(target, target)

        lines.append(
            f"| {name} | {target} | {affinity:.1f} "
            f"| Approved antileishmanial/antifungal | {desc} inhibitor |"
        )

    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _build_per_target_section(rows: list[dict[str, str]]) -> str:
    """Build per-target analysis subsections.

    Args:
        rows: All ADMET-filtered result rows.

    Returns:
        Markdown string with one subsection per target.
    """
    # Group rows by target.
    targets: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        gene = row.get("target_gene_name", "unknown")
        targets.setdefault(gene, []).append(row)

    lines: list[str] = ["## Per-Target Analysis", ""]

    for gene_name in sorted(targets):
        target_rows = targets[gene_name]
        desc = TARGET_DESCRIPTIONS.get(gene_name, gene_name)

        # Find best hit by binding affinity (most negative).
        best_row = min(
            target_rows,
            key=lambda r: _parse_float(r.get("binding_affinity", "0")),
        )
        best_compound = best_row.get(
            "compound_name",
            best_row.get("compound_id", "unknown"),
        )
        best_affinity = _parse_float(best_row.get("binding_affinity", "0"))

        compounds_tested = len(target_rows)
        hits_below_threshold = sum(
            1 for r in target_rows
            if _parse_float(r.get("binding_affinity", "0")) < AFFINITY_THRESHOLD
        )

        lines.extend([
            f"### {gene_name} ({desc})",
            f"- Best hit: {best_compound} at {best_affinity:.1f} kcal/mol",
            f"- Compounds tested: {compounds_tested}",
            f"- Hits with affinity < {AFFINITY_THRESHOLD}: {hits_below_threshold}",
            "",
        ])

    lines.extend(["---", ""])
    return "\n".join(lines)


def _build_validation_section() -> str:
    """Build the validation notes section with cross-references and limitations.

    Returns:
        Markdown string for the validation section.
    """
    return """## Validation Notes

### Cross-reference with Literature
- TryR: Known inhibitors should rank in top 10 if pipeline is calibrated
- SMT: Azole antifungals (ketoconazole, itraconazole) expected to bind
- Allopurinol: Known purine analog, should show affinity for HGPRT/XPRT

### Limitations
- Docking scores are estimates -- in vitro validation required
- AlphaFold structures may have low confidence in binding site regions
- Grid box positioning affects results significantly
- Binding affinity does not account for selectivity over human homolog

---
"""


def _build_next_steps_section() -> str:
    """Build the next steps section.

    Returns:
        Markdown string for the next steps section.
    """
    return """## Next Steps
1. Validate top 3 hits with molecular dynamics simulation (100 ns)
2. Test approved drug candidates in L. infantum promastigote assay
3. For novel compounds: assess synthetic accessibility
4. Compare binding poses between parasite and human homolog enzymes
"""


def _build_report(
    rows: list[dict[str, str]],
    top_hits: list[dict[str, str]],
    generated_date: str,
) -> str:
    """Assemble the full Markdown report from all sections.

    Args:
        rows: All ADMET-filtered result rows.
        top_hits: Ranked top N hit rows.
        generated_date: ISO date string for the header.

    Returns:
        Complete Markdown report as a single string.
    """
    sections: list[str] = [
        _build_summary_section(rows, top_hits, generated_date),
        _build_top_hits_table(top_hits),
        _build_repurposing_section(rows),
        _build_per_target_section(rows),
        _build_validation_section(),
        _build_next_steps_section(),
    ]
    return "\n".join(sections)


# ---------------------------------------------------------------------------
# PyMOL visualization scripts
# ---------------------------------------------------------------------------


def _generate_pymol_scripts(top_hits: list[dict[str, str]]) -> int:
    """Generate PyMOL .pml scripts for each top hit with an existing docked pose.

    For each hit, checks whether the receptor PDB and docked PDBQT exist
    on disk before writing the script.

    Args:
        top_hits: Ranked top N hit rows.

    Returns:
        Number of PyMOL scripts written.
    """
    PYMOL_DIR.mkdir(parents=True, exist_ok=True)
    written = 0

    for row in top_hits:
        target = row.get("target_gene_name", "")
        compound_id = row.get("compound_id", "")

        if not target or not compound_id:
            continue

        receptor_pdb = STRUCTURES_DIR / f"{target}_receptor.pdb"
        docked_pdbqt = DOCKING_DIR / target / f"{compound_id}_out.pdbqt"

        if not docked_pdbqt.exists():
            logger.debug(
                "Docked pose not found for %s/%s, skipping PyMOL script.",
                target,
                compound_id,
            )
            continue

        script_path = PYMOL_DIR / f"{target}_{compound_id}.pml"
        script_content = (
            f"load {receptor_pdb}\n"
            f"load {docked_pdbqt}\n"
            "hide everything\n"
            "show cartoon, *_receptor\n"
            "show sticks, *_out\n"
            "color cyan, *_receptor\n"
            "color yellow, *_out\n"
            "zoom *_out\n"
        )

        script_path.write_text(script_content, encoding="utf-8")
        written += 1
        logger.debug("Wrote PyMOL script: %s", script_path)

    if written:
        logger.info("Generated %d PyMOL visualization scripts in %s.", written, PYMOL_DIR)

    return written


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_docking_report(dry_run: bool = False) -> str:
    """Generate the docking report from ADMET-filtered results.

    Produces three outputs:

    1. ``results/docking_top_hits.csv`` -- top 10 compound-target pairs
       ranked by composite score.
    2. ``results/docking_report.md`` -- full Markdown analysis report with
       summary, top hits table, repurposing candidates, per-target analysis,
       validation notes, and next steps.
    3. ``results/docking_pymol/*.pml`` -- PyMOL visualization scripts for
       each top hit whose docked pose PDBQT exists on disk.

    Args:
        dry_run: If ``True``, log what would be generated without writing
                 any files to disk.

    Returns:
        Path to the generated Markdown report file
        (``results/docking_report.md``).

    Raises:
        FileNotFoundError: If the ADMET-filtered CSV does not exist.
        RuntimeError: If the CSV is empty or contains no data rows.
    """
    report_path = REPORT_OUTPUT

    # Step 1: Load ADMET-filtered results.
    rows = _load_admet_results(ADMET_INPUT)

    # Step 2: Rank top hits.
    top_hits = _rank_top_hits(rows)
    targets_in_top = sorted({r.get("target_gene_name", "") for r in top_hits})

    if dry_run:
        logger.info("[DRY RUN] Would write top %d hits to %s.", len(top_hits), TOP_HITS_OUTPUT)
        logger.info("[DRY RUN] Would write Markdown report to %s.", report_path)
        logger.info("[DRY RUN] Would generate PyMOL scripts in %s.", PYMOL_DIR)
        logger.info(
            "[DRY RUN] Generated docking report with %d top hits across %d targets.",
            len(top_hits),
            len(targets_in_top),
        )
        return str(report_path)

    # Step 3: Write top hits CSV.
    _write_top_hits_csv(top_hits, TOP_HITS_OUTPUT)

    # Step 4: Generate Markdown report.
    generated_date = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")
    report_content = _build_report(rows, top_hits, generated_date)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_content)

    logger.info("Wrote docking report to %s.", report_path)

    # Step 5: Generate PyMOL visualization scripts.
    _generate_pymol_scripts(top_hits)

    # Step 6: Summary log.
    logger.info(
        "Generated docking report with %d top hits across %d targets.",
        len(top_hits),
        len(targets_in_top),
    )

    return str(report_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Generate a comprehensive molecular docking report from "
            "ADMET-filtered results."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would be generated without writing files.",
    )
    args = parser.parse_args()

    result_path = generate_docking_report(dry_run=args.dry_run)
    logger.info("Complete: %s", result_path)
