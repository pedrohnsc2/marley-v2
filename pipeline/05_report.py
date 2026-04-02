"""Generate the final ranked report of vaccine candidates with scores and annotations.

Fetches all scored candidates from Supabase, ranks them by ``final_score``,
and produces two output artefacts:

1. **CSV** -- top candidates with all scores (``results/top_candidates.csv``).
2. **Markdown report** -- summary statistics, ranked table, methodology,
   and next steps (``results/report.md``).

Usage:
    python -m pipeline.05_report
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import pandas as pd

from core.db import get_all_candidates
from core.logger import get_logger
from core.models import Candidate, STATUS_APPROVED

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTPUT_CSV: str = "results/top_candidates.csv"
OUTPUT_REPORT: str = "results/report.md"
TOP_N: int = 20
TOP_DISPLAY: int = 10

logger = get_logger("report")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ensure_output_dir() -> None:
    """Create the results directory if it does not exist."""
    Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)


def _build_summary_stats(
    all_candidates: list[Candidate],
    scored: list[Candidate],
) -> dict[str, int]:
    """Compute summary statistics across the candidate pool.

    Args:
        all_candidates: Every candidate row from the database.
        scored: The subset that passed all filters and has a final score.

    Returns:
        Dictionary with keys: ``total``, ``surface_filter``,
        ``conservation_filter``, ``final_candidates``.
    """
    total = len(all_candidates)
    surface_filter = sum(
        1 for c in all_candidates if c.has_signal_peptide
    )
    conservation_filter = sum(
        1 for c in all_candidates
        if c.conservation_score > 0.0 and c.status == STATUS_APPROVED
    )
    final_candidates = len(scored)

    return {
        "total": total,
        "surface_filter": surface_filter,
        "conservation_filter": conservation_filter,
        "final_candidates": final_candidates,
    }


def _generate_csv(candidates: list[Candidate]) -> str:
    """Write the top candidates to a CSV file.

    Args:
        candidates: Sorted list of candidates (descending by final_score).

    Returns:
        Path to the written CSV file.
    """
    top = candidates[:TOP_N]

    df = pd.DataFrame(
        [
            {
                "gene_id": c.gene_id,
                "gene_name": c.gene_name,
                "has_signal_peptide": c.has_signal_peptide,
                "conservation_score": round(c.conservation_score, 4),
                "immunogenicity_score": round(c.immunogenicity_score, 4),
                "final_score": round(c.final_score, 4),
                "status": c.status,
            }
            for c in top
        ]
    )

    df.to_csv(OUTPUT_CSV, index=False, quoting=csv.QUOTE_NONNUMERIC)
    logger.info("CSV written to %s (%d rows).", OUTPUT_CSV, len(df))
    return OUTPUT_CSV


def _generate_markdown(
    candidates: list[Candidate],
    stats: dict[str, int],
) -> str:
    """Write a Markdown report summarising the pipeline results.

    Args:
        candidates: Sorted list of candidates (descending by final_score).
        stats: Summary statistics from :func:`_build_summary_stats`.

    Returns:
        Path to the written Markdown file.
    """
    today = datetime.now().strftime("%Y-%m-%d")
    top = candidates[:TOP_N]

    lines: list[str] = []

    # Title
    lines.append(f"# Marley -- Vaccine Candidate Report")
    lines.append("")
    lines.append(f"**Generated:** {today}")
    lines.append("")

    # Summary statistics
    lines.append("## Summary Statistics")
    lines.append("")
    lines.append(f"| Metric | Count |")
    lines.append(f"|---|---|")
    lines.append(f"| Total proteins analysed | {stats['total']} |")
    lines.append(
        f"| Passed surface filter (signal peptide) | {stats['surface_filter']} |"
    )
    lines.append(
        f"| Passed conservation filter | {stats['conservation_filter']} |"
    )
    lines.append(f"| Final scored candidates | {stats['final_candidates']} |")
    lines.append("")

    # Top candidates table
    lines.append(f"## Top {TOP_N} Candidates")
    lines.append("")
    lines.append(
        "| Rank | Gene ID | Gene Name | Conservation | Immunogenicity | Final Score |"
    )
    lines.append("|---|---|---|---|---|---|")

    for rank, c in enumerate(top, start=1):
        lines.append(
            f"| {rank} | {c.gene_id} | {c.gene_name} "
            f"| {c.conservation_score:.4f} | {c.immunogenicity_score:.4f} "
            f"| {c.final_score:.4f} |"
        )

    lines.append("")

    # Methodology
    lines.append("## Methodology")
    lines.append("")
    lines.append(
        "This report was generated by the **Marley** reverse vaccinology pipeline, "
        "which applies a multi-stage computational screen to identify promising "
        "vaccine candidates against *Leishmania infantum* (canine visceral leishmaniasis)."
    )
    lines.append("")
    lines.append("1. **Genome fetch** -- The complete *L. infantum* JPCM5 proteome was "
                 "downloaded from TriTrypDB.")
    lines.append("2. **Surface filter** -- Proteins were screened for Sec/SPI signal "
                 "peptides using SignalP 6.0, retaining surface-exposed and secreted proteins.")
    lines.append("3. **Conservation analysis** -- BLAST searches against related "
                 "*Leishmania* strains determined sequence conservation scores.")
    lines.append("4. **Immunogenicity scoring** -- MHC-I binding predictions (IEDB) for "
                 "canine DLA alleles estimated each candidate's immunogenic potential.")
    lines.append("5. **Final ranking** -- A weighted composite score "
                 "(40% conservation + 60% immunogenicity) produced the final ranking.")
    lines.append("")

    # Next steps
    lines.append("## Next Steps")
    lines.append("")
    lines.append("- **Experimental validation**: Express and purify top candidates for "
                 "in-vitro immunoassays (ELISPOT, cytokine profiling).")
    lines.append("- **B-cell epitope mapping**: Complement MHC-I predictions with B-cell "
                 "linear and conformational epitope analysis.")
    lines.append("- **In-vivo testing**: Evaluate protective efficacy in a murine or "
                 "canine challenge model.")
    lines.append("- **Adjuvant selection**: Screen adjuvant formulations to maximise "
                 "Th1-biased immune response.")
    lines.append("- **Multi-valent design**: Consider combining top candidates into a "
                 "multi-epitope or chimeric construct.")
    lines.append("")

    report_text = "\n".join(lines)
    Path(OUTPUT_REPORT).write_text(report_text, encoding="utf-8")
    logger.info("Markdown report written to %s.", OUTPUT_REPORT)
    return OUTPUT_REPORT


def _print_top_candidates(candidates: list[Candidate]) -> None:
    """Print the top candidates to the terminal with formatted output.

    Args:
        candidates: Sorted list of candidates (descending by final_score).
    """
    top = candidates[:TOP_DISPLAY]

    if not top:
        print("\n  No scored candidates to display.\n")
        return

    header = (
        f"{'Rank':<6} {'Gene ID':<25} {'Gene Name':<30} "
        f"{'Conserv.':<10} {'Immuno.':<10} {'Final':<10}"
    )
    separator = "-" * len(header)

    print(f"\n  Top {TOP_DISPLAY} Vaccine Candidates")
    print(f"  {separator}")
    print(f"  {header}")
    print(f"  {separator}")

    for rank, c in enumerate(top, start=1):
        row = (
            f"{rank:<6} {c.gene_id:<25} {c.gene_name:<30} "
            f"{c.conservation_score:<10.4f} {c.immunogenicity_score:<10.4f} "
            f"{c.final_score:<10.4f}"
        )
        print(f"  {row}")

    print(f"  {separator}\n")


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------


def generate_report() -> tuple[str, str]:
    """Generate the final pipeline report and CSV output.

    Fetches all candidates from Supabase, filters to those that are approved
    and have been scored, sorts by ``final_score`` descending, and produces:

    - ``results/top_candidates.csv`` -- top 20 candidates with all scores.
    - ``results/report.md`` -- full Markdown report with statistics,
      candidate table, methodology, and next steps.

    The top 10 candidates are also printed to the terminal.

    Returns:
        A tuple of ``(csv_path, report_path)``.
    """
    logger.info("Generating final report ...")

    _ensure_output_dir()

    # Fetch all candidates from Supabase.
    all_candidates = get_all_candidates()
    logger.info("Fetched %d total candidate(s) from database.", len(all_candidates))

    # Filter to approved/scored candidates only.
    scored = [
        c for c in all_candidates
        if c.status == STATUS_APPROVED and c.final_score > 0.0
    ]
    logger.info(
        "%d candidate(s) are approved with a final score.", len(scored)
    )

    # Sort by final_score descending.
    scored.sort(key=lambda c: c.final_score, reverse=True)

    # Summary statistics.
    stats = _build_summary_stats(all_candidates, scored)

    # Generate artefacts.
    csv_path = _generate_csv(scored)
    report_path = _generate_markdown(scored, stats)

    # Terminal output.
    _print_top_candidates(scored)

    logger.info("Report generation complete.")
    return csv_path, report_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    csv_out, report_out = generate_report()
    print(f"CSV:    {csv_out}")
    print(f"Report: {report_out}")
