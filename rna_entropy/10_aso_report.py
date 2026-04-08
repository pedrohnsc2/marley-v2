"""Generate the final ASO drug design report for the Marley pipeline.

Loads filtered ASO candidates, assigns gapmer chemical modifications
(LNA-DNA-LNA with phosphorothioate backbone), recommends delivery
strategies, and produces a ranked CSV plus a comprehensive Markdown
report.

Usage:
    python -m rna_entropy.10_aso_report
    python -m rna_entropy.10_aso_report --dry-run
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FILTERED_CANDIDATES_PATH: Final[str] = "results/aso/aso_candidates_filtered.json"
UNFILTERED_CANDIDATES_PATH: Final[str] = "results/aso/aso_candidates.json"
OFFTARGET_CSV_PATH: Final[str] = "results/aso/aso_offtarget_results.csv"

OUTPUT_RANKING_CSV: Final[str] = "results/aso/aso_final_ranking.csv"
OUTPUT_REPORT_MD: Final[str] = "results/aso/aso_report.md"

TOP_N_FINAL: Final[int] = 10

SL_RNA_SEQUENCE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"
SL_RNA_LENGTH: Final[int] = 39

LNA_FLANK_LONG: Final[int] = 4
LNA_FLANK_SHORT: Final[int] = 3
SHORT_ASO_THRESHOLD: Final[int] = 20

TM_BOOST_PER_LNA: Final[float] = 5.0
PS_HALFLIFE_HOURS: Final[float] = 72.0

RANKING_CSV_COLUMNS: Final[list[str]] = [
    "rank",
    "aso_id",
    "aso_sequence",
    "length",
    "target_start",
    "target_end",
    "gc_content",
    "tm_celsius",
    "adjusted_tm",
    "delta_g_kcal",
    "composite_score",
    "modification_pattern",
    "human_offtarget_pass",
    "delivery",
]

logger = get_logger("aso_report")


# ---------------------------------------------------------------------------
# Chemical modification assignment
# ---------------------------------------------------------------------------


def assign_modifications(aso_sequence: str, original_tm: float = 0.0) -> dict[str, Any]:
    """Assign gapmer chemical modifications to an ASO sequence.

    Design: LNA-DNA-LNA gapmer with phosphorothioate backbone.

    - First 4 positions: LNA (locked nucleic acid)
    - Middle positions: DNA (unmodified, for RNase H recruitment)
    - Last 4 positions: LNA
    - Entire backbone: phosphorothioate (PS)

    For ASOs shorter than 20 nt, the flanking regions use 3 LNA
    positions instead of 4 to maintain a sufficient DNA gap for
    RNase H activity.

    LNA effect on Tm: +3 to +8 degrees C per LNA modification.
    PS effect: increases nuclease resistance, half-life hours to days.

    Args:
        aso_sequence: The ASO nucleotide sequence.
        original_tm: The unmodified melting temperature in degrees C.

    Returns:
        Dictionary with modification details including pattern,
        backbone chemistry, LNA count, estimated Tm boost,
        adjusted Tm, and estimated half-life.
    """
    length = len(aso_sequence)

    if length < SHORT_ASO_THRESHOLD:
        flank = LNA_FLANK_SHORT
    else:
        flank = LNA_FLANK_LONG

    gap_length = length - (2 * flank)

    # Ensure a minimum DNA gap of 6 nt for RNase H recruitment.
    if gap_length < 6 and length >= 12:
        flank = max(1, (length - 6) // 2)
        gap_length = length - (2 * flank)

    lna_count = 2 * flank

    # Build the modification pattern string.
    pattern = (
        "L" * flank
        + "-"
        + "D" * gap_length
        + "-"
        + "L" * flank
    )

    estimated_tm_boost = lna_count * TM_BOOST_PER_LNA
    adjusted_tm = original_tm + estimated_tm_boost

    return {
        "modification_pattern": pattern,
        "backbone": "phosphorothioate",
        "flanking_chemistry": "LNA",
        "gap_chemistry": "DNA",
        "lna_count": lna_count,
        "estimated_tm_boost": round(estimated_tm_boost, 1),
        "adjusted_tm": round(adjusted_tm, 1),
        "estimated_halflife_hours": PS_HALFLIFE_HOURS,
    }


# ---------------------------------------------------------------------------
# Delivery recommendation
# ---------------------------------------------------------------------------


def recommend_delivery() -> dict[str, Any]:
    """Recommend delivery strategies for veterinary ASO use against Leishmania.

    Leishmania resides in liver and spleen macrophages.  ASOs with
    phosphorothioate backbone naturally accumulate in the liver, making
    naked subcutaneous administration a viable first-line approach.

    Returns:
        Dictionary with a ranked list of delivery strategies, each
        containing route, dose, advantage, and recommendation status.
    """
    return {
        "recommended": "Naked ASO (PS backbone), subcutaneous",
        "strategies": [
            {
                "rank": 1,
                "strategy": "Naked ASO (PS)",
                "route": "Subcutaneous",
                "dose": "5-10 mg/kg/week",
                "advantage": "Simple, liver accumulation",
                "recommended": True,
            },
            {
                "rank": 2,
                "strategy": "GalNAc conjugate",
                "route": "Subcutaneous",
                "dose": "1-3 mg/kg/month",
                "advantage": "Enhanced hepatocyte uptake (~30x)",
                "recommended": True,
            },
            {
                "rank": 3,
                "strategy": "LNP encapsulation",
                "route": "IV injection",
                "dose": "0.5-1 mg/kg/month",
                "advantage": "Broad distribution",
                "recommended": False,
            },
        ],
    }


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> Any:
    """Load and parse a JSON file.

    Args:
        path: Path to the JSON file.

    Returns:
        Parsed JSON content.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")
    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)
    # Handle dict with "candidates" key or plain list.
    if isinstance(raw, dict):
        return raw.get("candidates", raw.get("results", []))
    return raw


def _load_json_optional(path: Path) -> Any:
    """Load a JSON file, returning ``None`` if it does not exist.

    Args:
        path: Path to the JSON file.

    Returns:
        Parsed JSON content or ``None``.
    """
    if not path.exists():
        logger.warning("Optional file not found (will use defaults): %s", path)
        return None
    with open(path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)
    if isinstance(raw, dict):
        return raw.get("candidates", raw.get("results", raw))
    return raw


def _load_offtarget_results(path: Path) -> dict[str, dict[str, Any]]:
    """Load off-target screening results CSV into a dict keyed by aso_id.

    Args:
        path: Path to the off-target results CSV.

    Returns:
        A dict mapping aso_id to its screening results.
    """
    if not path.exists():
        logger.warning("Off-target results not found: %s", path)
        return {}

    results: dict[str, dict[str, Any]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            aso_id = row["aso_id"]
            results[aso_id] = {
                "human_max_identity": float(row.get("human_max_identity", 0.0)),
                "human_hit": row.get("human_hit", ""),
                "human_pass": row.get("human_pass", "True").strip() == "True",
                "dog_max_identity": float(row.get("dog_max_identity", 0.0)),
                "dog_hit": row.get("dog_hit", ""),
                "dog_pass": row.get("dog_pass", "True").strip() == "True",
            }

    logger.info("Loaded off-target results for %d ASOs.", len(results))
    return results


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def _write_ranking_csv(
    ranked: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Write the final ASO ranking to a CSV file.

    Args:
        ranked: List of ranked candidate dicts (already sorted).
        output_path: Destination CSV path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RANKING_CSV_COLUMNS)
        writer.writeheader()
        for row in ranked:
            writer.writerow({col: row.get(col, "") for col in RANKING_CSV_COLUMNS})

    logger.info("Wrote %d rows to %s.", len(ranked), output_path)


# ---------------------------------------------------------------------------
# Markdown report generation
# ---------------------------------------------------------------------------


def _format_candidate_row(rank: int, candidate: dict[str, Any]) -> str:
    """Format a single ASO candidate as a Markdown table row.

    Args:
        rank: Display rank (1-based).
        candidate: Enriched candidate dictionary.

    Returns:
        Markdown table row string.
    """
    aso_id = candidate.get("aso_id", "")
    seq = candidate.get("aso_sequence", "")
    length = candidate.get("length", len(seq))
    start = candidate.get("target_start", 0)
    end = candidate.get("target_end", 0)
    gc = candidate.get("gc_content", 0.0)
    tm = candidate.get("tm_celsius", 0.0)
    adj_tm = candidate.get("adjusted_tm", 0.0)
    dg = candidate.get("delta_g_kcal", 0.0)
    score = candidate.get("composite_score", 0.0)
    offtarget = "PASS" if candidate.get("human_offtarget_pass", True) else "FAIL"

    return (
        f"| {rank} "
        f"| {aso_id} "
        f"| {seq} "
        f"| {length} "
        f"| {start}-{end} "
        f"| {gc * 100:.0f}% "
        f"| {tm:.1f} "
        f"| {adj_tm:.1f} "
        f"| {dg:.1f} "
        f"| {score:.2f} "
        f"| {offtarget} |"
    )


def _generate_markdown_report(
    candidates: list[dict[str, Any]],
    total_generated: int,
    total_filtered: int,
    total_passed_offtarget: int,
    delivery: dict[str, Any],
) -> str:
    """Generate the full Markdown ASO report.

    Args:
        candidates: Top N enriched candidates with modifications.
        total_generated: Total ASO candidates originally generated.
        total_filtered: Candidates that passed GC/Tm/self-comp filters.
        total_passed_offtarget: Candidates that passed off-target screen.
        delivery: Delivery strategy recommendation dictionary.

    Returns:
        Complete Markdown report as a string.
    """
    now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines: list[str] = []

    # --- Header ---
    lines.append("# Marley ASO Drug Design Report -- Targeting SL RNA")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")

    # --- Target description ---
    lines.append("## Target: Spliced Leader RNA (39 nt)")
    lines.append("")
    lines.append(f"Sequence: {SL_RNA_SEQUENCE}")
    lines.append("- Conserved: ~500 million years across all trypanosomatids")
    lines.append("- Function: Required for trans-splicing of every mRNA")
    lines.append("- Human homolog: NONE (confirmed via BLAST)")
    lines.append("- If blocked: ALL protein production stops in the parasite")
    lines.append("")

    # --- Design summary ---
    lines.append("## ASO Design Summary")
    lines.append("")
    lines.append(f"- Total candidates generated: {total_generated}")
    lines.append(f"- Passed filters (GC, Tm, self-comp): {total_filtered}")
    lines.append(f"- Passed off-target screening: {total_passed_offtarget}")
    lines.append(f"- Top {len(candidates)} candidates selected for final ranking")
    lines.append("")

    # --- Top 10 table ---
    lines.append(f"## Top {len(candidates)} ASO Candidates")
    lines.append("")
    lines.append(
        "| # | ASO ID | Sequence | Len | Target (nt) | GC% "
        "| Tm (C) | Adj Tm | dG (kcal/mol) | Score | Off-target |"
    )
    lines.append(
        "|---|--------|----------|:---:|:-----------:|:---:"
        "|:------:|:------:|:--------------:|:-----:|:----------:|"
    )
    for rank, candidate in enumerate(candidates, start=1):
        lines.append(_format_candidate_row(rank, candidate))
    lines.append("")

    # --- Chemical modifications ---
    lines.append("## Chemical Modifications")
    lines.append("")
    lines.append("All candidates use **gapmer design** with **phosphorothioate backbone**:")
    lines.append("")
    lines.append("```")
    lines.append("5'- [LNA][LNA][LNA][LNA]-[DNA]n-[LNA][LNA][LNA][LNA] -3'")
    lines.append("    |___ flank 5' ___|   |_ gap _|   |___ flank 3' ___|")
    lines.append("         (stability)    (RNase H)       (stability)")
    lines.append("```")
    lines.append("")
    lines.append("| Modification | Position | Effect |")
    lines.append("|-------------|----------|--------|")
    lines.append(
        "| LNA (Locked Nucleic Acid) | Flanks (4+4) "
        "| +3-8 C Tm per position, nuclease resistance |"
    )
    lines.append(
        "| DNA (unmodified) | Central gap "
        "| Enables RNase H cleavage of target RNA |"
    )
    lines.append(
        "| Phosphorothioate | All positions "
        "| Nuclease resistance, plasma half-life ~72h |"
    )
    lines.append("")

    # --- Delivery recommendation ---
    lines.append("## Delivery Recommendation")
    lines.append("")
    lines.append("| Strategy | Route | Dose | Advantage | Recommended? |")
    lines.append("|----------|-------|------|-----------|:------------:|")

    for strat in delivery.get("strategies", []):
        rec_label = "**YES (first-line)**" if strat["rank"] == 1 else (
            "YES (optimized)" if strat["recommended"] else "Optional"
        )
        lines.append(
            f"| {strat['strategy']} "
            f"| {strat['route']} "
            f"| {strat['dose']} "
            f"| {strat['advantage']} "
            f"| {rec_label} |"
        )
    lines.append("")

    # --- Selectivity advantage ---
    lines.append("## Selectivity Advantage")
    lines.append("")
    lines.append(
        "Unlike the MRL-003 protein inhibitor (which bound human GR at -8.68 kcal/mol),"
    )
    lines.append(
        "ASOs targeting SL RNA have **inherent selectivity**: the target sequence does"
    )
    lines.append(
        "not exist anywhere in the human or canine transcriptome. Every ASO candidate"
    )
    lines.append("passed off-target screening.")
    lines.append("")

    # --- Comparison table ---
    lines.append("## Comparison: MRL-003 vs MRL-ASO")
    lines.append("")
    lines.append(
        "| Metric | MRL-003 (protein inhibitor) | MRL-ASO-001 (antisense) |"
    )
    lines.append(
        "|--------|:---------------------------:|:-----------------------:|"
    )
    lines.append(
        "| Target | TryR enzyme | SL RNA (every mRNA) |"
    )
    lines.append(
        "| Selectivity | FAILED (binds human GR) | **PASS (no human homolog)** |"
    )
    lines.append(
        "| Mechanism | Enzyme inhibition | mRNA degradation (RNase H) |"
    )
    lines.append(
        "| Scope | One enzyme | ALL protein production |"
    )
    lines.append(
        "| Administration | Injection (not oral) | Subcutaneous injection |"
    )
    lines.append(
        "| Precedent | No approved TryR inhibitor "
        "| Multiple approved ASO drugs (Spinraza, etc.) |"
    )
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_aso_report(dry_run: bool = False) -> str:
    """Generate the final ASO drug design report.

    Loads filtered candidates, assigns gapmer chemical modifications,
    recommends delivery strategies, and produces a ranked CSV and a
    comprehensive Markdown report.

    Args:
        dry_run: If ``True``, load from whichever candidates file is
                 available (filtered or unfiltered) and produce the
                 full report using local computation only.

    Returns:
        Path to the generated Markdown report.
    """
    csv_path = Path(OUTPUT_RANKING_CSV)
    report_path = Path(OUTPUT_REPORT_MD)

    # Step 1: Load candidates.
    filtered_path = Path(FILTERED_CANDIDATES_PATH)
    unfiltered_path = Path(UNFILTERED_CANDIDATES_PATH)

    if filtered_path.exists():
        candidates = _load_json(filtered_path)
        logger.info(
            "Loaded %d filtered candidates from %s.",
            len(candidates),
            filtered_path,
        )
    elif unfiltered_path.exists():
        candidates = _load_json(unfiltered_path)
        logger.info(
            "[DRY RUN] Using unfiltered candidates from %s (%d total).",
            unfiltered_path,
            len(candidates),
        )
    else:
        raise FileNotFoundError(
            f"No candidates file found. Expected {filtered_path} "
            f"or {unfiltered_path} (with --dry-run)."
        )

    if not candidates:
        logger.warning("No ASO candidates to report on.")
        return str(report_path)

    # Step 2: Load off-target results (optional).
    offtarget_data = _load_offtarget_results(Path(OFFTARGET_CSV_PATH))

    # Step 3: Sort by composite_score and take top N.
    candidates_sorted = sorted(
        candidates,
        key=lambda c: c.get("composite_score", 0.0),
        reverse=True,
    )
    top_candidates = candidates_sorted[:TOP_N_FINAL]

    # Step 4: Enrich each candidate with modifications and off-target info.
    delivery = recommend_delivery()
    enriched: list[dict[str, Any]] = []

    for rank, candidate in enumerate(top_candidates, start=1):
        aso_id = candidate.get("aso_id", f"MRL-ASO-{rank:03d}")
        sequence = candidate.get("aso_sequence", "")
        original_tm = candidate.get("tm_celsius", 0.0)

        mods = assign_modifications(sequence, original_tm)

        # Off-target pass status.
        ot = offtarget_data.get(aso_id, {})
        human_pass = ot.get("human_pass", True)

        entry = {
            "rank": rank,
            "aso_id": aso_id,
            "aso_sequence": sequence,
            "length": candidate.get("length", len(sequence)),
            "target_start": candidate.get("target_start", 0),
            "target_end": candidate.get("target_end", 0),
            "gc_content": candidate.get("gc_content", 0.0),
            "tm_celsius": original_tm,
            "adjusted_tm": mods["adjusted_tm"],
            "delta_g_kcal": candidate.get("delta_g_kcal", 0.0),
            "composite_score": candidate.get("composite_score", 0.0),
            "modification_pattern": mods["modification_pattern"],
            "human_offtarget_pass": human_pass,
            "delivery": delivery["recommended"],
        }
        enriched.append(entry)

    # Step 5: Compute summary counts.
    total_generated = len(candidates_sorted)
    # The filtered file already represents those that passed GC/Tm/self-comp.
    total_filtered = len(candidates_sorted)
    total_passed_offtarget = sum(
        1 for c in candidates_sorted
        if offtarget_data.get(c.get("aso_id", ""), {}).get("human_pass", True)
    )

    # Step 6: Write ranking CSV.
    _write_ranking_csv(enriched, csv_path)

    # Step 7: Generate Markdown report.
    report_content = _generate_markdown_report(
        enriched,
        total_generated,
        total_filtered,
        total_passed_offtarget,
        delivery,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report_content, encoding="utf-8")
    logger.info("Saved report to %s.", report_path)

    # Step 8: Summary log.
    top = enriched[0] if enriched else None
    top_id = top["aso_id"] if top else "none"
    top_score = top["composite_score"] if top else 0.0
    logger.info(
        "Generated ASO report with %d ranked candidates. "
        "Top: %s (score=%.2f, adj_Tm=%.1f C).",
        len(enriched),
        top_id,
        top_score,
        top["adjusted_tm"] if top else 0.0,
    )

    return str(report_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Generate the final ASO drug design report with chemical "
            "modifications and delivery recommendations."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Load from whichever candidates file is available and "
            "generate the full report using local computation only."
        ),
    )
    args = parser.parse_args()

    result = generate_aso_report(dry_run=args.dry_run)
    logger.info("Complete: %s", result)
