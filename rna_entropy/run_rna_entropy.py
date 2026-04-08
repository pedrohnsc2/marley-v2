"""Entry point for running the RNA entropy analysis pipeline end-to-end.

Executes all seven pipeline stages sequentially:

1. Fetch transcriptome (UniProt mRNA download)
2. Codon usage analysis (RSCU computation)
3. Shannon entropy (per-transcript information content)
4. SL RNA analysis (spliced leader detection)
5. Human comparison (entropy delta calculation)
6. Structure prediction (MFE via ViennaRNA or GC fallback)
7. Report generation (CSV + Markdown summary)

Usage:
    python rna_entropy/run_rna_entropy.py                  # full run
    python rna_entropy/run_rna_entropy.py --skip-fetch     # skip download
    python rna_entropy/run_rna_entropy.py --dry-run        # test mode
    python rna_entropy/run_rna_entropy.py --force          # force re-run
    python rna_entropy/run_rna_entropy.py --priority-only  # report only
    python rna_entropy/run_rna_entropy.py --top-n 30       # custom top-N
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from typing import Any, Callable

from core.logger import get_logger

logger = get_logger("rna_pipeline")

# ---------------------------------------------------------------------------
# Module imports via importlib (numbered module names)
# ---------------------------------------------------------------------------


def _load_pipeline_modules() -> dict[str, Any]:
    """Import numbered RNA entropy pipeline modules using importlib.

    Returns:
        Dictionary mapping short names to imported module objects.
    """
    modules = {
        "fetch": importlib.import_module("rna_entropy.01_fetch_transcriptome"),
        "codon": importlib.import_module("rna_entropy.02_codon_usage"),
        "entropy": importlib.import_module("rna_entropy.03_shannon_entropy"),
        "sl_rna": importlib.import_module("rna_entropy.04_sl_rna_analysis"),
        "human": importlib.import_module("rna_entropy.05_human_comparison"),
        "structure": importlib.import_module("rna_entropy.06_structure_prediction"),
        "report": importlib.import_module("rna_entropy.07_report"),
    }
    return modules


# ---------------------------------------------------------------------------
# Stage definitions
# ---------------------------------------------------------------------------


class StageResult:
    """Container for the outcome of a single pipeline stage."""

    def __init__(self, name: str) -> None:
        self.name: str = name
        self.status: str = "pending"
        self.duration_seconds: float = 0.0
        self.output: str = ""
        self.error: str = ""


def _confirm_continue(stage_name: str) -> bool:
    """Prompt the user for confirmation before running a stage.

    Args:
        stage_name: Human-readable name of the upcoming stage.

    Returns:
        ``True`` if the user wants to continue, ``False`` to quit.
    """
    try:
        response = input(
            f"\n  [Stage: {stage_name}] "
            f"Press Enter to continue or 'q' to quit: "
        )
    except (EOFError, KeyboardInterrupt):
        return False

    return response.strip().lower() != "q"


def _run_stage(
    name: str,
    func: Callable[..., Any],
    kwargs: dict[str, Any],
    dry_run: bool = False,
    skip_confirm: bool = False,
) -> StageResult:
    """Execute a single pipeline stage with timing and error handling.

    Args:
        name: Human-readable stage name.
        func: The callable to invoke for this stage.
        kwargs: Keyword arguments passed to *func*.
        dry_run: If ``True``, log what would be done without calling *func*.
        skip_confirm: If ``True``, do not prompt for user confirmation.

    Returns:
        A :class:`StageResult` with timing and status information.
    """
    result = StageResult(name)

    # Confirmation gate.
    if not skip_confirm:
        if not _confirm_continue(name):
            result.status = "skipped (user quit)"
            logger.info("User chose to quit before stage: %s", name)
            return result

    # Dry-run mode.
    if dry_run:
        logger.info("[DRY RUN] Would execute: %s", name)
        result.status = "dry-run"
        result.output = "(dry run -- no action taken)"
        return result

    # Execute the stage.
    logger.info("Starting stage: %s", name)
    start = time.time()

    try:
        output = func(**kwargs)
        elapsed = time.time() - start
        result.duration_seconds = elapsed
        result.status = "success"
        result.output = str(output) if output is not None else ""
        logger.info(
            "Stage '%s' completed in %.1f seconds.", name, elapsed
        )
    except Exception as exc:
        elapsed = time.time() - start
        result.duration_seconds = elapsed
        result.status = "failed"
        result.error = str(exc)
        logger.error(
            "Stage '%s' failed after %.1f seconds: %s",
            name,
            elapsed,
            exc,
        )

    return result


def _print_summary(results: list[StageResult]) -> None:
    """Print a summary table of all pipeline stages.

    Args:
        results: List of :class:`StageResult` objects in execution order.
    """
    header = f"{'#':<4} {'Stage':<35} {'Status':<18} {'Duration':<12}"
    separator = "-" * len(header)

    print(f"\n  RNA Entropy Pipeline Summary")
    print(f"  {separator}")
    print(f"  {header}")
    print(f"  {separator}")

    total_time = 0.0
    for idx, r in enumerate(results, start=1):
        duration_str = (
            f"{r.duration_seconds:.1f}s"
            if r.duration_seconds > 0
            else "--"
        )
        total_time += r.duration_seconds
        print(f"  {idx:<4} {r.name:<35} {r.status:<18} {duration_str:<12}")

    print(f"  {separator}")
    print(f"  {'Total':<39} {'':<18} {total_time:.1f}s")
    print(f"  {separator}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Parse arguments and run the RNA entropy pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley RNA Entropy -- Run the full information theory "
            "analysis pipeline for Leishmania infantum RNA targets."
        ),
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Skip step 1 (transcriptome fetch) if FASTA files already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Test mode: log what each stage would do without external calls.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of all stages even if output files already exist.",
    )
    parser.add_argument(
        "--priority-only",
        action="store_true",
        help="Skip stages 1-6 and generate report with priority targets only.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Number of top targets to include (passed to relevant modules).",
    )
    args = parser.parse_args()

    skip_confirm = args.dry_run
    force = args.force

    logger.info("=" * 60)
    logger.info("Marley RNA Entropy Pipeline -- Starting")
    logger.info(
        "Options: skip_fetch=%s  dry_run=%s  force=%s  priority_only=%s  top_n=%s",
        args.skip_fetch,
        args.dry_run,
        args.force,
        args.priority_only,
        args.top_n,
    )
    logger.info("=" * 60)

    # Load pipeline modules.
    try:
        mods = _load_pipeline_modules()
    except ImportError as exc:
        logger.error("Failed to import RNA entropy modules: %s", exc)
        sys.exit(1)

    results: list[StageResult] = []

    # ------------------------------------------------------------------
    # Priority-only mode: skip stages 1-6
    # ------------------------------------------------------------------
    if args.priority_only:
        logger.info("Priority-only mode: skipping stages 1-6.")
        for stage_name in [
            "1. Fetch Transcriptome",
            "2. Codon Usage Analysis",
            "3. Shannon Entropy",
            "4. SL RNA Analysis",
            "5. Human Comparison",
            "6. Structure Prediction",
        ]:
            skipped = StageResult(stage_name)
            skipped.status = "skipped (--priority-only)"
            results.append(skipped)

        r = _run_stage(
            name="7. Report Generation",
            func=mods["report"].generate_rna_report,
            kwargs={"dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)

        _print_summary(results)
        if r.status == "failed":
            sys.exit(1)
        logger.info("RNA entropy pipeline completed (priority-only).")
        return

    # ------------------------------------------------------------------
    # Stage 1: Fetch Transcriptome
    # ------------------------------------------------------------------
    if args.skip_fetch:
        skipped = StageResult("1. Fetch Transcriptome")
        skipped.status = "skipped (--skip-fetch)"
        results.append(skipped)
        logger.info("Skipping transcriptome fetch (--skip-fetch).")
    else:
        r = _run_stage(
            name="1. Fetch Transcriptome",
            func=mods["fetch"].fetch_transcriptome,
            kwargs={"force": force, "dry_run": args.dry_run},
            dry_run=False,
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

    # ------------------------------------------------------------------
    # Stage 2: Codon Usage Analysis
    # ------------------------------------------------------------------
    r = _run_stage(
        name="2. Codon Usage Analysis",
        func=mods["codon"].analyze_codon_usage,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 3: Shannon Entropy
    # ------------------------------------------------------------------
    r = _run_stage(
        name="3. Shannon Entropy",
        func=mods["entropy"].calculate_shannon_entropy,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 4: SL RNA Analysis
    # ------------------------------------------------------------------
    r = _run_stage(
        name="4. SL RNA Analysis",
        func=mods["sl_rna"].analyze_sl_rna,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 5: Human Comparison
    # ------------------------------------------------------------------
    r = _run_stage(
        name="5. Human Comparison",
        func=mods["human"].compare_with_human,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 6: Structure Prediction
    # ------------------------------------------------------------------
    r = _run_stage(
        name="6. Structure Prediction",
        func=mods["structure"].predict_rna_structures,
        kwargs={"force": force, "dry_run": args.dry_run},
        dry_run=False,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 7: Report Generation
    # ------------------------------------------------------------------
    r = _run_stage(
        name="7. Report Generation",
        func=mods["report"].generate_rna_report,
        kwargs={"dry_run": args.dry_run},
        skip_confirm=skip_confirm,
    )
    results.append(r)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    _print_summary(results)

    # Check for any failures.
    failures = [r for r in results if r.status == "failed"]
    if failures:
        logger.warning(
            "%d stage(s) failed. Review the logs above for details.",
            len(failures),
        )
        sys.exit(1)

    logger.info("RNA entropy pipeline completed successfully.")


if __name__ == "__main__":
    main()
