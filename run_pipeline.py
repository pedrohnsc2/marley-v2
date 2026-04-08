"""Entry point for running the full reverse vaccinology pipeline end-to-end.

Executes all five pipeline stages sequentially:

1. Fetch genome (proteome download from TriTrypDB)
2. Filter surface proteins (SignalP 6.0 signal peptide screen)
3. Conservation analysis (BLAST against Leishmania strains)
4. Immunogenicity scoring (IEDB MHC-I binding predictions)
5. Report generation (CSV + Markdown summary)

Usage:
    python run_pipeline.py                 # full run with confirmations
    python run_pipeline.py --skip-fetch    # skip proteome download
    python run_pipeline.py --dry-run       # test mode, no API calls
    python run_pipeline.py --force         # force re-run all stages
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from typing import Any, Callable

from core.logger import get_logger

logger = get_logger("pipeline")

# ---------------------------------------------------------------------------
# Module imports via importlib (numbered module names)
# ---------------------------------------------------------------------------


def _load_pipeline_modules() -> dict[str, Any]:
    """Import numbered pipeline modules using importlib.

    Returns:
        Dictionary mapping short names to imported module objects.
    """
    modules = {
        "fetch": importlib.import_module("pipeline.01_fetch_genome"),
        "filter": importlib.import_module("pipeline.02_filter_surface"),
        "conservation": importlib.import_module("pipeline.03_conservation"),
        "immunogenicity": importlib.import_module("pipeline.04_immunogenicity"),
        "report": importlib.import_module("pipeline.05_report"),
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

    print(f"\n  Pipeline Summary")
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
    """Parse arguments and run the full pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley -- Run the full reverse vaccinology pipeline "
            "for Leishmania infantum vaccine candidate discovery."
        ),
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Skip step 1 (genome fetch) if the FASTA file already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Test mode: log what each stage would do without calling external APIs.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of all stages even if output files already exist.",
    )
    args = parser.parse_args()

    skip_confirm = args.dry_run
    force = args.force

    logger.info("=" * 60)
    logger.info("Marley Pipeline -- Starting")
    logger.info(
        "Options: skip_fetch=%s  dry_run=%s  force=%s",
        args.skip_fetch,
        args.dry_run,
        args.force,
    )
    logger.info("=" * 60)

    # Load pipeline modules.
    try:
        mods = _load_pipeline_modules()
    except ImportError as exc:
        logger.error("Failed to import pipeline modules: %s", exc)
        sys.exit(1)

    results: list[StageResult] = []

    # ------------------------------------------------------------------
    # Stage 1: Fetch Genome
    # ------------------------------------------------------------------
    if args.skip_fetch:
        skipped = StageResult("1. Fetch Genome")
        skipped.status = "skipped (--skip-fetch)"
        results.append(skipped)
        logger.info("Skipping genome fetch (--skip-fetch).")
    else:
        r = _run_stage(
            name="1. Fetch Genome",
            func=mods["fetch"].fetch_genome,
            kwargs={"force": force},
            dry_run=args.dry_run,
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

    # ------------------------------------------------------------------
    # Stage 2: Filter Surface Proteins
    # ------------------------------------------------------------------
    r = _run_stage(
        name="2. Filter Surface Proteins",
        func=mods["filter"].filter_surface_proteins,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 3: Conservation Analysis
    # ------------------------------------------------------------------
    r = _run_stage(
        name="3. Conservation Analysis",
        func=mods["conservation"].analyze_conservation,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 4: Immunogenicity Scoring
    # ------------------------------------------------------------------
    r = _run_stage(
        name="4. Immunogenicity Scoring",
        func=mods["immunogenicity"].score_immunogenicity,
        kwargs={"force": force},
        dry_run=args.dry_run,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 5: Report Generation
    # ------------------------------------------------------------------
    r = _run_stage(
        name="5. Report Generation",
        func=mods["report"].generate_report,
        kwargs={},
        dry_run=args.dry_run,
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

    logger.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
