"""Marley v4 -- Vaccine Optimization Pipeline.

Optimizes the v1 mRNA vaccine construct through:
1. Epitope optimization via APL (Altered Peptide Ligand) design
2. Cross-reactivity safety check against dog proteome
3. Adjuvant screening (6 candidates)
4. Construct reconstruction with comparison report

Usage:
    python run_v4_pipeline.py
    python run_v4_pipeline.py --dry-run
    python run_v4_pipeline.py --force
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from typing import Any, Callable

from core.logger import get_logger

logger = get_logger("v4_pipeline")

# ---------------------------------------------------------------------------
# Module imports via importlib (numbered module names)
# ---------------------------------------------------------------------------


def _load_pipeline_modules() -> dict[str, Any]:
    """Import numbered v4 pipeline modules using importlib.

    Returns:
        Dictionary mapping short names to imported module objects.
    """
    modules = {
        "optimize": importlib.import_module("pipeline.08_epitope_optimize"),
        "safety": importlib.import_module("pipeline.09_safety_check"),
        "adjuvant": importlib.import_module("pipeline.10_adjuvant_screen"),
        "reconstruct": importlib.import_module("pipeline.11_reconstruct"),
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

    print(f"\n  v4 Pipeline Summary")
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
    """Parse arguments and run the v4 optimization pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley v4 -- Run the vaccine optimization pipeline to improve "
            "the v1 mRNA construct through epitope optimization, safety "
            "checks, adjuvant screening, and construct reconstruction."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Test mode: log what each stage would do without calling "
            "external APIs or writing final outputs."
        ),
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
    logger.info("Marley v4 Optimization Pipeline -- Starting")
    logger.info(
        "Options: dry_run=%s  force=%s",
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
    # Stage 1: Epitope Optimization
    # ------------------------------------------------------------------
    r = _run_stage(
        name="1. Epitope Optimization",
        func=mods["optimize"].optimize_epitopes,
        kwargs={"force": force, "dry_run": args.dry_run},
        dry_run=False,  # Module handles dry_run internally
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 2: Safety Check
    # ------------------------------------------------------------------
    r = _run_stage(
        name="2. Safety Check",
        func=mods["safety"].run_safety_check,
        kwargs={"force": force, "dry_run": args.dry_run},
        dry_run=False,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 3: Adjuvant Screening
    # ------------------------------------------------------------------
    r = _run_stage(
        name="3. Adjuvant Screening",
        func=mods["adjuvant"].screen_adjuvants,
        kwargs={"force": force, "dry_run": args.dry_run},
        dry_run=False,
        skip_confirm=skip_confirm,
    )
    results.append(r)
    if r.status == "skipped (user quit)":
        _print_summary(results)
        return

    # ------------------------------------------------------------------
    # Stage 4: Reconstruct & Compare
    # ------------------------------------------------------------------
    r = _run_stage(
        name="4. Reconstruct & Compare",
        func=mods["reconstruct"].reconstruct_optimized,
        kwargs={"force": force, "dry_run": args.dry_run},
        dry_run=False,
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

    logger.info("v4 Optimization Pipeline completed successfully.")


if __name__ == "__main__":
    main()
