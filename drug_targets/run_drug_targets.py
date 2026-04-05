"""Entry point for the Marley v2/v3 drug target discovery pipeline.

Executes pipeline stages sequentially:

v2 (drug target discovery):
1. Fetch enzymes (UniProt download by metabolic pathway)
2. Human comparison (BLAST against human proteome)
3. Essentiality check (DEG + TriTrypDB knockout data)
4. Druggability scoring (composite score + AlphaFold links)
5. Report generation (CSV + Markdown summary)

v3 (molecular docking, activated with --docking):
6. Fetch structures (AlphaFold PDB + PDBQT preparation)
7. Compound library (ChEMBL + repurposing library)
8. Docking (AutoDock Vina)
9. ADMET filter (Lipinski + pkCSM)
10. Docking report (top hits + PyMOL scripts)

Usage:
    python -m drug_targets.run_drug_targets
    python -m drug_targets.run_drug_targets --dry-run
    python -m drug_targets.run_drug_targets --priority-only
    python -m drug_targets.run_drug_targets --docking --top-n 5
    python -m drug_targets.run_drug_targets --docking --dry-run
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from typing import Any, Callable

from core.logger import get_logger

logger = get_logger("drug_targets_pipeline")

# ---------------------------------------------------------------------------
# Module imports via importlib (numbered module names)
# ---------------------------------------------------------------------------


def _load_modules() -> dict[str, Any]:
    """Import numbered drug target pipeline modules using importlib.

    Returns:
        Dictionary mapping short names to imported module objects.
    """
    modules = {
        "fetch": importlib.import_module("drug_targets.01_fetch_enzymes"),
        "comparison": importlib.import_module("drug_targets.02_human_comparison"),
        "essentiality": importlib.import_module("drug_targets.03_essentiality"),
        "druggability": importlib.import_module("drug_targets.04_druggability"),
        "report": importlib.import_module("drug_targets.05_report"),
    }
    return modules


def _load_docking_modules() -> dict[str, Any]:
    """Import v3 docking pipeline modules using importlib.

    Returns:
        Dictionary mapping short names to imported module objects.
    """
    modules = {
        "structures": importlib.import_module("drug_targets.06_fetch_structures"),
        "compounds": importlib.import_module("drug_targets.07_compound_library"),
        "docking": importlib.import_module("drug_targets.08_docking"),
        "admet": importlib.import_module("drug_targets.09_admet_filter"),
        "docking_report": importlib.import_module("drug_targets.10_docking_report"),
    }
    return modules


# ---------------------------------------------------------------------------
# Stage execution
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
        dry_run: If ``True``, pass dry_run=True to the stage function.
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

    # Execute the stage.
    logger.info("Starting stage: %s", name)
    start = time.time()

    try:
        if dry_run:
            kwargs["dry_run"] = True
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
    header = f"{'#':<4} {'Stage':<40} {'Status':<18} {'Duration':<12}"
    separator = "-" * len(header)

    print(f"\n  Drug Target Pipeline Summary")
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
        print(f"  {idx:<4} {r.name:<40} {r.status:<18} {duration_str:<12}")

    print(f"  {separator}")
    print(f"  {'Total':<44} {'':<18} {total_time:.1f}s")
    print(f"  {separator}\n")


def _print_pathway_summary() -> None:
    """Print a summary of targets per metabolic pathway."""
    try:
        from core.db import get_all_drug_targets

        targets = get_all_drug_targets()
        if not targets:
            return

        pathway_counts: dict[str, int] = {}
        for t in targets:
            pathway = t.pathway or "unknown"
            pathway_counts[pathway] = pathway_counts.get(pathway, 0) + 1

        print("\n  Targets by Metabolic Pathway")
        print("  " + "-" * 40)
        for pathway, count in sorted(
            pathway_counts.items(), key=lambda x: x[1], reverse=True
        ):
            print(f"    {pathway:<30} {count}")
        print("  " + "-" * 40)
        print(f"    {'Total':<30} {len(targets)}\n")

    except Exception:
        logger.warning("Could not generate pathway summary (database unavailable).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Parse arguments and run the drug target discovery pipeline."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley v2/v3 -- Drug target discovery + molecular docking "
            "pipeline for Leishmania infantum enzymatic targets."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Test mode: log what each stage would do without calling external APIs.",
    )
    parser.add_argument(
        "--priority-only",
        action="store_true",
        help="Load only pre-validated priority targets (skip API-dependent stages).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of all stages even if output files already exist.",
    )
    parser.add_argument(
        "--docking",
        action="store_true",
        help="Run v3 molecular docking stages (06-10) after drug target discovery.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top targets to use for docking (default: 5).",
    )
    parser.add_argument(
        "--max-compounds",
        type=int,
        default=100,
        help="Maximum compounds per target for docking (default: 100).",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=16,
        help="AutoDock Vina exhaustiveness parameter (default: 16).",
    )
    parser.add_argument(
        "--blind-docking",
        action="store_true",
        help="Use blind docking (large grid box covering entire protein).",
    )
    parser.add_argument(
        "--predict-admet",
        action="store_true",
        help="Run pkCSM ADMET predictions (slower, requires API access).",
    )
    args = parser.parse_args()

    skip_confirm = args.dry_run
    force = args.force

    logger.info("=" * 60)
    logger.info("Marley v2/v3 -- Drug Target Discovery Pipeline")
    logger.info(
        "Options: dry_run=%s  priority_only=%s  force=%s  docking=%s",
        args.dry_run,
        args.priority_only,
        args.force,
        args.docking,
    )
    logger.info("=" * 60)

    # Load modules.
    try:
        mods = _load_modules()
    except ImportError as exc:
        logger.error("Failed to import drug target modules: %s", exc)
        sys.exit(1)

    results: list[StageResult] = []

    if args.priority_only:
        # --priority-only: skip stages 1-3, run only druggability + report.
        logger.info("Priority-only mode: skipping fetch, comparison, and essentiality.")

        for stage_name in [
            "1. Fetch Enzymes",
            "2. Human Comparison",
            "3. Essentiality Check",
        ]:
            skipped = StageResult(stage_name)
            skipped.status = "skipped (--priority-only)"
            results.append(skipped)

        # Stage 4: Druggability (loads priority targets).
        r = _run_stage(
            name="4. Druggability Scoring",
            func=mods["druggability"].score_druggability,
            kwargs={"force": force, "dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 5: Report.
        r = _run_stage(
            name="5. Report Generation",
            func=mods["report"].generate_drug_targets_report,
            kwargs={"dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)

    else:
        # Full pipeline: all 5 stages.

        # Stage 1: Fetch Enzymes.
        r = _run_stage(
            name="1. Fetch Enzymes",
            func=mods["fetch"].fetch_enzymes,
            kwargs={"force": force, "dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 2: Human Comparison.
        r = _run_stage(
            name="2. Human Comparison",
            func=mods["comparison"].compare_with_human,
            kwargs={"force": force, "dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 3: Essentiality Check.
        r = _run_stage(
            name="3. Essentiality Check",
            func=mods["essentiality"].check_essentiality,
            kwargs={"force": force, "dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 4: Druggability Scoring.
        r = _run_stage(
            name="4. Druggability Scoring",
            func=mods["druggability"].score_druggability,
            kwargs={"force": force, "dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 5: Report Generation.
        r = _run_stage(
            name="5. Report Generation",
            func=mods["report"].generate_drug_targets_report,
            kwargs={"dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)

    # ------------------------------------------------------------------
    # v3: Molecular Docking (stages 6-10, if --docking)
    # ------------------------------------------------------------------
    if args.docking:
        logger.info("=" * 60)
        logger.info("Marley v3 -- Molecular Docking Pipeline")
        logger.info("=" * 60)

        try:
            docking_mods = _load_docking_modules()
        except ImportError as exc:
            logger.error("Failed to import docking modules: %s", exc)
            logger.error(
                "Install docking dependencies: pip install rdkit-pypi meeko vina openbabel-wheel chembl-webresource-client"
            )
            sys.exit(1)

        # Stage 6: Fetch Structures.
        r = _run_stage(
            name="6. Fetch Structures (AlphaFold)",
            func=docking_mods["structures"].fetch_and_prepare_structures,
            kwargs={"top_n": args.top_n, "force": force, "dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 7: Compound Library.
        r = _run_stage(
            name="7. Compound Library (ChEMBL)",
            func=docking_mods["compounds"].build_compound_library,
            kwargs={
                "top_n": args.top_n,
                "max_compounds": args.max_compounds,
                "force": force,
                "dry_run": args.dry_run,
            },
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 8: Molecular Docking.
        r = _run_stage(
            name="8. Molecular Docking (Vina)",
            func=docking_mods["docking"].run_docking_campaign,
            kwargs={
                "top_n": args.top_n,
                "exhaustiveness": args.exhaustiveness,
                "blind_docking": args.blind_docking,
                "force": force,
                "dry_run": args.dry_run,
            },
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 9: ADMET Filter.
        r = _run_stage(
            name="9. ADMET Filter (Lipinski + pkCSM)",
            func=docking_mods["admet"].filter_admet,
            kwargs={
                "predict_admet": args.predict_admet,
                "force": force,
                "dry_run": args.dry_run,
            },
            skip_confirm=skip_confirm,
        )
        results.append(r)
        if r.status == "skipped (user quit)":
            _print_summary(results)
            return

        # Stage 10: Docking Report.
        r = _run_stage(
            name="10. Docking Report",
            func=docking_mods["docking_report"].generate_docking_report,
            kwargs={"dry_run": args.dry_run},
            skip_confirm=skip_confirm,
        )
        results.append(r)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    _print_summary(results)

    if not args.dry_run:
        _print_pathway_summary()

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
