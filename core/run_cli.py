"""CLI unificado para execucao de pipelines do Marley.

Ponto de entrada unico que substitui os scripts individuais
(run_pipeline.py, drug_targets/run_drug_targets.py, etc.)
adicionando rastreamento de runs, parametrizacao e comparacao.

Uso:
    # Listar pipelines disponiveis
    python -m core.run_cli --list-pipelines

    # Listar presets de um pipeline
    python -m core.run_cli --list-presets vaccine

    # Executar com preset default
    python -m core.run_cli --pipeline vaccine

    # Executar com preset especifico
    python -m core.run_cli --pipeline vaccine --preset leishmania_infantum

    # Override de parametros
    python -m core.run_cli --pipeline vaccine --param conservation_threshold=0.7

    # Executar com config JSON custom
    python -m core.run_cli --pipeline vaccine --config path/to/config.json

    # Listar runs anteriores
    python -m core.run_cli --list-runs vaccine

    # Comparar dois runs
    python -m core.run_cli --compare RUN_ID_A RUN_ID_B

    # Dry-run (simula sem executar stages reais)
    python -m core.run_cli --pipeline vaccine --dry-run
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
import sys
import time
from pathlib import Path
from typing import Any

from core.logger import get_logger
from core.registry import PIPELINES, get_pipeline, list_pipelines
from core.run import RunManager, PipelineRun

logger = get_logger("run_cli")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PRESETS_DIR = PROJECT_ROOT / "configs" / "presets"

# Security: frozen set of allowed module paths from the registry.
# importlib.import_module() will only load modules in this set.
_ALLOWED_MODULES: frozenset[str] = frozenset(
    stage.module_path
    for pipe in PIPELINES.values()
    for stage in pipe.stages
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_preset(pipeline_id: str, preset_name: str) -> dict[str, Any]:
    """Carrega um arquivo de preset JSON."""
    preset_path = PRESETS_DIR / pipeline_id / f"{preset_name}.json"
    if not preset_path.exists():
        raise FileNotFoundError(
            f"Preset nao encontrado: {preset_path}\n"
            f"Presets disponiveis: {_list_preset_names(pipeline_id)}"
        )
    with open(preset_path, encoding="utf-8") as fh:
        return json.load(fh)


def _list_preset_names(pipeline_id: str) -> list[str]:
    """Lista nomes de presets disponiveis para um pipeline."""
    preset_dir = PRESETS_DIR / pipeline_id
    if not preset_dir.exists():
        return []
    return sorted(p.stem for p in preset_dir.glob("*.json"))


def _parse_param_override(param_str: str) -> tuple[str, Any]:
    """Parseia um override 'key=value' com inferencia de tipo."""
    if "=" not in param_str:
        raise ValueError(f"Formato invalido: '{param_str}'. Use 'key=value'.")
    key, val_str = param_str.split("=", 1)
    key = key.strip()
    val_str = val_str.strip()

    # Inferencia de tipo
    if val_str.lower() in ("true", "false"):
        return key, val_str.lower() == "true"
    try:
        return key, int(val_str)
    except ValueError:
        pass
    try:
        return key, float(val_str)
    except ValueError:
        pass
    return key, val_str


def _build_parameters(
    pipeline_id: str,
    preset: str | None,
    config_path: str | None,
    param_overrides: list[str] | None,
) -> dict[str, Any]:
    """Constroi o dicionario de parametros a partir de preset + overrides."""
    params: dict[str, Any] = {}

    if config_path:
        # Config JSON explicita tem prioridade sobre preset
        with open(config_path, encoding="utf-8") as fh:
            params = json.load(fh)
    elif preset:
        params = _load_preset(pipeline_id, preset)
    else:
        # Tenta carregar preset default
        pipe_def = get_pipeline(pipeline_id)
        try:
            params = _load_preset(pipeline_id, pipe_def.default_preset)
        except FileNotFoundError:
            pass

    # Aplica overrides
    if param_overrides:
        for override in param_overrides:
            key, val = _parse_param_override(override)
            params[key] = val

    return params


# ---------------------------------------------------------------------------
# Pipeline execution (split into create + execute for launcher reuse)
# ---------------------------------------------------------------------------

def create_pipeline_run(
    pipeline_id: str,
    parameters: dict[str, Any],
    tags: list[str] | None = None,
    notes: str = "",
    dry_run: bool = False,
    user_id: str | None = None,
    team_id: str | None = None,
) -> PipelineRun:
    """Cria um PipelineRun sem executar stages. Retorna o run persistido."""
    mgr = RunManager()

    run_tags = tags or []
    if dry_run:
        run_tags.append("dry-run")

    run = mgr.create_run(
        pipeline=pipeline_id,
        parameters=parameters,
        tags=run_tags,
        notes=notes,
        user_id=user_id,
        team_id=team_id,
    )
    return run


def execute_pipeline_run(
    run: PipelineRun,
    parameters: dict[str, Any],
    dry_run: bool = False,
) -> None:
    """Executa os stages de um PipelineRun ja criado."""
    pipe_def = get_pipeline(run.pipeline)
    mgr = RunManager()

    logger.info("=" * 60)
    logger.info("Marley Pipeline Runner")
    logger.info("  Pipeline: %s (%s)", pipe_def.display_name, run.pipeline)
    logger.info("  Run ID:   %s", run.run_id)
    logger.info("  Stages:   %d", len(pipe_def.stages))
    logger.info("  Params:   %s", json.dumps(parameters, indent=2, default=str))
    logger.info("=" * 60)

    mgr.start_run(run)

    for stage_def in pipe_def.stages:
        logger.info("-" * 60)
        logger.info("Stage: %s — %s", stage_def.stage_id, stage_def.name)
        logger.info("-" * 60)

        mgr.start_stage(run, stage_def.stage_id, stage_def.name)

        if dry_run:
            logger.info("[DRY RUN] Skipping: %s", stage_def.name)
            mgr.complete_stage(
                run, stage_def.stage_id,
                status="success",
                key_metrics={"dry_run": True},
            )
            continue

        # Security: only import modules from the trusted registry
        if stage_def.module_path not in _ALLOWED_MODULES:
            mgr.complete_stage(
                run, stage_def.stage_id,
                status="failed",
                error=f"Untrusted module: {stage_def.module_path}",
            )
            logger.error("BLOCKED untrusted module: %s", stage_def.module_path)
            continue

        start = time.time()
        try:
            mod = importlib.import_module(stage_def.module_path)
            func = getattr(mod, stage_def.entry_function)

            sig = inspect.signature(func)
            kwargs: dict[str, Any] = {}

            if "force" in sig.parameters:
                kwargs["force"] = parameters.get("force_rerun", False)
            if "dry_run" in sig.parameters:
                kwargs["dry_run"] = False
            if "config" in sig.parameters:
                # aso_math modules expect a TargetConfig object, not a dict
                if run.pipeline == "aso_math":
                    from aso_math.target_config import TargetConfig
                    valid_fields = {f.name for f in __import__("dataclasses").fields(TargetConfig)}
                    filtered = {k: v for k, v in parameters.items() if k in valid_fields}
                    kwargs["config"] = TargetConfig(**filtered)
                else:
                    kwargs["config"] = parameters
            if "top_n" in sig.parameters:
                kwargs["top_n"] = parameters.get("top_n_targets", 5)

            result = func(**kwargs)
            elapsed = round(time.time() - start, 2)

            key_metrics: dict[str, Any] = {"duration_s": elapsed}
            if isinstance(result, dict):
                key_metrics.update(result)
            elif isinstance(result, str) and result:
                key_metrics["output"] = result[:200]

            mgr.complete_stage(
                run, stage_def.stage_id,
                status="success",
                key_metrics=key_metrics,
            )
            logger.info("Stage %s: OK (%.2f s)", stage_def.stage_id, elapsed)

        except Exception as exc:
            elapsed = round(time.time() - start, 2)
            from core.errors import classify_error
            err_info = classify_error(exc, stage_id=stage_def.stage_id)
            mgr.complete_stage(
                run, stage_def.stage_id,
                status="failed",
                error=str(exc),
                error_info=err_info.to_dict(),
                key_metrics={"duration_s": elapsed},
            )
            logger.error("Stage %s FAILED [%s]: %s", stage_def.stage_id, err_info.code, exc)

    mgr.complete_run(run)
    _print_run_summary(run)


def _execute_pipeline(
    pipeline_id: str,
    parameters: dict[str, Any],
    tags: list[str] | None = None,
    notes: str = "",
    dry_run: bool = False,
) -> PipelineRun:
    """Executa um pipeline completo (create + execute). Backwards compatible."""
    run = create_pipeline_run(pipeline_id, parameters, tags, notes, dry_run)
    execute_pipeline_run(run, parameters, dry_run)
    return run


def _print_run_summary(run: PipelineRun) -> None:
    """Imprime resumo formatado de um run."""
    header = f"{'#':<4} {'Stage':<40} {'Status':<12} {'Duration':<10}"
    sep = "-" * len(header)

    print(f"\n  Run Summary: {run.run_id}")
    print(f"  Pipeline:    {run.pipeline}")
    print(f"  Status:      {run.status}")
    print(f"  {sep}")
    print(f"  {header}")
    print(f"  {sep}")

    for idx, stage in enumerate(run.stages, 1):
        dur = f"{stage.duration_s:.1f}s" if stage.duration_s > 0 else "--"
        status_icon = "OK" if stage.status == "success" else stage.status.upper()
        print(f"  {idx:<4} {stage.name:<40} {status_icon:<12} {dur:<10}")

    print(f"  {sep}")
    print(f"  {'Total':<44} {'':<12} {run.total_duration_s:.1f}s")
    print(f"  {sep}")
    print(f"  Output: {run.output_dir}")
    print()


def _print_runs_table(runs: list[PipelineRun]) -> None:
    """Imprime tabela de runs."""
    if not runs:
        print("  Nenhum run encontrado.")
        return

    header = f"{'Run ID':<45} {'Pipeline':<12} {'Status':<12} {'Duration':<10} {'Date':<20}"
    sep = "-" * len(header)

    print(f"\n  {sep}")
    print(f"  {header}")
    print(f"  {sep}")

    for run in runs:
        date = run.created_at[:19] if run.created_at else "--"
        dur = f"{run.total_duration_s:.1f}s" if run.total_duration_s > 0 else "--"
        print(f"  {run.run_id:<45} {run.pipeline:<12} {run.status:<12} {dur:<10} {date:<20}")

    print(f"  {sep}")
    print(f"  Total: {len(runs)} run(s)")
    print()


def _print_comparison(comparison: dict[str, Any]) -> None:
    """Imprime comparacao entre dois runs."""
    print(f"\n  Comparison: {comparison['run_a']} vs {comparison['run_b']}")
    print(f"  Pipeline: {comparison['pipeline']}")
    print()

    diffs = comparison["parameter_diff"]
    if diffs:
        print("  Parameter Differences:")
        print(f"  {'Parameter':<30} {'Run A':<20} {'Run B':<20}")
        print(f"  {'-' * 70}")
        for d in diffs:
            print(f"  {d['param']:<30} {str(d['run_a']):<20} {str(d['run_b']):<20}")
    else:
        print("  Parameters: identical")

    print()
    print(f"  Duration: {comparison['duration_a']:.1f}s vs {comparison['duration_b']:.1f}s")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Marley Pipeline Runner — execucao dinamica e parametrizada",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:
  python -m core.run_cli --list-pipelines
  python -m core.run_cli --pipeline vaccine --dry-run
  python -m core.run_cli --pipeline aso_math --preset trypanosoma_cruzi
  python -m core.run_cli --pipeline drug --param top_n_targets=30
  python -m core.run_cli --list-runs vaccine
  python -m core.run_cli --compare RUN_A RUN_B
""",
    )

    # Modo de operacao (mutuamente exclusivo de forma logica)
    parser.add_argument(
        "--pipeline", "-p",
        type=str,
        help="Pipeline a executar (vaccine, drug, docking, rna, aso_math)",
    )
    parser.add_argument(
        "--list-pipelines",
        action="store_true",
        help="Listar pipelines disponiveis",
    )
    parser.add_argument(
        "--list-presets",
        type=str,
        metavar="PIPELINE",
        help="Listar presets disponiveis para um pipeline",
    )
    parser.add_argument(
        "--list-runs",
        type=str,
        nargs="?",
        const="__all__",
        metavar="PIPELINE",
        help="Listar runs anteriores (opcional: filtrar por pipeline)",
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("RUN_A", "RUN_B"),
        help="Comparar dois runs",
    )

    # Parametrizacao
    parser.add_argument(
        "--preset",
        type=str,
        help="Nome do preset (ex: leishmania_infantum, trypanosoma_cruzi)",
    )
    parser.add_argument(
        "--config",
        type=str,
        help="Caminho para arquivo JSON de configuracao",
    )
    parser.add_argument(
        "--param",
        action="append",
        help="Override de parametro (key=value). Pode ser repetido.",
    )

    # Controle
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simula execucao sem rodar stages reais",
    )
    parser.add_argument(
        "--tag",
        action="append",
        help="Tag para o run (pode ser repetido)",
    )
    parser.add_argument(
        "--notes",
        type=str,
        default="",
        help="Notas descritivas para o run",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # List pipelines
    # ------------------------------------------------------------------
    if args.list_pipelines:
        print("\n  Available Pipelines:")
        print(f"  {'-' * 70}")
        for p in list_pipelines():
            print(f"  {p['id']:<12} {p['name']:<30} ({p['stages']} stages)")
            print(f"  {'':12} {p['description']}")
            print()
        return

    # ------------------------------------------------------------------
    # List presets
    # ------------------------------------------------------------------
    if args.list_presets:
        names = _list_preset_names(args.list_presets)
        if not names:
            print(f"  Nenhum preset encontrado para '{args.list_presets}'.")
            return
        print(f"\n  Presets for '{args.list_presets}':")
        for name in names:
            print(f"    - {name}")
        print()
        return

    # ------------------------------------------------------------------
    # List runs
    # ------------------------------------------------------------------
    if args.list_runs is not None:
        mgr = RunManager()
        pipeline_filter = None if args.list_runs == "__all__" else args.list_runs
        runs = mgr.list_runs(pipeline=pipeline_filter)
        _print_runs_table(runs)
        return

    # ------------------------------------------------------------------
    # Compare runs
    # ------------------------------------------------------------------
    if args.compare:
        mgr = RunManager()
        try:
            comparison = mgr.compare_runs(args.compare[0], args.compare[1])
            _print_comparison(comparison)
        except FileNotFoundError as exc:
            print(f"  Erro: {exc}")
            sys.exit(1)
        return

    # ------------------------------------------------------------------
    # Execute pipeline
    # ------------------------------------------------------------------
    if not args.pipeline:
        parser.print_help()
        sys.exit(1)

    try:
        parameters = _build_parameters(
            args.pipeline,
            preset=args.preset,
            config_path=args.config,
            param_overrides=args.param,
        )

        run = _execute_pipeline(
            pipeline_id=args.pipeline,
            parameters=parameters,
            tags=args.tag,
            notes=args.notes,
            dry_run=args.dry_run,
        )

        # Exit code baseado no status
        if run.status == "failed":
            sys.exit(1)

    except (ValueError, FileNotFoundError) as exc:
        logger.error(str(exc))
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Execucao cancelada pelo usuario.")
        sys.exit(130)


if __name__ == "__main__":
    main()
