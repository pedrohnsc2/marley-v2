"""Launcher para execucao de pipeline via web.

Le JSON do stdin, cria um PipelineRun, faz fork da execucao
em processo background, e retorna o run_id imediatamente.

Uso (chamado pelo Next.js API route):
    echo '{"pipeline":"vaccine","preset":"leishmania_infantum"}' | .venv/bin/python -m core.launcher

Modo de execucao (chamado internamente pelo fork):
    .venv/bin/python -m core.launcher --execute <run_id>
"""

from __future__ import annotations

import json
import subprocess
import sys
import traceback
from pathlib import Path

from core.logger import get_logger
from core.registry import PIPELINES
from core.run import RunManager

logger = get_logger("launcher")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYTHON_BIN = PROJECT_ROOT / ".venv" / "bin" / "python"


def _validate_input(data: dict) -> None:
    """Validate launcher input. Raises ValueError on invalid input."""
    pipeline = data.get("pipeline")
    if not pipeline or pipeline not in PIPELINES:
        available = ", ".join(sorted(PIPELINES.keys()))
        raise ValueError(f"Invalid pipeline: '{pipeline}'. Available: {available}")

    preset = data.get("preset", "")
    if preset and not all(c.isalnum() or c == "_" for c in preset):
        raise ValueError(f"Invalid preset name: '{preset}'")


def _launch(data: dict) -> dict:
    """Create a run and fork execution. Returns run info."""
    from core.run_cli import create_pipeline_run, _build_parameters

    pipeline_id = data["pipeline"]

    # Build parameters from preset + overrides
    # SECURITY: config_path is always None — never accept file paths from web
    parameters = _build_parameters(
        pipeline_id,
        preset=data.get("preset"),
        config_path=None,
        param_overrides=None,
    )
    # Merge inline parameter overrides from the web form
    if data.get("parameters") and isinstance(data["parameters"], dict):
        parameters.update(data["parameters"])

    tags = data.get("tags", [])
    if not isinstance(tags, list):
        tags = []
    notes = str(data.get("notes", ""))[:500]

    # Extract optional user/team ownership
    user_id = data.get("user_id") if isinstance(data.get("user_id"), str) else None
    team_id = data.get("team_id") if isinstance(data.get("team_id"), str) else None

    # Create the run (writes to disk + Supabase)
    run = create_pipeline_run(
        pipeline_id=pipeline_id,
        parameters=parameters,
        tags=tags,
        notes=notes,
        user_id=user_id,
        team_id=team_id,
    )

    # Fork execution into detached background process
    run_dir = Path(run.output_dir).parent
    log_path = run_dir / "output.log"
    log_file = open(log_path, "w")

    child = subprocess.Popen(
        [str(PYTHON_BIN), "-m", "core.launcher", "--execute", run.run_id],
        cwd=str(PROJECT_ROOT),
        stdout=log_file,
        stderr=log_file,
        start_new_session=True,
    )

    # Write PID file for crash detection
    pid_path = run_dir / "worker.pid"
    pid_path.write_text(str(child.pid))

    return {"run_id": run.run_id, "status": "created", "pid": child.pid}


def _execute(run_id: str) -> None:
    """Execute a pre-created run. Called by the forked process."""
    from core.run_cli import execute_pipeline_run

    mgr = RunManager()

    try:
        run = mgr.load_run(run_id)
        execute_pipeline_run(run, run.parameters)
    except Exception as exc:
        logger.error("Run %s crashed: %s", run_id, traceback.format_exc())
        try:
            from core.errors import classify_error
            err_info = classify_error(exc)
            run = mgr.load_run(run_id)
            mgr.fail_run(run, error=traceback.format_exc()[:2000])
            # Attach structured error info to run notes for diagnostics
            run.notes = (
                f"{run.notes}\n[{err_info.code}] {err_info.suggestion}".strip()
            )
            mgr._save(run)
        except Exception:
            pass
        raise
    finally:
        # Clean up PID file
        try:
            pid_path = Path(mgr._run_dir(run_id)) / "worker.pid"
            pid_path.unlink(missing_ok=True)
        except Exception:
            pass


def main() -> None:
    if len(sys.argv) >= 3 and sys.argv[1] == "--execute":
        run_id = sys.argv[2]
        # Validate run_id format (alphanumeric + underscore + hyphen)
        if not all(c.isalnum() or c in "_-" for c in run_id):
            logger.error("Invalid run_id: %s", run_id)
            sys.exit(1)
        logger.info("Executing run: %s", run_id)
        _execute(run_id)
        return

    # Launch mode: read JSON from stdin, create run, fork execution
    try:
        raw = sys.stdin.read()
        data = json.loads(raw)
        _validate_input(data)
        result = _launch(data)
        print(json.dumps(result))
    except (json.JSONDecodeError, ValueError) as exc:
        print(json.dumps({"error": str(exc)}))
        sys.exit(1)
    except Exception:
        print(json.dumps({"error": "Internal launcher error"}))
        sys.exit(1)


if __name__ == "__main__":
    main()
