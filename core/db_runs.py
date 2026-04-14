"""Supabase persistence for pipeline run tracking.

Best-effort sync layer: all writes are fire-and-forget.  JSON on disk
(managed by ``RunManager``) remains the source of truth.  If Supabase
is unreachable the pipeline execution is not affected.

Tables expected in Supabase:
    - ``pipeline_runs``  (PK: run_id)
    - ``pipeline_stages`` (PK: run_id + stage_id)
"""

from __future__ import annotations

from typing import Any

from core.logger import get_logger
from core.run import PipelineRun, StageRecord

logger = get_logger("db_runs")

RUNS_TABLE = "pipeline_runs"
STAGES_TABLE = "pipeline_stages"


def _get_client():
    """Reuse the lazily-initialised Supabase client from core.db."""
    from core.db import _get_client as get_db_client

    return get_db_client()


# ---------------------------------------------------------------------------
# Write helpers
# ---------------------------------------------------------------------------


def upsert_run(run: PipelineRun) -> None:
    """Insert or update a pipeline run row.

    Silently swallows errors so the caller never has to worry about
    Supabase availability.
    """
    try:
        client = _get_client()
        payload: dict[str, Any] = {
            "run_id": run.run_id,
            "pipeline": run.pipeline,
            "status": run.status,
            "created_at": run.created_at or None,
            "started_at": run.started_at,
            "completed_at": run.completed_at,
            "git_sha": run.git_sha,
            "parameters": run.parameters,
            "total_duration_s": run.total_duration_s,
            "output_dir": run.output_dir,
            "notes": run.notes,
            "tags": run.tags,
        }
        client.table(RUNS_TABLE).upsert(payload, on_conflict="run_id").execute()
        logger.info("Upserted run %s (status=%s)", run.run_id, run.status)
    except Exception:
        logger.warning(
            "Failed to sync run %s to Supabase (continuing with local JSON)",
            run.run_id,
        )


def upsert_stage(run_id: str, stage: StageRecord) -> None:
    """Insert or update a pipeline stage row.

    Silently swallows errors so the caller never has to worry about
    Supabase availability.
    """
    try:
        client = _get_client()
        payload: dict[str, Any] = {
            "run_id": run_id,
            "stage_id": stage.stage_id,
            "name": stage.name,
            "status": stage.status,
            "started_at": stage.started_at,
            "completed_at": stage.completed_at,
            "duration_s": stage.duration_s,
            "error": stage.error,
            "error_info": stage.error_info,
            "key_metrics": stage.key_metrics,
        }
        client.table(STAGES_TABLE).upsert(
            payload, on_conflict="run_id,stage_id"
        ).execute()
    except Exception:
        logger.warning(
            "Failed to sync stage %s/%s to Supabase", run_id, stage.stage_id
        )


# ---------------------------------------------------------------------------
# Read helpers
# ---------------------------------------------------------------------------


def get_run(run_id: str) -> dict[str, Any] | None:
    """Fetch a single run from Supabase by run_id.

    Returns ``None`` on any error or if the row does not exist.
    """
    try:
        client = _get_client()
        response = (
            client.table(RUNS_TABLE).select("*").eq("run_id", run_id).execute()
        )
        return response.data[0] if response.data else None
    except Exception:
        logger.warning("Failed to fetch run %s from Supabase", run_id)
        return None


def get_stages(run_id: str) -> list[dict[str, Any]]:
    """Fetch all stages for a run from Supabase.

    Returns an empty list on any error.
    """
    try:
        client = _get_client()
        response = (
            client.table(STAGES_TABLE)
            .select("*")
            .eq("run_id", run_id)
            .order("stage_id")
            .execute()
        )
        return response.data or []
    except Exception:
        logger.warning("Failed to fetch stages for run %s from Supabase", run_id)
        return []
