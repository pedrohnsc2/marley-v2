"""FastAPI service exposing pipeline execution to the Next.js frontend.

Replaces the execFileSync(python -m core.launcher) approach used in
local dev. In production, the Next.js container calls this API over
the internal Docker network instead of spawning a subprocess.

Endpoints:
    POST /api/runs/start     — enqueue a pipeline run
    GET  /api/runs            — list runs (delegates to RunManager)
    GET  /api/runs/{run_id}   — get single run metadata
    GET  /api/presets/{pipeline} — list available presets
    GET  /health              — health check

Depends on Redis for the job queue (via rq — simpler than Celery for
a small team, zero config, same Redis).
"""

from __future__ import annotations

import hmac
import os
import re
from typing import Any

from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, field_validator

load_dotenv()

app = FastAPI(title="Marley Worker API", version="1.0.0")

# CORS: allow the Next.js frontend origin
ALLOWED_ORIGINS = os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type", "X-API-Key", "Authorization"],
)

# ---------------------------------------------------------------------------
# Auth
# ---------------------------------------------------------------------------

API_SECRET_KEY = os.getenv("API_SECRET_KEY", "")


def _check_auth(request: Request) -> None:
    """Validate X-API-Key header. Uses constant-time comparison."""
    if not API_SECRET_KEY:
        return  # dev mode — skip auth
    provided = request.headers.get("x-api-key", "")
    if not hmac.compare_digest(provided, API_SECRET_KEY):
        raise HTTPException(status_code=401, detail="Unauthorized")


# ---------------------------------------------------------------------------
# Request/Response models
# ---------------------------------------------------------------------------

ALLOWED_PIPELINES = {"vaccine", "drug", "docking", "rna", "aso_math"}
PRESET_RE = re.compile(r"^[a-z0-9_]{1,50}$")
PARAM_KEY_RE = re.compile(r"^[a-zA-Z0-9_]+$")


class StartRunRequest(BaseModel):
    pipeline: str
    preset: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    tags: list[str] = Field(default_factory=list, max_length=10)
    notes: str = Field(default="", max_length=500)
    user_id: str | None = None
    team_id: str | None = None

    @field_validator("pipeline")
    @classmethod
    def validate_pipeline(cls, v: str) -> str:
        if v not in ALLOWED_PIPELINES:
            raise ValueError(f"Invalid pipeline. Allowed: {sorted(ALLOWED_PIPELINES)}")
        return v

    @field_validator("preset")
    @classmethod
    def validate_preset(cls, v: str | None) -> str | None:
        if v is not None and not PRESET_RE.match(v):
            raise ValueError("Invalid preset name")
        return v

    @field_validator("parameters")
    @classmethod
    def validate_parameters(cls, v: dict[str, Any]) -> dict[str, Any]:
        for key, value in v.items():
            if not PARAM_KEY_RE.match(key):
                raise ValueError(f"Invalid parameter key: {key}")
            if isinstance(value, str) and (".." in value or "/" in value or "\\" in value):
                raise ValueError("Invalid parameter value")
        return v


class StartRunResponse(BaseModel):
    run_id: str
    status: str
    job_id: str


# ---------------------------------------------------------------------------
# Job queue (rq + Redis)
# ---------------------------------------------------------------------------

import redis
from rq import Queue

REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")
redis_conn = redis.from_url(REDIS_URL)
pipeline_queue = Queue("pipelines", connection=redis_conn, default_timeout=1800)  # 30 min max


def _execute_run_job(run_id: str, parameters: dict[str, Any]) -> None:
    """Job function executed by rq worker in a separate process."""
    from core.run import RunManager
    from core.run_cli import execute_pipeline_run

    mgr = RunManager()
    run = mgr.load_run(run_id)
    execute_pipeline_run(run, parameters)


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------


@app.get("/health")
async def health() -> dict[str, str]:
    """Health check — verifies Redis connectivity."""
    try:
        redis_conn.ping()
    except Exception:
        raise HTTPException(status_code=503, detail="Redis unreachable")
    return {"status": "healthy"}


@app.post("/api/runs/start", response_model=StartRunResponse, status_code=201)
async def start_run(body: StartRunRequest, request: Request) -> StartRunResponse:
    """Create a pipeline run and enqueue it for background execution."""
    _check_auth(request)

    # Check queue depth (concurrency control)
    queued = len(pipeline_queue)
    if queued >= 10:
        raise HTTPException(status_code=429, detail="Too many queued runs. Try again later.")

    from core.run_cli import create_pipeline_run, _build_parameters

    # Build parameters
    parameters = _build_parameters(
        body.pipeline,
        preset=body.preset,
        config_path=None,  # never accept file paths from HTTP
        param_overrides=None,
    )
    if body.parameters:
        parameters.update(body.parameters)

    # Create the run record (persists to disk + Supabase)
    run = create_pipeline_run(
        pipeline_id=body.pipeline,
        parameters=parameters,
        tags=body.tags,
        notes=body.notes,
        user_id=body.user_id,
        team_id=body.team_id,
    )

    # Enqueue for background execution
    job = pipeline_queue.enqueue(
        _execute_run_job,
        run.run_id,
        parameters,
        job_id=f"run-{run.run_id}",
        job_timeout=1800,
    )

    return StartRunResponse(
        run_id=run.run_id,
        status="created",
        job_id=job.id,
    )


@app.get("/api/runs")
async def list_runs(request: Request, pipeline: str | None = None) -> list[dict[str, Any]]:
    """List pipeline runs, optionally filtered by pipeline name."""
    _check_auth(request)
    from core.run import RunManager

    mgr = RunManager()
    runs = mgr.list_runs(pipeline=pipeline)
    return [r.to_dict() for r in runs]


@app.get("/api/runs/{run_id}")
async def get_run(run_id: str, request: Request) -> dict[str, Any]:
    """Get metadata for a specific run."""
    _check_auth(request)
    if ".." in run_id or "/" in run_id or "\\" in run_id:
        raise HTTPException(status_code=400, detail="Invalid run_id")

    from core.run import RunManager

    mgr = RunManager()
    try:
        run = mgr.load_run(run_id)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Run not found")
    return run.to_dict()


@app.get("/api/presets/{pipeline}")
async def list_presets(pipeline: str, request: Request) -> list[str]:
    """List available preset names for a pipeline."""
    _check_auth(request)
    from core.run_cli import _list_preset_names

    return _list_preset_names(pipeline)


@app.get("/api/queue/status")
async def queue_status(request: Request) -> dict[str, int]:
    """Return current queue depth and worker count."""
    _check_auth(request)
    from rq import Worker

    workers = Worker.all(connection=redis_conn)
    return {
        "queued": len(pipeline_queue),
        "active_workers": len(workers),
        "busy_workers": sum(1 for w in workers if w.get_state() == "busy"),
    }
