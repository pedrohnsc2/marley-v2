"""Modelo de execucao de pipeline — PipelineRun + StageRecord.

Cada execucao do pipeline e registrada como um PipelineRun que contem
N StageRecords. Os dados sao persistidos em JSON no diretorio
results/runs/<run_id>/.

Uso:
    from core.run import RunManager

    mgr = RunManager()
    run = mgr.create_run("vaccine", parameters={"organism": "leishmania_infantum"})
    mgr.start_stage(run, "01_fetch_genome")
    # ... executa o stage ...
    mgr.complete_stage(run, "01_fetch_genome", key_metrics={"proteins": 8527})
    mgr.complete_run(run)
"""

from __future__ import annotations

import json
import os
import subprocess
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class StageRecord:
    """Resultado de um unico stage dentro de uma execucao."""

    stage_id: str
    name: str
    status: str = "pending"  # pending | running | success | failed | skipped
    started_at: str | None = None
    completed_at: str | None = None
    duration_s: float = 0.0
    error: str | None = None
    error_info: dict[str, Any] | None = None  # structured ErrorInfo
    output_file: str | None = None
    key_metrics: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage_id": self.stage_id,
            "name": self.name,
            "status": self.status,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "duration_s": self.duration_s,
            "error": self.error,
            "error_info": self.error_info,
            "output_file": self.output_file,
            "key_metrics": self.key_metrics,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StageRecord:
        return cls(
            stage_id=data["stage_id"],
            name=data["name"],
            status=data.get("status", "pending"),
            started_at=data.get("started_at"),
            completed_at=data.get("completed_at"),
            duration_s=data.get("duration_s", 0.0),
            error=data.get("error"),
            error_info=data.get("error_info"),
            output_file=data.get("output_file"),
            key_metrics=data.get("key_metrics", {}),
        )


@dataclass
class PipelineRun:
    """Uma execucao completa de um pipeline com parametros especificos."""

    run_id: str
    pipeline: str  # "vaccine" | "aso_math" | "drug" | "docking" | "rna"
    status: str = "created"  # created | running | completed | failed | cancelled
    created_at: str = ""
    started_at: str | None = None
    completed_at: str | None = None
    git_sha: str = ""
    parameters: dict[str, Any] = field(default_factory=dict)
    stages: list[StageRecord] = field(default_factory=list)
    total_duration_s: float = 0.0
    output_dir: str = ""
    notes: str = ""
    tags: list[str] = field(default_factory=list)
    user_id: str | None = None
    team_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "pipeline": self.pipeline,
            "status": self.status,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "git_sha": self.git_sha,
            "parameters": self.parameters,
            "stages": [s.to_dict() for s in self.stages],
            "total_duration_s": self.total_duration_s,
            "output_dir": self.output_dir,
            "notes": self.notes,
            "tags": self.tags,
            "user_id": self.user_id,
            "team_id": self.team_id,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PipelineRun:
        stages = [StageRecord.from_dict(s) for s in data.get("stages", [])]
        return cls(
            run_id=data["run_id"],
            pipeline=data["pipeline"],
            status=data.get("status", "created"),
            created_at=data.get("created_at", ""),
            started_at=data.get("started_at"),
            completed_at=data.get("completed_at"),
            git_sha=data.get("git_sha", ""),
            parameters=data.get("parameters", {}),
            stages=stages,
            total_duration_s=data.get("total_duration_s", 0.0),
            output_dir=data.get("output_dir", ""),
            notes=data.get("notes", ""),
            tags=data.get("tags", []),
            user_id=data.get("user_id"),
            team_id=data.get("team_id"),
        )


# ---------------------------------------------------------------------------
# Run Manager
# ---------------------------------------------------------------------------

def _now_iso() -> str:
    return datetime.now(tz=timezone.utc).isoformat()


def _git_sha() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


class RunManager:
    """Gerencia criacao, atualizacao e persistencia de PipelineRuns.

    Todos os dados ficam em ``results/runs/<run_id>/``.
    Um symlink ``results/latest/<pipeline>`` aponta para o output do
    ultimo run bem-sucedido de cada pipeline.
    """

    def __init__(self, base_dir: Path | None = None, use_supabase: bool = True) -> None:
        if base_dir is None:
            base_dir = Path(__file__).resolve().parent.parent / "results"
        self.base_dir = base_dir
        self.runs_dir = base_dir / "runs"
        self.latest_dir = base_dir / "latest"
        self.use_supabase = use_supabase

    def _run_dir(self, run_id: str) -> Path:
        return self.runs_dir / run_id

    def _metadata_path(self, run_id: str) -> Path:
        return self._run_dir(run_id) / "metadata.json"

    def _save(self, run: PipelineRun) -> Path:
        """Persiste o estado atual do run como metadata.json."""
        run_dir = self._run_dir(run.run_id)
        run_dir.mkdir(parents=True, exist_ok=True)
        meta_path = run_dir / "metadata.json"
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(run.to_dict(), fh, indent=2, ensure_ascii=False)

        # Tambem salva parameters.json separado (mais facil de comparar)
        params_path = run_dir / "parameters.json"
        with open(params_path, "w", encoding="utf-8") as fh:
            json.dump(run.parameters, fh, indent=2, ensure_ascii=False)

        # Sync to Supabase (best-effort — JSON is the source of truth)
        if self.use_supabase:
            try:
                from core.db_runs import upsert_run
                upsert_run(run)
            except Exception:
                pass  # Supabase is best-effort; JSON on disk is authoritative

        return meta_path

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def create_run(
        self,
        pipeline: str,
        parameters: dict[str, Any] | None = None,
        tags: list[str] | None = None,
        notes: str = "",
        user_id: str | None = None,
        team_id: str | None = None,
    ) -> PipelineRun:
        """Cria um novo PipelineRun e persiste no disco."""
        now = _now_iso()
        date_prefix = datetime.now(tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
        short_id = uuid.uuid4().hex[:8]
        run_id = f"{date_prefix}_{pipeline}_{short_id}"

        run = PipelineRun(
            run_id=run_id,
            pipeline=pipeline,
            status="created",
            created_at=now,
            git_sha=_git_sha(),
            parameters=parameters or {},
            output_dir=str(self._run_dir(run_id) / "outputs"),
            notes=notes,
            tags=tags or [],
            user_id=user_id,
            team_id=team_id,
        )

        # Cria diretorios
        run_dir = self._run_dir(run_id)
        (run_dir / "stages").mkdir(parents=True, exist_ok=True)
        (run_dir / "outputs").mkdir(parents=True, exist_ok=True)

        self._save(run)
        return run

    def start_run(self, run: PipelineRun) -> None:
        """Marca o run como em execucao."""
        run.status = "running"
        run.started_at = _now_iso()
        self._save(run)

    def start_stage(self, run: PipelineRun, stage_id: str, name: str = "") -> StageRecord:
        """Registra o inicio de um stage."""
        stage = StageRecord(
            stage_id=stage_id,
            name=name or stage_id,
            status="running",
            started_at=_now_iso(),
        )
        run.stages.append(stage)
        self._save(run)
        return stage

    def complete_stage(
        self,
        run: PipelineRun,
        stage_id: str,
        status: str = "success",
        error: str | None = None,
        error_info: dict[str, Any] | None = None,
        output_file: str | None = None,
        key_metrics: dict[str, Any] | None = None,
    ) -> None:
        """Marca um stage como completo."""
        for stage in run.stages:
            if stage.stage_id == stage_id:
                stage.status = status
                stage.completed_at = _now_iso()
                if stage.started_at:
                    start = datetime.fromisoformat(stage.started_at)
                    end = datetime.fromisoformat(stage.completed_at)
                    stage.duration_s = round((end - start).total_seconds(), 2)
                stage.error = error
                stage.error_info = error_info
                stage.output_file = output_file
                stage.key_metrics = key_metrics or {}

                # Salva resultado do stage individualmente
                stage_dir = self._run_dir(run.run_id) / "stages"
                stage_path = stage_dir / f"{stage_id}.json"
                with open(stage_path, "w", encoding="utf-8") as fh:
                    json.dump(stage.to_dict(), fh, indent=2, ensure_ascii=False)
                break

        self._save(run)

    def complete_run(self, run: PipelineRun, status: str = "completed") -> None:
        """Finaliza o run e atualiza symlink latest."""
        run.status = status
        run.completed_at = _now_iso()
        run.total_duration_s = sum(s.duration_s for s in run.stages)

        # Se algum stage falhou, marca o run como failed
        if any(s.status == "failed" for s in run.stages):
            run.status = "failed"

        self._save(run)

        # Atualiza symlink latest/ se sucesso
        if run.status == "completed":
            self._update_latest(run)

    def fail_run(self, run: PipelineRun, error: str = "") -> None:
        """Marca o run como falho."""
        run.status = "failed"
        run.completed_at = _now_iso()
        run.total_duration_s = sum(s.duration_s for s in run.stages)
        if error:
            run.notes = f"{run.notes}\nError: {error}".strip()
        self._save(run)

    def cancel_run(self, run: PipelineRun) -> None:
        """Marca o run como cancelado."""
        run.status = "cancelled"
        run.completed_at = _now_iso()
        run.total_duration_s = sum(s.duration_s for s in run.stages)
        self._save(run)

    # ------------------------------------------------------------------
    # Latest symlinks
    # ------------------------------------------------------------------

    def _update_latest(self, run: PipelineRun) -> None:
        """Atualiza o symlink results/latest/<pipeline> para o output do run."""
        self.latest_dir.mkdir(parents=True, exist_ok=True)
        link = self.latest_dir / run.pipeline
        target = Path(run.output_dir)

        # Remove symlink existente
        if link.is_symlink() or link.exists():
            link.unlink()

        # Cria novo symlink relativo
        try:
            rel_target = os.path.relpath(target, self.latest_dir)
            link.symlink_to(rel_target)
        except OSError:
            # Fallback para symlink absoluto
            link.symlink_to(target)

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def load_run(self, run_id: str) -> PipelineRun:
        """Carrega um run existente do disco."""
        meta_path = self._metadata_path(run_id)
        if not meta_path.exists():
            raise FileNotFoundError(f"Run nao encontrado: {run_id}")
        with open(meta_path, encoding="utf-8") as fh:
            return PipelineRun.from_dict(json.load(fh))

    def list_runs(self, pipeline: str | None = None) -> list[PipelineRun]:
        """Lista todos os runs, opcionalmente filtrados por pipeline."""
        runs: list[PipelineRun] = []
        if not self.runs_dir.exists():
            return runs

        for entry in sorted(self.runs_dir.iterdir(), reverse=True):
            meta = entry / "metadata.json"
            if not meta.exists():
                continue
            try:
                with open(meta, encoding="utf-8") as fh:
                    data = json.load(fh)
                if pipeline and data.get("pipeline") != pipeline:
                    continue
                runs.append(PipelineRun.from_dict(data))
            except (json.JSONDecodeError, KeyError):
                continue

        return runs

    def get_run_output_dir(self, run_id: str) -> Path:
        """Retorna o diretorio de outputs de um run."""
        return self._run_dir(run_id) / "outputs"

    def compare_runs(self, run_id_a: str, run_id_b: str) -> dict[str, Any]:
        """Compara dois runs retornando diffs de parametros e metricas."""
        run_a = self.load_run(run_id_a)
        run_b = self.load_run(run_id_b)

        # Diff de parametros
        all_keys = set(run_a.parameters) | set(run_b.parameters)
        param_diff: list[dict[str, Any]] = []
        for key in sorted(all_keys):
            val_a = run_a.parameters.get(key)
            val_b = run_b.parameters.get(key)
            if val_a != val_b:
                param_diff.append({"param": key, "run_a": val_a, "run_b": val_b})

        # Metricas por stage
        metrics_a = {s.stage_id: s.key_metrics for s in run_a.stages}
        metrics_b = {s.stage_id: s.key_metrics for s in run_b.stages}

        return {
            "run_a": run_id_a,
            "run_b": run_id_b,
            "pipeline": run_a.pipeline,
            "parameter_diff": param_diff,
            "metrics_a": metrics_a,
            "metrics_b": metrics_b,
            "duration_a": run_a.total_duration_s,
            "duration_b": run_b.total_duration_s,
        }
