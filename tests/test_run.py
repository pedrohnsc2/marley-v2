"""Tests for core.run — PipelineRun, StageRecord, and RunManager."""

from __future__ import annotations

import json

import pytest

from core.run import PipelineRun, RunManager, StageRecord


# ---------------------------------------------------------------------------
# StageRecord data model
# ---------------------------------------------------------------------------


class TestStageRecordRoundTrip:
    """StageRecord.to_dict / from_dict round-trip."""

    def test_round_trip_all_fields(self) -> None:
        stage = StageRecord(
            stage_id="01_fetch",
            name="Fetch Genome",
            status="success",
            started_at="2026-04-14T10:00:00+00:00",
            completed_at="2026-04-14T10:01:00+00:00",
            duration_s=60.0,
            error=None,
            error_info={"code": "network.timeout", "category": "network",
                        "message": "timed out", "stage_id": "01_fetch",
                        "suggestion": "retry"},
            output_file="/tmp/out.fasta",
            key_metrics={"proteins": 8527},
        )
        d = stage.to_dict()
        restored = StageRecord.from_dict(d)
        assert restored.stage_id == stage.stage_id
        assert restored.name == stage.name
        assert restored.status == stage.status
        assert restored.started_at == stage.started_at
        assert restored.completed_at == stage.completed_at
        assert restored.duration_s == stage.duration_s
        assert restored.error == stage.error
        assert restored.error_info == stage.error_info
        assert restored.output_file == stage.output_file
        assert restored.key_metrics == stage.key_metrics

    def test_round_trip_minimal(self) -> None:
        stage = StageRecord(stage_id="02_filter", name="Filter")
        d = stage.to_dict()
        restored = StageRecord.from_dict(d)
        assert restored.stage_id == "02_filter"
        assert restored.name == "Filter"
        assert restored.status == "pending"
        assert restored.error_info is None

    def test_from_dict_missing_optional_fields(self) -> None:
        """Minimal dict without optional keys should still deserialize."""
        data = {"stage_id": "03_score", "name": "Score"}
        stage = StageRecord.from_dict(data)
        assert stage.status == "pending"
        assert stage.duration_s == 0.0
        assert stage.key_metrics == {}
        assert stage.error_info is None

    def test_to_dict_keys(self) -> None:
        stage = StageRecord(stage_id="x", name="X")
        expected_keys = {
            "stage_id", "name", "status", "started_at", "completed_at",
            "duration_s", "error", "error_info", "output_file", "key_metrics",
        }
        assert set(stage.to_dict().keys()) == expected_keys


# ---------------------------------------------------------------------------
# PipelineRun data model
# ---------------------------------------------------------------------------


class TestPipelineRunRoundTrip:
    """PipelineRun.to_dict / from_dict round-trip."""

    def test_round_trip_full(self) -> None:
        stage = StageRecord(stage_id="01_fetch", name="Fetch")
        run = PipelineRun(
            run_id="20260414_100000_vaccine_abc12345",
            pipeline="vaccine",
            status="completed",
            created_at="2026-04-14T10:00:00+00:00",
            started_at="2026-04-14T10:00:01+00:00",
            completed_at="2026-04-14T10:05:00+00:00",
            git_sha="abc1234",
            parameters={"organism": "leishmania_infantum"},
            stages=[stage],
            total_duration_s=300.0,
            output_dir="/tmp/output",
            notes="test run",
            tags=["ci", "nightly"],
            user_id="user-123",
            team_id="team-456",
        )
        d = run.to_dict()
        restored = PipelineRun.from_dict(d)

        assert restored.run_id == run.run_id
        assert restored.pipeline == run.pipeline
        assert restored.status == run.status
        assert restored.created_at == run.created_at
        assert restored.started_at == run.started_at
        assert restored.completed_at == run.completed_at
        assert restored.git_sha == run.git_sha
        assert restored.parameters == run.parameters
        assert len(restored.stages) == 1
        assert restored.stages[0].stage_id == "01_fetch"
        assert restored.total_duration_s == run.total_duration_s
        assert restored.output_dir == run.output_dir
        assert restored.notes == run.notes
        assert restored.tags == run.tags
        assert restored.user_id == "user-123"
        assert restored.team_id == "team-456"

    def test_from_dict_missing_user_team_backward_compat(self) -> None:
        """Old JSON without user_id/team_id should deserialize with None."""
        data = {
            "run_id": "old_run_001",
            "pipeline": "vaccine",
            "status": "completed",
            "parameters": {},
        }
        run = PipelineRun.from_dict(data)
        assert run.user_id is None
        assert run.team_id is None
        assert run.run_id == "old_run_001"

    def test_from_dict_all_fields_populated(self) -> None:
        data = {
            "run_id": "run_full",
            "pipeline": "drug",
            "status": "running",
            "created_at": "2026-04-14T00:00:00+00:00",
            "started_at": "2026-04-14T00:00:01+00:00",
            "completed_at": None,
            "git_sha": "def5678",
            "parameters": {"top_n": 10},
            "stages": [],
            "total_duration_s": 0.0,
            "output_dir": "/results/run_full/outputs",
            "notes": "in progress",
            "tags": ["experiment"],
            "user_id": "usr_aaa",
            "team_id": "team_bbb",
        }
        run = PipelineRun.from_dict(data)
        assert run.user_id == "usr_aaa"
        assert run.team_id == "team_bbb"
        assert run.pipeline == "drug"
        assert run.tags == ["experiment"]

    def test_to_dict_includes_user_and_team(self) -> None:
        run = PipelineRun(
            run_id="r1", pipeline="vaccine",
            user_id="u1", team_id="t1",
        )
        d = run.to_dict()
        assert "user_id" in d
        assert "team_id" in d
        assert d["user_id"] == "u1"
        assert d["team_id"] == "t1"

    def test_to_dict_user_team_none_when_absent(self) -> None:
        run = PipelineRun(run_id="r2", pipeline="vaccine")
        d = run.to_dict()
        assert d["user_id"] is None
        assert d["team_id"] is None


# ---------------------------------------------------------------------------
# RunManager — create
# ---------------------------------------------------------------------------


class TestRunManagerCreate:
    """RunManager.create_run() writes the correct filesystem structure."""

    def test_create_run_creates_dirs_and_files(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine", parameters={"organism": "leishmania"})

        run_dir = mgr.runs_dir / run.run_id
        assert run_dir.exists()
        assert (run_dir / "metadata.json").exists()
        assert (run_dir / "parameters.json").exists()
        assert (run_dir / "stages").is_dir()
        assert (run_dir / "outputs").is_dir()

    def test_create_run_with_user_id_and_team_id(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run(
            "vaccine",
            user_id="user-abc",
            team_id="team-xyz",
        )

        meta_path = mgr.runs_dir / run.run_id / "metadata.json"
        with open(meta_path) as fh:
            meta = json.load(fh)

        assert meta["user_id"] == "user-abc"
        assert meta["team_id"] == "team-xyz"
        assert run.user_id == "user-abc"
        assert run.team_id == "team-xyz"

    def test_create_run_without_user_id_backward_compat(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")

        meta_path = mgr.runs_dir / run.run_id / "metadata.json"
        with open(meta_path) as fh:
            meta = json.load(fh)

        assert meta["user_id"] is None
        assert meta["team_id"] is None
        assert run.user_id is None

    def test_create_run_parameters_persisted(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        params = {"organism": "leishmania_infantum", "threshold": 0.7}
        run = mgr.create_run("vaccine", parameters=params)

        params_path = mgr.runs_dir / run.run_id / "parameters.json"
        with open(params_path) as fh:
            loaded = json.load(fh)
        assert loaded == params

    def test_create_run_status_is_created(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("drug")
        assert run.status == "created"

    def test_create_run_tags_and_notes(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine", tags=["ci"], notes="test note")
        assert run.tags == ["ci"]
        assert run.notes == "test note"

    def test_create_run_id_contains_pipeline_name(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("aso_math")
        assert "aso_math" in run.run_id


# ---------------------------------------------------------------------------
# RunManager — lifecycle
# ---------------------------------------------------------------------------


class TestRunManagerLifecycle:
    """RunManager lifecycle transitions: start, stage, complete, fail, cancel."""

    def test_start_run_sets_status_and_started_at(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        assert run.started_at is None

        mgr.start_run(run)
        assert run.status == "running"
        assert run.started_at is not None

    def test_start_and_complete_stage(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)

        stage = mgr.start_stage(run, "01_fetch", "Fetch Genome")
        assert stage.status == "running"
        assert stage.started_at is not None
        assert len(run.stages) == 1

        mgr.complete_stage(
            run, "01_fetch",
            status="success",
            key_metrics={"proteins": 100},
        )

        completed = run.stages[0]
        assert completed.status == "success"
        assert completed.completed_at is not None
        assert completed.key_metrics == {"proteins": 100}

        # Stage JSON file should exist
        stage_file = mgr.runs_dir / run.run_id / "stages" / "01_fetch.json"
        assert stage_file.exists()

    def test_complete_stage_with_error(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)
        mgr.start_stage(run, "02_filter", "Filter")

        error_info = {"code": "network.timeout", "category": "network",
                      "message": "timed out", "stage_id": "02_filter",
                      "suggestion": "retry"}
        mgr.complete_stage(
            run, "02_filter",
            status="failed",
            error="Connection timed out",
            error_info=error_info,
        )

        assert run.stages[0].status == "failed"
        assert run.stages[0].error == "Connection timed out"
        assert run.stages[0].error_info == error_info

    def test_complete_run_sets_status_and_duration(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)
        mgr.start_stage(run, "01_fetch", "Fetch")
        mgr.complete_stage(run, "01_fetch", status="success")
        mgr.complete_run(run)

        assert run.status == "completed"
        assert run.completed_at is not None
        assert run.total_duration_s >= 0

    def test_complete_run_failed_stage_marks_run_failed(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)
        mgr.start_stage(run, "01_fetch", "Fetch")
        mgr.complete_stage(run, "01_fetch", status="failed", error="boom")
        mgr.complete_run(run)

        assert run.status == "failed"

    def test_complete_run_creates_latest_symlink(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)
        mgr.complete_run(run, status="completed")

        link = mgr.latest_dir / "vaccine"
        assert link.is_symlink() or link.exists()

    def test_fail_run(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)
        mgr.fail_run(run, error="Fatal crash")

        assert run.status == "failed"
        assert run.completed_at is not None
        assert "Fatal crash" in run.notes

    def test_cancel_run(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run("vaccine")
        mgr.start_run(run)
        mgr.cancel_run(run)

        assert run.status == "cancelled"
        assert run.completed_at is not None


# ---------------------------------------------------------------------------
# RunManager — queries
# ---------------------------------------------------------------------------


class TestRunManagerLoad:
    """RunManager.load_run() restores a run from disk."""

    def test_load_run(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run = mgr.create_run(
            "vaccine",
            user_id="u1",
            team_id="t1",
            parameters={"organism": "leishmania"},
        )

        loaded = mgr.load_run(run.run_id)
        assert loaded.run_id == run.run_id
        assert loaded.pipeline == "vaccine"
        assert loaded.user_id == "u1"
        assert loaded.team_id == "t1"
        assert loaded.parameters == {"organism": "leishmania"}

    def test_load_run_nonexistent_raises(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        with pytest.raises(FileNotFoundError, match="Run nao encontrado"):
            mgr.load_run("nonexistent_run_id")


class TestRunManagerList:
    """RunManager.list_runs() returns runs from disk."""

    def test_list_runs_empty(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        assert mgr.list_runs() == []

    def test_list_runs_returns_all(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        mgr.create_run("vaccine")
        mgr.create_run("drug")
        mgr.create_run("vaccine")

        runs = mgr.list_runs()
        assert len(runs) == 3

    def test_list_runs_filter_by_pipeline(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        mgr.create_run("vaccine")
        mgr.create_run("drug")
        mgr.create_run("vaccine")

        vaccine_runs = mgr.list_runs(pipeline="vaccine")
        assert len(vaccine_runs) == 2
        assert all(r.pipeline == "vaccine" for r in vaccine_runs)

        drug_runs = mgr.list_runs(pipeline="drug")
        assert len(drug_runs) == 1

    def test_list_runs_nonexistent_pipeline_returns_empty(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        mgr.create_run("vaccine")
        assert mgr.list_runs(pipeline="nonexistent") == []


class TestRunManagerCompare:
    """RunManager.compare_runs() diffs parameters and metrics."""

    def test_compare_runs_parameter_diff(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run_a = mgr.create_run("vaccine", parameters={"threshold": 0.5, "organism": "leish"})
        run_b = mgr.create_run("vaccine", parameters={"threshold": 0.7, "organism": "leish"})

        comparison = mgr.compare_runs(run_a.run_id, run_b.run_id)
        assert comparison["run_a"] == run_a.run_id
        assert comparison["run_b"] == run_b.run_id
        assert comparison["pipeline"] == "vaccine"

        # threshold differs, organism does not
        diff_params = {d["param"] for d in comparison["parameter_diff"]}
        assert "threshold" in diff_params
        assert "organism" not in diff_params

    def test_compare_runs_identical_params(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run_a = mgr.create_run("vaccine", parameters={"x": 1})
        run_b = mgr.create_run("vaccine", parameters={"x": 1})

        comparison = mgr.compare_runs(run_a.run_id, run_b.run_id)
        assert comparison["parameter_diff"] == []

    def test_compare_runs_nonexistent_raises(self, tmp_path: object) -> None:
        mgr = RunManager(base_dir=tmp_path, use_supabase=False)  # type: ignore[arg-type]
        run_a = mgr.create_run("vaccine")
        with pytest.raises(FileNotFoundError):
            mgr.compare_runs(run_a.run_id, "nonexistent")
