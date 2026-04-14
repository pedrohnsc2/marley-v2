"""Tests for core.db_runs — Supabase persistence with mocked client."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from core.run import PipelineRun, StageRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_mock_client() -> MagicMock:
    """Create a mock Supabase client with chained table().upsert().execute()."""
    mock_client = MagicMock()
    mock_table = MagicMock()
    mock_client.table.return_value = mock_table
    mock_table.upsert.return_value = mock_table
    mock_table.select.return_value = mock_table
    mock_table.eq.return_value = mock_table
    mock_table.order.return_value = mock_table

    mock_response = MagicMock()
    mock_response.data = []
    mock_table.execute.return_value = mock_response

    return mock_client


def _sample_run(
    user_id: str | None = "usr-1",
    team_id: str | None = "team-2",
) -> PipelineRun:
    """Build a sample PipelineRun for testing."""
    return PipelineRun(
        run_id="20260414_100000_vaccine_abc12345",
        pipeline="vaccine",
        status="running",
        created_at="2026-04-14T10:00:00+00:00",
        started_at="2026-04-14T10:00:01+00:00",
        git_sha="abc1234",
        parameters={"organism": "leishmania_infantum"},
        total_duration_s=0.0,
        output_dir="/results/runs/test/outputs",
        notes="test",
        tags=["ci"],
        user_id=user_id,
        team_id=team_id,
    )


# ---------------------------------------------------------------------------
# upsert_run
# ---------------------------------------------------------------------------


class TestUpsertRun:
    """core.db_runs.upsert_run() sends the correct payload to Supabase."""

    def test_payload_includes_user_and_team(self) -> None:
        mock_client = _make_mock_client()
        run = _sample_run(user_id="usr-abc", team_id="team-xyz")

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_run
            upsert_run(run)

        # Extract the payload passed to .upsert()
        mock_client.table.assert_called_with("pipeline_runs")
        call_args = mock_client.table().upsert.call_args
        payload = call_args[0][0]

        assert payload["user_id"] == "usr-abc"
        assert payload["team_id"] == "team-xyz"
        assert payload["run_id"] == run.run_id
        assert payload["pipeline"] == "vaccine"

    def test_payload_user_team_none(self) -> None:
        """When user_id/team_id are None, payload should contain None."""
        mock_client = _make_mock_client()
        run = _sample_run(user_id=None, team_id=None)

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_run
            upsert_run(run)

        call_args = mock_client.table().upsert.call_args
        payload = call_args[0][0]
        assert payload["user_id"] is None
        assert payload["team_id"] is None

    def test_payload_has_all_expected_keys(self) -> None:
        mock_client = _make_mock_client()
        run = _sample_run()

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_run
            upsert_run(run)

        call_args = mock_client.table().upsert.call_args
        payload = call_args[0][0]

        expected_keys = {
            "run_id", "pipeline", "status", "created_at", "started_at",
            "completed_at", "git_sha", "parameters", "total_duration_s",
            "output_dir", "notes", "tags", "user_id", "team_id",
        }
        assert set(payload.keys()) == expected_keys

    def test_swallows_exceptions_silently(self) -> None:
        """upsert_run should not raise even if Supabase fails."""
        mock_client = MagicMock()
        mock_client.table.side_effect = Exception("Supabase down")
        run = _sample_run()

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_run
            # Should not raise
            upsert_run(run)

    def test_swallows_client_init_error(self) -> None:
        """If _get_client itself fails, upsert_run should still not raise."""
        run = _sample_run()

        with patch("core.db_runs._get_client", side_effect=Exception("no creds")):
            from core.db_runs import upsert_run
            upsert_run(run)


# ---------------------------------------------------------------------------
# upsert_stage
# ---------------------------------------------------------------------------


class TestUpsertStage:
    """core.db_runs.upsert_stage() sends stage data to Supabase."""

    def test_upsert_stage_payload(self) -> None:
        mock_client = _make_mock_client()
        stage = StageRecord(
            stage_id="01_fetch",
            name="Fetch Genome",
            status="success",
            started_at="2026-04-14T10:00:00+00:00",
            completed_at="2026-04-14T10:01:00+00:00",
            duration_s=60.0,
            error=None,
            error_info=None,
            key_metrics={"proteins": 100},
        )

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_stage
            upsert_stage("run-123", stage)

        mock_client.table.assert_called_with("pipeline_stages")
        call_args = mock_client.table().upsert.call_args
        payload = call_args[0][0]

        assert payload["run_id"] == "run-123"
        assert payload["stage_id"] == "01_fetch"
        assert payload["name"] == "Fetch Genome"
        assert payload["status"] == "success"
        assert payload["duration_s"] == 60.0
        assert payload["key_metrics"] == {"proteins": 100}

    def test_upsert_stage_with_error_info(self) -> None:
        mock_client = _make_mock_client()
        stage = StageRecord(
            stage_id="02_filter",
            name="Filter",
            status="failed",
            error="timed out",
            error_info={"code": "network.timeout"},
        )

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_stage
            upsert_stage("run-456", stage)

        call_args = mock_client.table().upsert.call_args
        payload = call_args[0][0]
        assert payload["error"] == "timed out"
        assert payload["error_info"] == {"code": "network.timeout"}

    def test_upsert_stage_swallows_exceptions(self) -> None:
        mock_client = MagicMock()
        mock_client.table.side_effect = Exception("DB error")
        stage = StageRecord(stage_id="x", name="X")

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import upsert_stage
            upsert_stage("run-789", stage)  # should not raise


# ---------------------------------------------------------------------------
# get_run
# ---------------------------------------------------------------------------


class TestGetRun:
    """core.db_runs.get_run() fetches a run or returns None."""

    def test_get_run_found(self) -> None:
        mock_client = _make_mock_client()
        expected_row = {"run_id": "run-1", "pipeline": "vaccine", "status": "completed"}
        mock_response = MagicMock()
        mock_response.data = [expected_row]
        mock_client.table().select().eq().execute.return_value = mock_response

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import get_run
            result = get_run("run-1")

        assert result == expected_row

    def test_get_run_not_found(self) -> None:
        mock_client = _make_mock_client()
        mock_response = MagicMock()
        mock_response.data = []
        mock_client.table().select().eq().execute.return_value = mock_response

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import get_run
            result = get_run("nonexistent")

        assert result is None

    def test_get_run_returns_none_on_error(self) -> None:
        """On any exception, get_run should return None."""
        with patch("core.db_runs._get_client", side_effect=Exception("DB down")):
            from core.db_runs import get_run
            result = get_run("any-id")

        assert result is None

    def test_get_run_returns_none_on_query_error(self) -> None:
        """If the query itself raises, get_run should still return None."""
        mock_client = MagicMock()
        mock_client.table.return_value.select.side_effect = Exception("query err")

        with patch("core.db_runs._get_client", return_value=mock_client):
            from core.db_runs import get_run
            result = get_run("err-id")

        assert result is None
