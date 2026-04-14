"""Tests for core.worker_api — FastAPI endpoints.

FastAPI + redis + rq may not be installed locally. All tests are skipped
gracefully via pytest.importorskip() if the dependencies are missing.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

# Skip entire module if fastapi/httpx/redis/rq are not installed
fastapi = pytest.importorskip("fastapi")
pytest.importorskip("httpx")
pytest.importorskip("redis")
pytest.importorskip("rq")


# ---------------------------------------------------------------------------
# Fixtures — mock Redis + rq before importing the app
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _mock_redis_and_queue(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace the module-level redis_conn and pipeline_queue with mocks.

    This must happen before any request hits the app, so it's autouse.
    """
    mock_redis = MagicMock()
    mock_redis.ping.return_value = True

    mock_queue = MagicMock()
    mock_queue.__len__ = MagicMock(return_value=0)

    import core.worker_api as wa
    monkeypatch.setattr(wa, "redis_conn", mock_redis)
    monkeypatch.setattr(wa, "pipeline_queue", mock_queue)


@pytest.fixture()
def client():
    """Create a TestClient bound to the worker_api app."""
    from starlette.testclient import TestClient

    from core.worker_api import app
    return TestClient(app)


# ---------------------------------------------------------------------------
# Health check
# ---------------------------------------------------------------------------


class TestHealthEndpoint:
    """GET /health returns service status."""

    def test_health_ok(self, client) -> None:
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json()["status"] == "healthy"

    def test_health_redis_down(self, client, monkeypatch: pytest.MonkeyPatch) -> None:
        """When Redis is unreachable, /health should return 503."""
        import core.worker_api as wa

        mock_redis = MagicMock()
        mock_redis.ping.side_effect = Exception("Connection refused")
        monkeypatch.setattr(wa, "redis_conn", mock_redis)

        resp = client.get("/health")
        assert resp.status_code == 503


# ---------------------------------------------------------------------------
# POST /api/runs/start
# ---------------------------------------------------------------------------


class TestStartRunEndpoint:
    """POST /api/runs/start creates and enqueues a pipeline run."""

    def test_invalid_pipeline_returns_422(self, client) -> None:
        resp = client.post("/api/runs/start", json={
            "pipeline": "nonexistent",
        })
        assert resp.status_code == 422

    def test_valid_pipeline_returns_201(self, client, tmp_path, monkeypatch) -> None:
        """A valid request should create a run and return 201."""
        from core.run import RunManager
        import core.worker_api as wa

        # Mock the queue enqueue to return a fake job
        mock_job = MagicMock()
        mock_job.id = "job-123"
        wa.pipeline_queue.enqueue.return_value = mock_job

        # Use a tmp-dir RunManager so no real files pollute the project
        tmp_mgr = RunManager(base_dir=tmp_path, use_supabase=False)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            resp = client.post("/api/runs/start", json={
                "pipeline": "vaccine",
                "user_id": "usr-1",
                "team_id": "team-2",
            })

        assert resp.status_code == 201
        data = resp.json()
        assert "run_id" in data
        assert data["status"] == "created"
        assert data["job_id"] == "job-123"

    def test_queue_full_returns_429(self, client, monkeypatch) -> None:
        """When 10+ jobs are queued, endpoint should reject with 429."""
        import core.worker_api as wa

        wa.pipeline_queue.__len__.return_value = 10

        resp = client.post("/api/runs/start", json={
            "pipeline": "vaccine",
        })
        assert resp.status_code == 429

    def test_invalid_preset_returns_422(self, client) -> None:
        resp = client.post("/api/runs/start", json={
            "pipeline": "vaccine",
            "preset": "../../../etc/passwd",
        })
        assert resp.status_code == 422

    def test_invalid_parameter_key_returns_422(self, client) -> None:
        resp = client.post("/api/runs/start", json={
            "pipeline": "vaccine",
            "parameters": {"bad-key!": "value"},
        })
        assert resp.status_code == 422

    def test_path_traversal_in_parameter_value_returns_422(self, client) -> None:
        resp = client.post("/api/runs/start", json={
            "pipeline": "vaccine",
            "parameters": {"config": "../../secrets"},
        })
        assert resp.status_code == 422


# ---------------------------------------------------------------------------
# GET /api/queue/status
# ---------------------------------------------------------------------------


class TestQueueStatusEndpoint:
    """GET /api/queue/status returns queue metrics."""

    def test_queue_status_returns_counts(self, client, monkeypatch) -> None:
        import core.worker_api as wa

        wa.pipeline_queue.__len__.return_value = 3

        # Mock rq.Worker.all() to return fake workers
        mock_worker_busy = MagicMock()
        mock_worker_busy.get_state.return_value = "busy"
        mock_worker_idle = MagicMock()
        mock_worker_idle.get_state.return_value = "idle"

        with patch("core.worker_api.Worker") as MockWorker:
            MockWorker.all.return_value = [mock_worker_busy, mock_worker_idle]
            resp = client.get("/api/queue/status")

        assert resp.status_code == 200
        data = resp.json()
        assert data["queued"] == 3
        assert data["active_workers"] == 2
        assert data["busy_workers"] == 1


# ---------------------------------------------------------------------------
# StartRunRequest model validation
# ---------------------------------------------------------------------------


class TestStartRunRequestModel:
    """Pydantic validation on the StartRunRequest model."""

    def test_user_id_and_team_id_optional(self) -> None:
        from core.worker_api import StartRunRequest

        req = StartRunRequest(pipeline="vaccine")
        assert req.user_id is None
        assert req.team_id is None

    def test_user_id_and_team_id_accepted(self) -> None:
        from core.worker_api import StartRunRequest

        req = StartRunRequest(
            pipeline="vaccine",
            user_id="usr-abc",
            team_id="team-xyz",
        )
        assert req.user_id == "usr-abc"
        assert req.team_id == "team-xyz"

    def test_invalid_pipeline_rejected(self) -> None:
        from pydantic import ValidationError

        from core.worker_api import StartRunRequest

        with pytest.raises(ValidationError):
            StartRunRequest(pipeline="invalid_pipeline")

    def test_valid_preset_accepted(self) -> None:
        from core.worker_api import StartRunRequest

        req = StartRunRequest(pipeline="vaccine", preset="leishmania_infantum")
        assert req.preset == "leishmania_infantum"

    def test_invalid_preset_rejected(self) -> None:
        from pydantic import ValidationError

        from core.worker_api import StartRunRequest

        with pytest.raises(ValidationError):
            StartRunRequest(pipeline="vaccine", preset="bad preset!")

    def test_notes_max_length(self) -> None:
        from core.worker_api import StartRunRequest

        req = StartRunRequest(pipeline="vaccine", notes="x" * 500)
        assert len(req.notes) == 500

    def test_default_parameters_empty(self) -> None:
        from core.worker_api import StartRunRequest

        req = StartRunRequest(pipeline="vaccine")
        assert req.parameters == {}
        assert req.tags == []
        assert req.notes == ""
