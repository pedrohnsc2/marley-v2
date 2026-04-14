"""Tests for core.run_cli — create_pipeline_run function."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from core.run import PipelineRun, RunManager


# ---------------------------------------------------------------------------
# create_pipeline_run
# ---------------------------------------------------------------------------


class TestCreatePipelineRun:
    """create_pipeline_run() creates a PipelineRun via RunManager."""

    def _make_tmp_run_manager(self, tmp_path: Path) -> RunManager:
        """Create a RunManager rooted in tmp_path with Supabase disabled."""
        return RunManager(base_dir=tmp_path, use_supabase=False)

    def test_create_with_user_and_team(self, tmp_path: Path) -> None:
        """user_id and team_id are passed through to the created run."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="vaccine",
                parameters={"organism": "leishmania_infantum"},
                user_id="usr-111",
                team_id="team-222",
            )

        assert run.user_id == "usr-111"
        assert run.team_id == "team-222"
        assert run.pipeline == "vaccine"

    def test_create_without_user_backward_compat(self, tmp_path: Path) -> None:
        """Omitting user_id/team_id should produce None values."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="vaccine",
                parameters={},
            )

        assert run.user_id is None
        assert run.team_id is None

    def test_run_persisted_to_disk(self, tmp_path: Path) -> None:
        """The run should be written as metadata.json on disk."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="drug",
                parameters={"top_n": 5},
                user_id="u1",
                team_id="t1",
                notes="disk test",
            )

        meta_path = tmp_mgr.runs_dir / run.run_id / "metadata.json"
        assert meta_path.exists()

        with open(meta_path) as fh:
            meta = json.load(fh)

        assert meta["pipeline"] == "drug"
        assert meta["user_id"] == "u1"
        assert meta["team_id"] == "t1"
        assert meta["notes"] == "disk test"
        assert meta["parameters"]["top_n"] == 5

    def test_dry_run_tag_added(self, tmp_path: Path) -> None:
        """When dry_run=True, 'dry-run' tag should be appended."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="vaccine",
                parameters={},
                dry_run=True,
            )

        assert "dry-run" in run.tags

    def test_dry_run_false_no_tag(self, tmp_path: Path) -> None:
        """When dry_run=False, 'dry-run' tag should not appear."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="vaccine",
                parameters={},
                dry_run=False,
            )

        assert "dry-run" not in run.tags

    def test_custom_tags_preserved(self, tmp_path: Path) -> None:
        """Tags passed in should be present on the resulting run."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="vaccine",
                parameters={},
                tags=["nightly", "ci"],
            )

        assert "nightly" in run.tags
        assert "ci" in run.tags

    def test_status_is_created(self, tmp_path: Path) -> None:
        """Newly created run should have status 'created'."""
        tmp_mgr = self._make_tmp_run_manager(tmp_path)

        with patch("core.run_cli.RunManager", return_value=tmp_mgr):
            from core.run_cli import create_pipeline_run

            run = create_pipeline_run(
                pipeline_id="vaccine",
                parameters={},
            )

        assert run.status == "created"
