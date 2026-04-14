"""Tests for core.launcher — input validation."""

from __future__ import annotations

import pytest

from core.launcher import _validate_input


# ---------------------------------------------------------------------------
# _validate_input
# ---------------------------------------------------------------------------


class TestValidateInputPipeline:
    """Validation of the 'pipeline' field in launcher input."""

    def test_valid_pipelines_accepted(self) -> None:
        """All registered pipelines should pass validation."""
        for pipeline in ("vaccine", "drug", "docking", "rna", "aso_math"):
            _validate_input({"pipeline": pipeline})  # should not raise

    def test_invalid_pipeline_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid pipeline"):
            _validate_input({"pipeline": "nonexistent_pipeline"})

    def test_empty_pipeline_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid pipeline"):
            _validate_input({"pipeline": ""})

    def test_missing_pipeline_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid pipeline"):
            _validate_input({})

    def test_none_pipeline_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid pipeline"):
            _validate_input({"pipeline": None})


class TestValidateInputPreset:
    """Validation of the optional 'preset' field in launcher input."""

    def test_valid_preset_accepted(self) -> None:
        _validate_input({"pipeline": "vaccine", "preset": "leishmania_infantum"})

    def test_alphanumeric_underscore_accepted(self) -> None:
        _validate_input({"pipeline": "vaccine", "preset": "trypanosoma_cruzi"})

    def test_empty_preset_accepted(self) -> None:
        """Empty string preset should pass (treated as no preset)."""
        _validate_input({"pipeline": "vaccine", "preset": ""})

    def test_preset_with_special_chars_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid preset name"):
            _validate_input({"pipeline": "vaccine", "preset": "../../etc/passwd"})

    def test_preset_with_spaces_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid preset name"):
            _validate_input({"pipeline": "vaccine", "preset": "bad preset"})

    def test_preset_with_semicolon_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid preset name"):
            _validate_input({"pipeline": "vaccine", "preset": "preset;rm -rf"})

    def test_preset_with_dots_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid preset name"):
            _validate_input({"pipeline": "vaccine", "preset": "foo.bar"})

    def test_no_preset_key_accepted(self) -> None:
        """Missing preset key should not raise (defaults to empty)."""
        _validate_input({"pipeline": "vaccine"})
