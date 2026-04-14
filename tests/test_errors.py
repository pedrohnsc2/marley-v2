"""Tests for core.errors — error classification system."""

from __future__ import annotations

import subprocess

import pytest

from core.errors import ErrorInfo, classify_error
from core.run import StageRecord


# ---------------------------------------------------------------------------
# ErrorInfo data model
# ---------------------------------------------------------------------------

class TestErrorInfoRoundTrip:
    """ErrorInfo.to_dict / from_dict round-trip."""

    def test_round_trip_full(self) -> None:
        info = ErrorInfo(
            category="network",
            code="network.timeout",
            message="Connection timed out",
            stage_id="01_fetch",
            suggestion="Check your internet connection.",
        )
        d = info.to_dict()
        restored = ErrorInfo.from_dict(d)
        assert restored == info

    def test_round_trip_minimal(self) -> None:
        info = ErrorInfo(category="internal", code="internal.unknown", message="oops")
        d = info.to_dict()
        restored = ErrorInfo.from_dict(d)
        assert restored == info

    def test_to_dict_keys(self) -> None:
        info = ErrorInfo(category="data", code="data.missing_input", message="nope")
        d = info.to_dict()
        assert set(d.keys()) == {"category", "code", "message", "stage_id", "suggestion"}

    def test_frozen(self) -> None:
        info = ErrorInfo(category="data", code="data.missing_input", message="x")
        with pytest.raises(AttributeError):
            info.category = "network"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Dependency errors
# ---------------------------------------------------------------------------

class TestDependencyClassification:
    def test_module_not_found(self) -> None:
        exc = ModuleNotFoundError("No module named 'biopython'")
        info = classify_error(exc)
        assert info.code == "dependency.module_missing"
        assert info.category == "dependency"
        assert info.suggestion  # non-empty

    def test_tool_missing_command_not_found(self) -> None:
        exc = FileNotFoundError("[Errno 2] command not found: 'blastp'")
        info = classify_error(exc)
        assert info.code == "dependency.tool_missing"

    def test_tool_missing_bin_path(self) -> None:
        exc = FileNotFoundError("No such file: /usr/local/bin/muscle")
        info = classify_error(exc)
        assert info.code == "dependency.tool_missing"


# ---------------------------------------------------------------------------
# Network errors
# ---------------------------------------------------------------------------

class TestNetworkClassification:
    def test_timeout_pattern(self) -> None:
        exc = Exception("Connection timed out after 30s")
        info = classify_error(exc)
        assert info.code == "network.timeout"

    def test_read_timeout(self) -> None:
        exc = Exception("ReadTimeout: read timed out")
        info = classify_error(exc)
        assert info.code == "network.timeout"

    def test_connection_refused(self) -> None:
        exc = ConnectionRefusedError("[Errno 111] Connection refused")
        info = classify_error(exc)
        assert info.code == "network.connection_refused"

    def test_connection_reset(self) -> None:
        exc = ConnectionResetError("Connection reset by peer")
        info = classify_error(exc)
        assert info.code == "network.connection_refused"

    def test_ssl_error(self) -> None:
        exc = Exception("SSLError: CERTIFICATE_VERIFY_FAILED")
        info = classify_error(exc)
        assert info.code == "network.ssl_error"

    def test_dns_failure(self) -> None:
        exc = Exception("Name or service not known")
        info = classify_error(exc)
        assert info.code == "network.dns_failure"

    def test_dns_getaddrinfo(self) -> None:
        exc = Exception("getaddrinfo failed for api.uniprot.org")
        info = classify_error(exc)
        assert info.code == "network.dns_failure"

    def test_http_500(self) -> None:
        exc = Exception("HTTP 500 Server Error")
        info = classify_error(exc)
        assert info.code == "network.api_error"

    def test_http_404(self) -> None:
        exc = Exception("HTTP 404 Not Found")
        info = classify_error(exc)
        assert info.code == "network.api_error"

    def test_http_error_class(self) -> None:
        exc = Exception("HTTPError: 502 Bad Gateway")
        info = classify_error(exc)
        assert info.code == "network.api_error"


# ---------------------------------------------------------------------------
# Data errors
# ---------------------------------------------------------------------------

class TestDataClassification:
    def test_file_not_found_general(self) -> None:
        exc = FileNotFoundError("No such file: results/proteins.fasta")
        info = classify_error(exc)
        assert info.code == "data.missing_input"

    def test_invalid_format_value_error(self) -> None:
        exc = ValueError("invalid format: expected FASTA header")
        info = classify_error(exc)
        assert info.code == "data.invalid_format"

    def test_invalid_format_parse(self) -> None:
        exc = KeyError("parse error in column 3")
        info = classify_error(exc)
        assert info.code == "data.invalid_format"

    def test_unexpected_column(self) -> None:
        exc = ValueError("unexpected column 'foo' in input CSV")
        info = classify_error(exc)
        assert info.code == "data.invalid_format"

    def test_empty_result(self) -> None:
        exc = Exception("empty result set after filtering")
        info = classify_error(exc)
        assert info.code == "data.empty_result"

    def test_no_candidates(self) -> None:
        exc = Exception("no candidates passed the threshold")
        info = classify_error(exc)
        assert info.code == "data.empty_result"

    def test_zero_sequences(self) -> None:
        exc = Exception("Found 0 sequences in input FASTA")
        info = classify_error(exc)
        assert info.code == "data.empty_result"

    def test_zero_records(self) -> None:
        exc = Exception("zero records returned from query")
        info = classify_error(exc)
        assert info.code == "data.empty_result"


# ---------------------------------------------------------------------------
# Compute errors
# ---------------------------------------------------------------------------

class TestComputeClassification:
    def test_memory_error(self) -> None:
        exc = MemoryError("unable to allocate array")
        info = classify_error(exc)
        assert info.code == "compute.out_of_memory"

    def test_overflow(self) -> None:
        exc = OverflowError("math range error")
        info = classify_error(exc)
        assert info.code == "compute.numerical_error"

    def test_zero_division(self) -> None:
        exc = ZeroDivisionError("division by zero")
        info = classify_error(exc)
        assert info.code == "compute.numerical_error"

    def test_floating_point(self) -> None:
        exc = FloatingPointError("underflow")
        info = classify_error(exc)
        assert info.code == "compute.numerical_error"

    def test_subprocess_returncode(self) -> None:
        exc = Exception("Process exited with returncode 1")
        info = classify_error(exc)
        assert info.code == "compute.subprocess_failed"

    def test_subprocess_exit_code(self) -> None:
        exc = Exception("exit code 137 (OOM killed)")
        info = classify_error(exc)
        assert info.code == "compute.subprocess_failed"

    def test_subprocess_nonzero(self) -> None:
        exc = Exception("non-zero exit status")
        info = classify_error(exc)
        assert info.code == "compute.subprocess_failed"

    def test_called_process_error(self) -> None:
        exc = subprocess.CalledProcessError(1, "blastp")
        info = classify_error(exc)
        assert info.code == "compute.subprocess_failed"


# ---------------------------------------------------------------------------
# Config errors
# ---------------------------------------------------------------------------

class TestConfigClassification:
    def test_preset_not_found(self) -> None:
        exc = Exception("preset 'exotic_species' not found")
        info = classify_error(exc)
        assert info.code == "config.preset_not_found"

    def test_invalid_parameter_unexpected_keyword(self) -> None:
        exc = TypeError("got an unexpected keyword argument 'foobar'")
        info = classify_error(exc)
        assert info.code == "config.invalid_parameter"

    def test_invalid_parameter_missing_argument(self) -> None:
        exc = TypeError("missing 1 required positional argument: 'config'")
        info = classify_error(exc)
        assert info.code == "config.invalid_parameter"


# ---------------------------------------------------------------------------
# Permission errors
# ---------------------------------------------------------------------------

class TestPermissionClassification:
    def test_permission_denied(self) -> None:
        exc = PermissionError("[Errno 13] Permission denied: '/root/data'")
        info = classify_error(exc)
        assert info.code == "permission.file_access"

    def test_eacces(self) -> None:
        exc = OSError("[Errno 13] EACCES: /var/log/pipeline.log")
        info = classify_error(exc)
        assert info.code == "permission.file_access"

    def test_api_401_unauthorized(self) -> None:
        exc = Exception("Unauthorized: invalid token")
        info = classify_error(exc)
        assert info.code == "permission.api_auth"

    def test_api_403_forbidden(self) -> None:
        exc = Exception("Forbidden — check your credentials")
        info = classify_error(exc)
        assert info.code == "permission.api_auth"

    def test_http_401_matches_api_error_first(self) -> None:
        """HTTP 401 matches the broader HTTP error rule first due to rule
        ordering. This is by design — the HTTP pattern fires before the
        auth keyword pattern."""
        exc = Exception("HTTP 401 Unauthorized")
        info = classify_error(exc)
        assert info.code == "network.api_error"

    def test_api_key_invalid(self) -> None:
        exc = Exception("API key is invalid or expired")
        info = classify_error(exc)
        assert info.code == "permission.api_auth"


# ---------------------------------------------------------------------------
# Internal errors
# ---------------------------------------------------------------------------

class TestInternalClassification:
    def test_assertion_error(self) -> None:
        exc = AssertionError("invariant violated")
        info = classify_error(exc)
        assert info.code == "internal.assertion"
        assert info.category == "internal"

    def test_fallback_unknown(self) -> None:
        exc = RuntimeError("something completely unexpected")
        info = classify_error(exc)
        assert info.code == "internal.unknown"
        assert info.category == "internal"
        assert info.suggestion  # non-empty


# ---------------------------------------------------------------------------
# Stage ID propagation
# ---------------------------------------------------------------------------

class TestStageIdPropagation:
    def test_stage_id_set(self) -> None:
        exc = FileNotFoundError("missing proteins.fasta")
        info = classify_error(exc, stage_id="03_filter")
        assert info.stage_id == "03_filter"
        assert info.code == "data.missing_input"

    def test_stage_id_default_empty(self) -> None:
        exc = RuntimeError("boom")
        info = classify_error(exc)
        assert info.stage_id == ""


# ---------------------------------------------------------------------------
# Rule ordering (first-match semantics)
# ---------------------------------------------------------------------------

class TestRuleOrdering:
    def test_file_not_found_with_bin_matches_tool_missing_not_data(self) -> None:
        """FileNotFoundError + 'bin/' should match dependency.tool_missing,
        not data.missing_input."""
        exc = FileNotFoundError("/usr/local/bin/clustalw not found")
        info = classify_error(exc)
        assert info.code == "dependency.tool_missing"

    def test_timeout_pattern_beats_generic(self) -> None:
        """A ValueError that mentions 'timed out' should match network.timeout
        (pattern-based) rather than data rules."""
        exc = ValueError("Request timed out after 60 seconds")
        info = classify_error(exc)
        assert info.code == "network.timeout"

    def test_permission_error_without_pattern_is_not_matched(self) -> None:
        """OSError without 'permission denied' or 'EACCES' should fall through
        to internal.unknown (or another matching rule)."""
        exc = OSError("Some other OS error")
        info = classify_error(exc)
        assert info.code != "permission.file_access"


# ---------------------------------------------------------------------------
# Backwards compatibility: StageRecord without error_info
# ---------------------------------------------------------------------------

class TestStageRecordBackwardsCompat:
    def test_from_dict_without_error_info(self) -> None:
        """Old metadata.json files without error_info should still load."""
        old_data = {
            "stage_id": "01_fetch",
            "name": "Fetch Genome",
            "status": "failed",
            "started_at": "2026-04-12T10:00:00+00:00",
            "completed_at": "2026-04-12T10:01:00+00:00",
            "duration_s": 60.0,
            "error": "Connection refused",
            "output_file": None,
            "key_metrics": {"duration_s": 60.0},
        }
        stage = StageRecord.from_dict(old_data)
        assert stage.error_info is None
        assert stage.error == "Connection refused"
        assert stage.stage_id == "01_fetch"

    def test_from_dict_with_error_info(self) -> None:
        """New metadata.json files with error_info should load correctly."""
        new_data = {
            "stage_id": "02_filter",
            "name": "Filter Proteins",
            "status": "failed",
            "error": "Connection timed out",
            "error_info": {
                "category": "network",
                "code": "network.timeout",
                "message": "Connection timed out",
                "stage_id": "02_filter",
                "suggestion": "Check your internet connection.",
            },
        }
        stage = StageRecord.from_dict(new_data)
        assert stage.error_info is not None
        assert stage.error_info["code"] == "network.timeout"

    def test_to_dict_includes_error_info(self) -> None:
        stage = StageRecord(
            stage_id="03_score",
            name="Score",
            error_info={"category": "data", "code": "data.empty_result",
                        "message": "empty", "stage_id": "03_score",
                        "suggestion": "Relax filters."},
        )
        d = stage.to_dict()
        assert "error_info" in d
        assert d["error_info"]["code"] == "data.empty_result"

    def test_to_dict_error_info_none(self) -> None:
        stage = StageRecord(stage_id="04_report", name="Report")
        d = stage.to_dict()
        assert "error_info" in d
        assert d["error_info"] is None
