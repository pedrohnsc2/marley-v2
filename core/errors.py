"""Error classification for pipeline execution.

Classifies exceptions into structured ErrorInfo objects with category,
code, human-readable suggestion, and the raw message. Used by the
pipeline runner to produce researcher-friendly error diagnostics.

Categories:
    data        — missing/malformed input files or empty results
    network     — connectivity, DNS, SSL, HTTP errors
    compute     — OOM, numerical failures, subprocess crashes
    config      — invalid parameters, missing presets
    dependency  — missing Python modules or external tools
    permission  — file system or API authentication errors
    internal    — assertion failures and uncategorised errors
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any


# ---------------------------------------------------------------------------
# ErrorInfo data model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ErrorInfo:
    """Structured representation of a classified pipeline error."""

    category: str        # e.g. "network"
    code: str            # e.g. "network.timeout"
    message: str         # raw exception message
    stage_id: str = ""
    suggestion: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "category": self.category,
            "code": self.code,
            "message": self.message,
            "stage_id": self.stage_id,
            "suggestion": self.suggestion,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ErrorInfo:
        return cls(
            category=data["category"],
            code=data["code"],
            message=data.get("message", ""),
            stage_id=data.get("stage_id", ""),
            suggestion=data.get("suggestion", ""),
        )


# ---------------------------------------------------------------------------
# Suggestions per code
# ---------------------------------------------------------------------------

_SUGGESTIONS: dict[str, str] = {
    "data.missing_input": (
        "A required input file was not found. Check that previous pipeline "
        "stages completed successfully and that file paths are correct."
    ),
    "data.invalid_format": (
        "The input data has an unexpected format. Verify the file is not "
        "corrupted and matches the expected schema (FASTA, CSV, JSON, etc.)."
    ),
    "data.empty_result": (
        "A processing step produced zero results. Try relaxing filter "
        "thresholds (e.g. conservation_threshold, e_value) or verify "
        "that the input data contains valid entries."
    ),
    "network.timeout": (
        "A network request timed out. Check your internet connection and "
        "retry. External services (UniProt, NCBI) may be temporarily slow."
    ),
    "network.connection_refused": (
        "Could not connect to a remote service. Verify the server is "
        "reachable and that no firewall is blocking the connection."
    ),
    "network.dns_failure": (
        "DNS resolution failed. Check your internet connection and DNS "
        "settings. The remote host may be temporarily unreachable."
    ),
    "network.api_error": (
        "An external API returned an HTTP error. Check the service status "
        "page and verify that your query parameters are valid."
    ),
    "network.ssl_error": (
        "SSL/TLS verification failed. Ensure your system certificates are "
        "up to date. If using a corporate proxy, configure SSL_CERT_FILE."
    ),
    "compute.out_of_memory": (
        "The process ran out of memory. Try reducing batch sizes, using "
        "a smaller dataset, or running on a machine with more RAM."
    ),
    "compute.numerical_error": (
        "A numerical computation failed (overflow, division by zero, etc.). "
        "Check input values for extremes or missing data that could cause "
        "mathematical errors."
    ),
    "compute.timeout": (
        "A computation exceeded the time limit. Try reducing the dataset "
        "size or increasing the timeout setting."
    ),
    "compute.subprocess_failed": (
        "An external tool exited with an error. Check that the tool is "
        "installed correctly and review its output log for details."
    ),
    "config.invalid_parameter": (
        "A configuration parameter is invalid. Review the parameter names "
        "and types in your preset or command-line overrides."
    ),
    "config.missing_parameter": (
        "A required configuration parameter is missing. Check the pipeline "
        "documentation for required parameters."
    ),
    "config.preset_not_found": (
        "The specified preset was not found. Run --list-presets to see "
        "available presets for this pipeline."
    ),
    "dependency.module_missing": (
        "A required Python module is not installed. Run "
        "'pip install -r requirements.txt' to install dependencies."
    ),
    "dependency.tool_missing": (
        "A required external tool is not installed or not on PATH. "
        "Check the installation guide for required system dependencies."
    ),
    "permission.file_access": (
        "File access was denied. Check file permissions and ensure the "
        "output directory is writable."
    ),
    "permission.api_auth": (
        "API authentication failed. Verify your API key or credentials "
        "are set correctly in the environment or .env file."
    ),
    "internal.assertion": (
        "An internal assertion failed. This is likely a bug — please "
        "report it with the full error traceback."
    ),
    "internal.unknown": (
        "An unexpected error occurred. Review the full traceback for "
        "details and consider reporting this as a bug."
    ),
}


# ---------------------------------------------------------------------------
# Classification rules
# ---------------------------------------------------------------------------

# Each rule is: (exception_types | None, regex_pattern | None, code)
# - exception_types: tuple of types to match via isinstance, or None for any
# - regex_pattern: compiled regex to match against str(exc), or None for any
# - code: the error code to assign on match
#
# First match wins. The list order matters.

_RULES: list[tuple[tuple[type, ...] | None, re.Pattern[str] | None, str]] = [
    # -- dependency --
    (
        (ModuleNotFoundError,),
        None,
        "dependency.module_missing",
    ),
    (
        (FileNotFoundError,),
        re.compile(r"(command not found|bin/)", re.IGNORECASE),
        "dependency.tool_missing",
    ),

    # -- network --
    (
        None,
        re.compile(r"(timed?\s*out|TimeoutError|ReadTimeout)", re.IGNORECASE),
        "network.timeout",
    ),
    (
        (ConnectionError, ConnectionRefusedError, ConnectionResetError),
        None,
        "network.connection_refused",
    ),
    (
        None,
        re.compile(r"(SSLError|CERTIFICATE_VERIFY_FAILED)", re.IGNORECASE),
        "network.ssl_error",
    ),
    (
        None,
        re.compile(
            r"(Name or service not known|getaddrinfo failed|DNS)",
            re.IGNORECASE,
        ),
        "network.dns_failure",
    ),
    (
        None,
        re.compile(r"(HTTP\s*(4\d\d|5\d\d)|HTTPError)", re.IGNORECASE),
        "network.api_error",
    ),

    # -- data --
    (
        (FileNotFoundError,),
        None,
        "data.missing_input",
    ),
    (
        (ValueError, KeyError),
        re.compile(
            r"(invalid.*format|parse|unexpected.*column)",
            re.IGNORECASE,
        ),
        "data.invalid_format",
    ),
    (
        None,
        re.compile(
            r"(empty.*result|no.*candidates|0 sequences|zero.*records)",
            re.IGNORECASE,
        ),
        "data.empty_result",
    ),

    # -- compute --
    (
        (MemoryError,),
        None,
        "compute.out_of_memory",
    ),
    (
        (OverflowError, ZeroDivisionError, FloatingPointError),
        None,
        "compute.numerical_error",
    ),
    (
        None,
        re.compile(
            r"(returncode|exit\s*code|non-?zero|CalledProcessError)",
            re.IGNORECASE,
        ),
        "compute.subprocess_failed",
    ),

    # -- config --
    (
        None,
        re.compile(r"(preset.*not found)", re.IGNORECASE),
        "config.preset_not_found",
    ),
    (
        (TypeError,),
        re.compile(
            r"(unexpected keyword|missing.*argument)",
            re.IGNORECASE,
        ),
        "config.invalid_parameter",
    ),

    # -- permission --
    (
        (PermissionError, OSError),
        re.compile(r"(permission denied|EACCES)", re.IGNORECASE),
        "permission.file_access",
    ),
    (
        None,
        re.compile(r"(401|403|Unauthorized|Forbidden|API.key)", re.IGNORECASE),
        "permission.api_auth",
    ),

    # -- internal --
    (
        (AssertionError,),
        None,
        "internal.assertion",
    ),
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def classify_error(exc: BaseException, stage_id: str = "") -> ErrorInfo:
    """Classify an exception into a structured ErrorInfo.

    Walks the rule chain in order; first match wins.  Falls back to
    ``internal.unknown`` if no rule matches.

    Args:
        exc: The exception to classify.
        stage_id: Optional pipeline stage where the error occurred.

    Returns:
        A frozen ErrorInfo dataclass with category, code, suggestion,
        and the raw exception message.
    """
    msg = str(exc)

    for exc_types, pattern, code in _RULES:
        # Check exception type constraint
        if exc_types is not None and not isinstance(exc, exc_types):
            continue
        # Check regex pattern constraint
        if pattern is not None and not pattern.search(msg):
            continue
        # Match found
        category = code.split(".")[0]
        return ErrorInfo(
            category=category,
            code=code,
            message=msg,
            stage_id=stage_id,
            suggestion=_SUGGESTIONS.get(code, ""),
        )

    # Fallback
    return ErrorInfo(
        category="internal",
        code="internal.unknown",
        message=msg,
        stage_id=stage_id,
        suggestion=_SUGGESTIONS["internal.unknown"],
    )
