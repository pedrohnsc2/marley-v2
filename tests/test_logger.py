"""Tests for the Marley logging module."""

from __future__ import annotations

import logging
import re

from core.logger import LOG_PREFIX, get_logger


def test_get_logger_returns_logger() -> None:
    """get_logger should return a logging.Logger instance."""
    logger = get_logger("test_module")
    assert isinstance(logger, logging.Logger)


def test_logger_format(capfd: object) -> None:
    """Log output should match the [MARLEY][timestamp] LEVEL - message format."""
    # Use a unique name to avoid handler reuse from other tests.
    logger = get_logger("format_check")
    logger.info("hello world")

    captured = capfd.readouterr()  # type: ignore[attr-defined]
    # Expected pattern: [MARLEY][YYYY-MM-DD HH:MM:SS] INFO - hello world
    pattern = r"\[MARLEY\]\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\] INFO - hello world"
    assert re.search(pattern, captured.out), f"Output did not match expected format: {captured.out!r}"


def test_logger_levels(capfd: object) -> None:
    """INFO, WARNING, and ERROR messages should all appear in output."""
    logger = get_logger("levels_check")

    logger.info("info message")
    logger.warning("warning message")
    logger.error("error message")

    captured = capfd.readouterr()  # type: ignore[attr-defined]
    assert "INFO - info message" in captured.out
    assert "WARNING - warning message" in captured.out
    assert "ERROR - error message" in captured.out
