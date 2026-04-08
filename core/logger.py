"""
Standardized logging for the Marley reverse vaccinology pipeline.

Provides a consistent log format across all modules:
    [MARLEY][YYYY-MM-DD HH:MM:SS] LEVEL - message
"""

import logging
import sys

LOG_PREFIX = "MARLEY"
LOG_FORMAT = f"[{LOG_PREFIX}][%(asctime)s] %(levelname)s - %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
DEFAULT_LOG_LEVEL = logging.INFO


class MarleyFormatter(logging.Formatter):
    """Custom formatter that applies the Marley log format."""

    def __init__(self) -> None:
        super().__init__(fmt=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)


def get_logger(name: str) -> logging.Logger:
    """Return a configured logger with the Marley format.

    If the logger already has handlers (i.e. it was previously configured),
    it is returned as-is to avoid duplicate output.

    Args:
        name: Logical name for the logger, typically the module name.

    Returns:
        A ``logging.Logger`` instance ready for use.
    """
    logger = logging.getLogger(f"{LOG_PREFIX}.{name}")

    if not logger.handlers:
        logger.setLevel(DEFAULT_LOG_LEVEL)

        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(DEFAULT_LOG_LEVEL)
        handler.setFormatter(MarleyFormatter())

        logger.addHandler(handler)
        logger.propagate = False

    return logger
