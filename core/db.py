"""
Supabase database helpers for the Marley pipeline.

Handles connection setup and CRUD operations on the ``candidates`` table.
All functions log through ``core.logger`` and raise on unrecoverable errors
after logging context.
"""

from __future__ import annotations

import os
from typing import Any

from dotenv import load_dotenv
from supabase import Client, create_client

from core.logger import get_logger
from core.models import Candidate

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

load_dotenv()

SUPABASE_URL: str = os.getenv("SUPABASE_URL", "")
SUPABASE_KEY: str = os.getenv("SUPABASE_KEY", "")
CANDIDATES_TABLE: str = "candidates"

logger = get_logger("db")

# ---------------------------------------------------------------------------
# Client initialisation
# ---------------------------------------------------------------------------

_client: Client | None = None


def _get_client() -> Client:
    """Return a lazily-initialised Supabase client.

    Raises:
        EnvironmentError: If ``SUPABASE_URL`` or ``SUPABASE_KEY`` are missing.
    """
    global _client  # noqa: PLW0603

    if _client is not None:
        return _client

    if not SUPABASE_URL or not SUPABASE_KEY:
        msg = (
            "SUPABASE_URL and SUPABASE_KEY must be set in the environment. "
            "Add them to your .env file."
        )
        logger.error(msg)
        raise EnvironmentError(msg)

    try:
        _client = create_client(SUPABASE_URL, SUPABASE_KEY)
        logger.info("Supabase client initialised successfully.")
    except Exception:
        logger.exception("Failed to create Supabase client.")
        raise

    return _client


# ---------------------------------------------------------------------------
# CRUD helpers
# ---------------------------------------------------------------------------


def upsert_candidate(candidate: Candidate) -> dict[str, Any]:
    """Insert or update a candidate row keyed by ``gene_id``.

    Args:
        candidate: The ``Candidate`` to persist.

    Returns:
        The Supabase response data for the upserted row.

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()
    payload = candidate.to_dict()

    try:
        response = (
            client.table(CANDIDATES_TABLE)
            .upsert(payload, on_conflict="gene_id")
            .execute()
        )
        logger.info("Upserted candidate %s.", candidate.gene_id)
        return response.data
    except Exception:
        logger.exception("Failed to upsert candidate %s.", candidate.gene_id)
        raise


def get_all_candidates() -> list[Candidate]:
    """Fetch every row from the candidates table.

    Returns:
        List of ``Candidate`` instances (may be empty).

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()

    try:
        response = client.table(CANDIDATES_TABLE).select("*").execute()
        candidates = [Candidate.from_dict(row) for row in response.data]
        logger.info("Fetched %d candidate(s).", len(candidates))
        return candidates
    except Exception:
        logger.exception("Failed to fetch candidates.")
        raise


def update_candidate(gene_id: str, data: dict[str, Any]) -> dict[str, Any]:
    """Apply a partial update to an existing candidate.

    Args:
        gene_id: The ``gene_id`` of the row to update.
        data: A dictionary of column-value pairs to set.

    Returns:
        The Supabase response data for the updated row.

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()

    try:
        response = (
            client.table(CANDIDATES_TABLE)
            .update(data)
            .eq("gene_id", gene_id)
            .execute()
        )
        logger.info("Updated candidate %s.", gene_id)
        return response.data
    except Exception:
        logger.exception("Failed to update candidate %s.", gene_id)
        raise
