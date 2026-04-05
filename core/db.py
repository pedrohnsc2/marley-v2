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
from core.models import Candidate, DrugTarget, Epitope, VaccineConstruct

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

load_dotenv()

SUPABASE_URL: str = os.getenv("SUPABASE_URL", "")
SUPABASE_KEY: str = os.getenv("SUPABASE_KEY", "")
CANDIDATES_TABLE: str = "candidates"
CONSTRUCTS_TABLE: str = "vaccine_constructs"
EPITOPES_TABLE: str = "construct_epitopes"
DRUG_TARGETS_TABLE: str = "drug_targets"

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


# ---------------------------------------------------------------------------
# Vaccine construct helpers (Module 06)
# ---------------------------------------------------------------------------


def upsert_construct(construct: VaccineConstruct) -> dict[str, Any]:
    """Insert or update a vaccine construct row keyed by ``construct_id``.

    Args:
        construct: The ``VaccineConstruct`` to persist.

    Returns:
        The Supabase response data for the upserted row.

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()
    payload = construct.to_dict()

    # Remove empty created_at so Supabase uses its DEFAULT now()
    if not payload.get("created_at"):
        payload.pop("created_at", None)

    try:
        response = (
            client.table(CONSTRUCTS_TABLE)
            .upsert(payload, on_conflict="construct_id")
            .execute()
        )
        logger.info("Upserted construct %s.", construct.construct_id)
        return response.data
    except Exception:
        logger.exception("Failed to upsert construct %s.", construct.construct_id)
        raise


def get_construct(construct_id: str) -> VaccineConstruct | None:
    """Fetch a vaccine construct by its ID.

    Args:
        construct_id: The unique construct identifier.

    Returns:
        A ``VaccineConstruct`` instance, or ``None`` if not found.

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()

    try:
        response = (
            client.table(CONSTRUCTS_TABLE)
            .select("*")
            .eq("construct_id", construct_id)
            .execute()
        )
        if not response.data:
            logger.info("Construct %s not found.", construct_id)
            return None

        construct = VaccineConstruct.from_dict(response.data[0])
        logger.info("Fetched construct %s.", construct_id)
        return construct
    except Exception:
        logger.exception("Failed to fetch construct %s.", construct_id)
        raise


def upsert_construct_epitopes(
    construct_id: str, epitopes: list[Epitope]
) -> None:
    """Insert epitopes for a construct, replacing any existing rows.

    Existing epitopes for the given *construct_id* are deleted first to
    avoid duplicates on re-runs.

    Args:
        construct_id: The construct these epitopes belong to.
        epitopes: List of ``Epitope`` instances to persist.

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()

    try:
        # Remove previous epitopes for this construct
        client.table(EPITOPES_TABLE).delete().eq(
            "construct_id", construct_id
        ).execute()

        if not epitopes:
            logger.info("Cleared epitopes for construct %s (none to insert).", construct_id)
            return

        rows = [
            {
                "construct_id": construct_id,
                "peptide": ep.sequence,
                "epitope_type": ep.epitope_type,
                "source_gene_id": ep.source_gene_id,
                "source_gene_name": ep.source_gene_name,
                "allele": ep.allele,
                "ic50": ep.ic50,
                "position_in_construct": ep.start_position,
            }
            for ep in epitopes
        ]

        client.table(EPITOPES_TABLE).insert(rows).execute()
        logger.info(
            "Inserted %d epitope(s) for construct %s.",
            len(epitopes),
            construct_id,
        )
    except Exception:
        logger.exception(
            "Failed to upsert epitopes for construct %s.", construct_id
        )
        raise


# ---------------------------------------------------------------------------
# Drug target helpers (Module v2)
# ---------------------------------------------------------------------------


def upsert_drug_target(target: DrugTarget) -> dict[str, Any]:
    """Insert or update a drug target row keyed by ``gene_id``.

    Args:
        target: The ``DrugTarget`` to persist.

    Returns:
        The Supabase response data for the upserted row.

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()
    payload = target.to_dict()

    try:
        response = (
            client.table(DRUG_TARGETS_TABLE)
            .upsert(payload, on_conflict="gene_id")
            .execute()
        )
        logger.info("Upserted drug target %s.", target.gene_id)
        return response.data
    except Exception:
        logger.exception("Failed to upsert drug target %s.", target.gene_id)
        raise


def get_all_drug_targets() -> list[DrugTarget]:
    """Fetch every row from the drug_targets table.

    Returns:
        List of ``DrugTarget`` instances (may be empty).

    Raises:
        Exception: On any Supabase communication failure.
    """
    client = _get_client()

    try:
        response = client.table(DRUG_TARGETS_TABLE).select("*").execute()
        targets = [DrugTarget.from_dict(row) for row in response.data]
        logger.info("Fetched %d drug target(s).", len(targets))
        return targets
    except Exception:
        logger.exception("Failed to fetch drug targets.")
        raise


def update_drug_target(gene_id: str, data: dict[str, Any]) -> dict[str, Any]:
    """Apply a partial update to an existing drug target.

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
            client.table(DRUG_TARGETS_TABLE)
            .update(data)
            .eq("gene_id", gene_id)
            .execute()
        )
        logger.info("Updated drug target %s.", gene_id)
        return response.data
    except Exception:
        logger.exception("Failed to update drug target %s.", gene_id)
        raise
