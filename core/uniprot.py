"""
UniProt REST API client for fetching protein sequences.

Retrieves amino-acid sequences for validated Leishmania vaccine antigens
using the UniProt REST API (https://rest.uniprot.org).  Supports both
direct accession lookup and free-text search.
"""

from __future__ import annotations

import time

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_BASE_URL: str = "https://rest.uniprot.org/uniprotkb"

REQUEST_TIMEOUT: int = 30  # seconds
API_DELAY_SECONDS: float = 1.0  # polite delay between requests

# Leishmania infantum organism ID in UniProt taxonomy.
LEISHMANIA_INFANTUM_ORGANISM_ID: int = 5671

# Verified UniProt accessions for validated antigens.
# ``None`` means the accession is unknown and must be resolved via search.
VALIDATED_ACCESSIONS: dict[str, str | None] = {
    "A2": "A4HZU7",
    "LACK": "J9XRQ3",
    "HSP70": "A0A6L0XI79",  # HSP70 for L. infantum
    "HSP83": None,  # needs search
    "KMP-11": None,  # needs search
    "LiHyp1": None,  # needs search
    "LiESP_Q": None,  # needs search
    "LBSap_antigens": None,  # composite, skip
    "Lutzomyia_longipalpis_proteins": None,  # vector protein, skip
}

# Antigens to skip during sequence fetching (composite or non-Leishmania).
SKIP_ANTIGENS: set[str] = {"LBSap_antigens", "Lutzomyia_longipalpis_proteins"}

logger = get_logger("uniprot")

# ---------------------------------------------------------------------------
# FASTA fetching by accession
# ---------------------------------------------------------------------------


def fetch_fasta_by_accession(accession: str) -> str | None:
    """Fetch a protein sequence from UniProt by accession number.

    Requests the FASTA representation and parses out the bare amino-acid
    sequence (header line removed, newlines stripped).

    Args:
        accession: UniProt accession (e.g. ``"A4HZU7"``).

    Returns:
        Amino-acid sequence string, or ``None`` on failure.
    """
    url = f"{UNIPROT_BASE_URL}/{accession}.fasta"
    logger.info("Fetching FASTA for accession %s ...", accession)

    try:
        response = requests.get(url, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("Failed to fetch %s: %s", accession, exc)
        return None

    return _parse_fasta_sequence(response.text, accession)


def _parse_fasta_sequence(fasta_text: str, label: str) -> str | None:
    """Extract the amino-acid sequence from raw FASTA text.

    Args:
        fasta_text: Full FASTA content (header + sequence lines).
        label: Human-readable label for log messages.

    Returns:
        Concatenated sequence string, or ``None`` if parsing fails.
    """
    lines = fasta_text.strip().splitlines()
    sequence_lines = [line.strip() for line in lines if line and not line.startswith(">")]

    if not sequence_lines:
        logger.warning("Empty sequence returned for %s.", label)
        return None

    sequence = "".join(sequence_lines)
    logger.info("Got %d-aa sequence for %s.", len(sequence), label)
    return sequence


# ---------------------------------------------------------------------------
# Free-text search
# ---------------------------------------------------------------------------


def search_uniprot(
    query: str,
    organism_id: int = LEISHMANIA_INFANTUM_ORGANISM_ID,
) -> str | None:
    """Search UniProt for a protein by name and return the top hit sequence.

    Constructs a query combining the free-text *query* with an organism
    filter and returns the amino-acid sequence of the first result.

    Args:
        query: Protein / gene name to search for.
        organism_id: UniProt taxonomy ID to restrict results (default:
                     *Leishmania infantum*, 5671).

    Returns:
        Amino-acid sequence of the top hit, or ``None`` if no results.
    """
    search_url = f"{UNIPROT_BASE_URL}/search"
    params = {
        "query": f"{query} AND organism_id:{organism_id}",
        "format": "json",
        "size": "1",
    }

    logger.info("Searching UniProt: %s (organism %d) ...", query, organism_id)

    try:
        response = requests.get(search_url, params=params, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("UniProt search failed for '%s': %s", query, exc)
        return None

    return _extract_sequence_from_json(response.json(), query)


def _extract_sequence_from_json(data: dict, label: str) -> str | None:
    """Pull the amino-acid sequence from a UniProt JSON search response.

    Args:
        data: Parsed JSON dict from the UniProt search endpoint.
        label: Human-readable label for log messages.

    Returns:
        Sequence string, or ``None`` if extraction fails.
    """
    results = data.get("results", [])
    if not results:
        logger.warning("No UniProt results for '%s'.", label)
        return None

    top_hit = results[0]
    sequence_obj = top_hit.get("sequence", {})
    sequence = sequence_obj.get("value")

    if not sequence:
        logger.warning("No sequence in top UniProt result for '%s'.", label)
        return None

    accession = top_hit.get("primaryAccession", "unknown")
    logger.info(
        "Found %s -> %s (%d aa).",
        label,
        accession,
        len(sequence),
    )
    return sequence


# ---------------------------------------------------------------------------
# Bulk fetch for all validated antigens
# ---------------------------------------------------------------------------


def fetch_all_validated_sequences() -> dict[str, str]:
    """Fetch sequences for all validated antigens from UniProt.

    Iterates over ``VALIDATED_ACCESSIONS``, using the known accession when
    available and falling back to a name-based search otherwise.  Antigens
    in ``SKIP_ANTIGENS`` are silently skipped.

    Returns:
        Dictionary mapping ``gene_id`` to its amino-acid sequence.
        Only antigens whose sequences were successfully retrieved are
        included.
    """
    sequences: dict[str, str] = {}

    for gene_id, accession in VALIDATED_ACCESSIONS.items():
        if gene_id in SKIP_ANTIGENS:
            logger.info("Skipping %s (composite / vector protein).", gene_id)
            continue

        # Polite delay between API calls.
        if sequences:
            time.sleep(API_DELAY_SECONDS)

        if accession:
            seq = fetch_fasta_by_accession(accession)
        else:
            seq = search_uniprot(gene_id)

        if seq:
            sequences[gene_id] = seq
        else:
            logger.warning("Could not retrieve sequence for %s.", gene_id)

    logger.info(
        "Retrieved %d/%d validated antigen sequence(s).",
        len(sequences),
        len(VALIDATED_ACCESSIONS) - len(SKIP_ANTIGENS),
    )
    return sequences
