"""Fetch mRNA transcriptome sequences for codon usage and entropy analysis.

Downloads *Leishmania infantum* JPCM5 and human reference mRNA sequences
from the UniProt REST API.  Saves FASTA files to ``data/raw/`` for
downstream codon usage comparison and Shannon entropy profiling.

Usage:
    python -m rna_entropy.01_fetch_transcriptome
    python -m rna_entropy.01_fetch_transcriptome --force
    python -m rna_entropy.01_fetch_transcriptome --dry-run
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_SEARCH_URL: Final[str] = "https://rest.uniprot.org/uniprotkb/search"

LINF_QUERY: Final[str] = "organism_id:5671"
LINF_SIZE: Final[int] = 500
LINF_OUTPUT: Final[str] = "data/raw/transcriptome_linf.fasta"

HUMAN_QUERY: Final[str] = "organism_id:9606+keyword:reference+proteome"
HUMAN_SIZE: Final[int] = 100
HUMAN_OUTPUT: Final[str] = "data/raw/transcriptome_human.fasta"

MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 120
USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (rna-entropy-analysis)"

logger = get_logger("fetch_transcriptome")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _fetch_with_retry(url: str, params: dict | None = None) -> requests.Response:
    """Perform an HTTP GET with exponential backoff retry.

    Args:
        url: The URL to fetch.
        params: Optional query parameters.

    Returns:
        The successful ``requests.Response``.

    Raises:
        requests.RequestException: After all retries are exhausted.
    """
    headers = {"User-Agent": USER_AGENT, "Accept": "text/plain"}

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(
                url, params=params, headers=headers, timeout=REQUEST_TIMEOUT
            )
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            if attempt < MAX_RETRIES - 1:
                delay = RETRY_DELAYS[attempt]
                logger.warning(
                    "Request failed (attempt %d/%d): %s — retrying in %ds",
                    attempt + 1,
                    MAX_RETRIES,
                    exc,
                    delay,
                )
                time.sleep(delay)
            else:
                logger.error("All %d attempts failed: %s", MAX_RETRIES, exc)
                raise

    raise requests.RequestException("Unexpected retry loop exit")


def count_sequences(fasta_path: str) -> int:
    """Count the number of sequences in a FASTA file.

    Each sequence starts with a ``>`` character at the beginning of a line.

    Args:
        fasta_path: Absolute or relative path to the FASTA file.

    Returns:
        Number of sequences found.

    Raises:
        FileNotFoundError: If *fasta_path* does not exist.
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    count = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def _download_fasta(query: str, size: int, output_path: Path) -> int:
    """Download FASTA sequences from UniProt and save to disk.

    Args:
        query: UniProt search query string.
        size: Maximum number of records to fetch.
        output_path: Destination path for the FASTA file.

    Returns:
        Number of sequences in the downloaded file.

    Raises:
        RuntimeError: If the download fails or produces no sequences.
    """
    params = {
        "query": query,
        "format": "fasta",
        "size": str(size),
    }

    logger.info("Fetching from UniProt: query=%s, size=%d", query, size)
    response = _fetch_with_retry(UNIPROT_SEARCH_URL, params=params)

    content = response.text
    if not content.strip() or not content.lstrip().startswith(">"):
        msg = f"UniProt returned no valid FASTA data for query: {query}"
        logger.error(msg)
        raise RuntimeError(msg)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content, encoding="utf-8")

    seq_count = count_sequences(str(output_path))
    logger.info("Saved %d sequences to %s", seq_count, output_path)
    return seq_count


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def fetch_transcriptome(force: bool = False, dry_run: bool = False) -> str:
    """Download L. infantum and human mRNA transcriptome sequences.

    Fetches protein/mRNA sequences from UniProt for both organisms and
    saves them as FASTA files for downstream codon usage and entropy
    analysis.

    Args:
        force: Re-download even if output files already exist.
        dry_run: Log what would be downloaded without making API calls.

    Returns:
        Path to the L. infantum FASTA file.

    Raises:
        RuntimeError: If downloads fail after all retry attempts.
    """
    linf_path = Path(LINF_OUTPUT)
    human_path = Path(HUMAN_OUTPUT)

    # --- Dry run mode ---------------------------------------------------------
    if dry_run:
        logger.info("[DRY RUN] Would download L. infantum transcripts:")
        logger.info("[DRY RUN]   URL: %s?query=%s&size=%d", UNIPROT_SEARCH_URL, LINF_QUERY, LINF_SIZE)
        logger.info("[DRY RUN]   Output: %s", linf_path)
        logger.info("[DRY RUN] Would download human reference transcripts:")
        logger.info("[DRY RUN]   URL: %s?query=%s&size=%d", UNIPROT_SEARCH_URL, HUMAN_QUERY, HUMAN_SIZE)
        logger.info("[DRY RUN]   Output: %s", human_path)
        return str(linf_path)

    # --- Check for existing files ---------------------------------------------
    if linf_path.exists() and human_path.exists() and not force:
        linf_count = count_sequences(str(linf_path))
        human_count = count_sequences(str(human_path))
        logger.info(
            "Transcriptome files already exist: %d L. infantum and %d human "
            "transcripts. Use --force to re-download.",
            linf_count,
            human_count,
        )
        return str(linf_path)

    # --- Download L. infantum transcriptome -----------------------------------
    linf_count = _download_fasta(LINF_QUERY, LINF_SIZE, linf_path)

    # --- Download human reference transcriptome -------------------------------
    human_count = _download_fasta(HUMAN_QUERY, HUMAN_SIZE, human_path)

    # --- Summary --------------------------------------------------------------
    logger.info(
        "%d L. infantum and %d human transcripts downloaded.",
        linf_count,
        human_count,
    )

    return str(linf_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Fetch L. infantum and human transcriptome from UniProt.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if FASTA files already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would be done without making API calls.",
    )
    args = parser.parse_args()

    result = fetch_transcriptome(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
