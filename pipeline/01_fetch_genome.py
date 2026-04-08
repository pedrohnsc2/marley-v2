"""Fetch the *Leishmania infantum* JPCM5 proteome from TriTrypDB.

Downloads all protein sequences for the organism and saves them as a
single FASTA file at ``data/raw/proteins.fasta``.  Implements retry with
exponential back-off and a graceful fallback with manual-download
instructions when the remote API is unreachable.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Final

import requests
from tqdm import tqdm

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TRITRYPDB_BASE_URL: Final[str] = "https://tritrypdb.org/tritrypdb/service"
ORGANISM: Final[str] = "linfantumJPCM5"
OUTPUT_DIR: Final[str] = "data/raw"
OUTPUT_FILE: Final[str] = "data/raw/proteins.fasta"
MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]  # exponential backoff in seconds

# Alternative direct download URL used when the REST search endpoint fails.
FALLBACK_DOWNLOAD_URL: Final[str] = (
    "https://tritrypdb.org/common/downloads/release-68/"
    "LinfantumJPCM5/fasta/data/"
    "TriTrypDB-68_LinfantumJPCM5_AnnotatedProteins.fasta"
)

# HTTP settings
REQUEST_TIMEOUT: Final[int] = 120  # seconds
DOWNLOAD_CHUNK_SIZE: Final[int] = 8192  # bytes
USER_AGENT: Final[str] = "Marley-Pipeline/1.0 (reverse-vaccinology research)"

MANUAL_DOWNLOAD_URL: Final[str] = (
    "https://tritrypdb.org/tritrypdb/app/downloads/"
    "release-current/LinfantumJPCM5/fasta/data/"
)

logger = get_logger("fetch_genome")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Download helpers (private)
# ---------------------------------------------------------------------------


def _build_search_url() -> str:
    """Return the TriTrypDB REST endpoint for fetching gene records."""
    return (
        f"{TRITRYPDB_BASE_URL}/record-types/transcript/searches/"
        "GenesByTaxonGene/reports/standard"
    )


def _common_headers() -> dict[str, str]:
    """Return HTTP headers shared across all requests."""
    return {
        "User-Agent": USER_AGENT,
        "Accept": "application/json",
    }


def _try_api_download(session: requests.Session) -> str | None:
    """Attempt to download protein FASTA via the TriTrypDB REST search API.

    The endpoint accepts a POST with organism filter parameters and returns
    gene records. We request the FASTA report attachment.

    Args:
        session: A ``requests.Session`` pre-configured with headers.

    Returns:
        The FASTA content as a string, or ``None`` on failure.
    """
    url = _build_search_url()
    payload = {
        "searchConfig": {
            "parameters": {
                "organism": [f"[{ORGANISM}]"],
            },
            "wdkWeight": 10,
        },
        "reportConfig": {
            "attributes": ["primary_key", "product"],
            "tables": [],
            "attachmentType": "plain",
        },
    }

    logger.info("Trying TriTrypDB search API: %s", url)
    try:
        resp = session.post(url, json=payload, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        # If we get JSON back, we need the FASTA report instead.
        # The search API may not directly return FASTA; fall through to
        # the direct download in that case.
        content_type = resp.headers.get("Content-Type", "")
        if "fasta" in content_type or resp.text.lstrip().startswith(">"):
            logger.info("Received FASTA from search API.")
            return resp.text
        logger.warning(
            "Search API returned non-FASTA content (%s). "
            "Falling back to direct download.",
            content_type,
        )
    except requests.RequestException as exc:
        logger.warning("Search API request failed: %s", exc)

    return None


def _try_direct_download(session: requests.Session, dest: Path) -> bool:
    """Stream-download the proteome FASTA from the TriTrypDB file server.

    Uses ``tqdm`` to display a progress bar when the server provides a
    ``Content-Length`` header.

    Args:
        session: A ``requests.Session`` pre-configured with headers.
        dest: Destination ``Path`` for the downloaded file.

    Returns:
        ``True`` on success, ``False`` on failure.
    """
    logger.info("Trying direct download: %s", FALLBACK_DOWNLOAD_URL)
    try:
        resp = session.get(
            FALLBACK_DOWNLOAD_URL,
            timeout=REQUEST_TIMEOUT,
            stream=True,
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("Direct download request failed: %s", exc)
        return False

    total_size = int(resp.headers.get("Content-Length", 0))
    progress = tqdm(
        total=total_size or None,
        unit="B",
        unit_scale=True,
        desc="Downloading proteome",
        disable=total_size == 0,
    )

    try:
        with open(dest, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=DOWNLOAD_CHUNK_SIZE):
                if chunk:
                    fh.write(chunk)
                    progress.update(len(chunk))
    except Exception:
        logger.exception("Error while writing downloaded data.")
        # Clean up partial file
        if dest.exists():
            dest.unlink()
        return False
    finally:
        progress.close()

    # Quick sanity check: the file should start with ">"
    with open(dest, "r", encoding="utf-8", errors="replace") as fh:
        first_char = fh.read(1)
    if first_char != ">":
        logger.error(
            "Downloaded file does not look like FASTA (first char: %r). "
            "Removing corrupt file.",
            first_char,
        )
        dest.unlink()
        return False

    return True


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def fetch_genome(force: bool = False) -> str:
    """Download the *L. infantum* JPCM5 proteome and save it locally.

    The function first checks whether the output file already exists.
    Unless *force* is ``True``, an existing file is reused.  Two download
    strategies are attempted with retry and exponential back-off:

    1. **TriTrypDB search REST API** -- structured query for the organism.
    2. **Direct file download** -- pre-built FASTA from the TriTrypDB
       release server.

    If both strategies fail after all retries, a ``RuntimeError`` is raised
    with instructions for manual download.

    Args:
        force: Re-download even if the local file already exists.

    Returns:
        Path to the downloaded FASTA file (``data/raw/proteins.fasta``).

    Raises:
        RuntimeError: If the proteome could not be downloaded after all
            retry attempts.
    """
    output_path = Path(OUTPUT_FILE)

    # --- Check for existing file -------------------------------------------
    if output_path.exists() and not force:
        seq_count = count_sequences(str(output_path))
        logger.warning(
            "Output file already exists: %s (%d sequences). "
            "Use force=True to re-download.",
            output_path,
            seq_count,
        )
        return str(output_path)

    # --- Ensure output directory exists ------------------------------------
    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Output directory ready: %s", output_dir)

    # --- Attempt download with retries -------------------------------------
    session = requests.Session()
    session.headers.update(_common_headers())

    for attempt in range(1, MAX_RETRIES + 1):
        logger.info("Download attempt %d/%d", attempt, MAX_RETRIES)

        # Strategy 1: REST search API
        fasta_text = _try_api_download(session)
        if fasta_text is not None:
            output_path.write_text(fasta_text, encoding="utf-8")
            seq_count = count_sequences(str(output_path))
            logger.info(
                "Proteome saved via API: %s (%d sequences)",
                output_path,
                seq_count,
            )
            return str(output_path)

        # Strategy 2: direct file download
        if _try_direct_download(session, output_path):
            seq_count = count_sequences(str(output_path))
            logger.info(
                "Proteome saved via direct download: %s (%d sequences)",
                output_path,
                seq_count,
            )
            return str(output_path)

        # Both strategies failed -- back off before retrying
        if attempt < MAX_RETRIES:
            delay = RETRY_DELAYS[attempt - 1]
            logger.warning(
                "Both download strategies failed. "
                "Retrying in %d seconds...",
                delay,
            )
            time.sleep(delay)

    # --- All retries exhausted ---------------------------------------------
    msg = (
        "Failed to download the L. infantum JPCM5 proteome after "
        f"{MAX_RETRIES} attempts.\n\n"
        "You can download it manually:\n"
        f"  1. Visit {MANUAL_DOWNLOAD_URL}\n"
        "  2. Download the *AnnotatedProteins.fasta* file\n"
        f"  3. Save it as {OUTPUT_FILE}\n"
        "  4. Re-run this pipeline step.\n"
    )
    logger.error(msg)
    raise RuntimeError(msg)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """Run the fetch step as a standalone script.

    Accepts an optional ``--force`` flag to re-download even when the
    output file already exists.
    """
    force = "--force" in sys.argv
    try:
        path = fetch_genome(force=force)
        logger.info("Fetch complete: %s", path)
    except RuntimeError:
        logger.error("Fetch step failed. See messages above for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
