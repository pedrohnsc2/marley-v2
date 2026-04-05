"""Fetch AlphaFold structures and prepare receptor PDBQT files for docking.

Downloads AlphaFold PDB models for top-ranked drug targets, then converts
them to PDBQT format suitable for molecular docking with AutoDock Vina.

Usage:
    python -m drug_targets.06_fetch_structures
    python -m drug_targets.06_fetch_structures --top-n 10
    python -m drug_targets.06_fetch_structures --force
    python -m drug_targets.06_fetch_structures --dry-run
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path
from typing import Final

import requests

from core.db import get_all_drug_targets
from core.logger import get_logger
from core.models import DrugTarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STRUCTURES_DIR: Final[Path] = Path("data/structures")

ALPHAFOLD_PDB_URL: Final[str] = (
    "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
)

MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]

KNOWN_UNIPROT_IDS: Final[dict[str, str]] = {
    "TryS": "Q4QJ30",
    "TryR": "Q4Q457",
    "ADL": "A4I7R2",
    "GMPS": "A4IBM8",
    "SMT": "A4I4D2",
    "6PGDH": "Q18L04",
    "XPRT": "E9AGS8",
    "HGPRT": "Q1PC47",
    "APRT": "A4HY89",
    "SDM": "A4ICL7",
    "SQS": "A4I3T7",
    "TPI": "A4I0S4",
    "PGK": "A4HSP0",
    "AK": "A4I9S5",
}

# Regex to extract a UniProt accession from an AlphaFold URL.
_UNIPROT_FROM_URL_RE: Final[re.Pattern[str]] = re.compile(
    r"[/=]([A-Z0-9]{6,10})(?:[/\-?&#]|$)", re.IGNORECASE
)

logger = get_logger("fetch_structures")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _resolve_uniprot_id(target: DrugTarget) -> str | None:
    """Determine the UniProt accession for a drug target.

    Resolution order:
    1. Extract from ``alphafold_url`` if present and parseable.
    2. Look up ``gene_name`` in ``KNOWN_UNIPROT_IDS``.

    Args:
        target: The drug target to resolve.

    Returns:
        A UniProt accession string, or ``None`` if unresolvable.
    """
    if target.alphafold_url:
        match = _UNIPROT_FROM_URL_RE.search(target.alphafold_url)
        if match:
            return match.group(1)

    return KNOWN_UNIPROT_IDS.get(target.gene_name)


def _download_with_retry(url: str, dest: Path) -> bool:
    """Download a file from *url* to *dest* with exponential backoff.

    Args:
        url: The URL to download.
        dest: Local file path to save the downloaded content.

    Returns:
        ``True`` if the download succeeded, ``False`` otherwise.
    """
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(response.content)
            logger.info("Downloaded %s -> %s", url, dest)
            return True
        except requests.RequestException as exc:
            if attempt < MAX_RETRIES - 1:
                delay = RETRY_DELAYS[attempt]
                logger.warning(
                    "Attempt %d/%d failed for %s: %s. Retrying in %ds.",
                    attempt + 1,
                    MAX_RETRIES,
                    url,
                    exc,
                    delay,
                )
                time.sleep(delay)
            else:
                logger.error(
                    "All %d attempts failed for %s: %s",
                    MAX_RETRIES,
                    url,
                    exc,
                )
    return False


def _convert_pdb_to_pdbqt_meeko(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Convert PDB to PDBQT using meeko (primary method).

    Args:
        pdb_path: Input PDB file path.
        pdbqt_path: Output PDBQT file path.

    Returns:
        ``True`` if conversion succeeded, ``False`` otherwise.
    """
    try:
        from meeko import MoleculePreparation, PDBMoleculeSetup  # type: ignore[import-untyped]

        setup = PDBMoleculeSetup.from_pdb_file(str(pdb_path))
        preparator = MoleculePreparation()
        preparator.prepare(setup)

        pdbqt_string = preparator.write_pdbqt_string()
        pdbqt_path.write_text(pdbqt_string, encoding="utf-8")
        logger.info("Converted %s -> %s via meeko.", pdb_path.name, pdbqt_path.name)
        return True
    except ImportError:
        logger.debug("meeko not available, trying fallback methods.")
        return False
    except Exception as exc:
        logger.warning("meeko conversion failed for %s: %s", pdb_path.name, exc)
        return False


def _convert_pdb_to_pdbqt_obabel(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Convert PDB to PDBQT using Open Babel CLI (fallback method).

    Args:
        pdb_path: Input PDB file path.
        pdbqt_path: Output PDBQT file path.

    Returns:
        ``True`` if conversion succeeded, ``False`` otherwise.
    """
    try:
        result = subprocess.run(
            [
                "obabel",
                "-ipdb", str(pdb_path),
                "-opdbqt",
                "-O", str(pdbqt_path),
                "-xrh",
            ],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        if result.returncode == 0 and pdbqt_path.exists():
            logger.info(
                "Converted %s -> %s via obabel.", pdb_path.name, pdbqt_path.name
            )
            return True

        logger.warning(
            "obabel returned code %d for %s: %s",
            result.returncode,
            pdb_path.name,
            result.stderr.strip(),
        )
        return False
    except FileNotFoundError:
        logger.debug("obabel not found on PATH, trying fallback.")
        return False
    except subprocess.TimeoutExpired:
        logger.warning("obabel timed out for %s.", pdb_path.name)
        return False
    except Exception as exc:
        logger.warning("obabel conversion failed for %s: %s", pdb_path.name, exc)
        return False


def _convert_pdb_to_pdbqt_simple(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Convert PDB to PDBQT with a minimal Python-only approach (last resort).

    Reads ATOM/HETATM records from the PDB, assigns Gasteiger-like neutral
    charges (0.000) and AutoDock atom types based on element, then writes a
    basic PDBQT file.  This is functional enough for initial docking setup,
    though production runs should use meeko or Open Babel.

    Args:
        pdb_path: Input PDB file path.
        pdbqt_path: Output PDBQT file path.

    Returns:
        ``True`` if conversion succeeded, ``False`` otherwise.
    """
    # Map common PDB elements to AutoDock atom types.
    ad_type_map: dict[str, str] = {
        "C": "C",
        "N": "N",
        "O": "OA",
        "S": "SA",
        "H": "HD",
        "P": "P",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
        "FE": "Fe",
        "ZN": "Zn",
        "MG": "Mg",
        "CA": "Ca",
        "MN": "Mn",
    }

    try:
        pdb_lines = pdb_path.read_text(encoding="utf-8").splitlines()
        pdbqt_lines: list[str] = [
            "REMARK  PDBQT generated by Marley pipeline (simple converter)",
            "REMARK  For production use, install meeko or Open Babel",
        ]

        for line in pdb_lines:
            record = line[:6].strip()
            if record not in ("ATOM", "HETATM"):
                # Pass through non-coordinate records (MODEL, END, etc.).
                if record in ("MODEL", "ENDMDL", "END", "TER"):
                    pdbqt_lines.append(line)
                continue

            # Extract element from columns 77-78 (standard PDB format).
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                # Fall back to atom name (columns 13-16), strip digits.
                atom_name = line[12:16].strip()
                element = "".join(c for c in atom_name if c.isalpha())[:2].upper()
                if len(element) > 1 and element not in ad_type_map:
                    element = element[0]

            ad_type = ad_type_map.get(element, "C")

            # Pad line to 54 chars (coordinates end), then rebuild.
            padded = line.ljust(78)
            # PDBQT format: columns 1-54 same as PDB, then occupancy(55-60),
            # temp factor(61-66), partial charge(67-76), atom type(77-78).
            coord_part = padded[:54]
            occupancy = padded[54:60] if padded[54:60].strip() else "  1.00"
            temp_factor = padded[60:66] if padded[60:66].strip() else "  0.00"

            pdbqt_line = (
                f"{coord_part}{occupancy}{temp_factor}    +0.000 {ad_type:<2s}"
            )
            pdbqt_lines.append(pdbqt_line)

        pdbqt_path.write_text("\n".join(pdbqt_lines) + "\n", encoding="utf-8")
        logger.info(
            "Converted %s -> %s via simple converter.", pdb_path.name, pdbqt_path.name
        )
        return True
    except Exception as exc:
        logger.error(
            "Simple PDBQT conversion failed for %s: %s", pdb_path.name, exc
        )
        return False


def _convert_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Convert a PDB file to PDBQT, trying multiple methods.

    Attempts meeko first, then Open Babel CLI, then a simple Python-based
    converter as a last resort.

    Args:
        pdb_path: Input PDB file path.
        pdbqt_path: Output PDBQT file path.

    Returns:
        ``True`` if any method succeeded, ``False`` if all failed.
    """
    if _convert_pdb_to_pdbqt_meeko(pdb_path, pdbqt_path):
        return True
    if _convert_pdb_to_pdbqt_obabel(pdb_path, pdbqt_path):
        return True
    return _convert_pdb_to_pdbqt_simple(pdb_path, pdbqt_path)


def _get_top_targets(top_n: int) -> list[DrugTarget]:
    """Load the top-N drug targets sorted by druggability score.

    Only targets with a resolvable UniProt ID are included.

    Args:
        top_n: Maximum number of targets to return.

    Returns:
        List of ``DrugTarget`` instances with resolvable UniProt IDs,
        sorted by druggability_score descending.
    """
    all_targets = get_all_drug_targets()

    # Sort by druggability_score descending.
    all_targets.sort(key=lambda t: t.druggability_score, reverse=True)

    # Filter to targets with resolvable UniProt IDs.
    resolvable: list[DrugTarget] = []
    for target in all_targets:
        if _resolve_uniprot_id(target) is not None:
            resolvable.append(target)
            if len(resolvable) >= top_n:
                break

    return resolvable


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def fetch_and_prepare_structures(
    top_n: int = 5,
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Download AlphaFold PDB structures and convert to PDBQT for docking.

    Fetches structures for the top-N drug targets (sorted by druggability
    score) from the AlphaFold database. Each PDB is then converted to
    PDBQT format using meeko, Open Babel, or a simple fallback converter.

    Args:
        top_n: Number of top-scoring drug targets to process.
        force: Re-download and re-convert even if files already exist.
        dry_run: Log planned actions without downloading or converting.

    Returns:
        Path to the structures output directory.

    Raises:
        RuntimeError: If no targets with resolvable UniProt IDs are found.
    """
    structures_dir = STRUCTURES_DIR
    structures_dir.mkdir(parents=True, exist_ok=True)

    targets = _get_top_targets(top_n)

    if not targets:
        raise RuntimeError(
            "No drug targets with resolvable UniProt IDs found. "
            "Run 04_druggability.py first to populate the database."
        )

    logger.info(
        "Selected %d target(s) for structure retrieval (top_n=%d).",
        len(targets),
        top_n,
    )

    downloaded_count = 0

    for target in targets:
        uniprot_id = _resolve_uniprot_id(target)
        if uniprot_id is None:
            # Should not happen after filtering, but guard anyway.
            continue

        gene_name = target.gene_name
        pdb_path = structures_dir / f"{gene_name}_{uniprot_id}.pdb"
        pdbqt_path = structures_dir / f"{gene_name}_{uniprot_id}_receptor.pdbqt"

        download_url = ALPHAFOLD_PDB_URL.format(uniprot_id=uniprot_id)

        if dry_run:
            logger.info(
                "[DRY RUN] Would download %s -> %s", download_url, pdb_path
            )
            logger.info(
                "[DRY RUN] Would convert %s -> %s", pdb_path.name, pdbqt_path.name
            )
            downloaded_count += 1
            continue

        # Download PDB.
        if pdb_path.exists() and not force:
            logger.info(
                "PDB already exists: %s. Use --force to re-download.", pdb_path
            )
        else:
            if not _download_with_retry(download_url, pdb_path):
                logger.error(
                    "Skipping %s (%s): PDB download failed.",
                    gene_name,
                    uniprot_id,
                )
                continue

        # Convert to PDBQT.
        if pdbqt_path.exists() and not force:
            logger.info(
                "PDBQT already exists: %s. Use --force to re-convert.",
                pdbqt_path,
            )
        else:
            if not _convert_pdb_to_pdbqt(pdb_path, pdbqt_path):
                logger.error(
                    "Failed to convert %s to PDBQT. PDB is still available at %s.",
                    gene_name,
                    pdb_path,
                )
                # PDB was downloaded successfully; count it even if conversion failed.

        downloaded_count += 1

    action = "Would download and prepare" if dry_run else "Downloaded and prepared"
    logger.info(
        "%s %d structures for docking.", action, downloaded_count
    )

    return str(structures_dir)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Fetch AlphaFold structures and prepare PDBQT for docking.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top drug targets to fetch structures for (default: 5).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download and re-convert even if files already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would be done without downloading or converting.",
    )
    args = parser.parse_args()

    result = fetch_and_prepare_structures(
        top_n=args.top_n,
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Complete: %s", result)
