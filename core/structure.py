"""3D structure prediction and visualization for vaccine constructs.

Provides an ESMFold API client for protein structure prediction, region
mapping for construct annotation, pLDDT quality extraction, and script
generators for PyMOL and UCSF ChimeraX visualization.

Usage:
    from core.structure import predict_structure, map_construct_regions

    pdb_path = predict_structure(protein_seq, "results/construct/predicted.pdb")
    regions = map_construct_regions(protein_seq, "tPA", "L7L12", epitopes)
"""

from __future__ import annotations

import importlib
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from core.logger import get_logger

logger = get_logger("structure")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ESMFOLD_URL: str = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_SEQUENCE_LENGTH: int = 800
REQUEST_TIMEOUT: int = 120
MAX_RETRIES: int = 2
RETRY_DELAY: int = 5

COLOR_SCHEME: dict[str, str] = {
    "signal_peptide": "#3B82F6",
    "adjuvant": "#22C55E",
    "linker": "#9CA3AF",
    "ctl_epitope": "#EF4444",
    "htl_epitope": "#F59E0B",
}


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class ConstructRegion:
    """A contiguous region within the multi-epitope vaccine construct.

    Attributes:
        name: Unique identifier, e.g. ``"signal_peptide"``, ``"CTL_1"``.
        region_type: Category key into :data:`COLOR_SCHEME`.
        start: 1-based residue start position (inclusive).
        end: 1-based residue end position (inclusive).
        color: Hex color string for visualization.
        label: Human-readable display label.
    """

    name: str
    region_type: str
    start: int
    end: int
    color: str
    label: str


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _load_construct_constants() -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Import SIGNAL_PEPTIDES, ADJUVANTS, and LINKERS from the construct module.

    Uses importlib so that this module does not create a hard import-time
    dependency on the pipeline package.

    Returns:
        Tuple of (signal_peptides, adjuvants, linkers) dictionaries.

    Raises:
        ImportError: If ``pipeline.06_construct`` cannot be imported.
    """
    mod = importlib.import_module("pipeline.06_construct")
    signal_peptides: dict[str, str] = getattr(mod, "SIGNAL_PEPTIDES")
    adjuvants: dict[str, str] = getattr(mod, "ADJUVANTS")
    linkers: dict[str, str] = getattr(mod, "LINKERS")
    return signal_peptides, adjuvants, linkers


def _validate_sequence(protein_seq: str) -> bool:
    """Check that a protein sequence is non-empty and within length limits.

    Args:
        protein_seq: Amino-acid sequence string.

    Returns:
        ``True`` if the sequence passes validation.
    """
    if not protein_seq or not protein_seq.strip():
        logger.error("Protein sequence is empty")
        return False
    if len(protein_seq) > MAX_SEQUENCE_LENGTH:
        logger.error(
            "Sequence length %d exceeds maximum %d residues",
            len(protein_seq),
            MAX_SEQUENCE_LENGTH,
        )
        return False
    return True


# ---------------------------------------------------------------------------
# 1. ESMFold API client
# ---------------------------------------------------------------------------


def predict_structure(
    protein_seq: str,
    output_path: str | Path,
    timeout: int = REQUEST_TIMEOUT,
) -> Path | None:
    """Predict 3D structure via the ESMFold API and save the resulting PDB.

    Sends the raw amino-acid sequence to the ESMFold REST endpoint. On
    transient failure the request is retried up to :data:`MAX_RETRIES`
    times with a :data:`RETRY_DELAY` second backoff between attempts.

    Args:
        protein_seq: Single-letter amino-acid sequence (no FASTA header).
        output_path: Destination file path for the PDB output.
        timeout: HTTP request timeout in seconds.

    Returns:
        :class:`~pathlib.Path` to the saved PDB file on success, or
        ``None`` if prediction failed (non-blocking).
    """
    if not _validate_sequence(protein_seq):
        return None

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    clean_seq = protein_seq.strip()
    attempts = 0

    while attempts <= MAX_RETRIES:
        attempts += 1
        try:
            logger.info(
                "ESMFold request attempt %d/%d (%d residues)",
                attempts,
                MAX_RETRIES + 1,
                len(clean_seq),
            )
            response = requests.post(
                ESMFOLD_URL,
                data=clean_seq,
                headers={"Content-Type": "text/plain"},
                timeout=timeout,
            )
            response.raise_for_status()

            pdb_content = response.text
            if "ATOM" not in pdb_content:
                logger.error("ESMFold response does not contain ATOM records")
                return None

            output_path.write_text(pdb_content, encoding="utf-8")
            logger.info("PDB saved to %s", output_path)
            return output_path

        except requests.exceptions.RequestException as exc:
            logger.warning("ESMFold request failed (attempt %d): %s", attempts, exc)
            if attempts <= MAX_RETRIES:
                logger.info("Retrying in %d seconds...", RETRY_DELAY)
                time.sleep(RETRY_DELAY)

    logger.error(
        "ESMFold prediction failed after %d attempts", MAX_RETRIES + 1
    )
    return None


# ---------------------------------------------------------------------------
# 2. Region mapper
# ---------------------------------------------------------------------------


def map_construct_regions(
    protein_seq: str,
    signal_peptide_name: str,
    adjuvant_name: str,
    epitopes: list[dict],
) -> list[ConstructRegion]:
    """Walk the construct sequence and assign region boundaries.

    The expected construct layout is::

        signal_peptide -> adjuvant -> EAAAK -> [CTL1-AAY-CTL2-AAY-...]
            -> GPGPG -> [HTL1-GPGPG-HTL2-GPGPG-...]

    Region boundaries are computed from the known component lengths
    defined in ``pipeline.06_construct``.

    Args:
        protein_seq: Full protein sequence of the construct.
        signal_peptide_name: Key into ``SIGNAL_PEPTIDES`` (e.g. ``"tPA"``).
        adjuvant_name: Key into ``ADJUVANTS`` (e.g. ``"L7L12"``).
        epitopes: List of epitope dicts from the construct card JSON.
            Each dict must contain ``"sequence"`` and ``"epitope_type"``
            (``"CTL"`` or ``"HTL"``).

    Returns:
        Ordered list of :class:`ConstructRegion` covering the full
        construct.
    """
    try:
        signal_peptides, adjuvants, linkers = _load_construct_constants()
    except ImportError:
        logger.error("Cannot import construct constants from pipeline.06_construct")
        return []

    sp_seq = signal_peptides.get(signal_peptide_name, "")
    adj_seq = adjuvants.get(adjuvant_name, "")
    eaaak = linkers.get("adjuvant_to_epitopes", "EAAAK")
    ctl_linker = linkers.get("ctl", "AAY")
    htl_linker = linkers.get("htl", "GPGPG")

    if not sp_seq:
        logger.warning(
            "Unknown signal peptide '%s'; region mapping may be inaccurate",
            signal_peptide_name,
        )
    if not adj_seq:
        logger.warning(
            "Unknown adjuvant '%s'; region mapping may be inaccurate",
            adjuvant_name,
        )

    regions: list[ConstructRegion] = []
    pos = 1  # 1-based position cursor

    # --- Signal peptide ---
    if sp_seq:
        sp_len = len(sp_seq)
        regions.append(ConstructRegion(
            name="signal_peptide",
            region_type="signal_peptide",
            start=pos,
            end=pos + sp_len - 1,
            color=COLOR_SCHEME["signal_peptide"],
            label=f"Signal peptide ({signal_peptide_name})",
        ))
        pos += sp_len

    # --- Adjuvant ---
    if adj_seq:
        adj_len = len(adj_seq)
        regions.append(ConstructRegion(
            name="adjuvant",
            region_type="adjuvant",
            start=pos,
            end=pos + adj_len - 1,
            color=COLOR_SCHEME["adjuvant"],
            label=f"Adjuvant ({adjuvant_name})",
        ))
        pos += adj_len

    # --- EAAAK linker (adjuvant -> epitopes) ---
    eaaak_len = len(eaaak)
    regions.append(ConstructRegion(
        name="linker_EAAAK",
        region_type="linker",
        start=pos,
        end=pos + eaaak_len - 1,
        color=COLOR_SCHEME["linker"],
        label="Linker (EAAAK)",
    ))
    pos += eaaak_len

    # --- Separate epitopes by type ---
    ctl_epitopes = [e for e in epitopes if e.get("epitope_type") == "CTL"]
    htl_epitopes = [e for e in epitopes if e.get("epitope_type") == "HTL"]

    # --- CTL epitopes joined by AAY ---
    for idx, ep in enumerate(ctl_epitopes):
        ep_seq = ep.get("sequence", "")
        ep_len = len(ep_seq)
        if ep_len == 0:
            continue

        regions.append(ConstructRegion(
            name=f"CTL_{idx + 1}",
            region_type="ctl_epitope",
            start=pos,
            end=pos + ep_len - 1,
            color=COLOR_SCHEME["ctl_epitope"],
            label=f"CTL epitope {idx + 1}",
        ))
        pos += ep_len

        # Add AAY linker between CTL epitopes (not after the last one)
        if idx < len(ctl_epitopes) - 1:
            link_len = len(ctl_linker)
            regions.append(ConstructRegion(
                name=f"linker_AAY_{idx + 1}",
                region_type="linker",
                start=pos,
                end=pos + link_len - 1,
                color=COLOR_SCHEME["linker"],
                label="Linker (AAY)",
            ))
            pos += link_len

    # --- GPGPG linker between CTL and HTL blocks ---
    if ctl_epitopes and htl_epitopes:
        gpgpg_len = len(htl_linker)
        regions.append(ConstructRegion(
            name="linker_GPGPG_bridge",
            region_type="linker",
            start=pos,
            end=pos + gpgpg_len - 1,
            color=COLOR_SCHEME["linker"],
            label="Linker (GPGPG)",
        ))
        pos += gpgpg_len

    # --- HTL epitopes joined by GPGPG ---
    for idx, ep in enumerate(htl_epitopes):
        ep_seq = ep.get("sequence", "")
        ep_len = len(ep_seq)
        if ep_len == 0:
            continue

        regions.append(ConstructRegion(
            name=f"HTL_{idx + 1}",
            region_type="htl_epitope",
            start=pos,
            end=pos + ep_len - 1,
            color=COLOR_SCHEME["htl_epitope"],
            label=f"HTL epitope {idx + 1}",
        ))
        pos += ep_len

        # Add GPGPG linker between HTL epitopes (not after the last one)
        if idx < len(htl_epitopes) - 1:
            gpgpg_len = len(htl_linker)
            regions.append(ConstructRegion(
                name=f"linker_GPGPG_{idx + 1}",
                region_type="linker",
                start=pos,
                end=pos + gpgpg_len - 1,
                color=COLOR_SCHEME["linker"],
                label="Linker (GPGPG)",
            ))
            pos += gpgpg_len

    # Sanity check: mapped length vs. actual sequence length
    mapped_length = pos - 1
    actual_length = len(protein_seq)
    if mapped_length != actual_length:
        logger.warning(
            "Region mapping covers %d residues but sequence has %d; "
            "boundaries may be misaligned",
            mapped_length,
            actual_length,
        )

    return regions


# ---------------------------------------------------------------------------
# 3. pLDDT extraction
# ---------------------------------------------------------------------------


def extract_plddt_scores(pdb_content: str) -> list[tuple[int, float]]:
    """Extract per-residue pLDDT scores from an ESMFold PDB file.

    ESMFold stores the predicted Local Distance Difference Test (pLDDT)
    confidence score in the B-factor column of each ATOM record.  This
    function reads only CA (alpha-carbon) atoms to produce one score per
    residue.

    Args:
        pdb_content: Raw PDB file content as a string.

    Returns:
        List of ``(residue_number, plddt_score)`` tuples sorted by
        residue number.  Returns an empty list if no CA atoms are found.
    """
    scores: list[tuple[int, float]] = []

    for line in pdb_content.splitlines():
        if not line.startswith("ATOM"):
            continue

        # PDB fixed-width format:
        #   columns 13-16: atom name
        #   columns 23-26: residue sequence number
        #   columns 61-66: B-factor (temperature factor / pLDDT)
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        try:
            residue_num = int(line[22:26].strip())
            b_factor = float(line[60:66].strip())
            scores.append((residue_num, b_factor))
        except (ValueError, IndexError) as exc:
            logger.warning("Failed to parse ATOM line: %s (%s)", line.rstrip(), exc)
            continue

    if not scores:
        logger.warning("No CA atoms found in PDB content")

    return scores


# ---------------------------------------------------------------------------
# 4. PyMOL script generator
# ---------------------------------------------------------------------------


def generate_pymol_script(
    pdb_filename: str,
    regions: list[ConstructRegion],
    output_path: str | Path,
) -> Path:
    """Generate a PyMOL ``.pml`` script that colors the structure by region.

    The script loads the PDB file, sets a white background, applies a
    default gray color, then selects and colors each construct region.

    Args:
        pdb_filename: Name of the PDB file to load (relative path).
        regions: Ordered list of construct regions from
            :func:`map_construct_regions`.
        output_path: Destination file path for the ``.pml`` script.

    Returns:
        :class:`~pathlib.Path` to the saved script file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = [
        f"load {pdb_filename}",
        "bg_color white",
        "color gray, all",
        "cartoon loop",
    ]

    for region in regions:
        selection_name = region.name
        resi_range = f"{region.start}-{region.end}"
        hex_color = region.color.lstrip("#")
        lines.append(f"select {selection_name}, resi {resi_range}")
        lines.append(f"color 0x{hex_color}, {selection_name}")

    lines.append("deselect")
    lines.append("zoom")
    lines.append("")  # trailing newline

    script_content = "\n".join(lines)
    output_path.write_text(script_content, encoding="utf-8")
    logger.info("PyMOL script saved to %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# 5. ChimeraX script generator
# ---------------------------------------------------------------------------


def generate_chimerax_script(
    pdb_filename: str,
    regions: list[ConstructRegion],
    output_path: str | Path,
) -> Path:
    """Generate a UCSF ChimeraX ``.cxc`` command script.

    The script opens the PDB file, applies a default gray color, then
    colors each construct region by its residue range.

    Args:
        pdb_filename: Name of the PDB file to open (relative path).
        regions: Ordered list of construct regions from
            :func:`map_construct_regions`.
        output_path: Destination file path for the ``.cxc`` script.

    Returns:
        :class:`~pathlib.Path` to the saved script file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = [
        f"open {pdb_filename}",
        "color gray",
    ]

    for region in regions:
        resi_range = f"{region.start}-{region.end}"
        lines.append(f"color {region.color} :{resi_range}")

    lines.append("")  # trailing newline

    script_content = "\n".join(lines)
    output_path.write_text(script_content, encoding="utf-8")
    logger.info("ChimeraX script saved to %s", output_path)
    return output_path
