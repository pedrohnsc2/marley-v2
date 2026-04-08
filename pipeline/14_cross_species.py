"""Cross-species analysis of TryR drug target and vaccine epitope conservation.

Fetches trypanothione reductase (TryR) sequences from multiple Leishmania
and Trypanosoma species, computes pairwise identity against L. infantum,
checks active-site residue conservation, optionally docks MRL-003 into
predicted structures of high-identity orthologs, and evaluates vaccine
epitope conservation across species.

Usage:
    python -m pipeline.14_cross_species
    python -m pipeline.14_cross_species --force
    python -m pipeline.14_cross_species --dry-run
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("cross_species")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTPUT_DIR: Final[Path] = Path("results/cross_species")

CONSTRUCT_CARD_PATH: Final[str] = "results/construct/construct_card.json"

ESMFOLD_URL: Final[str] = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_TIMEOUT: Final[int] = 180
ESMFOLD_MAX_RETRIES: Final[int] = 2
ESMFOLD_RETRY_DELAY: Final[int] = 10

UNIPROT_SEARCH_URL: Final[str] = (
    "https://rest.uniprot.org/uniprotkb/search"
)
UNIPROT_TIMEOUT: Final[int] = 30
API_DELAY: Final[float] = 1.0

# Tool paths inside the project virtualenv.
OBABEL_PATH: Final[str] = ".venv/bin/obabel"
VINA_PATH: Final[str] = ".venv/bin/vina"

# MRL-003 ligand PDBQT (pre-prepared by the drug target pipeline).
MRL003_PDBQT: Final[Path] = Path(
    "data/compounds/TryR/MRL-003_amide_tail.pdbqt"
)

# L. infantum validated docking affinity from README / de-novo design.
LINFANTUM_MRL003_AFFINITY: Final[float] = -7.74

# Blind docking box size (Angstroms per axis).
BLIND_BOX_SIZE: Final[float] = 40.0
VINA_EXHAUSTIVENESS: Final[int] = 16

# Species to analyze with their NCBI/UniProt taxonomy IDs.
SPECIES: Final[dict[str, str]] = {
    "L. infantum": "5671",
    "L. donovani": "5661",
    "L. major": "5664",
    "L. braziliensis": "5660",
    "L. mexicana": "5665",
    "L. tropica": "5669",
    "T. cruzi": "5693",
    "T. brucei": "5691",
}

# Identity threshold (%) above which ESMFold + docking is attempted.
IDENTITY_THRESHOLD_FOR_DOCKING: Final[float] = 90.0

# Epitope conservation threshold (fraction) for "conserved" call.
EPITOPE_IDENTITY_THRESHOLD: Final[float] = 0.80

# Active site residues in L. infantum TryR (1-based positions).
# Catalytic residues + substrate-binding pocket from docking contact map.
ACTIVE_SITE_RESIDUES: Final[dict[str, int]] = {
    "Cys52": 52,
    "Cys57": 57,
    "His461": 461,
    "Glu466": 466,
    "Tyr198": 198,
    "Ile199": 199,
    "Phe396": 396,
    "Leu399": 399,
    "Pro398": 398,
    "Ser464": 464,
}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> dict | list:
    """Load and return a JSON file."""
    with open(path) as fh:
        return json.load(fh)


def _write_json(data: object, path: Path) -> None:
    """Write *data* as pretty-printed JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)
    logger.info("Wrote %s", path)


def _write_csv(
    rows: list[dict],
    columns: list[str],
    path: Path,
) -> None:
    """Write a list of dicts as a CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    logger.info("Wrote %d rows to %s", len(rows), path)


def _write_text(text: str, path: Path) -> None:
    """Write plain text to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    logger.info("Wrote %s", path)


def _read_fasta_sequence(fasta_path: Path) -> str:
    """Read a single-record FASTA file and return the plain sequence."""
    lines: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    return "".join(lines)


# ---------------------------------------------------------------------------
# 1. Fetch TryR sequences from UniProt
# ---------------------------------------------------------------------------


def _fetch_tryr_fasta(taxid: str) -> str | None:
    """Fetch the top trypanothione reductase FASTA for a given taxon ID.

    Queries UniProt's REST search endpoint with a protein name + organism
    filter and returns the amino-acid sequence of the first hit.

    Args:
        taxid: NCBI taxonomy ID as a string.

    Returns:
        Amino-acid sequence string, or ``None`` on failure.
    """
    url = (
        f"{UNIPROT_SEARCH_URL}?query=protein_name:trypanothione+reductase"
        f"+organism_id:{taxid}&format=fasta&size=1"
    )
    try:
        resp = requests.get(url, timeout=UNIPROT_TIMEOUT)
        resp.raise_for_status()
    except requests.RequestException as exc:
        logger.warning("UniProt fetch failed for taxid %s: %s", taxid, exc)
        return None

    fasta_text = resp.text.strip()
    if not fasta_text:
        logger.warning("Empty FASTA response for taxid %s.", taxid)
        return None

    seq_lines = [
        line.strip()
        for line in fasta_text.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    if not seq_lines:
        logger.warning("No sequence lines in FASTA for taxid %s.", taxid)
        return None

    sequence = "".join(seq_lines)
    logger.info("Got %d-aa TryR for taxid %s.", len(sequence), taxid)
    return sequence


def fetch_all_tryr_sequences() -> dict[str, str]:
    """Fetch TryR sequences for all species in :data:`SPECIES`.

    Returns:
        Mapping of species name to amino-acid sequence.
    """
    sequences: dict[str, str] = {}
    for species_name, taxid in SPECIES.items():
        if sequences:
            time.sleep(API_DELAY)
        seq = _fetch_tryr_fasta(taxid)
        if seq:
            sequences[species_name] = seq
        else:
            logger.warning("Could not retrieve TryR for %s.", species_name)
    logger.info(
        "Retrieved TryR sequences for %d/%d species.",
        len(sequences),
        len(SPECIES),
    )
    return sequences


# ---------------------------------------------------------------------------
# 2. Pairwise alignment and identity calculation
# ---------------------------------------------------------------------------


def _pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Compute percent identity between two sequences via global alignment.

    Tries BioPython pairwise2 first; falls back to a simple
    length-normalized matching-residue count when BioPython is unavailable.

    Args:
        seq_a: First amino-acid sequence.
        seq_b: Second amino-acid sequence.

    Returns:
        Percent identity as a float in [0, 100].
    """
    try:
        from Bio import pairwise2  # type: ignore[import-untyped]

        alignments = pairwise2.align.globalxx(seq_a, seq_b, one_alignment_only=True)
        if not alignments:
            return 0.0
        aligned_a, aligned_b, _score, _begin, _end = alignments[0]
        matches = sum(
            1 for a, b in zip(aligned_a, aligned_b) if a == b and a != "-"
        )
        aligned_length = max(len(aligned_a), 1)
        return (matches / aligned_length) * 100.0
    except ImportError:
        pass

    # Fallback: simple matching residues over the longer sequence length.
    min_len = min(len(seq_a), len(seq_b))
    max_len = max(len(seq_a), len(seq_b))
    if max_len == 0:
        return 0.0
    matches = sum(1 for i in range(min_len) if seq_a[i] == seq_b[i])
    return (matches / max_len) * 100.0


def compute_identities(
    reference: str,
    sequences: dict[str, str],
) -> dict[str, float]:
    """Compute pairwise identity of each species TryR against the reference.

    Args:
        reference: L. infantum TryR sequence.
        sequences: Mapping of species name to TryR sequence.

    Returns:
        Mapping of species name to percent identity.
    """
    identities: dict[str, float] = {}
    for species_name, seq in sequences.items():
        if species_name == "L. infantum":
            identities[species_name] = 100.0
        else:
            identity = _pairwise_identity(reference, seq)
            identities[species_name] = round(identity, 1)
            logger.info(
                "TryR identity %s vs L. infantum: %.1f%%",
                species_name,
                identity,
            )
    return identities


# ---------------------------------------------------------------------------
# 3. Active site conservation
# ---------------------------------------------------------------------------


def _aligned_residue_at(
    ref_seq: str,
    query_seq: str,
    ref_position: int,
) -> tuple[str | None, str | None]:
    """Return the reference and query residues at a 1-based reference position.

    Performs a global alignment (BioPython preferred) and maps the reference
    position to the aligned coordinates.

    Args:
        ref_seq: Reference (L. infantum) sequence.
        query_seq: Query (other species) sequence.
        ref_position: 1-based position in the reference.

    Returns:
        Tuple of (reference_residue, query_residue), or (None, None) on
        failure.
    """
    try:
        from Bio import pairwise2  # type: ignore[import-untyped]

        alignments = pairwise2.align.globalxx(
            ref_seq, query_seq, one_alignment_only=True,
        )
        if not alignments:
            return None, None
        aligned_ref, aligned_query, *_ = alignments[0]

        ref_idx = 0  # counts non-gap positions in the reference
        for i, (r, q) in enumerate(zip(aligned_ref, aligned_query)):
            if r != "-":
                ref_idx += 1
            if ref_idx == ref_position:
                ref_res = r if r != "-" else None
                query_res = q if q != "-" else None
                return ref_res, query_res
        return None, None
    except ImportError:
        pass

    # Fallback: direct positional comparison (no gaps assumed).
    idx = ref_position - 1
    ref_res = ref_seq[idx] if idx < len(ref_seq) else None
    query_res = query_seq[idx] if idx < len(query_seq) else None
    return ref_res, query_res


def check_active_site_conservation(
    ref_seq: str,
    query_seq: str,
) -> dict[str, bool]:
    """Check whether each active-site residue is conserved in the query.

    Args:
        ref_seq: L. infantum TryR sequence.
        query_seq: Other species TryR sequence.

    Returns:
        Mapping of residue label to conservation boolean.
    """
    result: dict[str, bool] = {}
    for label, pos in ACTIVE_SITE_RESIDUES.items():
        ref_res, query_res = _aligned_residue_at(ref_seq, query_seq, pos)
        conserved = (
            ref_res is not None
            and query_res is not None
            and ref_res == query_res
        )
        result[label] = conserved
    return result


# ---------------------------------------------------------------------------
# 4. Cross-species docking (ESMFold + Vina)
# ---------------------------------------------------------------------------


def _predict_structure_esmfold(sequence: str, output_pdb: Path) -> bool:
    """Predict 3D structure via ESMFold and save to a PDB file.

    Args:
        sequence: Amino-acid sequence (single-letter).
        output_pdb: Destination PDB path.

    Returns:
        True on success, False on failure.
    """
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    clean_seq = sequence.strip()

    for attempt in range(1, ESMFOLD_MAX_RETRIES + 2):
        try:
            logger.info(
                "ESMFold request attempt %d/%d (%d residues)",
                attempt,
                ESMFOLD_MAX_RETRIES + 1,
                len(clean_seq),
            )
            resp = requests.post(
                ESMFOLD_URL,
                data=clean_seq,
                headers={"Content-Type": "text/plain"},
                timeout=ESMFOLD_TIMEOUT,
            )
            resp.raise_for_status()

            pdb_content = resp.text
            if "ATOM" not in pdb_content:
                logger.error("ESMFold response lacks ATOM records.")
                return False

            output_pdb.write_text(pdb_content, encoding="utf-8")
            logger.info("PDB saved to %s", output_pdb)
            return True

        except requests.RequestException as exc:
            logger.warning(
                "ESMFold attempt %d failed: %s", attempt, exc,
            )
            if attempt <= ESMFOLD_MAX_RETRIES:
                time.sleep(ESMFOLD_RETRY_DELAY)

    logger.error(
        "ESMFold prediction failed after %d attempts.",
        ESMFOLD_MAX_RETRIES + 1,
    )
    return False


def _compute_centroid(pdbqt_path: Path) -> dict[str, float]:
    """Parse PDBQT ATOM records and compute the protein centroid.

    Returns a grid box dictionary centered on the centroid with a uniform
    size of :data:`BLIND_BOX_SIZE` Angstroms per axis.

    Args:
        pdbqt_path: Path to a receptor PDBQT file.

    Returns:
        Dictionary with center_x/y/z and size_x/y/z.

    Raises:
        ValueError: If no ATOM records are found.
    """
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []

    with open(pdbqt_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except (ValueError, IndexError):
                    continue

    if not xs:
        msg = f"No ATOM records found in {pdbqt_path}"
        raise ValueError(msg)

    return {
        "center_x": sum(xs) / len(xs),
        "center_y": sum(ys) / len(ys),
        "center_z": sum(zs) / len(zs),
        "size_x": BLIND_BOX_SIZE,
        "size_y": BLIND_BOX_SIZE,
        "size_z": BLIND_BOX_SIZE,
    }


def _pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path) -> bool:
    """Convert a PDB file to PDBQT using Open Babel.

    Args:
        pdb_path: Input PDB path.
        pdbqt_path: Output PDBQT path.

    Returns:
        True on success, False on failure.
    """
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        OBABEL_PATH,
        "-ipdb", str(pdb_path),
        "-opdbqt",
        "-O", str(pdbqt_path),
        "-xr",  # receptor mode: add partial charges
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=120,
        )
        logger.info("Converted %s -> %s", pdb_path.name, pdbqt_path.name)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        logger.warning("obabel conversion failed: %s", exc)
        return False
    except subprocess.TimeoutExpired:
        logger.warning("obabel timed out converting %s.", pdb_path.name)
        return False


def _parse_vina_output(stdout: str) -> float | None:
    """Extract best-pose binding affinity from Vina CLI output.

    Args:
        stdout: Standard output captured from the vina command.

    Returns:
        Binding affinity in kcal/mol, or None if parsing fails.
    """
    pattern = re.compile(
        r"^\s*1\s+([-+]?\d+\.?\d*)\s+",
    )
    for line in stdout.splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1))
    return None


def _dock_mrl003(receptor_pdbqt: Path, output_pdbqt: Path) -> float | None:
    """Dock MRL-003 into a receptor using blind docking with AutoDock Vina.

    Args:
        receptor_pdbqt: Path to the receptor PDBQT file.
        output_pdbqt: Path to write the docked pose.

    Returns:
        Best binding affinity in kcal/mol, or None on failure.
    """
    if not MRL003_PDBQT.exists():
        logger.warning("MRL-003 PDBQT not found at %s.", MRL003_PDBQT)
        return None

    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    try:
        grid_box = _compute_centroid(receptor_pdbqt)
    except ValueError as exc:
        logger.warning("Cannot compute centroid: %s", exc)
        return None

    cmd = [
        VINA_PATH,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(MRL003_PDBQT),
        "--center_x", str(round(grid_box["center_x"], 3)),
        "--center_y", str(round(grid_box["center_y"], 3)),
        "--center_z", str(round(grid_box["center_z"], 3)),
        "--size_x", str(grid_box["size_x"]),
        "--size_y", str(grid_box["size_y"]),
        "--size_z", str(grid_box["size_z"]),
        "--exhaustiveness", str(VINA_EXHAUSTIVENESS),
        "--num_modes", "5",
        "--out", str(output_pdbqt),
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=600,
        )
        affinity = _parse_vina_output(result.stdout)
        if affinity is not None:
            logger.info("MRL-003 affinity: %.2f kcal/mol", affinity)
        return affinity
    except FileNotFoundError:
        logger.warning("Vina not found at %s.", VINA_PATH)
        return None
    except subprocess.CalledProcessError as exc:
        logger.warning("Vina docking failed: %s", exc.stderr[:200])
        return None
    except subprocess.TimeoutExpired:
        logger.warning("Vina timed out (600s limit).")
        return None


def dock_cross_species(
    species_name: str,
    sequence: str,
    work_dir: Path,
) -> float | None:
    """Run ESMFold prediction and MRL-003 docking for one species.

    Args:
        species_name: Display name (e.g. "L. donovani").
        sequence: TryR amino-acid sequence.
        work_dir: Directory to store intermediate files.

    Returns:
        Best binding affinity, or None on failure.
    """
    safe_name = species_name.replace(" ", "_").replace(".", "")
    pdb_path = work_dir / f"{safe_name}_tryr.pdb"
    pdbqt_path = work_dir / f"{safe_name}_tryr_receptor.pdbqt"
    output_path = work_dir / f"{safe_name}_MRL003_out.pdbqt"

    logger.info("Predicting structure for %s TryR (%d aa)...", species_name, len(sequence))

    if not _predict_structure_esmfold(sequence, pdb_path):
        return None

    if not _pdb_to_pdbqt(pdb_path, pdbqt_path):
        return None

    return _dock_mrl003(pdbqt_path, output_path)


# ---------------------------------------------------------------------------
# 5. Vaccine epitope conservation
# ---------------------------------------------------------------------------


def _fetch_protein_by_gene_id(gene_id: str, taxid: str) -> str | None:
    """Fetch a protein sequence from UniProt by gene ID and taxon.

    Extracts the gene name from the full gene_id string (e.g.
    ``"LINF_240013900-T1-p1"`` -> ``"LINF_240013900"``) and searches
    UniProt for it in the given organism.

    Args:
        gene_id: Gene identifier from the construct card.
        taxid: NCBI taxonomy ID for the species.

    Returns:
        Amino-acid sequence, or None on failure.
    """
    # Extract base gene name: strip transcript/protein suffixes.
    base_name = gene_id.split("-")[0] if "-" in gene_id else gene_id

    url = UNIPROT_SEARCH_URL
    params = {
        "query": f"gene:{base_name} AND organism_id:{taxid}",
        "format": "json",
        "size": "1",
    }

    try:
        resp = requests.get(url, params=params, timeout=UNIPROT_TIMEOUT)
        resp.raise_for_status()
    except requests.RequestException as exc:
        logger.debug(
            "UniProt search failed for %s in taxid %s: %s",
            base_name, taxid, exc,
        )
        return None

    data = resp.json()
    results = data.get("results", [])
    if not results:
        return None

    sequence = results[0].get("sequence", {}).get("value")
    return sequence


def _local_epitope_identity(epitope: str, protein: str) -> float:
    """Compute the best local identity of an epitope within a protein.

    Slides the epitope across the protein and returns the best match
    fraction.  Uses BioPython local alignment when available.

    Args:
        epitope: Short peptide sequence.
        protein: Full-length protein sequence.

    Returns:
        Best identity fraction in [0.0, 1.0].
    """
    if not protein or not epitope:
        return 0.0

    try:
        from Bio import pairwise2  # type: ignore[import-untyped]

        alignments = pairwise2.align.localxx(
            protein, epitope, one_alignment_only=True,
        )
        if alignments:
            _aligned_prot, _aligned_ep, score, _begin, _end = alignments[0]
            return score / len(epitope)
        return 0.0
    except ImportError:
        pass

    # Fallback: sliding window.
    ep_len = len(epitope)
    if ep_len > len(protein):
        return 0.0

    best = 0
    for i in range(len(protein) - ep_len + 1):
        window = protein[i : i + ep_len]
        matches = sum(1 for a, b in zip(epitope, window) if a == b)
        best = max(best, matches)

    return best / ep_len


def check_epitope_conservation(
    epitopes: list[dict],
    species_name: str,
    taxid: str,
) -> list[dict]:
    """Check conservation of each epitope in a given species.

    Args:
        epitopes: List of epitope dicts from the construct card (must
            contain ``"peptide"`` and ``"gene_id"``).
        species_name: Display name (e.g. "L. donovani").
        taxid: NCBI taxonomy ID.

    Returns:
        List of dicts with keys: species, epitope, gene_id, conserved,
        identity_pct.
    """
    results: list[dict] = []

    # Cache fetched proteins by gene_id to avoid redundant requests.
    protein_cache: dict[str, str | None] = {}

    for ep in epitopes:
        peptide = ep.get("peptide", "")
        gene_id = ep.get("gene_id", "")

        if gene_id not in protein_cache:
            time.sleep(API_DELAY)
            protein_cache[gene_id] = _fetch_protein_by_gene_id(gene_id, taxid)

        protein = protein_cache[gene_id]

        if protein is None:
            identity = 0.0
        else:
            identity = _local_epitope_identity(peptide, protein)

        conserved = identity >= EPITOPE_IDENTITY_THRESHOLD

        results.append({
            "species": species_name,
            "epitope": peptide,
            "gene_id": gene_id,
            "conserved": "yes" if conserved else "no",
            "identity_pct": round(identity * 100.0, 1),
        })

    conserved_count = sum(1 for r in results if r["conserved"] == "yes")
    logger.info(
        "%s: %d/%d epitopes conserved.",
        species_name,
        conserved_count,
        len(results),
    )
    return results


# ---------------------------------------------------------------------------
# 6. Report generation
# ---------------------------------------------------------------------------


def _build_report(
    tryr_rows: list[dict],
    vaccine_rows: list[dict],
    epitope_count: int,
) -> str:
    """Generate the Markdown cross-species report.

    Args:
        tryr_rows: Per-species TryR analysis results.
        vaccine_rows: Per-species epitope conservation results.
        epitope_count: Total number of epitopes evaluated.

    Returns:
        Report as a Markdown string.
    """
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines: list[str] = []

    lines.append("# Cross-Species Analysis")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")

    # -- TryR Conservation & MRL-003 Spectrum --
    lines.append("## TryR Conservation & MRL-003 Spectrum")
    lines.append("")
    lines.append(
        "| Species | TryR Identity | Active Site Conserved "
        "| MRL-003 Affinity | Would work? |"
    )
    lines.append(
        "|---------|:------------:|:---------------------:"
        "|:----------------:|:-----------:|"
    )
    for row in tryr_rows:
        species = row["species"]
        identity = f"{row['identity_pct']:.1f}%"
        conserved = row["active_site_conserved"]
        affinity_val = row.get("mrl003_affinity")
        if affinity_val is not None:
            affinity = f"{affinity_val:.2f} kcal/mol"
        else:
            affinity = "N/A"

        # Determine if it would work: identity >= 90 and active site >= 80%.
        total_sites = len(ACTIVE_SITE_RESIDUES)
        conserved_parts = conserved.split("/")
        conserved_count = int(conserved_parts[0]) if conserved_parts[0].isdigit() else 0
        conserved_frac = conserved_count / total_sites if total_sites else 0

        if species == "L. infantum":
            would_work = "YES (validated)"
        elif row["identity_pct"] >= 90.0 and conserved_frac >= 0.8:
            if affinity_val is not None and affinity_val <= -6.0:
                would_work = "YES"
            elif affinity_val is not None:
                would_work = "MAYBE"
            else:
                would_work = "LIKELY"
        elif row["identity_pct"] >= 70.0:
            would_work = "MAYBE"
        else:
            would_work = "NO"

        lines.append(
            f"| {species} | {identity} | {conserved} "
            f"| {affinity} | {would_work} |"
        )

    lines.append("")

    # -- Vaccine Epitope Conservation --
    lines.append("## Vaccine Epitope Conservation")
    lines.append("")
    lines.append("| Species | Epitopes Conserved | Coverage |")
    lines.append("|---------|:-----------------:|:--------:|")

    # Aggregate by species.
    species_epitope_counts: dict[str, tuple[int, int]] = {}
    for row in vaccine_rows:
        sp = row["species"]
        if sp not in species_epitope_counts:
            species_epitope_counts[sp] = (0, 0)
        total = species_epitope_counts[sp][1] + 1
        conserved = species_epitope_counts[sp][0] + (
            1 if row["conserved"] == "yes" else 0
        )
        species_epitope_counts[sp] = (conserved, total)

    for sp_name in SPECIES:
        if sp_name in species_epitope_counts:
            c, t = species_epitope_counts[sp_name]
            coverage = f"{c / t * 100:.0f}%" if t > 0 else "N/A"
            lines.append(f"| {sp_name} | {c}/{t} | {coverage} |")
        else:
            lines.append(f"| {sp_name} | N/A | N/A |")

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def analyze_cross_species(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the full cross-species analysis pipeline.

    Steps:
        1. Fetch TryR sequences from UniProt for all target species.
        2. Compute pairwise identity against L. infantum TryR.
        3. Check active-site residue conservation.
        4. For high-identity species (>90%), predict structure and dock
           MRL-003 (skipped in dry-run mode).
        5. Load vaccine epitopes and check conservation across species.
        6. Generate CSV and Markdown reports.

    Args:
        force: If True, overwrite existing output even if present.
        dry_run: If True, skip ESMFold and docking; generate reports
            with alignment and conservation data only.

    Returns:
        Path to the output Markdown report.
    """
    report_path = OUTPUT_DIR / "cross_species_report.md"

    if report_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            report_path,
        )
        return str(report_path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    work_dir = OUTPUT_DIR / "structures"
    work_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Fetch TryR sequences
    # ------------------------------------------------------------------
    logger.info("Fetching TryR sequences from UniProt...")
    tryr_sequences = fetch_all_tryr_sequences()

    if "L. infantum" not in tryr_sequences:
        logger.error("Could not retrieve L. infantum TryR. Cannot proceed.")
        sys.exit(1)

    ref_seq = tryr_sequences["L. infantum"]

    # ------------------------------------------------------------------
    # 2. Compute pairwise identities
    # ------------------------------------------------------------------
    logger.info("Computing pairwise identities...")
    identities = compute_identities(ref_seq, tryr_sequences)

    # ------------------------------------------------------------------
    # 3. Active site conservation
    # ------------------------------------------------------------------
    logger.info("Checking active site conservation...")
    active_site_results: dict[str, dict[str, bool]] = {}
    for species_name, seq in tryr_sequences.items():
        if species_name == "L. infantum":
            # All sites trivially conserved.
            active_site_results[species_name] = {
                label: True for label in ACTIVE_SITE_RESIDUES
            }
        else:
            active_site_results[species_name] = check_active_site_conservation(
                ref_seq, seq,
            )

    # ------------------------------------------------------------------
    # 4. Cross-species docking (only for >90% identity, skip in dry run)
    # ------------------------------------------------------------------
    affinities: dict[str, float | None] = {}
    affinities["L. infantum"] = LINFANTUM_MRL003_AFFINITY

    for species_name in tryr_sequences:
        if species_name == "L. infantum":
            continue

        identity = identities.get(species_name, 0.0)
        if identity < IDENTITY_THRESHOLD_FOR_DOCKING:
            logger.info(
                "Skipping docking for %s (%.1f%% identity < %.0f%% threshold).",
                species_name,
                identity,
                IDENTITY_THRESHOLD_FOR_DOCKING,
            )
            affinities[species_name] = None
            continue

        if dry_run:
            logger.info(
                "[DRY RUN] Would dock MRL-003 into %s TryR (%.1f%% identity).",
                species_name,
                identity,
            )
            affinities[species_name] = None
            continue

        affinity = dock_cross_species(
            species_name,
            tryr_sequences[species_name],
            work_dir,
        )
        affinities[species_name] = affinity

    # ------------------------------------------------------------------
    # Build TryR rows
    # ------------------------------------------------------------------
    total_sites = len(ACTIVE_SITE_RESIDUES)
    tryr_rows: list[dict] = []

    for species_name in SPECIES:
        if species_name not in tryr_sequences:
            tryr_rows.append({
                "species": species_name,
                "taxid": SPECIES[species_name],
                "identity_pct": 0.0,
                "active_site_conserved": f"0/{total_sites}",
                "mrl003_affinity": None,
            })
            continue

        identity = identities.get(species_name, 0.0)
        as_result = active_site_results.get(species_name, {})
        conserved_count = sum(1 for v in as_result.values() if v)

        tryr_rows.append({
            "species": species_name,
            "taxid": SPECIES[species_name],
            "identity_pct": identity,
            "active_site_conserved": f"{conserved_count}/{total_sites}",
            "mrl003_affinity": affinities.get(species_name),
        })

    # ------------------------------------------------------------------
    # 5. Vaccine epitope conservation
    # ------------------------------------------------------------------
    logger.info("Checking vaccine epitope conservation...")
    card_path = Path(CONSTRUCT_CARD_PATH)
    vaccine_rows: list[dict] = []

    if not card_path.exists():
        logger.warning(
            "Construct card not found at %s. Skipping epitope conservation.",
            card_path,
        )
        epitope_count = 0
    else:
        card = _load_json(card_path)
        epitopes = card.get("epitopes", []) if isinstance(card, dict) else []
        epitope_count = len(epitopes)

        # Deduplicate epitopes by peptide sequence.
        seen_peptides: set[str] = set()
        unique_epitopes: list[dict] = []
        for ep in epitopes:
            pep = ep.get("peptide", "")
            if pep and pep not in seen_peptides:
                seen_peptides.add(pep)
                unique_epitopes.append(ep)

        epitope_count = len(unique_epitopes)
        logger.info("Evaluating %d unique epitopes across species.", epitope_count)

        for species_name, taxid in SPECIES.items():
            if species_name == "L. infantum":
                # All epitopes are trivially conserved in the source species.
                for ep in unique_epitopes:
                    vaccine_rows.append({
                        "species": species_name,
                        "epitope": ep.get("peptide", ""),
                        "gene_id": ep.get("gene_id", ""),
                        "conserved": "yes",
                        "identity_pct": 100.0,
                    })
            else:
                rows = check_epitope_conservation(
                    unique_epitopes, species_name, taxid,
                )
                vaccine_rows.extend(rows)

    # ------------------------------------------------------------------
    # 6. Generate outputs
    # ------------------------------------------------------------------
    logger.info("Writing output files...")

    # TryR CSV
    tryr_csv_columns = [
        "species", "taxid", "identity_pct",
        "active_site_conserved", "mrl003_affinity",
    ]
    _write_csv(tryr_rows, tryr_csv_columns, OUTPUT_DIR / "cross_species_tryr.csv")

    # Vaccine CSV
    vaccine_csv_columns = [
        "species", "epitope", "gene_id", "conserved", "identity_pct",
    ]
    if vaccine_rows:
        _write_csv(
            vaccine_rows,
            vaccine_csv_columns,
            OUTPUT_DIR / "cross_species_vaccine.csv",
        )

    # Markdown report
    report = _build_report(tryr_rows, vaccine_rows, epitope_count)
    _write_text(report, report_path)

    # Summary JSON
    summary = {
        "generated_at": datetime.now().isoformat(),
        "dry_run": dry_run,
        "species_analyzed": len(tryr_sequences),
        "tryr_identity_range": {
            "min": min(identities.values()) if identities else 0.0,
            "max": max(identities.values()) if identities else 0.0,
        },
        "species_docked": sum(
            1 for v in affinities.values() if v is not None
        ),
        "epitopes_evaluated": epitope_count,
        "report_path": str(report_path),
    }
    _write_json(summary, OUTPUT_DIR / "cross_species_summary.json")

    logger.info("Cross-species analysis complete. Report: %s", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Cross-species analysis of TryR conservation, MRL-003 docking "
            "spectrum, and vaccine epitope conservation across Leishmania "
            "and Trypanosoma species."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing outputs even if present.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Fetch sequences and compute alignments/conservation but skip "
            "ESMFold structure prediction and docking."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = analyze_cross_species(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
