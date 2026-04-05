"""Build a compound library for molecular docking against top drug targets.

Queries ChEMBL for known inhibitors, merges a curated repurposing library
of approved antileishmanial drugs, deduplicates by SMILES, and prepares
ligand PDBQT files for downstream virtual screening.

Usage:
    python -m drug_targets.07_compound_library
    python -m drug_targets.07_compound_library --top-n 3 --max-compounds 50
    python -m drug_targets.07_compound_library --dry-run
    python -m drug_targets.07_compound_library --force
"""

from __future__ import annotations

import csv
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger
from core.models import DockingCompound

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RESULTS_DIR: Final[Path] = Path("results")
OUTPUT_FILE: Final[Path] = RESULTS_DIR / "compound_library.csv"
COMPOUNDS_DIR: Final[Path] = Path("data/compounds")

CHEMBL_ACTIVITY_URL: Final[str] = (
    "https://www.ebi.ac.uk/chembl/api/data/activity.json"
)
MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 60
USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (compound-library)"

CSV_COLUMNS: Final[list[str]] = [
    "compound_id",
    "name",
    "smiles",
    "source",
    "is_approved_drug",
    "target_gene_name",
]

logger = get_logger("compound_library")

# ---------------------------------------------------------------------------
# Curated repurposing library
# ---------------------------------------------------------------------------

REPURPOSING_LIBRARY: Final[list[dict[str, str]]] = [
    {
        "compound_id": "CHEMBL95",
        "name": "Miltefosine",
        "smiles": "CCCCCCCCCCCCCCCCOP(=O)([O-])OCC[N+](C)(C)C",
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL75",
        "name": "Amphotericin B",
        "smiles": (
            "O[C@H]1C=CC=CC=CC=CC=CC=CC=C[C@@H](C[C@@H]2O[C@@H]"
            "(O[C@@H]3[C@@H](O)[C@@H](N)[C@@H](O)[C@@H](C)O3)"
            "[C@H](C)[C@@H](O)[C@@H]2O)OC(=O)C=CC=C[C@@H](O)"
            "C[C@@H](O)CC(=O)O[C@@H]1C"
        ),
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL36",
        "name": "Ketoconazole",
        "smiles": (
            "ClC1=CC=C(C(OC[C@@H]2CO[C@@](CN3C=CN=C3)(O2)"
            "C2=CC=C(Cl)C=C2)C2=CC=C(C=C2)OC)C=C1"
        ),
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL607",
        "name": "Itraconazole",
        "smiles": (
            "CCC(C)N1N=CN(C1=O)C1=CC=C(OC[C@@H]2CO[C@@]"
            "(CN3C=CN=C3)(O2)C2=CC=C(Cl)C=C2Cl)C=C1"
        ),
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL185",
        "name": "Fluconazole",
        "smiles": "OC(CN1C=NC=N1)(CN1C=NC=N1)C1=CC=C(F)C=C1F",
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL194",
        "name": "Pentamidine",
        "smiles": "N=C(N)C1=CC=C(OCCCCCOC2=CC=C(C(=N)N)C=C2)C=C1",
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL118",
        "name": "Paromomycin",
        "smiles": (
            "NC[C@H]1OC(O[C@@H]2[C@@H](O)[C@H]"
            "(O[C@@H]3[C@@H](O)[C@H](N)[C@@H](O)[C@H](CO)O3)"
            "O[C@@H]2CO)[C@H](N)[C@@H](O)[C@@H]1O"
        ),
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL611",
        "name": "Sodium stibogluconate",
        "smiles": (
            "[Na+].[Na+].[Na+].OC(=O)[C@H](O)[C@@H](O)"
            "[C@H](O[Sb]([O-])[O-])[C@H](O)CO."
            "[O-][Sb]1O[C@@H]([C@H](O)CO)[C@H](O)[C@H](C(O)=O)O1"
        ),
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL1542",
        "name": "Allopurinol",
        "smiles": "O=C1N=CN=C2NNC=C12",
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL92",
        "name": "Methotrexate",
        "smiles": (
            "CN(CC1=NC2=C(N=C1)N=C(N)N=C2N)C1=CC=C"
            "(C(=O)N[C@@H](CCC(O)=O)C(O)=O)C=C1"
        ),
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL428647",
        "name": "Sitamaquine",
        "smiles": "COC1=CC2=C(NC(CCCCCN(CC)CC)C)C=CC2=CC1=C",
        "source": "repurposing_hub",
    },
    {
        "compound_id": "CHEMBL1200986",
        "name": "Sinefungin",
        "smiles": (
            "N[C@@H](CC[C@H](N)C(O)=O)C[C@@H]1OC"
            "([C@@H](O)[C@H]1O)N1C=NC2=C(N)N=CN=C12"
        ),
        "source": "repurposing_hub",
    },
]


# ---------------------------------------------------------------------------
# Data container for compound entries before persistence
# ---------------------------------------------------------------------------


@dataclass
class _CompoundEntry:
    """Internal container for a compound associated with a specific target."""

    compound_id: str
    name: str
    smiles: str
    source: str
    is_approved_drug: bool
    target_gene_name: str


# ---------------------------------------------------------------------------
# ChEMBL API helpers
# ---------------------------------------------------------------------------


def _chembl_request(params: dict[str, str]) -> dict | None:
    """Send a GET request to the ChEMBL activity endpoint with retry/backoff.

    Args:
        params: Query parameters for the ChEMBL activity endpoint.

    Returns:
        Parsed JSON response dict on success, or ``None`` after all retries
        are exhausted.
    """
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}

    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(
                CHEMBL_ACTIVITY_URL,
                params=params,
                headers=headers,
                timeout=REQUEST_TIMEOUT,
            )
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            delay = RETRY_DELAYS[attempt] if attempt < len(RETRY_DELAYS) else 8
            logger.warning(
                "ChEMBL request failed (attempt %d/%d): %s. "
                "Retrying in %ds.",
                attempt + 1,
                MAX_RETRIES,
                exc,
                delay,
            )
            time.sleep(delay)

    logger.error("ChEMBL request failed after %d retries.", MAX_RETRIES)
    return None


def _search_chembl_for_target(
    gene_name: str,
    enzyme_class: str,
    max_compounds: int,
) -> list[_CompoundEntry]:
    """Search ChEMBL for known inhibitors of a given target.

    Queries against the target name and enzyme class, filtering for
    activities with IC50 or Ki measurements and reported SMILES.

    Args:
        gene_name: Gene/enzyme name (e.g. "TryR").
        enzyme_class: Enzyme class descriptor (e.g. "trypanothione_reductase").
        max_compounds: Maximum number of compounds to return per target.

    Returns:
        List of ``_CompoundEntry`` objects found for this target.
    """
    entries: list[_CompoundEntry] = []
    seen_smiles: set[str] = set()

    # Try both gene name and enzyme class as search terms.
    search_terms = [gene_name]
    if enzyme_class:
        # Convert underscores to spaces for the API query.
        search_terms.append(enzyme_class.replace("_", " "))

    for term in search_terms:
        if len(entries) >= max_compounds:
            break

        params = {
            "target_organism": "Leishmania",
            "target_pref_name__icontains": term,
            "standard_type__in": "IC50,Ki",
            "limit": str(min(max_compounds - len(entries), 50)),
            "format": "json",
        }

        data = _chembl_request(params)
        if data is None:
            continue

        activities = data.get("activities", [])
        logger.info(
            "ChEMBL returned %d activities for term '%s'.",
            len(activities),
            term,
        )

        for act in activities:
            smiles = act.get("canonical_smiles", "")
            chembl_id = act.get("molecule_chembl_id", "")
            mol_name = act.get("molecule_pref_name") or chembl_id

            if not smiles or not chembl_id:
                continue
            if smiles in seen_smiles:
                continue

            seen_smiles.add(smiles)

            # Compounds at max_phase >= 3 are considered approved/late-stage.
            max_phase = 0
            try:
                max_phase = int(act.get("molecule_max_phase") or 0)
            except (ValueError, TypeError):
                pass

            entries.append(
                _CompoundEntry(
                    compound_id=chembl_id,
                    name=mol_name,
                    smiles=smiles,
                    source="chembl",
                    is_approved_drug=max_phase >= 3,
                    target_gene_name=gene_name,
                )
            )

            if len(entries) >= max_compounds:
                break

    logger.info(
        "Found %d unique compounds from ChEMBL for target %s.",
        len(entries),
        gene_name,
    )
    return entries


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------


def _deduplicate_compounds(
    compounds: list[_CompoundEntry],
) -> list[_CompoundEntry]:
    """Remove duplicate compounds by canonical SMILES.

    When duplicates are found, the first occurrence is kept.  Approved
    drugs are preferred over non-approved duplicates.

    Args:
        compounds: List of compound entries, possibly with duplicates.

    Returns:
        Deduplicated list preserving insertion order.
    """
    seen: dict[str, int] = {}
    result: list[_CompoundEntry] = []

    for entry in compounds:
        canonical = entry.smiles.strip()
        if canonical in seen:
            idx = seen[canonical]
            # Prefer approved drugs over non-approved.
            if entry.is_approved_drug and not result[idx].is_approved_drug:
                result[idx] = entry
            continue

        seen[canonical] = len(result)
        result.append(entry)

    removed = len(compounds) - len(result)
    if removed > 0:
        logger.info("Deduplicated %d compounds (removed %d duplicates).", len(result), removed)

    return result


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def _write_csv(compounds: list[_CompoundEntry], output_path: Path) -> None:
    """Write the compound library to a CSV file.

    Args:
        compounds: List of compound entries to write.
        output_path: Destination file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        for entry in compounds:
            writer.writerow({
                "compound_id": entry.compound_id,
                "name": entry.name,
                "smiles": entry.smiles,
                "source": entry.source,
                "is_approved_drug": entry.is_approved_drug,
                "target_gene_name": entry.target_gene_name,
            })

    logger.info("Wrote %d compounds to %s.", len(compounds), output_path)


# ---------------------------------------------------------------------------
# Supabase persistence
# ---------------------------------------------------------------------------


def _persist_compounds(compounds: list[_CompoundEntry]) -> None:
    """Persist compounds to Supabase via upsert_docking_compound.

    Imports the database function inside this function to avoid circular
    imports.  Handles missing table errors gracefully.

    Args:
        compounds: List of compound entries to persist.
    """
    try:
        from core.db import upsert_docking_compound
    except ImportError:
        logger.warning(
            "Could not import upsert_docking_compound. "
            "Skipping Supabase persistence."
        )
        return

    persisted = 0
    for entry in compounds:
        dc = DockingCompound(
            compound_id=entry.compound_id,
            name=entry.name,
            smiles=entry.smiles,
            source=entry.source,
            is_approved_drug=entry.is_approved_drug,
        )
        try:
            upsert_docking_compound(dc)
            persisted += 1
        except Exception:
            logger.warning(
                "Failed to persist compound %s (%s). "
                "Table may not exist yet.",
                entry.compound_id,
                entry.name,
            )

    logger.info("Persisted %d / %d compounds to Supabase.", persisted, len(compounds))


# ---------------------------------------------------------------------------
# PDBQT generation
# ---------------------------------------------------------------------------


def _smiles_to_pdbqt_rdkit(smiles: str, output_path: Path) -> bool:
    """Convert SMILES to PDBQT using RDKit and Meeko.

    Generates 3D coordinates via RDKit's ETKDG, optimises geometry with
    MMFF94, then writes a PDBQT file via Meeko.

    Args:
        smiles: Canonical SMILES string.
        output_path: Destination PDBQT file path.

    Returns:
        True if conversion succeeded, False otherwise.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import meeko

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            return False

        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

        preparator = meeko.MoleculePreparation()
        mol_setup = preparator.prepare(mol)

        # Meeko >= 0.5 returns a list of setups.
        if isinstance(mol_setup, list):
            mol_setup = mol_setup[0]

        pdbqt_string = meeko.PDBQTWriterLegacy.write_string(mol_setup)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(pdbqt_string, encoding="utf-8")
        return True

    except Exception as exc:
        logger.debug("RDKit/Meeko conversion failed: %s", exc)
        return False


def _smiles_to_pdbqt_obabel(smiles: str, output_path: Path) -> bool:
    """Convert SMILES to PDBQT using Open Babel as a fallback.

    Args:
        smiles: Canonical SMILES string.
        output_path: Destination PDBQT file path.

    Returns:
        True if conversion succeeded, False otherwise.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        result = subprocess.run(
            [
                "obabel",
                f"-:{smiles}",
                "--gen3d",
                "-opdbqt",
                "-O",
                str(output_path),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode == 0 and output_path.exists():
            return True

        logger.debug("Open Babel stderr: %s", result.stderr.strip())
        return False

    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        logger.debug("Open Babel conversion failed: %s", exc)
        return False


def _save_smiles_fallback(smiles: str, compound_id: str, output_dir: Path) -> None:
    """Save SMILES string to a text file as a last resort.

    Args:
        smiles: Canonical SMILES string.
        compound_id: Compound identifier for the filename.
        output_dir: Directory in which to write the file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    smi_path = output_dir / f"{compound_id}.smi"
    smi_path.write_text(f"{smiles} {compound_id}\n", encoding="utf-8")
    logger.warning(
        "Saved raw SMILES for %s to %s (manual PDBQT conversion needed).",
        compound_id,
        smi_path,
    )


def _prepare_pdbqt_files(compounds: list[_CompoundEntry]) -> int:
    """Generate PDBQT ligand files for all compounds.

    Tries RDKit + Meeko first, then Open Babel, then saves raw SMILES
    as a last resort.

    Args:
        compounds: List of compound entries with SMILES.

    Returns:
        Number of compounds successfully converted to PDBQT.
    """
    converted = 0

    for entry in compounds:
        gene_dir = COMPOUNDS_DIR / entry.target_gene_name
        pdbqt_path = gene_dir / f"{entry.compound_id}.pdbqt"

        if pdbqt_path.exists():
            logger.debug(
                "PDBQT already exists for %s, skipping.", entry.compound_id
            )
            converted += 1
            continue

        # Primary: RDKit + Meeko.
        if _smiles_to_pdbqt_rdkit(entry.smiles, pdbqt_path):
            logger.debug("Converted %s via RDKit/Meeko.", entry.compound_id)
            converted += 1
            continue

        # Fallback: Open Babel.
        if _smiles_to_pdbqt_obabel(entry.smiles, pdbqt_path):
            logger.debug("Converted %s via Open Babel.", entry.compound_id)
            converted += 1
            continue

        # Last resort: save SMILES for manual processing.
        _save_smiles_fallback(entry.smiles, entry.compound_id, gene_dir)

    logger.info(
        "Prepared PDBQT files for %d / %d compounds.", converted, len(compounds)
    )
    return converted


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def build_compound_library(
    top_n: int = 5,
    max_compounds: int = 100,
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Build a compound library for molecular docking against top drug targets.

    Loads the highest-scoring drug targets from Supabase, searches ChEMBL
    for known inhibitors, merges a curated repurposing library of approved
    antileishmanial drugs, deduplicates by SMILES, writes a CSV, persists
    compounds to Supabase, and generates PDBQT ligand files.

    Args:
        top_n: Number of top drug targets to build libraries for.
        max_compounds: Maximum number of ChEMBL compounds per target.
        force: Overwrite the output CSV even if it already exists.
        dry_run: Use only the repurposing library, skip ChEMBL queries
                 and PDBQT generation.

    Returns:
        Path to the output compound library CSV.

    Raises:
        RuntimeError: If no drug targets are available.
    """
    output_path = OUTPUT_FILE

    if output_path.exists() and not force:
        logger.info(
            "Output file already exists at %s. Use --force to rebuild.",
            output_path,
        )
        return str(output_path)

    # ------------------------------------------------------------------
    # Step 1: Load top-N drug targets from Supabase
    # ------------------------------------------------------------------

    from core.db import get_all_drug_targets

    all_targets = get_all_drug_targets()
    if not all_targets:
        raise RuntimeError(
            "No drug targets found in Supabase. "
            "Run the druggability pipeline first."
        )

    # Sort by druggability score descending, take top N.
    all_targets.sort(key=lambda t: t.druggability_score, reverse=True)
    targets = all_targets[:top_n]

    logger.info(
        "Selected top %d targets: %s",
        len(targets),
        ", ".join(f"{t.gene_name} ({t.druggability_score:.2f})" for t in targets),
    )

    # ------------------------------------------------------------------
    # Step 2: Collect compounds
    # ------------------------------------------------------------------

    all_compounds: list[_CompoundEntry] = []

    if dry_run:
        logger.info("[DRY RUN] Skipping ChEMBL queries.")
    else:
        # 2a: Search ChEMBL for known inhibitors of each target.
        for target in targets:
            chembl_compounds = _search_chembl_for_target(
                gene_name=target.gene_name,
                enzyme_class=target.enzyme_class,
                max_compounds=max_compounds,
            )
            all_compounds.extend(chembl_compounds)

    # 2b: Add repurposing library for each target.
    for target in targets:
        for drug in REPURPOSING_LIBRARY:
            all_compounds.append(
                _CompoundEntry(
                    compound_id=drug["compound_id"],
                    name=drug["name"],
                    smiles=drug["smiles"],
                    source=drug["source"],
                    is_approved_drug=True,
                    target_gene_name=target.gene_name,
                )
            )

    # 2c: Deduplicate by SMILES.
    all_compounds = _deduplicate_compounds(all_compounds)

    if not all_compounds:
        logger.warning("No compounds found for any target.")
        return str(output_path)

    # ------------------------------------------------------------------
    # Step 3: Save compound library CSV and persist to Supabase
    # ------------------------------------------------------------------

    _write_csv(all_compounds, output_path)

    if not dry_run:
        _persist_compounds(all_compounds)

    # ------------------------------------------------------------------
    # Step 4: Prepare ligand PDBQT files
    # ------------------------------------------------------------------

    if dry_run:
        logger.info("[DRY RUN] Skipping PDBQT generation.")
    else:
        _prepare_pdbqt_files(all_compounds)

    # ------------------------------------------------------------------
    # Step 5: Summary
    # ------------------------------------------------------------------

    target_names = {t.gene_name for t in targets}
    logger.info(
        "Built library of %d compounds for %d targets.",
        len(all_compounds),
        len(target_names),
    )

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Build a compound library for molecular docking against "
            "top L. infantum drug targets."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top drug targets to use (default: 5).",
    )
    parser.add_argument(
        "--max-compounds",
        type=int,
        default=100,
        help="Maximum ChEMBL compounds per target (default: 100).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild even if compound_library.csv already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Use repurposing library only, skip ChEMBL queries "
            "and PDBQT generation."
        ),
    )
    args = parser.parse_args()

    result = build_compound_library(
        top_n=args.top_n,
        max_compounds=args.max_compounds,
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Complete: %s", result)
