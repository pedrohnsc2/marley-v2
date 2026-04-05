"""Fetch annotated enzymes from *Leishmania infantum* via UniProt.

Downloads all protein sequences with EC number annotations for priority
metabolic pathways.  Saves a FASTA file and persists each enzyme in the
``drug_targets`` Supabase table with ``status='pending'``.

Primary source: UniProt REST API (organism_id:5671 for *L. infantum*).

Usage:
    python -m drug_targets.01_fetch_enzymes
    python -m drug_targets.01_fetch_enzymes --force
    python -m drug_targets.01_fetch_enzymes --dry-run
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Final

import requests

from core.logger import get_logger
from core.models import DrugTarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

UNIPROT_SEARCH_URL: Final[str] = "https://rest.uniprot.org/uniprotkb/search"
ORGANISM_ID: Final[str] = "5671"  # Leishmania infantum (taxonomy)
OUTPUT_DIR: Final[str] = "data/raw"
OUTPUT_FILE: Final[str] = "data/raw/enzymes.fasta"
MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 120
USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (drug-target-discovery)"

logger = get_logger("fetch_enzymes")

# ---------------------------------------------------------------------------
# Priority pathway definitions
# ---------------------------------------------------------------------------

PRIORITY_PATHWAYS: Final[dict[str, list[str]]] = {
    "purine_salvage": ["HGPRT", "XPRT", "APRT", "ADL", "GMPS", "AK"],
    "sterol_biosynthesis": ["SMT", "SDM", "SQS"],
    "trypanothione": ["TryS", "TryR"],
    "pentose_phosphate": ["6PGDH"],
    "glycolysis": ["TPI", "PGK"],
}

# Map enzyme name -> pathway for reverse lookup.
ENZYME_TO_PATHWAY: Final[dict[str, str]] = {
    enzyme: pathway
    for pathway, enzymes in PRIORITY_PATHWAYS.items()
    for enzyme in enzymes
}

# EC numbers associated with priority enzymes.
PRIORITY_EC_NUMBERS: Final[dict[str, str]] = {
    "HGPRT": "2.4.2.8",
    "XPRT": "2.4.2.22",
    "APRT": "2.4.2.7",
    "ADL": "4.3.2.2",
    "GMPS": "6.3.5.2",
    "AK": "2.7.4.3",
    "SMT": "2.1.1.41",
    "SDM": "1.14.13.70",
    "SQS": "2.5.1.21",
    "TryS": "6.3.1.9",
    "TryR": "1.8.1.12",
    "6PGDH": "1.1.1.44",
    "TPI": "5.3.1.1",
    "PGK": "2.7.2.3",
}

# Enzyme class descriptions.
ENZYME_CLASS_MAP: Final[dict[str, str]] = {
    "HGPRT": "hypoxanthine_guanine_phosphoribosyltransferase",
    "XPRT": "xanthine_phosphoribosyltransferase",
    "APRT": "adenine_phosphoribosyltransferase",
    "ADL": "adenylosuccinate_lyase",
    "GMPS": "GMP_synthetase",
    "AK": "adenylate_kinase",
    "SMT": "sterol_24C_methyltransferase",
    "SDM": "sterol_14alpha_demethylase",
    "SQS": "squalene_synthase",
    "TryS": "trypanothione_synthetase",
    "TryR": "trypanothione_reductase",
    "6PGDH": "6_phosphogluconate_dehydrogenase",
    "TPI": "triosephosphate_isomerase",
    "PGK": "phosphoglycerate_kinase",
}


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
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}

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


def _search_uniprot(ec_number: str, enzyme_name: str) -> list[dict[str, str]]:
    """Search UniProt for *L. infantum* enzymes by EC number.

    Args:
        ec_number: Enzyme Commission number (e.g. "2.4.2.8").
        enzyme_name: Short enzyme name for annotation.

    Returns:
        List of dicts with keys: gene_id, gene_name, product, sequence.
    """
    query = f"organism_id:{ORGANISM_ID} AND ec:{ec_number}"
    params = {
        "query": query,
        "format": "json",
        "fields": "accession,gene_names,protein_name,ec,sequence",
        "size": "50",
    }

    response = _fetch_with_retry(UNIPROT_SEARCH_URL, params=params)
    data = response.json()

    records: list[dict[str, str]] = []
    results = data.get("results", [])

    for entry in results:
        accession = entry.get("primaryAccession", "")
        sequence_data = entry.get("sequence", {})
        sequence = sequence_data.get("value", "")

        if not accession or not sequence:
            continue

        # Extract protein name.
        protein_desc = entry.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {}).get("value", "")
        if not full_name:
            sub_names = protein_desc.get("submissionNames", [])
            if sub_names:
                full_name = sub_names[0].get("fullName", {}).get("value", "")

        records.append({
            "gene_id": accession,
            "gene_name": enzyme_name,
            "product": full_name,
            "sequence": sequence,
        })

    return records


def _write_fasta(records: list[dict[str, str]], output_path: Path) -> None:
    """Write enzyme records to a FASTA file.

    Args:
        records: List of dicts with gene_id, gene_name, product, sequence.
        output_path: Path to the output FASTA file.
    """
    with open(output_path, "w", encoding="utf-8") as fh:
        for rec in records:
            header = f">{rec['gene_id']} | {rec['gene_name']} | {rec['product']}"
            fh.write(f"{header}\n{rec['sequence']}\n")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def fetch_enzymes(force: bool = False, dry_run: bool = False) -> str:
    """Fetch all priority enzymes from UniProt and persist them.

    Queries UniProt for each priority enzyme EC number filtered by
    *L. infantum* taxonomy, saves results to ``data/raw/enzymes.fasta``,
    and upserts each enzyme into the ``drug_targets`` Supabase table.

    Args:
        force: Re-fetch even if the output file already exists.
        dry_run: Log actions without making API calls or writing files.

    Returns:
        Path to the output FASTA file.

    Raises:
        RuntimeError: If no enzymes are found after querying all pathways.
    """
    output_path = Path(OUTPUT_FILE)

    if output_path.exists() and not force:
        logger.info(
            "Enzyme file already exists at %s. Use --force to re-fetch.",
            output_path,
        )
        return str(output_path)

    if dry_run:
        logger.info("[DRY RUN] Would fetch enzymes for %d pathways.", len(PRIORITY_PATHWAYS))
        for pathway, enzymes in PRIORITY_PATHWAYS.items():
            logger.info("[DRY RUN]   %s: %s", pathway, ", ".join(enzymes))
        return str(output_path)

    # Ensure output directory exists.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_records: list[dict[str, str]] = []
    pathway_counts: dict[str, int] = {}

    for enzyme_name, ec_number in PRIORITY_EC_NUMBERS.items():
        pathway = ENZYME_TO_PATHWAY[enzyme_name]
        logger.info(
            "Fetching %s (EC %s) from pathway '%s' via UniProt...",
            enzyme_name,
            ec_number,
            pathway,
        )

        try:
            records = _search_uniprot(ec_number, enzyme_name)

            if records:
                all_records.extend(records)
                pathway_counts[pathway] = pathway_counts.get(pathway, 0) + len(records)
                logger.info(
                    "  Found %d record(s) for %s.", len(records), enzyme_name
                )
            else:
                logger.warning("  No records found for %s (EC %s).", enzyme_name, ec_number)

        except requests.RequestException:
            logger.error(
                "Failed to fetch %s — continuing with remaining enzymes.",
                enzyme_name,
            )

    if not all_records:
        msg = "No enzymes found across any pathway. Check UniProt availability."
        logger.error(msg)
        raise RuntimeError(msg)

    # Write FASTA.
    _write_fasta(all_records, output_path)
    logger.info("Wrote %d enzymes to %s.", len(all_records), output_path)

    # Persist to Supabase.
    _persist_to_database(all_records)

    # Summary.
    pathways_found = len(pathway_counts)
    logger.info(
        "%d enzymes found in %d metabolic pathways.",
        len(all_records),
        pathways_found,
    )
    for pathway, count in sorted(pathway_counts.items()):
        logger.info("  %s: %d enzyme(s)", pathway, count)

    return str(output_path)


def _persist_to_database(records: list[dict[str, str]]) -> None:
    """Upsert fetched enzyme records into the drug_targets table.

    Args:
        records: List of enzyme dicts from UniProt parsing.
    """
    from core.db import upsert_drug_target

    for rec in records:
        enzyme_name = rec["gene_name"]
        target = DrugTarget(
            gene_id=rec["gene_id"],
            gene_name=enzyme_name,
            sequence=rec["sequence"],
            enzyme_class=ENZYME_CLASS_MAP.get(enzyme_name, ""),
            pathway=ENZYME_TO_PATHWAY.get(enzyme_name, ""),
            status="pending",
        )

        try:
            upsert_drug_target(target)
        except Exception:
            logger.error(
                "Failed to persist %s to database — continuing.", rec["gene_id"]
            )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Fetch annotated enzymes from L. infantum (UniProt).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-fetch even if enzymes.fasta already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would be done without making API calls.",
    )
    args = parser.parse_args()

    result = fetch_enzymes(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
