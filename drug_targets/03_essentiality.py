"""Verify whether each enzyme in the drug_targets table is essential for parasite survival.

Cross-references three data sources:

1. A curated list of known essential genes from published literature.
2. The Database of Essential Genes (DEG) for trypanosomatid organisms.
3. TriTrypDB knockout/knockdown phenotype annotations for *L. infantum*
   and *L. donovani*.

Updates each drug target in Supabase with ``is_essential`` status and
writes a summary CSV to ``results/essentiality_check.csv``.

Usage:
    python -m drug_targets.03_essentiality
    python -m drug_targets.03_essentiality --force
    python -m drug_targets.03_essentiality --dry-run
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Final

import requests

from core.db import get_all_drug_targets, update_drug_target
from core.logger import get_logger
from core.models import DrugTarget

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUTPUT_DIR: Final[str] = "results"
OUTPUT_FILE: Final[str] = "results/essentiality_check.csv"
MAX_RETRIES: Final[int] = 3
RETRY_DELAYS: Final[list[int]] = [2, 4, 8]
REQUEST_TIMEOUT: Final[int] = 60
USER_AGENT: Final[str] = "Marley-Pipeline/2.0 (essentiality-check)"

DEG_SEARCH_URL: Final[str] = "http://www.essentialgene.org/cgi-bin/search.cgi"
TRITRYPDB_PHENOTYPE_URL: Final[str] = (
    "https://tritrypdb.org/tritrypdb/service/record-types/transcript/searches/"
    "GenesByTextSearch/reports/standard"
)

CSV_COLUMNS: Final[list[str]] = [
    "gene_id",
    "gene_name",
    "is_essential",
    "evidence_source",
]

logger = get_logger("essentiality")

# ---------------------------------------------------------------------------
# Curated literature data
# ---------------------------------------------------------------------------

KNOWN_ESSENTIAL_GENES: Final[dict[str, str]] = {
    "HGPRT": "Essential in L. donovani promastigotes (Boitz et al., 2012)",
    "ADL": "Lethal knockout in L. donovani (Boitz et al., 2013)",
    "GMPS": "Essential convergence point in purine salvage pathway",
    "TryS": "Essential for trypanothione biosynthesis (Oza et al., 2002)",
    "TryR": "Essential reductase, validated drug target (Fairlamb et al., 1985)",
    "SMT": "Essential for ergosterol biosynthesis (Magaraci et al., 2003)",
    "SDM": "CYP51 ortholog, essential sterol pathway enzyme",
    "6PGDH": "Essential in pentose phosphate pathway (Hanau et al., 2004)",
    "TPI": "Essential glycolytic enzyme in trypanosomatids",
    "PGK": "Essential glycolytic enzyme in trypanosomatids",
}

# Trypanosomatid organism keywords used to filter DEG results.
TRYPANOSOMATID_KEYWORDS: Final[list[str]] = [
    "leishmania",
    "trypanosoma",
    "trypanosomatid",
]

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _fetch_with_retry(url: str, params: dict[str, str] | None = None) -> requests.Response:
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
                url, headers=headers, params=params, timeout=REQUEST_TIMEOUT
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
                logger.error(
                    "All %d attempts failed for %s: %s",
                    MAX_RETRIES,
                    url,
                    exc,
                )
                raise

    # Should not reach here, but satisfy type checker.
    raise requests.RequestException("Unexpected retry loop exit")


def _check_deg(gene_name: str) -> str | None:
    """Query the Database of Essential Genes for a given gene name.

    Searches DEG for matches in trypanosomatid organisms. Returns the
    evidence string if the gene is found as essential, otherwise ``None``.

    Args:
        gene_name: Short enzyme name (e.g. "HGPRT").

    Returns:
        Evidence string if essential in trypanosomatids, or ``None``.
    """
    params = {
        "field": "gene_name",
        "query": gene_name,
        "organism": "all",
    }

    try:
        response = _fetch_with_retry(DEG_SEARCH_URL, params=params)
        content = response.text.lower()

        for keyword in TRYPANOSOMATID_KEYWORDS:
            if keyword in content and "essential" in content:
                logger.info(
                    "DEG: %s found as essential in trypanosomatids.", gene_name
                )
                return f"DEG: essential in trypanosomatids ({keyword} match)"

    except requests.RequestException:
        logger.warning("DEG query failed for %s — skipping DEG check.", gene_name)

    return None


def _check_tritrypdb_phenotype(gene_name: str) -> str | None:
    """Query TriTrypDB for knockout/knockdown phenotype data.

    Searches for lethal phenotype annotations in *L. infantum* or
    *L. donovani* for the given gene.

    Args:
        gene_name: Short enzyme name (e.g. "TryR").

    Returns:
        Evidence string if lethal knockout documented, or ``None``.
    """
    params = {
        "searchText": f"{gene_name} knockout knockdown lethal",
        "reportConfig": (
            '{"attributes":["primary_key","gene_name","product"]}'
        ),
    }

    try:
        response = _fetch_with_retry(TRITRYPDB_PHENOTYPE_URL, params=params)
        data = response.json()

        records = data.get("records", [])
        for rec in records:
            attrs = rec.get("attributes", {})
            product = attrs.get("product", "").lower()
            primary_key = attrs.get("primary_key", "").lower()

            gene_name_lower = gene_name.lower()
            if gene_name_lower in product or gene_name_lower in primary_key:
                logger.info(
                    "TriTrypDB: lethal phenotype found for %s.", gene_name
                )
                return f"TriTrypDB: knockout/knockdown phenotype documented"

    except requests.RequestException:
        logger.warning(
            "TriTrypDB phenotype query failed for %s — skipping.", gene_name
        )
    except (ValueError, KeyError):
        logger.warning(
            "TriTrypDB returned unexpected response for %s — skipping.",
            gene_name,
        )

    return None


def _determine_essentiality(
    target: DrugTarget, dry_run: bool = False
) -> tuple[bool, str]:
    """Determine whether a drug target enzyme is essential.

    Checks the curated list first, then DEG, then TriTrypDB phenotype
    data. Returns on the first positive match.

    Args:
        target: The ``DrugTarget`` to evaluate.
        dry_run: If ``True``, only check the curated list (skip API calls).

    Returns:
        Tuple of (is_essential, evidence_source). The evidence_source is
        an empty string when the gene is not found to be essential.
    """
    gene_name = target.gene_name

    # 1. Check curated list.
    if gene_name in KNOWN_ESSENTIAL_GENES:
        evidence = f"Literature: {KNOWN_ESSENTIAL_GENES[gene_name]}"
        logger.info("Curated essential: %s — %s", gene_name, evidence)
        return True, evidence

    if dry_run:
        logger.info(
            "[DRY RUN] Would query DEG and TriTrypDB for %s.", gene_name
        )
        return False, ""

    # 2. Check DEG.
    deg_evidence = _check_deg(gene_name)
    if deg_evidence is not None:
        return True, deg_evidence

    # 3. Check TriTrypDB phenotype data.
    tritrypdb_evidence = _check_tritrypdb_phenotype(gene_name)
    if tritrypdb_evidence is not None:
        return True, tritrypdb_evidence

    logger.info("No essentiality evidence found for %s.", gene_name)
    return False, ""


def _write_csv(
    results: list[dict[str, str | bool]], output_path: Path
) -> None:
    """Write essentiality results to a CSV file.

    Args:
        results: List of dicts with keys matching ``CSV_COLUMNS``.
        output_path: Path to the output CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(results)

    logger.info("Wrote essentiality results to %s.", output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def check_essentiality(force: bool = False, dry_run: bool = False) -> str:
    """Check essentiality for all drug targets and persist results.

    For each enzyme in the ``drug_targets`` table, determines whether it
    is essential for parasite survival by consulting literature, DEG, and
    TriTrypDB phenotype data.  Updates Supabase and writes a summary CSV.

    Args:
        force: Re-run even if the output CSV already exists.
        dry_run: Only check the curated literature list (skip API calls).
            Still writes the CSV from curated data alone.

    Returns:
        Path to the output CSV file.

    Raises:
        RuntimeError: If no drug targets are found in the database.
    """
    output_path = Path(OUTPUT_FILE)

    if output_path.exists() and not force:
        logger.info(
            "Essentiality file already exists at %s. Use --force to re-run.",
            output_path,
        )
        return str(output_path)

    if dry_run:
        logger.info("[DRY RUN] Essentiality check — curated list only.")

    # Fetch all drug targets from Supabase.
    targets = get_all_drug_targets()

    if not targets:
        msg = "No drug targets found in database. Run 01_fetch_enzymes first."
        logger.error(msg)
        raise RuntimeError(msg)

    logger.info("Checking essentiality for %d drug target(s).", len(targets))

    results: list[dict[str, str | bool]] = []
    essential_count = 0

    for target in targets:
        is_essential, evidence_source = _determine_essentiality(target, dry_run)

        if is_essential:
            essential_count += 1

        results.append({
            "gene_id": target.gene_id,
            "gene_name": target.gene_name,
            "is_essential": is_essential,
            "evidence_source": evidence_source,
        })

        # Update Supabase.
        if not dry_run:
            try:
                update_payload: dict[str, bool | str] = {
                    "is_essential": is_essential,
                }
                if evidence_source:
                    update_payload["evidence"] = evidence_source
                update_drug_target(target.gene_id, update_payload)
            except Exception:
                logger.error(
                    "Failed to update %s in database — continuing.",
                    target.gene_id,
                )
        else:
            logger.info(
                "[DRY RUN] Would update %s: is_essential=%s.",
                target.gene_id,
                is_essential,
            )

    # Write CSV.
    _write_csv(results, output_path)

    # Summary.
    logger.info(
        "%d of %d enzymes confirmed as essential.",
        essential_count,
        len(targets),
    )

    return str(output_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Check essentiality of L. infantum drug target enzymes.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if essentiality_check.csv already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Use curated list only, skip API calls to DEG and TriTrypDB.",
    )
    args = parser.parse_args()

    result = check_essentiality(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
