"""Score candidate antigens for B-cell and T-cell epitope immunogenicity using prediction tools.

Predicts MHC-I binding affinity for canine (DLA) alleles via the IEDB API,
then combines immunogenicity with conservation into a weighted final score.
"""

from __future__ import annotations

import csv
import io
import time
from pathlib import Path

import requests
from tqdm import tqdm

from core.db import update_candidate
from core.logger import get_logger
from core.models import CONSERVATION_WEIGHT, IMMUNOGENICITY_WEIGHT

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

IEDB_MHC_URL: str = "https://tools.iedb.org/mhci/"
IEDB_API_URL: str = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

INPUT_FILE: str = "results/conserved_candidates.csv"
OUTPUT_FILE: str = "results/scored_candidates.csv"

MHC_ALLELES: list[str] = ["DLA-88*50101", "DLA-88*50801"]
PEPTIDE_LENGTH: int = 9  # 9-mer peptides
IC50_THRESHOLD: float = 500  # nM -- below this = good binder
PREDICTION_METHOD: str = "recommended"

# Polite delay between IEDB API calls to avoid hammering the service.
API_DELAY_SECONDS: float = 1.0

logger = get_logger("immunogenicity")

# ---------------------------------------------------------------------------
# Helper: FASTA parsing
# ---------------------------------------------------------------------------


def _load_sequences_from_fasta(fasta_path: str) -> dict[str, str]:
    """Parse a FASTA file into a mapping of gene_id -> amino-acid sequence.

    The gene_id is extracted as the first whitespace-delimited token on each
    header line (the part after '>').

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Dictionary mapping gene identifiers to their full sequences.

    Raises:
        FileNotFoundError: If *fasta_path* does not exist.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(parts)

    return sequences


# ---------------------------------------------------------------------------
# Peptide generation
# ---------------------------------------------------------------------------


def generate_peptides(sequence: str, length: int = PEPTIDE_LENGTH) -> list[str]:
    """Generate all overlapping peptides of *length* from *sequence*.

    Args:
        sequence: Full amino-acid sequence string.
        length: Desired peptide length (default ``PEPTIDE_LENGTH``).

    Returns:
        List of peptide substrings.  Empty if the sequence is shorter
        than *length*.
    """
    if len(sequence) < length:
        return []
    return [sequence[i : i + length] for i in range(len(sequence) - length + 1)]


# ---------------------------------------------------------------------------
# IEDB API interaction
# ---------------------------------------------------------------------------


def predict_binding(
    sequence: str,
    allele: str,
    length: int = PEPTIDE_LENGTH,
) -> list[dict]:
    """Submit a sequence to the IEDB MHC-I binding prediction API.

    The IEDB tools API accepts form-encoded POST data and returns
    tab-separated text.

    Args:
        sequence: Amino-acid sequence to analyse.
        allele: MHC allele name (e.g. ``"DLA-88*50101"``).
        length: Peptide length for the prediction.

    Returns:
        List of dicts, each containing:
        - ``peptide`` (str): The predicted peptide.
        - ``allele`` (str): The MHC allele used.
        - ``ic50`` (float): Predicted IC50 in nM.
        - ``rank`` (float): Percentile rank.

        Returns an empty list when the API call fails or returns
        unparsable data.
    """
    payload = {
        "method": PREDICTION_METHOD,
        "sequence_text": sequence,
        "allele": allele,
        "length": str(length),
    }

    try:
        response = requests.post(IEDB_API_URL, data=payload, timeout=120)
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.warning(
            "IEDB API request failed for allele %s: %s", allele, exc,
        )
        return []

    return _parse_iedb_response(response.text, allele)


def _parse_iedb_response(text: str, allele: str) -> list[dict]:
    """Parse tab-separated IEDB MHC-I prediction output.

    Args:
        text: Raw response body from the IEDB API.
        allele: Allele string to attach to each record (fall-back if the
                response does not include it).

    Returns:
        Parsed list of prediction dicts.  Records with missing or
        unparsable IC50 values are silently skipped.
    """
    results: list[dict] = []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")

    for row in reader:
        try:
            peptide = row.get("peptide", "").strip()
            ic50_raw = row.get("ic50", row.get("ic50_score", "")).strip()
            rank_raw = row.get("percentile_rank", row.get("rank", "0")).strip()

            if not peptide or not ic50_raw:
                continue

            results.append(
                {
                    "peptide": peptide,
                    "allele": row.get("allele", allele).strip(),
                    "ic50": float(ic50_raw),
                    "rank": float(rank_raw) if rank_raw else 0.0,
                }
            )
        except (ValueError, KeyError) as exc:
            logger.debug("Skipping unparsable IEDB row: %s (%s)", row, exc)
            continue

    return results


# ---------------------------------------------------------------------------
# Immunogenicity scoring
# ---------------------------------------------------------------------------


def calculate_immunogenicity(predictions: list[dict]) -> float:
    """Derive an immunogenicity score from MHC-I binding predictions.

    The score is the fraction of unique peptides whose best (lowest) IC50
    across all alleles falls below ``IC50_THRESHOLD``.

    Args:
        predictions: List of prediction dicts (as returned by
                     ``predict_binding``).

    Returns:
        Score in [0.0, 1.0].  Returns 0.0 when *predictions* is empty.
    """
    if not predictions:
        return 0.0

    # Collect the best (lowest) IC50 per peptide across alleles.
    best_ic50: dict[str, float] = {}
    for pred in predictions:
        peptide = pred["peptide"]
        ic50 = pred["ic50"]
        if peptide not in best_ic50 or ic50 < best_ic50[peptide]:
            best_ic50[peptide] = ic50

    total_peptides = len(best_ic50)
    if total_peptides == 0:
        return 0.0

    good_binders = sum(1 for ic50 in best_ic50.values() if ic50 < IC50_THRESHOLD)
    return good_binders / total_peptides


# ---------------------------------------------------------------------------
# Main pipeline step
# ---------------------------------------------------------------------------


def score_immunogenicity(force: bool = False) -> str:
    """Run the immunogenicity scoring pipeline step.

    Reads conserved candidates from ``INPUT_FILE``, queries the IEDB API
    for MHC-I binding predictions, calculates immunogenicity and final
    scores, persists results to Supabase, and writes ``OUTPUT_FILE``.

    Args:
        force: If ``True``, re-score even when ``OUTPUT_FILE`` already
               exists.

    Returns:
        Path to the output CSV (``OUTPUT_FILE``).

    Raises:
        FileNotFoundError: If ``INPUT_FILE`` or the FASTA sequence file
                           is missing.
    """
    output_path = Path(OUTPUT_FILE)

    if output_path.exists() and not force:
        logger.info("Output file %s already exists. Use force=True to re-run.", OUTPUT_FILE)
        return OUTPUT_FILE

    # Ensure output directory exists.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load input candidates
    # ------------------------------------------------------------------
    input_path = Path(INPUT_FILE)
    if not input_path.exists():
        msg = f"Input file not found: {INPUT_FILE}"
        logger.error(msg)
        raise FileNotFoundError(msg)

    candidates: list[dict] = []
    with open(input_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            candidates.append(row)

    logger.info("Loaded %d conserved candidate(s) from %s.", len(candidates), INPUT_FILE)

    if not candidates:
        logger.warning("No candidates to score. Writing empty output.")
        _write_output([], output_path)
        return OUTPUT_FILE

    # ------------------------------------------------------------------
    # 2. Load sequences from FASTA
    # ------------------------------------------------------------------
    fasta_path = "data/raw/surface_proteins.fasta"
    sequences: dict[str, str] = {}
    if Path(fasta_path).exists():
        sequences = _load_sequences_from_fasta(fasta_path)
        logger.info("Loaded %d sequence(s) from %s.", len(sequences), fasta_path)
    else:
        logger.warning(
            "FASTA file %s not found. Will attempt to use sequence from DB or skip.",
            fasta_path,
        )

    # ------------------------------------------------------------------
    # 3. Score each candidate
    # ------------------------------------------------------------------
    scored: list[dict] = []

    for candidate in tqdm(candidates, desc="Immunogenicity scoring", unit="gene"):
        gene_id: str = candidate.get("gene_id", "")
        gene_name: str = candidate.get("gene_name", "")
        conservation_score: float = float(candidate.get("conservation_score", 0.0))

        sequence = sequences.get(gene_id, candidate.get("sequence", ""))

        if not sequence:
            logger.warning("No sequence for %s (%s). Setting immunogenicity to 0.", gene_id, gene_name)
            immunogenicity_score = 0.0
        else:
            all_predictions: list[dict] = []
            for allele in MHC_ALLELES:
                preds = predict_binding(sequence, allele, PEPTIDE_LENGTH)
                all_predictions.extend(preds)
                # Be polite to the IEDB service.
                time.sleep(API_DELAY_SECONDS)

            immunogenicity_score = calculate_immunogenicity(all_predictions)

        final_score = (
            CONSERVATION_WEIGHT * conservation_score
            + IMMUNOGENICITY_WEIGHT * immunogenicity_score
        )

        logger.info(
            "%s (%s): conservation=%.3f  immunogenicity=%.3f  final=%.3f",
            gene_id,
            gene_name,
            conservation_score,
            immunogenicity_score,
            final_score,
        )

        # Persist to Supabase.
        try:
            update_candidate(
                gene_id,
                {
                    "immunogenicity_score": round(immunogenicity_score, 4),
                    "final_score": round(final_score, 4),
                },
            )
        except Exception:
            logger.warning("Could not update Supabase for %s. Continuing.", gene_id)

        scored.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "conservation_score": round(conservation_score, 4),
                "immunogenicity_score": round(immunogenicity_score, 4),
                "final_score": round(final_score, 4),
            }
        )

    # ------------------------------------------------------------------
    # 4. Write output CSV
    # ------------------------------------------------------------------
    _write_output(scored, output_path)

    # ------------------------------------------------------------------
    # 5. Summary
    # ------------------------------------------------------------------
    if scored:
        best = max(scored, key=lambda r: r["final_score"])
        avg_final = sum(r["final_score"] for r in scored) / len(scored)
        logger.info(
            "Scoring complete. %d candidate(s) processed. "
            "Best: %s (%s) with final_score=%.3f. Average final_score=%.3f.",
            len(scored),
            best["gene_id"],
            best["gene_name"],
            best["final_score"],
            avg_final,
        )
    else:
        logger.info("Scoring complete. No candidates were processed.")

    return OUTPUT_FILE


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _write_output(rows: list[dict], path: Path) -> None:
    """Write scored candidate rows to a CSV file.

    Args:
        rows: List of dicts to write.
        path: Destination file path.
    """
    fieldnames = [
        "gene_id",
        "gene_name",
        "conservation_score",
        "immunogenicity_score",
        "final_score",
    ]

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Results written to %s (%d row(s)).", path, len(rows))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Step 4: Score candidate immunogenicity via IEDB MHC-I predictions.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run scoring even if the output file already exists.",
    )
    args = parser.parse_args()

    result_path = score_immunogenicity(force=args.force)
    print(f"Output: {result_path}")
