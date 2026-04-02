"""Run the full Marley pipeline on all 8527 proteins.

This script runs each stage sequentially, saving progress between stages
so it can be resumed if interrupted.
"""

import importlib
import time
import csv
from pathlib import Path
from Bio import SeqIO

from core.logger import get_logger
from core.models import Candidate, STATUS_APPROVED, STATUS_REJECTED
from core.models import CONSERVATION_WEIGHT, IMMUNOGENICITY_WEIGHT
from core.db import upsert_candidate, update_candidate, get_all_candidates

logger = get_logger("full_pipeline")

# Import numbered modules
mod_filter = importlib.import_module("pipeline.02_filter_surface")
mod_conservation = importlib.import_module("pipeline.03_conservation")
mod_immuno = importlib.import_module("pipeline.04_immunogenicity")


def stage_signalp():
    """Stage 2: Run SignalP on all proteins."""
    input_file = "data/raw/proteins.fasta"
    output_file = "data/raw/surface_proteins.fasta"

    if Path(output_file).exists():
        count = sum(1 for _ in SeqIO.parse(output_file, "fasta"))
        if count > 0:
            logger.info("SignalP already done: %d surface proteins found. Skipping.", count)
            return count

    proteins = mod_filter.parse_fasta(input_file)
    logger.info("Starting SignalP on %d proteins ...", len(proteins))

    chunk_size = 500
    all_hits = set()
    all_results = {}

    for i in range(0, len(proteins), chunk_size):
        chunk = proteins[i:i + chunk_size]
        chunk_num = i // chunk_size + 1
        total_chunks = (len(proteins) + chunk_size - 1) // chunk_size
        logger.info("SignalP chunk %d/%d (%d proteins) ...", chunk_num, total_chunks, len(chunk))

        sequences = [(p[0], p[2]) for p in chunk]
        temp_fasta = mod_filter._write_temp_fasta(sequences)

        try:
            results = mod_filter.run_signalp(temp_fasta)
            all_results.update(results)
            hits = mod_filter._extract_signal_hits(results)
            all_hits.update(hits)
            logger.info("Chunk %d: %d hits (total so far: %d)", chunk_num, len(hits), len(all_hits))
        except Exception as exc:
            logger.error("Chunk %d failed: %s. Continuing...", chunk_num, exc)
        finally:
            Path(temp_fasta).unlink(missing_ok=True)

    # Save filtered proteins
    records = list(SeqIO.parse(input_file, "fasta"))
    filtered = [r for r in records if r.id in all_hits]
    SeqIO.write(filtered, output_file, "fasta")
    logger.info("SignalP complete: %d of %d proteins have signal peptide.", len(filtered), len(records))

    # Upsert to Supabase
    for r in filtered:
        c = Candidate(
            gene_id=r.id,
            gene_name=r.description.replace(r.id, "").strip(),
            sequence=str(r.seq),
            has_signal_peptide=True,
            filters_passed=["surface_filter"],
            status=STATUS_APPROVED,
        )
        try:
            upsert_candidate(c)
        except Exception:
            logger.warning("Failed to upsert %s", r.id)

    return len(filtered)


def stage_blast():
    """Stage 3: Run BLAST conservation analysis on surface proteins."""
    input_file = "data/raw/surface_proteins.fasta"
    output_file = "results/conserved_candidates.csv"
    threshold = mod_conservation.CONSERVATION_THRESHOLD

    # Check for existing progress
    done_ids = set()
    if Path(output_file).exists():
        with open(output_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                done_ids.add(row["gene_id"])
        logger.info("Resuming BLAST: %d already done.", len(done_ids))

    records = list(SeqIO.parse(input_file, "fasta"))
    remaining = [r for r in records if r.id not in done_ids]
    logger.info("Starting BLAST on %d proteins (%d remaining) ...", len(records), len(remaining))

    Path("results").mkdir(exist_ok=True)

    # Open CSV in append mode
    write_header = not Path(output_file).exists() or len(done_ids) == 0
    csvfile = open(output_file, "a", newline="")
    writer = csv.DictWriter(csvfile, fieldnames=["gene_id", "gene_name", "conservation_score", "status"])
    if write_header:
        writer.writeheader()

    approved = len([d for d in done_ids])  # rough count
    for i, record in enumerate(remaining, 1):
        gene_id = record.id
        gene_name = record.description.replace(record.id, "").strip()
        sequence = str(record.seq)

        logger.info("[%d/%d] BLAST: %s (%d aa) ...", i, len(remaining), gene_id, len(sequence))
        start = time.time()

        try:
            hits = mod_conservation.run_blast(sequence, gene_id)
            score = mod_conservation.calculate_conservation(hits)
        except Exception as exc:
            logger.error("BLAST failed for %s: %s. Setting score=0.", gene_id, exc)
            score = 0.0

        elapsed = time.time() - start
        status = STATUS_APPROVED if score >= threshold else STATUS_REJECTED

        writer.writerow({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "conservation_score": round(score, 4),
            "status": status,
        })
        csvfile.flush()

        try:
            update_candidate(gene_id, {
                "conservation_score": round(score, 4),
                "filters_passed": ["surface_filter", "conservation"] if status == STATUS_APPROVED else ["surface_filter"],
                "status": status,
            })
        except Exception:
            logger.warning("Failed to update Supabase for %s", gene_id)

        logger.info("  %s: score=%.4f (%s) [%.1fs]", gene_id, score, status, elapsed)
        if status == STATUS_APPROVED:
            approved += 1

    csvfile.close()
    logger.info("BLAST complete: %d approved candidates.", approved)
    return approved


def stage_iedb():
    """Stage 4: Run IEDB immunogenicity scoring on conserved candidates."""
    input_file = "results/conserved_candidates.csv"
    output_file = "results/scored_candidates.csv"
    fasta_file = "data/raw/surface_proteins.fasta"

    # Load sequences
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    # Load conserved candidates
    candidates = []
    with open(input_file) as f:
        for row in csv.DictReader(f):
            if row["status"] == STATUS_APPROVED:
                candidates.append(row)

    logger.info("Starting IEDB on %d conserved candidates ...", len(candidates))

    Path("results").mkdir(exist_ok=True)
    results = []

    for i, cand in enumerate(candidates, 1):
        gene_id = cand["gene_id"]
        gene_name = cand["gene_name"]
        conservation = float(cand["conservation_score"])
        sequence = sequences.get(gene_id, "")

        if not sequence:
            logger.warning("No sequence for %s, skipping.", gene_id)
            continue

        logger.info("[%d/%d] IEDB: %s ...", i, len(candidates), gene_id)
        start = time.time()

        all_preds = []
        for allele in mod_immuno.MHC_ALLELES:
            preds = mod_immuno.predict_binding(sequence, allele)
            all_preds.extend(preds)
            time.sleep(1)  # Rate limiting

        immunogenicity = mod_immuno.calculate_immunogenicity(all_preds)
        final = CONSERVATION_WEIGHT * conservation + IMMUNOGENICITY_WEIGHT * immunogenicity
        elapsed = time.time() - start

        logger.info("  %s: immuno=%.4f final=%.4f [%.1fs]", gene_id, immunogenicity, final, elapsed)

        try:
            update_candidate(gene_id, {
                "immunogenicity_score": round(immunogenicity, 4),
                "final_score": round(final, 4),
            })
        except Exception:
            logger.warning("Failed to update Supabase for %s", gene_id)

        results.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "conservation_score": round(conservation, 4),
            "immunogenicity_score": round(immunogenicity, 4),
            "final_score": round(final, 4),
        })

    # Write output CSV
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene_id", "gene_name", "conservation_score", "immunogenicity_score", "final_score"])
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda r: r["final_score"], reverse=True))

    logger.info("IEDB complete: %d candidates scored.", len(results))
    return len(results)


def stage_report():
    """Stage 5: Generate final report."""
    mod_report = importlib.import_module("pipeline.05_report")
    # Load validated antigens first
    mod_immuno.load_validated_antigens()
    csv_path, report_path = mod_report.generate_report()
    logger.info("Report generated: %s, %s", csv_path, report_path)


if __name__ == "__main__":
    overall_start = time.time()

    print("=" * 60)
    print("  MARLEY — Full Pipeline Run")
    print("=" * 60)

    # Stage 2: SignalP
    print("\n[STAGE 2] SignalP — filtering surface proteins ...")
    t = time.time()
    n_surface = stage_signalp()
    print(f"  Result: {n_surface} surface proteins ({time.time()-t:.0f}s)")

    # Stage 3: BLAST
    print(f"\n[STAGE 3] BLAST — conservation analysis on {n_surface} proteins ...")
    t = time.time()
    n_conserved = stage_blast()
    print(f"  Result: {n_conserved} conserved candidates ({time.time()-t:.0f}s)")

    # Stage 4: IEDB
    print(f"\n[STAGE 4] IEDB — immunogenicity scoring on {n_conserved} candidates ...")
    t = time.time()
    n_scored = stage_iedb()
    print(f"  Result: {n_scored} candidates scored ({time.time()-t:.0f}s)")

    # Stage 5: Report
    print("\n[STAGE 5] Generating report ...")
    t = time.time()
    stage_report()
    print(f"  Report generated ({time.time()-t:.0f}s)")

    total = time.time() - overall_start
    print("\n" + "=" * 60)
    print(f"  PIPELINE COMPLETE — Total time: {total/60:.1f} minutes")
    print("=" * 60)
