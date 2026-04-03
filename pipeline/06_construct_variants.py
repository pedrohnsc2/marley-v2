"""Generate multiple mRNA vaccine construct variants for comparison.

Orchestrates the creation of three distinct construct designs (A, B, C)
varying in epitope sources, epitope types, and adjuvant choice.  Produces
per-variant output directories and a unified comparison report.

Usage:
    python -m pipeline.06_construct_variants
    python -m pipeline.06_construct_variants --force
"""

from __future__ import annotations

import csv
import importlib
import json
from datetime import datetime
from pathlib import Path
from typing import Any

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Import functions from the numbered module via importlib
# ---------------------------------------------------------------------------

_mod06 = importlib.import_module("pipeline.06_construct")

select_epitopes = _mod06.select_epitopes
assemble_construct = _mod06.assemble_construct
reverse_translate_optimized = _mod06.reverse_translate_optimized
assemble_mrna = _mod06.assemble_mrna
compute_physicochemical = _mod06.compute_physicochemical
predict_antigenicity = _mod06.predict_antigenicity
predict_allergenicity = _mod06.predict_allergenicity
calculate_gc_content = _mod06.calculate_gc_content
write_outputs = _mod06.write_outputs

SelectedEpitope = _mod06.SelectedEpitope
ConstructCard = _mod06.ConstructCard

SIGNAL_PEPTIDES = _mod06.SIGNAL_PEPTIDES
ADJUVANTS = _mod06.ADJUVANTS
LINKERS = _mod06.LINKERS
RESTRICTION_SITES = _mod06.RESTRICTION_SITES
get_optimal_codon = _mod06.get_optimal_codon

_load_sequences_from_fasta = _mod06._load_sequences_from_fasta

# HTL epitope selection -- the other agents are adding this to 06_construct.
# Import it if available; otherwise define a stub that raises at call time.
select_htl_epitopes = getattr(_mod06, "select_htl_epitopes", None)

logger = get_logger("construct_variants")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INPUT_FILE: str = "results/scored_candidates.csv"
FASTA_FILE: str = "data/raw/surface_proteins.fasta"
VALIDATED_FASTA_FILE: str = "data/validated_sequences.fasta"
OUTPUT_DIR: str = "results/construct"
COMPARISON_REPORT: str = "variant_comparison.md"

# ---------------------------------------------------------------------------
# Validated antigens (same list used by module 04)
# ---------------------------------------------------------------------------

VALIDATED_ANTIGENS: list[dict] = [
    {
        "gene_id": "LiHyp1",
        "gene_name": "LiHyp1",
        "final_score": 0.95,
    },
    {
        "gene_id": "A2",
        "gene_name": "A2",
        "final_score": 0.92,
    },
    {
        "gene_id": "LBSap_antigens",
        "gene_name": "LBSap antigens",
        "final_score": 0.90,
    },
    {
        "gene_id": "Lutzomyia_longipalpis_proteins",
        "gene_name": "Lutzomyia longipalpis proteins",
        "final_score": 0.88,
    },
    {
        "gene_id": "KMP-11",
        "gene_name": "KMP-11",
        "final_score": 0.85,
    },
    {
        "gene_id": "LiESP_Q",
        "gene_name": "LiESP/Q",
        "final_score": 0.83,
    },
    {
        "gene_id": "LACK",
        "gene_name": "LACK",
        "final_score": 0.82,
    },
    {
        "gene_id": "HSP70_HSP83",
        "gene_name": "HSP70/HSP83",
        "final_score": 0.80,
    },
]

VALIDATED_GENE_IDS: set[str] = {a["gene_id"] for a in VALIDATED_ANTIGENS}

# ---------------------------------------------------------------------------
# Variant definitions
# ---------------------------------------------------------------------------

VARIANTS: list[dict[str, Any]] = [
    {
        "variant_id": "A",
        "description": "Computational candidates only, CTL epitopes, L7L12 adjuvant",
        "include_validated": False,
        "include_htl": False,
        "signal_peptide": "tPA",
        "adjuvant": "L7L12",
    },
    {
        "variant_id": "B",
        "description": "Validated + computational, CTL + HTL epitopes, L7L12 adjuvant",
        "include_validated": True,
        "include_htl": True,
        "signal_peptide": "tPA",
        "adjuvant": "L7L12",
    },
    {
        "variant_id": "C",
        "description": "Validated + computational, CTL + HTL epitopes, RS09 adjuvant",
        "include_validated": True,
        "include_htl": True,
        "signal_peptide": "tPA",
        "adjuvant": "RS09",
    },
]


# ---------------------------------------------------------------------------
# Helper: load candidates with optional validated-antigen merging
# ---------------------------------------------------------------------------


def _load_candidates(include_validated: bool) -> list[dict]:
    """Load scored candidates from CSV, optionally merging validated antigens.

    Args:
        include_validated: When ``True``, append validated antigens from the
            curated list so they are considered during epitope selection.

    Returns:
        List of candidate dicts with at least ``gene_id``, ``gene_name``,
        and ``final_score`` keys.

    Raises:
        FileNotFoundError: If the scored-candidates CSV does not exist.
    """
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

    logger.info("Loaded %d scored candidate(s) from %s.", len(candidates), INPUT_FILE)

    if include_validated:
        existing_ids = {c["gene_id"] for c in candidates}
        added = 0
        for antigen in VALIDATED_ANTIGENS:
            if antigen["gene_id"] not in existing_ids:
                candidates.append(antigen)
                added += 1
        if added:
            logger.info("Merged %d validated antigen(s) into candidate pool.", added)

    if not candidates:
        logger.error("No candidates available. Cannot design construct.")
        raise ValueError("No candidates available for construct design.")

    return candidates


def _load_sequences(include_validated: bool) -> dict[str, str]:
    """Load protein sequences from FASTA files.

    Always loads the main surface-proteins FASTA.  When
    *include_validated* is ``True``, also loads the validated-sequences
    FASTA if it exists.

    Args:
        include_validated: Whether to include validated-sequence FASTA.

    Returns:
        Mapping of gene_id to amino-acid sequence.
    """
    sequences: dict[str, str] = {}

    if Path(FASTA_FILE).exists():
        sequences = _load_sequences_from_fasta(FASTA_FILE)
        logger.info("Loaded %d sequence(s) from %s.", len(sequences), FASTA_FILE)
    else:
        logger.warning(
            "FASTA file %s not found. Epitope prediction may be limited.",
            FASTA_FILE,
        )

    if include_validated and Path(VALIDATED_FASTA_FILE).exists():
        validated_seqs = _load_sequences_from_fasta(VALIDATED_FASTA_FILE)
        sequences.update(validated_seqs)
        logger.info(
            "Loaded %d validated sequence(s) from %s.",
            len(validated_seqs),
            VALIDATED_FASTA_FILE,
        )

    return sequences


# ---------------------------------------------------------------------------
# Single-variant pipeline
# ---------------------------------------------------------------------------


def _design_single_variant(
    variant: dict[str, Any],
    ctl_epitopes: list[SelectedEpitope],
    htl_epitopes: list[SelectedEpitope] | None,
) -> dict[str, Any]:
    """Run the construct pipeline for a single variant definition.

    Uses pre-computed CTL and HTL epitopes (shared across variants that
    use the same candidate/sequence pool) to avoid redundant IEDB calls.

    Args:
        variant: One entry from :data:`VARIANTS`.
        ctl_epitopes: Pre-selected CTL epitopes for this candidate pool.
        htl_epitopes: Pre-selected HTL epitopes, or ``None`` if HTL is
            not included in this variant.

    Returns:
        Dictionary with variant metadata and computed properties.
    """
    variant_id = variant["variant_id"]
    adjuvant = variant["adjuvant"]
    signal_peptide = variant["signal_peptide"]
    include_htl = variant["include_htl"]

    variant_dir = str(Path(OUTPUT_DIR) / f"variant_{variant_id}")
    Path(variant_dir).mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("VARIANT %s: %s", variant_id, variant["description"])
    logger.info("=" * 60)

    # -- Combine epitopes for assembly ------------------------------------
    combined_epitopes = list(ctl_epitopes)
    if include_htl and htl_epitopes:
        combined_epitopes.extend(htl_epitopes)

    if not combined_epitopes:
        logger.error("No epitopes available for variant %s.", variant_id)
        raise RuntimeError(
            f"No epitopes available for variant {variant_id}."
        )

    # -- Assemble protein construct ----------------------------------------
    # assemble_construct expects SelectedEpitope objects with a .peptide
    # attribute.  HTL epitopes use GPGPG linker; the current assemble_construct
    # uses CTL linker (AAY).  We pass htl_epitopes separately if the function
    # supports it, otherwise we assemble manually.
    if include_htl and htl_epitopes and hasattr(_mod06, "assemble_construct"):
        # Check if assemble_construct accepts htl_epitopes parameter
        import inspect

        sig = inspect.signature(assemble_construct)
        if "htl_epitopes" in sig.parameters:
            construct = assemble_construct(
                ctl_epitopes,
                signal_peptide=signal_peptide,
                adjuvant=adjuvant,
                htl_epitopes=htl_epitopes,
            )
        else:
            construct = _assemble_with_htl(
                ctl_epitopes, htl_epitopes, signal_peptide, adjuvant
            )
    else:
        construct = assemble_construct(
            ctl_epitopes,
            signal_peptide=signal_peptide,
            adjuvant=adjuvant,
        )

    # -- Physicochemical analysis ------------------------------------------
    logger.info("Computing physicochemical properties for variant %s ...", variant_id)
    physicochemical = compute_physicochemical(construct)

    # -- Codon optimisation ------------------------------------------------
    logger.info("Running codon optimisation for variant %s ...", variant_id)
    cds = reverse_translate_optimized(construct)
    gc = calculate_gc_content(cds)

    # -- mRNA assembly -----------------------------------------------------
    mrna = assemble_mrna(cds)

    # -- Safety predictions ------------------------------------------------
    logger.info("Predicting antigenicity for variant %s ...", variant_id)
    vaxijen_score = predict_antigenicity(construct)

    logger.info("Predicting allergenicity for variant %s ...", variant_id)
    allertop_result = predict_allergenicity(construct)

    # -- Restriction site accounting ---------------------------------------
    naive_cds = "".join(get_optimal_codon(aa) for aa in construct)
    naive_site_count = sum(naive_cds.count(site) for site in RESTRICTION_SITES)
    final_site_count = sum(cds.count(site) for site in RESTRICTION_SITES)
    sites_removed = max(0, naive_site_count - final_site_count)

    # -- Build identity card -----------------------------------------------
    ctl_count = len(ctl_epitopes)
    htl_count = len(htl_epitopes) if (include_htl and htl_epitopes) else 0

    card = ConstructCard(
        signal_peptide_name=signal_peptide,
        adjuvant_name=adjuvant,
        epitope_count=ctl_count + htl_count,
        protein_length=len(construct),
        mrna_length=len(mrna),
        gc_content=gc,
        molecular_weight=physicochemical.get("molecular_weight", 0.0),
        isoelectric_point=physicochemical.get("isoelectric_point", 0.0),
        instability_index=physicochemical.get("instability_index", 0.0),
        gravy=physicochemical.get("gravy", 0.0),
        aromaticity=physicochemical.get("aromaticity", 0.0),
        vaxijen_score=vaxijen_score,
        allertop_result=allertop_result,
        restriction_sites_removed=sites_removed,
        epitopes=[
            {
                "peptide": ep.peptide,
                "gene_id": ep.gene_id,
                "gene_name": ep.gene_name,
                "allele": ep.allele,
                "ic50": ep.ic50,
                "rank": ep.rank,
                "position": ep.position,
            }
            for ep in combined_epitopes
        ],
    )

    # -- Write per-variant outputs -----------------------------------------
    write_outputs(construct, mrna, combined_epitopes, physicochemical, card, variant_dir)

    # -- Summary log -------------------------------------------------------
    logger.info("Variant %s complete: %d aa protein, %d nt mRNA, %d CTL + %d HTL epitopes.",
                variant_id, len(construct), len(mrna), ctl_count, htl_count)

    return {
        "variant_id": variant_id,
        "description": variant["description"],
        "include_validated": variant["include_validated"],
        "adjuvant": adjuvant,
        "signal_peptide": signal_peptide,
        "ctl_count": ctl_count,
        "htl_count": htl_count,
        "protein_length": len(construct),
        "mrna_length": len(mrna),
        "gc_content": gc,
        "molecular_weight": physicochemical.get("molecular_weight", 0.0),
        "isoelectric_point": physicochemical.get("isoelectric_point", 0.0),
        "instability_index": physicochemical.get("instability_index", 0.0),
        "gravy": physicochemical.get("gravy", 0.0),
        "vaxijen_score": vaxijen_score,
        "allertop_result": allertop_result,
        "output_dir": variant_dir,
    }


def _assemble_with_htl(
    ctl_epitopes: list[SelectedEpitope],
    htl_epitopes: list[SelectedEpitope],
    signal_peptide: str,
    adjuvant: str,
) -> str:
    """Assemble a construct with both CTL and HTL epitope cassettes.

    Fallback used when :func:`assemble_construct` does not yet accept
    an ``htl_epitopes`` parameter.

    Layout::

        [Signal peptide] + [Adjuvant] + EAAAK +
        [CTL1] + AAY + [CTL2] + ... +
        EAAAK +
        [HTL1] + GPGPG + [HTL2] + ...

    Args:
        ctl_epitopes: Selected CTL epitopes.
        htl_epitopes: Selected HTL epitopes.
        signal_peptide: Key into :data:`SIGNAL_PEPTIDES`.
        adjuvant: Key into :data:`ADJUVANTS`.

    Returns:
        Full amino-acid sequence of the vaccine construct.
    """
    parts: list[str] = [
        SIGNAL_PEPTIDES[signal_peptide],
        ADJUVANTS[adjuvant],
        LINKERS["adjuvant_to_epitopes"],
    ]

    # CTL cassette with AAY linkers
    for i, ep in enumerate(ctl_epitopes):
        if i > 0:
            parts.append(LINKERS["ctl"])
        parts.append(ep.peptide)

    # Bridge between CTL and HTL cassettes
    if htl_epitopes:
        parts.append(LINKERS["adjuvant_to_epitopes"])
        for i, ep in enumerate(htl_epitopes):
            if i > 0:
                parts.append(LINKERS["htl"])
            parts.append(ep.peptide)

    construct = "".join(parts)
    logger.info(
        "Assembled construct with HTL: %d aa (%s signal + %s adjuvant + %d CTL + %d HTL).",
        len(construct),
        signal_peptide,
        adjuvant,
        len(ctl_epitopes),
        len(htl_epitopes),
    )
    return construct


# ---------------------------------------------------------------------------
# Comparison report generation
# ---------------------------------------------------------------------------


def _generate_comparison_report(
    results: list[dict[str, Any]],
    output_dir: str = OUTPUT_DIR,
) -> str:
    """Generate a Markdown comparison report for all variants.

    Args:
        results: List of result dicts from :func:`_design_single_variant`.
        output_dir: Directory where the report is written.

    Returns:
        Path to the written comparison report.
    """
    today = datetime.now().strftime("%Y-%m-%d")
    lines: list[str] = []

    lines.append("# Marley -- Vaccine Construct Variant Comparison")
    lines.append("")
    lines.append(f"**Generated:** {today}")
    lines.append("")

    # -- Summary table -----------------------------------------------------
    header = ["Property"] + [f"Variant {r['variant_id']}" for r in results]
    separator = ["-" * max(8, len(h)) for h in header]

    def _row(label: str, values: list[str]) -> str:
        cells = [label] + values
        return "| " + " | ".join(cells) + " |"

    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(separator) + " |")

    lines.append(_row("Description", [r["description"] for r in results]))
    lines.append(_row("CTL epitopes", [str(r["ctl_count"]) for r in results]))
    lines.append(_row("HTL epitopes", [str(r["htl_count"]) for r in results]))
    lines.append(_row(
        "Protein length (aa)",
        [str(r["protein_length"]) for r in results],
    ))
    lines.append(_row(
        "mRNA length (nt)",
        [str(r["mrna_length"]) for r in results],
    ))
    lines.append(_row(
        "MW (kDa)",
        [f"{r['molecular_weight'] / 1000:.1f}" for r in results],
    ))
    lines.append(_row(
        "pI",
        [f"{r['isoelectric_point']:.2f}" for r in results],
    ))
    lines.append(_row(
        "Instability index",
        [f"{r['instability_index']:.2f}" for r in results],
    ))
    lines.append(_row(
        "GC content (%)",
        [f"{r['gc_content'] * 100:.1f}" for r in results],
    ))
    lines.append(_row(
        "Adjuvant",
        [
            f"{r['adjuvant']} ({len(ADJUVANTS[r['adjuvant']])} aa)"
            for r in results
        ],
    ))
    lines.append(_row(
        "Validated antigens",
        ["Yes" if r["include_validated"] else "No" for r in results],
    ))

    # VaxiJen row
    vaxijen_values: list[str] = []
    for r in results:
        if r["vaxijen_score"] is not None:
            vaxijen_values.append(f"{r['vaxijen_score']:.4f}")
        else:
            vaxijen_values.append("N/A")
    lines.append(_row("VaxiJen score", vaxijen_values))

    # AllerTOP row
    allertop_values: list[str] = []
    for r in results:
        if r["allertop_result"] is not None:
            allertop_values.append(r["allertop_result"])
        else:
            allertop_values.append("N/A")
    lines.append(_row("AllerTOP", allertop_values))

    lines.append("")

    # -- Stability assessment ----------------------------------------------
    lines.append("## Stability Assessment")
    lines.append("")
    for r in results:
        ii = r["instability_index"]
        stability = "stable" if ii < 40 else "unstable"
        lines.append(
            f"- **Variant {r['variant_id']}:** instability index = "
            f"{ii:.2f} ({stability}, threshold < 40)"
        )
    lines.append("")

    # -- Output locations --------------------------------------------------
    lines.append("## Output Locations")
    lines.append("")
    for r in results:
        lines.append(f"- **Variant {r['variant_id']}:** `{r['output_dir']}/`")
    lines.append("")

    # -- Disclaimer --------------------------------------------------------
    lines.append("## Disclaimer")
    lines.append("")
    lines.append(
        "This is a **computational vaccine design comparison** generated by "
        "the Marley reverse vaccinology pipeline. All variants are *in silico* "
        "predictions and must be validated experimentally before any "
        "translational application."
    )
    lines.append("")

    # -- Write report ------------------------------------------------------
    report_path = Path(output_dir) / COMPARISON_REPORT
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as fh:
        fh.write("\n".join(lines))

    logger.info("Comparison report written to %s.", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------


def design_all_variants(force: bool = False) -> str:
    """Design all construct variants and produce comparison report.

    For each variant:
    1. Load candidates (filter by include_validated)
    2. Load sequences from FASTA + validated_sequences.fasta
    3. Select CTL epitopes (IEDB MHC-I)
    4. Select HTL epitopes if include_htl (IEDB MHC-II)
    5. Assemble construct with specified adjuvant
    6. Codon optimize, compute physicochemical, check safety
    7. Write outputs to results/construct/variant_{id}/

    Then generate comparison report.

    Args:
        force: Re-run even when output directories already contain results.

    Returns:
        Path to comparison report.
    """
    comparison_path = Path(OUTPUT_DIR) / COMPARISON_REPORT
    if comparison_path.exists() and not force:
        logger.info(
            "Comparison report already exists at %s. Use force=True to re-run.",
            comparison_path,
        )
        return str(comparison_path)

    logger.info("=" * 60)
    logger.info("MULTI-VARIANT CONSTRUCT GENERATION")
    logger.info("=" * 60)

    # ------------------------------------------------------------------
    # Pre-compute epitopes to avoid redundant IEDB API calls.
    # We need two candidate pools: computational-only and computational+validated.
    # ------------------------------------------------------------------

    # Pool 1: computational only
    candidates_comp = _load_candidates(include_validated=False)
    sequences_comp = _load_sequences(include_validated=False)

    # Pool 2: computational + validated
    candidates_full = _load_candidates(include_validated=True)
    sequences_full = _load_sequences(include_validated=True)

    # -- CTL epitopes (cached per pool) ------------------------------------
    logger.info("Selecting CTL epitopes for computational-only pool ...")
    ctl_comp = select_epitopes(candidates_comp, sequences_comp)

    logger.info("Selecting CTL epitopes for full (validated + computational) pool ...")
    ctl_full = select_epitopes(candidates_full, sequences_full)

    # -- HTL epitopes (only needed for variants B and C) -------------------
    htl_full: list[SelectedEpitope] | None = None
    if any(v["include_htl"] for v in VARIANTS):
        if select_htl_epitopes is not None:
            logger.info("Selecting HTL epitopes for full pool ...")
            htl_full = select_htl_epitopes(candidates_full, sequences_full)
        else:
            logger.warning(
                "select_htl_epitopes not yet available in 06_construct. "
                "HTL epitopes will be empty for this run."
            )
            htl_full = []

    # ------------------------------------------------------------------
    # Build each variant
    # ------------------------------------------------------------------
    results: list[dict[str, Any]] = []

    for variant in VARIANTS:
        variant_id = variant["variant_id"]
        include_validated = variant["include_validated"]
        include_htl = variant["include_htl"]

        # Pick the right pre-computed epitopes
        ctl = ctl_full if include_validated else ctl_comp
        htl = htl_full if include_htl else None

        if not ctl:
            logger.error(
                "No CTL epitopes for variant %s. Skipping.", variant_id
            )
            continue

        result = _design_single_variant(variant, ctl, htl)
        results.append(result)

    if not results:
        logger.error("No variants could be designed.")
        raise RuntimeError("All variant designs failed.")

    # ------------------------------------------------------------------
    # Generate comparison report
    # ------------------------------------------------------------------
    report_path = _generate_comparison_report(results)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("MULTI-VARIANT GENERATION COMPLETE")
    logger.info("=" * 60)
    logger.info("  Variants designed: %d / %d", len(results), len(VARIANTS))
    for r in results:
        logger.info(
            "  Variant %s: %d aa, %d epitopes (%d CTL + %d HTL), %s adjuvant",
            r["variant_id"],
            r["protein_length"],
            r["ctl_count"] + r["htl_count"],
            r["ctl_count"],
            r["htl_count"],
            r["adjuvant"],
        )
    logger.info("  Comparison report: %s", report_path)
    logger.info("=" * 60)

    return report_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Step 6b: Generate multiple mRNA vaccine construct variants "
            "for comparison."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if the output directory already contains results.",
    )
    args = parser.parse_args()

    report = design_all_variants(force=args.force)
    print(f"Comparison report: {report}")
