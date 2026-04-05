"""Reconstruct the optimized mRNA vaccine with multiple ordering strategies.

Loads safe epitopes from the v4 optimization pipeline, selects the best
protein adjuvant, builds four construct variants (one per epitope ordering
strategy), computes physicochemical properties, and generates a comparison
report against the v1 baseline.

Usage:
    python -m pipeline.11_reconstruct
    python -m pipeline.11_reconstruct --force
    python -m pipeline.11_reconstruct --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("reconstruct")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SAFETY_REPORT_PATH: Final[str] = "results/v4_optimization/safety_report.json"
ADJUVANT_RANKING_PATH: Final[str] = (
    "results/v4_optimization/adjuvant_ranking.json"
)
V1_CONSTRUCT_CARD_PATH: Final[str] = "results/construct/construct_card.json"
V1_CONSTRUCT_CARD_FALLBACKS: Final[list[str]] = [
    "results/construct/variant_A/construct_card.json",
    "results/construct/variant_B/construct_card.json",
    "results/construct/variant_C/construct_card.json",
]
OUTPUT_DIR: Final[str] = "results/v4_optimization"

ORDERING_STRATEGIES: Final[list[str]] = [
    "ic50_ascending",
    "ic50_descending",
    "gene_clustered",
    "alternating_type",
]

IMPROVEMENT_RATIO_THRESHOLD: Final[float] = 1.1

# tPA signal peptide (same as v1 construct)
TPA_SIGNAL_PEPTIDE: Final[str] = "MDAMKRGLCCVLLLCGAVFVSAS"

# Linkers
LINKER_ADJUVANT: Final[str] = "EAAAK"
LINKER_CTL: Final[str] = "AAY"
LINKER_HTL: Final[str] = "GPGPG"

# Dog codon table -- most frequent codon per amino acid for Canis lupus
# familiaris, derived from the Codon Usage Database (Kazusa).
DOG_CODONS: Final[dict[str, str]] = {
    "M": "ATG", "F": "TTC", "L": "CTG", "S": "AGC", "Y": "TAC",
    "C": "TGC", "W": "TGG", "P": "CCC", "H": "CAC", "Q": "CAG",
    "R": "AGG", "I": "ATC", "T": "ACC", "N": "AAC", "K": "AAG",
    "V": "GTG", "A": "GCC", "D": "GAC", "E": "GAG", "G": "GGC",
    "*": "TGA",
}

# mRNA cassette components
FIVE_PRIME_UTR: Final[str] = (
    "GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGAGCCACC"
)
THREE_PRIME_UTR_SINGLE: Final[str] = (
    "GCUGGAGCCUCGGUGGCCAUGCUUCUUGCCCCUUGGGCCUCCCCCCAG"
    "CCCCUCCUCCCCUUCCUGCACCCGUACCCCCGUGGUCUUUGAAUAAAGUCUGA"
)
THREE_PRIME_UTR: Final[str] = THREE_PRIME_UTR_SINGLE * 2
POLY_A_TAIL: Final[str] = "A" * 120

# Approximate amino acid molecular weights (Da), including water loss
# during peptide-bond formation.
AA_WEIGHTS: Final[dict[str, float]] = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "E": 147.13, "Q": 146.15, "G": 75.03, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}
WATER_WEIGHT: Final[float] = 18.015


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> dict | list:
    """Load and return a JSON file."""
    with open(path) as fh:
        return json.load(fh)


def _write_json(data: object, path: Path) -> None:
    """Write *data* as pretty-printed JSON to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
    logger.info("Wrote %s", path)


def _write_fasta(header: str, sequence: str, path: Path) -> None:
    """Write a single-record FASTA file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            fh.write(sequence[i : i + 80] + "\n")
    logger.info("Wrote %s", path)


def _write_text(text: str, path: Path) -> None:
    """Write plain text to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    logger.info("Wrote %s", path)


def _find_v1_construct_card() -> Path | None:
    """Locate the v1 construct_card.json, returning None if not found."""
    primary = Path(V1_CONSTRUCT_CARD_PATH)
    if primary.exists():
        return primary
    for fallback in V1_CONSTRUCT_CARD_FALLBACKS:
        p = Path(fallback)
        if p.exists():
            return p
    return None


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------


def _load_safe_epitopes(safety_path: Path) -> list[dict]:
    """Load epitopes from the safety report and keep only safe ones.

    For each safe epitope the optimized sequence is used when the
    improvement ratio exceeds the threshold; otherwise the original
    sequence is kept.

    Returns:
        List of epitope dicts with a resolved ``sequence`` key.
    """
    data = _load_json(safety_path)
    epitopes_raw: list[dict] = data if isinstance(data, list) else data.get("epitopes", [])

    safe: list[dict] = []
    for ep in epitopes_raw:
        if not ep.get("is_safe", False):
            continue

        ratio = ep.get("improvement_ratio", 1.0)
        if ratio > IMPROVEMENT_RATIO_THRESHOLD:
            resolved_seq = ep.get("optimized_sequence", ep.get("original_sequence", ep.get("peptide", "")))
        else:
            resolved_seq = ep.get("original_sequence", ep.get("peptide", ""))

        safe.append({
            **ep,
            "sequence": resolved_seq,
        })

    return safe


def _select_best_adjuvant(adjuvant_path: Path) -> dict:
    """Select the protein adjuvant with the highest VaxiJen score.

    Returns:
        Dict with at least ``name``, ``sequence``, and ``vaxijen_score``.

    Raises:
        ValueError: If no protein adjuvant is found in the ranking.
    """
    data = _load_json(adjuvant_path)
    candidates: list[dict] = data if isinstance(data, list) else data.get("adjuvants", [])

    protein_candidates = [c for c in candidates if c.get("is_protein", False)]
    if not protein_candidates:
        raise ValueError(
            "No protein adjuvant with is_protein=True found in "
            f"{adjuvant_path}"
        )

    best = max(protein_candidates, key=lambda c: c.get("vaxijen_score", 0.0))
    return best


# ---------------------------------------------------------------------------
# Epitope ordering strategies
# ---------------------------------------------------------------------------


def _order_ic50_ascending(epitopes: list[dict]) -> list[dict]:
    """Sort epitopes by IC50 ascending (best binders first)."""
    return sorted(epitopes, key=lambda e: e.get("ic50", 9999.0))


def _order_ic50_descending(epitopes: list[dict]) -> list[dict]:
    """Sort epitopes by IC50 descending (best binders last)."""
    return sorted(epitopes, key=lambda e: e.get("ic50", 9999.0), reverse=True)


def _order_gene_clustered(epitopes: list[dict]) -> list[dict]:
    """Group epitopes by gene, then sort within each group by IC50."""
    return sorted(
        epitopes,
        key=lambda e: (e.get("gene_id", ""), e.get("ic50", 9999.0)),
    )


def _order_alternating_type(epitopes: list[dict]) -> list[dict]:
    """Interleave CTL and HTL epitopes in alternating order."""
    ctl = [e for e in epitopes if e.get("epitope_type", "CTL") == "CTL"]
    htl = [e for e in epitopes if e.get("epitope_type", "CTL") == "HTL"]

    ctl.sort(key=lambda e: e.get("ic50", 9999.0))
    htl.sort(key=lambda e: e.get("ic50", 9999.0))

    merged: list[dict] = []
    ci, hi = 0, 0
    while ci < len(ctl) or hi < len(htl):
        if ci < len(ctl):
            merged.append(ctl[ci])
            ci += 1
        if hi < len(htl):
            merged.append(htl[hi])
            hi += 1

    return merged


ORDERING_FUNCTIONS: Final[dict[str, object]] = {
    "ic50_ascending": _order_ic50_ascending,
    "ic50_descending": _order_ic50_descending,
    "gene_clustered": _order_gene_clustered,
    "alternating_type": _order_alternating_type,
}

ORDERING_NOTES: Final[dict[str, str]] = {
    "ic50_ascending": "Best binders first",
    "ic50_descending": "Best binders last",
    "gene_clustered": "Same-gene epitopes together",
    "alternating_type": "CTL-HTL alternating",
}


# ---------------------------------------------------------------------------
# Construct assembly
# ---------------------------------------------------------------------------


def _split_epitopes_by_type(
    epitopes: list[dict],
) -> tuple[list[dict], list[dict]]:
    """Separate epitopes into CTL and HTL lists, preserving order."""
    ctl = [e for e in epitopes if e.get("epitope_type", "CTL") == "CTL"]
    htl = [e for e in epitopes if e.get("epitope_type", "CTL") == "HTL"]
    return ctl, htl


def _assemble_protein(
    adjuvant_seq: str,
    ordered_epitopes: list[dict],
) -> str:
    """Assemble the multi-epitope protein construct.

    Structure:
        [tPA] + [adjuvant] + EAAAK + [CTL epitopes with AAY] + GPGPG +
        [HTL epitopes with GPGPG]
    """
    ctl, htl = _split_epitopes_by_type(ordered_epitopes)

    parts: list[str] = [TPA_SIGNAL_PEPTIDE, adjuvant_seq, LINKER_ADJUVANT]

    # CTL block: epitopes joined by AAY linker
    ctl_seqs = [e["sequence"] for e in ctl]
    if ctl_seqs:
        parts.append(LINKER_CTL.join(ctl_seqs))

    # Transition linker between CTL and HTL blocks
    parts.append(LINKER_HTL)

    # HTL block: epitopes joined by GPGPG linker
    htl_seqs = [e["sequence"] for e in htl]
    if htl_seqs:
        parts.append(LINKER_HTL.join(htl_seqs))

    return "".join(parts)


def _codon_optimize(protein: str) -> str:
    """Translate protein to DNA using the dog-preferred codon table.

    Args:
        protein: Amino acid sequence (single-letter codes).

    Returns:
        DNA coding sequence (no stop codon appended).
    """
    codons: list[str] = []
    for aa in protein:
        codon = DOG_CODONS.get(aa)
        if codon is None:
            raise ValueError(f"Unknown amino acid '{aa}' in protein sequence.")
        codons.append(codon)
    return "".join(codons)


def _build_mrna(orf_dna: str) -> str:
    """Assemble the full mRNA cassette from the coding DNA.

    Converts DNA to RNA (T -> U) and wraps in UTRs + poly(A).

    Structure:
        5'UTR + Kozak(in UTR) + AUG + ORF(minus start) + stop + 3'UTR + poly(A)
    """
    # The ORF already starts with ATG from the tPA signal peptide.
    orf_rna = orf_dna.replace("T", "U")
    stop_rna = DOG_CODONS["*"].replace("T", "U")

    return FIVE_PRIME_UTR + orf_rna + stop_rna + THREE_PRIME_UTR + POLY_A_TAIL


# ---------------------------------------------------------------------------
# Physicochemical properties
# ---------------------------------------------------------------------------


def _compute_molecular_weight(protein: str) -> float:
    """Approximate molecular weight in Daltons.

    Sum of individual AA weights minus (n-1) water molecules lost during
    peptide bond formation.
    """
    if not protein:
        return 0.0
    total = sum(AA_WEIGHTS.get(aa, 0.0) for aa in protein)
    water_loss = (len(protein) - 1) * WATER_WEIGHT
    return total - water_loss


def _compute_instability_index(protein: str) -> float:
    """Simplified instability index based on DKGS dipeptide motifs.

    This is a rough heuristic; a full Guruprasad implementation would use
    the 400-element DIWV weight table.  Here we count occurrences of the
    four destabilizing dipeptides (DK, GS, DG, KS) and scale by length.
    """
    if len(protein) < 2:
        return 0.0

    destabilizing = {"DK", "GS", "DG", "KS"}
    count = sum(
        1
        for i in range(len(protein) - 1)
        if protein[i : i + 2] in destabilizing
    )
    return (10.0 / len(protein)) * count * 100.0


def _compute_gc_content(mrna: str) -> float:
    """GC content of an RNA sequence as a fraction in [0, 1]."""
    if not mrna:
        return 0.0
    gc = sum(1 for nt in mrna if nt in "GCgc")
    return gc / len(mrna)


def _compute_properties(protein: str, mrna: str) -> dict:
    """Compute a dict of physicochemical and sequence-level properties."""
    mw = _compute_molecular_weight(protein)
    return {
        "protein_length_aa": len(protein),
        "mrna_length_nt": len(mrna),
        "molecular_weight_da": round(mw, 2),
        "molecular_weight_kda": round(mw / 1000.0, 2),
        "gc_content": round(_compute_gc_content(mrna), 4),
        "gc_content_pct": round(_compute_gc_content(mrna) * 100.0, 2),
        "instability_index": round(_compute_instability_index(protein), 2),
    }


# ---------------------------------------------------------------------------
# Comparison report
# ---------------------------------------------------------------------------


def _load_v1_properties() -> dict | None:
    """Try to load the v1 construct card and extract comparable properties."""
    card_path = _find_v1_construct_card()
    if card_path is None:
        logger.warning("v1 construct card not found; comparison will be partial.")
        return None

    card = _load_json(card_path)
    return {
        "epitope_count": card.get("epitope_count", 0),
        "protein_length_aa": card.get("protein_length", 0),
        "mrna_length_nt": card.get("mrna_length", 0),
        "gc_content_pct": round(card.get("gc_content", 0.0) * 100.0, 2),
        "molecular_weight_kda": round(
            card.get("molecular_weight", 0.0) / 1000.0, 2
        ),
    }


def _format_change(old: float | int, new: float | int) -> str:
    """Format a change value with sign prefix."""
    diff = new - old
    if isinstance(diff, float):
        return f"{diff:+.2f}"
    return f"{diff:+d}"


def _build_comparison_report(
    safe_epitopes: list[dict],
    all_epitopes_raw: list[dict],
    adjuvant: dict,
    strategy_results: dict[str, dict],
    best_strategy: str,
) -> str:
    """Generate the Markdown comparison report.

    Args:
        safe_epitopes: Epitopes that passed safety and were included.
        all_epitopes_raw: All epitopes from the safety report (before filter).
        adjuvant: Selected adjuvant dict.
        strategy_results: Mapping of strategy name to properties dict.
        best_strategy: Name of the recommended strategy.

    Returns:
        The full Markdown report as a string.
    """
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    total_raw = len(all_epitopes_raw)
    safe_count = len(safe_epitopes)
    unsafe_count = total_raw - safe_count
    optimized_count = sum(
        1 for e in safe_epitopes
        if e.get("improvement_ratio", 1.0) > IMPROVEMENT_RATIO_THRESHOLD
    )

    v1 = _load_v1_properties()
    best_props = strategy_results[best_strategy]

    lines: list[str] = []
    lines.append("# Marley v4 -- Vaccine Optimization Report")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")

    # -- Epitope Optimization Summary --
    lines.append("## Epitope Optimization Summary")
    lines.append(f"- Original epitopes: {total_raw}")
    lines.append(f"- Optimized (improved IC50): {optimized_count}")
    lines.append(f"- Removed (unsafe): {unsafe_count}")
    lines.append(f"- Safe and optimized: {safe_count}")
    lines.append("")

    # -- Adjuvant Selection --
    adj_name = adjuvant.get("name", "Unknown")
    adj_score = adjuvant.get("vaxijen_score", 0.0)
    lines.append("## Adjuvant Selection")
    lines.append(
        f"- Best protein adjuvant: {adj_name} (VaxiJen={adj_score:.4f})"
    )
    lines.append(
        "- Recommended co-adjuvant: CpG ODN 2006 (Th1 bias for dogs)"
    )
    lines.append("")

    # -- Construct Comparison --
    lines.append("## Construct Comparison")
    lines.append("")

    if v1 is not None:
        lines.append(
            "| Metric | Original (v1) | Optimized (v4 best) | Change |"
        )
        lines.append("|--------|--------------|--------------------|---------| ")
        lines.append(
            f"| Epitopes | {v1['epitope_count']} | {safe_count} "
            f"| {_format_change(v1['epitope_count'], safe_count)} |"
        )
        lines.append(
            f"| Protein length | {v1['protein_length_aa']} aa "
            f"| {best_props['protein_length_aa']} aa "
            f"| {_format_change(v1['protein_length_aa'], best_props['protein_length_aa'])} aa |"
        )
        lines.append(
            f"| mRNA length | {v1['mrna_length_nt']} nt "
            f"| {best_props['mrna_length_nt']} nt "
            f"| {_format_change(v1['mrna_length_nt'], best_props['mrna_length_nt'])} nt |"
        )
        lines.append(
            f"| GC content | {v1['gc_content_pct']}% "
            f"| {best_props['gc_content_pct']}% "
            f"| {_format_change(v1['gc_content_pct'], best_props['gc_content_pct'])}% |"
        )
        lines.append(
            f"| Molecular weight | {v1['molecular_weight_kda']} kDa "
            f"| {best_props['molecular_weight_kda']} kDa "
            f"| {_format_change(v1['molecular_weight_kda'], best_props['molecular_weight_kda'])} kDa |"
        )
    else:
        lines.append(
            "| Metric | Optimized (v4 best) |"
        )
        lines.append("|--------|---------------------|")
        lines.append(f"| Epitopes | {safe_count} |")
        lines.append(
            f"| Protein length | {best_props['protein_length_aa']} aa |"
        )
        lines.append(
            f"| mRNA length | {best_props['mrna_length_nt']} nt |"
        )
        lines.append(f"| GC content | {best_props['gc_content_pct']}% |")
        lines.append(
            f"| Molecular weight | {best_props['molecular_weight_kda']} kDa |"
        )

    lines.append("")

    # -- Epitope-by-Epitope Comparison --
    lines.append("## Epitope-by-Epitope Comparison")
    lines.append("")
    lines.append(
        "| # | Original | Optimized | IC50 orig | IC50 opt "
        "| Improvement | Safe? |"
    )
    lines.append(
        "|---|----------|-----------|-----------|----------|"
        "-------------|-------|"
    )
    for idx, ep in enumerate(all_epitopes_raw, start=1):
        orig_seq = ep.get("original_sequence", ep.get("peptide", ""))
        opt_seq = ep.get("optimized_sequence", orig_seq)
        ic50_orig = ep.get("original_ic50", ep.get("ic50", 0.0))
        ic50_opt = ep.get("optimized_ic50", ic50_orig)
        ratio = ep.get("improvement_ratio", 1.0)
        is_safe = ep.get("is_safe", False)

        lines.append(
            f"| {idx} | {orig_seq} | {opt_seq} "
            f"| {ic50_orig:.1f} | {ic50_opt:.1f} "
            f"| {ratio:.2f}x | {'Yes' if is_safe else 'No'} |"
        )

    lines.append("")

    # -- Ordering Strategy Comparison --
    lines.append("## Ordering Strategy Comparison")
    lines.append("")
    lines.append("| Strategy | Length | GC% | Notes |")
    lines.append("|----------|--------|-----|-------|")
    for strat in ORDERING_STRATEGIES:
        props = strategy_results[strat]
        note = ORDERING_NOTES.get(strat, "")
        lines.append(
            f"| {strat} | {props['protein_length_aa']} aa "
            f"| {props['gc_content_pct']}% | {note} |"
        )

    lines.append("")

    # -- Recommended Formulation --
    lines.append("## Recommended Formulation")
    lines.append(f"- mRNA construct: {best_strategy} variant")
    lines.append(f"- Adjuvant (fused): {adj_name}")
    lines.append("- Co-adjuvant: CpG ODN 2006 (co-administered)")
    lines.append("- Delivery: Lipid nanoparticle (LNP)")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def reconstruct_optimized(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the construct reconstruction pipeline.

    1. Load safe epitopes and best adjuvant from v4 optimization outputs.
    2. Build four construct variants (one per ordering strategy).
    3. Compute physicochemical properties for each variant.
    4. Generate a comparison report against the v1 baseline.
    5. Save all outputs to ``results/v4_optimization/``.

    Args:
        force: If True, overwrite existing output even if present.
        dry_run: If True, build constructs with simplified sequences
            (skip codon optimization) and write the comparison report
            with available data.

    Returns:
        Path to the comparison report Markdown file.
    """
    output_base = Path(OUTPUT_DIR)
    report_path = output_base / "comparison_report.md"

    if report_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            report_path,
        )
        return str(report_path)

    # ------------------------------------------------------------------
    # 1. Load inputs
    # ------------------------------------------------------------------
    safety_path = Path(SAFETY_REPORT_PATH)
    if not safety_path.exists():
        logger.error("Safety report not found at %s.", safety_path)
        sys.exit(1)

    adjuvant_path = Path(ADJUVANT_RANKING_PATH)
    if not adjuvant_path.exists():
        logger.error("Adjuvant ranking not found at %s.", adjuvant_path)
        sys.exit(1)

    all_epitopes_raw = _load_json(safety_path)
    if isinstance(all_epitopes_raw, dict):
        all_epitopes_raw = all_epitopes_raw.get("epitopes", [])

    safe_epitopes = _load_safe_epitopes(safety_path)
    if not safe_epitopes:
        logger.error("No safe epitopes found in %s.", safety_path)
        sys.exit(1)

    logger.info(
        "Loaded %d safe epitopes out of %d total.",
        len(safe_epitopes),
        len(all_epitopes_raw),
    )

    adjuvant = _select_best_adjuvant(adjuvant_path)
    adjuvant_seq = adjuvant.get("sequence", "")
    logger.info(
        "Selected adjuvant: %s (VaxiJen=%.4f)",
        adjuvant.get("name", "Unknown"),
        adjuvant.get("vaxijen_score", 0.0),
    )

    # ------------------------------------------------------------------
    # 2. Build constructs for each ordering strategy
    # ------------------------------------------------------------------
    strategy_results: dict[str, dict] = {}
    strategy_proteins: dict[str, str] = {}
    strategy_mrnas: dict[str, str] = {}

    for strategy in ORDERING_STRATEGIES:
        logger.info("Building construct with ordering: %s", strategy)

        order_fn = ORDERING_FUNCTIONS[strategy]
        ordered = order_fn(safe_epitopes)

        protein = _assemble_protein(adjuvant_seq, ordered)

        if dry_run:
            # In dry-run mode, skip codon optimization; use protein only.
            mrna = f"(dry-run: {len(protein) * 3} nt estimated)"
            props = {
                "protein_length_aa": len(protein),
                "mrna_length_nt": len(protein) * 3,
                "molecular_weight_da": round(
                    _compute_molecular_weight(protein), 2
                ),
                "molecular_weight_kda": round(
                    _compute_molecular_weight(protein) / 1000.0, 2
                ),
                "gc_content": 0.0,
                "gc_content_pct": 0.0,
                "instability_index": round(
                    _compute_instability_index(protein), 2
                ),
            }
        else:
            orf_dna = _codon_optimize(protein)
            mrna = _build_mrna(orf_dna)
            props = _compute_properties(protein, mrna)

        strategy_results[strategy] = props
        strategy_proteins[strategy] = protein
        strategy_mrnas[strategy] = mrna

        logger.info(
            "  %s: protein=%d aa, mRNA=%s nt, GC=%.1f%%",
            strategy,
            props["protein_length_aa"],
            props["mrna_length_nt"],
            props.get("gc_content_pct", 0.0),
        )

    # ------------------------------------------------------------------
    # 3. Select best strategy (lowest instability, then highest GC)
    # ------------------------------------------------------------------
    best_strategy = min(
        ORDERING_STRATEGIES,
        key=lambda s: (
            strategy_results[s]["instability_index"],
            -strategy_results[s].get("gc_content_pct", 0.0),
        ),
    )
    logger.info("Recommended ordering strategy: %s", best_strategy)

    # ------------------------------------------------------------------
    # 4. Save outputs for each strategy
    # ------------------------------------------------------------------
    for strategy in ORDERING_STRATEGIES:
        strat_dir = output_base / strategy
        strat_dir.mkdir(parents=True, exist_ok=True)

        props = strategy_results[strategy]
        protein = strategy_proteins[strategy]
        mrna = strategy_mrnas[strategy]

        # Construct card JSON
        card = {
            "strategy": strategy,
            "signal_peptide": "tPA",
            "adjuvant_name": adjuvant.get("name", "Unknown"),
            "adjuvant_vaxijen_score": adjuvant.get("vaxijen_score", 0.0),
            "epitope_count": len(safe_epitopes),
            "ctl_count": sum(
                1 for e in safe_epitopes
                if e.get("epitope_type", "CTL") == "CTL"
            ),
            "htl_count": sum(
                1 for e in safe_epitopes
                if e.get("epitope_type", "CTL") == "HTL"
            ),
            **props,
            "is_recommended": strategy == best_strategy,
            "generated_at": datetime.now().isoformat(),
        }
        _write_json(card, strat_dir / "construct_card.json")

        # Protein FASTA
        _write_fasta(
            header=(
                f"Marley_v4_{strategy} | "
                f"{props['protein_length_aa']} aa | "
                f"MW={props['molecular_weight_kda']} kDa"
            ),
            sequence=protein,
            path=strat_dir / "vaccine_construct.fasta",
        )

        # mRNA FASTA (skip in dry run since it is a placeholder string)
        if not dry_run:
            _write_fasta(
                header=(
                    f"Marley_v4_{strategy}_mRNA | "
                    f"{props['mrna_length_nt']} nt | "
                    f"GC={props['gc_content_pct']}%"
                ),
                sequence=mrna,
                path=strat_dir / "vaccine_mrna.fasta",
            )

    # ------------------------------------------------------------------
    # 5. Generate and save comparison report
    # ------------------------------------------------------------------
    report = _build_comparison_report(
        safe_epitopes=safe_epitopes,
        all_epitopes_raw=all_epitopes_raw,
        adjuvant=adjuvant,
        strategy_results=strategy_results,
        best_strategy=best_strategy,
    )
    _write_text(report, report_path)

    logger.info("Comparison report written to %s", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley v4 -- Reconstruct the optimized vaccine construct "
            "with multiple epitope ordering strategies and generate a "
            "comparison report."
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
            "Build constructs with simplified sequences (skip codon "
            "optimization) and write the report with available data."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = reconstruct_optimized(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
