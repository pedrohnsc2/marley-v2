"""Optimize mRNA delivery: structure analysis, CAI, LNP formulation.

Loads the vaccine mRNA from the construct stage, analyzes secondary
structure (via ViennaRNA when available, otherwise basic heuristics),
calculates the Codon Adaptation Index for Canis lupus familiaris,
identifies suboptimal codons, and recommends an LNP formulation with
dosing schedule.

Usage:
    python -m pipeline.15_delivery
    python -m pipeline.15_delivery --force
    python -m pipeline.15_delivery --dry-run
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

logger = get_logger("delivery")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MRNA_FASTA_PATH: Final[str] = "results/construct/vaccine_mrna.fasta"
OUTPUT_DIR: Final[str] = "results/delivery"

AVG_NUCLEOTIDE_MW: Final[float] = 330.0  # Da per nucleotide (average)
PHOSPHATE_CHARGE: Final[int] = -1  # One negative charge per phosphate

NP_RATIO: Final[str] = "6:1"
TARGET_PARTICLE_SIZE: Final[str] = "80-100 nm"
ENCAPSULATION_EFFICIENCY: Final[str] = ">90%"
DOSE_PER_DOG_UG: Final[str] = "50-100"
DOG_WEIGHT_KG: Final[int] = 30

CAI_TARGET: Final[float] = 0.80
WORST_CODONS_TO_REPORT: Final[int] = 10

HAIRPIN_MIN_LEN: Final[int] = 6
HAIRPIN_WINDOW: Final[int] = 50

LNP_COMPOSITION: Final[list[dict[str, str]]] = [
    {
        "component": "SM-102 (ionizable lipid)",
        "mol_pct": "50.0%",
        "role": "mRNA complexation + endosomal escape",
    },
    {
        "component": "Cholesterol",
        "mol_pct": "38.5%",
        "role": "Membrane stability",
    },
    {
        "component": "DSPC",
        "mol_pct": "10.0%",
        "role": "Structural lipid",
    },
    {
        "component": "DMG-PEG2000",
        "mol_pct": "1.5%",
        "role": "Steric stabilization, prevents aggregation",
    },
]

# Reference codon usage for dog (Canis lupus familiaris) from Kazusa.
# Values are relative synonymous codon usage fractions per amino acid.
DOG_CODON_USAGE: Final[dict[str, dict[str, float]]] = {
    "F": {"TTT": 0.43, "TTC": 0.57},
    "L": {
        "TTA": 0.07, "TTG": 0.12, "CTT": 0.12,
        "CTC": 0.20, "CTA": 0.07, "CTG": 0.42,
    },
    "I": {"ATT": 0.34, "ATC": 0.51, "ATA": 0.15},
    "V": {"GTT": 0.17, "GTC": 0.25, "GTA": 0.10, "GTG": 0.48},
    "S": {
        "TCT": 0.18, "TCC": 0.23, "TCA": 0.14,
        "TCG": 0.06, "AGT": 0.14, "AGC": 0.25,
    },
    "P": {"CCT": 0.29, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11},
    "T": {"ACT": 0.23, "ACC": 0.37, "ACA": 0.27, "ACG": 0.13},
    "A": {"GCT": 0.26, "GCC": 0.41, "GCA": 0.22, "GCG": 0.11},
    "Y": {"TAT": 0.42, "TAC": 0.58},
    "H": {"CAT": 0.40, "CAC": 0.60},
    "Q": {"CAA": 0.25, "CAG": 0.75},
    "N": {"AAT": 0.44, "AAC": 0.56},
    "K": {"AAA": 0.41, "AAG": 0.59},
    "D": {"GAT": 0.44, "GAC": 0.56},
    "E": {"GAA": 0.41, "GAG": 0.59},
    "C": {"TGT": 0.43, "TGC": 0.57},
    "R": {
        "CGT": 0.09, "CGC": 0.19, "CGA": 0.11,
        "CGG": 0.21, "AGA": 0.20, "AGG": 0.20,
    },
    "G": {"GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25},
    "W": {"TGG": 1.00},
    "M": {"ATG": 1.00},
    "*": {"TAA": 0.28, "TAG": 0.20, "TGA": 0.52},
}

# Reverse lookup: codon -> (amino_acid, frequency)
_CODON_TO_AA: Final[dict[str, tuple[str, float]]] = {}
_MAX_FREQ: Final[dict[str, float]] = {}

# Populated at module load time (below the dict definitions).


# ---------------------------------------------------------------------------
# Initialise reverse tables
# ---------------------------------------------------------------------------

def _init_codon_tables() -> None:
    """Build reverse-lookup tables from DOG_CODON_USAGE."""
    for aa, codons in DOG_CODON_USAGE.items():
        max_f = max(codons.values())
        _MAX_FREQ[aa] = max_f
        for codon, freq in codons.items():
            _CODON_TO_AA[codon] = (aa, freq)


_init_codon_tables()


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _load_fasta_sequence(path: Path) -> str:
    """Load a single-record FASTA file and return the sequence string.

    Strips whitespace and joins continuation lines.  Returns the raw
    sequence (RNA with U or DNA with T, depending on the file).
    """
    lines: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            lines.append(line)
    if not lines:
        raise ValueError(f"No sequence data found in {path}")
    return "".join(lines)


def _write_json(data: object, path: Path) -> None:
    """Write *data* as pretty-printed JSON to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
    logger.info("Wrote %s", path)


def _write_text(text: str, path: Path) -> None:
    """Write plain text to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    logger.info("Wrote %s", path)


# ---------------------------------------------------------------------------
# mRNA secondary structure analysis
# ---------------------------------------------------------------------------


def _compute_gc_content(sequence: str) -> float:
    """Return GC fraction for a nucleotide sequence."""
    if not sequence:
        return 0.0
    gc = sum(1 for nt in sequence.upper() if nt in "GC")
    return gc / len(sequence)


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA/RNA-like sequence (uppercase)."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G",
                  "U": "A"}
    return "".join(complement.get(nt, nt) for nt in reversed(seq.upper()))


def _count_hairpin_regions(sequence: str) -> int:
    """Count potential hairpin-forming palindromic regions.

    Scans with a sliding window and checks whether any sub-sequence of
    length >= HAIRPIN_MIN_LEN has its reverse complement within the same
    window.  This is a rough heuristic, not a thermodynamic prediction.
    """
    seq = sequence.upper()
    count = 0
    for start in range(0, len(seq) - HAIRPIN_WINDOW + 1, HAIRPIN_WINDOW // 2):
        window = seq[start : start + HAIRPIN_WINDOW]
        for k in range(HAIRPIN_MIN_LEN, len(window) // 2 + 1):
            left = window[:k]
            rc = _reverse_complement(left)
            if rc in window[k:]:
                count += 1
                break  # one hit per window is enough
    return count


def _gc_in_first_n(sequence: str, n: int = 50) -> float:
    """Return GC fraction in the first *n* nucleotides of a sequence."""
    segment = sequence[:n]
    return _compute_gc_content(segment)


def _analyze_structure(sequence: str) -> dict:
    """Analyse mRNA secondary structure.

    Attempts to use ViennaRNA (RNA.fold) for thermodynamic predictions.
    Falls back to heuristic metrics when ViennaRNA is not installed.

    Returns a dict of structure-related metrics.
    """
    result: dict = {
        "viennarna_available": False,
        "gc_content": round(_compute_gc_content(sequence), 4),
        "gc_content_pct": round(_compute_gc_content(sequence) * 100.0, 2),
    }

    try:
        import RNA  # ViennaRNA Python bindings  # noqa: F811

        result["viennarna_available"] = True

        # Fold 5'UTR + first 50 codons (150 nt of coding region)
        # The mRNA starts with the 5'UTR; approximate the UTR as first ~46 nt
        # (matching the known FIVE_PRIME_UTR length from module 11).
        utr_plus_start = sequence[:196]  # ~46 UTR + 150 coding
        utr_struct, utr_mfe = RNA.fold(utr_plus_start)
        result["five_prime_structure"] = utr_struct
        result["five_prime_mfe_kcal"] = round(utr_mfe, 2)
        result["five_prime_assessment"] = (
            "free" if utr_mfe > -15.0 else "structured"
        )

        # Fold full mRNA
        full_struct, full_mfe = RNA.fold(sequence)
        result["full_structure_dot_bracket"] = full_struct
        result["full_mfe_kcal"] = round(full_mfe, 2)

        logger.info(
            "ViennaRNA: 5'UTR MFE=%.2f kcal/mol (%s), full MFE=%.2f kcal/mol",
            utr_mfe,
            result["five_prime_assessment"],
            full_mfe,
        )

    except ImportError:
        logger.info(
            "ViennaRNA not available; using heuristic structure metrics."
        )
        result["hairpin_regions"] = _count_hairpin_regions(sequence)
        result["gc_first_50nt_coding"] = round(
            _gc_in_first_n(sequence[46:], 50), 4
        )
        result["gc_first_50nt_coding_pct"] = round(
            result["gc_first_50nt_coding"] * 100.0, 2
        )
        # Low GC in the first 50 nt of coding region favours ribosome entry.
        result["five_prime_assessment"] = (
            "favourable (low GC)"
            if result["gc_first_50nt_coding"] < 0.50
            else "potentially structured (high GC)"
        )

    return result


# ---------------------------------------------------------------------------
# Codon Adaptation Index
# ---------------------------------------------------------------------------


def _rna_to_dna(sequence: str) -> str:
    """Convert RNA (U) to DNA (T)."""
    return sequence.replace("U", "T").replace("u", "t")


def _extract_orf_codons(mrna_sequence: str) -> list[str]:
    """Extract codons from the ORF within the mRNA sequence.

    Finds the first AUG (or ATG) and reads triplets until a stop codon
    or end of sequence.
    """
    seq = _rna_to_dna(mrna_sequence).upper()

    start = seq.find("ATG")
    if start == -1:
        raise ValueError("No start codon (ATG/AUG) found in mRNA sequence.")

    codons: list[str] = []
    stop_codons = {"TAA", "TAG", "TGA"}
    for i in range(start, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) < 3:
            break
        codons.append(codon)
        if codon in stop_codons:
            break

    return codons


def _calculate_cai(codons: list[str]) -> float:
    """Calculate the Codon Adaptation Index for a list of codons.

    CAI = geometric mean of (w_i) for each codon, where
    w_i = freq(codon) / max_freq(amino_acid).

    Stop codons and unknown codons are excluded from the calculation.
    """
    log_sum = 0.0
    count = 0

    for codon in codons:
        if codon not in _CODON_TO_AA:
            continue
        aa, freq = _CODON_TO_AA[codon]
        if aa == "*":
            continue  # skip stop codons
        max_f = _MAX_FREQ[aa]
        if max_f == 0:
            continue
        w = freq / max_f
        if w > 0:
            log_sum += math.log(w)
            count += 1

    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


def _identify_worst_codons(
    codons: list[str],
    limit: int = WORST_CODONS_TO_REPORT,
) -> list[dict]:
    """Identify codons with the lowest relative adaptiveness.

    Returns up to *limit* entries sorted by ascending w-value, each
    containing position, codon, amino acid, current frequency, and
    the best synonym.
    """
    scored: list[dict] = []

    for idx, codon in enumerate(codons):
        if codon not in _CODON_TO_AA:
            continue
        aa, freq = _CODON_TO_AA[codon]
        if aa == "*":
            continue
        max_f = _MAX_FREQ[aa]
        if max_f == 0:
            continue
        w = freq / max_f

        # Already optimal
        if w >= 0.99:
            continue

        # Find best synonym
        best_codon = max(
            DOG_CODON_USAGE[aa], key=DOG_CODON_USAGE[aa].get,  # type: ignore[arg-type]
        )

        scored.append({
            "position": idx + 1,
            "codon": codon,
            "amino_acid": aa,
            "current_frequency": freq,
            "relative_adaptiveness": round(w, 4),
            "best_synonym": best_codon,
            "best_frequency": DOG_CODON_USAGE[aa][best_codon],
        })

    scored.sort(key=lambda x: x["relative_adaptiveness"])
    return scored[:limit]


def _suggest_reoptimization(
    codons: list[str],
    cai: float,
) -> dict:
    """Suggest codon re-optimization if CAI is below target.

    Returns a dict with the assessment and worst codons to fix.
    """
    needs_optimization = cai < CAI_TARGET
    worst = _identify_worst_codons(codons) if needs_optimization else []

    return {
        "needs_optimization": needs_optimization,
        "current_cai": round(cai, 4),
        "target_cai": CAI_TARGET,
        "worst_codons": worst,
        "codons_to_fix": len(worst),
    }


# ---------------------------------------------------------------------------
# LNP formulation
# ---------------------------------------------------------------------------


def _compute_lnp_formulation(mrna_length: int) -> dict:
    """Compute LNP formulation parameters based on mRNA length.

    Returns a dict with molecular weight, charge, composition, and
    dosing recommendations.
    """
    mw_da = mrna_length * AVG_NUCLEOTIDE_MW
    mw_kda = mw_da / 1000.0
    charge = mrna_length * PHOSPHATE_CHARGE

    return {
        "mrna_molecular_weight_da": round(mw_da, 2),
        "mrna_molecular_weight_kda": round(mw_kda, 2),
        "total_charge": charge,
        "np_ratio": NP_RATIO,
        "lipid_composition": [
            {
                "component": entry["component"],
                "mol_pct": entry["mol_pct"],
                "role": entry["role"],
            }
            for entry in LNP_COMPOSITION
        ],
        "target_particle_size_nm": TARGET_PARTICLE_SIZE,
        "estimated_encapsulation_efficiency": ENCAPSULATION_EFFICIENCY,
        "dosing": {
            "dose_per_dog_ug": DOSE_PER_DOG_UG,
            "dog_weight_kg": DOG_WEIGHT_KG,
            "schedule": [
                {"dose": 1, "day": 0, "label": "Day 0"},
                {"dose": 2, "day": 28, "label": "Day 28"},
                {"dose": 3, "day": 56, "label": "Day 56 (booster)"},
            ],
            "annual_booster": "Recommended for endemic areas",
        },
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def _build_delivery_report(
    mrna_length: int,
    structure: dict,
    cai: float,
    optimization: dict,
    lnp: dict,
) -> str:
    """Generate the Markdown delivery report."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines: list[str] = []
    lines.append("# mRNA Delivery Optimization Report")
    lines.append("")
    lines.append(f"Generated: {now}")
    lines.append("")

    # -- mRNA Characteristics --
    lines.append("## mRNA Characteristics")
    lines.append(f"- Length: {mrna_length} nt")
    lines.append(f"- GC content: {structure['gc_content_pct']}%")
    lines.append(f"- CAI (Canis lupus familiaris): {cai:.2f}")

    if structure.get("viennarna_available"):
        lines.append(
            f"- Minimum Free Energy: {structure['full_mfe_kcal']} kcal/mol"
        )
        lines.append(
            f"- 5'UTR + start region MFE: "
            f"{structure['five_prime_mfe_kcal']} kcal/mol"
        )
    else:
        lines.append("- Minimum Free Energy: N/A (ViennaRNA not installed)")
        lines.append(
            f"- Potential hairpin regions: {structure.get('hairpin_regions', 'N/A')}"
        )
        lines.append(
            f"- GC in first 50 nt of coding region: "
            f"{structure.get('gc_first_50nt_coding_pct', 'N/A')}%"
        )

    lines.append(f"- 5'UTR structure: {structure['five_prime_assessment']}")
    lines.append("")

    # -- Codon Optimization --
    lines.append("## Codon Optimization")
    lines.append(f"- Current CAI: {cai:.2f}")
    lines.append(f"- Target CAI: >= {CAI_TARGET:.2f}")

    if optimization["needs_optimization"]:
        worst = optimization["worst_codons"]
        lines.append(f"- Codons to optimize: {len(worst)}")
        lines.append("")
        lines.append(
            "| # | Position | Codon | AA | w | Best synonym | Best freq |"
        )
        lines.append(
            "|---|----------|-------|----|---|-------------|-----------|"
        )
        for i, entry in enumerate(worst, 1):
            lines.append(
                f"| {i} | {entry['position']} | {entry['codon']} "
                f"| {entry['amino_acid']} | {entry['relative_adaptiveness']:.4f} "
                f"| {entry['best_synonym']} | {entry['best_frequency']:.2f} |"
            )
    else:
        lines.append(f"- Status: CAI is already >= {CAI_TARGET:.2f}; no changes needed")

    lines.append("")

    # -- LNP Formulation --
    lines.append("## Recommended LNP Formulation")
    lines.append("")
    lines.append("| Component | Mol% | Role |")
    lines.append("|-----------|:----:|------|")
    for entry in lnp["lipid_composition"]:
        lines.append(
            f"| {entry['component']} | {entry['mol_pct']} | {entry['role']} |"
        )
    lines.append("")
    lines.append(f"- N/P ratio: {lnp['np_ratio']}")
    lines.append(f"- Target particle size: {lnp['target_particle_size_nm']}")
    lines.append(f"- mRNA MW: {lnp['mrna_molecular_weight_kda']:.1f} kDa")
    lines.append(
        f"- Estimated encapsulation efficiency: "
        f"{lnp['estimated_encapsulation_efficiency']}"
    )
    lines.append(
        f"- Estimated dose per dog ({DOG_WEIGHT_KG} kg): "
        f"{DOSE_PER_DOG_UG} ug mRNA"
    )
    lines.append("")

    # -- Dosing Schedule --
    lines.append("## Dosing Schedule")
    for step in lnp["dosing"]["schedule"]:
        lines.append(f"- Dose {step['dose']}: {step['label']}")
    lines.append(f"- {lnp['dosing']['annual_booster']}")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def optimize_delivery(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Run the mRNA delivery optimization analysis.

    1. Load mRNA sequence from the construct stage.
    2. Analyse secondary structure (ViennaRNA or fallback heuristics).
    3. Calculate Codon Adaptation Index for Canis lupus familiaris.
    4. Suggest codon re-optimization if CAI < 0.80.
    5. Compute LNP formulation recommendation.
    6. Generate delivery report and JSON analysis.

    Args:
        force: If True, overwrite existing outputs even if present.
        dry_run: If True, compute all metrics (local computation only)
            and generate the full report.  ViennaRNA is used if installed,
            otherwise the fallback heuristics are used.

    Returns:
        Path to the delivery report Markdown file.
    """
    output_base = Path(OUTPUT_DIR)
    report_path = output_base / "delivery_report.md"
    json_path = output_base / "delivery_analysis.json"

    if report_path.exists() and not force and not dry_run:
        logger.info(
            "Output already exists at %s. Use --force to overwrite.",
            report_path,
        )
        return str(report_path)

    # ------------------------------------------------------------------
    # 1. Load mRNA sequence
    # ------------------------------------------------------------------
    mrna_path = Path(MRNA_FASTA_PATH)
    if not mrna_path.exists():
        logger.error("mRNA FASTA not found at %s.", mrna_path)
        sys.exit(1)

    mrna_sequence = _load_fasta_sequence(mrna_path)
    mrna_length = len(mrna_sequence)
    logger.info("Loaded mRNA sequence: %d nt from %s", mrna_length, mrna_path)

    # ------------------------------------------------------------------
    # 2. mRNA secondary structure analysis
    # ------------------------------------------------------------------
    logger.info("Analysing mRNA secondary structure...")
    structure = _analyze_structure(mrna_sequence)
    logger.info(
        "GC content: %.1f%%, 5'UTR assessment: %s",
        structure["gc_content_pct"],
        structure["five_prime_assessment"],
    )

    # ------------------------------------------------------------------
    # 3. Codon Adaptation Index
    # ------------------------------------------------------------------
    logger.info("Calculating Codon Adaptation Index (CAI)...")
    codons = _extract_orf_codons(mrna_sequence)
    cai = _calculate_cai(codons)
    logger.info("CAI (Canis lupus familiaris): %.4f (%d codons)", cai, len(codons))

    # ------------------------------------------------------------------
    # 4. Suggest codon re-optimization if CAI < 0.80
    # ------------------------------------------------------------------
    optimization = _suggest_reoptimization(codons, cai)
    if optimization["needs_optimization"]:
        logger.info(
            "CAI %.4f < %.2f target; %d codons flagged for re-optimization.",
            cai,
            CAI_TARGET,
            optimization["codons_to_fix"],
        )
    else:
        logger.info("CAI %.4f meets target (>= %.2f).", cai, CAI_TARGET)

    # ------------------------------------------------------------------
    # 5. LNP formulation recommendation
    # ------------------------------------------------------------------
    logger.info("Computing LNP formulation parameters...")
    lnp = _compute_lnp_formulation(mrna_length)
    logger.info(
        "mRNA MW: %.1f kDa, charge: %d, N/P ratio: %s",
        lnp["mrna_molecular_weight_kda"],
        lnp["total_charge"],
        lnp["np_ratio"],
    )

    # ------------------------------------------------------------------
    # 6. Generate delivery report
    # ------------------------------------------------------------------
    logger.info("Generating delivery report...")

    analysis_data = {
        "mrna_source": str(mrna_path),
        "mrna_length_nt": mrna_length,
        "structure": structure,
        "cai": round(cai, 4),
        "codon_optimization": optimization,
        "lnp_formulation": lnp,
        "generated_at": datetime.now().isoformat(),
    }
    _write_json(analysis_data, json_path)

    report = _build_delivery_report(
        mrna_length=mrna_length,
        structure=structure,
        cai=cai,
        optimization=optimization,
        lnp=lnp,
    )
    _write_text(report, report_path)

    logger.info("Delivery analysis complete. Report: %s", report_path)
    return str(report_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Marley -- mRNA delivery optimization: secondary structure "
            "analysis, Codon Adaptation Index, LNP formulation, and "
            "dosing schedule."
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
            "Compute all metrics (local computation only) and generate "
            "the full report. ViennaRNA is used if installed."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    result_path = optimize_delivery(
        force=args.force,
        dry_run=args.dry_run,
    )
    logger.info("Done. Output: %s", result_path)
