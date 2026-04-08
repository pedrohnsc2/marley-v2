"""Design antisense oligonucleotides (ASOs) targeting the Leishmania SL RNA.

The Spliced Leader RNA is a conserved 39 nt sequence added to the 5' end of
every trypanosomatid mRNA via trans-splicing.  Because it is absent in
mammalian hosts, it represents a selective therapeutic target.

This module:

1. Predicts the secondary structure of the SL RNA (ViennaRNA or fallback).
2. Enumerates all possible ASO candidates of length 18--25 nt.
3. Computes thermodynamic properties (Tm, delta_G) using the SantaLucia 1998
   nearest-neighbor model.
4. Scores candidates on accessibility, Tm, binding energy, and
   self-complementarity.
5. Outputs ranked candidates as CSV and JSON.

Usage:
    python -m rna_entropy.08_aso_design
    python -m rna_entropy.08_aso_design --force
    python -m rna_entropy.08_aso_design --dry-run
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# The 39 nt Spliced Leader RNA sequence of Leishmania infantum.
SL_SEQUENCE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"

ASO_MIN_LENGTH: Final[int] = 18
ASO_MAX_LENGTH: Final[int] = 25

# GC content bounds — relaxed lower bound because SL has GC = 28.2%
GC_MIN: Final[float] = 0.25
GC_MAX: Final[float] = 0.65

# Melting temperature bounds (Celsius)
TM_MIN: Final[float] = 40.0
TM_MAX: Final[float] = 70.0

# Maximum consecutive self-complementary base pairs allowed
SELF_COMP_MAX: Final[int] = 4

# Standard ASO concentration for Tm calculation (250 nM)
ASO_CONCENTRATION: Final[float] = 250e-9  # mol/L

# Gas constant in cal / (mol * K)
R_GAS: Final[float] = 1.987

# SantaLucia 1998 unified nearest-neighbor parameters for DNA/DNA duplexes.
# Each key is a dinucleotide pair XY/X'Y' where X'Y' is the complement.
# Values: (delta_H in kcal/mol, delta_S in cal/(mol*K))
NN_PARAMS: Final[dict[str, tuple[float, float]]] = {
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}

# Initiation parameters: added once per duplex (delta_H kcal/mol, delta_S cal/(mol*K))
NN_INIT: Final[tuple[float, float]] = (0.1, -2.8)

# DNA complement mapping (includes U for RNA input)
COMPLEMENT: Final[dict[str, str]] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "U": "A",
}

# Output paths
OUTPUT_DIR: Final[str] = "results/aso"
STRUCTURE_JSON: Final[str] = "results/aso/sl_secondary_structure.json"
CANDIDATES_CSV: Final[str] = "results/aso/aso_candidates.csv"
CANDIDATES_JSON: Final[str] = "results/aso/aso_candidates.json"

# CSV column order
CSV_COLUMNS: Final[list[str]] = [
    "rank",
    "aso_id",
    "aso_sequence",
    "length",
    "target_start",
    "target_end",
    "gc_content",
    "tm_celsius",
    "delta_g_kcal",
    "accessibility_score",
    "self_comp_pairs",
    "composite_score",
]

# Composite score weights
WEIGHT_ACCESSIBILITY: Final[float] = 0.35
WEIGHT_TM: Final[float] = 0.25
WEIGHT_DG: Final[float] = 0.25
WEIGHT_SELF_COMP: Final[float] = 0.15

logger = get_logger("aso_design")


# ---------------------------------------------------------------------------
# SL RNA secondary structure prediction
# ---------------------------------------------------------------------------


def predict_sl_structure(sequence: str) -> dict[str, Any]:
    """Predict secondary structure of the SL RNA.

    Attempts to use ViennaRNA (RNA.fold) for accurate minimum free energy
    (MFE) structure prediction.  Falls back to a simple GC-based estimation
    when ViennaRNA is not installed.

    Args:
        sequence: RNA or DNA sequence to fold.

    Returns:
        Dictionary with keys:
            - sequence: the input sequence
            - dot_bracket: MFE structure in dot-bracket notation
            - mfe: minimum free energy in kcal/mol
            - accessible_regions: list of (start, end) tuples for unpaired runs
    """
    gc_count = sum(1 for b in sequence.upper() if b in ("G", "C"))
    gc_content = gc_count / len(sequence) if sequence else 0.0

    try:
        import RNA  # type: ignore[import-untyped]

        structure, mfe = RNA.fold(sequence)
        logger.info("ViennaRNA fold: MFE = %.2f kcal/mol", mfe)
    except ImportError:
        # Fallback: for a low-GC RNA (~28%), most positions are unpaired.
        # Estimate all positions as accessible (dot-bracket = all dots).
        # MFE approximation: -0.3 * length * GC_fraction
        structure = "." * len(sequence)
        mfe = -0.3 * len(sequence) * gc_content
        logger.warning(
            "ViennaRNA not available. Using fallback structure estimation "
            "(MFE ~ %.2f kcal/mol).",
            mfe,
        )

    # Extract accessible (unpaired) regions: contiguous runs of '.' in the
    # dot-bracket notation.
    accessible_regions: list[tuple[int, int]] = []
    i = 0
    while i < len(structure):
        if structure[i] == ".":
            start = i
            while i < len(structure) and structure[i] == ".":
                i += 1
            # Region spans [start, i) — store as inclusive (start, end)
            accessible_regions.append((start, i - 1))
        else:
            i += 1

    return {
        "sequence": sequence,
        "dot_bracket": structure,
        "mfe": round(mfe, 2),
        "accessible_regions": accessible_regions,
    }


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA/RNA sequence.

    Uses the COMPLEMENT mapping: A<->T, G<->C, U->A.
    Unknown bases are kept as-is.

    Args:
        seq: Input nucleotide sequence (case-insensitive).

    Returns:
        Reverse complement string in uppercase.
    """
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq.upper()))


def _compute_gc_content(sequence: str) -> float:
    """Calculate GC content as a fraction in [0.0, 1.0].

    Args:
        sequence: Nucleotide sequence (case-insensitive).

    Returns:
        GC fraction.  Returns 0.0 for an empty sequence.
    """
    if not sequence:
        return 0.0
    upper = sequence.upper()
    gc_count = upper.count("G") + upper.count("C")
    return gc_count / len(upper)


# ---------------------------------------------------------------------------
# Thermodynamic calculations (SantaLucia 1998 nearest-neighbor model)
# ---------------------------------------------------------------------------


def _get_nn_params(dinucleotide: str) -> tuple[float, float]:
    """Look up nearest-neighbor parameters for a dinucleotide pair.

    The dinucleotide is converted to its canonical NN_PARAMS key by trying
    both the forward representation (XY/X'Y') and the reverse complement.

    Args:
        dinucleotide: Two-character DNA string (e.g., "AT").

    Returns:
        Tuple of (delta_H in kcal/mol, delta_S in cal/(mol*K)).
        Returns (0.0, 0.0) if no match is found.
    """
    base1, base2 = dinucleotide[0], dinucleotide[1]
    comp1 = COMPLEMENT.get(base1, base1)
    comp2 = COMPLEMENT.get(base2, base2)

    # Forward key: XY/X'Y'  (e.g., "AT" -> "AT/TA")
    forward_key = f"{base1}{base2}/{comp1}{comp2}"
    if forward_key in NN_PARAMS:
        return NN_PARAMS[forward_key]

    # Reverse complement key: Y'X'/YX  (e.g., "AT" -> "TA/AT")
    reverse_key = f"{comp2}{comp1}/{base2}{base1}"
    if reverse_key in NN_PARAMS:
        return NN_PARAMS[reverse_key]

    # No match — log once and return zero contribution
    logger.debug("No NN parameters for dinucleotide %s", dinucleotide)
    return (0.0, 0.0)


def compute_tm(sequence: str) -> float:
    """Calculate melting temperature using the SantaLucia 1998 nearest-neighbor model.

    For sequences >= 14 nt (nearest-neighbor method):
        Tm = (delta_H * 1000) / (delta_S + R * ln(Ct / 4)) - 273.15

    Where:
        delta_H = sum of nearest-neighbor enthalpy contributions + initiation
        delta_S = sum of nearest-neighbor entropy contributions + initiation
        R = 1.987 cal/(mol*K)  (gas constant)
        Ct = total strand concentration (250 nM for ASOs)

    For sequences < 14 nt (Wallace rule):
        Tm = 2 * (count_A + count_T) + 4 * (count_G + count_C)

    Args:
        sequence: DNA oligonucleotide sequence (case-insensitive).

    Returns:
        Predicted melting temperature in degrees Celsius.
    """
    seq = sequence.upper()

    if len(seq) < 14:
        # Wallace rule: Tm = 2*(A+T) + 4*(G+C)
        count_at = seq.count("A") + seq.count("T")
        count_gc = seq.count("G") + seq.count("C")
        return float(2 * count_at + 4 * count_gc)

    # Sum nearest-neighbor parameters over all consecutive dinucleotides
    total_dh = NN_INIT[0]  # kcal/mol — initiation enthalpy
    total_ds = NN_INIT[1]  # cal/(mol*K) — initiation entropy

    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        dh, ds = _get_nn_params(dinuc)
        total_dh += dh  # accumulate enthalpy (kcal/mol)
        total_ds += ds  # accumulate entropy (cal/(mol*K))

    # Tm = (delta_H * 1000) / (delta_S + R * ln(Ct/4)) - 273.15
    # delta_H is in kcal/mol, multiply by 1000 to convert to cal/mol
    # delta_S is already in cal/(mol*K)
    # Ct = ASO_CONCENTRATION (250 nM = 250e-9 M)
    ln_ct_over_4 = math.log(ASO_CONCENTRATION / 4.0)
    denominator = total_ds + R_GAS * ln_ct_over_4

    if abs(denominator) < 1e-12:
        # Avoid division by zero — return midpoint of valid range
        return (TM_MIN + TM_MAX) / 2.0

    # Convert delta_H from kcal to cal (multiply by 1000), then divide
    tm_kelvin = (total_dh * 1000.0) / denominator
    tm_celsius = tm_kelvin - 273.15

    return round(tm_celsius, 2)


def compute_binding_energy(aso_seq: str, target_seq: str) -> float:
    """Calculate delta_G of ASO-target hybridization at 37 C (310.15 K).

    Uses the same SantaLucia 1998 nearest-neighbor parameters.

    delta_G_37 = delta_H - (310.15 * delta_S / 1000)

    Where:
        delta_H = sum of NN enthalpies + initiation  (kcal/mol)
        delta_S = sum of NN entropies + initiation   (cal/(mol*K))
        310.15 K = 37 C in Kelvin
        Division by 1000 converts delta_S from cal to kcal

    A more negative delta_G indicates more favorable (stronger) binding.

    Args:
        aso_seq: Antisense oligonucleotide sequence (5'->3').
        target_seq: Target RNA/DNA sequence (5'->3').

    Returns:
        delta_G in kcal/mol (negative = favorable binding).
    """
    seq = aso_seq.upper()

    # Sum nearest-neighbor parameters along the ASO sequence
    total_dh = NN_INIT[0]  # kcal/mol — initiation enthalpy
    total_ds = NN_INIT[1]  # cal/(mol*K) — initiation entropy

    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        dh, ds = _get_nn_params(dinuc)
        total_dh += dh  # accumulate enthalpy (kcal/mol)
        total_ds += ds  # accumulate entropy (cal/(mol*K))

    # delta_G = delta_H - T * delta_S
    # T = 310.15 K (37 C)
    # delta_S is in cal/(mol*K), divide by 1000 to get kcal/(mol*K)
    temperature_k = 310.15  # 37 C in Kelvin
    delta_g = total_dh - temperature_k * (total_ds / 1000.0)

    return round(delta_g, 2)


# ---------------------------------------------------------------------------
# Self-complementarity detection
# ---------------------------------------------------------------------------


def check_self_complementarity(
    sequence: str,
) -> tuple[int, list[tuple[int, int]]]:
    """Detect potential hairpin formation in an ASO sequence.

    For each pair of positions (i, j) where j > i + 3 (minimum loop size = 3),
    checks if the bases at i and j are complementary.  Tracks the longest run
    of consecutive complementary pairs that could form a hairpin stem.

    Args:
        sequence: ASO DNA sequence (case-insensitive).

    Returns:
        Tuple of (max_consecutive_pairs, hairpin_positions) where
        hairpin_positions is a list of (stem_start, stem_end) tuples for
        each detected hairpin region with >= 3 consecutive pairs.
    """
    seq = sequence.upper()
    n = len(seq)
    max_consecutive = 0
    hairpin_positions: list[tuple[int, int]] = []

    # Slide a potential stem start along the sequence
    for i in range(n):
        # The partner must be at least 4 positions away (3 nt loop minimum)
        for j in range(i + 4, n):
            # Count consecutive complementary pairs starting at (i, j)
            # moving inward: i+k pairs with j-k
            consecutive = 0
            k = 0
            while (i + k) < (j - k) and (j - k) >= 0:
                base_left = seq[i + k]
                base_right = seq[j - k]
                # Check if bases are complementary
                if COMPLEMENT.get(base_left) == base_right:
                    consecutive += 1
                    k += 1
                else:
                    break

            if consecutive > max_consecutive:
                max_consecutive = consecutive

            # Record hairpin positions if stem is >= 3 pairs
            if consecutive >= 3:
                hairpin_positions.append((i, i + consecutive - 1))

    return max_consecutive, hairpin_positions


# ---------------------------------------------------------------------------
# Candidate generation and scoring
# ---------------------------------------------------------------------------


def _compute_accessibility(
    dot_bracket: str, start: int, end: int
) -> float:
    """Compute the fraction of unpaired positions in a target region.

    Accessibility is defined as the fraction of '.' characters in the
    dot-bracket notation over the region [start, end).

    Args:
        dot_bracket: Full dot-bracket structure string.
        start: Start index (inclusive).
        end: End index (exclusive).

    Returns:
        Fraction of unpaired (accessible) positions in [0.0, 1.0].
    """
    region = dot_bracket[start:end]
    if not region:
        return 0.0
    # Count unpaired positions (dots) divided by region length
    unpaired = region.count(".")
    return unpaired / len(region)


def _normalize_tm(tm: float) -> float:
    """Normalize Tm to [0, 1] based on the acceptable range.

    norm_tm = (Tm - TM_MIN) / (TM_MAX - TM_MIN), clamped to [0, 1]

    Args:
        tm: Melting temperature in Celsius.

    Returns:
        Normalized Tm score in [0.0, 1.0].
    """
    if TM_MAX == TM_MIN:
        return 0.5
    normalized = (tm - TM_MIN) / (TM_MAX - TM_MIN)
    return max(0.0, min(1.0, normalized))


def _normalize_dg(delta_g: float) -> float:
    """Normalize delta_G to [0, 1].

    norm_dg = min(1.0, abs(delta_G) / 50.0)

    A delta_G of -50 kcal/mol or more negative maps to 1.0.

    Args:
        delta_g: Binding free energy in kcal/mol.

    Returns:
        Normalized binding energy score in [0.0, 1.0].
    """
    return min(1.0, abs(delta_g) / 50.0)


def _compute_composite_score(
    accessibility: float,
    norm_tm: float,
    norm_dg: float,
    self_comp: int,
    max_self_comp: int,
) -> float:
    """Compute the composite ASO score.

    composite = 0.35 * accessibility
              + 0.25 * norm_tm
              + 0.25 * norm_dg
              + 0.15 * (1 - self_comp / max_self_comp)

    Where:
        accessibility = fraction of unpaired target positions
        norm_tm = normalized melting temperature [0, 1]
        norm_dg = normalized binding energy [0, 1]
        self_comp / max_self_comp = penalty for self-complementarity

    Args:
        accessibility: Target region accessibility [0, 1].
        norm_tm: Normalized Tm [0, 1].
        norm_dg: Normalized delta_G [0, 1].
        self_comp: Number of max consecutive self-complementary pairs.
        max_self_comp: Maximum self-comp observed across all candidates
                       (used for normalization; minimum 1 to avoid division by zero).

    Returns:
        Composite score in [0.0, 1.0].
    """
    # Self-complementarity penalty: lower self-comp is better, so we use
    # (1 - self_comp / max_self_comp).  Clamped to [0, 1].
    safe_max = max(max_self_comp, 1)
    self_comp_score = 1.0 - min(self_comp / safe_max, 1.0)

    score = (
        WEIGHT_ACCESSIBILITY * accessibility
        + WEIGHT_TM * norm_tm
        + WEIGHT_DG * norm_dg
        + WEIGHT_SELF_COMP * self_comp_score
    )

    return round(score, 4)


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------


def _write_json(data: Any, output_path: Path) -> None:
    """Write data as pretty-printed JSON.

    Args:
        data: Serializable data (dict, list, etc.).
        output_path: Destination file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    logger.info("Wrote %s", output_path)


def _write_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Write candidate rows to a CSV file.

    Args:
        rows: List of dicts with keys matching CSV_COLUMNS.
        output_path: Destination CSV file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d row(s) to %s", len(rows), output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_aso_candidates(
    force: bool = False, dry_run: bool = False
) -> str:
    """Generate all ASO candidates targeting the Leishmania SL RNA.

    Algorithm:
    1. Predict SL secondary structure (ViennaRNA or fallback).
    2. For each length L in [ASO_MIN_LENGTH, ASO_MAX_LENGTH]:
       For each position i in [0, len(SL) - L]:
         - Extract target window: SL[i : i + L]
         - Generate ASO: reverse_complement(target)
         - Compute GC content, Tm, delta_G, self-complementarity
         - Compute accessibility from secondary structure
    3. Filter candidates: GC >= GC_MIN, Tm >= TM_MIN, self_comp <= SELF_COMP_MAX
    4. Score: composite = 0.35*accessibility + 0.25*norm_tm
                         + 0.25*norm_dg + 0.15*(1 - self_comp/max)
    5. Rank and assign IDs: MRL-ASO-001, MRL-ASO-002, ...
    6. Save to results/aso/ as CSV and JSON.

    This module performs no API calls — all calculations are local.
    In dry-run mode, extra diagnostic information is logged.

    Args:
        force: Re-run even if output files already exist.
        dry_run: Log extra diagnostic information (all computation is
                 local regardless).

    Returns:
        Path to the output JSON file.
    """
    json_path = Path(CANDIDATES_JSON)
    csv_path = Path(CANDIDATES_CSV)
    structure_path = Path(STRUCTURE_JSON)

    if json_path.exists() and not force:
        logger.info(
            "ASO candidates already exist at %s. Use --force to re-run.",
            json_path,
        )
        return str(json_path)

    if dry_run:
        logger.info(
            "[DRY RUN] All ASO computations are local — "
            "dry-run logs extra diagnostics."
        )

    # ------------------------------------------------------------------
    # Step 1: Predict SL RNA secondary structure
    # ------------------------------------------------------------------
    logger.info(
        "Predicting SL RNA secondary structure (length=%d, GC=%.1f%%)...",
        len(SL_SEQUENCE),
        _compute_gc_content(SL_SEQUENCE) * 100,
    )
    sl_structure = predict_sl_structure(SL_SEQUENCE)
    _write_json(sl_structure, structure_path)
    logger.info(
        "SL structure: MFE=%.2f kcal/mol, %d accessible region(s)",
        sl_structure["mfe"],
        len(sl_structure["accessible_regions"]),
    )

    dot_bracket: str = sl_structure["dot_bracket"]

    if dry_run:
        logger.info("[DRY RUN] Dot-bracket: %s", dot_bracket)

    # ------------------------------------------------------------------
    # Step 2: Enumerate and evaluate all ASO candidates
    # ------------------------------------------------------------------
    sl_len = len(SL_SEQUENCE)
    raw_candidates: list[dict[str, Any]] = []

    for length in range(ASO_MIN_LENGTH, ASO_MAX_LENGTH + 1):
        for start in range(0, sl_len - length + 1):
            end = start + length

            # Extract target region from SL RNA
            target = SL_SEQUENCE[start:end]

            # ASO is the reverse complement of the target
            aso_seq = reverse_complement(target)

            # GC content of the ASO
            gc = _compute_gc_content(aso_seq)

            # Melting temperature via nearest-neighbor model
            tm = compute_tm(aso_seq)

            # Binding free energy at 37 C
            delta_g = compute_binding_energy(aso_seq, target)

            # Self-complementarity check for hairpin potential
            self_comp, hairpin_pos = check_self_complementarity(aso_seq)

            # Accessibility: fraction of unpaired bases in the target region
            accessibility = _compute_accessibility(dot_bracket, start, end)

            raw_candidates.append({
                "aso_sequence": aso_seq,
                "length": length,
                "target_start": start,
                "target_end": end,
                "target_sequence": target,
                "gc_content": round(gc, 4),
                "tm_celsius": tm,
                "delta_g_kcal": delta_g,
                "accessibility_score": round(accessibility, 4),
                "self_comp_pairs": self_comp,
                "hairpin_positions": hairpin_pos,
            })

    logger.info(
        "Enumerated %d raw ASO candidates (lengths %d--%d).",
        len(raw_candidates),
        ASO_MIN_LENGTH,
        ASO_MAX_LENGTH,
    )

    if dry_run:
        logger.info(
            "[DRY RUN] GC range: %.2f -- %.2f",
            min(c["gc_content"] for c in raw_candidates),
            max(c["gc_content"] for c in raw_candidates),
        )
        logger.info(
            "[DRY RUN] Tm range: %.1f -- %.1f C",
            min(c["tm_celsius"] for c in raw_candidates),
            max(c["tm_celsius"] for c in raw_candidates),
        )

    # ------------------------------------------------------------------
    # Step 3: Filter candidates
    # ------------------------------------------------------------------
    filtered: list[dict[str, Any]] = []
    for cand in raw_candidates:
        # GC content must be at least GC_MIN (no upper filter — SL is low-GC)
        if cand["gc_content"] < GC_MIN:
            continue
        # Tm must be at least TM_MIN for stable binding
        if cand["tm_celsius"] < TM_MIN:
            continue
        # Self-complementarity must not exceed SELF_COMP_MAX
        if cand["self_comp_pairs"] > SELF_COMP_MAX:
            continue
        filtered.append(cand)

    logger.info(
        "After filtering: %d of %d candidates passed "
        "(GC >= %.0f%%, Tm >= %.0f C, self-comp <= %d).",
        len(filtered),
        len(raw_candidates),
        GC_MIN * 100,
        TM_MIN,
        SELF_COMP_MAX,
    )

    if not filtered:
        logger.warning(
            "No candidates passed filters. Relaxing: keeping all %d candidates.",
            len(raw_candidates),
        )
        filtered = raw_candidates

    # ------------------------------------------------------------------
    # Step 4: Score candidates
    # ------------------------------------------------------------------
    # Find the maximum self-complementarity for normalization
    max_self_comp = max(
        (c["self_comp_pairs"] for c in filtered), default=1
    )

    for cand in filtered:
        n_tm = _normalize_tm(cand["tm_celsius"])
        n_dg = _normalize_dg(cand["delta_g_kcal"])
        cand["composite_score"] = _compute_composite_score(
            accessibility=cand["accessibility_score"],
            norm_tm=n_tm,
            norm_dg=n_dg,
            self_comp=cand["self_comp_pairs"],
            max_self_comp=max_self_comp,
        )

    # Sort by composite score descending
    filtered.sort(key=lambda c: c["composite_score"], reverse=True)

    # ------------------------------------------------------------------
    # Step 5: Assign IDs and ranks
    # ------------------------------------------------------------------
    for rank, cand in enumerate(filtered, start=1):
        cand["rank"] = rank
        cand["aso_id"] = f"MRL-ASO-{rank:03d}"

    logger.info(
        "Top candidate: %s (score=%.4f, Tm=%.1f C, dG=%.1f kcal/mol)",
        filtered[0]["aso_id"] if filtered else "none",
        filtered[0]["composite_score"] if filtered else 0.0,
        filtered[0]["tm_celsius"] if filtered else 0.0,
        filtered[0]["delta_g_kcal"] if filtered else 0.0,
    )

    # ------------------------------------------------------------------
    # Step 6: Save outputs
    # ------------------------------------------------------------------
    # CSV rows (subset of fields, no hairpin_positions or target_sequence)
    csv_rows: list[dict[str, Any]] = []
    for cand in filtered:
        csv_rows.append({
            "rank": cand["rank"],
            "aso_id": cand["aso_id"],
            "aso_sequence": cand["aso_sequence"],
            "length": cand["length"],
            "target_start": cand["target_start"],
            "target_end": cand["target_end"],
            "gc_content": cand["gc_content"],
            "tm_celsius": cand["tm_celsius"],
            "delta_g_kcal": cand["delta_g_kcal"],
            "accessibility_score": cand["accessibility_score"],
            "self_comp_pairs": cand["self_comp_pairs"],
            "composite_score": cand["composite_score"],
        })
    _write_csv(csv_rows, csv_path)

    # JSON output (full data including metadata)
    json_output: dict[str, Any] = {
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "sl_sequence": SL_SEQUENCE,
        "sl_length": len(SL_SEQUENCE),
        "sl_gc_content": round(_compute_gc_content(SL_SEQUENCE), 4),
        "sl_mfe": sl_structure["mfe"],
        "aso_length_range": [ASO_MIN_LENGTH, ASO_MAX_LENGTH],
        "filters": {
            "gc_min": GC_MIN,
            "gc_max": GC_MAX,
            "tm_min": TM_MIN,
            "tm_max": TM_MAX,
            "self_comp_max": SELF_COMP_MAX,
        },
        "scoring_weights": {
            "accessibility": WEIGHT_ACCESSIBILITY,
            "tm": WEIGHT_TM,
            "delta_g": WEIGHT_DG,
            "self_complementarity": WEIGHT_SELF_COMP,
        },
        "total_enumerated": len(raw_candidates),
        "total_after_filter": len(filtered),
        "candidates": [
            {
                "rank": c["rank"],
                "aso_id": c["aso_id"],
                "aso_sequence": c["aso_sequence"],
                "length": c["length"],
                "target_start": c["target_start"],
                "target_end": c["target_end"],
                "target_sequence": c["target_sequence"],
                "gc_content": c["gc_content"],
                "tm_celsius": c["tm_celsius"],
                "delta_g_kcal": c["delta_g_kcal"],
                "accessibility_score": c["accessibility_score"],
                "self_comp_pairs": c["self_comp_pairs"],
                "composite_score": c["composite_score"],
            }
            for c in filtered
        ],
    }
    _write_json(json_output, json_path)

    logger.info(
        "ASO design complete: %d candidates saved to %s and %s.",
        len(filtered),
        csv_path,
        json_path,
    )

    return str(json_path)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Design antisense oligonucleotides (ASOs) targeting the "
            "Leishmania Spliced Leader RNA."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run design even if results already exist.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run with extra diagnostic logging (all computation is local).",
    )
    args = parser.parse_args()

    result = generate_aso_candidates(force=args.force, dry_run=args.dry_run)
    logger.info("Complete: %s", result)
