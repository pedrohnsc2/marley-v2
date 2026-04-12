"""Modulo 08 -- Predicao de estrutura e acessibilidade do SL RNA."""

from __future__ import annotations

from typing import Any

from aso_math.config import (
    ASO_SEQUENCE, ASO_TARGET_END, ASO_TARGET_START, SL_SEQUENCE,
)
from aso_math.envelope import Timer, create_envelope, write_result

from .accessibility import analyze_accessibility
from .structure import predict_structure


MODULE_NAME = "08_sl_structure"


def main() -> dict[str, Any]:
    """Executa analise de estrutura e acessibilidade do SL RNA."""
    envelope = create_envelope(MODULE_NAME)
    with Timer() as timer:
        free_result = predict_structure(SL_SEQUENCE)
        aso_blocked = set(range(ASO_TARGET_START, ASO_TARGET_END))
        bound_result = predict_structure(SL_SEQUENCE, blocked_positions=aso_blocked)
        accessibility = analyze_accessibility(
            free_result.dot_bracket, ASO_TARGET_START, ASO_TARGET_END,
        )
        free_pairs_set = set(free_result.pairs)
        bound_pairs_set = set(bound_result.pairs)
        disrupted_pairs = free_pairs_set - bound_pairs_set
        new_pairs = bound_pairs_set - free_pairs_set
        maintained_pairs = free_pairs_set & bound_pairs_set
        disrupted_in_aso = [
            (i, j) for i, j in disrupted_pairs
            if i in aso_blocked or j in aso_blocked
        ]
        if accessibility.classification == "accessible":
            conclusion = (
                f"A regiao-alvo do ASO (pos {ASO_TARGET_START}-{ASO_TARGET_END - 1}) "
                f"e predominantemente acessivel (score {accessibility.mean_accessibility:.2f})."
            )
        elif accessibility.classification == "partially_blocked":
            conclusion = (
                f"A regiao-alvo (pos {ASO_TARGET_START}-{ASO_TARGET_END - 1}) "
                f"esta parcialmente bloqueada (score {accessibility.mean_accessibility:.2f})."
            )
        else:
            conclusion = (
                f"A regiao-alvo (pos {ASO_TARGET_START}-{ASO_TARGET_END - 1}) "
                f"esta bloqueada (score {accessibility.mean_accessibility:.2f})."
            )
    envelope["runtime_seconds"] = timer.elapsed
    envelope["status"] = "success"
    envelope["summary"]["conclusion"] = conclusion
    envelope["summary"]["key_metrics"] = {
        "sl_length": len(SL_SEQUENCE),
        "free_n_pairs": free_result.n_pairs,
        "free_mfe_estimate": free_result.mfe_estimate,
        "free_paired_fraction": free_result.paired_fraction,
        "bound_n_pairs": bound_result.n_pairs,
        "bound_mfe_estimate": bound_result.mfe_estimate,
        "target_accessibility": accessibility.mean_accessibility,
        "target_classification": accessibility.classification,
        "pairs_disrupted": len(disrupted_pairs),
        "pairs_new": len(new_pairs),
        "pairs_maintained": len(maintained_pairs),
    }
    envelope["data"] = {
        "sl_sequence": SL_SEQUENCE,
        "aso_sequence": ASO_SEQUENCE,
        "aso_binding_region": f"pos {ASO_TARGET_START}-{ASO_TARGET_END - 1}",
        "free_structure": {
            "dot_bracket": free_result.dot_bracket,
            "n_pairs": free_result.n_pairs,
            "score": free_result.score,
            "mfe_estimate": free_result.mfe_estimate,
            "paired_fraction": free_result.paired_fraction,
        },
        "bound_structure": {
            "dot_bracket": bound_result.dot_bracket,
            "n_pairs": bound_result.n_pairs,
            "score": bound_result.score,
            "mfe_estimate": bound_result.mfe_estimate,
            "paired_fraction": bound_result.paired_fraction,
        },
        "accessibility": {
            "target_start": accessibility.target_start,
            "target_end": accessibility.target_end,
            "region_dot_bracket": accessibility.region_dot_bracket,
            "per_position_scores": accessibility.per_position_scores,
            "mean_accessibility": accessibility.mean_accessibility,
            "classification": accessibility.classification,
            "n_accessible": accessibility.n_accessible,
            "n_blocked": accessibility.n_blocked,
            "alternative_windows": accessibility.alternative_windows,
        },
        "structural_impact": {
            "pairs_disrupted": len(disrupted_pairs),
            "pairs_disrupted_in_aso_region": len(disrupted_in_aso),
            "pairs_new_rearrangement": len(new_pairs),
            "pairs_maintained": len(maintained_pairs),
            "score_change": round(bound_result.score - free_result.score, 4),
            "mfe_change": round(bound_result.mfe_estimate - free_result.mfe_estimate, 2),
        },
    }
    write_result(envelope)
    return envelope


if __name__ == "__main__":
    main()
