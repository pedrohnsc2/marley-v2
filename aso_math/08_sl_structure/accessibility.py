"""Analise de acessibilidade da regiao-alvo do ASO no SL RNA.

Dado o dot-bracket de uma estrutura secundaria, avalia quais posicoes
estao disponiveis para ligacao de um oligonucleotideo antisense.

Logica:
    - '.' (nao-pareada) = acessivel (score 1.0)
    - '(' ou ')' (pareada) = bloqueada (score 0.0)

Ref: Ding Y & Lawrence CE (2003) Nucleic Acids Res 31(24):7280-7301
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final


ACCESSIBLE_THRESHOLD: Final[float] = 0.7
PARTIAL_THRESHOLD: Final[float] = 0.5
MIN_WINDOW_SIZE: Final[int] = 15


@dataclass(frozen=True)
class AccessibilityReport:
    """Relatorio de acessibilidade de uma regiao-alvo."""
    target_start: int
    target_end: int
    region_dot_bracket: str
    per_position_scores: list[float] = field(default_factory=list)
    mean_accessibility: float = 0.0
    classification: str = "blocked"
    n_accessible: int = 0
    n_blocked: int = 0
    alternative_windows: list[dict[str, object]] = field(default_factory=list)


def _score_position(char: str) -> float:
    """1.0 se nao-pareado ('.'), 0.0 se pareado."""
    return 1.0 if char == "." else 0.0


def _classify_accessibility(mean_score: float) -> str:
    """Classifica: accessible (>0.7), partially_blocked (0.5-0.7), blocked (<0.5)."""
    if mean_score > ACCESSIBLE_THRESHOLD:
        return "accessible"
    if mean_score >= PARTIAL_THRESHOLD:
        return "partially_blocked"
    return "blocked"


def _find_alternative_windows(
    dot_bracket: str, current_start: int, current_end: int,
    window_size: int | None = None,
) -> list[dict[str, object]]:
    """Busca janelas alternativas com maior acessibilidade."""
    n = len(dot_bracket)
    win = window_size or (current_end - current_start)
    if win < MIN_WINDOW_SIZE:
        win = MIN_WINDOW_SIZE
    if win > n:
        return []
    current_scores = [_score_position(dot_bracket[p]) for p in range(current_start, min(current_end, n))]
    current_mean = sum(current_scores) / len(current_scores) if current_scores else 0.0
    alternatives: list[dict[str, object]] = []
    for start in range(n - win + 1):
        end = start + win
        if start == current_start and end == current_end:
            continue
        scores = [_score_position(dot_bracket[p]) for p in range(start, end)]
        mean_acc = sum(scores) / len(scores)
        if mean_acc > current_mean:
            alternatives.append({
                "start": start, "end": end,
                "mean_accessibility": round(mean_acc, 4),
                "classification": _classify_accessibility(mean_acc),
            })
    alternatives.sort(key=lambda w: w["mean_accessibility"], reverse=True)
    return alternatives[:5]


def analyze_accessibility(
    dot_bracket: str, target_start: int, target_end: int,
) -> AccessibilityReport:
    """Analisa a acessibilidade de uma regiao-alvo no dot-bracket."""
    n = len(dot_bracket)
    if target_start < 0 or target_end > n or target_start >= target_end:
        raise ValueError(
            f"Regiao invalida: start={target_start}, end={target_end}, "
            f"comprimento do dot-bracket={n}"
        )
    region_db = dot_bracket[target_start:target_end]
    per_position = [_score_position(c) for c in region_db]
    region_len = target_end - target_start
    n_accessible = sum(1 for s in per_position if s == 1.0)
    n_blocked = region_len - n_accessible
    mean_acc = sum(per_position) / region_len if region_len > 0 else 0.0
    classification = _classify_accessibility(mean_acc)
    alternatives: list[dict[str, object]] = []
    if classification != "accessible":
        alternatives = _find_alternative_windows(
            dot_bracket, target_start, target_end,
        )
    return AccessibilityReport(
        target_start=target_start, target_end=target_end,
        region_dot_bracket=region_db, per_position_scores=per_position,
        mean_accessibility=round(mean_acc, 4),
        classification=classification,
        n_accessible=n_accessible, n_blocked=n_blocked,
        alternative_windows=alternatives,
    )
