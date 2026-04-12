"""Predicao de estrutura secundaria de RNA via algoritmo de Nussinov.

Adaptado de marley_ai/04_rna_fm/structure.py para uso standalone no
suite aso_math, sem dependencias externas (numpy removido).

Ref: Nussinov R et al. (1978) SIAM J Appl Math 35(1):68-82
Ref: Zuker M & Stiegler P (1981) Nucleic Acids Res 9(1):133-148
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final


BASE_PAIR_RULES: Final[dict[tuple[str, str], float]] = {
    ("A", "U"): 1.0,
    ("U", "A"): 1.0,
    ("G", "C"): 1.0,
    ("C", "G"): 1.0,
    ("G", "U"): 0.5,
    ("U", "G"): 0.5,
}

MIN_LOOP_LENGTH: Final[int] = 4
MFE_SCALE_KCAL: Final[float] = -1.5


@dataclass(frozen=True)
class StructureResult:
    """Resultado da predicao de estrutura secundaria."""
    sequence: str
    dot_bracket: str
    pairs: list[tuple[int, int]] = field(default_factory=list)
    n_pairs: int = 0
    score: float = 0.0
    mfe_estimate: float = 0.0
    paired_fraction: float = 0.0


def can_pair(nuc_i: str, nuc_j: str) -> bool:
    """Verifica se dois nucleotideos podem parear (WC ou wobble)."""
    return (nuc_i, nuc_j) in BASE_PAIR_RULES


def pair_score(nuc_i: str, nuc_j: str) -> float:
    """Score de pareamento: WC=1.0, wobble=0.5, invalido=0.0."""
    return BASE_PAIR_RULES.get((nuc_i, nuc_j), 0.0)


def nussinov_fill(
    sequence: str,
    blocked_positions: set[int] | None = None,
) -> list[list[float]]:
    """Preenche a matriz DP do Nussinov."""
    n = len(sequence)
    dp: list[list[float]] = [[0.0] * n for _ in range(n)]
    blocked = blocked_positions or set()
    for length in range(MIN_LOOP_LENGTH + 1, n):
        for i in range(n - length):
            j = i + length
            dp[i][j] = dp[i + 1][j]
            dp[i][j] = max(dp[i][j], dp[i][j - 1])
            if (
                i not in blocked
                and j not in blocked
                and can_pair(sequence[i], sequence[j])
            ):
                score = pair_score(sequence[i], sequence[j])
                if i + 1 <= j - 1:
                    dp[i][j] = max(dp[i][j], dp[i + 1][j - 1] + score)
                else:
                    dp[i][j] = max(dp[i][j], score)
            for k in range(i + 1, j):
                dp[i][j] = max(dp[i][j], dp[i][k] + dp[k + 1][j])
    return dp


def nussinov_traceback(
    dp: list[list[float]],
    sequence: str,
    blocked_positions: set[int] | None = None,
) -> list[tuple[int, int]]:
    """Traceback na matriz de Nussinov para recuperar pares de bases."""
    n = len(sequence)
    blocked = blocked_positions or set()
    pairs: list[tuple[int, int]] = []
    stack: list[tuple[int, int]] = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j or j - i < MIN_LOOP_LENGTH:
            continue
        if dp[i][j] == dp[i + 1][j]:
            stack.append((i + 1, j))
            continue
        if dp[i][j] == dp[i][j - 1]:
            stack.append((i, j - 1))
            continue
        if (
            i not in blocked
            and j not in blocked
            and can_pair(sequence[i], sequence[j])
        ):
            score = pair_score(sequence[i], sequence[j])
            prev = dp[i + 1][j - 1] if i + 1 <= j - 1 else 0.0
            if abs(dp[i][j] - (prev + score)) < 1e-9:
                pairs.append((i, j))
                if i + 1 <= j - 1:
                    stack.append((i + 1, j - 1))
                continue
        for k in range(i + 1, j):
            if abs(dp[i][j] - (dp[i][k] + dp[k + 1][j])) < 1e-9:
                stack.append((i, k))
                stack.append((k + 1, j))
                break
    return sorted(pairs)


def pairs_to_dotbracket(pairs: list[tuple[int, int]], length: int) -> str:
    """Converte lista de pares de bases para notacao dot-bracket."""
    structure = ["."] * length
    for i, j in pairs:
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure)


def predict_structure(
    sequence: str,
    blocked_positions: set[int] | None = None,
) -> StructureResult:
    """Prediz estrutura secundaria de uma sequencia de RNA."""
    if not sequence:
        return StructureResult(
            sequence="", dot_bracket="", pairs=[],
            n_pairs=0, score=0.0, mfe_estimate=0.0, paired_fraction=0.0,
        )
    seq_rna = sequence.upper().replace("T", "U")
    dp = nussinov_fill(seq_rna, blocked_positions)
    pairs = nussinov_traceback(dp, seq_rna, blocked_positions)
    dot_bracket = pairs_to_dotbracket(pairs, len(seq_rna))
    score = dp[0][len(seq_rna) - 1] if len(seq_rna) > 1 else 0.0
    mfe_estimate = round(score * MFE_SCALE_KCAL, 2)
    paired_positions: set[int] = set()
    for i, j in pairs:
        paired_positions.add(i)
        paired_positions.add(j)
    paired_fraction = len(paired_positions) / len(seq_rna)
    return StructureResult(
        sequence=seq_rna, dot_bracket=dot_bracket, pairs=pairs,
        n_pairs=len(pairs), score=score, mfe_estimate=mfe_estimate,
        paired_fraction=round(paired_fraction, 4),
    )
