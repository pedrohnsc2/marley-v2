"""Predicao de estrutura secundaria de RNA via algoritmo de Nussinov.

Implementa o algoritmo classico de programacao dinamica para maximizar
o numero de pares de bases em uma sequencia de RNA, respeitando regras
de pareamento Watson-Crick (AU, GC) e wobble (GU).

A predicao de estrutura e usada para:
    1. Determinar a conformacao do SL RNA de L. infantum em estado livre
    2. Simular o efeito da ligacao do MRL-ASO-001 (bloqueio de posicoes)
    3. Quantificar a disrupcao estrutural causada pelo ASO

Ref: Nussinov R et al. (1978) SIAM J Appl Math 35(1):68-82
Ref: Zuker M & Stiegler P (1981) Nucleic Acids Res 9(1):133-148
"""

from __future__ import annotations

import numpy as np

from .config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    BASE_PAIR_RULES,
    MIN_LOOP_LENGTH,
    SL_SEQUENCES,
)


# ---------------------------------------------------------------------------
# Algoritmo de Nussinov
# ---------------------------------------------------------------------------

def can_pair(nuc_i: str, nuc_j: str) -> bool:
    """Verifica se dois nucleotideos podem parear.

    Args:
        nuc_i: Nucleotideo na posicao i (A, U, G, C).
        nuc_j: Nucleotideo na posicao j (A, U, G, C).

    Returns:
        True se o par e valido (Watson-Crick ou wobble).
    """
    return (nuc_i, nuc_j) in BASE_PAIR_RULES


def pair_score(nuc_i: str, nuc_j: str) -> float:
    """Retorna o score de pareamento entre dois nucleotideos.

    Watson-Crick (AU, GC) = 1.0, wobble (GU) = 0.5, invalido = 0.0.

    Args:
        nuc_i: Nucleotideo na posicao i.
        nuc_j: Nucleotideo na posicao j.

    Returns:
        Score de pareamento (0.0, 0.5 ou 1.0).
    """
    return BASE_PAIR_RULES.get((nuc_i, nuc_j), 0.0)


def nussinov_fill(
    sequence: str,
    blocked_positions: set[int] | None = None,
) -> np.ndarray:
    """Preenche a matriz de programacao dinamica do Nussinov.

    Calcula o numero maximo de pares de bases para toda subsequencia [i,j].
    Posicoes bloqueadas (ex: ligadas ao ASO) nao podem formar pares.

    Recorrencia:
        dp[i][j] = max(
            dp[i+1][j],                          # i nao pareado
            dp[i][j-1],                          # j nao pareado
            dp[i+1][j-1] + score(i,j),           # i-j pareados
            max_k(dp[i][k] + dp[k+1][j])         # bifurcacao
        )

    Args:
        sequence: Sequencia de RNA (ex: "AACUAACGCUAUAUAAG...").
        blocked_positions: Conjunto de posicoes que nao podem parear
                          (simulacao de ligacao do ASO).

    Returns:
        Matriz dp (n x n) com scores otimos por subsequencia.
    """
    n = len(sequence)
    dp = np.zeros((n, n), dtype=float)
    blocked = blocked_positions or set()

    # Preenche diagonais crescentes (comprimento da subsequencia)
    for length in range(MIN_LOOP_LENGTH + 1, n):
        for i in range(n - length):
            j = i + length

            # Opcao 1: posicao i nao pareada
            dp[i][j] = dp[i + 1][j]

            # Opcao 2: posicao j nao pareada
            dp[i][j] = max(dp[i][j], dp[i][j - 1])

            # Opcao 3: i e j pareados (se permitido)
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

            # Opcao 4: bifurcacao (dividir em duas sub-estruturas)
            for k in range(i + 1, j):
                dp[i][j] = max(dp[i][j], dp[i][k] + dp[k + 1][j])

    return dp


def nussinov_traceback(
    dp: np.ndarray,
    sequence: str,
    blocked_positions: set[int] | None = None,
) -> list[tuple[int, int]]:
    """Traceback na matriz de Nussinov para recuperar pares de bases.

    Reconstroi a estrutura otima seguindo as decisoes registradas na
    matriz de programacao dinamica.

    Args:
        dp: Matriz preenchida pelo nussinov_fill.
        sequence: Sequencia de RNA original.
        blocked_positions: Posicoes bloqueadas (devem ser as mesmas do fill).

    Returns:
        Lista de tuplas (i, j) representando pares de bases.
    """
    n = len(sequence)
    blocked = blocked_positions or set()
    pairs: list[tuple[int, int]] = []

    # Stack para traceback iterativo (evita recursao profunda)
    stack: list[tuple[int, int]] = [(0, n - 1)]

    while stack:
        i, j = stack.pop()

        if i >= j or j - i < MIN_LOOP_LENGTH:
            continue

        # Caso 1: i nao pareado
        if dp[i][j] == dp[i + 1][j]:
            stack.append((i + 1, j))
            continue

        # Caso 2: j nao pareado
        if dp[i][j] == dp[i][j - 1]:
            stack.append((i, j - 1))
            continue

        # Caso 3: i-j pareados
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

        # Caso 4: bifurcacao
        for k in range(i + 1, j):
            if abs(dp[i][j] - (dp[i][k] + dp[k + 1][j])) < 1e-9:
                stack.append((i, k))
                stack.append((k + 1, j))
                break

    return sorted(pairs)


def pairs_to_dotbracket(pairs: list[tuple[int, int]], length: int) -> str:
    """Converte lista de pares de bases para notacao dot-bracket.

    Na notacao dot-bracket:
        - '(' indica base pareada (abertura)
        - ')' indica base pareada (fechamento)
        - '.' indica base nao-pareada

    Args:
        pairs: Lista de tuplas (i, j) com i < j.
        length: Comprimento total da sequencia.

    Returns:
        String em notacao dot-bracket (ex: "..((....))..").
    """
    structure = ["."] * length
    for i, j in pairs:
        structure[i] = "("
        structure[j] = ")"
    return "".join(structure)


# ---------------------------------------------------------------------------
# Funcoes de analise de estrutura do SL RNA
# ---------------------------------------------------------------------------

def predict_structure(
    sequence: str,
    blocked_positions: set[int] | None = None,
) -> dict[str, object]:
    """Prediz estrutura secundaria de uma sequencia de RNA.

    Executa o algoritmo de Nussinov completo: fill + traceback + dot-bracket.

    Args:
        sequence: Sequencia de RNA.
        blocked_positions: Posicoes bloqueadas (opcional).

    Returns:
        Dicionario com:
            - "sequence": sequencia original
            - "dot_bracket": notacao dot-bracket
            - "pairs": lista de pares de bases
            - "n_pairs": numero total de pares
            - "score": score otimo da matriz dp
            - "paired_fraction": fracao de posicoes pareadas
    """
    # Normaliza: DNA -> RNA
    seq_rna = sequence.upper().replace("T", "U")

    dp = nussinov_fill(seq_rna, blocked_positions)
    pairs = nussinov_traceback(dp, seq_rna, blocked_positions)
    dot_bracket = pairs_to_dotbracket(pairs, len(seq_rna))

    # Calcula fracao de posicoes pareadas
    paired_positions = set()
    for i, j in pairs:
        paired_positions.add(i)
        paired_positions.add(j)
    paired_fraction = len(paired_positions) / len(seq_rna) if seq_rna else 0.0

    return {
        "sequence": seq_rna,
        "dot_bracket": dot_bracket,
        "pairs": pairs,
        "n_pairs": len(pairs),
        "score": float(dp[0][len(seq_rna) - 1]) if seq_rna else 0.0,
        "paired_fraction": round(paired_fraction, 4),
    }


def analyze_aso_structural_impact() -> dict[str, object]:
    """Analisa o impacto estrutural da ligacao do MRL-ASO-001 ao SL RNA.

    Compara a estrutura secundaria do SL RNA de L. infantum em dois estados:
        1. Livre (sem ASO) — conformacao nativa
        2. Ligado ao ASO — posicoes 5-29 bloqueadas

    Quantifica a disrupcao estrutural: pares perdidos, ganhos, e alteracao
    da forca de dobragem.

    Returns:
        Dicionario com resultados comparativos.
    """
    sl_seq = SL_SEQUENCES["L_infantum"]

    # Estado 1: SL RNA livre (conformacao nativa)
    free_structure = predict_structure(sl_seq)

    # Estado 2: SL RNA com ASO ligado (posicoes 5-29 bloqueadas)
    # O ASO se liga por complementaridade, impedindo pareamento intramolecular
    aso_blocked = set(range(ASO_TARGET_START, ASO_TARGET_END))
    bound_structure = predict_structure(sl_seq, blocked_positions=aso_blocked)

    # Analise comparativa: pares perdidos e ganhos
    free_pairs_set = set(free_structure["pairs"])
    bound_pairs_set = set(bound_structure["pairs"])
    disrupted_pairs = free_pairs_set - bound_pairs_set  # Pares perdidos
    new_pairs = bound_pairs_set - free_pairs_set         # Pares novos (rearranjo)
    maintained_pairs = free_pairs_set & bound_pairs_set  # Pares mantidos

    # Posicoes cujo pareamento foi diretamente afetado pelo ASO
    disrupted_in_aso_region = [
        (i, j) for i, j in disrupted_pairs
        if i in aso_blocked or j in aso_blocked
    ]

    return {
        "sl_sequence": sl_seq,
        "aso_sequence": ASO_SEQUENCE,
        "aso_binding_region": f"pos {ASO_TARGET_START}-{ASO_TARGET_END - 1}",
        "free_structure": {
            "dot_bracket": free_structure["dot_bracket"],
            "n_pairs": free_structure["n_pairs"],
            "score": free_structure["score"],
            "paired_fraction": free_structure["paired_fraction"],
        },
        "bound_structure": {
            "dot_bracket": bound_structure["dot_bracket"],
            "n_pairs": bound_structure["n_pairs"],
            "score": bound_structure["score"],
            "paired_fraction": bound_structure["paired_fraction"],
        },
        "structural_impact": {
            "pairs_disrupted": len(disrupted_pairs),
            "pairs_disrupted_in_aso_region": len(disrupted_in_aso_region),
            "pairs_new_rearrangement": len(new_pairs),
            "pairs_maintained": len(maintained_pairs),
            "score_change": round(
                bound_structure["score"] - free_structure["score"], 4
            ),
            "paired_fraction_change": round(
                bound_structure["paired_fraction"] - free_structure["paired_fraction"], 4
            ),
        },
        "interpretation": _interpret_structural_impact(
            free_structure, bound_structure, len(disrupted_pairs)
        ),
    }


def _interpret_structural_impact(
    free: dict, bound: dict, n_disrupted: int
) -> str:
    """Gera interpretacao textual do impacto estrutural do ASO.

    Args:
        free: Resultado da predicao de estrutura livre.
        bound: Resultado da predicao de estrutura com ASO.
        n_disrupted: Numero de pares disrupted.

    Returns:
        String com interpretacao biologica.
    """
    if n_disrupted == 0:
        return (
            "O ASO nao disrupta pares de bases previstos na regiao de ligacao. "
            "A regiao alvo ja e predominantemente single-stranded (nao pareada), "
            "o que favorece a acessibilidade do ASO."
        )

    total_free = free["n_pairs"]
    pct_disrupted = (n_disrupted / total_free * 100) if total_free > 0 else 0

    if pct_disrupted > 50:
        severity = "massiva"
    elif pct_disrupted > 25:
        severity = "significativa"
    else:
        severity = "moderada"

    return (
        f"O ASO causa disrupcao {severity} da estrutura secundaria: "
        f"{n_disrupted}/{total_free} pares ({pct_disrupted:.1f}%) sao perdidos. "
        f"Isso indica que o ASO compete diretamente com a estrutura intramolecular "
        f"do SL RNA, desestabilizando a conformacao nativa e potencialmente "
        f"expondo o RNA para degradacao por RNase H."
    )
