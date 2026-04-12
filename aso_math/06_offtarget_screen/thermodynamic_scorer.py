"""Scoring termodinamico de hits off-target.

Para cada hit de complementaridade encontrado pelo aligner, calcula
dG e Tm da ligacao parcial e classifica o nivel de risco.

Niveis de risco (baseados em Crooke et al. 2017):
  - safe:    dG > -15 kcal/mol (ligacao fraca, sem atividade RNase H)
  - monitor: -20 < dG <= -15 kcal/mol (ligacao moderada, requer monitoramento)
  - danger:  dG <= -20 kcal/mol (ligacao forte, off-target provavel)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from aso_math.thermo import compute_dg, compute_tm


# ---------------------------------------------------------------------------
# Constantes de classificacao de risco
# ---------------------------------------------------------------------------

_DG_SAFE_THRESHOLD: float = -15.0      # acima disso: seguro
_DG_MONITOR_THRESHOLD: float = -20.0   # entre -20 e -15: monitorar


# ---------------------------------------------------------------------------
# Tipos
# ---------------------------------------------------------------------------


@dataclass
class ScoredHit:
    """Hit off-target com score termodinamico e classificacao de risco."""

    transcript_id: str
    match_length: int
    position: int
    matched_region: str
    dg_kcal: float
    tm_celsius: float
    risk_level: str  # "safe", "monitor", "danger"

    def to_dict(self) -> dict[str, Any]:
        """Serializa para dicionario JSON-compativel."""
        return {
            "transcript_id": self.transcript_id,
            "match_length": self.match_length,
            "position": self.position,
            "matched_region": self.matched_region,
            "dg_kcal": self.dg_kcal,
            "tm_celsius": self.tm_celsius,
            "risk_level": self.risk_level,
        }


# ---------------------------------------------------------------------------
# Classificacao de risco
# ---------------------------------------------------------------------------


def classify_risk(dg: float) -> str:
    """Classifica risco de um hit off-target baseado em dG.

    Args:
        dg: delta G de ligacao em kcal/mol (valores mais negativos = mais forte).

    Returns:
        "safe", "monitor", ou "danger".
    """
    if dg > _DG_SAFE_THRESHOLD:
        return "safe"
    if dg > _DG_MONITOR_THRESHOLD:
        return "monitor"
    return "danger"


# ---------------------------------------------------------------------------
# Scoring de hits
# ---------------------------------------------------------------------------


def _compute_best_subseq_dg(aso_seq: str, match_length: int) -> float:
    """Calcula o dG mais negativo entre todas as substrings do ASO com dado comprimento.

    Isso e uma estimativa conservadora: assume que o off-target alinha
    com a regiao do ASO que produz a ligacao mais forte.

    Args:
        aso_seq: Sequencia completa do ASO.
        match_length: Comprimento do trecho complementar.

    Returns:
        delta_G em kcal/mol (valor mais negativo = pior caso).
    """
    if match_length <= 0:
        return 0.0

    aso_upper = aso_seq.upper()
    aso_len = len(aso_upper)

    if match_length >= aso_len:
        return compute_dg(aso_upper)

    best_dg = 0.0
    for start in range(aso_len - match_length + 1):
        subseq = aso_upper[start : start + match_length]
        dg = compute_dg(subseq)
        if dg < best_dg:
            best_dg = dg

    return best_dg


def _compute_best_subseq_tm(aso_seq: str, match_length: int) -> float:
    """Calcula o Tm mais alto entre todas as substrings do ASO com dado comprimento.

    Args:
        aso_seq: Sequencia completa do ASO.
        match_length: Comprimento do trecho complementar.

    Returns:
        Tm em graus Celsius (valor mais alto = pior caso).
    """
    if match_length <= 0:
        return 0.0

    aso_upper = aso_seq.upper()
    aso_len = len(aso_upper)

    if match_length >= aso_len:
        return compute_tm(aso_upper)

    best_tm = 0.0
    for start in range(aso_len - match_length + 1):
        subseq = aso_upper[start : start + match_length]
        tm = compute_tm(subseq)
        if tm > best_tm:
            best_tm = tm

    return best_tm


def score_hits(
    hits: list[Any],
    aso_seq: str,
) -> list[ScoredHit]:
    """Score e classifica uma lista de hits off-target.

    Para cada hit, calcula dG e Tm da ligacao parcial usando o pior caso
    (substring do ASO com a ligacao mais forte para o comprimento do match).

    Args:
        hits: Lista de OffTargetHit (ou objetos com transcript_id,
              match_length, position, matched_region).
        aso_seq: Sequencia completa do ASO (5'->3').

    Returns:
        Lista de ScoredHit ordenada por dG crescente (mais perigoso primeiro).
    """
    scored: list[ScoredHit] = []

    for hit in hits:
        dg = _compute_best_subseq_dg(aso_seq, hit.match_length)
        tm = _compute_best_subseq_tm(aso_seq, hit.match_length)
        risk = classify_risk(dg)

        scored.append(
            ScoredHit(
                transcript_id=hit.transcript_id,
                match_length=hit.match_length,
                position=hit.position,
                matched_region=hit.matched_region,
                dg_kcal=round(dg, 2),
                tm_celsius=round(tm, 2),
                risk_level=risk,
            )
        )

    scored.sort(key=lambda s: s.dg_kcal)
    return scored
