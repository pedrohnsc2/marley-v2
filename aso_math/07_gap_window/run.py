"""Modulo 07 -- Otimizacao de janela de gap para design gapmer LNA-DNA-LNA.

Objetivo: enumerar TODAS as configuracoes validas de gapmer para MRL-ASO-001
(25 nt), modelar a eficiencia de RNase H para cada gap, e rankear por
score composto integrando termodinamica, clivagem e resistencia a nucleases.

Design atual de MRL-ASO-001: LNA 4 + DNA 17 + LNA 4 (4-17-4).
"""

from __future__ import annotations

from typing import Any

from aso_math.config import (
    ASO_SEQUENCE,
    SL_SEQUENCE,
)
from aso_math.envelope import Timer, create_envelope, write_result
from core.logger import get_logger

from .gap_enumerator import GapWindowConfig, enumerate_gap_configs
from .rnase_h_model import ScoredConfig, score_config

logger = get_logger("07_gap_window")

# Configuracao LNA conhecida de MRL-ASO-001
MRL_LNA_5P: int = 4
MRL_DNA_GAP: int = 17
MRL_LNA_3P: int = 4


def _scored_to_dict(scored: ScoredConfig, rank: int) -> dict[str, Any]:
    """Converte ScoredConfig para dicionario serializavel."""
    return {
        "rank": rank,
        "notation": scored.config.notation(),
        "lna_5p": scored.config.lna_5p,
        "dna_gap": scored.config.dna_gap,
        "lna_3p": scored.config.lna_3p,
        "total_lna": scored.config.total_lna,
        "dg_binding_kcal": scored.dg_binding,
        "tm_base_celsius": scored.tm_base,
        "tm_adjusted_celsius": scored.tm_adjusted,
        "rnase_h_efficiency": scored.rnase_h_efficiency,
        "nuclease_resistance": scored.nuclease_resistance,
        "composite_score": scored.composite_score,
    }


def _find_current_design(
    scored_list: list[ScoredConfig],
    lna_5p: int = MRL_LNA_5P,
    dna_gap: int = MRL_DNA_GAP,
    lna_3p: int = MRL_LNA_3P,
) -> tuple[ScoredConfig | None, int]:
    """Localiza o design atual nos resultados rankeados.

    Returns:
        Tupla (ScoredConfig, rank_1_indexed) ou (None, -1) se nao encontrado.
    """
    for i, scored in enumerate(scored_list):
        if (
            scored.config.lna_5p == lna_5p
            and scored.config.dna_gap == dna_gap
            and scored.config.lna_3p == lna_3p
        ):
            return scored, i + 1
    return None, -1


def main(
    aso_seq: str | None = None,
    sl_seq: str | None = None,
) -> dict[str, Any]:
    """Executa a otimizacao de janela de gap para o ASO especificado.

    Fluxo:
        1. Enumerar todas as configuracoes LNA-DNA-LNA validas para 25 nt
        2. Calcular score composto para cada configuracao
        3. Rankear por score decrescente
        4. Localizar o design atual (4-17-4) e reportar status
        5. Gravar resultados em JSON

    Args:
        aso_seq: Sequencia do ASO. Se None, usa ASO_SEQUENCE do config.
        sl_seq: Sequencia do SL RNA alvo. Se None, usa SL_SEQUENCE do config.

    Returns:
        Envelope JSON completo com resultados.
    """
    _aso_seq = (aso_seq or ASO_SEQUENCE).upper()
    _sl_seq = (sl_seq or SL_SEQUENCE).upper()
    aso_length = len(_aso_seq)

    logger.info("=" * 60)
    logger.info("MODULO 07: Otimizacao de Janela de Gap (LNA-DNA-LNA)")
    logger.info("=" * 60)
    logger.info("ASO: %s (%d nt)", _aso_seq, aso_length)

    envelope = create_envelope("07_gap_window")

    with Timer() as timer:
        # --- Passo 1: Enumerar configuracoes validas ---
        logger.info(
            "Passo 1/4: Enumerando configuracoes gapmer para %d nt...",
            aso_length,
        )
        configs = enumerate_gap_configs(aso_length=aso_length)
        logger.info("  %d configuracoes validas encontradas.", len(configs))

        if not configs:
            envelope["status"] = "failed"
            envelope["summary"]["conclusion"] = (
                f"Nenhuma configuracao gapmer valida para ASO de {aso_length} nt."
            )
            output_path = write_result(envelope)
            logger.warning("Nenhuma configuracao valida. Resultado: %s", output_path)
            return envelope

        # --- Passo 2: Calcular scores ---
        logger.info("Passo 2/4: Calculando scores compostos...")
        scored_list: list[ScoredConfig] = []
        for config in configs:
            scored = score_config(config, _aso_seq, _sl_seq)
            scored_list.append(scored)

        # Rankear por score decrescente
        scored_list.sort(key=lambda s: -s.composite_score)

        logger.info(
            "  Top 3: %s",
            [
                f"{s.config.notation()} (score={s.composite_score:.4f})"
                for s in scored_list[:3]
            ],
        )

        # --- Passo 3: Construir ranking table ---
        logger.info("Passo 3/4: Construindo tabela de ranking...")
        ranking_table: list[dict[str, Any]] = []
        for rank, scored in enumerate(scored_list, start=1):
            ranking_table.append(_scored_to_dict(scored, rank))

        # --- Passo 4: Comparar com design atual ---
        logger.info("Passo 4/4: Comparando com design atual (4-17-4)...")
        current_design, current_rank = _find_current_design(scored_list)

        if current_design is not None:
            best = scored_list[0]
            if current_rank == 1:
                conclusion = (
                    f"Design atual ({current_design.config.notation()}) "
                    f"CONFIRMADO como configuracao otima entre {len(configs)} "
                    f"gapmers validos. Score: {current_design.composite_score:.4f}."
                )
            else:
                score_diff = round(best.composite_score - current_design.composite_score, 4)
                conclusion = (
                    f"Design atual ({current_design.config.notation()}) "
                    f"esta no rank #{current_rank} de {len(configs)} "
                    f"(score {current_design.composite_score:.4f}). "
                    f"Melhor configuracao: {best.config.notation()} "
                    f"com score {best.composite_score:.4f} "
                    f"(+{score_diff} sobre o atual)."
                )

            logger.info(
                "  Design atual: rank #%d, score %.4f, RNase H %.2f",
                current_rank,
                current_design.composite_score,
                current_design.rnase_h_efficiency,
            )
        else:
            best = scored_list[0]
            conclusion = (
                f"Design atual (4-17-4) NAO encontrado no espaco de busca. "
                f"Melhor configuracao: {best.config.notation()} "
                f"com score {best.composite_score:.4f}."
            )
            logger.warning("Design atual nao encontrado no espaco de busca!")

        logger.info("CONCLUSAO: %s", conclusion)

    # --- Montar envelope ---
    envelope["runtime_seconds"] = timer.elapsed
    envelope["status"] = "success"
    envelope["summary"]["conclusion"] = conclusion
    envelope["summary"]["key_metrics"] = {
        "total_configs_evaluated": len(configs),
        "best_config": scored_list[0].config.notation() if scored_list else "N/A",
        "best_score": scored_list[0].composite_score if scored_list else 0.0,
        "current_design_rank": current_rank,
        "current_design_score": (
            current_design.composite_score if current_design else 0.0
        ),
    }

    # Dados de comparacao do design atual
    current_comparison: dict[str, Any] = {}
    if current_design is not None:
        current_comparison = {
            "found": True,
            "rank": current_rank,
            **_scored_to_dict(current_design, current_rank),
        }
    else:
        current_comparison = {"found": False, "rank": -1}

    envelope["data"] = {
        "aso_sequence": _aso_seq,
        "aso_length": aso_length,
        "sl_sequence": _sl_seq,
        "enumeration": {
            "total_valid_configs": len(configs),
            "constraints": {
                "lna_5p_range": [2, 5],
                "lna_3p_range": [2, 5],
                "dna_gap_min": 8,
                "total_lna_max": 10,
            },
        },
        "ranking": ranking_table,
        "current_design_comparison": current_comparison,
        "scoring_weights": {
            "w_dg": 0.30,
            "w_tm": 0.20,
            "w_rnase_h": 0.30,
            "w_nuclease": 0.20,
        },
    }

    # Gravar resultado
    output_path = write_result(envelope)
    logger.info("Resultado gravado em: %s", output_path)
    logger.info("Tempo de execucao: %.2f segundos", timer.elapsed)

    return envelope


if __name__ == "__main__":
    main()
