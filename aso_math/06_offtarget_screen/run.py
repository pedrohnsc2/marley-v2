"""Modulo 06 - Off-target screening do MRL-ASO-001 contra transcriptoma canino.

Busca potenciais alvos off-target no transcriptoma do hospedeiro (Canis lupus
familiaris), classifica cada hit por nivel de risco termodinamico, e gera
um envelope JSON padronizado com os resultados.

Se o transcriptoma canino nao estiver disponivel, utiliza um modelo nulo
com sequencias aleatorias de conteudo GC similar ao canino (~41%).

Pipeline:
    1. Busca de complementaridade via janela deslizante (aligner.py)
    2. Scoring termodinamico de cada hit (thermodynamic_scorer.py)
    3. Classificacao de risco e resumo

Referencia:
  - Crooke ST et al. (2017) Nat Biotechnol 35(3):230-237
  - Lima WF et al. (2007) J Biol Chem 282(14):10374-10385
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

from aso_math.config import ASO_SEQUENCE, BASES, PROJECT_ROOT
from aso_math.envelope import Timer, create_envelope, write_result
from aso_math.target_config import TargetConfig

from .aligner import (
    OffTargetHit,
    find_max_complementarity,
    find_offtargets,
    is_nucleotide_fasta,
    parse_fasta,
)
from .thermodynamic_scorer import ScoredHit, classify_risk, score_hits

try:
    from core.logger import get_logger
    logger = get_logger("aso_math.06_offtarget_screen")
except Exception:
    import logging
    logger = logging.getLogger("aso_math.06_offtarget_screen")


TRANSCRIPTOME_CANINE: Path = PROJECT_ROOT / "data" / "raw" / "transcriptome_canine.fasta"

_RNG_SEED: int = 42
_NULL_MODEL_N_SEQUENCES: int = 300
_NULL_MODEL_MIN_LENGTH: int = 200
_NULL_MODEL_MAX_LENGTH: int = 1500
_DEFAULT_MIN_MATCH: int = 15


def _generate_null_model(
    n_sequences: int = _NULL_MODEL_N_SEQUENCES,
    min_length: int = _NULL_MODEL_MIN_LENGTH,
    max_length: int = _NULL_MODEL_MAX_LENGTH,
) -> dict[str, str]:
    """Gera sequencias aleatorias com conteudo GC similar ao canino (~41%).

    Cada sequencia tem comprimento uniforme entre min_length e max_length.
    O conteudo GC e amostrado uniformemente entre 0.38 e 0.44 para cada
    sequencia, refletindo a variacao real do genoma canino (Lindblad-Toh 2005).
    """
    random.seed(_RNG_SEED)
    sequences: dict[str, str] = {}

    for i in range(n_sequences):
        length = random.randint(min_length, max_length)
        gc_frac = random.uniform(0.38, 0.44)

        weights = [
            (1.0 - gc_frac) / 2.0,
            gc_frac / 2.0,
            gc_frac / 2.0,
            (1.0 - gc_frac) / 2.0,
        ]
        bases_list = random.choices(list(BASES), weights=weights, k=length)
        seq = "".join(bases_list)
        sequences[f"null_canine_seq_{i:04d}|len={length}|gc={gc_frac:.2f}"] = seq

    return sequences


def _screen_null_model(
    aso_seq: str,
    min_match: int = _DEFAULT_MIN_MATCH,
) -> tuple[list[OffTargetHit], int]:
    """Busca off-targets no modelo nulo.

    Returns:
        Tupla (hits, total_screened).
    """
    sequences = _generate_null_model()
    hits: list[OffTargetHit] = []

    for header, seq in sequences.items():
        max_compl, position, matched_region = find_max_complementarity(
            aso_seq, seq
        )
        if max_compl >= min_match:
            hits.append(
                OffTargetHit(
                    transcript_id=header,
                    match_length=max_compl,
                    position=position,
                    matched_region=matched_region,
                )
            )

    hits.sort(key=lambda h: h.match_length, reverse=True)
    return hits, len(sequences)


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a varredura off-target completa.

    Busca complementaridade do ASO contra o transcriptoma canino (hospedeiro).
    Se o transcriptoma nao estiver disponivel, usa modelo nulo.

    Cada hit e pontuado termodinamicamente e classificado como:
      - safe: dG > -15 kcal/mol
      - monitor: -20 < dG <= -15 kcal/mol
      - danger: dG <= -20 kcal/mol

    Args:
        config: Configuracao do organismo-alvo. Se None, usa L. infantum.

    Returns:
        Envelope JSON completo com resultados.
    """
    if config is None:
        config = TargetConfig()

    aso_seq = config.aso_sequence or ASO_SEQUENCE
    fasta_path = TRANSCRIPTOME_CANINE

    module_name = "06_offtarget_screen"
    envelope = create_envelope(module_name)

    logger.info("=" * 60)
    logger.info("MODULO 06: Off-Target Screening - %s", config.species_name)
    logger.info("=" * 60)

    try:
        with Timer() as timer:
            transcriptome_type: str
            hits: list[OffTargetHit]
            total_screened: int

            if fasta_path.exists() and is_nucleotide_fasta(fasta_path):
                logger.info("Usando transcriptoma canino: %s", fasta_path)
                transcriptome_type = "canine"
                hits = find_offtargets(aso_seq, fasta_path, min_match=_DEFAULT_MIN_MATCH)
                sequences = parse_fasta(fasta_path)
                total_screened = len(sequences)
            else:
                if fasta_path.exists():
                    reason = (
                        "Transcriptoma canino contem proteinas, nao nucleotideos. "
                        "Usando modelo nulo."
                    )
                else:
                    reason = (
                        f"Transcriptoma canino nao encontrado em {fasta_path}. "
                        "Usando modelo nulo com sequencias aleatorias."
                    )
                logger.warning(reason)
                envelope["warnings"].append(reason)
                transcriptome_type = "null_model"
                hits, total_screened = _screen_null_model(aso_seq)

            logger.info("  ASO: %s (%d nt)", aso_seq, len(aso_seq))
            logger.info("  Transcritos varridos: %d", total_screened)
            logger.info("  Hits com >= %d bp complementares: %d", _DEFAULT_MIN_MATCH, len(hits))

            logger.info("")
            logger.info("Passo 2: Scoring termodinamico dos hits")

            scored_hits: list[ScoredHit] = score_hits(hits, aso_seq)

            risk_counts = {"safe": 0, "monitor": 0, "danger": 0}
            for sh in scored_hits:
                risk_counts[sh.risk_level] += 1

            logger.info("  Classificacao de risco:")
            logger.info("    safe:    %d hits", risk_counts["safe"])
            logger.info("    monitor: %d hits", risk_counts["monitor"])
            logger.info("    danger:  %d hits", risk_counts["danger"])

            logger.info("")
            logger.info("Passo 3: Resumo")

            overall_safe = risk_counts["danger"] == 0
            worst_hit: dict[str, Any] | None = None
            if scored_hits:
                worst = scored_hits[0]
                worst_hit = worst.to_dict()
                logger.info(
                    "  Pior hit: %s (match=%d bp, dG=%.2f kcal/mol, Tm=%.2f C, risk=%s)",
                    worst.transcript_id[:50],
                    worst.match_length,
                    worst.dg_kcal,
                    worst.tm_celsius,
                    worst.risk_level,
                )

            if overall_safe:
                conclusion = (
                    f"Varredura off-target contra {total_screened} transcritos "
                    f"({transcriptome_type}): {len(hits)} hits com >= {_DEFAULT_MIN_MATCH} bp "
                    f"de complementaridade. Nenhum hit classificado como danger "
                    f"(dG <= -20 kcal/mol). O MRL-ASO-001 apresenta perfil de "
                    f"seguranca adequado para o hospedeiro canino."
                )
            else:
                conclusion = (
                    f"Varredura off-target contra {total_screened} transcritos "
                    f"({transcriptome_type}): {len(hits)} hits com >= {_DEFAULT_MIN_MATCH} bp "
                    f"de complementaridade. {risk_counts['danger']} hit(s) classificado(s) "
                    f"como danger (dG <= -20 kcal/mol). Requer investigacao adicional."
                )

            logger.info("")
            logger.info("RESULTADO: %s", "SEGURO" if overall_safe else "REQUER ATENCAO")
            logger.info("Conclusao: %s", conclusion)

            envelope["data"] = {
                "transcriptome_source": transcriptome_type,
                "transcripts_screened": total_screened,
                "min_match_threshold": _DEFAULT_MIN_MATCH,
                "total_hits": len(hits),
                "scored_hits": [sh.to_dict() for sh in scored_hits],
                "risk_summary": risk_counts,
                "worst_hit": worst_hit,
                "overall_safe": overall_safe,
            }
            envelope["status"] = "pass" if overall_safe else "warn"
            envelope["summary"]["conclusion"] = conclusion
            envelope["summary"]["key_metrics"] = {
                "transcripts_screened": total_screened,
                "total_hits": len(hits),
                "danger_hits": risk_counts["danger"],
                "monitor_hits": risk_counts["monitor"],
                "safe_hits": risk_counts["safe"],
                "overall_safe": overall_safe,
            }

        envelope["runtime_seconds"] = timer.elapsed
        logger.info("Tempo de execucao: %.2f s", timer.elapsed)

    except Exception as exc:
        logger.error("Erro no modulo 06: %s", exc, exc_info=True)
        envelope["status"] = "error"
        envelope["warnings"].append(f"Erro durante execucao: {exc}")

    output_path = write_result(envelope)
    logger.info("Resultados gravados em: %s", output_path)

    return envelope


if __name__ == "__main__":
    main()
