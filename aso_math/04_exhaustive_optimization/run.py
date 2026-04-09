"""Modulo 04 — Otimizacao exaustiva: enumeracao completa do espaco de design ASO.

Objetivo: enumerar TODAS as combinacoes possiveis de (posicao_inicio x comprimento
x configuracao_LNA x tamanho_gap) e provar que a configuracao de MRL-ASO-001 e
otima — ou encontrar honestamente uma melhor.

Espaco de busca:
    ~170 sequencias candidatas (janelas no SL RNA de 18-27 nt)
    x ~10 configuracoes LNA validas por comprimento (gapmer 5'-LNA / DNA gap / 3'-LNA)
    = ~1700 combinacoes totais

Por que enumeracao exaustiva e nao Bayesian Optimization:
    O espaco de ~1700 combinacoes e pequeno o suficiente para avaliacao exata
    em segundos. Enumeracao exaustiva fornece PROVA MATEMATICA de otimalidade
    (ou contra-exemplo), enquanto metodos estocasticos so fornecem estimativas.

Design gapmer LNA:
    LNA (Locked Nucleic Acid) nos flancos 5' e 3' aumentam Tm (~3 C por posicao,
    McTigue et al. 2004) e resistencia a nucleases. O gap central de DNA e
    necessario para ativacao de RNase H (minimo 8 nt, Crooke et al. 2017).

MRL-ASO-001 atual: 25 nt, 4 LNA 5' + 17 DNA gap + 4 LNA 3'.
"""

from __future__ import annotations

from typing import Any

from aso_math.config import (
    ASO_SEQUENCE,
    LENGTH_SCAN_MAX,
    LENGTH_SCAN_MIN,
    SL_LENGTH,
    SL_SEQUENCE,
)
from aso_math.envelope import Timer, create_envelope, write_result
from aso_math.target_config import TargetConfig
from aso_math.thermo import (
    compute_dg,
    compute_hairpin_dg,
    compute_self_dimer_dg,
    compute_tm,
    gc_content,
    reverse_complement,
)
from core.logger import get_logger

logger = get_logger("04_exhaustive_optimization")

# ---------------------------------------------------------------------------
# Constantes de design LNA gapmer
# ---------------------------------------------------------------------------

# Faixa de LNA nos flancos (2-5 posicoes em cada extremidade)
LNA_FLANK_MIN: int = 2
LNA_FLANK_MAX: int = 5

# Gap minimo de DNA para ativacao de RNase H
# Ref: Crooke ST et al. (2017) Nat Rev Drug Discov — minimo 8-10 nt
RNASE_H_MIN_GAP: int = 8

# Boost de Tm por posicao LNA
# Ref: McTigue PM et al. (2004) Biochemistry 43(18):5388-5405
# Media observada: +2-4 C por LNA, usamos +3 C como media conservadora
TM_BOOST_PER_LNA: float = 3.0

# Limiar de LNA total acima do qual penalizamos (custo de sintese + toxicidade)
LNA_TOXICITY_THRESHOLD: int = 10

# ---------------------------------------------------------------------------
# Pesos da funcao de score unificada
# ---------------------------------------------------------------------------

W_DG: float = 0.35       # Energia livre de ligacao (quanto mais negativa, melhor)
W_TM: float = 0.25       # Temperatura de melting ajustada por LNA
W_GC: float = 0.15       # Conteudo GC proximo do ideal (~45%)
W_STRUCT: float = 0.15   # Ausencia de estruturas secundarias problematicas
W_LNA: float = 0.10      # Qualidade do gap para RNase H

# Limiares de estrutura secundaria
HAIRPIN_DG_THRESHOLD: float = -3.0      # kcal/mol — hairpin problematico
SELF_DIMER_DG_THRESHOLD: float = -5.0   # kcal/mol — self-dimer problematico

# Configuracao LNA conhecida de MRL-ASO-001
MRL_LNA_5PRIME: int = 4
MRL_LNA_3PRIME: int = 4


# ---------------------------------------------------------------------------
# 1. Enumeracao do espaco de sequencias
# ---------------------------------------------------------------------------


def enumerate_sequence_space(
    sl_sequence: str = SL_SEQUENCE,
    sl_length: int = SL_LENGTH,
    length_scan_min: int = LENGTH_SCAN_MIN,
    length_scan_max: int = LENGTH_SCAN_MAX,
) -> list[dict[str, Any]]:
    """Enumera todas as sequencias ASO validas sobre o SL RNA.

    Para cada comprimento L em [length_scan_min, length_scan_max] e cada
    posicao de inicio s em [0, sl_length-L], gera o ASO como
    reverse_complement(SL_RNA[s:s+L]) e computa todas as propriedades
    termodinamicas.

    Args:
        sl_sequence: Sequencia do SL RNA alvo. Default: L. infantum.
        sl_length: Comprimento do SL RNA. Default: 39.
        length_scan_min: Comprimento minimo do ASO a testar. Default: 18.
        length_scan_max: Comprimento maximo do ASO a testar. Default: 27.

    Returns:
        Lista de candidatos ordenada por dG (mais negativo primeiro).
    """
    candidates: list[dict[str, Any]] = []

    for length in range(length_scan_min, length_scan_max + 1):
        max_start = sl_length - length
        for start in range(max_start + 1):
            target_window = sl_sequence[start:start + length]
            aso_seq = reverse_complement(target_window)

            dg = compute_dg(aso_seq)
            tm = compute_tm(aso_seq)
            gc = gc_content(aso_seq)
            hairpin_dg = compute_hairpin_dg(aso_seq)
            self_dimer_dg = compute_self_dimer_dg(aso_seq)

            candidates.append({
                "sequence": aso_seq,
                "length": length,
                "target_start": start,
                "target_end": start + length,
                "dg_binding_kcal": dg,
                "tm_celsius": tm,
                "gc_content": round(gc, 4),
                "hairpin_dg_kcal": hairpin_dg,
                "self_dimer_dg_kcal": self_dimer_dg,
            })

    # Ordenar por dG (mais negativo = melhor)
    candidates.sort(key=lambda c: c["dg_binding_kcal"])

    logger.info(
        "Espaco de sequencias: %d candidatos enumerados, dG range [%.2f, %.2f]",
        len(candidates),
        candidates[0]["dg_binding_kcal"] if candidates else 0.0,
        candidates[-1]["dg_binding_kcal"] if candidates else 0.0,
    )

    return candidates


# ---------------------------------------------------------------------------
# 2. Enumeracao de configuracoes LNA gapmer
# ---------------------------------------------------------------------------


def enumerate_lna_configurations(aso_length: int) -> list[dict[str, Any]]:
    """Enumera todas as configuracoes LNA gapmer validas para um dado comprimento.

    Para cada combinacao de (lna_5prime, lna_3prime) em [2, 5] x [2, 5],
    verifica se o gap de DNA resultante >= 8 nt (requisito RNase H).

    Para cada configuracao valida, calcula:
        - gap_size: comprimento do segmento de DNA central
        - tm_boost: aumento estimado de Tm por LNA (+3 C por posicao)
        - total_lna: numero total de nucleotideos LNA
        - nuclease_resistance: fator proporcional ao total de LNA
        - toxicity_penalty: penalidade se total_lna > 10

    Args:
        aso_length: Comprimento do ASO candidato.

    Returns:
        Lista de configuracoes validas, ordenadas por gap_size decrescente.
    """
    configs: list[dict[str, Any]] = []

    for lna_5 in range(LNA_FLANK_MIN, LNA_FLANK_MAX + 1):
        for lna_3 in range(LNA_FLANK_MIN, LNA_FLANK_MAX + 1):
            gap_size = aso_length - lna_5 - lna_3

            # Restricao biologica: gap minimo para RNase H
            if gap_size < RNASE_H_MIN_GAP:
                continue

            total_lna = lna_5 + lna_3
            tm_boost = total_lna * TM_BOOST_PER_LNA

            # Resistencia a nucleases: proporcional ao total de LNA
            # Normalizada em [0, 1] com base no maximo possivel (10 LNA)
            nuclease_resistance = min(1.0, total_lna / 10.0)

            # Penalidade de toxicidade: se total_lna > 10, penalizamos
            # proporcionalmente ao excesso
            toxicity_penalty = max(0.0, (total_lna - LNA_TOXICITY_THRESHOLD) * 0.1)

            configs.append({
                "lna_5prime": lna_5,
                "lna_3prime": lna_3,
                "gap_size": gap_size,
                "total_lna": total_lna,
                "tm_boost_celsius": tm_boost,
                "nuclease_resistance": round(nuclease_resistance, 4),
                "toxicity_penalty": round(toxicity_penalty, 4),
            })

    # Ordenar por gap_size decrescente (maior gap = melhor para RNase H)
    configs.sort(key=lambda c: -c["gap_size"])

    return configs


# ---------------------------------------------------------------------------
# 3. Funcao de score unificada
# ---------------------------------------------------------------------------


def compute_design_score(
    candidate: dict[str, Any],
    lna_config: dict[str, Any],
) -> float:
    """Calcula score composto para uma combinacao (sequencia, config LNA).

    O score integra 5 dimensoes com pesos pre-definidos:
        score = w_dg * norm_dg + w_tm * norm_tm + w_gc * gc_score
              + w_struct * struct_score + w_lna * lna_score

    Cada componente e normalizado em [0, 1] para comparabilidade.
    Score mais alto = melhor candidato.

    Componentes:
        norm_dg:     abs(dG) / 50.0, capped em 1.0. Mais negativo = melhor.
        norm_tm:     Funcao triangular centrada em 75 C (faixa ideal 60-90 C).
                     Tm muito alta causa ligacao off-target e toxicidade;
                     Tm muito baixa causa ligacao fraca in vivo.
                     Ref: Crooke ST (2007) — ASO Tm otimo entre 65-85 C.
        gc_score:    1.0 - 2 * |GC - 0.45|. GC ideal ~45%.
        struct_score: 1.0 se hairpin e self-dimer sao fracos, penalidade gradual.
        lna_score:   Combina qualidade do gap para RNase H com penalidade de
                     toxicidade e bonus por gap otimo (10-15 nt).

    Args:
        candidate: Dicionario com propriedades da sequencia candidata.
        lna_config: Dicionario com configuracao LNA gapmer.

    Returns:
        Score composto em [0, 1], arredondado a 4 casas.
    """
    # --- Componente 1: Energia livre de ligacao ---
    dg = candidate["dg_binding_kcal"]
    norm_dg = min(1.0, abs(dg) / 50.0)

    # --- Componente 2: Tm ajustada por LNA ---
    # Funcao triangular: pico em 75 C, zero abaixo de 55 C e acima de 95 C.
    # Penaliza Tm excessiva (off-target binding, toxicidade hepatica).
    # Ref: Monia BP et al. (1993) JBC — janela termica ideal para ASOs
    adjusted_tm = candidate["tm_celsius"] + lna_config["tm_boost_celsius"]
    tm_optimal = 75.0
    tm_half_width = 20.0  # Faixa de 55-95 C
    norm_tm = max(0.0, 1.0 - abs(adjusted_tm - tm_optimal) / tm_half_width)

    # --- Componente 3: Conteudo GC ---
    gc = candidate["gc_content"]
    gc_score = max(0.0, 1.0 - 2.0 * abs(gc - 0.45))

    # --- Componente 4: Estrutura secundaria ---
    hairpin_dg = candidate["hairpin_dg_kcal"]
    self_dimer_dg = candidate["self_dimer_dg_kcal"]

    # Score de 1.0 se ambos acima dos limiares, penalidade gradual caso contrario
    hairpin_ok = hairpin_dg > HAIRPIN_DG_THRESHOLD
    dimer_ok = self_dimer_dg > SELF_DIMER_DG_THRESHOLD

    if hairpin_ok and dimer_ok:
        struct_score = 1.0
    elif hairpin_ok or dimer_ok:
        # Um problematico: penalidade parcial proporcional a severidade
        penalty = 0.0
        if not hairpin_ok:
            penalty += min(0.5, abs(hairpin_dg - HAIRPIN_DG_THRESHOLD) * 0.1)
        if not dimer_ok:
            penalty += min(0.5, abs(self_dimer_dg - SELF_DIMER_DG_THRESHOLD) * 0.1)
        struct_score = max(0.0, 1.0 - penalty)
    else:
        # Ambos problematicos: penalidade severa
        penalty_h = min(0.5, abs(hairpin_dg - HAIRPIN_DG_THRESHOLD) * 0.1)
        penalty_d = min(0.5, abs(self_dimer_dg - SELF_DIMER_DG_THRESHOLD) * 0.1)
        struct_score = max(0.0, 1.0 - penalty_h - penalty_d)

    # --- Componente 5: Qualidade do design LNA para RNase H ---
    gap_size = lna_config["gap_size"]
    # Gap otimo: 10-15 nt. Gaps menores reduzem eficiencia de RNase H,
    # gaps muito grandes reduzem a protecao dos flancos LNA.
    # Ref: Kurreck J (2003) Eur J Biochem — gap 10-14 nt ideal para gapmers
    if gap_size <= 15:
        lna_score = max(0.0, min(1.0, gap_size / 10.0))
    else:
        # Penalidade suave para gaps excessivos (> 15 nt)
        # DNA desprotegido no meio e vulneravel a nucleases
        lna_score = max(0.0, 1.0 - (gap_size - 15) * 0.05)

    # Aplicar penalidade de toxicidade por excesso de LNA
    lna_score = max(0.0, lna_score - lna_config["toxicity_penalty"])

    # --- Score composto ---
    score = (
        W_DG * norm_dg
        + W_TM * norm_tm
        + W_GC * gc_score
        + W_STRUCT * struct_score
        + W_LNA * lna_score
    )

    return round(score, 4)


# ---------------------------------------------------------------------------
# 4. Enumeracao exaustiva completa
# ---------------------------------------------------------------------------


def run_full_enumeration(
    candidates: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Enumera todas as combinacoes (sequencia x configuracao LNA) e rankeia.

    Para cada candidato de sequencia, gera todas as configuracoes LNA validas
    para seu comprimento, calcula o score composto e rankeia globalmente.

    Args:
        candidates: Lista de candidatos de sequencia (output de enumerate_sequence_space).

    Returns:
        Lista de todas as combinacoes, ordenada por score decrescente.
    """
    # Cache de configs LNA por comprimento (evita recalculo)
    lna_cache: dict[int, list[dict[str, Any]]] = {}
    all_designs: list[dict[str, Any]] = []

    for candidate in candidates:
        length = candidate["length"]

        # Buscar configs LNA do cache ou gerar
        if length not in lna_cache:
            lna_cache[length] = enumerate_lna_configurations(length)

        configs = lna_cache[length]

        for lna_config in configs:
            score = compute_design_score(candidate, lna_config)

            adjusted_tm = candidate["tm_celsius"] + lna_config["tm_boost_celsius"]

            all_designs.append({
                "sequence": candidate["sequence"],
                "length": candidate["length"],
                "target_start": candidate["target_start"],
                "target_end": candidate["target_end"],
                "dg_binding_kcal": candidate["dg_binding_kcal"],
                "tm_base_celsius": candidate["tm_celsius"],
                "tm_adjusted_celsius": round(adjusted_tm, 2),
                "gc_content": candidate["gc_content"],
                "hairpin_dg_kcal": candidate["hairpin_dg_kcal"],
                "self_dimer_dg_kcal": candidate["self_dimer_dg_kcal"],
                "lna_5prime": lna_config["lna_5prime"],
                "lna_3prime": lna_config["lna_3prime"],
                "gap_size": lna_config["gap_size"],
                "total_lna": lna_config["total_lna"],
                "nuclease_resistance": lna_config["nuclease_resistance"],
                "score": score,
            })

    # Ordenar por score decrescente (melhor primeiro)
    all_designs.sort(key=lambda d: -d["score"])

    # Atribuir rank
    for rank, design in enumerate(all_designs, start=1):
        design["rank"] = rank

    logger.info(
        "Enumeracao completa: %d combinacoes avaliadas, "
        "score range [%.4f, %.4f]",
        len(all_designs),
        all_designs[-1]["score"] if all_designs else 0.0,
        all_designs[0]["score"] if all_designs else 0.0,
    )

    return all_designs


# ---------------------------------------------------------------------------
# 5. Comparacao com MRL-ASO-001
# ---------------------------------------------------------------------------


def compare_with_mrl_aso_001(
    ranked_results: list[dict[str, Any]],
    aso_sequence: str = ASO_SEQUENCE,
    lna_5prime: int = MRL_LNA_5PRIME,
    lna_3prime: int = MRL_LNA_3PRIME,
) -> dict[str, Any]:
    """Localiza o ASO de referencia nos resultados rankeados e reporta status.

    Busca pela sequencia exata do ASO candidato com a configuracao LNA
    especificada. Reporta honestamente o rank, score e margem.

    Args:
        ranked_results: Lista ordenada por score decrescente.
        aso_sequence: Sequencia do ASO a localizar. Default: MRL-ASO-001.
        lna_5prime: Posicoes LNA no flanco 5'. Default: 4.
        lna_3prime: Posicoes LNA no flanco 3'. Default: 4.

    Returns:
        Dicionario com rank, score, margem e conclusao.
    """
    # Calcular gap esperado para o log de conclusao
    expected_gap = len(aso_sequence) - lna_5prime - lna_3prime

    # Buscar ASO com config LNA exata
    mrl_entry: dict[str, Any] | None = None
    mrl_rank: int = -1

    for design in ranked_results:
        if (
            design["sequence"] == aso_sequence
            and design["lna_5prime"] == lna_5prime
            and design["lna_3prime"] == lna_3prime
        ):
            mrl_entry = design
            mrl_rank = design["rank"]
            break

    if mrl_entry is None:
        logger.warning(
            "ASO referencia (seq=%s, LNA %d+%d) NAO encontrado nos resultados!",
            aso_sequence,
            lna_5prime,
            lna_3prime,
        )
        return {
            "found": False,
            "mrl_aso_001_rank": -1,
            "mrl_aso_001_score": 0.0,
            "is_optimal": False,
            "conclusion": "ASO de referencia nao encontrado no espaco de busca.",
        }

    # Score do primeiro e segundo colocados
    top_score = ranked_results[0]["score"]
    second_score = ranked_results[1]["score"] if len(ranked_results) > 1 else 0.0
    margin_to_second = round(top_score - second_score, 4)

    mrl_score = mrl_entry["score"]
    is_optimal = mrl_rank == 1

    if is_optimal:
        # Calcular margem do ASO sobre o segundo
        conclusion = (
            f"ASO referencia (LNA {lna_5prime}+{lna_3prime}, gap {expected_gap}) "
            f"CONFIRMADO como design otimo. "
            f"Score: {mrl_score:.4f}. "
            f"Margem sobre o segundo colocado: {margin_to_second:.4f} "
            f"({margin_to_second / top_score * 100:.1f}% relativa)."
        )
    else:
        # Reportar honestamente quem e melhor
        winner = ranked_results[0]
        score_diff = round(winner["score"] - mrl_score, 4)
        pct_better = round(score_diff / mrl_score * 100, 1) if mrl_score > 0 else 0.0

        conclusion = (
            f"ASO referencia (LNA {lna_5prime}+{lna_3prime}, gap {expected_gap}) "
            f"esta no rank #{mrl_rank} "
            f"com score {mrl_score:.4f}. "
            f"O design otimo e: seq={winner['sequence']}, "
            f"LNA {winner['lna_5prime']}+{winner['lna_3prime']}, "
            f"gap {winner['gap_size']}, score {winner['score']:.4f} "
            f"({pct_better}% melhor)."
        )

    logger.info("Comparacao ASO referencia: rank #%d, score %.4f", mrl_rank, mrl_score)

    return {
        "found": True,
        "mrl_aso_001_rank": mrl_rank,
        "mrl_aso_001_score": mrl_score,
        "mrl_aso_001_entry": mrl_entry,
        "is_optimal": is_optimal,
        "margin_to_second": margin_to_second,
        "conclusion": conclusion,
    }


# ---------------------------------------------------------------------------
# 6. Analise de frente de Pareto
# ---------------------------------------------------------------------------


def compute_pareto_front(
    ranked_results: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Identifica candidatos Pareto-otimos em 4 dimensoes com trade-offs reais.

    Um candidato e Pareto-otimo se nenhum outro candidato e estritamente
    melhor em TODAS as dimensoes simultaneamente. A escolha de objetivos
    garante trade-offs genuinos (nao e possivel maximizar todos ao mesmo tempo).

    Dimensoes (todas maximizadas):
        1. Afinidade de ligacao: abs(dG) — favorece sequencias mais longas e GC-ricas
        2. Proximidade ao Tm ideal: 1 - |Tm_adj - 75|/20 — penaliza extremos
        3. Qualidade estrutural: ausencia de hairpin/self-dimer
        4. Eficiencia de RNase H: gap no range 10-15 nt (gap muito grande perde
           protecao dos flancos LNA, muito pequeno reduz clivagem)

    Trade-offs biologicos capturados:
        - dG vs Tm: mais LNA aumenta Tm mas pode ultrapassar o ideal
        - dG vs gap: sequencias longas com pouco LNA tem gap grande mas baixa afinidade LNA
        - Tm vs gap: mais LNA melhora Tm mas reduz gap

    Args:
        ranked_results: Lista de todos os designs avaliados.

    Returns:
        Lista de candidatos Pareto-otimos, ordenada por score composto.
    """
    # Extrair vetores de objetivos (todos maximizados)
    designs_with_objectives: list[tuple[dict[str, Any], list[float]]] = []

    for design in ranked_results:
        # Objetivo 1: Afinidade de ligacao
        obj_dg = abs(design["dg_binding_kcal"])

        # Objetivo 2: Proximidade ao Tm ideal (75 C)
        adjusted_tm = design["tm_adjusted_celsius"]
        obj_tm = max(0.0, 1.0 - abs(adjusted_tm - 75.0) / 20.0)

        # Objetivo 3: Qualidade estrutural
        hairpin_ok = design["hairpin_dg_kcal"] > HAIRPIN_DG_THRESHOLD
        dimer_ok = design["self_dimer_dg_kcal"] > SELF_DIMER_DG_THRESHOLD
        if hairpin_ok and dimer_ok:
            obj_struct = 1.0
        elif hairpin_ok or dimer_ok:
            obj_struct = 0.5
        else:
            obj_struct = 0.0

        # Objetivo 4: Eficiencia do gap para RNase H (pico em 10-15 nt)
        gap = design["gap_size"]
        if gap <= 15:
            obj_gap = min(1.0, gap / 10.0)
        else:
            obj_gap = max(0.0, 1.0 - (gap - 15) * 0.05)

        designs_with_objectives.append(
            (design, [obj_dg, obj_tm, obj_struct, obj_gap])
        )

    # Filtro de Pareto: candidato i e dominado se existe j tal que
    # j >= i em todas as dimensoes e j > i em pelo menos uma
    num_objectives = 4
    pareto_front: list[dict[str, Any]] = []

    n = len(designs_with_objectives)
    for i in range(n):
        design_i, obj_i = designs_with_objectives[i]
        is_dominated = False

        for j in range(n):
            if i == j:
                continue
            _, obj_j = designs_with_objectives[j]

            # Checar se j domina i
            all_geq = all(obj_j[k] >= obj_i[k] for k in range(num_objectives))
            any_greater = any(obj_j[k] > obj_i[k] for k in range(num_objectives))

            if all_geq and any_greater:
                is_dominated = True
                break

        if not is_dominated:
            pareto_front.append(design_i)

    # Ordenar por score composto
    pareto_front.sort(key=lambda d: -d["score"])

    logger.info(
        "Frente de Pareto: %d candidatos nao-dominados de %d total",
        len(pareto_front),
        n,
    )

    return pareto_front


# ---------------------------------------------------------------------------
# 7. Funcoes auxiliares para construcao do output
# ---------------------------------------------------------------------------


def _build_length_distribution(candidates: list[dict[str, Any]]) -> dict[str, int]:
    """Conta candidatos por comprimento para o relatorio."""
    dist: dict[str, int] = {}
    for c in candidates:
        key = str(c["length"])
        dist[key] = dist.get(key, 0) + 1
    return dist


def _build_best_per_length(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Encontra o melhor candidato (por dG) para cada comprimento."""
    best: dict[int, dict[str, Any]] = {}
    for c in candidates:
        length = c["length"]
        if length not in best or c["dg_binding_kcal"] < best[length]["dg_binding_kcal"]:
            best[length] = c

    result = list(best.values())
    result.sort(key=lambda c: c["length"])
    return result


# ---------------------------------------------------------------------------
# 8. Orquestrador principal
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a otimizacao exaustiva completa do espaco de design ASO.

    Fluxo:
        1. Enumerar todas as sequencias candidatas sobre o SL RNA
        2. Para cada sequencia, enumerar todas as configuracoes LNA validas
        3. Calcular score composto para cada combinacao
        4. Rankear globalmente
        5. Localizar ASO de referencia e reportar status de otimalidade
        6. Calcular frente de Pareto multi-objetivo
        7. Gravar resultados em JSON

    Args:
        config: Configuracao parametrizavel do organismo-alvo.
                Se None, usa defaults de L. infantum (retrocompativel).

    Retorna o envelope completo para uso programatico.
    """
    # Retrocompatibilidade: se nenhum config fornecido, usar defaults de L. infantum
    if config is None:
        config = TargetConfig()

    # Extrair parametros do config para uso local
    sl_sequence = config.sl_sequence
    sl_length = config.sl_length
    aso_sequence = config.aso_sequence if config.aso_sequence else ASO_SEQUENCE
    length_scan_min = config.length_scan_min
    length_scan_max = config.length_scan_max
    lna_5prime = config.lna_5prime
    lna_3prime = config.lna_3prime

    logger.info("=" * 60)
    logger.info("MODULO 04: Otimizacao Exaustiva do Espaco de Design ASO")
    logger.info("=" * 60)

    envelope = create_envelope("04_exhaustive_optimization")

    with Timer() as timer:
        # --- Passo 1: Enumerar espaco de sequencias ---
        logger.info(
            "Passo 1/5: Enumerando sequencias candidatas "
            "(comprimento %d-%d nt sobre SL RNA de %d nt)...",
            length_scan_min, length_scan_max, sl_length,
        )
        candidates = enumerate_sequence_space(
            sl_sequence=sl_sequence,
            sl_length=sl_length,
            length_scan_min=length_scan_min,
            length_scan_max=length_scan_max,
        )

        length_distribution = _build_length_distribution(candidates)
        best_per_length = _build_best_per_length(candidates)

        logger.info(
            "  Total de sequencias: %d, distribuicao por comprimento: %s",
            len(candidates),
            length_distribution,
        )

        # --- Passo 2: Verificar configs LNA por comprimento ---
        logger.info("Passo 2/5: Enumerando configuracoes LNA gapmer por comprimento...")

        total_lna_configs = 0
        lna_summary: dict[int, int] = {}
        for length in range(length_scan_min, length_scan_max + 1):
            configs = enumerate_lna_configurations(length)
            lna_summary[length] = len(configs)
            # Contar candidatos com este comprimento
            count_at_length = sum(1 for c in candidates if c["length"] == length)
            total_lna_configs += len(configs) * count_at_length

        logger.info(
            "  Configs LNA por comprimento: %s",
            {str(k): v for k, v in sorted(lna_summary.items())},
        )
        logger.info(
            "  Total de combinacoes (sequencia x LNA): %d",
            total_lna_configs,
        )

        # --- Passo 3: Enumeracao exaustiva ---
        logger.info(
            "Passo 3/5: Avaliando %d sequencias x configs LNA...",
            len(candidates),
        )
        ranked_results = run_full_enumeration(candidates)

        logger.info(
            "  Top 3: %s",
            [
                f"score={d['score']:.4f}, seq={d['sequence'][:12]}..., "
                f"LNA {d['lna_5prime']}+{d['lna_3prime']}"
                for d in ranked_results[:3]
            ],
        )

        # --- Passo 4: Comparar com ASO de referencia ---
        logger.info("Passo 4/5: Localizando ASO de referencia no ranking...")
        comparison = compare_with_mrl_aso_001(
            ranked_results,
            aso_sequence=aso_sequence,
            lna_5prime=lna_5prime,
            lna_3prime=lna_3prime,
        )

        logger.info("  %s", comparison["conclusion"])

        # --- Passo 5: Frente de Pareto ---
        logger.info("Passo 5/5: Calculando frente de Pareto (dG x Tm x estrutura)...")
        pareto_front = compute_pareto_front(ranked_results)

        # Verificar se ASO de referencia esta na frente de Pareto
        mrl_in_pareto = any(
            d["sequence"] == aso_sequence
            and d["lna_5prime"] == lna_5prime
            and d["lna_3prime"] == lna_3prime
            for d in pareto_front
        )

        logger.info(
            "  Pareto: %d candidatos nao-dominados, ASO referencia na frente: %s",
            len(pareto_front),
            "SIM" if mrl_in_pareto else "NAO",
        )

        # --- Construir conclusao ---
        is_optimal = comparison.get("is_optimal", False)
        mrl_rank = comparison.get("mrl_aso_001_rank", -1)
        mrl_score = comparison.get("mrl_aso_001_score", 0.0)

        if is_optimal:
            conclusion = (
                f"Enumeracao exaustiva de {len(ranked_results)} combinacoes "
                f"(sequencia x LNA gapmer) CONFIRMA ASO referencia como design otimo. "
                f"Score: {mrl_score:.4f}, margem sobre #2: "
                f"{comparison.get('margin_to_second', 0):.4f}."
            )
        else:
            winner = ranked_results[0]
            conclusion = (
                f"Enumeracao exaustiva de {len(ranked_results)} combinacoes "
                f"posiciona ASO referencia no rank #{mrl_rank} (score {mrl_score:.4f}). "
                f"Design otimo: seq={winner['sequence']}, "
                f"LNA {winner['lna_5prime']}+{winner['lna_3prime']}, "
                f"score {winner['score']:.4f}."
            )

        logger.info("CONCLUSAO: %s", conclusion)

    # --- Montar envelope ---
    envelope["runtime_seconds"] = timer.elapsed
    envelope["status"] = "success"
    envelope["summary"]["conclusion"] = conclusion
    envelope["summary"]["key_metrics"] = {
        "total_sequences": len(candidates),
        "total_combinations": len(ranked_results),
        "mrl_aso_001_rank": mrl_rank,
        "mrl_aso_001_score": mrl_score,
        "is_optimal": is_optimal,
        "pareto_front_size": len(pareto_front),
        "mrl_in_pareto": mrl_in_pareto,
    }

    # Top 10 para o relatorio (sem incluir todos os ~1700)
    top_10 = []
    for d in ranked_results[:10]:
        top_10.append({
            "rank": d["rank"],
            "sequence": d["sequence"],
            "length": d["length"],
            "target_start": d["target_start"],
            "target_end": d["target_end"],
            "dg_binding_kcal": d["dg_binding_kcal"],
            "tm_base_celsius": d["tm_base_celsius"],
            "tm_adjusted_celsius": d["tm_adjusted_celsius"],
            "gc_content": d["gc_content"],
            "lna_5prime": d["lna_5prime"],
            "lna_3prime": d["lna_3prime"],
            "gap_size": d["gap_size"],
            "score": d["score"],
        })

    # Pareto front resumido
    pareto_summary = []
    for d in pareto_front[:20]:  # Limitar a 20 para o JSON
        pareto_summary.append({
            "rank": d["rank"],
            "sequence": d["sequence"],
            "length": d["length"],
            "dg_binding_kcal": d["dg_binding_kcal"],
            "tm_adjusted_celsius": d["tm_adjusted_celsius"],
            "gc_content": d["gc_content"],
            "lna_5prime": d["lna_5prime"],
            "lna_3prime": d["lna_3prime"],
            "gap_size": d["gap_size"],
            "score": d["score"],
        })

    envelope["data"] = {
        "sequence_space": {
            "total_candidates": len(candidates),
            "length_distribution": length_distribution,
            "best_per_length": [
                {
                    "length": c["length"],
                    "sequence": c["sequence"],
                    "dg_binding_kcal": c["dg_binding_kcal"],
                    "tm_celsius": c["tm_celsius"],
                    "gc_content": c["gc_content"],
                }
                for c in best_per_length
            ],
        },
        "lna_configurations": {
            "total_valid": len(ranked_results),
            "lna_range_5prime": [LNA_FLANK_MIN, LNA_FLANK_MAX],
            "lna_range_3prime": [LNA_FLANK_MIN, LNA_FLANK_MAX],
            "min_gap_size": RNASE_H_MIN_GAP,
            "configs_per_length": {
                str(k): v for k, v in sorted(lna_summary.items())
            },
        },
        "full_enumeration": {
            "total_evaluated": len(ranked_results),
            "top_10": top_10,
            "mrl_aso_001_rank": mrl_rank,
            "mrl_aso_001_score": mrl_score,
            "margin_to_second": comparison.get("margin_to_second", 0.0),
            "is_optimal": is_optimal,
            "scoring_weights": {
                "w_dg": W_DG,
                "w_tm": W_TM,
                "w_gc": W_GC,
                "w_struct": W_STRUCT,
                "w_lna": W_LNA,
            },
        },
        "pareto_front": {
            "description": (
                "Sequencias Pareto-otimas considerando dG, Tm ajustada "
                "e qualidade estrutural. Nenhum outro candidato domina "
                "simultaneamente em todas as dimensoes."
            ),
            "total_non_dominated": len(pareto_front),
            "mrl_in_pareto": mrl_in_pareto,
            "candidates": pareto_summary,
        },
        "mrl_aso_001_comparison": comparison,
    }

    # Gravar resultado
    output_path = write_result(envelope)
    logger.info("Resultado gravado em: %s", output_path)
    logger.info("Tempo de execucao: %.2f segundos", timer.elapsed)

    return envelope


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
