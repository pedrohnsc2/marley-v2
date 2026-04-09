"""Modulo 01 — Mapeamento completo da paisagem termodinamica ao redor de MRL-ASO-001.

Objetivo: provar que MRL-ASO-001 se encontra no minimo global (ou local profundo)
de energia livre de ligacao dentro da vizinhanca de Hamming-1, Hamming-2, e
sobre todas as janelas de comprimento variavel no SL RNA alvo.

A logica matematica e simples:
    Se dG(wildtype) <= dG(mutante) para TODOS os mutantes na vizinhanca,
    entao o wildtype e um otimo local. Se alem disso ele tambem vence todas
    as janelas de comprimento variavel, e um otimo global sobre o espaco de
    busca factivel.

O fitness score incorpora penalidades por estruturas secundarias indesejaveis
(hairpin e auto-dimerizacao), pois um ASO com dG muito negativo mas que forma
hairpin nao funciona in vivo.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any

from aso_math.config import (
    ASO_KNOWN_DG,
    ASO_KNOWN_TM,
    ASO_LENGTH,
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    BASES,
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

logger = get_logger("01_thermodynamic_landscape")

# ---------------------------------------------------------------------------
# Constantes do fitness score
# Penalizamos hairpin apenas se dG_hairpin < -2.0 kcal/mol (estrutura estavel)
# e self-dimer apenas se dG_self_dimer < -3.0 kcal/mol.
# O peso 0.3 foi calibrado para que penalidades severas (~5 kcal/mol)
# desqualifiquem o candidato sem dominar candidatos sem problemas.
# ---------------------------------------------------------------------------

HAIRPIN_THRESHOLD: float = -2.0
SELF_DIMER_THRESHOLD: float = -3.0
PENALTY_WEIGHT: float = 0.3


# ---------------------------------------------------------------------------
# Funcoes auxiliares
# ---------------------------------------------------------------------------


def _compute_fitness(dg_binding: float, hairpin_dg: float, self_dimer_dg: float) -> float:
    """Calcula o fitness score composto.

    fitness = dG_binding - 0.3 * max(0, hairpin_penalty) - 0.3 * max(0, dimer_penalty)

    Onde hairpin_penalty ativa apenas se o hairpin for mais estavel que o limiar,
    e analogamente para self-dimer. Valor mais negativo = melhor.

    A intuicao: queremos o dG mais negativo possivel, mas descontamos
    candidatos que formam estruturas secundarias problematicas.
    """
    hairpin_penalty = max(0.0, hairpin_dg - HAIRPIN_THRESHOLD)
    dimer_penalty = max(0.0, self_dimer_dg - SELF_DIMER_THRESHOLD)
    return round(dg_binding - PENALTY_WEIGHT * hairpin_penalty - PENALTY_WEIGHT * dimer_penalty, 4)


def _characterize_aso(sequence: str) -> dict[str, Any]:
    """Calcula todas as propriedades termodinamicas de um ASO candidato.

    Retorna dicionario com dG, Tm, hairpin, self-dimer, GC e fitness.
    Funcao pura — sem efeitos colaterais.
    """
    dg = compute_dg(sequence)
    tm = compute_tm(sequence)
    hairpin = compute_hairpin_dg(sequence)
    self_dimer = compute_self_dimer_dg(sequence)
    gc = gc_content(sequence)
    fitness = _compute_fitness(dg, hairpin, self_dimer)

    return {
        "sequence": sequence,
        "length": len(sequence),
        "dg_binding_kcal": dg,
        "tm_celsius": tm,
        "dg_hairpin_kcal": hairpin,
        "dg_self_dimer_kcal": self_dimer,
        "gc_content": round(gc, 4),
        "fitness_score": fitness,
    }


# ---------------------------------------------------------------------------
# 1. Mutantes pontuais (Hamming-1)
# Para cada posicao i em [0, 24], testamos as 3 bases alternativas.
# Total: 25 * 3 = 75 mutantes.
# Se nenhum tiver dG mais negativo que o wildtype, provamos otimalidade local.
# ---------------------------------------------------------------------------


def generate_point_mutants(aso_sequence: str | None = None) -> list[dict[str, Any]]:
    """Gera e avalia todos os mutantes pontuais de um ASO candidato.

    Para cada posicao, substitui pela 3 bases alternativas e calcula
    o perfil termodinamico completo. Retorna lista ordenada por dG.

    Args:
        aso_sequence: Sequencia do ASO. Se None, usa ASO_SEQUENCE do config.
    """
    # Sequencia do wildtype — usa parametro ou fallback para constante global
    wt = (aso_sequence or ASO_SEQUENCE).upper()
    aso_length = len(wt)
    mutants: list[dict[str, Any]] = []

    for pos in range(aso_length):
        original_base = wt[pos]
        for base in BASES:
            if base == original_base:
                continue

            # Construir mutante: substituir base na posicao pos
            mutant_seq = wt[:pos] + base + wt[pos + 1:]
            profile = _characterize_aso(mutant_seq)
            profile["position"] = pos
            profile["original_base"] = original_base
            profile["mutant_base"] = base
            profile["mutation"] = f"{original_base}{pos + 1}{base}"
            mutants.append(profile)

    # Ordenar por dG de ligacao (mais negativo = melhor = primeiro)
    mutants.sort(key=lambda m: m["dg_binding_kcal"])

    logger.info(
        "Mutantes pontuais: %d avaliados, melhor dG = %.2f, pior dG = %.2f",
        len(mutants),
        mutants[0]["dg_binding_kcal"],
        mutants[-1]["dg_binding_kcal"],
    )

    return mutants


# ---------------------------------------------------------------------------
# 2. Varredura de comprimento (length scan)
# Para cada comprimento L em [18, 27] e cada posicao de inicio valida
# no SL RNA, geramos o ASO complementar (reverse complement da janela alvo).
# Isso cobre o espaco completo de ASOs anti-SL de comprimento factivel.
# ---------------------------------------------------------------------------


def scan_length_variants(
    sl_sequence: str | None = None,
    sl_length: int | None = None,
    length_scan_min: int | None = None,
    length_scan_max: int | None = None,
) -> list[dict[str, Any]]:
    """Varre todas as janelas de comprimento variavel sobre o SL RNA.

    Para cada comprimento L em [length_scan_min, length_scan_max] e cada
    posicao de inicio s em [0, sl_length - L], gera o ASO complementar
    e calcula dG, Tm e GC.

    Args:
        sl_sequence: Sequencia do SL RNA. Se None, usa SL_SEQUENCE do config.
        sl_length: Comprimento do SL RNA. Se None, calcula de sl_sequence.
        length_scan_min: Comprimento minimo da varredura. Se None, usa LENGTH_SCAN_MIN.
        length_scan_max: Comprimento maximo da varredura. Se None, usa LENGTH_SCAN_MAX.
    """
    # Parametros — usa argumentos ou fallback para constantes globais
    _sl_seq = (sl_sequence or SL_SEQUENCE).upper()
    _sl_len = sl_length if sl_length is not None else len(_sl_seq)
    _len_min = length_scan_min if length_scan_min is not None else LENGTH_SCAN_MIN
    _len_max = length_scan_max if length_scan_max is not None else LENGTH_SCAN_MAX

    variants: list[dict[str, Any]] = []

    for length in range(_len_min, _len_max + 1):
        max_start = _sl_len - length
        for start in range(max_start + 1):
            target_window = _sl_seq[start:start + length]
            aso_candidate = reverse_complement(target_window)

            dg = compute_dg(aso_candidate)
            tm = compute_tm(aso_candidate)
            gc = gc_content(aso_candidate)

            variants.append({
                "sequence": aso_candidate,
                "length": length,
                "target_start": start,
                "target_end": start + length,
                "target_window": target_window,
                "dg_binding_kcal": dg,
                "tm_celsius": tm,
                "gc_content": round(gc, 4),
            })

    # Ordenar por dG (mais negativo primeiro)
    variants.sort(key=lambda v: v["dg_binding_kcal"])

    total = len(variants)
    expected = sum(_sl_len - l + 1 for l in range(_len_min, _len_max + 1))
    logger.info(
        "Varredura de comprimento: %d candidatos avaliados (esperado: %d), melhor dG = %.2f",
        total,
        expected,
        variants[0]["dg_binding_kcal"] if variants else 0.0,
    )

    return variants


# ---------------------------------------------------------------------------
# 3. Dados para heatmap (25 posicoes x 4 bases)
# Matriz de dG para visualizacao. Na posicao do wildtype, usamos o dG real;
# nas alternativas, usamos o dG do mutante correspondente.
# ---------------------------------------------------------------------------


def build_heatmap_data(
    wildtype_dg: float,
    point_mutants: list[dict[str, Any]],
    aso_sequence: str | None = None,
) -> dict[str, Any]:
    """Constroi matriz NxB de valores de dG para heatmap.

    Cada celula (i, j) contem o dG do ASO com base BASES[j] na posicao i.
    Para a base wildtype, usamos wildtype_dg. Para mutantes, buscamos
    na lista de mutantes pontuais.

    A visualizacao permite identificar rapidamente quais posicoes sao
    mais sensiveis a mutacao (colunas com grande variacao) e quais
    sao robustas (colunas uniformes).

    Args:
        wildtype_dg: dG do ASO wildtype.
        point_mutants: Lista de mutantes pontuais com dG.
        aso_sequence: Sequencia do ASO. Se None, usa ASO_SEQUENCE do config.
    """
    wt = (aso_sequence or ASO_SEQUENCE).upper()
    aso_length = len(wt)
    bases_list = list(BASES)  # ["A", "C", "G", "T"]

    # Indexar mutantes por (posicao, base) para lookup O(1)
    mutant_index: dict[tuple[int, str], float] = {}
    for m in point_mutants:
        mutant_index[(m["position"], m["mutant_base"])] = m["dg_binding_kcal"]

    positions = list(range(aso_length))
    dg_matrix: list[list[float]] = []

    for pos in positions:
        row: list[float] = []
        for base in bases_list:
            if base == wt[pos]:
                # Base wildtype — dG do ASO original
                row.append(wildtype_dg)
            else:
                # Base mutante — buscar no indice
                row.append(mutant_index.get((pos, base), 0.0))
        dg_matrix.append(row)

    return {
        "positions": positions,
        "bases": bases_list,
        "dg_matrix": dg_matrix,
    }


# ---------------------------------------------------------------------------
# 4. Mutantes duplos (Hamming-2)
# C(25, 2) = 300 pares de posicoes. Para cada par, 3 x 3 = 9 combinacoes
# de bases alternativas. Total: 300 * 9 = 2700 mutantes.
# Isso estende a prova de otimalidade de Hamming-1 para Hamming-2,
# refutando a critica de que mutacoes pontuais so provam otimalidade local.
# ---------------------------------------------------------------------------


def generate_double_mutants(aso_sequence: str | None = None) -> list[dict[str, Any]]:
    """Gera e avalia todos os mutantes duplos de um ASO candidato.

    Para cada par de posicoes (i, j) com i < j, testa todas as 9
    combinacoes de bases alternativas e calcula o perfil completo.

    A intuicao matematica: se o wildtype e um minimo local em Hamming-1
    mas nao em Hamming-2, existiria um "vale" acessivel apenas por
    mutacao dupla simultanea — biologicamente improvavel mas
    matematicamente possivel. Esta varredura descarta essa possibilidade.

    Args:
        aso_sequence: Sequencia do ASO. Se None, usa ASO_SEQUENCE do config.
    """
    wt = (aso_sequence or ASO_SEQUENCE).upper()
    aso_length = len(wt)
    mutants: list[dict[str, Any]] = []

    # Gerar bases alternativas para cada posicao (excluindo wildtype)
    alt_bases: list[list[str]] = []
    for pos in range(aso_length):
        alt_bases.append([b for b in BASES if b != wt[pos]])

    # Iterar sobre todos os pares C(N, 2)
    for pos1, pos2 in combinations(range(aso_length), 2):
        for base1 in alt_bases[pos1]:
            for base2 in alt_bases[pos2]:
                # Construir duplo mutante
                seq_list = list(wt)
                seq_list[pos1] = base1
                seq_list[pos2] = base2
                mutant_seq = "".join(seq_list)

                profile = _characterize_aso(mutant_seq)
                profile["positions"] = [pos1, pos2]
                profile["original_bases"] = [wt[pos1], wt[pos2]]
                profile["mutant_bases"] = [base1, base2]
                profile["mutation"] = (
                    f"{wt[pos1]}{pos1 + 1}{base1}+"
                    f"{wt[pos2]}{pos2 + 1}{base2}"
                )
                mutants.append(profile)

    # Ordenar por dG (mais negativo = melhor)
    mutants.sort(key=lambda m: m["dg_binding_kcal"])

    logger.info(
        "Mutantes duplos: %d avaliados, melhor dG = %.2f, pior dG = %.2f",
        len(mutants),
        mutants[0]["dg_binding_kcal"],
        mutants[-1]["dg_binding_kcal"],
    )

    return mutants


# ---------------------------------------------------------------------------
# 5. Orquestrador principal
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a analise completa da paisagem termodinamica.

    Fluxo:
        1. Caracterizar o wildtype ASO
        2. Gerar e avaliar mutantes pontuais (Hamming-1)
        3. Gerar e avaliar mutantes duplos (Hamming-2)
        4. Varrer todas as janelas de comprimento sobre o SL RNA
        5. Construir dados de heatmap para visualizacao
        6. Determinar se o wildtype e o otimo global
        7. Gravar resultados em JSON

    Args:
        config: Configuracao do organismo alvo. Se None, usa defaults de
                L. infantum (retrocompativel com o pipeline original).

    Retorna o envelope completo para uso programatico.
    """
    # --- Construir config padrao se nenhum foi fornecido (retrocompatibilidade) ---
    if config is None:
        config = TargetConfig(
            aso_sequence=ASO_SEQUENCE,
            sl_sequence=SL_SEQUENCE,
            aso_target_start=ASO_TARGET_START,
            aso_target_end=ASO_TARGET_END,
            length_scan_min=LENGTH_SCAN_MIN,
            length_scan_max=LENGTH_SCAN_MAX,
            known_dg=ASO_KNOWN_DG,
            known_tm=ASO_KNOWN_TM,
        )

    # Extrair valores do config para uso local
    aso_seq = config.aso_sequence.upper()
    sl_seq = config.sl_sequence.upper()
    sl_len = config.sl_length
    target_start = config.aso_target_start
    target_end = config.aso_target_end
    len_min = config.length_scan_min
    len_max = config.length_scan_max
    aso_name = config.aso_name or "ASO"

    # Se aso_sequence esta vazio, auto-design via varredura de comprimento
    if not aso_seq:
        logger.info("ASO nao definido — executando auto-design via varredura de comprimento...")
        auto_variants = scan_length_variants(
            sl_sequence=sl_seq,
            sl_length=sl_len,
            length_scan_min=len_min,
            length_scan_max=len_max,
        )
        if not auto_variants:
            raise ValueError(
                f"Nenhum ASO candidato encontrado para SL '{sl_seq}' "
                f"com comprimentos {len_min}-{len_max}"
            )
        # Selecionar o melhor candidato (menor dG = ligacao mais forte)
        best_auto = auto_variants[0]
        aso_seq = best_auto["sequence"]
        target_start = best_auto["target_start"]
        target_end = best_auto["target_end"]
        logger.info(
            "Auto-design: selecionado %s (dG=%.2f, L=%d, pos %d-%d)",
            aso_seq, best_auto["dg_binding_kcal"], best_auto["length"],
            target_start, target_end,
        )
        # Atualizar config com o ASO auto-desenhado (frozen, entao recria)
        config = TargetConfig(
            organism_slug=config.organism_slug,
            species_name=config.species_name,
            sl_sequence=config.sl_sequence,
            aso_sequence=aso_seq,
            aso_name=config.aso_name or f"auto_{config.organism_slug}",
            aso_target_start=target_start,
            aso_target_end=target_end,
            length_scan_min=len_min,
            length_scan_max=len_max,
            known_dg=config.known_dg,
            known_tm=config.known_tm,
            mutation_rate=config.mutation_rate,
            generation_time_hours=config.generation_time_hours,
            sl_copy_number=config.sl_copy_number,
            dg_functional_threshold=config.dg_functional_threshold,
            aso_concentration=config.aso_concentration,
            lna_5prime=config.lna_5prime,
            lna_3prime=config.lna_3prime,
            related_sl_sequences=config.related_sl_sequences,
            divergence_time_mya=config.divergence_time_mya,
            host_transcriptome_path=config.host_transcriptome_path,
            output_dir=config.output_dir,
        )

    aso_length = len(aso_seq)

    logger.info("=" * 60)
    logger.info("MODULO 01: Paisagem Termodinamica — %s (%s)", aso_name, config.species_name)
    logger.info("=" * 60)

    envelope = create_envelope("01_thermodynamic_landscape")

    with Timer() as timer:
        # --- Passo 1: Caracterizar wildtype ---
        logger.info("Passo 1/5: Caracterizando wildtype %s...", aso_name)
        wt_profile = _characterize_aso(aso_seq)
        wt_profile["target_start"] = target_start
        wt_profile["target_end"] = target_end

        wt_dg = wt_profile["dg_binding_kcal"]
        wt_fitness = wt_profile["fitness_score"]

        logger.info(
            "  Wildtype: dG = %.2f kcal/mol, Tm = %.2f C, GC = %.2f, fitness = %.4f",
            wt_dg,
            wt_profile["tm_celsius"],
            wt_profile["gc_content"],
            wt_fitness,
        )

        # Validacao: conferir que os valores batem com os conhecidos (se disponiveis)
        if config.known_dg is not None and config.known_tm is not None:
            dg_diff = abs(wt_dg - config.known_dg)
            tm_diff = abs(wt_profile["tm_celsius"] - config.known_tm)
            if dg_diff > 0.1 or tm_diff > 0.5:
                msg = (
                    f"AVISO: valores calculados divergem dos conhecidos! "
                    f"dG diff = {dg_diff:.2f}, Tm diff = {tm_diff:.2f}"
                )
                logger.warning(msg)
                envelope["warnings"].append(msg)

        # --- Passo 2: Mutantes pontuais (Hamming-1) ---
        n_point = aso_length * 3
        logger.info("Passo 2/5: Gerando %d mutantes pontuais (Hamming-1)...", n_point)
        point_mutants = generate_point_mutants(aso_sequence=aso_seq)

        # Contar quantos sao melhores que o wildtype (por dG e por fitness)
        better_dg_single = [m for m in point_mutants if m["dg_binding_kcal"] < wt_dg]
        better_fitness_single = [m for m in point_mutants if m["fitness_score"] < wt_fitness]

        logger.info(
            "  Mutantes com dG melhor que wildtype: %d / %d",
            len(better_dg_single),
            len(point_mutants),
        )
        logger.info(
            "  Mutantes com fitness melhor que wildtype: %d / %d",
            len(better_fitness_single),
            len(point_mutants),
        )

        # --- Passo 3: Mutantes duplos (Hamming-2) ---
        from math import comb
        n_double = comb(aso_length, 2) * 9
        logger.info("Passo 3/5: Gerando %d mutantes duplos (Hamming-2)...", n_double)
        double_mutants = generate_double_mutants(aso_sequence=aso_seq)

        better_dg_double = [m for m in double_mutants if m["dg_binding_kcal"] < wt_dg]
        better_fitness_double = [m for m in double_mutants if m["fitness_score"] < wt_fitness]

        logger.info(
            "  Duplos com dG melhor que wildtype: %d / %d",
            len(better_dg_double),
            len(double_mutants),
        )

        # --- Passo 4: Varredura de comprimento ---
        logger.info(
            "Passo 4/5: Varrendo janelas de comprimento %d-%d nt...",
            len_min, len_max,
        )
        length_variants = scan_length_variants(
            sl_sequence=sl_seq,
            sl_length=sl_len,
            length_scan_min=len_min,
            length_scan_max=len_max,
        )

        # Melhor por comprimento
        best_per_length: list[dict[str, Any]] = []
        for length in range(len_min, len_max + 1):
            candidates = [v for v in length_variants if v["length"] == length]
            if candidates:
                best = min(candidates, key=lambda v: v["dg_binding_kcal"])
                best_per_length.append(best)

        global_best_variant = length_variants[0] if length_variants else None

        # Verificar se o wildtype e o melhor entre todas as janelas
        variant_is_better = (
            global_best_variant is not None
            and global_best_variant["dg_binding_kcal"] < wt_dg
        )

        logger.info(
            "  Melhor variante de comprimento: dG = %.2f (seq = %s, L = %d)",
            global_best_variant["dg_binding_kcal"] if global_best_variant else 0.0,
            global_best_variant["sequence"][:15] + "..." if global_best_variant else "N/A",
            global_best_variant["length"] if global_best_variant else 0,
        )

        # --- Passo 5: Heatmap ---
        logger.info("Passo 5/5: Construindo dados de heatmap %dx4...", aso_length)
        heatmap = build_heatmap_data(wt_dg, point_mutants, aso_sequence=aso_seq)

        # --- Conclusao ---
        is_optimal_hamming1 = len(better_dg_single) == 0
        is_optimal_hamming2 = len(better_dg_double) == 0
        is_optimal_length = not variant_is_better
        is_global_optimal = is_optimal_hamming1 and is_optimal_hamming2 and is_optimal_length

        if is_global_optimal:
            conclusion = (
                f"{aso_name} e o OTIMO GLOBAL verificado: nenhum dos {len(point_mutants)} "
                f"mutantes pontuais (Hamming-1), nenhum dos {len(double_mutants)} mutantes "
                f"duplos (Hamming-2), e nenhuma das "
                f"{len(length_variants)} janelas de comprimento variavel ({len_min}-"
                f"{len_max} nt) apresentou dG de ligacao mais favoravel que "
                f"{wt_dg:.2f} kcal/mol."
            )
        else:
            # Relatar honestamente o que encontramos
            parts: list[str] = []
            if not is_optimal_hamming1:
                best_single = point_mutants[0]
                parts.append(
                    f"{len(better_dg_single)} mutante(s) pontual(is) com dG melhor "
                    f"(melhor: {best_single['mutation']} com dG = "
                    f"{best_single['dg_binding_kcal']:.2f} kcal/mol)"
                )
            if not is_optimal_hamming2:
                best_double = double_mutants[0]
                parts.append(
                    f"{len(better_dg_double)} mutante(s) duplo(s) com dG melhor "
                    f"(melhor: {best_double['mutation']} com dG = "
                    f"{best_double['dg_binding_kcal']:.2f} kcal/mol)"
                )
            if not is_optimal_length:
                parts.append(
                    f"variante de comprimento com dG melhor: "
                    f"{global_best_variant['sequence']} "  # type: ignore[index]
                    f"(L={global_best_variant['length']}, "  # type: ignore[index]
                    f"dG={global_best_variant['dg_binding_kcal']:.2f} kcal/mol)"  # type: ignore[index]
                )

            conclusion = (
                f"{aso_name} NAO e o otimo global. Encontrado(s): "
                + "; ".join(parts)
                + f". O wildtype tem dG = {wt_dg:.2f} kcal/mol."
            )

        logger.info("CONCLUSAO: %s", conclusion)

    # --- Montar envelope ---
    envelope["runtime_seconds"] = timer.elapsed
    envelope["status"] = "success"
    envelope["summary"]["conclusion"] = conclusion
    envelope["summary"]["key_metrics"] = {
        "wildtype_dg_kcal": wt_dg,
        "wildtype_tm_celsius": wt_profile["tm_celsius"],
        "wildtype_fitness": wt_fitness,
        "point_mutants_evaluated": len(point_mutants),
        "point_mutants_better": len(better_dg_single),
        "double_mutants_evaluated": len(double_mutants),
        "double_mutants_better": len(better_dg_double),
        "length_variants_evaluated": len(length_variants),
        "is_global_optimal": is_global_optimal,
    }

    envelope["data"] = {
        "wildtype": wt_profile,
        "point_mutants": {
            "total": len(point_mutants),
            "better_than_wildtype": len(better_dg_single),
            "better_fitness_than_wildtype": len(better_fitness_single),
            "best_mutant": point_mutants[0] if point_mutants else None,
            "worst_mutant": point_mutants[-1] if point_mutants else None,
            "mutants": point_mutants,
        },
        "length_variants": {
            "total_evaluated": len(length_variants),
            "lengths_scanned": list(range(len_min, len_max + 1)),
            "best_per_length": best_per_length,
            "global_best": global_best_variant,
            "is_optimal": is_optimal_length,
        },
        "heatmap": heatmap,
        "double_mutants": {
            "total": len(double_mutants),
            "sampled": len(double_mutants),  # varredura exaustiva, nao amostrada
            "better_than_wildtype": len(better_dg_double),
            "better_fitness_than_wildtype": len(better_fitness_double),
            "best_double_mutant": double_mutants[0] if double_mutants else None,
            "worst_double_mutant": double_mutants[-1] if double_mutants else None,
        },
    }

    # Gravar resultado — caminho depende do organismo
    if config.organism_slug != "leishmania_infantum":
        # Organismo diferente: gravar com prefixo do slug no output_dir do config
        output_name = f"{config.organism_slug}_01_thermodynamic_landscape"
        config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = write_result(envelope, module_name=output_name)
    else:
        # L. infantum: comportamento original (retrocompativel)
        output_path = write_result(envelope)

    logger.info("Resultado gravado em: %s", output_path)
    logger.info("Tempo de execucao: %.2f segundos", timer.elapsed)

    return envelope


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
