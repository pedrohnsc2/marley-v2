"""Funcoes termodinamicas compartilhadas — modelo nearest-neighbor SantaLucia 1998.

Reimplementacao independente das funcoes de rna_entropy/08_aso_design.py
para manter aso_math/ desacoplado. As tabelas NN sao constantes cientificas
publicadas, nao logica de negocio — duplicar e seguro.

Validacao: compute_tm("ACAGAAACTGATACTTATATAGCGT") deve retornar 68.48
           compute_dg("ACAGAAACTGATACTTATATAGCGT") deve retornar -27.97
"""

from __future__ import annotations

import math

from aso_math.config import (
    ASO_CONCENTRATION,
    COMPLEMENT,
    NN_INIT,
    NN_PARAMS,
    R_GAS,
    T_PHYSIOLOGICAL,
)


# ---------------------------------------------------------------------------
# Utilidades de sequencia
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    """Retorna o complemento reverso de uma sequencia DNA/RNA.

    Usa o mapeamento COMPLEMENT: A<->T, G<->C, U->A.
    Bases desconhecidas sao mantidas.
    """
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq.upper()))


def gc_content(seq: str) -> float:
    """Calcula fracao de GC em [0.0, 1.0]."""
    if not seq:
        return 0.0
    upper = seq.upper()
    return (upper.count("G") + upper.count("C")) / len(upper)


# ---------------------------------------------------------------------------
# Lookup nearest-neighbor
# ---------------------------------------------------------------------------


def _get_nn_params(dinucleotide: str) -> tuple[float, float]:
    """Busca parametros nearest-neighbor para um dinucleotideo.

    Tenta a representacao forward (XY/X'Y') e depois o complemento reverso.
    Retorna (0.0, 0.0) se nao encontrar — isso significa que o dinucleotideo
    nao tem par Watson-Crick valido (ex: mismatch).

    NOTA: Para mutantes com mismatches, (0.0, 0.0) SUBestima o custo real
    do mismatch (mismatches tem dG positivo). Isso faz mutantes parecerem
    MAIS favoraveis do que realmente sao, o que torna a prova de otimalidade
    do wildtype MAIS FORTE (conservadora).
    """
    base1, base2 = dinucleotide[0], dinucleotide[1]
    comp1 = COMPLEMENT.get(base1, base1)
    comp2 = COMPLEMENT.get(base2, base2)

    # Forward: XY/X'Y'
    forward_key = f"{base1}{base2}/{comp1}{comp2}"
    if forward_key in NN_PARAMS:
        return NN_PARAMS[forward_key]

    # Reverse complement: Y'X'/YX
    reverse_key = f"{comp2}{comp1}/{base2}{base1}"
    if reverse_key in NN_PARAMS:
        return NN_PARAMS[reverse_key]

    return (0.0, 0.0)


# ---------------------------------------------------------------------------
# Calculo de Tm (temperatura de melting)
# ---------------------------------------------------------------------------


def compute_tm(sequence: str) -> float:
    """Calcula Tm pelo modelo nearest-neighbor (SantaLucia 1998).

    Para sequencias >= 14 nt:
        Tm = (dH * 1000) / (dS + R * ln(Ct/4)) - 273.15

    Para sequencias < 14 nt (regra de Wallace):
        Tm = 2*(A+T) + 4*(G+C)

    Args:
        sequence: Oligonucleotideo DNA (case-insensitive).

    Returns:
        Tm em graus Celsius, arredondado a 2 casas decimais.
    """
    seq = sequence.upper()

    if len(seq) < 14:
        count_at = seq.count("A") + seq.count("T")
        count_gc = seq.count("G") + seq.count("C")
        return float(2 * count_at + 4 * count_gc)

    # Soma parametros NN ao longo dos dinucleotideos consecutivos
    total_dh = NN_INIT[0]  # kcal/mol — entalpia de iniciacao
    total_ds = NN_INIT[1]  # cal/(mol*K) — entropia de iniciacao

    for i in range(len(seq) - 1):
        dh, ds = _get_nn_params(seq[i:i + 2])
        total_dh += dh
        total_ds += ds

    # Tm = (dH * 1000) / (dS + R * ln(Ct/4)) - 273.15
    ln_ct_over_4 = math.log(ASO_CONCENTRATION / 4.0)
    denominator = total_ds + R_GAS * ln_ct_over_4

    if abs(denominator) < 1e-12:
        return 55.0  # fallback seguro

    tm_kelvin = (total_dh * 1000.0) / denominator
    return round(tm_kelvin - 273.15, 2)


# ---------------------------------------------------------------------------
# Calculo de dG (energia livre de binding)
# ---------------------------------------------------------------------------


def compute_dg(sequence: str) -> float:
    """Calcula delta_G de hibridizacao ASO-alvo a 37 C (310.15 K).

    dG = dH - T * dS/1000

    Onde:
        dH = soma de entalpias NN + iniciacao (kcal/mol)
        dS = soma de entropias NN + iniciacao (cal/(mol*K))
        T = 310.15 K (37 C)

    Valor mais negativo = ligacao mais forte.

    Args:
        sequence: Sequencia do ASO (5'->3').

    Returns:
        delta_G em kcal/mol, arredondado a 2 casas decimais.
    """
    seq = sequence.upper()

    total_dh = NN_INIT[0]
    total_ds = NN_INIT[1]

    for i in range(len(seq) - 1):
        dh, ds = _get_nn_params(seq[i:i + 2])
        total_dh += dh
        total_ds += ds

    delta_g = total_dh - T_PHYSIOLOGICAL * (total_ds / 1000.0)
    return round(delta_g, 2)


# ---------------------------------------------------------------------------
# Deteccao de hairpin e self-dimer
# ---------------------------------------------------------------------------


def compute_hairpin_dg(sequence: str) -> float:
    """Calcula delta_G da formacao de hairpin mais estavel.

    Para cada par (i, j) com j > i + 3 (loop minimo de 3 nt),
    encontra o stem mais longo e calcula dG do stem usando NN params.
    Penalidade de loop: +3.0 kcal/mol para loops de 3-4 nt (Mathews 1999).

    Retorna o dG mais negativo (hairpin mais estavel) encontrado.
    Valor mais negativo = hairpin mais problematico.
    """
    seq = sequence.upper()
    n = len(seq)
    best_dg = 0.0  # sem hairpin = 0.0

    for i in range(n):
        for j in range(i + 4, n):
            # Contar pares complementares consecutivos: i+k com j-k
            stem_pairs: list[str] = []
            k = 0
            while (i + k) < (j - k):
                if COMPLEMENT.get(seq[i + k]) == seq[j - k]:
                    stem_pairs.append(seq[i + k : i + k + 1])
                    k += 1
                else:
                    break

            if len(stem_pairs) < 2:
                continue

            # Calcular dG do stem usando NN params
            stem_dh = NN_INIT[0]
            stem_ds = NN_INIT[1]
            for p in range(len(stem_pairs) - 1):
                dinuc = stem_pairs[p] + stem_pairs[p + 1]
                dh, ds = _get_nn_params(dinuc)
                stem_dh += dh
                stem_ds += ds

            stem_dg = stem_dh - T_PHYSIOLOGICAL * (stem_ds / 1000.0)

            # Penalidade de loop (3-4 nt: +3.0, 5-6: +3.5, 7+: +4.0)
            loop_size = (j - k) - (i + k) + 1
            if loop_size <= 4:
                loop_penalty = 3.0
            elif loop_size <= 6:
                loop_penalty = 3.5
            else:
                loop_penalty = 4.0

            hairpin_dg = stem_dg + loop_penalty

            if hairpin_dg < best_dg:
                best_dg = hairpin_dg

    return round(best_dg, 2)


def compute_self_dimer_dg(sequence: str) -> float:
    """Calcula delta_G de auto-dimerizacao (melhor alinhamento parcial).

    Desliza a sequencia contra si mesma em todos os offsets,
    encontra o alinhamento com mais pares complementares consecutivos,
    e calcula dG desse alinhamento.

    Retorna o dG mais negativo encontrado.
    """
    seq = sequence.upper()
    n = len(seq)
    rc = reverse_complement(seq)
    best_dg = 0.0

    # Testar todos os offsets de alinhamento
    for offset in range(-(n - 2), n - 1):
        # Encontrar regiao de overlap
        if offset >= 0:
            start_seq = offset
            start_rc = 0
        else:
            start_seq = 0
            start_rc = -offset

        overlap_len = min(n - start_seq, n - start_rc)
        if overlap_len < 2:
            continue

        # Encontrar maior run de matches consecutivos
        best_run_start = -1
        best_run_len = 0
        current_start = -1
        current_len = 0

        for i in range(overlap_len):
            if seq[start_seq + i] == rc[start_rc + i]:
                if current_start < 0:
                    current_start = i
                current_len += 1
            else:
                if current_len > best_run_len:
                    best_run_len = current_len
                    best_run_start = current_start
                current_start = -1
                current_len = 0

        if current_len > best_run_len:
            best_run_len = current_len
            best_run_start = current_start

        if best_run_len < 2:
            continue

        # Calcular dG do run de matches
        run_dh = NN_INIT[0]
        run_ds = NN_INIT[1]
        for i in range(best_run_len - 1):
            pos = start_seq + best_run_start + i
            dinuc = seq[pos:pos + 2]
            dh, ds = _get_nn_params(dinuc)
            run_dh += dh
            run_ds += ds

        run_dg = run_dh - T_PHYSIOLOGICAL * (run_ds / 1000.0)
        if run_dg < best_dg:
            best_dg = run_dg

    return round(best_dg, 2)
