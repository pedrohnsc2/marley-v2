"""C0 -- Verificacao de ortogonalidade do SL RNA entre L. infantum e L. tarentolae.

VERIFICACAO MANDATORIA antes de qualquer resultado da Plataforma C.

O ASO MRL-ASO-001 foi desenhado para se ligar ao SL RNA de L. infantum,
inibindo o trans-splicing e matando o parasita patogenico. Porem,
L. tarentolae (nosso vetor de expressao) tambem possui SL RNA conservado.

Se o ASO tiver afinidade similar pelo SL RNA de L. tarentolae, ele pode
matar o vetor de expressao -- tornando a terapia combinada (ASO + vacina)
inviavel sem modificacoes.

Este modulo calcula:
    1. Alinhamento local e complementaridade ASO vs L. tarentolae SL RNA
    2. DeltaG de ligacao ASO:L. infantum vs ASO:L. tarentolae (SantaLucia 1998)
    3. Diferenca de DeltaG como indicador de seletividade
    4. Veredito: PASS se diferenca >= 2 kcal/mol, WARNING caso contrario

Referencia termodinamica: SantaLucia J Jr. (1998) PNAS 95(4):1460-1465
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

from aso_math.config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    COMPLEMENT,
    NN_INIT,
    NN_PARAMS,
    R_GAS,
    SL_SEQUENCE,
    T_PHYSIOLOGICAL,
    ASO_CONCENTRATION,
)

# ---------------------------------------------------------------------------
# Sequencia do SL RNA de L. tarentolae
# ---------------------------------------------------------------------------

# Mini-exon do Spliced Leader RNA de L. tarentolae (39 nt).
# Difere de L. infantum nas posicoes 35-36: "TT" vs "AT".
# Fonte: Breitling et al. (2002) Protein Expr Purif 25(2):209-218
# e GenBank accession M10983 (SL RNA gene cluster de L. tarentolae).
SL_TARENTOLAE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG"
SL_TARENTOLAE_LENGTH: Final[int] = len(SL_TARENTOLAE)

# Limiar minimo de diferenca de DeltaG para considerar o ASO seletivo
# para L. infantum sobre L. tarentolae. Baseado na sensibilidade
# termodinamica: 2 kcal/mol corresponde a ~20x diferenca em Kd a 37C
# (DeltaDeltaG = -RT ln(Kd1/Kd2)).
DG_SELECTIVITY_THRESHOLD: Final[float] = 2.0  # kcal/mol


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass
class OrthogonalityResult:
    """Resultado da verificacao de ortogonalidade ASO vs SL RNAs."""

    # Sequencias de entrada
    aso_sequence: str
    sl_infantum: str
    sl_tarentolae: str

    # Regiao alvo no SL RNA de L. infantum (onde o ASO foi desenhado)
    target_start: int
    target_end: int
    target_region_infantum: str

    # Regiao correspondente em L. tarentolae
    target_region_tarentolae: str

    # Alinhamento posicao-a-posicao
    alignment_length: int
    matches: int
    mismatches: int
    mismatch_positions: list[int]
    complementarity_infantum: float
    complementarity_tarentolae: float

    # Termodinamica (SantaLucia 1998, DNA/DNA)
    dg_infantum: float          # kcal/mol -- ASO vs alvo em L. infantum
    dg_tarentolae: float        # kcal/mol -- ASO vs regiao em L. tarentolae
    dh_infantum: float          # kcal/mol
    dh_tarentolae: float        # kcal/mol
    ds_infantum: float          # cal/(mol*K)
    ds_tarentolae: float        # cal/(mol*K)
    tm_infantum: float          # Celsius
    tm_tarentolae: float        # Celsius
    delta_dg: float             # |DG_infantum - DG_tarentolae| em kcal/mol
    kd_ratio: float             # Kd_tarentolae / Kd_infantum (fold selectivity)

    # Veredito
    passes_threshold: bool
    threshold_kcal: float
    verdict: str
    warnings: list[str] = field(default_factory=list)
    clinical_implications: str = ""


# ---------------------------------------------------------------------------
# Funcoes termodinamicas (espelhando aso_math para consistencia)
# ---------------------------------------------------------------------------


def _reverse_complement(seq: str) -> str:
    """Retorna o complemento reverso de uma sequencia de DNA.

    O ASO (DNA) se liga na direcao 3'->5' ao RNA alvo, entao
    precisamos do complemento reverso para calcular os pares.

    Args:
        seq: Sequencia de DNA (aceita T e U).

    Returns:
        Complemento reverso em DNA (sempre T, nunca U).
    """
    return "".join(COMPLEMENT.get(nt, "N") for nt in reversed(seq.upper()))


def _get_complement(nt: str) -> str:
    """Retorna o complemento Watson-Crick de um nucleotideo.

    Args:
        nt: Nucleotideo unico (A, T, G, C, U).

    Returns:
        Complemento em DNA.
    """
    return COMPLEMENT.get(nt.upper(), "N")


def _calculate_nn_thermodynamics(
    seq1: str,
    seq2: str,
) -> tuple[float, float, float, float]:
    """Calcula DeltaH, DeltaS, Tm e DeltaG pelo modelo nearest-neighbor.

    Usa os parametros de SantaLucia (1998) para duplexes DNA/DNA.
    Para cada par de dinucleotideos adjacentes, soma as contribuicoes
    entalpicas e entropicas. Inclui correcao de iniciacao.

    A Tm e calculada para a concentracao de ASO configurada (250 nM):
        Tm = DeltaH / (DeltaS + R * ln(Ct/4)) - 273.15

    O DeltaG e calculado a temperatura fisiologica (37C):
        DeltaG = DeltaH - T * DeltaS

    Args:
        seq1: Fita 5'->3' (ASO ou seu complemento reverso).
        seq2: Fita 3'->5' complementar (alvo).
              Deve ter mesmo comprimento que seq1.

    Returns:
        Tupla (DeltaH_kcal, DeltaS_cal, Tm_celsius, DeltaG_kcal).
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            f"Sequencias devem ter mesmo comprimento: {len(seq1)} vs {len(seq2)}"
        )

    total_h, total_s = NN_INIT  # Correcao de iniciacao

    for i in range(len(seq1) - 1):
        # Par de dinucleotideos: XY/X'Y' onde X' e Y' sao complementos
        dn_fwd = seq1[i] + seq1[i + 1]
        dn_rev = seq2[i] + seq2[i + 1]
        key = f"{dn_fwd}/{dn_rev}"

        if key in NN_PARAMS:
            dh, ds = NN_PARAMS[key]
        else:
            # Tentar o complemento reverso do par
            rc_key = f"{_reverse_complement(dn_rev)}/{_reverse_complement(dn_fwd)}"
            if rc_key in NN_PARAMS:
                dh, ds = NN_PARAMS[rc_key]
            else:
                # Par nao encontrado -- usar media conservadora
                # (mismatch reduz estabilidade)
                dh, ds = -5.0, -14.0

        total_h += dh
        total_s += ds

    # DeltaS em cal/(mol*K) -- converter para kcal para Tm
    # Tm = DeltaH / (DeltaS + R * ln(Ct/4)) - 273.15
    ds_for_tm = total_s + R_GAS * math.log(ASO_CONCENTRATION / 4.0)
    if abs(ds_for_tm) < 1e-10:
        tm_celsius = 0.0
    else:
        tm_kelvin = (total_h * 1000.0) / ds_for_tm  # H em kcal, S em cal
        tm_celsius = tm_kelvin - 273.15

    # DeltaG a 37C (temperatura fisiologica)
    dg = total_h - (T_PHYSIOLOGICAL * total_s / 1000.0)

    return total_h, total_s, round(tm_celsius, 2), round(dg, 2)


def _align_aso_to_target(aso: str, target_region: str) -> tuple[int, int, list[int]]:
    """Alinha o ASO (complemento reverso) contra a regiao alvo.

    O ASO e antisenso -- se liga ao RNA na direcao antiparalela.
    Para avaliar complementaridade, comparamos o complemento reverso
    do ASO com a fita sense do alvo.

    Args:
        aso: Sequencia do ASO 5'->3'.
        target_region: Regiao do SL RNA (fita sense).

    Returns:
        Tupla (matches, mismatches, lista de posicoes com mismatch).
    """
    # O complemento reverso do ASO deve se alinhar com o alvo
    aso_rc = _reverse_complement(aso)

    # Truncar ao comprimento minimo (podem diferir ligeiramente)
    align_len = min(len(aso_rc), len(target_region))
    aso_rc = aso_rc[:align_len]
    target_sub = target_region[:align_len]

    matches = 0
    mismatches = 0
    mismatch_pos: list[int] = []

    for i in range(align_len):
        if aso_rc[i] == target_sub[i]:
            matches += 1
        else:
            mismatches += 1
            mismatch_pos.append(i)

    return matches, mismatches, mismatch_pos


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def verify_sl_rna_orthogonality() -> OrthogonalityResult:
    """Executa a verificacao completa de ortogonalidade ASO vs SL RNAs.

    Etapas:
        1. Extrair a regiao alvo do ASO em ambos os SL RNAs
        2. Alinhar e contar matches/mismatches
        3. Calcular termodinamica de ligacao para ambos
        4. Avaliar seletividade pela diferenca de DeltaG
        5. Emitir veredito

    Returns:
        OrthogonalityResult com todos os dados calculados.

    Nota:
        Esta funcao DEVE ser chamada antes de qualquer outro calculo
        da Plataforma C. Se o resultado for FAIL, toda a plataforma
        deve reportar o problema de forma proeminente.
    """
    # Regiao alvo do ASO no SL RNA de L. infantum
    target_infantum = SL_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END]

    # Regiao correspondente em L. tarentolae (mesmas posicoes)
    target_tarentolae = SL_TARENTOLAE[ASO_TARGET_START:ASO_TARGET_END]

    # --- Alinhamento e complementaridade ---
    # Contra L. infantum (deve ser perfeito -- e o alvo original)
    matches_inf, mismatches_inf, mismatch_pos_inf = _align_aso_to_target(
        ASO_SEQUENCE, target_infantum
    )
    complementarity_inf = matches_inf / len(target_infantum) if target_infantum else 0.0

    # Contra L. tarentolae
    matches_tar, mismatches_tar, mismatch_pos_tar = _align_aso_to_target(
        ASO_SEQUENCE, target_tarentolae
    )
    align_len = min(len(ASO_SEQUENCE), len(target_tarentolae))
    complementarity_tar = matches_tar / align_len if align_len > 0 else 0.0

    # --- Termodinamica (SantaLucia 1998) ---
    # ASO vs L. infantum
    aso_rc = _reverse_complement(ASO_SEQUENCE)
    # Para os calculos NN, usamos a regiao complementar
    # Fita 1: ASO_rc (5'->3'), Fita 2: complemento do ASO_rc = alvo sense

    # Calcular para L. infantum (match perfeito esperado)
    dh_inf, ds_inf, tm_inf, dg_inf = _calculate_nn_thermodynamics(
        aso_rc[:len(target_infantum)],
        target_infantum,
    )

    # Calcular para L. tarentolae (com mismatches)
    dh_tar, ds_tar, tm_tar, dg_tar = _calculate_nn_thermodynamics(
        aso_rc[:len(target_tarentolae)],
        target_tarentolae,
    )

    # --- Seletividade ---
    # Delta DeltaG: diferenca absoluta de afinidade
    # DeltaG mais negativo = ligacao mais forte
    # Se DG_infantum << DG_tarentolae, o ASO prefere L. infantum (bom!)
    delta_dg = abs(dg_inf - dg_tar)

    # Razao de Kd: Kd = exp(DG / RT)
    # Kd_ratio = Kd_tar / Kd_inf = exp((DG_tar - DG_inf) / (R * T))
    dg_diff = dg_tar - dg_inf  # Positivo se ASO liga mais forte a L. infantum
    rt_kcal = (R_GAS * T_PHYSIOLOGICAL) / 1000.0  # R em kcal/(mol*K)
    kd_ratio = math.exp(dg_diff / rt_kcal) if rt_kcal > 0 else 1.0

    # --- Veredito ---
    passes = delta_dg >= DG_SELECTIVITY_THRESHOLD
    warnings: list[str] = []

    if passes:
        verdict = (
            f"PASS -- ASO MRL-ASO-001 e seletivo para L. infantum sobre "
            f"L. tarentolae (DeltaDeltaG = {delta_dg:.2f} kcal/mol >= "
            f"{DG_SELECTIVITY_THRESHOLD} kcal/mol). "
            f"A terapia combinada ASO + vacina L. tarentolae e viavel."
        )
    else:
        verdict = (
            f"WARNING -- Diferenca de DeltaG insuficiente "
            f"({delta_dg:.2f} kcal/mol < {DG_SELECTIVITY_THRESHOLD} kcal/mol). "
            f"O ASO MRL-ASO-001 pode afetar a viabilidade de L. tarentolae "
            f"como vetor de expressao. REQUER AVALIACAO EXPERIMENTAL."
        )
        warnings.append(
            f"CRITICO: DeltaDeltaG = {delta_dg:.2f} kcal/mol esta abaixo do "
            f"limiar de {DG_SELECTIVITY_THRESHOLD} kcal/mol. "
            f"Isso significa que o ASO pode ter afinidade significativa pelo "
            f"SL RNA de L. tarentolae, potencialmente inibindo o trans-splicing "
            f"do vetor de expressao."
        )

    # Adicionar avisos sobre mismatches encontrados
    if mismatches_tar > 0:
        diff_positions = []
        for pos in mismatch_pos_tar:
            abs_pos = ASO_TARGET_START + pos
            diff_positions.append(
                f"pos {abs_pos}: infantum={target_infantum[pos]}, "
                f"tarentolae={target_tarentolae[pos]}"
            )
        warnings.append(
            f"Mismatches ASO:L.tarentolae em {mismatches_tar} posicao(oes): "
            + "; ".join(diff_positions)
        )
    else:
        warnings.append(
            "ALERTA: Zero mismatches detectados -- o ASO tem complementaridade "
            "perfeita com ambos os SL RNAs. A seletividade depende inteiramente "
            "de fatores cineticos e de acessibilidade estrutural."
        )

    # --- Implicacoes clinicas ---
    if passes:
        clinical = (
            "COMPATIVEL: A terapia combinada ASO + vacina L. tarentolae e viavel. "
            "Protocolo sugerido: (1) administrar ASO para reduzir carga parasitaria "
            "de L. infantum, (2) apos clearance do ASO (~5 meias-vidas, ~48h para "
            "ASOs fosforotioados), administrar vacina L. tarentolae. "
            "A diferenca termodinamica garante que concentracoes residuais de ASO "
            "nao comprometem significativamente o vetor vacinal."
        )
    else:
        clinical = (
            "INCOMPATIVEL sem modificacoes. Opcoes: "
            "(1) Administrar vacina L. tarentolae ANTES do tratamento com ASO; "
            "(2) Usar intervalo longo (>1 semana) entre ASO e vacinacao; "
            "(3) Modificar o ASO com LNA/2'-OMe nas posicoes de mismatch para "
            "aumentar seletividade; "
            "(4) Considerar vetor de expressao alternativo (E. coli, Platform B)."
        )

    return OrthogonalityResult(
        aso_sequence=ASO_SEQUENCE,
        sl_infantum=SL_SEQUENCE,
        sl_tarentolae=SL_TARENTOLAE,
        target_start=ASO_TARGET_START,
        target_end=ASO_TARGET_END,
        target_region_infantum=target_infantum,
        target_region_tarentolae=target_tarentolae,
        alignment_length=align_len,
        matches=matches_tar,
        mismatches=mismatches_tar,
        mismatch_positions=mismatch_pos_tar,
        complementarity_infantum=round(complementarity_inf, 4),
        complementarity_tarentolae=round(complementarity_tar, 4),
        dg_infantum=dg_inf,
        dg_tarentolae=dg_tar,
        dh_infantum=round(dh_inf, 2),
        dh_tarentolae=round(dh_tar, 2),
        ds_infantum=round(ds_inf, 2),
        ds_tarentolae=round(ds_tar, 2),
        tm_infantum=tm_inf,
        tm_tarentolae=tm_tar,
        delta_dg=round(delta_dg, 2),
        kd_ratio=round(kd_ratio, 2),
        passes_threshold=passes,
        threshold_kcal=DG_SELECTIVITY_THRESHOLD,
        verdict=verdict,
        warnings=warnings,
        clinical_implications=clinical,
    )


def to_dict(result: OrthogonalityResult) -> dict:
    """Serializa o resultado de ortogonalidade para JSON.

    Args:
        result: Resultado da verificacao de ortogonalidade.

    Returns:
        Dicionario serializavel para JSON.
    """
    return {
        "sequences": {
            "aso_mrl_001": result.aso_sequence,
            "sl_rna_l_infantum": result.sl_infantum,
            "sl_rna_l_tarentolae": result.sl_tarentolae,
        },
        "target_region": {
            "start": result.target_start,
            "end": result.target_end,
            "infantum": result.target_region_infantum,
            "tarentolae": result.target_region_tarentolae,
        },
        "alignment": {
            "length": result.alignment_length,
            "matches": result.matches,
            "mismatches": result.mismatches,
            "mismatch_positions": result.mismatch_positions,
            "complementarity_infantum": result.complementarity_infantum,
            "complementarity_tarentolae": result.complementarity_tarentolae,
        },
        "thermodynamics": {
            "l_infantum": {
                "dH_kcal_mol": result.dh_infantum,
                "dS_cal_mol_K": result.ds_infantum,
                "Tm_celsius": result.tm_infantum,
                "dG_37C_kcal_mol": result.dg_infantum,
            },
            "l_tarentolae": {
                "dH_kcal_mol": result.dh_tarentolae,
                "dS_cal_mol_K": result.ds_tarentolae,
                "Tm_celsius": result.tm_tarentolae,
                "dG_37C_kcal_mol": result.dg_tarentolae,
            },
            "selectivity": {
                "delta_dG_kcal_mol": result.delta_dg,
                "Kd_ratio_tarentolae_over_infantum": result.kd_ratio,
                "threshold_kcal_mol": result.threshold_kcal,
            },
        },
        "verdict": {
            "passes_threshold": result.passes_threshold,
            "summary": result.verdict,
            "warnings": result.warnings,
            "clinical_implications": result.clinical_implications,
        },
    }
