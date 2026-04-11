"""Modulo A — Estabilidade de MRL-ASO-001 em pH acido do fagolisossomo.

Pergunta central: MRL-ASO-001 permanece estavel e funcional no
ambiente hostil do fagolisossomo de macrofagos caninos (pH 4.5)?

O fagolisossomo e o compartimento onde amastigotas de L. infantum
residem e se replicam. Este ambiente apresenta:
- pH 4.5-5.0 (vs pH 7.4 no sangue)
- Nucleases acidas (catepsinas, DNase II)
- Proteases e hidrolases lisossomais
- Alta concentracao de ions metalicos (Fe2+, Zn2+)

O ASO MRL-ASO-001 precisa FUNCIONAR neste ambiente porque o alvo
(SL RNA) esta no MESMO compartimento — nao e necessario escape
endossomal, mas e necessario resistir a degradacao e manter
afinidade de ligacao suficiente.

Este modulo executa quatro analises independentes:
1. Perfil de estabilidade termodinamica em gradiente de pH
2. Resistencia a nucleases por quimica de backbone
3. Estabilidade conformacional de nucleotideos LNA
4. Comparacao com ASOs clinicamente aprovados pelo FDA

Saida: JSON + relatorio Markdown em aso_delivery/module_a_stability/results/

Referencias principais:
- SantaLucia J Jr (1998) PNAS 95(4):1460-1465
- Siegfried NA et al. (2010) Biochemistry 49(15):3225-3236
- Eckstein F (2014) Nucleic Acids Res 42(6):3777-3788
- Vester B, Wengel J (2004) Biochemistry 43(42):13233-13241
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from aso_delivery.module_a_stability.comparisons import build_comparison_table
from aso_math.config import ASO_SEQUENCE as _ASO_SEQ
from aso_delivery.module_a_stability.models import (
    K_PO_PH45,
    K_PS_PH45,
    compute_dg_at_ph,
    compute_fraction_bound,
    compute_fraction_remaining,
    compute_gapmer_degradation,
    compute_lna_conformation,
    compute_nuclease_resistance,
    compute_protonation_correction,
    compute_tm_at_ph,
)
from core.logger import get_logger

logger = get_logger("module_a_stability")

# ---------------------------------------------------------------------------
# Constantes do modulo
# ---------------------------------------------------------------------------

MODULE_NAME: Final[str] = "module_a_stability"
MODULE_VERSION: Final[str] = "1.0.0"

# Diretorio de resultados
RESULTS_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# Perfil de pH a testar — representa a jornada do ASO do sangue ate o fagolisossomo
# Cada pH corresponde a um compartimento celular especifico:
PH_PROFILE: Final[dict[str, float]] = {
    "blood_plasma": 7.4,          # sangue circulante
    "early_endosome": 6.5,        # endossomo precoce (apos captacao)
    "late_endosome": 5.0,         # endossomo tardio (fusao com lisossomo)
    "phagolysosome": 4.5,         # fagolisossomo (destino final)
}

# Especificacoes do ASO — parametros termicos conhecidos do pipeline
# (a sequencia nao e incluida neste arquivo por design)
ASO_DG_NEUTRAL: Final[float] = -28.0     # kcal/mol a pH 7.4
ASO_TM_DNA: Final[float] = 61.4          # Celsius (DNA-DNA sem LNA)
ASO_TM_LNA: Final[float] = 108.5         # Celsius (com contribuicao LNA)
ASO_LENGTH: Final[int] = 25
ASO_N_BP: Final[int] = 25                # pares de bases no duplex

# Composicao de bases do ASO (derivada da sequencia no config, sem duplica-la)
ASO_N_CYTOSINES: Final[int] = _ASO_SEQ.upper().count("C")
ASO_N_ADENINES: Final[int] = _ASO_SEQ.upper().count("A")

# Arquitetura do gapmer
LNA_5PRIME: Final[int] = 5
DNA_GAP: Final[int] = 15
LNA_3PRIME: Final[int] = 5
TOTAL_LNA: Final[int] = LNA_5PRIME + LNA_3PRIME

# Limiar de funcionalidade termodinamica
# ASO perde funcao se dG > -15 kcal/mol (insuficiente para recrutar RNase H)
DG_FUNCTIONAL_THRESHOLD: Final[float] = -15.0

# Boost de Tm por LNA: cada residuo LNA adiciona ~3-8 C ao Tm
# Ref: McTigue PM et al. (2004) Biochemistry 43(18):5388-5405
# Para MRL-ASO-001, o boost total e:
# Tm_LNA - Tm_DNA = 108.5 - 61.4 = 47.1 C para 10 LNA = ~4.7 C/LNA
LNA_TM_BOOST_PER_RESIDUE: Final[float] = (ASO_TM_LNA - ASO_TM_DNA) / TOTAL_LNA


# ---------------------------------------------------------------------------
# 1. Perfil de estabilidade termodinamica por pH
# ---------------------------------------------------------------------------


def analyze_ph_stability() -> dict[str, Any]:
    """Calcula perfil completo de estabilidade termodinamica em gradiente de pH.

    Para cada pH na jornada intracelular (7.4 -> 6.5 -> 5.0 -> 4.5),
    calcula:
    - Correcao de dG por protonacao de C e A
    - dG corrigido
    - Tm estimado
    - Fracao de ASO ligado ao alvo

    A logica biologica: conforme o ASO desce pelo pathway endossomal,
    o pH diminui e a protonacao de bases aumenta, desestabilizando
    o duplex. A pergunta e se o dG se mantem negativo o suficiente
    para funcao terapeutica no destino final (pH 4.5).

    Returns:
        Dicionario com perfil de estabilidade por pH.
    """
    logger.info("Analise 1/4: Perfil de estabilidade termodinamica por pH")
    logger.info("  Composicao: %d citosinas, %d adeninas em %d nt",
                ASO_N_CYTOSINES, ASO_N_ADENINES, ASO_LENGTH)

    ph_results: dict[str, Any] = {}

    for compartment, ph in PH_PROFILE.items():
        # Correcao de protonacao
        correction = compute_protonation_correction(
            ph=ph,
            n_cytosines=ASO_N_CYTOSINES,
            n_adenines=ASO_N_ADENINES,
        )

        # dG corrigido
        dg_corrected = compute_dg_at_ph(
            dg_neutral=ASO_DG_NEUTRAL,
            ph=ph,
            n_cytosines=ASO_N_CYTOSINES,
            n_adenines=ASO_N_ADENINES,
        )

        # Tm estimado (usando Tm com LNA como base)
        tm_corrected = compute_tm_at_ph(
            tm_neutral=ASO_TM_LNA,
            dg_neutral=ASO_DG_NEUTRAL,
            dg_at_ph=dg_corrected,
            n_bp=ASO_N_BP,
        )

        # Fracao ligada
        f_bound = compute_fraction_bound(dg_corrected)

        # Funcionalidade mantida?
        functional = dg_corrected <= DG_FUNCTIONAL_THRESHOLD

        entry = {
            "ph": ph,
            "compartment": compartment,
            "delta_g_kcal": dg_corrected,
            "delta_g_correction_kcal": correction.ddg_total,
            "tm_celsius": tm_corrected,
            "fraction_bound": f_bound,
            "functional": functional,
            "protonation_detail": {
                "fraction_protonated_C": correction.fraction_protonated_c,
                "fraction_protonated_A": correction.fraction_protonated_a,
                "ddg_cytosines": correction.ddg_cytosines,
                "ddg_adenines": correction.ddg_adenines,
            },
        }

        ph_results[f"pH_{ph}"] = entry

        logger.info(
            "  pH %.1f (%s): dG = %.2f kcal/mol (correcao +%.4f), "
            "Tm = %.1f C, f_bound = %.4f, funcional = %s",
            ph, compartment, dg_corrected, correction.ddg_total,
            tm_corrected, f_bound, functional,
        )

    # Verificacao global: funcional no fagolisossomo?
    phagolysosome = ph_results["pH_4.5"]
    overall_functional = phagolysosome["functional"]

    logger.info(
        "  CONCLUSAO pH: dG no fagolisossomo = %.2f kcal/mol (limiar: %.1f) -> %s",
        phagolysosome["delta_g_kcal"],
        DG_FUNCTIONAL_THRESHOLD,
        "FUNCIONAL" if overall_functional else "NAO FUNCIONAL",
    )

    return {
        "profile": ph_results,
        "functional_at_target_ph": overall_functional,
        "dg_threshold_kcal": DG_FUNCTIONAL_THRESHOLD,
    }


# ---------------------------------------------------------------------------
# 2. Resistencia a nucleases
# ---------------------------------------------------------------------------


def analyze_nuclease_resistance() -> dict[str, Any]:
    """Compara resistencia a nucleases entre diferentes quimicas de backbone.

    Modela degradacao por cinetica de primeira ordem para tres cenarios:
    1. PO (fosfodiester nativo) — baseline, rapidamente degradado
    2. PS (fosforotioato puro) — 500x mais resistente
    3. LNA-DNA-LNA gapmer PS — protecao adicional nos flancos

    A motivacao: o fagolisossomo contem nucleases acidas (DNase II,
    catepsinas) que degradam DNA/RNA. O backbone PS substitui um oxigenio
    nao-ponte por enxofre, impedindo a coordenacao do ion metalico
    necessario para catalise. LNA adiciona protecao mecanica adicional
    (o acucar travado impede acesso ao backbone).

    Returns:
        Dicionario com perfis de resistencia e comparacao.
    """
    logger.info("Analise 2/4: Resistencia a nucleases por quimica de backbone")

    # PO — backbone nativo (controle negativo)
    po_profile = compute_nuclease_resistance("PO (phosphodiester)", K_PO_PH45)
    logger.info(
        "  PO: t1/2 = %.2f h, f(24h) = %.6f -> %s",
        po_profile.half_life_hours,
        po_profile.fraction_at_24h,
        "INSUFICIENTE" if not po_profile.therapeutic_window_met else "OK",
    )

    # PS — fosforotioato puro
    ps_profile = compute_nuclease_resistance("PS (phosphorothioate)", K_PS_PH45)
    logger.info(
        "  PS: t1/2 = %.2f h, f(24h) = %.6f -> %s",
        ps_profile.half_life_hours,
        ps_profile.fraction_at_24h,
        "SUFICIENTE" if ps_profile.therapeutic_window_met else "INSUFICIENTE",
    )

    # LNA-DNA-LNA gapmer
    gapmer = compute_gapmer_degradation(
        lna_5prime=LNA_5PRIME,
        dna_gap=DNA_GAP,
        lna_3prime=LNA_3PRIME,
        k_ps=K_PS_PH45,
    )

    gapmer_profile = compute_nuclease_resistance(
        "LNA-DNA-LNA gapmer (PS backbone)",
        gapmer.k_effective,
    )
    logger.info(
        "  Gapmer: k_eff = %.6f/h, t1/2 = %.2f h, f(24h) = %.6f -> %s",
        gapmer.k_effective,
        gapmer_profile.half_life_hours,
        gapmer_profile.fraction_at_24h,
        "SUFICIENTE" if gapmer_profile.therapeutic_window_met else "INSUFICIENTE",
    )

    # Curva de degradacao temporal (0-96 horas) para todos os backbones
    time_points = [0, 1, 2, 4, 6, 12, 24, 48, 72, 96]
    degradation_curves: dict[str, list[dict[str, float]]] = {
        "PO": [],
        "PS": [],
        "LNA_gapmer": [],
    }

    for t in time_points:
        degradation_curves["PO"].append({
            "time_hours": t,
            "fraction_remaining": compute_fraction_remaining(K_PO_PH45, t),
        })
        degradation_curves["PS"].append({
            "time_hours": t,
            "fraction_remaining": compute_fraction_remaining(K_PS_PH45, t),
        })
        degradation_curves["LNA_gapmer"].append({
            "time_hours": t,
            "fraction_remaining": compute_fraction_remaining(gapmer.k_effective, t),
        })

    # PS/PO ratio de resistencia
    ps_po_ratio = K_PO_PH45 / K_PS_PH45 if K_PS_PH45 > 0 else float("inf")
    gapmer_po_ratio = K_PO_PH45 / gapmer.k_effective if gapmer.k_effective > 0 else float("inf")

    return {
        "PO_halflife_hours": po_profile.half_life_hours,
        "PS_halflife_hours": ps_profile.half_life_hours,
        "LNA_gapmer_halflife_hours": gapmer_profile.half_life_hours,
        "therapeutic_window_met": gapmer_profile.therapeutic_window_met,
        "profiles": {
            "PO": asdict(po_profile),
            "PS": asdict(ps_profile),
            "LNA_gapmer": asdict(gapmer_profile),
        },
        "gapmer_model": asdict(gapmer),
        "resistance_ratios": {
            "PS_vs_PO": round(ps_po_ratio, 1),
            "gapmer_vs_PO": round(gapmer_po_ratio, 1),
        },
        "degradation_curves": degradation_curves,
    }


# ---------------------------------------------------------------------------
# 3. Estabilidade conformacional LNA
# ---------------------------------------------------------------------------


def analyze_lna_stability() -> dict[str, Any]:
    """Avalia estabilidade conformacional dos nucleotideos LNA em pH acido.

    LNA (Locked Nucleic Acid) tem ponte metileno C2'-C4' que trava
    o acucar na conformacao C3'-endo. Esta conformacao e essencial
    para hibridizacao de alta afinidade (geometria tipo RNA-A).

    Em pH acido, a protonacao de bases pode perturbar o equilibrio
    conformacional do acucar em DNA convencional, mas a ponte
    covalente do LNA e INSENSIVEL a pH — mantem C3'-endo por
    restricao mecanica, nao por equilibrio termodinamico.

    Returns:
        Dicionario com perfis conformacionais por pH.
    """
    logger.info("Analise 3/4: Estabilidade conformacional LNA por pH")

    lna_results: dict[str, Any] = {}

    for compartment, ph in PH_PROFILE.items():
        profile = compute_lna_conformation(
            ph=ph,
            n_lna=TOTAL_LNA,
            n_dna=DNA_GAP,
        )

        entry = {
            "ph": ph,
            "compartment": compartment,
            "c3_endo_fraction_lna": profile.c3_endo_fraction_lna,
            "c3_endo_fraction_dna": profile.c3_endo_fraction_dna,
            "geometry_maintained": profile.geometry_maintained,
            "free_energy_penalty_kcal": profile.free_energy_penalty,
            "n_lna_residues": profile.n_lna_residues,
            "n_dna_residues": profile.n_dna_residues,
        }

        lna_results[f"pH_{ph}"] = entry

        logger.info(
            "  pH %.1f: LNA C3'-endo = %.1f%%, DNA C3'-endo = %.1f%%, "
            "geometria %s, penalidade = %.4f kcal/mol",
            ph,
            profile.c3_endo_fraction_lna * 100,
            profile.c3_endo_fraction_dna * 100,
            "MANTIDA" if profile.geometry_maintained else "PERDIDA",
            profile.free_energy_penalty,
        )

    # Comparacao pH 7.4 vs pH 4.5
    neutral = lna_results["pH_7.4"]
    acidic = lna_results["pH_4.5"]

    logger.info(
        "  CONCLUSAO LNA: C3'-endo LNA cai de %.1f%% (pH 7.4) para %.1f%% (pH 4.5) -> %s",
        neutral["c3_endo_fraction_lna"] * 100,
        acidic["c3_endo_fraction_lna"] * 100,
        "GEOMETRIA MANTIDA" if acidic["geometry_maintained"] else "GEOMETRIA COMPROMETIDA",
    )

    return {
        "profile": lna_results,
        "c3_endo_fraction_pH_7_4": neutral["c3_endo_fraction_lna"],
        "c3_endo_fraction_pH_4_5": acidic["c3_endo_fraction_lna"],
        "geometry_maintained": acidic["geometry_maintained"],
        "lna_tm_boost_per_residue_celsius": round(LNA_TM_BOOST_PER_RESIDUE, 2),
        "total_lna_residues": TOTAL_LNA,
        "reference": "Vester B, Wengel J (2004) Biochemistry 43(42):13233-13241",
    }


# ---------------------------------------------------------------------------
# 4. Comparacao com ASOs FDA
# ---------------------------------------------------------------------------


def analyze_fda_comparison(gapmer_half_life_hours: float) -> dict[str, Any]:
    """Compara MRL-ASO-001 com ASOs aprovados pelo FDA.

    Contextualizacao clinica: embora MRL-ASO-001 seja pre-clinico,
    a comparacao com ASOs aprovados valida que a quimica escolhida
    tem precedente regulatorio e que a estabilidade estimada esta
    dentro da faixa de compostos clinicamente viáveis.

    O comparador mais relevante e mipomersen (Kynamro), que:
    - Usa gapmer 5-10-5 com backbone PS (similar ao MRL-ASO-001)
    - Enfrenta pH 4.5 em lisossomos hepaticos
    - Tem meia-vida em tecido de ~30 dias

    Args:
        gapmer_half_life_hours: Meia-vida estimada do gapmer MRL-ASO-001.

    Returns:
        Tabela comparativa com todos os ASOs.
    """
    logger.info("Analise 4/4: Comparacao com ASOs aprovados pelo FDA")

    table = build_comparison_table(gapmer_half_life_hours)

    for name, data in table.items():
        logger.info(
            "  %s: %d nt, %s, pH %.1f, t1/2 = %.1f h",
            name,
            data["length_nt"],
            data["chemistry"],
            data["target_ph"],
            data["half_life_hours"],
        )

    return table


# ---------------------------------------------------------------------------
# Geracao de relatorio Markdown
# ---------------------------------------------------------------------------


def generate_markdown_report(results: dict[str, Any]) -> str:
    """Gera relatorio Markdown com os resultados completos.

    O relatorio e estruturado para leitura tecnica, com tabelas
    e conclusoes por secao. Inclui citacoes e contexto biologico.

    Args:
        results: Dicionario com todos os resultados das 4 analises.

    Returns:
        String com conteudo Markdown do relatorio.
    """
    lines: list[str] = []

    lines.append("# Module A: Phagolysosomal Stability Analysis of MRL-ASO-001")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"**Module version:** {MODULE_VERSION}")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(results["executive_summary"])
    lines.append("")

    # ---- Secao 1: pH Stability ----
    lines.append("## 1. Thermodynamic Stability Across pH Gradient")
    lines.append("")
    lines.append("The ASO traverses a pH gradient from blood (7.4) to phagolysosome (4.5).")
    lines.append("Protonation of cytosine N3 (pKa 4.2) and adenine N1 (pKa 3.5)")
    lines.append("destabilizes Watson-Crick base pairs at low pH.")
    lines.append("")

    ph_data = results["ph_stability_profile"]
    lines.append("| Compartment | pH | dG (kcal/mol) | Correction | Tm (C) | Fraction Bound | Functional |")
    lines.append("|---|---|---|---|---|---|---|")

    for key, entry in ph_data["profile"].items():
        lines.append(
            f"| {entry['compartment']} | {entry['ph']} | {entry['delta_g_kcal']:.2f} | "
            f"+{entry['delta_g_correction_kcal']:.4f} | {entry['tm_celsius']:.1f} | "
            f"{entry['fraction_bound']:.4f} | {'Yes' if entry['functional'] else 'No'} |"
        )

    lines.append("")
    functional_status = "PASSES" if ph_data["functional_at_target_ph"] else "FAILS"
    lines.append(f"**Verdict:** MRL-ASO-001 **{functional_status}** the pH stability test.")
    lines.append(f"dG at pH 4.5 remains well below the functional threshold of {ph_data['dg_threshold_kcal']} kcal/mol.")
    lines.append("")

    # ---- Secao 2: Nuclease Resistance ----
    nr = results["nuclease_resistance"]
    lines.append("## 2. Nuclease Resistance at pH 4.5")
    lines.append("")
    lines.append("First-order degradation kinetics: [ASO](t) = [ASO]_0 * exp(-k*t)")
    lines.append("")
    lines.append("| Backbone | k (1/h) | Half-life (h) | f(24h) | f(48h) | Therapeutic Window |")
    lines.append("|---|---|---|---|---|---|")

    for btype in ["PO", "PS", "LNA_gapmer"]:
        p = nr["profiles"][btype]
        lines.append(
            f"| {p['backbone_type']} | {p['k_degradation']:.6f} | {p['half_life_hours']:.2f} | "
            f"{p['fraction_at_24h']:.6f} | {p['fraction_at_48h']:.6f} | "
            f"{'Met' if p['therapeutic_window_met'] else 'Not met'} |"
        )

    lines.append("")
    lines.append(f"**PS/PO resistance ratio:** {nr['resistance_ratios']['PS_vs_PO']}x")
    lines.append(f"**Gapmer/PO resistance ratio:** {nr['resistance_ratios']['gapmer_vs_PO']}x")
    lines.append("")

    gapmer = nr["gapmer_model"]
    lines.append(f"**Gapmer architecture:** {gapmer['lna_5prime_length']}-{gapmer['dna_gap_length']}-{gapmer['lna_3prime_length']} "
                 f"(LNA-DNA-LNA)")
    lines.append(f"**Effective k:** {gapmer['k_effective']:.8f}/h")
    lines.append(f"**Estimated half-life:** {gapmer['half_life_hours']:.2f} hours")
    lines.append("")
    window_status = "PASSES" if nr["therapeutic_window_met"] else "FAILS"
    lines.append(f"**Verdict:** MRL-ASO-001 **{window_status}** the nuclease resistance test.")
    lines.append("")

    # ---- Secao 3: LNA Conformational Stability ----
    lna = results["lna_stability"]
    lines.append("## 3. LNA Conformational Stability")
    lines.append("")
    lines.append("LNA nucleotides maintain C3'-endo sugar pucker via covalent methylene bridge,")
    lines.append("independent of pH. This is critical for high-affinity hybridization.")
    lines.append("")
    lines.append("| pH | Compartment | LNA C3'-endo | DNA C3'-endo | Geometry | Penalty (kcal/mol) |")
    lines.append("|---|---|---|---|---|---|")

    for key, entry in lna["profile"].items():
        lines.append(
            f"| {entry['ph']} | {entry['compartment']} | "
            f"{entry['c3_endo_fraction_lna']*100:.1f}% | "
            f"{entry['c3_endo_fraction_dna']*100:.1f}% | "
            f"{'Maintained' if entry['geometry_maintained'] else 'Compromised'} | "
            f"{entry['free_energy_penalty_kcal']:.4f} |"
        )

    lines.append("")
    lines.append(f"**LNA Tm boost:** {lna['lna_tm_boost_per_residue_celsius']:.2f} C per LNA residue "
                 f"({lna['total_lna_residues']} residues total)")
    geom_status = "PASSES" if lna["geometry_maintained"] else "FAILS"
    lines.append(f"**Verdict:** MRL-ASO-001 **{geom_status}** the LNA conformational stability test.")
    lines.append("")

    # ---- Secao 4: FDA Comparison ----
    fda = results["fda_comparison"]
    lines.append("## 4. Comparison with FDA-Approved ASOs")
    lines.append("")
    lines.append("| Property | Nusinersen | Mipomersen | Inotersen | MRL-ASO-001 |")
    lines.append("|---|---|---|---|---|")

    props = ["length_nt", "backbone", "gapmer_design", "route", "target_ph",
             "half_life_hours", "mechanism"]
    prop_labels = {
        "length_nt": "Length (nt)",
        "backbone": "Backbone",
        "gapmer_design": "Gapmer design",
        "route": "Route",
        "target_ph": "Target pH",
        "half_life_hours": "Half-life (h)",
        "mechanism": "Mechanism",
    }

    for prop in props:
        label = prop_labels[prop]
        vals = [str(fda[name].get(prop, "N/A")) for name in ["nusinersen", "mipomersen", "inotersen", "mrl_aso_001"]]
        lines.append(f"| {label} | {' | '.join(vals)} |")

    lines.append("")
    lines.append("**Key insight:** MRL-ASO-001 has a unique advantage over subcutaneous ASOs like")
    lines.append("mipomersen and inotersen: the target RNA (SL RNA) and the ASO colocalize in the")
    lines.append("same acidic compartment (phagolysosome), eliminating the need for endosomal escape.")
    lines.append("")

    # ---- Overall Conclusion ----
    lines.append("## Overall Conclusion")
    lines.append("")
    lines.append(results["overall_conclusion"])
    lines.append("")
    lines.append("## References")
    lines.append("")
    lines.append("1. SantaLucia J Jr (1998) PNAS 95(4):1460-1465 — Nearest-neighbor parameters")
    lines.append("2. Siegfried NA et al. (2010) Biochemistry 49(15):3225-3236 — Base protonation at low pH")
    lines.append("3. Eckstein F (2014) Nucleic Acids Res 42(6):3777-3788 — Phosphorothioate stability")
    lines.append("4. Vester B, Wengel J (2004) Biochemistry 43(42):13233-13241 — LNA conformational locking")
    lines.append("5. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review")
    lines.append("6. Kurreck J (2003) Eur J Biochem 270(8):1628-1644 — ASO design principles")
    lines.append("7. Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51 — ASO pharmacokinetics")
    lines.append("8. McTigue PM et al. (2004) Biochemistry 43(18):5388-5405 — LNA Tm enhancement")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Orquestrador principal
# ---------------------------------------------------------------------------


def main() -> dict[str, Any]:
    """Executa todas as analises de estabilidade em pH acido.

    Fluxo:
        1. Perfil termodinamico por pH (4 compartimentos)
        2. Resistencia a nucleases (3 quimicas de backbone)
        3. Estabilidade conformacional LNA (4 pHs)
        4. Comparacao com ASOs FDA (3 compostos)
        5. Gerar relatorio Markdown
        6. Gravar resultados JSON + Markdown

    Returns:
        Envelope completo com todos os resultados.
    """
    logger.info("=" * 70)
    logger.info("MODULO A: Estabilidade no Fagolisossomo — MRL-ASO-001")
    logger.info("=" * 70)
    logger.info("Pergunta: MRL-ASO-001 permanece estavel e funcional em pH 4.5?")
    logger.info("")

    start_time = time.time()

    # --- Analise 1: Estabilidade termodinamica por pH ---
    ph_stability = analyze_ph_stability()
    logger.info("")

    # --- Analise 2: Resistencia a nucleases ---
    nuclease_resistance = analyze_nuclease_resistance()
    logger.info("")

    # --- Analise 3: Estabilidade LNA ---
    lna_stability = analyze_lna_stability()
    logger.info("")

    # --- Analise 4: Comparacao FDA ---
    gapmer_half_life = nuclease_resistance["LNA_gapmer_halflife_hours"]
    fda_comparison = analyze_fda_comparison(gapmer_half_life)
    logger.info("")

    # --- Compilar resultados ---
    elapsed = round(time.time() - start_time, 2)

    # Determinar status global
    # MRL-ASO-001 passa se TODAS as condicoes forem satisfeitas:
    all_pass = (
        ph_stability["functional_at_target_ph"]
        and nuclease_resistance["therapeutic_window_met"]
        and lna_stability["geometry_maintained"]
    )

    # Sumario executivo
    phagolysosome_dg = ph_stability["profile"]["pH_4.5"]["delta_g_kcal"]
    phagolysosome_fbound = ph_stability["profile"]["pH_4.5"]["fraction_bound"]

    if all_pass:
        executive_summary = (
            f"MRL-ASO-001 PASSES all phagolysosomal stability tests. "
            f"At pH 4.5: dG = {phagolysosome_dg:.2f} kcal/mol "
            f"(threshold: {DG_FUNCTIONAL_THRESHOLD} kcal/mol), "
            f"fraction bound = {phagolysosome_fbound:.4f}, "
            f"gapmer half-life = {gapmer_half_life:.2f} hours, "
            f"LNA geometry maintained at {lna_stability['c3_endo_fraction_pH_4_5']*100:.1f}% C3'-endo. "
            f"The LNA-DNA-LNA gapmer with full PS backbone is well-suited for "
            f"the hostile phagolysosomal environment of L. infantum amastigotes."
        )
    else:
        failures: list[str] = []
        if not ph_stability["functional_at_target_ph"]:
            failures.append("pH stability (dG above functional threshold)")
        if not nuclease_resistance["therapeutic_window_met"]:
            failures.append("nuclease resistance (half-life below 24h)")
        if not lna_stability["geometry_maintained"]:
            failures.append("LNA conformation (geometry compromised)")
        executive_summary = (
            f"MRL-ASO-001 FAILS {len(failures)} stability test(s): "
            f"{', '.join(failures)}. Further optimization may be required."
        )

    # Conclusao geral
    if all_pass:
        overall_conclusion = (
            f"MRL-ASO-001 demonstrates robust stability across all four analyses:\n\n"
            f"1. **pH stability:** dG decreases from {ASO_DG_NEUTRAL:.1f} kcal/mol (pH 7.4) "
            f"to {phagolysosome_dg:.2f} kcal/mol (pH 4.5), a loss of only "
            f"{abs(phagolysosome_dg - ASO_DG_NEUTRAL):.2f} kcal/mol. "
            f"The duplex remains far below the functional threshold "
            f"({DG_FUNCTIONAL_THRESHOLD} kcal/mol).\n\n"
            f"2. **Nuclease resistance:** The LNA-DNA-LNA gapmer with PS backbone "
            f"has an estimated half-life of {gapmer_half_life:.2f} hours at pH 4.5, "
            f"far exceeding the 24-hour therapeutic window requirement.\n\n"
            f"3. **LNA conformation:** The methylene bridge maintains C3'-endo sugar pucker "
            f"at {lna_stability['c3_endo_fraction_pH_4_5']*100:.1f}% even at pH 4.5, "
            f"ensuring optimal hybridization geometry.\n\n"
            f"4. **FDA comparison:** MRL-ASO-001's chemistry (LNA gapmer, full PS) "
            f"matches or exceeds the stability profiles of clinically approved ASOs "
            f"like mipomersen that operate in similar acidic compartments.\n\n"
            f"**The ASO is predicted to remain stable and functional in the "
            f"phagolysosomal environment of canine macrophages.**"
        )
    else:
        overall_conclusion = (
            f"MRL-ASO-001 does not pass all stability criteria. "
            f"Failed tests: {', '.join(failures)}. "
            f"Consider alternative modifications to improve stability at pH 4.5."
        )

    # Montar envelope final
    results: dict[str, Any] = {
        "module": MODULE_NAME,
        "version": MODULE_VERSION,
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": elapsed,
        "status": "success" if all_pass else "partial_failure",
        "all_tests_passed": all_pass,
        "executive_summary": executive_summary,
        "overall_conclusion": overall_conclusion,
        "aso_specifications": {
            "name": "MRL-ASO-001",
            "length_nt": ASO_LENGTH,
            "gapmer_design": f"{LNA_5PRIME}-{DNA_GAP}-{LNA_3PRIME} (LNA-DNA-LNA)",
            "backbone": "full phosphorothioate (PS)",
            "dg_neutral_kcal": ASO_DG_NEUTRAL,
            "tm_dna_celsius": ASO_TM_DNA,
            "tm_lna_celsius": ASO_TM_LNA,
        },
        "ph_stability_profile": ph_stability,
        "nuclease_resistance": nuclease_resistance,
        "lna_stability": lna_stability,
        "fda_comparison": fda_comparison,
    }

    # --- Gravar resultados ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = RESULTS_DIR / "module_a_stability.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    logger.info("Resultados JSON gravados em: %s", json_path)

    # Markdown
    report_md = generate_markdown_report(results)
    md_path = RESULTS_DIR / "module_a_stability_report.md"
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Relatorio Markdown gravado em: %s", md_path)

    # --- Sumario final ---
    logger.info("")
    logger.info("=" * 70)
    logger.info("RESULTADO FINAL: %s", "APROVADO" if all_pass else "REPROVADO")
    logger.info("=" * 70)
    logger.info("  dG no fagolisossomo (pH 4.5): %.2f kcal/mol", phagolysosome_dg)
    logger.info("  Fracao ligada: %.4f", phagolysosome_fbound)
    logger.info("  Meia-vida do gapmer: %.2f horas", gapmer_half_life)
    logger.info("  LNA C3'-endo em pH 4.5: %.1f%%", lna_stability["c3_endo_fraction_pH_4_5"] * 100)
    logger.info("  Tempo de execucao: %.2f segundos", elapsed)

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
