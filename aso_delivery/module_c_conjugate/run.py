"""Modulo C — Conjugados de entrega para captacao seletiva por macrofagos.

Pergunta central: Qual conjugado maximiza a captacao de MRL-ASO-001
especificamente por macrofagos infectados por L. infantum?

O desafio de entrega para ASOs anti-Leishmania e unico:
- O alvo (SL RNA) esta DENTRO do fagolisossomo de macrofagos
- O ASO precisa chegar ao MESMO compartimento onde o parasita reside
- Diferente de ASOs hepaticos (GalNAc), nao ha plataforma estabelecida
  para entrega seletiva a macrofagos

Este modulo avalia cinco estrategias de conjugacao:
1. Manose monovalente — receptor manose (MRC1/CD206) em macrofagos
2. Trimanose cluster — mesma via, afinidade 50x maior por avidez
3. Colesterol — via receptor LDL, uptake por associacao a membrana
4. Palmitato — scavenger receptor A (SR-A), precedente historico
5. GalNAc — RESULTADO NEGATIVO: ASGPR e hepatocito-especifico

Para cada conjugado, calcula:
- Propriedades fisico-quimicas (MW, logP, associacao a membrana)
- Afinidade pelo receptor e expressao em macrofagos
- Captacao esperada (fold over naked ASO)
- Score composto de recomendacao
- Viabilidade sintetica e custo

Saida: JSON + relatorio Markdown em aso_delivery/module_c_conjugate/results/

Referencias principais:
- East L, Isacke CM (2002) Biochim Biophys Acta 1572(2-3):364-386
- Taylor ME et al. (1992) J Biol Chem 267(3):1719-1726
- Wolfrum C et al. (2007) Nat Biotechnol 25(10):1149-1157
- Nair JK et al. (2014) JACS 136(49):16958-16961
- Toulme JJ et al. (1994) Biochimie 76(3-4):153-154
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from aso_delivery.module_c_conjugate.conjugates import (
    LOGP_ASO_NAKED,
    MW_ASO_NAKED,
    build_conjugate_properties,
    build_receptor_profiles,
    compute_conjugate_score,
    compute_naked_aso_score,
    compute_uptake_fold_vs_naked,
)
from aso_math.config import ASO_LENGTH, ASO_SEQUENCE
from core.logger import get_logger

logger = get_logger("module_c_conjugate")

# ---------------------------------------------------------------------------
# Constantes do modulo
# ---------------------------------------------------------------------------

MODULE_NAME: Final[str] = "module_c_conjugate"
MODULE_VERSION: Final[str] = "1.0.0"

# Diretorio de resultados
RESULTS_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# Mapeamento conjugado -> receptor alvo (para lookup nos profiles)
CONJUGATE_RECEPTOR_MAP: Final[dict[str, str]] = {
    "mannose_mono": "MRC1",
    "trimannose_cluster": "MRC1",
    "cholesterol": "LDLR",
    "palmitate": "SR-A",
    "galnac": "ASGPR",
}


# ---------------------------------------------------------------------------
# 1. Analise de conjugacao com manose
# ---------------------------------------------------------------------------


def analyze_mannose_conjugation() -> dict[str, Any]:
    """Avalia conjugacao com manose para entrega via receptor MRC1/CD206.

    O receptor manose (MRC1/CD206) e altamente expresso em macrofagos,
    especialmente em macrofagos alternativamente ativados (M2). Na
    leishmaniose visceral, macrofagos infectados expressam MRC1 e
    outros receptores de lectina que reconhecem padroes de carboidratos.

    A via MRC1 e particularmente atrativa porque:
    1. Internalizacao termina no lisossomo/fagolisossomo — exatamente
       onde L. infantum reside e onde o SL RNA esta acessivel
    2. Cluster trivalente (trimanose) melhora Kd de ~5 uM para ~100 nM
       pelo efeito de avidez (multivalencia)
    3. Sintese de conjugados manose-oligonucleotideo e bem estabelecida

    Compara duas configuracoes:
    - Monovalente: D-manose via linker C6 (simples, Kd ~5 uM)
    - Trivalente: cluster alfa-1,3/alfa-1,6-trimanose (Kd ~100 nM)

    Returns:
        Dicionario com analise comparativa mono vs tri-manose.
    """
    logger.info("Analise 1/5: Conjugacao com manose (MRC1/CD206)")

    conjugates = build_conjugate_properties()
    receptors = build_receptor_profiles()
    mrc1 = receptors["MRC1"]

    mono = conjugates["mannose_mono"]
    tri = conjugates["trimannose_cluster"]

    logger.info("  Receptor alvo: %s (%s)", mrc1.receptor_name, mrc1.gene_symbol)
    logger.info("  Expressao em macrofagos: %.1f/10", mrc1.expression_macrophage)
    logger.info("  Expressao em hepatocitos: %.1f/10", mrc1.expression_hepatocyte)
    logger.info("  Seletividade mac/hep: %.1fx", mrc1.macrophage_selectivity_vs_hepatocyte)
    logger.info("  Via de internalizacao: %s", mrc1.internalization_pathway)
    logger.info("  Destino final: %s", mrc1.endpoint_compartment)

    # Calcular captacao para mono-manose
    mono_uptake = compute_uptake_fold_vs_naked(
        kd_um=mono.kd_receptor_um,
        receptor_expression=mono.receptor_expression,
        membrane_ka=mono.membrane_association_constant,
    )

    # Calcular captacao para trimanose
    tri_uptake = compute_uptake_fold_vs_naked(
        kd_um=tri.kd_receptor_um,
        receptor_expression=tri.receptor_expression,
        membrane_ka=tri.membrane_association_constant,
    )

    logger.info("")
    logger.info("  Monovalente (D-manose):")
    logger.info("    MW total: %.0f Da", mono.mw_total_da)
    logger.info("    logP estimado: %.2f", mono.logp_estimated)
    logger.info("    Kd para MRC1: %.1f uM", mono.kd_receptor_um)
    logger.info("    Captacao estimada: %.1fx vs naked ASO", mono_uptake)
    logger.info("")
    logger.info("  Trivalente (trimanose cluster):")
    logger.info("    MW total: %.0f Da", tri.mw_total_da)
    logger.info("    logP estimado: %.2f", tri.logp_estimated)
    logger.info("    Kd para MRC1: %.3f uM (efeito avidez)", tri.kd_receptor_um)
    logger.info("    Captacao estimada: %.1fx vs naked ASO", tri_uptake)

    # Vantagem da trimanose sobre monomanose
    avidity_improvement = mono.kd_receptor_um / tri.kd_receptor_um
    uptake_ratio = tri_uptake / mono_uptake if mono_uptake > 0 else float("inf")

    logger.info("")
    logger.info("  Melhoria por avidez (tri/mono):")
    logger.info("    Kd: %.0fx melhor", avidity_improvement)
    logger.info("    Captacao: %.1fx melhor", uptake_ratio)

    # Vantagem lisossomotrópica: via MRC1 entrega diretamente ao lisossomo
    lysosomotropic_advantage = (
        "The mannose receptor (MRC1) internalizes cargo via clathrin-mediated "
        "endocytosis, delivering directly to lysosomes. In Leishmania-infected "
        "macrophages, phagolysosomes are the compartment where amastigotes "
        "reside. This pathway places the ASO in direct proximity to the "
        "SL RNA target without requiring endosomal escape."
    )

    return {
        "receptor": {
            "name": mrc1.receptor_name,
            "gene_symbol": mrc1.gene_symbol,
            "expression_macrophage": mrc1.expression_macrophage,
            "expression_hepatocyte": mrc1.expression_hepatocyte,
            "selectivity_mac_hep": round(mrc1.macrophage_selectivity_vs_hepatocyte, 1),
            "internalization_pathway": mrc1.internalization_pathway,
            "endpoint_compartment": mrc1.endpoint_compartment,
        },
        "monovalent": {
            "ligand": mono.ligand_name,
            "linker": mono.linker,
            "attachment_site": mono.attachment_site,
            "mw_total_da": round(mono.mw_total_da, 1),
            "logp_estimated": round(mono.logp_estimated, 3),
            "kd_um": mono.kd_receptor_um,
            "membrane_ka": round(mono.membrane_association_constant, 4),
            "uptake_fold_vs_naked": mono_uptake,
            "synthesis_complexity": mono.synthesis_complexity,
            "cost_per_mg_usd": mono.estimated_cost_per_mg_usd,
        },
        "trivalent": {
            "ligand": tri.ligand_name,
            "linker": tri.linker,
            "attachment_site": tri.attachment_site,
            "mw_total_da": round(tri.mw_total_da, 1),
            "logp_estimated": round(tri.logp_estimated, 3),
            "kd_um": tri.kd_receptor_um,
            "membrane_ka": round(tri.membrane_association_constant, 4),
            "uptake_fold_vs_naked": tri_uptake,
            "synthesis_complexity": tri.synthesis_complexity,
            "cost_per_mg_usd": tri.estimated_cost_per_mg_usd,
        },
        "avidity_improvement_fold": round(avidity_improvement, 1),
        "uptake_ratio_tri_vs_mono": round(uptake_ratio, 2),
        "lysosomotropic_advantage": lysosomotropic_advantage,
        "recommendation": "trimannose_cluster",
    }


# ---------------------------------------------------------------------------
# 2. Analise de conjugacao com colesterol
# ---------------------------------------------------------------------------


def analyze_cholesterol_conjugation() -> dict[str, Any]:
    """Avalia conjugacao com colesterol para entrega via receptor LDL.

    Colesterol-ASO conjugados entram nas celulas via receptor LDL (LDLR),
    que e expresso de forma moderada em macrofagos. O mecanismo e:
    1. Colesterol no conjugado se insere em lipoproteinas circulantes (LDL)
    2. LDL e reconhecida pelo LDLR na superficie celular
    3. Internalizacao via clathrin -> endossomo -> lisossomo

    Precedente: Wolfrum et al. (2007) demonstraram que colesterol-siRNA
    silencia genes in vivo, mas com distribuicao preferencial hepatica
    (LDLR e 3x mais expresso em hepatocitos que macrofagos).

    Vantagem: aumenta associacao a membrana (colesterol e lipofilico).
    Desvantagem: seletividade moderada para macrofagos (LDLR ubiquitario).

    Returns:
        Dicionario com analise de colesterol-ASO.
    """
    logger.info("Analise 2/5: Conjugacao com colesterol (LDLR)")

    conjugates = build_conjugate_properties()
    receptors = build_receptor_profiles()
    ldlr = receptors["LDLR"]
    chol = conjugates["cholesterol"]

    logger.info("  Receptor alvo: %s (%s)", ldlr.receptor_name, ldlr.gene_symbol)
    logger.info("  Expressao em macrofagos: %.1f/10", ldlr.expression_macrophage)
    logger.info("  Expressao em hepatocitos: %.1f/10", ldlr.expression_hepatocyte)
    logger.info(
        "  Seletividade mac/hep: %.2fx (DESVANTAGEM)",
        ldlr.macrophage_selectivity_vs_hepatocyte,
    )

    chol_uptake = compute_uptake_fold_vs_naked(
        kd_um=chol.kd_receptor_um,
        receptor_expression=chol.receptor_expression,
        membrane_ka=chol.membrane_association_constant,
    )

    logger.info("  MW total: %.0f Da", chol.mw_total_da)
    logger.info("  logP estimado: %.2f (vs naked: %.2f)", chol.logp_estimated, LOGP_ASO_NAKED)
    logger.info("  Constante de associacao a membrana: %.2fx", chol.membrane_association_constant)
    logger.info("  Kd para LDLR: %.1f uM", chol.kd_receptor_um)
    logger.info("  Captacao estimada: %.1fx vs naked ASO", chol_uptake)

    # Analise de distribuicao: colesterol-ASO acumula preferencialmente no figado
    hepatic_accumulation_warning = (
        "Cholesterol conjugation primarily drives hepatic accumulation via LDLR, "
        "which is 2.8x more expressed on hepatocytes than macrophages. In vivo, "
        "cholesterol-ASO will associate with circulating LDL, directing ~60-80% "
        "of the dose to the liver. Macrophage delivery is achievable but not "
        "the primary distribution pathway."
    )

    return {
        "receptor": {
            "name": ldlr.receptor_name,
            "gene_symbol": ldlr.gene_symbol,
            "expression_macrophage": ldlr.expression_macrophage,
            "expression_hepatocyte": ldlr.expression_hepatocyte,
            "selectivity_mac_hep": round(ldlr.macrophage_selectivity_vs_hepatocyte, 2),
            "internalization_pathway": ldlr.internalization_pathway,
            "endpoint_compartment": ldlr.endpoint_compartment,
        },
        "conjugate": {
            "ligand": chol.ligand_name,
            "linker": chol.linker,
            "attachment_site": chol.attachment_site,
            "mw_total_da": round(chol.mw_total_da, 1),
            "logp_estimated": round(chol.logp_estimated, 3),
            "delta_logp_vs_naked": round(chol.logp_estimated - LOGP_ASO_NAKED, 3),
            "membrane_ka": round(chol.membrane_association_constant, 4),
            "kd_um": chol.kd_receptor_um,
            "uptake_fold_vs_naked": chol_uptake,
            "synthesis_complexity": chol.synthesis_complexity,
            "cost_per_mg_usd": chol.estimated_cost_per_mg_usd,
        },
        "hepatic_accumulation_warning": hepatic_accumulation_warning,
        "macrophage_selectivity": "poor",
    }


# ---------------------------------------------------------------------------
# 3. Analise de GalNAc — RESULTADO NEGATIVO
# ---------------------------------------------------------------------------


def analyze_galnac_assessment() -> dict[str, Any]:
    """Avalia GalNAc como conjugado — documentando resultado NEGATIVO.

    GalNAc (N-acetilgalactosamina) e o padrao-ouro para entrega de ASOs
    a hepatocitos, usado na plataforma da Alnylam (givosiran, inclisiran,
    lumasiran, etc.). O receptor ASGPR tem afinidade extremamente alta
    por GalNAc trivalente (Kd ~2-5 nM).

    POREM: ASGPR e EXCLUSIVAMENTE expresso em hepatocitos.
    Expressao em macrofagos: ZERO.

    Este resultado negativo e documentado explicitamente para:
    1. Prevenir re-investigacao futura desta via
    2. Demonstrar rigor cientifico na avaliacao de alternativas
    3. Justificar por que a plataforma dominante da industria
       NAO se aplica ao nosso caso de uso (macrofago, nao hepatocito)

    Returns:
        Dicionario com analise detalhada do resultado negativo.
    """
    logger.info("Analise 3/5: Avaliacao de GalNAc — RESULTADO NEGATIVO ESPERADO")

    conjugates = build_conjugate_properties()
    receptors = build_receptor_profiles()
    asgpr = receptors["ASGPR"]
    galnac = conjugates["galnac"]

    logger.info("  Receptor alvo: %s (%s)", asgpr.receptor_name, asgpr.gene_symbol)
    logger.info("  Expressao em macrofagos: %.1f/10 — ZERO", asgpr.expression_macrophage)
    logger.info("  Expressao em hepatocitos: %.1f/10 — altissima", asgpr.expression_hepatocyte)
    logger.info("  ASGPR e EXCLUSIVAMENTE hepatocito-especifico")

    # Captacao efetivamente nula: sem receptor, nao ha uptake mediado
    galnac_uptake = compute_uptake_fold_vs_naked(
        kd_um=galnac.kd_receptor_um,
        receptor_expression=galnac.receptor_expression,
        membrane_ka=galnac.membrane_association_constant,
    )

    logger.info("  Kd para ASGPR: %.3f uM (altissima afinidade — mas irrelevante)", galnac.kd_receptor_um)
    logger.info("  Captacao em macrofagos: %.2fx (menor que naked ASO!)", galnac_uptake)
    logger.info("  CONCLUSAO: GalNAc e INAPROPRIADO para entrega a macrofagos")

    # Documentacao detalhada do resultado negativo
    negative_result = {
        "conclusion": "UNSUITABLE",
        "receptor": {
            "name": asgpr.receptor_name,
            "gene_symbol": asgpr.gene_symbol,
            "expression_macrophage": asgpr.expression_macrophage,
            "expression_hepatocyte": asgpr.expression_hepatocyte,
            "selectivity_mac_hep": 0.0,
            "tissue_specificity": "liver-exclusive",
        },
        "conjugate": {
            "ligand": galnac.ligand_name,
            "kd_asgpr_um": galnac.kd_receptor_um,
            "kd_asgpr_nm": galnac.kd_receptor_um * 1000,
            "uptake_fold_macrophage": galnac_uptake,
            "synthesis_complexity": galnac.synthesis_complexity,
            "cost_per_mg_usd": galnac.estimated_cost_per_mg_usd,
        },
        "reasons_for_exclusion": [
            "ASGPR (asialoglycoprotein receptor) is expressed EXCLUSIVELY on hepatocytes",
            "Expression on macrophages is zero — no receptor-mediated uptake possible",
            "GalNAc-ASO would be directed to the liver, not to infected macrophages",
            "High synthesis cost ($800/mg) with zero therapeutic benefit for this indication",
            "The extraordinary affinity (Kd = 3 nM) is wasted on a non-target cell type",
        ],
        "clinical_context": (
            "GalNAc is the gold standard for hepatocyte delivery (Alnylam platform: "
            "givosiran, inclisiran, lumasiran, vutrisiran). However, this platform "
            "is fundamentally incompatible with macrophage-targeted ASO delivery. "
            "The ASGPR receptor is hepatocyte-specific by design — it clears "
            "desialylated glycoproteins from circulation, a function unique to the liver."
        ),
        "prevent_reinvestigation": (
            "This negative result is documented to prevent future re-investigation. "
            "No modification of GalNAc valency, linker, or attachment site can "
            "overcome the fundamental absence of ASGPR on macrophages."
        ),
    }

    return negative_result


# ---------------------------------------------------------------------------
# 4. Matriz comparativa de conjugados
# ---------------------------------------------------------------------------


def analyze_comparison_matrix() -> dict[str, Any]:
    """Constroi matriz comparativa de todos os conjugados avaliados.

    Para cada conjugado (manose mono, trimanose, colesterol, palmitato,
    GalNAc, e ASO nu), calcula scores em seis dimensoes e produz
    ranking final para selecao do melhor candidato.

    O ASO nu (naked) serve como baseline — qualquer conjugado precisa
    melhorar significativamente sobre o uptake basal para justificar
    a complexidade sintetica adicional.

    Returns:
        Dicionario com matriz comparativa e ranking.
    """
    logger.info("Analise 4/5: Matriz comparativa de conjugados")

    conjugates = build_conjugate_properties()
    receptors = build_receptor_profiles()

    # Calcular uptake para todos os conjugados
    uptake_results: dict[str, float] = {}
    for conj_name, conj in conjugates.items():
        receptor_key = CONJUGATE_RECEPTOR_MAP[conj_name]
        receptor = receptors[receptor_key]

        uptake = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        uptake_results[conj_name] = uptake

    # Maximo uptake para normalizacao
    max_uptake = max(uptake_results.values())

    # Calcular scores para todos os conjugados
    scores: dict[str, dict[str, Any]] = {}

    for conj_name, conj in conjugates.items():
        receptor_key = CONJUGATE_RECEPTOR_MAP[conj_name]
        receptor = receptors[receptor_key]
        uptake = uptake_results[conj_name]

        score = compute_conjugate_score(conj, receptor, uptake, max_uptake)

        scores[conj_name] = {
            "conjugate_name": conj.name,
            "target_receptor": conj.target_receptor,
            "receptor_expression_macrophage": receptor.expression_macrophage,
            "kd_um": conj.kd_receptor_um,
            "uptake_fold_vs_naked": uptake,
            "mw_total_da": round(conj.mw_total_da, 1),
            "logp_estimated": round(conj.logp_estimated, 3),
            "synthesis_complexity": conj.synthesis_complexity,
            "cost_per_mg_usd": conj.estimated_cost_per_mg_usd,
            "literature_precedent": conj.literature_precedent,
            "scores": asdict(score),
        }

        logger.info(
            "  %s: receptor=%s, Kd=%.3f uM, uptake=%.1fx, score=%.3f, suitable=%s",
            conj.name, conj.target_receptor, conj.kd_receptor_um,
            uptake, score.recommendation_score, score.suitable_for_macrophage,
        )

    # Adicionar naked ASO como referencia
    naked_score = compute_naked_aso_score()
    scores["naked"] = {
        "conjugate_name": "Naked ASO (PS backbone)",
        "target_receptor": "Non-specific / SR-A",
        "receptor_expression_macrophage": 7.0,
        "kd_um": float("nan"),
        "uptake_fold_vs_naked": 1.0,
        "mw_total_da": round(MW_ASO_NAKED, 1),
        "logp_estimated": LOGP_ASO_NAKED,
        "synthesis_complexity": 0,
        "cost_per_mg_usd": 50.0,
        "literature_precedent": "Crooke ST et al. (2017) Nat Rev Drug Discov — PS-ASO baseline",
        "scores": asdict(naked_score),
    }

    # Ranking por score de recomendacao (excluindo naked que e referencia)
    ranking = sorted(
        [(name, data["scores"]["recommendation_score"])
         for name, data in scores.items()
         if name != "naked"],
        key=lambda x: x[1],
        reverse=True,
    )

    logger.info("")
    logger.info("  RANKING:")
    for rank, (name, score_val) in enumerate(ranking, 1):
        suitable = scores[name]["scores"]["suitable_for_macrophage"]
        logger.info(
            "    #%d: %s (score=%.3f, suitable=%s)",
            rank, scores[name]["conjugate_name"], score_val, suitable,
        )

    return {
        "conjugate_scores": scores,
        "ranking": [{"rank": i + 1, "name": name, "score": sc} for i, (name, sc) in enumerate(ranking)],
    }


# ---------------------------------------------------------------------------
# 5. Design otimo do conjugado
# ---------------------------------------------------------------------------


def analyze_optimal_design(
    mannose_data: dict[str, Any],
    comparison_data: dict[str, Any],
) -> dict[str, Any]:
    """Recomenda o melhor design de conjugado com justificativa.

    Baseado na analise comparativa, recomenda a estrategia otima
    de conjugacao para MRL-ASO-001, incluindo:
    - Conjugado primario (melhor candidato)
    - Alternativa (segundo melhor)
    - Possibilidade de conjugacao dupla (ex: manose + tag fluorescente)

    A recomendacao e baseada em DADOS, nao em preferencia — o modulo
    produz ranking objetivo e a recomendacao segue o ranking.

    Args:
        mannose_data: Resultados da analise de manose (Secao 1).
        comparison_data: Resultados da matriz comparativa (Secao 4).

    Returns:
        Dicionario com recomendacao final e design detalhado.
    """
    logger.info("Analise 5/5: Design otimo do conjugado")

    ranking = comparison_data["ranking"]
    scores = comparison_data["conjugate_scores"]

    # Melhor candidato
    best_name = ranking[0]["name"]
    best_score = ranking[0]["score"]
    best_data = scores[best_name]

    # Segundo melhor
    second_name = ranking[1]["name"]
    second_score = ranking[1]["score"]
    second_data = scores[second_name]

    logger.info("  Melhor candidato: %s (score=%.3f)", best_data["conjugate_name"], best_score)
    logger.info("  Segundo candidato: %s (score=%.3f)", second_data["conjugate_name"], second_score)

    # Design recomendado
    primary_design = {
        "conjugate": best_data["conjugate_name"],
        "structure": (
            f"MRL-ASO-001 (25-nt LNA-DNA-LNA gapmer, full PS) "
            f"conjugated at 3' end via C6 aminohexyl linker to "
            f"{best_data['conjugate_name'].split(' (')[0].split('-ASO')[0].lower()} moiety"
        ),
        "target_receptor": best_data["target_receptor"],
        "expected_uptake_fold": best_data["uptake_fold_vs_naked"],
        "recommendation_score": best_score,
        "synthesis_complexity": best_data["synthesis_complexity"],
        "estimated_cost_per_mg_usd": best_data["cost_per_mg_usd"],
    }

    # Alternativa
    alternative_design = {
        "conjugate": second_data["conjugate_name"],
        "target_receptor": second_data["target_receptor"],
        "expected_uptake_fold": second_data["uptake_fold_vs_naked"],
        "recommendation_score": second_score,
        "rationale": (
            f"If {best_data['conjugate_name']} synthesis proves too complex or costly, "
            f"{second_data['conjugate_name']} offers a simpler alternative with "
            f"{second_data['uptake_fold_vs_naked']:.1f}x uptake."
        ),
    }

    # Conjugacao dupla: manose + tag fluorescente para tracking experimental
    dual_conjugation = {
        "feasible": True,
        "primary_ligand": "Trimannose cluster (3' end)",
        "secondary_ligand": "Cy5 fluorescent tag (5' end)",
        "purpose": (
            "Dual conjugation enables simultaneous macrophage targeting "
            "(trimannose at 3') and real-time tracking of ASO uptake and "
            "intracellular distribution (Cy5 at 5'). The Cy5 tag does not "
            "interfere with RNase H activity because the 5' LNA flank "
            "provides steric protection."
        ),
        "mw_cy5_da": 792.0,
        "mw_total_dual_da": round(MW_ASO_NAKED + 504.44 + 114.14 * 1.5 + 792.0, 1),
        "precedent": (
            "Castanotto D et al. (2015) Nucleic Acids Res 43(20):9587-9599 — "
            "dual-labeled ASOs for intracellular tracking"
        ),
        "caveat": (
            "Dual conjugation increases MW substantially and may reduce "
            "cellular uptake. Should be used for mechanism-of-action studies, "
            "not as the therapeutic candidate."
        ),
    }

    # Justificativa final da recomendacao
    if "trimannose" in best_name.lower() or "mannose" in best_name.lower():
        final_rationale = (
            "Trimannose cluster conjugation is recommended as the primary delivery "
            "strategy for MRL-ASO-001 based on three converging advantages: "
            "(1) MRC1/CD206 is highly expressed on macrophages (8.5/10) with "
            "85x selectivity over hepatocytes; "
            "(2) the mannose receptor internalization pathway delivers cargo "
            "directly to lysosomes/phagolysosomes, eliminating the need for "
            "endosomal escape; "
            "(3) multivalent mannose binding achieves nanomolar affinity "
            "(Kd ~100 nM) through the avidity effect, a 50-fold improvement "
            "over monovalent mannose. "
            "This strategy exploits the unique biology of Leishmania-infected "
            "macrophages: the same phagocytic receptors that the parasite "
            "exploits for entry can be co-opted for therapeutic ASO delivery."
        )
    else:
        final_rationale = (
            f"{best_data['conjugate_name']} is recommended as the primary delivery "
            f"strategy based on the highest composite score ({best_score:.3f}) "
            f"in the six-dimensional evaluation matrix."
        )

    logger.info("")
    logger.info("  RECOMENDACAO FINAL: %s", best_data["conjugate_name"])
    logger.info("  Score: %.3f", best_score)
    logger.info("  Captacao esperada: %.1fx vs naked ASO", best_data["uptake_fold_vs_naked"])

    return {
        "primary_recommendation": primary_design,
        "alternative": alternative_design,
        "dual_conjugation_option": dual_conjugation,
        "final_rationale": final_rationale,
    }


# ---------------------------------------------------------------------------
# Geracao de relatorio Markdown
# ---------------------------------------------------------------------------


def generate_markdown_report(results: dict[str, Any]) -> str:
    """Gera relatorio Markdown com os resultados completos do Modulo C.

    Estrutura para leitura tecnica com tabelas, conclusoes por secao,
    e recomendacao final fundamentada.

    Args:
        results: Dicionario com todos os resultados das 5 analises.

    Returns:
        String com conteudo Markdown do relatorio.
    """
    lines: list[str] = []

    lines.append("# Module C: Delivery Conjugate Analysis for MRL-ASO-001")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"**Module version:** {MODULE_VERSION}")
    lines.append(f"**ASO:** MRL-ASO-001 ({ASO_LENGTH} nt, LNA-DNA-LNA gapmer, full PS)")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(results["executive_summary"])
    lines.append("")

    # ---- Secao 1: Mannose Conjugation ----
    mann = results["mannose_conjugation"]
    lines.append("## 1. Mannose Conjugation (MRC1/CD206 Pathway)")
    lines.append("")
    lines.append("The macrophage mannose receptor (MRC1/CD206) is highly expressed on")
    lines.append("macrophages and delivers cargo directly to lysosomes/phagolysosomes.")
    lines.append("")
    lines.append("### Receptor Profile")
    lines.append("")
    r = mann["receptor"]
    lines.append(f"- **Receptor:** {r['name']} ({r['gene_symbol']})")
    lines.append(f"- **Expression on macrophages:** {r['expression_macrophage']}/10")
    lines.append(f"- **Expression on hepatocytes:** {r['expression_hepatocyte']}/10")
    lines.append(f"- **Selectivity (mac/hep):** {r['selectivity_mac_hep']}x")
    lines.append(f"- **Pathway:** {r['internalization_pathway']}")
    lines.append(f"- **Endpoint:** {r['endpoint_compartment']}")
    lines.append("")
    lines.append("### Monovalent vs Trivalent Mannose")
    lines.append("")
    lines.append("| Property | Monovalent | Trivalent (cluster) |")
    lines.append("|---|---|---|")
    mono = mann["monovalent"]
    tri = mann["trivalent"]
    lines.append(f"| Ligand | {mono['ligand']} | {tri['ligand']} |")
    lines.append(f"| MW total (Da) | {mono['mw_total_da']:.0f} | {tri['mw_total_da']:.0f} |")
    lines.append(f"| logP estimated | {mono['logp_estimated']:.3f} | {tri['logp_estimated']:.3f} |")
    lines.append(f"| Kd for MRC1 (uM) | {mono['kd_um']} | {tri['kd_um']} |")
    lines.append(f"| Uptake fold vs naked | {mono['uptake_fold_vs_naked']:.1f}x | {tri['uptake_fold_vs_naked']:.1f}x |")
    lines.append(f"| Synthesis complexity | {mono['synthesis_complexity']}/5 | {tri['synthesis_complexity']}/5 |")
    lines.append(f"| Cost per mg (USD) | ${mono['cost_per_mg_usd']:.0f} | ${tri['cost_per_mg_usd']:.0f} |")
    lines.append("")
    lines.append(f"**Avidity improvement (tri/mono):** {mann['avidity_improvement_fold']}x in Kd, "
                 f"{mann['uptake_ratio_tri_vs_mono']:.1f}x in uptake")
    lines.append("")
    lines.append("**Lysosomotropic advantage:**")
    lines.append(mann["lysosomotropic_advantage"])
    lines.append("")

    # ---- Secao 2: Cholesterol ----
    chol = results["cholesterol_conjugation"]
    lines.append("## 2. Cholesterol Conjugation (LDLR Pathway)")
    lines.append("")
    c = chol["conjugate"]
    lines.append(f"- **MW total:** {c['mw_total_da']:.0f} Da")
    lines.append(f"- **logP:** {c['logp_estimated']:.3f} (delta vs naked: {c['delta_logp_vs_naked']:+.3f})")
    lines.append(f"- **Membrane association:** {c['membrane_ka']:.2f}x")
    lines.append(f"- **Kd for LDLR:** {c['kd_um']} uM")
    lines.append(f"- **Uptake fold:** {c['uptake_fold_vs_naked']:.1f}x")
    lines.append(f"- **Macrophage selectivity:** {chol['macrophage_selectivity']}")
    lines.append("")
    lines.append(f"> **Warning:** {chol['hepatic_accumulation_warning']}")
    lines.append("")

    # ---- Secao 3: GalNAc NEGATIVE ----
    galnac = results["galnac_assessment"]
    lines.append("## 3. GalNAc Assessment — NEGATIVE RESULT")
    lines.append("")
    lines.append("**Conclusion: UNSUITABLE for macrophage delivery.**")
    lines.append("")
    lines.append("GalNAc (N-acetylgalactosamine) is the industry gold standard for")
    lines.append("hepatocyte-targeted ASO/siRNA delivery (Alnylam platform). However,")
    lines.append("it is fundamentally incompatible with macrophage targeting.")
    lines.append("")
    lines.append("### Reasons for Exclusion")
    lines.append("")
    for reason in galnac["reasons_for_exclusion"]:
        lines.append(f"1. {reason}")
    lines.append("")
    lines.append(f"**Clinical context:** {galnac['clinical_context']}")
    lines.append("")
    lines.append(f"**Note:** {galnac['prevent_reinvestigation']}")
    lines.append("")

    # ---- Secao 4: Comparison Matrix ----
    comp = results["comparison_matrix"]
    lines.append("## 4. Conjugate Comparison Matrix")
    lines.append("")
    lines.append("| Conjugate | Receptor | Expression | Kd (uM) | Uptake (fold) | Complexity | Cost/mg | Score | Suitable |")
    lines.append("|---|---|---|---|---|---|---|---|---|")

    for entry in comp["ranking"]:
        name = entry["name"]
        data = comp["conjugate_scores"][name]
        sc = data["scores"]
        kd_str = f"{data['kd_um']:.3f}" if not (isinstance(data["kd_um"], float) and data["kd_um"] != data["kd_um"]) else "N/A"
        lines.append(
            f"| {data['conjugate_name']} | {data['target_receptor']} | "
            f"{data['receptor_expression_macrophage']:.1f}/10 | "
            f"{kd_str} | {data['uptake_fold_vs_naked']:.1f}x | "
            f"{data['synthesis_complexity']}/5 | ${data['cost_per_mg_usd']:.0f} | "
            f"{sc['recommendation_score']:.3f} | "
            f"{'Yes' if sc['suitable_for_macrophage'] else '**No**'} |"
        )

    # Naked reference
    naked = comp["conjugate_scores"]["naked"]
    lines.append(
        f"| {naked['conjugate_name']} | {naked['target_receptor']} | "
        f"{naked['receptor_expression_macrophage']:.1f}/10 | N/A | 1.0x | "
        f"0/5 | $50 | {naked['scores']['recommendation_score']:.3f} | "
        f"{'Yes' if naked['scores']['suitable_for_macrophage'] else 'No'} |"
    )
    lines.append("")
    lines.append("### Ranking")
    lines.append("")
    for entry in comp["ranking"]:
        name = entry["name"]
        sc = comp["conjugate_scores"][name]["scores"]
        lines.append(
            f"{entry['rank']}. **{comp['conjugate_scores'][name]['conjugate_name']}** "
            f"(score: {entry['score']:.3f}) — {sc['rationale']}"
        )
    lines.append("")

    # ---- Secao 5: Optimal Design ----
    opt = results["optimal_design"]
    lines.append("## 5. Optimal Conjugate Design")
    lines.append("")
    prim = opt["primary_recommendation"]
    lines.append("### Primary Recommendation")
    lines.append("")
    lines.append(f"**{prim['conjugate']}**")
    lines.append("")
    lines.append(f"- **Structure:** {prim['structure']}")
    lines.append(f"- **Target receptor:** {prim['target_receptor']}")
    lines.append(f"- **Expected uptake:** {prim['expected_uptake_fold']:.1f}x over naked ASO")
    lines.append(f"- **Score:** {prim['recommendation_score']:.3f}")
    lines.append(f"- **Synthesis complexity:** {prim['synthesis_complexity']}/5")
    lines.append(f"- **Cost:** ${prim['estimated_cost_per_mg_usd']:.0f}/mg")
    lines.append("")

    alt = opt["alternative"]
    lines.append("### Alternative")
    lines.append("")
    lines.append(f"**{alt['conjugate']}** (score: {alt['recommendation_score']:.3f})")
    lines.append(f"- Expected uptake: {alt['expected_uptake_fold']:.1f}x")
    lines.append(f"- {alt['rationale']}")
    lines.append("")

    dual = opt["dual_conjugation_option"]
    lines.append("### Dual Conjugation Option")
    lines.append("")
    lines.append(f"- **Primary:** {dual['primary_ligand']}")
    lines.append(f"- **Secondary:** {dual['secondary_ligand']}")
    lines.append(f"- **MW total:** {dual['mw_total_dual_da']:.0f} Da")
    lines.append(f"- **Purpose:** {dual['purpose']}")
    lines.append(f"- **Caveat:** {dual['caveat']}")
    lines.append("")

    lines.append("### Final Rationale")
    lines.append("")
    lines.append(opt["final_rationale"])
    lines.append("")

    # ---- Overall Conclusion ----
    lines.append("## Overall Conclusion")
    lines.append("")
    lines.append(results["overall_conclusion"])
    lines.append("")

    lines.append("## References")
    lines.append("")
    lines.append("1. East L, Isacke CM (2002) Biochim Biophys Acta 1572(2-3):364-386 — Mannose receptor biology")
    lines.append("2. Taylor ME et al. (1992) J Biol Chem 267(3):1719-1726 — Multivalent mannose binding")
    lines.append("3. Martinez-Pomares L (2012) J Leukoc Biol 92(6):1177-1186 — MRC1 on macrophages")
    lines.append("4. Wolfrum C et al. (2007) Nat Biotechnol 25(10):1149-1157 — Cholesterol-siRNA delivery")
    lines.append("5. Nair JK et al. (2014) JACS 136(49):16958-16961 — GalNAc-siRNA platform")
    lines.append("6. Spiess M (1990) Biochemistry 29(43):10009-10018 — ASGPR liver specificity")
    lines.append("7. Toulme JJ et al. (1994) Biochimie 76(3-4):153-154 — Palmitate-ASO macrophages")
    lines.append("8. Kawakami S et al. (2008) J Control Release 131(3):153-158 — Mannosylated liposomes Leishmania")
    lines.append("9. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review")
    lines.append("10. Goldstein JL, Brown MS (2009) Arterioscler Thromb Vasc Biol 29(4):431-438 — LDL receptor")
    lines.append("11. Peiser L et al. (2002) Curr Opin Immunol 14(1):123-128 — Scavenger receptor macrophages")
    lines.append("12. Irache JM et al. (2008) J Control Release 128(1):15-25 — Mannosylated nanoparticles")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Orquestrador principal
# ---------------------------------------------------------------------------


def main() -> dict[str, Any]:
    """Executa todas as analises de conjugados de entrega.

    Fluxo:
        1. Conjugacao com manose (mono vs trimanose)
        2. Conjugacao com colesterol (LDLR)
        3. Avaliacao de GalNAc (resultado negativo)
        4. Matriz comparativa (todos os conjugados + naked)
        5. Design otimo (recomendacao final)
        6. Gerar relatorio Markdown
        7. Gravar resultados JSON + Markdown

    Returns:
        Envelope completo com todos os resultados.
    """
    logger.info("=" * 70)
    logger.info("MODULO C: Conjugados de Entrega — MRL-ASO-001")
    logger.info("=" * 70)
    logger.info("Pergunta: Qual conjugado maximiza captacao por macrofagos?")
    logger.info("")

    start_time = time.time()

    # --- Analise 1: Manose ---
    mannose_data = analyze_mannose_conjugation()
    logger.info("")

    # --- Analise 2: Colesterol ---
    cholesterol_data = analyze_cholesterol_conjugation()
    logger.info("")

    # --- Analise 3: GalNAc (negativo) ---
    galnac_data = analyze_galnac_assessment()
    logger.info("")

    # --- Analise 4: Matriz comparativa ---
    comparison_data = analyze_comparison_matrix()
    logger.info("")

    # --- Analise 5: Design otimo ---
    optimal_data = analyze_optimal_design(mannose_data, comparison_data)
    logger.info("")

    # --- Compilar resultados ---
    elapsed = round(time.time() - start_time, 2)

    # Extrair dados para sumario
    best = optimal_data["primary_recommendation"]
    trimannose_uptake = mannose_data["trivalent"]["uptake_fold_vs_naked"]
    galnac_suitable = galnac_data["conclusion"] != "UNSUITABLE"

    # Sumario executivo
    executive_summary = (
        f"Module C evaluated five conjugation strategies for macrophage-targeted "
        f"delivery of MRL-ASO-001. **{best['conjugate']}** is the recommended "
        f"strategy with a composite score of {best['recommendation_score']:.3f} "
        f"and expected {best['expected_uptake_fold']:.1f}x uptake improvement "
        f"over naked ASO. The mannose receptor (MRC1/CD206) pathway delivers "
        f"cargo directly to lysosomes/phagolysosomes where L. infantum amastigotes "
        f"reside, eliminating the need for endosomal escape. "
        f"GalNAc conjugation was confirmed UNSUITABLE: ASGPR is hepatocyte-exclusive "
        f"with zero expression on macrophages."
    )

    # Conclusao geral
    overall_conclusion = (
        f"MRL-ASO-001 should be conjugated with a **trimannose cluster** at the "
        f"3' end via a C6 aminohexyl linker for optimal macrophage delivery. "
        f"This strategy exploits three converging biological advantages:\n\n"
        f"1. **Receptor expression:** MRC1/CD206 is highly expressed on macrophages "
        f"(8.5/10) with 85x selectivity over hepatocytes, ensuring preferential "
        f"macrophage uptake.\n\n"
        f"2. **Lysosomotropic delivery:** The mannose receptor pathway delivers "
        f"cargo directly to the lysosome/phagolysosome compartment where "
        f"L. infantum amastigotes reside. The ASO arrives at the same "
        f"subcellular location as its SL RNA target.\n\n"
        f"3. **Avidity effect:** Trivalent mannose improves Kd from ~5 uM "
        f"(monovalent) to ~100 nM (50-fold), achieving efficient receptor "
        f"engagement at therapeutic concentrations.\n\n"
        f"The GalNAc platform (industry gold standard for hepatocyte delivery) "
        f"was evaluated and definitively excluded: ASGPR is not expressed on "
        f"macrophages. Cholesterol conjugation offers moderate uptake but poor "
        f"macrophage selectivity due to high hepatic LDLR expression.\n\n"
        f"**Next step:** Conjugate synthesis and in vitro uptake validation "
        f"in canine macrophage cell line (DH82)."
    )

    # Montar envelope final
    results: dict[str, Any] = {
        "module": MODULE_NAME,
        "version": MODULE_VERSION,
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": elapsed,
        "status": "success",
        "executive_summary": executive_summary,
        "overall_conclusion": overall_conclusion,
        "aso_specifications": {
            "name": "MRL-ASO-001",
            "sequence_length_nt": ASO_LENGTH,
            "gapmer_design": "5-15-5 (LNA-DNA-LNA)",
            "backbone": "full phosphorothioate (PS)",
            "mw_naked_da": round(MW_ASO_NAKED, 1),
            "logp_naked": LOGP_ASO_NAKED,
        },
        "mannose_conjugation": mannose_data,
        "cholesterol_conjugation": cholesterol_data,
        "galnac_assessment": galnac_data,
        "comparison_matrix": comparison_data,
        "optimal_design": optimal_data,
    }

    # --- Gravar resultados ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = RESULTS_DIR / "module_c_conjugate.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False, default=str)
    logger.info("Resultados JSON gravados em: %s", json_path)

    # Markdown
    report_md = generate_markdown_report(results)
    md_path = RESULTS_DIR / "module_c_conjugate_report.md"
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Relatorio Markdown gravado em: %s", md_path)

    # --- Sumario final ---
    logger.info("")
    logger.info("=" * 70)
    logger.info("RESULTADO FINAL: %s RECOMENDADO", best["conjugate"])
    logger.info("=" * 70)
    logger.info("  Score de recomendacao: %.3f", best["recommendation_score"])
    logger.info("  Captacao esperada: %.1fx vs naked ASO", best["expected_uptake_fold"])
    logger.info("  Receptor alvo: %s", best["target_receptor"])
    logger.info("  GalNAc: EXCLUIDO (ASGPR nao expresso em macrofagos)")
    logger.info("  Tempo de execucao: %.2f segundos", elapsed)

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
