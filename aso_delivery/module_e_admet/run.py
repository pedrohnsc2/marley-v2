"""Modulo E — Perfil ADMET completo de MRL-ASO-001 em Canis lupus familiaris.

Pergunta central: qual e o perfil farmacocinetico e de seguranca de
MRL-ASO-001 apos injecao subcutanea em caes com leishmaniose visceral?

MRL-ASO-001 e um gapmer LNA-DNA-LNA de 25 nt com backbone PS completo
projetado para degradar o Spliced Leader RNA de Leishmania infantum
via recrutamento de RNase H1.

O perfil ADMET de PS-ASOs e fundamentalmente diferente de small molecules:
- Absorcao: SC com biodisponibilidade ~87% (vs oral irregular)
- Distribuicao: ligacao a albumina >90%, tropismo hepatoesplenico
- Metabolismo: degradacao por nucleases, NAO por CYP450
- Excrecao: renal de metabolitos curtos, meia-vida tecidual de semanas
- Toxicidade: efeitos de classe documentados (trombocitopenia, hepato, nefro)

NOTA ESPECIAL sobre TLR9: a ativacao imune por CpG e simultaneamente
um efeito de classe (toxicidade) e um mecanismo terapeutico (imunidade
anti-Leishmania). Este modulo avalia ambos os aspectos honestamente.

Este modulo executa seis analises independentes:
1. Absorcao (SC, biodisponibilidade, Tmax, Cmax)
2. Distribuicao (Vd, proteinas plasmaticas, tecidos-alvo)
3. Metabolismo (nucleases, mecanismo, DDI)
4. Excrecao (clearance, meia-vida terminal)
5. Toxicidade (efeitos de classe, orgao-especifica, TLR9)
6. Regime posologico (loading + manutencao, comparacao com miltefosine)

Saida: JSON + relatorio Markdown em aso_delivery/module_e_admet/results/

Referencias principais:
- Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
- Henry SP et al. (2008) Toxicology 252(1-3):97-105
- Frazier KS (2015) Toxicol Pathol 43(1):78-89
- Bennett CF (2019) Annu Rev Pharmacol Toxicol 59:447-464
- Yu RZ et al. (2007) J Pharmacol Exp Ther 320(1):108-116
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from aso_math.config import ASO_SEQUENCE as _ASO_SEQ
from aso_delivery.module_e_admet.pharmacokinetics import (
    BIOAVAILABILITY_SC_MEAN,
    BODY_WEIGHT_KG,
    CL_TOTAL_ML_MIN_KG,
    FORMULATION_CONC_MG_ML,
    FRACTION_UNBOUND,
    PLASMA_HALF_LIFE_HOURS,
    PROTEIN_BINDING_FRACTION,
    TISSUE_HALF_LIFE_DAYS_MEAN,
    TISSUE_PARTITION_COEFFICIENTS,
    VD_SS_L_PER_KG,
    compute_absorption,
    compute_distribution,
    compute_dosing_regimen,
    compute_excretion,
    compute_metabolism,
    compute_pk_simulation,
)
from aso_delivery.module_e_admet.toxicity import (
    compute_safety_profile,
)
from core.logger import get_logger

logger = get_logger("module_e_admet")

# ---------------------------------------------------------------------------
# Constantes do modulo
# ---------------------------------------------------------------------------

MODULE_NAME: Final[str] = "module_e_admet"
MODULE_VERSION: Final[str] = "1.0.0"

# Diretorio de resultados
RESULTS_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# Parametros do ASO (derivados de aso_math.config sem duplicar a sequencia)
ASO_LENGTH: Final[int] = len(_ASO_SEQ)
ASO_MW_DA: Final[float] = 8500.0
N_PS_LINKAGES: Final[int] = ASO_LENGTH - 1  # 24 para 25-mer

# Arquitetura do gapmer
LNA_5PRIME: Final[int] = 5
DNA_GAP: Final[int] = 15
LNA_3PRIME: Final[int] = 5


# ---------------------------------------------------------------------------
# 1. Analise de absorcao
# ---------------------------------------------------------------------------


def analyze_absorption() -> dict[str, Any]:
    """Calcula perfil de absorcao SC para doses de 5 e 10 mg/kg.

    PS-ASOs tem excelente biodisponibilidade SC (~87%) porque:
    - Moleculas de ~8500 Da sao grandes demais para absorcao oral
    - Mas pequenas o suficiente para drenagem linfatica/capilar SC
    - Backbone PS confere carga negativa -> solubilidade em agua

    A via SC e padrao para ASOs veterinarios (vs IV para emergencias).
    O deposito SC cria um reservatorio de liberacao sustentada,
    evitando picos plasmaticos que causam efeitos adversos.

    Returns:
        Perfil de absorcao para ambas as doses.
    """
    logger.info("Analise 1/6: Absorcao subcutanea")

    # Dose de manutencao: 5 mg/kg
    abs_5 = compute_absorption(dose_mg_kg=5.0)
    logger.info(
        "  Dose 5 mg/kg: abs = %.1f mg, inj vol = %.2f mL, "
        "F = %.0f%%, Cmax ~ %.1f ng/mL",
        abs_5.dose_absolute_mg,
        abs_5.injection_volume_ml,
        abs_5.bioavailability_percent,
        abs_5.cmax_estimated_ng_ml,
    )

    # Dose de carga: 10 mg/kg
    abs_10 = compute_absorption(dose_mg_kg=10.0)
    logger.info(
        "  Dose 10 mg/kg: abs = %.1f mg, inj vol = %.2f mL, "
        "F = %.0f%%, Cmax ~ %.1f ng/mL",
        abs_10.dose_absolute_mg,
        abs_10.injection_volume_ml,
        abs_10.bioavailability_percent,
        abs_10.cmax_estimated_ng_ml,
    )

    return {
        "route": "subcutaneous injection",
        "bioavailability_percent": abs_5.bioavailability_percent,
        "formulation_concentration_mg_ml": FORMULATION_CONC_MG_ML,
        "dose_5_mg_kg": asdict(abs_5),
        "dose_10_mg_kg": asdict(abs_10),
        "ka_per_hour": abs_5.ka_per_hour,
        "tmax_hours_range": list(abs_5.tmax_hours_range),
        "reference": "Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51",
    }


# ---------------------------------------------------------------------------
# 2. Analise de distribuicao
# ---------------------------------------------------------------------------


def analyze_distribution() -> dict[str, Any]:
    """Calcula perfil de distribuicao em tecidos caninos.

    O tropismo hepatoesplenico dos PS-ASOs e uma vantagem natural
    para leishmaniose visceral: L. infantum coloniza EXATAMENTE
    os tecidos onde PS-ASOs se acumulam em concentracoes elevadas.

    Figado: Kp = 30x -> celulas de Kupffer infectadas por amastigotas
    Baco: Kp = 15x -> macrofagos da zona marginal infectados
    Rim: Kp = 25x -> nao terapeuticamente relevante, mas monitorar

    Returns:
        Perfil de distribuicao completo.
    """
    logger.info("Analise 2/6: Distribuicao em tecidos caninos")

    # Dose de referencia: 5 mg/kg, F = 87%
    dose_absorbed = 5.0 * BODY_WEIGHT_KG * BIOAVAILABILITY_SC_MEAN
    dist = compute_distribution(dose_absorbed_mg=dose_absorbed)

    logger.info("  Vd,ss = %.2f L/kg (%.2f L total)", dist.vd_ss_l_per_kg, dist.vd_ss_l_total)
    logger.info("  Ligacao proteica: %.0f%% (fu = %.2f)", dist.protein_binding_percent, dist.fraction_unbound)
    logger.info("  Distribuicao tecidual (Kp x C_plasma):")

    for tissue, conc in dist.tissue_concentrations.items():
        kp = dist.tissue_partition_coefficients[tissue]
        relevance = ""
        if tissue in ("liver", "spleen", "lymph_nodes"):
            relevance = " <- TERAPEUTICAMENTE RELEVANTE (sitio de infeccao)"
        elif tissue == "kidney_cortex":
            relevance = " <- monitorar nefrotoxicidade"
        logger.info("    %s: Kp = %.1f, C = %.1f ng/mL%s", tissue, kp, conc, relevance)

    # Modelo de dois compartimentos
    two_comp = dist.two_compartment_params
    logger.info(
        "  Modelo 2-comp: alpha = %.4f/h, beta = %.4f/h",
        two_comp["alpha_per_hour"],
        two_comp["beta_per_hour"],
    )

    return {
        "vd_ss_l_per_kg": dist.vd_ss_l_per_kg,
        "vd_ss_l_total": dist.vd_ss_l_total,
        "protein_binding_percent": dist.protein_binding_percent,
        "fraction_unbound": dist.fraction_unbound,
        "tissue_concentrations_ng_ml": dist.tissue_concentrations,
        "tissue_partition_coefficients": dist.tissue_partition_coefficients,
        "two_compartment_model": dist.two_compartment_params,
        "therapeutic_relevance": {
            "liver": "PRIMARY site of L. infantum infection — Kp 30x ensures high ASO exposure",
            "spleen": "Major parasite reservoir — Kp 15x provides therapeutic concentrations",
            "kidney_cortex": "Not therapeutically relevant but monitor for nephrotoxicity",
            "lymph_nodes": "Secondary infection site — Kp 10x via lymphatic drainage",
        },
        "reference": "Yu RZ et al. (2007) J Pharmacol Exp Ther 320(1):108-116",
    }


# ---------------------------------------------------------------------------
# 3. Analise de metabolismo
# ---------------------------------------------------------------------------


def analyze_metabolism() -> dict[str, Any]:
    """Calcula perfil metabolico do ASO.

    Vantagem clinica dos PS-ASOs: sem metabolismo CYP450.
    Caes com leishmaniose frequentemente recebem multiplos farmacos:
    - Alopurinol (inibidor de xantina oxidase)
    - Antimoniais (meglumine antimoniate)
    - Miltefosine (alquilfosfolipidio)
    - Antibioticos (infeccoes secundarias)

    MRL-ASO-001 pode ser combinado com TODOS esses farmacos sem
    risco de interacao farmacocinetica, porque nao compete por CYP450.

    Returns:
        Perfil metabolico completo.
    """
    logger.info("Analise 3/6: Metabolismo")

    metab = compute_metabolism()

    logger.info("  Metabolismo CYP450: %s", "NAO" if not metab.cyp450_metabolism else "SIM")
    logger.info("  Degradacao: %s", metab.primary_degradation_pathway[:80] + "...")
    logger.info("  Meia-vida tecidual: %.0f dias", metab.tissue_half_life_days)
    logger.info("  Meia-vida plasmatica: %.1f horas", metab.plasma_half_life_hours)
    logger.info("  DDI: %s", metab.drug_drug_interactions[:80] + "...")

    return {
        "cyp450_metabolism": metab.cyp450_metabolism,
        "primary_degradation_pathway": metab.primary_degradation_pathway,
        "nuclease_protection": metab.nuclease_protection,
        "tissue_half_life_days": metab.tissue_half_life_days,
        "plasma_half_life_hours": metab.plasma_half_life_hours,
        "metabolites": metab.metabolites,
        "drug_drug_interactions": metab.drug_drug_interactions,
        "reference": "Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563",
    }


# ---------------------------------------------------------------------------
# 4. Analise de excrecao
# ---------------------------------------------------------------------------


def analyze_excretion() -> dict[str, Any]:
    """Calcula perfil de excrecao do ASO.

    A excrecao de PS-ASOs e dominada pela eliminacao renal de
    metabolitos encurtados. O ASO intacto e retido por ligacao
    a proteinas plasmaticas e acumulacao tecidual.

    A meia-vida terminal longa (3 semanas) e uma vantagem para
    leishmaniose: permite dosagem semanal (vs diaria para miltefosine).

    Returns:
        Perfil de excrecao completo.
    """
    logger.info("Analise 4/6: Excrecao")

    exc = compute_excretion()

    logger.info("  CL total: %.1f mL/min (%.1f mL/min/kg)", exc.cl_total_ml_min, exc.cl_total_ml_min_kg)
    logger.info("  CL renal: %.1f mL/min (%.0f%% do total)", exc.cl_renal_ml_min, exc.renal_fraction * 100)
    logger.info("  Meia-vida terminal: %.0f dias (%.0f horas)", exc.terminal_half_life_days, exc.terminal_half_life_hours)
    logger.info("  Steady state em: %.0f dias (~%.0f semanas)", exc.time_to_steady_state_days, exc.time_to_steady_state_days / 7)

    return {
        "cl_total_ml_min_kg": exc.cl_total_ml_min_kg,
        "cl_total_ml_min": exc.cl_total_ml_min,
        "cl_renal_ml_min": exc.cl_renal_ml_min,
        "cl_hepatic_ml_min": exc.cl_hepatic_ml_min,
        "renal_fraction_percent": exc.renal_fraction * 100,
        "terminal_half_life_days": exc.terminal_half_life_days,
        "terminal_half_life_hours": exc.terminal_half_life_hours,
        "time_to_steady_state_days": exc.time_to_steady_state_days,
        "reference": "Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51",
    }


# ---------------------------------------------------------------------------
# 5. Analise de toxicidade
# ---------------------------------------------------------------------------


def analyze_toxicity() -> dict[str, Any]:
    """Avalia perfil de seguranca completo do ASO.

    AVALIACAO HONESTA: este modulo reporta todos os riscos conhecidos
    de PS-ASOs sem minimiza-los. Os efeitos adversos sao reais e
    documentados em humanos (mipomersen, inotersen — ambos FDA-approved).

    Contexto: leishmaniose visceral canina e FATAL em >90% dos casos
    sem tratamento. O risco-beneficio de um tratamento com efeitos
    adversos manejáveis e claramente favoravel.

    Returns:
        Perfil de seguranca integrado.
    """
    logger.info("Analise 5/6: Toxicidade e seguranca")

    safety = compute_safety_profile(_ASO_SEQ)

    logger.info("  Indice terapeutico: %.1f (NOAEL %.0f / dose %.0f mg/kg/semana)",
                safety.therapeutic_index, safety.noael_mg_kg_week, safety.dose_mg_kg_week)
    logger.info("  Margem de seguranca: %.1fx (MTD / dose)", safety.safety_margin)
    logger.info("  Classificacao: %s", safety.overall_risk_classification)

    # Reportar cada endpoint
    for ep in safety.endpoints:
        logger.info(
            "  %s: severidade=%s, risk_score=%.1f/10, reversivel=%s",
            ep.name, ep.severity, ep.risk_score,
            "sim" if ep.reversible else "NAO",
        )

    # TLR9 — caso especial
    tlr9 = safety.tlr9_assessment
    logger.info("  TLR9: %d motivo(s) CpG, ativacao %s",
                tlr9.cpg_motifs_in_sequence,
                "ESPERADA (dual function)" if tlr9.tlr9_activation_expected else "minima")

    # Montar resultado
    endpoints_list = [asdict(ep) for ep in safety.endpoints]
    organ_list = [asdict(op) for op in safety.organ_profiles]

    return {
        "therapeutic_index": safety.therapeutic_index,
        "noael_mg_kg_week": safety.noael_mg_kg_week,
        "mtd_mg_kg_week": safety.mtd_mg_kg_week,
        "dose_mg_kg_week": safety.dose_mg_kg_week,
        "safety_margin": safety.safety_margin,
        "overall_risk_classification": safety.overall_risk_classification,
        "endpoints": endpoints_list,
        "organ_toxicity_profiles": organ_list,
        "tlr9_assessment": asdict(tlr9),
        "monitoring_schedule": safety.monitoring_schedule,
        "contraindications": safety.contraindications,
        "drug_interactions": safety.drug_interactions,
        "reference": "Henry SP et al. (2008) Toxicology 252(1-3):97-105",
    }


# ---------------------------------------------------------------------------
# 6. Regime posologico
# ---------------------------------------------------------------------------


def analyze_dosing() -> dict[str, Any]:
    """Calcula regime posologico recomendado.

    Baseado em escalonamento alometrico de mipomersen/inotersen
    e adaptado para leishmaniose visceral canina:
    - Loading dose para saturar tecidos rapidamente (parasitemia alta)
    - Manutencao semanal para manter concentracoes terapeuticas
    - Duracao de 12 semanas (baseada em clearance parasitario)

    Inclui comparacao com miltefosine (tratamento padrao).

    Returns:
        Regime posologico completo.
    """
    logger.info("Analise 6/6: Regime posologico recomendado")

    regimen = compute_dosing_regimen()

    logger.info("  Loading: %.0f mg/kg %s por %d semanas",
                regimen.loading_dose_mg_kg, regimen.loading_frequency,
                regimen.loading_duration_weeks)
    logger.info("  Manutencao: %.0f mg/kg %s por %d semanas",
                regimen.maintenance_dose_mg_kg, regimen.maintenance_frequency,
                regimen.treatment_duration_weeks)
    logger.info("  Duracao total: %d semanas", regimen.total_treatment_weeks)
    logger.info("  Vol injecao (loading): %.2f mL", regimen.injection_volume_loading_ml)
    logger.info("  Vol injecao (manut): %.2f mL", regimen.injection_volume_maintenance_ml)
    logger.info("  Ctrough estimado: %.1f ng/mL", regimen.estimated_trough_ng_ml)
    logger.info("  Indice terapeutico: %.1f", regimen.therapeutic_index)

    return asdict(regimen)


# ---------------------------------------------------------------------------
# Simulacao PK temporal
# ---------------------------------------------------------------------------


def run_pk_simulation() -> list[dict[str, float]]:
    """Executa simulacao PK de 168 horas (1 semana) apos dose unica.

    Modelo de dois compartimentos com absorcao SC de primeira ordem.
    Gera curva concentracao-tempo para plasma e tecidos.

    Returns:
        Lista de pontos temporais com concentracoes.
    """
    logger.info("Simulacao PK: modelo 2-compartimentos, 168h pos-dose")

    dose_mg = 5.0 * BODY_WEIGHT_KG  # 75 mg para cao de 15 kg
    sim = compute_pk_simulation(dose_mg=dose_mg, duration_hours=168.0)

    # Reportar pontos chave
    if sim:
        # Encontrar Cmax e Tmax
        cmax_point = max(sim, key=lambda p: p.plasma_conc_ng_ml)
        logger.info(
            "  Cmax plasma: %.1f ng/mL a t = %.1f h",
            cmax_point.plasma_conc_ng_ml, cmax_point.time_hours,
        )
        # Concentracao em 24h e 168h
        for target_t in [24.0, 168.0]:
            pts = [p for p in sim if abs(p.time_hours - target_t) < 1.0]
            if pts:
                p = pts[0]
                logger.info(
                    "  t = %.0f h: plasma = %.1f ng/mL, tecido = %.1f ng/mL",
                    p.time_hours, p.plasma_conc_ng_ml, p.tissue_conc_ng_ml,
                )

    return [asdict(p) for p in sim]


# ---------------------------------------------------------------------------
# Geracao de relatorio Markdown
# ---------------------------------------------------------------------------


def generate_markdown_report(results: dict[str, Any]) -> str:
    """Gera relatorio Markdown completo do perfil ADMET.

    Estrutura:
    1. Executive Summary
    2. Absorption (SC)
    3. Distribution (tissues)
    4. Metabolism (nucleases)
    5. Excretion (renal)
    6. Toxicity (class effects + TLR9)
    7. Dosing Regimen (loading + maintenance)
    8. PK Simulation
    9. References

    Args:
        results: Dicionario com todos os resultados das 6 analises.

    Returns:
        String com conteudo Markdown do relatorio.
    """
    lines: list[str] = []

    lines.append("# Module E: Full ADMET Profile of MRL-ASO-001 in Canis lupus familiaris")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"**Module version:** {MODULE_VERSION}")
    lines.append(f"**Target animal:** Canis lupus familiaris, {BODY_WEIGHT_KG} kg reference weight")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(results["executive_summary"])
    lines.append("")

    # ---- Secao 1: Absorption ----
    lines.append("## 1. Absorption")
    lines.append("")
    lines.append("**Route:** Subcutaneous injection (standard for veterinary ASOs)")
    lines.append("")

    abs_data = results["absorption"]
    lines.append(f"**Bioavailability:** {abs_data['bioavailability_percent']}%")
    lines.append(f"**Formulation:** {abs_data['formulation_concentration_mg_ml']} mg/mL aqueous solution")
    lines.append(f"**Absorption rate constant (ka):** {abs_data['ka_per_hour']} h^-1")
    tmax = abs_data["tmax_hours_range"]
    lines.append(f"**Tmax:** {tmax[0]}-{tmax[1]} hours post-injection")
    lines.append("")

    lines.append("| Parameter | 5 mg/kg (maintenance) | 10 mg/kg (loading) |")
    lines.append("|---|---|---|")
    d5 = abs_data["dose_5_mg_kg"]
    d10 = abs_data["dose_10_mg_kg"]
    lines.append(f"| Absolute dose (mg) | {d5['dose_absolute_mg']} | {d10['dose_absolute_mg']} |")
    lines.append(f"| Injection volume (mL) | {d5['injection_volume_ml']} | {d10['injection_volume_ml']} |")
    lines.append(f"| Dose absorbed (mg) | {d5['dose_absorbed_mg']} | {d10['dose_absorbed_mg']} |")
    lines.append(f"| Estimated Cmax (ng/mL) | {d5['cmax_estimated_ng_ml']} | {d10['cmax_estimated_ng_ml']} |")
    lines.append(f"| Estimated AUC (ng*h/mL) | {d5['auc_estimated_ng_h_ml']} | {d10['auc_estimated_ng_h_ml']} |")
    lines.append("")

    # ---- Secao 2: Distribution ----
    dist = results["distribution"]
    lines.append("## 2. Distribution")
    lines.append("")
    lines.append(f"**Vd,ss:** {dist['vd_ss_l_per_kg']} L/kg ({dist['vd_ss_l_total']} L total)")
    lines.append(f"**Protein binding:** {dist['protein_binding_percent']}% (primarily albumin)")
    lines.append(f"**Fraction unbound (fu):** {dist['fraction_unbound']}")
    lines.append("")

    lines.append("### Tissue Distribution (Kp and estimated concentrations after 5 mg/kg SC)")
    lines.append("")
    lines.append("| Tissue | Kp (tissue:plasma) | Estimated C (ng/mL) | Therapeutic Relevance |")
    lines.append("|---|---|---|---|")

    relevance_map = dist.get("therapeutic_relevance", {})
    for tissue, kp in dist["tissue_partition_coefficients"].items():
        conc = dist["tissue_concentrations_ng_ml"].get(tissue, "N/A")
        rel = relevance_map.get(tissue, "")
        lines.append(f"| {tissue} | {kp} | {conc} | {rel[:60]}{'...' if len(rel) > 60 else ''} |")
    lines.append("")

    lines.append("### Two-Compartment Model Parameters")
    lines.append("")
    two_comp = dist["two_compartment_model"]
    lines.append(f"- Vc (central volume): {two_comp['vc_l_per_kg']} L/kg")
    lines.append(f"- Vt (peripheral volume): {two_comp['vt_l_per_kg']} L/kg")
    lines.append(f"- k12 (central -> peripheral): {two_comp['k12_per_hour']} h^-1")
    lines.append(f"- k21 (peripheral -> central): {two_comp['k21_per_hour']} h^-1")
    lines.append(f"- ke (elimination from central): {two_comp['ke_per_hour']} h^-1")
    lines.append(f"- alpha (distribution phase): {two_comp['alpha_per_hour']} h^-1")
    lines.append(f"- beta (elimination phase): {two_comp['beta_per_hour']} h^-1")
    lines.append("")

    lines.append("**Key insight:** PS-ASO hepatosplenic tropism is a natural advantage for visceral")
    lines.append("leishmaniasis. The same tissue distribution pattern that limits PS-ASO utility in")
    lines.append("CNS diseases makes MRL-ASO-001 inherently suited for targeting L. infantum")
    lines.append("amastigotes in liver and spleen macrophages.")
    lines.append("")

    # ---- Secao 3: Metabolism ----
    met = results["metabolism"]
    lines.append("## 3. Metabolism")
    lines.append("")
    lines.append(f"**CYP450 metabolism:** {'Yes' if met['cyp450_metabolism'] else 'No'}")
    lines.append(f"**Primary pathway:** {met['primary_degradation_pathway']}")
    lines.append(f"**Nuclease protection:** {met['nuclease_protection']}")
    lines.append(f"**Tissue half-life:** {met['tissue_half_life_days']} days")
    lines.append(f"**Plasma half-life:** {met['plasma_half_life_hours']} hours")
    lines.append("")

    lines.append("### Metabolites")
    lines.append("")
    for m in met["metabolites"]:
        lines.append(f"- {m}")
    lines.append("")

    lines.append("### Drug-Drug Interactions")
    lines.append("")
    lines.append(met["drug_drug_interactions"])
    lines.append("")

    # ---- Secao 4: Excretion ----
    exc = results["excretion"]
    lines.append("## 4. Excretion")
    lines.append("")
    lines.append(f"**Total clearance:** {exc['cl_total_ml_min']} mL/min ({exc['cl_total_ml_min_kg']} mL/min/kg)")
    lines.append(f"**Renal clearance:** {exc['cl_renal_ml_min']} mL/min ({exc['renal_fraction_percent']}% of total)")
    lines.append(f"**Hepatic clearance:** {exc['cl_hepatic_ml_min']} mL/min")
    lines.append(f"**Terminal half-life:** {exc['terminal_half_life_days']} days ({exc['terminal_half_life_hours']} hours)")
    lines.append(f"**Time to steady state:** {exc['time_to_steady_state_days']} days (~{exc['time_to_steady_state_days']/7:.0f} weeks)")
    lines.append("")

    # ---- Secao 5: Toxicity ----
    tox = results["toxicity"]
    lines.append("## 5. Toxicity Assessment")
    lines.append("")
    lines.append(f"**Therapeutic index:** {tox['therapeutic_index']} (NOAEL {tox['noael_mg_kg_week']} / dose {tox['dose_mg_kg_week']} mg/kg/week)")
    lines.append(f"**Safety margin:** {tox['safety_margin']}x (MTD / dose)")
    lines.append(f"**Overall risk:** {tox['overall_risk_classification']}")
    lines.append("")

    lines.append("### Class Effects")
    lines.append("")
    lines.append("| Adverse Effect | Severity | Risk Score | Dose-Dependent | Reversible |")
    lines.append("|---|---|---|---|---|")
    for ep in tox["endpoints"]:
        lines.append(
            f"| {ep['name']} | {ep['severity']} | {ep['risk_score']}/10 | "
            f"{'Yes' if ep['dose_dependent'] else 'No'} | "
            f"{'Yes' if ep['reversible'] else 'No'} |"
        )
    lines.append("")

    # TLR9 — secao especial
    tlr9 = tox["tlr9_assessment"]
    lines.append("### TLR9 Activation: Dual Function Assessment")
    lines.append("")
    lines.append(f"**CpG motifs in MRL-ASO-001:** {tlr9['cpg_motifs_in_sequence']}")
    lines.append(f"**TLR9 activation expected:** {'Yes' if tlr9['tlr9_activation_expected'] else 'No'}")
    lines.append("")
    lines.append("**Therapeutic benefit:**")
    lines.append(tlr9["therapeutic_benefit"])
    lines.append("")
    lines.append("**Toxicity risk:**")
    lines.append(tlr9["toxicity_risk"])
    lines.append("")
    lines.append("**Dose window:**")
    lines.append(tlr9["dose_window_description"])
    lines.append("")
    lines.append("**Dual function assessment:**")
    lines.append(tlr9["dual_function_assessment"])
    lines.append("")

    # Orgao-especifica
    lines.append("### Organ-Specific Toxicity")
    lines.append("")
    lines.append("| Organ | Kp | Primary Cells | Threshold (mg/kg) | Expected at 5 mg/kg |")
    lines.append("|---|---|---|---|---|")
    for op in tox["organ_toxicity_profiles"]:
        lines.append(
            f"| {op['organ']} | {op['accumulation_factor']} | "
            f"{op['primary_cell_type_affected'][:30]}... | "
            f"{op['threshold_dose_mg_kg']} | "
            f"{op['expected_at_5_mg_kg'][:40]}... |"
        )
    lines.append("")

    # Monitoramento
    lines.append("### Monitoring Schedule")
    lines.append("")
    for phase, tests in tox["monitoring_schedule"].items():
        lines.append(f"**{phase.replace('_', ' ').title()}:** {tests}")
        lines.append("")

    # Contraindicacoes
    lines.append("### Contraindications")
    lines.append("")
    for ci in tox["contraindications"]:
        lines.append(f"- {ci}")
    lines.append("")

    # ---- Secao 6: Dosing Regimen ----
    dos = results["dosing_regimen"]
    lines.append("## 6. Recommended Dosing Regimen")
    lines.append("")
    lines.append("### Loading Phase")
    lines.append("")
    lines.append(f"- **Dose:** {dos['loading_dose_mg_kg']} mg/kg ({dos['loading_dose_mg']} mg absolute)")
    lines.append(f"- **Frequency:** {dos['loading_frequency']}")
    lines.append(f"- **Duration:** {dos['loading_duration_weeks']} weeks")
    lines.append(f"- **Injection volume:** {dos['injection_volume_loading_ml']} mL")
    lines.append("")
    lines.append("### Maintenance Phase")
    lines.append("")
    lines.append(f"- **Dose:** {dos['maintenance_dose_mg_kg']} mg/kg ({dos['maintenance_dose_mg']} mg absolute)")
    lines.append(f"- **Frequency:** {dos['maintenance_frequency']}")
    lines.append(f"- **Duration:** {dos['treatment_duration_weeks']} weeks")
    lines.append(f"- **Injection volume:** {dos['injection_volume_maintenance_ml']} mL")
    lines.append("")
    lines.append(f"**Total treatment duration:** {dos['total_treatment_weeks']} weeks")
    lines.append(f"**Estimated trough concentration:** {dos['estimated_trough_ng_ml']} ng/mL")
    lines.append(f"**Therapeutic index:** {dos['therapeutic_index']}")
    lines.append("")

    # Comparacao com miltefosine
    comp = dos["comparison_miltefosine"]
    lines.append("### Comparison with Miltefosine (Current Standard of Care)")
    lines.append("")
    lines.append("| Parameter | Miltefosine | MRL-ASO-001 |")
    lines.append("|---|---|---|")
    lines.append(f"| Dose | {comp['miltefosine_dose']} | {comp['mrl_aso_001_dose']} |")
    lines.append(f"| Efficacy | {comp['miltefosine_efficacy']} | {comp['mrl_aso_001_efficacy']} |")
    lines.append(f"| Toxicity | {comp['miltefosine_toxicity']} | {comp['mrl_aso_001_toxicity']} |")
    lines.append(f"| Resistance | {comp['miltefosine_resistance']} | {comp['mrl_aso_001_resistance']} |")
    lines.append(f"| Cost | {comp['miltefosine_cost']} | {comp['mrl_aso_001_cost']} |")
    lines.append("")
    lines.append(f"**Key advantage of MRL-ASO-001:** {comp['key_advantage']}")
    lines.append("")

    # ---- PK Simulation ----
    if "pk_simulation" in results and results["pk_simulation"]:
        sim = results["pk_simulation"]
        lines.append("## 7. PK Simulation (5 mg/kg SC, Single Dose, 168h)")
        lines.append("")
        lines.append("| Time (h) | Plasma (ng/mL) | Tissue (ng/mL) | Absorbed (mg) | Eliminated (mg) |")
        lines.append("|---|---|---|---|---|")

        # Selecionar pontos chave para a tabela (nao todos)
        key_times = {0, 1, 2, 3, 4, 6, 8, 12, 24, 48, 72, 96, 120, 144, 168}
        for p in sim:
            t = p["time_hours"]
            if t in key_times or (t <= 12 and t % 2 == 0):
                lines.append(
                    f"| {t:.0f} | {p['plasma_conc_ng_ml']:.1f} | "
                    f"{p['tissue_conc_ng_ml']:.1f} | {p['amount_absorbed_mg']:.2f} | "
                    f"{p['amount_eliminated_mg']:.2f} |"
                )
        lines.append("")

    # ---- Overall Conclusion ----
    lines.append("## Overall Conclusion")
    lines.append("")
    lines.append(results["overall_conclusion"])
    lines.append("")

    # ---- References ----
    lines.append("## References")
    lines.append("")
    lines.append("1. Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51 — ASO pharmacokinetics in animals")
    lines.append("2. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review")
    lines.append("3. Yu RZ et al. (2007) J Pharmacol Exp Ther 320(1):108-116 — PS-ASO PK in monkeys/dogs")
    lines.append("4. Henry SP et al. (2008) Toxicology 252(1-3):97-105 — ASO toxicology in animals")
    lines.append("5. Frazier KS (2015) Toxicol Pathol 43(1):78-89 — ASO-related organ toxicity")
    lines.append("6. Bennett CF (2019) Annu Rev Pharmacol Toxicol 59:447-464 — ASO pharmacology update")
    lines.append("7. Volpi S et al. (2012) J Immunol 188(12):5890-5897 — TLR9 activation by CpG in canines")
    lines.append("8. Mipomersen FDA label (2013) — SC bioavailability, hepatotoxicity, ISR data")
    lines.append("9. Inotersen FDA label (2018) — thrombocytopenia black box warning, PK data")
    lines.append("10. Raal FJ et al. (2010) Lancet 375(9719):998-1006 — mipomersen clinical efficacy")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Orquestrador principal
# ---------------------------------------------------------------------------


def main() -> dict[str, Any]:
    """Executa perfil ADMET completo de MRL-ASO-001.

    Fluxo:
        1. Absorcao SC (biodisponibilidade, Cmax, Tmax)
        2. Distribuicao (Vd, tecidos-alvo, modelo 2-compartimentos)
        3. Metabolismo (nucleases, DDI, meia-vida)
        4. Excrecao (clearance renal, meia-vida terminal)
        5. Toxicidade (efeitos de classe, TLR9, orgao-especifica)
        6. Regime posologico (loading + manutencao)
        7. Simulacao PK temporal
        8. Gerar relatorio Markdown
        9. Gravar resultados JSON + Markdown

    Returns:
        Envelope completo com todos os resultados.
    """
    logger.info("=" * 70)
    logger.info("MODULO E: Perfil ADMET Completo — MRL-ASO-001")
    logger.info("=" * 70)
    logger.info("Animal: Canis lupus familiaris, %.0f kg", BODY_WEIGHT_KG)
    logger.info("ASO: 25-nt LNA-DNA-LNA gapmer, PS backbone, MW ~8500 Da")
    logger.info("")

    start_time = time.time()

    # --- Analise 1: Absorcao ---
    absorption = analyze_absorption()
    logger.info("")

    # --- Analise 2: Distribuicao ---
    distribution = analyze_distribution()
    logger.info("")

    # --- Analise 3: Metabolismo ---
    metabolism = analyze_metabolism()
    logger.info("")

    # --- Analise 4: Excrecao ---
    excretion = analyze_excretion()
    logger.info("")

    # --- Analise 5: Toxicidade ---
    toxicity = analyze_toxicity()
    logger.info("")

    # --- Analise 6: Regime posologico ---
    dosing = analyze_dosing()
    logger.info("")

    # --- Simulacao PK ---
    pk_sim = run_pk_simulation()
    logger.info("")

    # --- Compilar resultados ---
    elapsed = round(time.time() - start_time, 2)

    # Sumario executivo
    ti = toxicity["therapeutic_index"]
    t_half = excretion["terminal_half_life_days"]
    bioav = absorption["bioavailability_percent"]

    executive_summary = (
        f"MRL-ASO-001 has a favorable ADMET profile for treating canine visceral "
        f"leishmaniasis via subcutaneous injection. "
        f"**Absorption:** {bioav}% bioavailability after SC injection, Tmax 2-4 hours. "
        f"**Distribution:** Vd,ss = {VD_SS_L_PER_KG} L/kg with preferential accumulation "
        f"in liver (Kp 30x) and spleen (Kp 15x) — the primary sites of L. infantum "
        f"infection. **Metabolism:** No CYP450 involvement; degraded by nucleases with "
        f"LNA flanks providing exonuclease resistance. **Excretion:** Renal elimination "
        f"of shortened metabolites; terminal half-life {t_half:.0f} days enabling weekly "
        f"dosing. **Toxicity:** Therapeutic index {ti:.1f}x (NOAEL/dose) with manageable "
        f"class effects (ISR, mild thrombocytopenia risk). TLR9 activation by CpG motifs "
        f"provides DUAL therapeutic function: direct SL RNA knockdown AND innate immune "
        f"activation against Leishmania."
    )

    # Conclusao geral
    overall_conclusion = (
        f"MRL-ASO-001 demonstrates pharmacokinetic properties well-suited for canine "
        f"visceral leishmaniasis treatment:\n\n"
        f"1. **Natural tissue tropism:** PS-ASO hepatosplenic accumulation (Kp 15-30x) "
        f"delivers the drug directly to L. infantum-infected macrophages without "
        f"need for specialized delivery systems.\n\n"
        f"2. **Long half-life:** Terminal t1/2 of {t_half:.0f} days enables convenient "
        f"once-weekly SC dosing during maintenance (vs daily oral miltefosine).\n\n"
        f"3. **No drug interactions:** Absence of CYP450 metabolism allows safe "
        f"combination with standard leishmaniasis drugs (allopurinol, antimonials).\n\n"
        f"4. **Favorable therapeutic index:** TI of {ti:.1f}x provides adequate safety "
        f"margin for chronic treatment. The most significant risk (thrombocytopenia) "
        f"requires monitoring but is manageable with dose adjustment.\n\n"
        f"5. **Dual mechanism:** TLR9 activation by CpG motifs in the PS backbone "
        f"provides built-in immune adjuvant activity — a unique therapeutic advantage "
        f"not available with any current anti-leishmanial treatment.\n\n"
        f"**HONEST RISKS:** Thrombocytopenia (5-10% clinically significant), injection "
        f"site reactions (70-90%), mild transaminase elevation (5-10%). These are "
        f"REAL risks that require monitoring, but are ACCEPTABLE in the context of "
        f"a disease with >90% mortality without treatment.\n\n"
        f"**Recommended regimen:** Loading phase (10 mg/kg 2x/week for 2 weeks) "
        f"followed by maintenance (5 mg/kg 1x/week for 10 weeks). "
        f"Total treatment: 12 weeks with weekly CBC monitoring."
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
            "molecular_weight_da": ASO_MW_DA,
            "gapmer_design": f"{LNA_5PRIME}-{DNA_GAP}-{LNA_3PRIME} (LNA-DNA-LNA)",
            "backbone": "full phosphorothioate (PS)",
            "ps_linkages": N_PS_LINKAGES,
        },
        "target_animal": {
            "species": "Canis lupus familiaris",
            "body_weight_kg": BODY_WEIGHT_KG,
            "disease": "visceral leishmaniasis (L. infantum)",
        },
        "absorption": absorption,
        "distribution": distribution,
        "metabolism": metabolism,
        "excretion": excretion,
        "toxicity": toxicity,
        "dosing_regimen": dosing,
        "pk_simulation": pk_sim,
    }

    # --- Gravar resultados ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = RESULTS_DIR / "module_e_admet.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    logger.info("Resultados JSON gravados em: %s", json_path)

    # Markdown
    report_md = generate_markdown_report(results)
    md_path = RESULTS_DIR / "module_e_admet_report.md"
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Relatorio Markdown gravado em: %s", md_path)

    # --- Sumario final ---
    logger.info("")
    logger.info("=" * 70)
    logger.info("MODULO E: PERFIL ADMET COMPLETO")
    logger.info("=" * 70)
    logger.info("  Biodisponibilidade SC: %.0f%%", bioav)
    logger.info("  Vd,ss: %.2f L/kg", VD_SS_L_PER_KG)
    logger.info("  Meia-vida terminal: %.0f dias", t_half)
    logger.info("  Indice terapeutico: %.1f", ti)
    logger.info("  Classificacao de risco: %s", toxicity["overall_risk_classification"])
    logger.info("  Regime: loading 10 mg/kg 2x/sem (2 sem) + manut 5 mg/kg 1x/sem (10 sem)")
    logger.info("  Tempo de execucao: %.2f segundos", elapsed)

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
