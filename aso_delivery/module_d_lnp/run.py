"""Modulo D — Nanoparticulas lipidicas (LNP) para entrega de MRL-ASO-001.

Pergunta central: encapsulacao em LNP pode melhorar a entrega de
MRL-ASO-001 a macrofagos infectados por L. infantum?

O macrofago infectado e o alvo celular primario. Amastigotas de
L. infantum residem e se replicam DENTRO do fagolisossomo (pH 4.5),
onde o SL RNA alvo esta localizado. A LNP ideal deve:

1. Proteger o ASO durante circulacao (pH 7.4, nucleases sericas)
2. Ser captada preferencialmente por macrofagos (targeting por manose)
3. Liberar a carga no fagolisossomo (liberacao pH-responsiva)
4. Ser estavel em armazenamento (cadeia fria ou liofilizada)

Este modulo executa seis analises:
1. Composicao lipidica da LNP (4 componentes, razoes molares)
2. Otimizacao de razao N/P (eficiencia de encapsulacao)
3. Cinetica de liberacao pH-responsiva (4 compartimentos)
4. Modificacoes para targeting de macrofagos (manose-PEG)
5. Estabilidade e armazenamento (3 temperaturas)
6. Comparacao: LNP vs ASO livre

Saida: JSON + relatorio Markdown em aso_delivery/module_d_lnp/results/

Referencias principais:
- Cullis PR, Hope MJ (2017) Mol Ther 25(7):1467-1475 — LNP design
- Jayaraman M et al. (2012) Angew Chem 51(34):8529-8533 — DLin-MC3-DMA
- Semple SC et al. (2010) Nat Biotechnol 28(2):172-176 — LNP optimization
- Patel S et al. (2020) J Control Release 327:146-160 — macrophage targeting
- Kulkarni JA et al. (2018) Nanoscale 10(10):4567-4573 — release mechanism
- Sahay G et al. (2013) Nat Biotechnol 31(7):653-658 — endosomal escape
- Schoenmaker L et al. (2021) Int J Pharm 601:120586 — LNP stability
- Ball RL et al. (2017) Drug Deliv Transl Res 7(1):89-103 — lyophilization
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from aso_delivery.module_d_lnp.formulation import (
    MANNOSE_UPTAKE_FACTOR_MEAN,
    MW_DLIN_MC3_DMA,
    MW_DSPC,
    MW_CHOLESTEROL,
    MW_PEG_LIPID,
    MW_MANNOSE_PEG,
    compute_lnp_composition,
    compute_macrophage_targeting,
    compute_storage_profiles,
    find_optimal_np,
    scan_np_ratios,
)
from aso_delivery.module_d_lnp.release import (
    compute_release_profile,
)
from aso_math.config import ASO_LENGTH, ASO_SEQUENCE
from core.logger import get_logger

logger = get_logger("module_d_lnp")

# ---------------------------------------------------------------------------
# Constantes do modulo
# ---------------------------------------------------------------------------

MODULE_NAME: Final[str] = "module_d_lnp"
MODULE_VERSION: Final[str] = "1.0.0"

# Diretorio de resultados
RESULTS_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# Cargas negativas no backbone PS do ASO (1 por ligacao internucleotidica)
ASO_PS_CHARGES: Final[int] = ASO_LENGTH - 1  # 24 para 25 nt

# Custos estimados por dose (USD) — ordem de grandeza para comparacao
# Ref: estimativas da literatura, nao valores exatos de fabricacao
COST_NAKED_ASO_USD: Final[float] = 5.0       # ASO livre sintetico (custo de GMP)
COST_LNP_ASO_USD: Final[float] = 15.0        # ASO em LNP (formulacao + lipidios)

# Complexidade relativa de fabricacao (escala 1-5)
COMPLEXITY_NAKED: Final[int] = 2    # sintese + purificacao
COMPLEXITY_LNP: Final[int] = 4      # sintese + formulacao + caracterizacao


# ---------------------------------------------------------------------------
# 1. Composicao lipidica
# ---------------------------------------------------------------------------


def analyze_lnp_composition() -> dict[str, Any]:
    """Calcula e reporta composicao da formulacao LNP.

    A formulacao padrao (Onpattro-like) usa quatro componentes com
    funcoes distintas:
    - Lipidio ionizavel (50%): complexacao com ASO + liberacao pH-responsiva
    - DSPC (10%): estrutura da bicamada
    - Colesterol (38.5%): estabilidade mecanica
    - PEG-lipidio (1.5%): stealth (evita clearance)

    Returns:
        Dicionario com composicao e propriedades.
    """
    logger.info("Analise 1/6: Composicao lipidica da LNP")

    comp = compute_lnp_composition()

    logger.info("  Formulacao: %s / %s / %s / %s",
                comp.ionizable_lipid, comp.helper_lipid,
                comp.sterol, comp.peg_lipid)
    logger.info("  Razoes molares: %.1f%% / %.1f%% / %.1f%% / %.1f%%",
                comp.mol_fraction_ionizable * 100,
                comp.mol_fraction_helper * 100,
                comp.mol_fraction_sterol * 100,
                comp.mol_fraction_peg * 100)
    logger.info("  MW medio ponderado: %.2f g/mol", comp.weighted_average_mw)
    logger.info("  Carga a pH 7.4: %.6f (stealth)", comp.charge_ph_7_4)
    logger.info("  Carga a pH 4.5: %.6f (cationico)", comp.charge_ph_4_5)
    logger.info("  pKa do ionizavel: %.2f", comp.pka_ionizable)

    return {
        "composition": asdict(comp),
        "components": {
            "ionizable_lipid": {
                "name": "DLin-MC3-DMA",
                "mw_g_per_mol": MW_DLIN_MC3_DMA,
                "mol_fraction": comp.mol_fraction_ionizable,
                "role": "pH-responsive encapsulation and release",
                "pka": comp.pka_ionizable,
            },
            "helper_lipid": {
                "name": "DSPC",
                "mw_g_per_mol": MW_DSPC,
                "mol_fraction": comp.mol_fraction_helper,
                "role": "Bilayer structural integrity",
            },
            "sterol": {
                "name": "Cholesterol",
                "mw_g_per_mol": MW_CHOLESTEROL,
                "mol_fraction": comp.mol_fraction_sterol,
                "role": "Membrane stability and fusogenicity",
            },
            "peg_lipid": {
                "name": "DMG-PEG2000",
                "mw_g_per_mol": MW_PEG_LIPID,
                "mol_fraction": comp.mol_fraction_peg,
                "role": "Stealth layer (prevents opsonization and RES clearance)",
            },
        },
    }


# ---------------------------------------------------------------------------
# 2. Otimizacao de razao N/P
# ---------------------------------------------------------------------------


def analyze_np_optimization() -> dict[str, Any]:
    """Varre razoes N/P de 1 a 20 e identifica o otimo para PS-ASOs.

    N/P = moles de amina ionizavel / moles de fosfato no ASO.
    Para MRL-ASO-001 (25 nt, backbone PS): 24 cargas negativas.

    A faixa otima para PS-ASOs (4-8) difere de mRNA (6-12) porque:
    - PS tem carga negativa PERMANENTE (vs PO que e -1 por fosfato)
    - PS tem maior afinidade eletrostatica por lipidios cationicos
    - Menor N/P necessario -> menos lipidio -> menor toxicidade

    Returns:
        Dicionario com varredura completa e N/P otimo.
    """
    logger.info("Analise 2/6: Otimizacao de razao N/P")
    logger.info("  ASO: %d nt, %d cargas PS negativas", ASO_LENGTH, ASO_PS_CHARGES)

    # Varredura com passo de 1
    scan = scan_np_ratios(ASO_LENGTH, np_min=1.0, np_max=20.0, step=1.0)
    optimal = find_optimal_np(scan)

    logger.info("  Varredura: N/P 1 a 20 (20 pontos)")
    logger.info("  N/P otimo: %.0f", optimal.np_ratio)
    logger.info("  Encapsulacao: %.1f%%", optimal.encapsulation_efficiency * 100)
    logger.info("  Diametro: %.1f nm", optimal.particle_diameter_nm)
    logger.info("  Zeta: %.2f mV", optimal.zeta_potential_mv)
    logger.info("  Fagocitavel: %s", optimal.suitable_for_phagocytosis)

    # Resultados detalhados da varredura
    scan_data = [asdict(r) for r in scan]

    return {
        "aso_length_nt": ASO_LENGTH,
        "ps_charges": ASO_PS_CHARGES,
        "scan_range": {"min": 1.0, "max": 20.0, "step": 1.0},
        "scan_results": scan_data,
        "optimal": asdict(optimal),
        "optimal_np_ratio": optimal.np_ratio,
        "optimal_encapsulation": optimal.encapsulation_efficiency,
        "optimal_diameter_nm": optimal.particle_diameter_nm,
        "ps_aso_note": (
            "PS-ASOs require lower N/P (4-8) than mRNA (6-12) due to "
            "permanent negative charge on phosphorothioate backbone."
        ),
    }


# ---------------------------------------------------------------------------
# 3. Cinetica de liberacao pH-responsiva
# ---------------------------------------------------------------------------


def analyze_ph_release() -> dict[str, Any]:
    """Calcula perfil de liberacao pH-responsiva ao longo do pathway intracelular.

    Modelo sigmoidal: f_released(pH) = 1 / (1 + exp(k * (pH - pKa)))
    pKa do DLin-MC3-DMA = 6.44

    Compartimentos:
    - pH 7.4 (sangue): liberacao minima -> LNP estavel em circulacao
    - pH 6.5 (endossomo precoce): liberacao parcial -> inicio da transicao
    - pH 5.0 (endossomo tardio): liberacao quase completa
    - pH 4.5 (fagolisossomo): liberacao maxima -> ASO disponivel para acao

    KEY INSIGHT: a liberacao no fagolisossomo e VANTAJOSA porque
    o alvo (SL RNA de L. infantum) esta no MESMO compartimento.
    Diferente de siRNA/mRNA onde escape endossomal e necessario,
    aqui o ASO PRECISA estar no compartimento acido.

    Returns:
        Dicionario com perfil completo de liberacao.
    """
    logger.info("Analise 3/6: Cinetica de liberacao pH-responsiva")

    profile = compute_release_profile()

    for compartment, release in profile.compartment_releases.items():
        logger.info(
            "  pH %.1f (%s): MC3 protonado %.1f%%, ASO liberado %.1f%% [%s]",
            release.ph,
            compartment,
            release.fraction_protonated * 100,
            release.fraction_released * 100,
            release.release_status,
        )

    logger.info("  INSIGHT: Liberacao no fagolisossomo = %.1f%% -> VANTAJOSO",
                profile.release_at_target * 100)

    # Serializar compartment_releases
    releases_dict = {
        name: asdict(rel)
        for name, rel in profile.compartment_releases.items()
    }

    return {
        "pka_ionizable_lipid": profile.pka_ionizable,
        "hill_slope": profile.hill_slope,
        "compartment_releases": releases_dict,
        "release_at_target_ph": profile.release_at_target,
        "target_compartment": profile.target_compartment,
        "advantage_for_macrophage_delivery": profile.advantage_for_macrophage,
        "mechanism": (
            "DLin-MC3-DMA is neutral at pH 7.4 (stealth) and protonated at pH <6.5. "
            "Protonation triggers interaction with endosomal phospholipids, "
            "forming inverted hexagonal (H_II) phase that disrupts the LNP "
            "and releases the ASO cargo."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Targeting de macrofagos
# ---------------------------------------------------------------------------


def analyze_macrophage_targeting(base_diameter_nm: float) -> dict[str, Any]:
    """Modela efeito de manose-PEG no targeting de macrofagos.

    Macrofagos infectados por L. infantum expressam CD206 (receptor
    de manose) em nivel elevado, especialmente macrofagos M2
    (perfil anti-inflamatorio induzido pelo parasita).

    Substituir parcialmente DMG-PEG2000 por manose-PEG-DSPE
    funcionaliza a superficie da LNP para reconhecimento ativo
    pelo CD206, aumentando captacao seletiva por macrofagos 2-5x.

    Args:
        base_diameter_nm: Diametro da LNP antes da modificacao.

    Returns:
        Dicionario com resultados do targeting.
    """
    logger.info("Analise 4/6: Targeting de macrofagos com manose-PEG")

    # Testar diferentes fracoes de substituicao
    fractions = [0.0, 0.25, 0.50, 0.75, 1.0]
    targeting_results: list[dict[str, Any]] = []

    for frac in fractions:
        result = compute_macrophage_targeting(base_diameter_nm, frac)
        targeting_results.append(asdict(result))
        logger.info(
            "  Manose-PEG %.0f%%: uptake %.1fx, diametro %.1f nm, fagocitavel: %s",
            frac * 100,
            result.uptake_fold_increase,
            result.particle_diameter_nm,
            result.suitable_for_phagocytosis,
        )

    # Recomendacao: 50% substituicao (equilibrio entre targeting e stealth)
    recommended = compute_macrophage_targeting(base_diameter_nm, 0.50)

    logger.info("  RECOMENDADO: 50%% manose-PEG -> %.1fx uptake em macrofagos",
                recommended.uptake_fold_increase)

    return {
        "receptor_target": "CD206 (mannose receptor / MRC1)",
        "biological_rationale": (
            "L. infantum-infected macrophages upregulate CD206 expression. "
            "Mannose-decorated LNPs are recognized by CD206 and internalized "
            "via receptor-mediated endocytosis, directing the ASO payload "
            "to the phagolysosome where the parasite resides."
        ),
        "scan_results": targeting_results,
        "recommended": asdict(recommended),
        "recommended_mannose_fraction": 0.50,
        "expected_uptake_fold": recommended.uptake_fold_increase,
        "mannose_peg_lipid": {
            "name": "Mannose-PEG2000-DSPE",
            "mw_g_per_mol": MW_MANNOSE_PEG,
            "replaces": "DMG-PEG2000",
        },
    }


# ---------------------------------------------------------------------------
# 5. Estabilidade e armazenamento
# ---------------------------------------------------------------------------


def analyze_storage_stability() -> dict[str, Any]:
    """Calcula estabilidade da LNP em tres temperaturas.

    Temperaturas testadas:
    - -20 C: armazenamento ultra-frio (padrao para reagentes)
    - 4 C: refrigerado (cadeia fria padrao)
    - 25 C: ambiente (relevante para Brasil — regioes endemicas de LV
      no Nordeste frequentemente nao tem infraestrutura de cadeia fria)

    Liofilizacao e avaliada como estrategia para eliminar a dependencia
    de cadeia fria em regioes tropicais.

    Returns:
        Dicionario com perfis de estabilidade.
    """
    logger.info("Analise 5/6: Estabilidade em armazenamento")

    profiles = compute_storage_profiles()
    stability_data: list[dict[str, Any]] = []

    for s in profiles:
        stability_data.append(asdict(s))
        logger.info(
            "  %+.0f C: k_agg = %.8f/dia, t1/2 = %.1f dias, "
            "f(90d) = %.2f%%, liofilizacao: %s",
            s.temperature_celsius,
            s.k_aggregation_per_day,
            s.half_life_days,
            s.fraction_intact_90_days * 100,
            s.lyophilization_feasible,
        )

    # Avaliar viabilidade para contexto brasileiro
    storage_25 = profiles[2]  # 25 C
    shelf_life_25c_adequate = storage_25.fraction_intact_90_days >= 0.80

    logger.info(
        "  Brasil: f(90d) a 25 C = %.2f%% -> %s",
        storage_25.fraction_intact_90_days * 100,
        "ADEQUADO" if shelf_life_25c_adequate else "REQUER LIOFILIZACAO",
    )

    return {
        "profiles": stability_data,
        "temperatures_tested": [-20.0, 4.0, 25.0],
        "shelf_life_25c_adequate": shelf_life_25c_adequate,
        "lyophilization": {
            "feasible": True,
            "cryoprotectants": ["sucrose (10% w/v)", "trehalose (10% w/v)"],
            "reconstitution": "5 min in sterile water, vortex 30 sec",
            "post_reconstitution_size_increase": "<20%",
            "post_reconstitution_ee": ">80%",
            "advantage_for_brazil": (
                "Eliminates cold chain dependency for deployment in "
                "Northeast Brazil (Bahia, Ceara, Maranhao) where "
                "visceral leishmaniasis is endemic and cold chain "
                "infrastructure is limited."
            ),
        },
    }


# ---------------------------------------------------------------------------
# 6. Comparacao LNP vs ASO livre
# ---------------------------------------------------------------------------


def analyze_lnp_vs_naked(
    encapsulation_efficiency: float,
    uptake_fold: float,
    release_at_target: float,
) -> dict[str, Any]:
    """Compara LNP-ASO com ASO livre (naked) em multiplas dimensoes.

    A comparacao cobre: encapsulacao, uptake celular, custo, armazenamento,
    e complexidade de fabricacao. O objetivo e avaliar se o beneficio
    terapeutico da LNP justifica a complexidade adicional.

    Args:
        encapsulation_efficiency: EE da LNP otima (0-1).
        uptake_fold: Fator de aumento de uptake por manose.
        release_at_target: Fracao liberada no fagolisossomo.

    Returns:
        Dicionario com comparacao ponto-a-ponto.
    """
    logger.info("Analise 6/6: Comparacao LNP vs ASO livre")

    comparison = {
        "naked_aso": {
            "encapsulation_efficiency": None,
            "cellular_uptake_relative": 1.0,
            "target_delivery": "passive (gymnosis or electroporation)",
            "nuclease_protection": "PS backbone only",
            "estimated_cost_per_dose_usd": COST_NAKED_ASO_USD,
            "storage_requirement": "Lyophilized powder, room temperature",
            "manufacturing_complexity": COMPLEXITY_NAKED,
            "manufacturing_description": "Solid-phase synthesis + HPLC purification",
            "advantages": [
                "Lower cost per dose",
                "Simpler manufacturing",
                "No cold chain needed (lyophilized)",
                "PS backbone provides nuclease resistance",
            ],
            "disadvantages": [
                "Low cellular uptake without delivery vehicle",
                "No macrophage targeting",
                "Requires high dose to compensate for poor uptake",
                "Limited biodistribution control",
            ],
        },
        "lnp_aso": {
            "encapsulation_efficiency": encapsulation_efficiency,
            "cellular_uptake_relative": uptake_fold,
            "target_delivery": "mannose receptor-mediated endocytosis",
            "nuclease_protection": "LNP encapsulation + PS backbone",
            "estimated_cost_per_dose_usd": COST_LNP_ASO_USD,
            "storage_requirement": "4 C (liquid) or lyophilized",
            "manufacturing_complexity": COMPLEXITY_LNP,
            "manufacturing_description": (
                "ASO synthesis + microfluidic LNP formulation + "
                "characterization (DLS, EE assay, zeta potential)"
            ),
            "advantages": [
                f"{uptake_fold:.1f}x higher macrophage uptake (mannose-targeted)",
                f"{encapsulation_efficiency*100:.0f}% encapsulation efficiency",
                f"{release_at_target*100:.0f}% release at phagolysosomal pH",
                "Dual protection: LNP shell + PS backbone",
                "Lower effective dose (higher delivery efficiency)",
                "Active targeting reduces off-target effects",
            ],
            "disadvantages": [
                "Higher cost per dose (3x vs naked)",
                "More complex manufacturing (microfluidics)",
                "Requires cold chain or lyophilization",
                "Additional QC steps (size, PDI, EE, zeta)",
            ],
        },
    }

    # Calcular beneficio liquido
    effective_dose_reduction = uptake_fold * encapsulation_efficiency
    cost_ratio = COST_LNP_ASO_USD / COST_NAKED_ASO_USD

    net_benefit = {
        "uptake_improvement_fold": uptake_fold,
        "encapsulation_efficiency": encapsulation_efficiency,
        "effective_dose_reduction_fold": round(effective_dose_reduction, 2),
        "cost_ratio_lnp_vs_naked": cost_ratio,
        "cost_effectiveness_ratio": round(effective_dose_reduction / cost_ratio, 2),
        "recommendation": (
            "LNP formulation is RECOMMENDED for MRL-ASO-001 because: "
            f"(1) {uptake_fold:.1f}x higher macrophage uptake via mannose-CD206, "
            f"(2) {release_at_target*100:.0f}% release at phagolysosomal pH 4.5, "
            f"(3) effective dose reduction of {effective_dose_reduction:.1f}x "
            f"justifies the {cost_ratio:.0f}x cost increase. "
            f"Cost-effectiveness ratio: {effective_dose_reduction/cost_ratio:.2f}."
        ),
    }

    logger.info("  ASO livre: uptake = 1.0x, custo = $%.2f", COST_NAKED_ASO_USD)
    logger.info("  LNP-ASO: uptake = %.1fx, custo = $%.2f", uptake_fold, COST_LNP_ASO_USD)
    logger.info("  Reducao efetiva de dose: %.1fx", effective_dose_reduction)
    logger.info("  Razao custo-eficacia: %.2f", effective_dose_reduction / cost_ratio)

    return {
        "comparison": comparison,
        "net_benefit": net_benefit,
    }


# ---------------------------------------------------------------------------
# Geracao de relatorio Markdown
# ---------------------------------------------------------------------------


def generate_markdown_report(results: dict[str, Any]) -> str:
    """Gera relatorio Markdown com os resultados completos.

    Estrutura: sumario executivo + 6 secoes de analise + conclusao +
    referencias. Tabelas formatadas para leitura tecnica.

    Args:
        results: Envelope completo de resultados.

    Returns:
        String com conteudo Markdown do relatorio.
    """
    lines: list[str] = []

    lines.append("# Module D: LNP Formulation Analysis for MRL-ASO-001")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"**Module version:** {MODULE_VERSION}")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(results["executive_summary"])
    lines.append("")

    # ---- Secao 1: LNP Composition ----
    comp = results["lnp_composition"]
    lines.append("## 1. LNP Composition")
    lines.append("")
    lines.append("Standard 4-component formulation based on Onpattro (patisiran).")
    lines.append("")
    lines.append("| Component | Name | Mol% | MW (g/mol) | Role |")
    lines.append("|---|---|---|---|---|")
    for key in ["ionizable_lipid", "helper_lipid", "sterol", "peg_lipid"]:
        c = comp["components"][key]
        lines.append(
            f"| {key.replace('_', ' ').title()} | {c['name']} | "
            f"{c['mol_fraction']*100:.1f}% | {c['mw_g_per_mol']:.1f} | {c['role']} |"
        )
    lines.append("")
    c_data = comp["composition"]
    lines.append(f"**Weighted average MW:** {c_data['weighted_average_mw']:.2f} g/mol")
    lines.append(f"**Charge at pH 7.4:** {c_data['charge_ph_7_4']:.6f} (near-neutral, stealth)")
    lines.append(f"**Charge at pH 4.5:** {c_data['charge_ph_4_5']:.6f} (cationic, triggers release)")
    lines.append(f"**pKa (DLin-MC3-DMA):** {c_data['pka_ionizable']}")
    lines.append("")

    # ---- Secao 2: N/P Optimization ----
    np_data = results["np_optimization"]
    lines.append("## 2. N/P Ratio Optimization")
    lines.append("")
    lines.append(f"ASO: {np_data['aso_length_nt']} nt, {np_data['ps_charges']} PS charges (negative).")
    lines.append(f"Optimal N/P range for PS-ASOs: 4-8 (vs 6-12 for mRNA).")
    lines.append("")
    lines.append("| N/P | EE (%) | Diameter (nm) | Zeta (mV) | Phagocytable | Optimal |")
    lines.append("|---|---|---|---|---|---|")
    for r in np_data["scan_results"]:
        lines.append(
            f"| {r['np_ratio']:.0f} | {r['encapsulation_efficiency']*100:.1f} | "
            f"{r['particle_diameter_nm']:.1f} | {r['zeta_potential_mv']:.1f} | "
            f"{'Yes' if r['suitable_for_phagocytosis'] else 'No'} | "
            f"{'**Yes**' if r['optimal'] else 'No'} |"
        )
    lines.append("")
    opt = np_data["optimal"]
    lines.append(f"**Selected N/P:** {opt['np_ratio']:.0f}")
    lines.append(f"**Encapsulation efficiency:** {opt['encapsulation_efficiency']*100:.1f}%")
    lines.append(f"**Particle diameter:** {opt['particle_diameter_nm']:.1f} nm")
    lines.append(f"**Zeta potential:** {opt['zeta_potential_mv']:.1f} mV")
    lines.append("")

    # ---- Secao 3: pH-Responsive Release ----
    rel = results["ph_release"]
    lines.append("## 3. pH-Responsive Release Kinetics")
    lines.append("")
    lines.append(f"Model: f_released(pH) = 1 / (1 + exp({rel['hill_slope']:.0f} * (pH - {rel['pka_ionizable_lipid']:.2f})))")
    lines.append("")
    lines.append("| Compartment | pH | MC3 Protonated (%) | ASO Released (%) | Status | Transit (min) |")
    lines.append("|---|---|---|---|---|---|")
    for name, r in rel["compartment_releases"].items():
        lines.append(
            f"| {r['compartment']} | {r['ph']} | {r['fraction_protonated']*100:.1f} | "
            f"{r['fraction_released']*100:.1f} | {r['release_status']} | {r['transit_time_min']:.0f} |"
        )
    lines.append("")
    lines.append(f"**Release at target (phagolysosome):** {rel['release_at_target_ph']*100:.1f}%")
    lines.append("")
    lines.append(f"**Key insight:** {rel['advantage_for_macrophage_delivery']}")
    lines.append("")

    # ---- Secao 4: Macrophage Targeting ----
    tgt = results["macrophage_targeting"]
    lines.append("## 4. Macrophage-Targeted LNP (Mannose-PEG)")
    lines.append("")
    lines.append(f"**Target receptor:** {tgt['receptor_target']}")
    lines.append("")
    lines.append(f"{tgt['biological_rationale']}")
    lines.append("")
    lines.append("| Mannose-PEG (%) | Uptake Fold | Diameter (nm) | Phagocytable |")
    lines.append("|---|---|---|---|")
    for r in tgt["scan_results"]:
        lines.append(
            f"| {r['mannose_peg_fraction']*100:.0f} | {r['uptake_fold_increase']:.1f}x | "
            f"{r['particle_diameter_nm']:.1f} | {'Yes' if r['suitable_for_phagocytosis'] else 'No'} |"
        )
    lines.append("")
    rec = tgt["recommended"]
    lines.append(f"**Recommended:** {tgt['recommended_mannose_fraction']*100:.0f}% mannose-PEG substitution")
    lines.append(f"**Expected uptake increase:** {rec['uptake_fold_increase']:.1f}x")
    lines.append("")

    # ---- Secao 5: Storage Stability ----
    stab = results["storage_stability"]
    lines.append("## 5. Storage Stability")
    lines.append("")
    lines.append("Arrhenius model: k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T))")
    lines.append("")
    lines.append("| Temp (C) | k_agg (1/day) | Half-life (days) | f(30d) | f(90d) | f(365d) |")
    lines.append("|---|---|---|---|---|---|")
    for p in stab["profiles"]:
        lines.append(
            f"| {p['temperature_celsius']:+.0f} | {p['k_aggregation_per_day']:.6f} | "
            f"{p['half_life_days']:.1f} | {p['fraction_intact_30_days']*100:.1f}% | "
            f"{p['fraction_intact_90_days']*100:.1f}% | {p['fraction_intact_365_days']*100:.1f}% |"
        )
    lines.append("")
    lines.append(f"**25 C shelf life adequate (>80% at 90 days):** "
                 f"{'Yes' if stab['shelf_life_25c_adequate'] else 'No'}")
    lines.append("")
    lyo = stab["lyophilization"]
    lines.append("### Lyophilization for Tropical Deployment")
    lines.append("")
    lines.append(f"- **Feasible:** {'Yes' if lyo['feasible'] else 'No'}")
    lines.append(f"- **Cryoprotectants:** {', '.join(lyo['cryoprotectants'])}")
    lines.append(f"- **Reconstitution:** {lyo['reconstitution']}")
    lines.append(f"- **Size increase post-reconstitution:** {lyo['post_reconstitution_size_increase']}")
    lines.append(f"- **EE post-reconstitution:** {lyo['post_reconstitution_ee']}")
    lines.append(f"- **Brazil advantage:** {lyo['advantage_for_brazil']}")
    lines.append("")

    # ---- Secao 6: LNP vs Naked ASO ----
    cmp = results["lnp_vs_naked"]
    lines.append("## 6. LNP vs Naked ASO Comparison")
    lines.append("")
    nb = cmp["net_benefit"]
    lines.append("| Property | Naked ASO | LNP-ASO |")
    lines.append("|---|---|---|")

    naked = cmp["comparison"]["naked_aso"]
    lnp = cmp["comparison"]["lnp_aso"]

    rows = [
        ("Cellular uptake (relative)", f"{naked['cellular_uptake_relative']:.1f}x", f"{lnp['cellular_uptake_relative']:.1f}x"),
        ("Encapsulation efficiency", "N/A", f"{lnp['encapsulation_efficiency']*100:.0f}%"),
        ("Target delivery", naked["target_delivery"], lnp["target_delivery"]),
        ("Nuclease protection", naked["nuclease_protection"], lnp["nuclease_protection"]),
        ("Cost per dose (USD)", f"${naked['estimated_cost_per_dose_usd']:.2f}", f"${lnp['estimated_cost_per_dose_usd']:.2f}"),
        ("Storage", naked["storage_requirement"], lnp["storage_requirement"]),
        ("Manufacturing complexity", f"{naked['manufacturing_complexity']}/5", f"{lnp['manufacturing_complexity']}/5"),
    ]

    for label, v1, v2 in rows:
        lines.append(f"| {label} | {v1} | {v2} |")

    lines.append("")
    lines.append(f"**Effective dose reduction:** {nb['effective_dose_reduction_fold']:.1f}x")
    lines.append(f"**Cost ratio (LNP/naked):** {nb['cost_ratio_lnp_vs_naked']:.0f}x")
    lines.append(f"**Cost-effectiveness ratio:** {nb['cost_effectiveness_ratio']:.2f}")
    lines.append("")
    lines.append(f"**Recommendation:** {nb['recommendation']}")
    lines.append("")

    # ---- Overall Conclusion ----
    lines.append("## Overall Conclusion")
    lines.append("")
    lines.append(results["overall_conclusion"])
    lines.append("")
    lines.append("## References")
    lines.append("")
    lines.append("1. Cullis PR, Hope MJ (2017) Mol Ther 25(7):1467-1475 — LNP design principles")
    lines.append("2. Jayaraman M et al. (2012) Angew Chem 51(34):8529-8533 — DLin-MC3-DMA optimization")
    lines.append("3. Semple SC et al. (2010) Nat Biotechnol 28(2):172-176 — LNP encapsulation")
    lines.append("4. Patel S et al. (2020) J Control Release 327:146-160 — Macrophage targeting with mannose")
    lines.append("5. Kulkarni JA et al. (2018) Nanoscale 10(10):4567-4573 — pH-responsive release mechanism")
    lines.append("6. Sahay G et al. (2013) Nat Biotechnol 31(7):653-658 — Endosomal escape and trafficking")
    lines.append("7. Schoenmaker L et al. (2021) Int J Pharm 601:120586 — LNP storage stability")
    lines.append("8. Ball RL et al. (2017) Drug Deliv Transl Res 7(1):89-103 — LNP lyophilization")
    lines.append("9. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Orquestrador principal
# ---------------------------------------------------------------------------


def main() -> dict[str, Any]:
    """Executa todas as analises de formulacao LNP.

    Fluxo:
        1. Composicao lipidica (4 componentes)
        2. Otimizacao N/P (varredura 1-20)
        3. Liberacao pH-responsiva (4 compartimentos)
        4. Targeting de macrofagos (manose-PEG)
        5. Estabilidade em armazenamento (3 temperaturas)
        6. Comparacao LNP vs ASO livre
        7. Gerar relatorio Markdown
        8. Gravar JSON + Markdown

    Returns:
        Envelope completo com todos os resultados.
    """
    logger.info("=" * 70)
    logger.info("MODULO D: Nanoparticulas Lipidicas (LNP) — MRL-ASO-001")
    logger.info("=" * 70)
    logger.info("Pergunta: LNP melhora a entrega de MRL-ASO-001 a macrofagos?")
    logger.info("")

    start_time = time.time()

    # --- Analise 1: Composicao ---
    lnp_composition = analyze_lnp_composition()
    logger.info("")

    # --- Analise 2: N/P ---
    np_optimization = analyze_np_optimization()
    logger.info("")

    # --- Analise 3: Liberacao pH ---
    ph_release = analyze_ph_release()
    logger.info("")

    # --- Analise 4: Targeting ---
    optimal_diameter = np_optimization["optimal_diameter_nm"]
    macrophage_targeting = analyze_macrophage_targeting(optimal_diameter)
    logger.info("")

    # --- Analise 5: Armazenamento ---
    storage_stability = analyze_storage_stability()
    logger.info("")

    # --- Analise 6: Comparacao ---
    lnp_vs_naked = analyze_lnp_vs_naked(
        encapsulation_efficiency=np_optimization["optimal_encapsulation"],
        uptake_fold=macrophage_targeting["expected_uptake_fold"],
        release_at_target=ph_release["release_at_target_ph"],
    )
    logger.info("")

    # --- Compilar resultados ---
    elapsed = round(time.time() - start_time, 2)

    # Criterios de sucesso:
    # 1. Encapsulacao >80%
    # 2. Liberacao >90% no fagolisossomo
    # 3. Tamanho 50-200 nm (fagocitavel)
    # 4. Uptake de macrofagos >2x
    ee = np_optimization["optimal_encapsulation"]
    release = ph_release["release_at_target_ph"]
    diameter = macrophage_targeting["recommended"]["particle_diameter_nm"]
    uptake = macrophage_targeting["expected_uptake_fold"]

    ee_pass = ee >= 0.80
    release_pass = release >= 0.90
    size_pass = 50.0 <= diameter <= 200.0
    uptake_pass = uptake >= 2.0
    all_pass = ee_pass and release_pass and size_pass and uptake_pass

    # Sumario executivo
    if all_pass:
        executive_summary = (
            f"LNP encapsulation is RECOMMENDED for MRL-ASO-001 delivery to macrophages. "
            f"Key metrics: encapsulation efficiency = {ee*100:.0f}% (N/P {np_optimization['optimal_np_ratio']:.0f}), "
            f"phagolysosomal release = {release*100:.0f}%, "
            f"particle diameter = {diameter:.0f} nm (phagocytable), "
            f"mannose-targeted uptake = {uptake:.1f}x vs untargeted. "
            f"The pH-responsive release of DLin-MC3-DMA (pKa 6.44) "
            f"delivers the ASO cargo directly to the phagolysosome "
            f"where L. infantum amastigotes and their SL RNA reside. "
            f"Lyophilization enables deployment in endemic regions of Brazil "
            f"without cold chain infrastructure."
        )
    else:
        failures: list[str] = []
        if not ee_pass:
            failures.append(f"encapsulation ({ee*100:.0f}% < 80%)")
        if not release_pass:
            failures.append(f"release ({release*100:.0f}% < 90%)")
        if not size_pass:
            failures.append(f"size ({diameter:.0f} nm outside 50-200 nm)")
        if not uptake_pass:
            failures.append(f"uptake ({uptake:.1f}x < 2.0x)")
        executive_summary = (
            f"LNP formulation does not meet all criteria: "
            f"{', '.join(failures)}. Optimization needed."
        )

    # Conclusao geral
    if all_pass:
        overall_conclusion = (
            f"LNP encapsulation provides significant advantages for MRL-ASO-001 delivery:\n\n"
            f"1. **Encapsulation:** {ee*100:.0f}% efficiency at N/P {np_optimization['optimal_np_ratio']:.0f}, "
            f"protecting the ASO during circulation.\n\n"
            f"2. **pH-responsive release:** {release*100:.0f}% of cargo released at phagolysosomal pH 4.5. "
            f"Unlike mRNA/siRNA therapeutics that require endosomal escape, MRL-ASO-001 benefits from "
            f"phagolysosomal release because the target (SL RNA) is in the SAME compartment.\n\n"
            f"3. **Macrophage targeting:** Mannose-PEG-DSPE decoration provides {uptake:.1f}x higher "
            f"uptake by infected macrophages via CD206 receptor-mediated endocytosis.\n\n"
            f"4. **Size:** {diameter:.0f} nm diameter is within the optimal range for "
            f"phagocytic uptake (50-200 nm).\n\n"
            f"5. **Storage:** Lyophilization with sucrose/trehalose enables room-temperature "
            f"storage, critical for deployment in endemic areas of Northeast Brazil.\n\n"
            f"6. **Cost-effectiveness:** Despite 3x higher cost per dose, the {uptake:.1f}x "
            f"uptake improvement and {ee*100:.0f}% encapsulation yield a cost-effectiveness "
            f"ratio of {lnp_vs_naked['net_benefit']['cost_effectiveness_ratio']:.2f}, "
            f"justifying the LNP formulation.\n\n"
            f"**The LNP-formulated MRL-ASO-001 is predicted to be a viable delivery strategy "
            f"for targeting L. infantum in canine macrophages.**"
        )
    else:
        overall_conclusion = (
            f"LNP formulation does not meet all criteria. "
            f"Failed: {', '.join(failures)}. "
            f"Further optimization of formulation parameters is required."
        )

    # Montar envelope final
    results: dict[str, Any] = {
        "module": MODULE_NAME,
        "version": MODULE_VERSION,
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": elapsed,
        "status": "success" if all_pass else "partial_failure",
        "all_criteria_met": all_pass,
        "criteria": {
            "encapsulation_ge_80pct": ee_pass,
            "release_ge_90pct": release_pass,
            "size_50_200nm": size_pass,
            "uptake_ge_2x": uptake_pass,
        },
        "executive_summary": executive_summary,
        "overall_conclusion": overall_conclusion,
        "aso_specifications": {
            "name": "MRL-ASO-001",
            "sequence_length_nt": ASO_LENGTH,
            "gapmer_design": "5-15-5 (LNA-DNA-LNA)",
            "backbone": "full phosphorothioate (PS)",
            "ps_charges": ASO_PS_CHARGES,
        },
        "lnp_composition": lnp_composition,
        "np_optimization": np_optimization,
        "ph_release": ph_release,
        "macrophage_targeting": macrophage_targeting,
        "storage_stability": storage_stability,
        "lnp_vs_naked": lnp_vs_naked,
    }

    # --- Gravar resultados ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = RESULTS_DIR / "module_d_lnp.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    logger.info("Resultados JSON gravados em: %s", json_path)

    # Markdown
    report_md = generate_markdown_report(results)
    md_path = RESULTS_DIR / "module_d_lnp_report.md"
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Relatorio Markdown gravado em: %s", md_path)

    # --- Sumario final ---
    logger.info("")
    logger.info("=" * 70)
    logger.info("RESULTADO FINAL: %s", "APROVADO" if all_pass else "REPROVADO")
    logger.info("=" * 70)
    logger.info("  Encapsulacao: %.1f%% (limiar: 80%%)", ee * 100)
    logger.info("  Liberacao no fagolisossomo: %.1f%% (limiar: 90%%)", release * 100)
    logger.info("  Diametro: %.0f nm (faixa: 50-200 nm)", diameter)
    logger.info("  Uptake macrofago: %.1fx (limiar: 2.0x)", uptake)
    logger.info("  Tempo de execucao: %.2f segundos", elapsed)

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
