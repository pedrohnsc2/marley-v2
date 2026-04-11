"""Plataforma A -- Otimizacao da vacina de mRNA contra LVC.

Ponto de entrada principal que orquestra os quatro sub-modulos:
    A1. Analise de redundancia de epitopos (epitope_analysis.py)
    A2. Otimizacao de pseudouridina (pseudouridine.py)
    A3. Formulacao de LNP veterinaria (lnp_formulation.py)
    A4. Analise estrategica de financiamento (strategic.py)

A plataforma A recebe o construto vacinal de mRNA existente (335 aa,
11 epitopos unicos) e otimiza para producao e entrega em caes.

Saidas:
    - results/platform_a_results.json -- dados estruturados completos
    - results/platform_a_report.md -- relatorio legivel em Markdown

Usage:
    python -m vaccine_platforms.platform_a_mrna.run
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Final

from core.logger import get_logger

from . import epitope_analysis as A1
from . import pseudouridine as A2
from . import lnp_formulation as A3
from . import strategic as A4

logger = get_logger("platform_a_mrna")

# ---------------------------------------------------------------------------
# Caminhos
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent.parent
CONSTRUCT_CARD_PATH: Final[Path] = (
    PROJECT_ROOT / "results" / "construct" / "construct_card.json"
)
OUTPUT_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_construct_card() -> dict:
    """Le o construct_card.json com dados do construto mRNA.

    Este arquivo contem os epitopos selecionados pelo pipeline de
    vacinologia reversa e os parametros biofisicos do construto original.

    Returns:
        Dicionario com os dados do construto.

    Raises:
        FileNotFoundError: Se o arquivo nao existe.
    """
    if not CONSTRUCT_CARD_PATH.exists():
        raise FileNotFoundError(
            f"construct_card.json nao encontrado em {CONSTRUCT_CARD_PATH}. "
            f"Execute o pipeline de construcao mRNA primeiro."
        )

    with open(CONSTRUCT_CARD_PATH, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _format_sequence(seq: str, width: int = 70) -> str:
    """Formata uma sequencia em linhas de largura fixa.

    Args:
        seq: Sequencia (proteica ou DNA).
        width: Largura maxima por linha.

    Returns:
        String formatada com quebras de linha.
    """
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def _generate_markdown_report(results: dict) -> str:
    """Gera o relatorio completo em Markdown para a plataforma A.

    O relatorio integra os resultados dos quatro sub-modulos (A1-A4)
    em um documento unico, formatado para leitura e apresentacao.

    Args:
        results: Dicionario completo de resultados.

    Returns:
        String com o relatorio formatado em Markdown.
    """
    a1 = results["epitope_redundancy"]
    a2 = results["pseudouridine_optimization"]
    a3 = results["lnp_formulation"]
    a4 = results["strategic_analysis"]

    report = f"""# Platform A -- mRNA Vaccine Optimization

**Generated:** {results['metadata']['generated']}
**Pipeline version:** {results['metadata']['pipeline_version']}
**Construct:** {results['metadata']['source_construct_length_aa']} aa, {results['metadata']['source_mrna_length_nt']} nt

---

## 1. Epitope Redundancy Analysis (A1)

### Overview

| Metric | Value |
|---|---|
| Total canonical epitopes | {a1['total_epitopes']} |
| Unique source genes | {a1['unique_source_genes']} |
| Unique DLA alleles | {a1['unique_alleles']} |
| Pairwise comparisons | {a1['pairwise_comparisons_count']} |
| Redundant pairs (>70% sim + same allele) | {a1['redundant_pairs_count']} |
| Minimum epitopes for 95% coverage | {a1['minimum_epitopes_95pct_coverage']} |

### Top 10 Most Similar Pairs

| Epitope A | Epitope B | Similarity | Same Gene | Same Allele |
|---|---|---|---|---|
"""
    for pair in a1["top_similar_pairs"][:10]:
        report += (
            f"| `{pair['epitope_a']}` | `{pair['epitope_b']}` | "
            f"{pair['hamming_similarity']:.2%} | "
            f"{'Yes' if pair['same_gene'] else 'No'} | "
            f"{'Yes' if pair['same_allele'] else 'No'} |\n"
        )

    report += "\n### Construct Variants\n\n"
    report += "| Epitopes | Length (aa) | MW (kDa) | II | Stable | Antigenicity | Coverage |\n"
    report += "|---|---|---|---|---|---|---|\n"

    for v in a1["construct_variants"]:
        report += (
            f"| {v['epitope_count']} | {v['protein_length_aa']} | "
            f"{v['molecular_weight_da']/1000:.1f} | "
            f"{v['instability_index']:.1f} | "
            f"{'Yes' if v['is_stable'] else 'No'} | "
            f"{v['antigenicity_proxy']:.4f} | "
            f"{v['theoretical_coverage_pct']:.1f}% |\n"
        )

    report += f"\n**Recommendation:** {a1['recommendation']}\n"

    # --- A2: Pseudouridina ---
    report += f"""
---

## 2. Pseudouridine Optimization (A2)

### mRNA Characteristics

| Metric | Value |
|---|---|
| mRNA length | {a2['mrna_length_nt']} nt |
| Total uridines | {a2['total_uridines']} |
| Uridine fraction | {a2['uridine_fraction']:.1%} |
| GC content | {a2['gc_content']:.1%} |

### Uridine Structural Context

| Category | Count |
|---|---|
| High base-pair score (>0.5) | {a2['uridine_context_summary']['high_bp_score_count']} |
| Medium base-pair score (0.25-0.5) | {a2['uridine_context_summary']['medium_bp_score_count']} |
| Low base-pair score (<0.25) | {a2['uridine_context_summary']['low_bp_score_count']} |
| Average base-pair score | {a2['uridine_context_summary']['avg_base_pair_score']:.4f} |

### Substitution Strategies

| Strategy | U Substituted | Translation Fold | Immunogenicity Reduction | Cost Factor |
|---|---|---|---|---|
"""
    for s in a2["strategies"]:
        report += (
            f"| {s['label']} | {s['uridines_substituted']}/{s['uridines_total']} | "
            f"{s['estimated_translation_fold']:.2f}x | "
            f"{s['estimated_immunogenicity_reduction_pct']:.1f}% | "
            f"{s['estimated_cost_factor']:.2f}x |\n"
        )

    report += f"\n**Optimal fraction:** {a2['optimal_fraction']:.0%}\n\n"
    report += f"**Rationale:** {a2['optimal_rationale']}\n"

    # --- A3: LNP ---
    dose = a3["dose"]
    comp = a3["lnp_composition"]
    report += f"""
---

## 3. Veterinary LNP Formulation (A3)

### Dose Calculation

| Parameter | Value |
|---|---|
| Target species | {dose['target_species']} |
| Body weight | {dose['body_weight_kg']} kg |
| Dose per kg | {dose['dose_ug_per_kg']} ug/kg |
| Total dose per injection | {dose['total_dose_ug']} ug |
| Doses per animal | {dose['doses_per_animal']} |
| Total mRNA per animal | {dose['total_mrna_per_animal_ug']} ug |

### LNP Composition (per dose)

| Component | Mol% | Amount |
|---|---|---|
| SM-102 (ionizable lipid) | {comp['molar_ratios_pct']['ionizable_lipid_SM102']}% | {comp['ionizable_lipid_nmol']:.2f} nmol |
| DSPC | {comp['molar_ratios_pct']['DSPC']}% | {comp['dspc_nmol']:.2f} nmol |
| Cholesterol | {comp['molar_ratios_pct']['cholesterol']}% | {comp['cholesterol_nmol']:.2f} nmol |
| PEG-DMG 2000 | {comp['molar_ratios_pct']['PEG_DMG_2000']}% | {comp['peg_lipid_nmol']:.2f} nmol |

| Metric | Value |
|---|---|
| N/P ratio | {comp['np_ratio']} |
| Total lipid mass | {comp['total_lipid_mass_ug']:.2f} ug |
| Lipid:mRNA ratio (w/w) | {comp['lipid_to_mrna_ratio_w_w']} |
| Encapsulation efficiency | {comp['encapsulation_efficiency_pct']}% |

### N/P Ratio Optimization

| N/P | Lipid Mass (ug) | Encapsulation | Transfection Score | Optimal |
|---|---|---|---|---|
"""
    for p in a3["np_ratio_optimization"]:
        report += (
            f"| {p['np_ratio']:.0f} | {p['total_lipid_mass_ug']:.2f} | "
            f"{p['estimated_encapsulation_pct']:.1f}% | "
            f"{p['estimated_transfection_score']:.3f} | "
            f"{'**Yes**' if p['is_optimal'] else ''} |\n"
        )

    report += "\n### Storage Options\n\n"
    report += "| Condition | Temperature | Shelf Life | Feasibility |\n"
    report += "|---|---|---|---|\n"

    for st in a3["storage_options"]:
        report += (
            f"| {st['label']} | {st['temperature_c']} C | "
            f"{st['shelf_life_months']} months | {st['feasibility']} |\n"
        )

    report += f"""
### Cost Estimate

| Component | Laboratory | Industrial |
|---|---|---|
| IVT (mRNA production) | ${a3['cost_laboratory']['ivt_usd']:.2f} | ${a3['cost_industrial']['ivt_usd']:.4f} |
| Lipids | ${a3['cost_laboratory']['lipids_usd']:.4f} | ${a3['cost_industrial']['lipids_usd']:.6f} |
| Formulation | ${a3['cost_laboratory']['formulation_usd']:.2f} | ${a3['cost_industrial']['formulation_usd']:.2f} |
| QC | ${a3['cost_laboratory']['qc_usd']:.2f} | ${a3['cost_industrial']['qc_usd']:.2f} |
| **Total per dose** | **${a3['cost_laboratory']['total_per_dose_usd']:.2f}** | **${a3['cost_industrial']['total_per_dose_usd']:.4f}** |
| **Total per animal** | **${a3['cost_laboratory']['total_per_animal_usd']:.2f}** | **${a3['cost_industrial']['total_per_animal_usd']:.4f}** |

**Recommendation:** {a3['recommendation']}
"""

    # --- A4: Estrategia ---
    mkt = a4["market"]
    report += f"""
---

## 4. Strategic Funding & Market Analysis (A4)

### Market Size (Brazil)

| Metric | Value |
|---|---|
| Total dogs in Brazil | {mkt['total_dogs_brazil']:,} |
| Dogs at risk (endemic areas) | {mkt['dogs_at_risk']:,} ({mkt['dogs_at_risk_fraction_pct']}%) |
| Mean seroprevalence | {mkt['seroprevalence_mean_pct']}% |
| Estimated infected dogs | {mkt['estimated_infected_dogs']:,} |
| Addressable market | {mkt['addressable_market_dogs']:,} dogs/year |
| Market value (BRL/year) | R$ {mkt['market_value_brl_year']:,.0f} |
| Market value (USD/year) | $ {mkt['market_value_usd_year']:,.0f} |
| Human cases/year | {mkt['human_cases_per_year']:,} (CFR: {mkt['human_fatality_rate_pct']}%) |
| Government control cost | R$ {mkt['annual_control_cost_brl']:,.0f}/year |

### Vaccine Sovereignty Argument

**{a4['sovereignty_argument']['title']}**

{a4['sovereignty_argument']['context']}

**COVID-19 precedent:** {a4['sovereignty_argument']['covid19_precedent']}

**One Health:** {a4['sovereignty_argument']['one_health']}

### Funding Sources (ranked by score)

| Rank | Source | Score | Probability | Amount (BRL) | Timeline |
|---|---|---|---|---|---|
"""
    for i, name in enumerate(a4["funding_ranked_by_score"], 1):
        fs = next(f for f in a4["funding_sources"] if f["name"] == name)
        report += (
            f"| {i} | {fs['name']} | {fs['score']:.1f}/10 | "
            f"{fs['probability_pct']:.0f}% | "
            f"R$ {fs['typical_amount_brl_min']/1e6:.1f}-{fs['typical_amount_brl_max']/1e6:.1f}M | "
            f"{fs['timeline_months_min']}-{fs['timeline_months_max']} months |\n"
        )

    report += f"""
### Regulatory Pathways

**Brazil (MAPA):** {a4['regulatory_brazil']['total_estimated_timeline']}
**USA (USDA-CVB):** {a4['regulatory_usda']['total_estimated_timeline']}

### Strategic Recommendation

{a4['recommendation']}

---

*Report generated by the Marley reverse vaccinology pipeline -- Platform A (mRNA optimization)*
"""
    return report


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> dict:
    """Executa o pipeline completo da plataforma A (mRNA optimization).

    Etapas:
        1. Le o construct_card.json com dados do construto
        2. Analise de redundancia de epitopos (A1)
        3. Otimizacao de pseudouridina (A2)
        4. Formulacao de LNP veterinaria (A3)
        5. Analise estrategica (A4)
        6. Gera relatorio e salva resultados

    Returns:
        Dicionario completo com todos os resultados.
    """
    logger.info("=" * 60)
    logger.info("Platform A -- mRNA Vaccine Optimization")
    logger.info("=" * 60)

    # ----- Carregar dados de entrada -----
    logger.info("Loading construct_card.json...")
    card = _load_construct_card()
    logger.info(
        "Construct: %d aa, %d nt, %d epitope entries",
        card.get("protein_length_aa", 0),
        card.get("mrna_length_nt", 0),
        len(card.get("epitopes", [])),
    )

    # ----- A1: Redundancia de epitopos -----
    logger.info("--- A1: Epitope Redundancy Analysis ---")
    redundancy = A1.run_redundancy_analysis()

    logger.info(
        "Epitopes: %d total, %d unique genes, %d unique alleles",
        redundancy.total_epitopes,
        redundancy.unique_genes,
        redundancy.unique_alleles,
    )
    logger.info(
        "Redundant pairs: %d, minimum for 95%% coverage: %d epitopes",
        len(redundancy.redundant_pairs),
        redundancy.minimum_epitopes_95pct,
    )

    # ----- A2: Pseudouridina -----
    logger.info("--- A2: Pseudouridine Optimization ---")
    psi_u = A2.run_pseudouridine_analysis()

    logger.info(
        "mRNA: %d nt, %d uridines (%.1f%%), GC=%.1f%%",
        psi_u.mrna_length_nt,
        psi_u.total_uridines,
        psi_u.uridine_fraction * 100,
        psi_u.gc_content * 100,
    )
    logger.info("Optimal psiU fraction: %.0f%%", psi_u.optimal_fraction * 100)

    # ----- A3: LNP veterinaria -----
    logger.info("--- A3: Veterinary LNP Formulation ---")
    lnp = A3.run_lnp_formulation(mrna_length_nt=psi_u.mrna_length_nt)

    logger.info(
        "Dose: %.0f ug (%.0f ug/kg x %.0f kg), N/P=%.0f",
        lnp.dose.total_dose_ug,
        lnp.dose.dose_ug_per_kg,
        lnp.dose.body_weight_kg,
        lnp.composition.np_ratio,
    )
    logger.info(
        "Cost/animal: $%.2f (lab), $%.4f (industrial)",
        lnp.cost_lab.total_per_animal_usd,
        lnp.cost_industrial.total_per_animal_usd,
    )

    # ----- A4: Analise estrategica -----
    logger.info("--- A4: Strategic Funding & Market Analysis ---")
    strategy = A4.run_strategic_analysis()

    logger.info(
        "Market: %s dogs at risk, R$ %.0fM/year addressable",
        f"{strategy.market.dogs_at_risk:,}",
        strategy.market.market_value_brl_year / 1e6,
    )
    logger.info(
        "Top funding: %s (score %.1f/10)",
        strategy.funding_ranked[0],
        max(fs["score"] for fs in strategy.funding_sources),
    )

    # ----- Montar resultado final -----
    now = datetime.now(timezone.utc).isoformat()

    results = {
        "metadata": {
            "generated": now,
            "pipeline_version": "0.1.0",
            "platform": "A_mrna_optimization",
            "source_construct": str(CONSTRUCT_CARD_PATH),
            "source_construct_length_aa": card.get("protein_length_aa", 0),
            "source_mrna_length_nt": card.get("mrna_length_nt", 0),
        },
        "epitope_redundancy": A1.to_dict(redundancy),
        "pseudouridine_optimization": A2.to_dict(psi_u),
        "lnp_formulation": A3.to_dict(lnp),
        "strategic_analysis": A4.to_dict(strategy),
    }

    # ----- Salvar resultados -----
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    json_path = OUTPUT_DIR / "platform_a_results.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    logger.info("Results saved: %s", json_path)

    report_path = OUTPUT_DIR / "platform_a_report.md"
    report_md = _generate_markdown_report(results)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Report saved: %s", report_path)

    logger.info("=" * 60)
    logger.info("Platform A pipeline complete.")
    logger.info("=" * 60)

    return results


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("Pipeline failed: %s", exc, exc_info=True)
        sys.exit(1)
