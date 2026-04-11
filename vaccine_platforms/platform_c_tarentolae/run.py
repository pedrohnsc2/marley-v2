"""Plataforma C -- Vacina em Leishmania tarentolae (LEXSY).

Ponto de entrada principal que orquestra os sub-modulos:
    C0. Verificacao de ortogonalidade SL RNA (MANDATORIA)
    C1. Justificativa da plataforma L. tarentolae
    C2. Redesenho do construto (construct.py)
    C3. Otimizacao de codons para trypanosomatideos (codon_optimizer.py)
    C4. Modelo de custos (cost_model.py)
    C5. Avaliacao de compatibilidade ASO + vacina

A verificacao de ortogonalidade (C0) DEVE passar antes que qualquer
resultado da plataforma seja considerado valido. Se falhar, o relatorio
documenta o problema de forma proeminente.

Saidas:
    - results/platform_c_results.json -- dados estruturados completos
    - results/platform_c_report.md   -- relatorio legivel em Markdown

Usage:
    python -m vaccine_platforms.platform_c_tarentolae.run
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Final

from core.logger import get_logger

from . import orthogonality as C0
from . import construct as C2
from . import codon_optimizer as C3
from . import cost_model as C4

logger = get_logger("platform_c_tarentolae")

# ---------------------------------------------------------------------------
# Caminhos
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _format_sequence(seq: str, width: int = 70) -> str:
    """Formata uma sequencia em linhas de largura fixa.

    Args:
        seq: Sequencia (proteica ou DNA).
        width: Largura maxima por linha.

    Returns:
        String formatada com quebras de linha.
    """
    lines = []
    for i in range(0, len(seq), width):
        lines.append(seq[i:i + width])
    return "\n".join(lines)


def _generate_markdown_report(results: dict) -> str:
    """Gera o relatorio em Markdown para a plataforma C.

    Inclui todas as secoes: ortogonalidade, justificativa, construto,
    codons, custos, e avaliacao de compatibilidade ASO.

    Args:
        results: Dicionario completo de resultados.

    Returns:
        String com o relatorio formatado em Markdown.
    """
    orth = results["orthogonality"]
    c = results["construct"]
    co = results["codon_optimization"]
    cost = results["cost"]
    compat = results["aso_compatibility"]

    # Status geral baseado na ortogonalidade
    orth_status = "PASS" if orth["verdict"]["passes_threshold"] else "WARNING"
    orth_icon = "[PASS]" if orth_status == "PASS" else "[WARNING]"

    report = f"""# Platform C -- L. tarentolae Live Vaccine (LEXSY)

**Generated:** {results['metadata']['generated']}
**Pipeline version:** {results['metadata']['pipeline_version']}
**Orthogonality check:** {orth_icon} {orth_status}

---

## 0. MANDATORY: SL RNA Orthogonality Check {orth_icon}

**Purpose:** Verify that ASO MRL-ASO-001 (designed to kill L. infantum) does
not also kill L. tarentolae (our expression vector).

### Sequence Comparison

| Property | L. infantum | L. tarentolae |
|---|---|---|
| SL RNA (39 nt) | `{orth['sequences']['sl_rna_l_infantum']}` | `{orth['sequences']['sl_rna_l_tarentolae']}` |
| ASO target region | `{orth['target_region']['infantum']}` | `{orth['target_region']['tarentolae']}` |

### Alignment

| Metric | Value |
|---|---|
| Alignment length | {orth['alignment']['length']} nt |
| Matches (ASO vs L. tarentolae) | {orth['alignment']['matches']} |
| Mismatches | {orth['alignment']['mismatches']} |
| Mismatch positions | {orth['alignment']['mismatch_positions']} |
| Complementarity to L. infantum | {orth['alignment']['complementarity_infantum']:.1%} |
| Complementarity to L. tarentolae | {orth['alignment']['complementarity_tarentolae']:.1%} |

### Thermodynamics (SantaLucia 1998)

| Parameter | L. infantum | L. tarentolae | Difference |
|---|---|---|---|
| DeltaH (kcal/mol) | {orth['thermodynamics']['l_infantum']['dH_kcal_mol']:.2f} | {orth['thermodynamics']['l_tarentolae']['dH_kcal_mol']:.2f} | -- |
| DeltaS (cal/mol/K) | {orth['thermodynamics']['l_infantum']['dS_cal_mol_K']:.2f} | {orth['thermodynamics']['l_tarentolae']['dS_cal_mol_K']:.2f} | -- |
| Tm (C) | {orth['thermodynamics']['l_infantum']['Tm_celsius']:.2f} | {orth['thermodynamics']['l_tarentolae']['Tm_celsius']:.2f} | {orth['thermodynamics']['l_infantum']['Tm_celsius'] - orth['thermodynamics']['l_tarentolae']['Tm_celsius']:.2f} |
| DeltaG at 37C (kcal/mol) | {orth['thermodynamics']['l_infantum']['dG_37C_kcal_mol']:.2f} | {orth['thermodynamics']['l_tarentolae']['dG_37C_kcal_mol']:.2f} | {orth['thermodynamics']['selectivity']['delta_dG_kcal_mol']:.2f} |

**Selectivity:** DeltaDeltaG = {orth['thermodynamics']['selectivity']['delta_dG_kcal_mol']:.2f} kcal/mol (threshold: >= {orth['thermodynamics']['selectivity']['threshold_kcal_mol']:.1f} kcal/mol)
**Kd ratio:** {orth['thermodynamics']['selectivity']['Kd_ratio_tarentolae_over_infantum']:.1f}x (L. tarentolae Kd / L. infantum Kd)

### Verdict

> {orth['verdict']['summary']}

"""

    # Adicionar avisos
    if orth["verdict"]["warnings"]:
        report += "**Warnings:**\n"
        for w in orth["verdict"]["warnings"]:
            report += f"- {w}\n"
        report += "\n"

    # --- C1: Justificativa ---
    rationale = c.get("rationale", {})
    report += f"""---

## 1. Why L. tarentolae? (C1)

| Aspect | Detail |
|---|---|
| Biosafety | {rationale.get('biosafety', 'N/A')} |
| Growth conditions | {rationale.get('growth_conditions', 'N/A')} |
| PTM capability | {rationale.get('ptm_capability', 'N/A')} |
| Intrinsic adjuvant | {rationale.get('intrinsic_adjuvant', 'N/A')} |
| Literature precedent | {rationale.get('literature_precedent', 'N/A')} |
| Commercial system | {rationale.get('commercial_system', 'N/A')} |
| Vector | {rationale.get('vector', 'N/A')} |
| Selection | {rationale.get('selection', 'N/A')} |

---

## 2. Construct Redesign (C2)

### Architecture

```
[SAP1 signal] - [L7/L12 adjuvant] - EAAAK - [11 epitopes + AAY] - GSGSGS - [His6] - [Strep-tag II]
```

Genomic integration context:
```
[SSU 5' flank] - [5'UTR] - ATG[CDS]TGA - [3'UTR] - [SSU 3' flank]
```

| Property | Value |
|---|---|
| Total protein length | {c['length_aa']} aa |
| Molecular weight | {c['molecular_weight_da']:,.2f} Da ({c['molecular_weight_da']/1000:.1f} kDa) |
| Isoelectric point (pI) | {c['isoelectric_point']:.2f} |
| Instability index | {c['instability_index']:.2f} ({'stable' if c['is_stable'] else 'UNSTABLE'}) |
| GRAVY | {c['gravy']:.4f} |
| Aromaticity | {c['aromaticity']:.4f} |
| Unique epitopes | {c['unique_epitope_count']} |
| Vector | {c['vector']} |
| Selection | {c['selection']} |

### Modifications from mRNA Platform

| Component | mRNA (Platform A) | L. tarentolae (Platform C) |
|---|---|---|
"""
    for mod in c.get("modifications_from_mrna", []):
        report += f"| {mod['component']} | {mod['mrna_platform'][:50]} | {mod['tarentolae_platform'][:50]} |\n"

    # Epitopos
    epitopes = c.get("unique_epitopes", [])
    report += """
### Epitope Cassette

| # | Peptide | Linker |
|---|---|---|
"""
    for i, ep in enumerate(epitopes):
        linker = "AAY" if i < len(epitopes) - 1 else "(C-terminal -> tags)"
        report += f"| {i+1} | `{ep}` | {linker} |\n"

    report += f"""
### Protein Sequence

```
{_format_sequence(c['protein_sequence'], width=70)}
```

---

## 3. Codon Optimization for L. tarentolae (C3)

| Metric | Value | Target |
|---|---|---|
| DNA insert length | {co['length_nt']} nt | -- |
| GC content | {co['gc_content']:.1%} | {co['gc_target']} |
| GC in range | {'Yes' if co['gc_in_range'] else 'NO'} | Yes |
| CAI | {co['cai']:.4f} | {co['cai_target']} |
| CAI above threshold | {'Yes' if co['cai_above_threshold'] else 'NO'} | Yes |
| Rare codons (< 10%) | {co['rare_codon_count']} | 0 |
| Restriction sites | {len(co['restriction_sites_found'])} | 0 |
| Homopolymer runs | {len(co['homopolymer_runs'])} | 0 |
| Poly-A signal | {'Found' if co['has_polya_signal'] else 'None'} | None |
| Poly-T signal | {'Found' if co['has_polyt_signal'] else 'None'} | None |
"""

    if co["warnings"]:
        report += "\n**Warnings:**\n"
        for w in co["warnings"]:
            report += f"- {w}\n"

    report += f"""
### DNA Insert Sequence (CDS only)

```
{_format_sequence(co['dna_sequence'], width=70)}
```

---

## 4. Cost Model (C4)

### Modalidade A: Vacina Viva (Organismo Inteiro)

| Item | Value |
|---|---|
| Growth medium cost | ${cost['live_vaccine']['growth_cost_per_l']:.2f}/L |
| Doses per liter | {cost['live_vaccine']['doses_per_l']:,.0f} |
| **Cost per dose** | **${cost['live_vaccine']['cost_per_dose_usd']:.6f}** |
| **Cost per animal (3 doses)** | **${cost['live_vaccine']['cost_per_animal_usd']:.6f}** |
| Industrial cost per dose | ${cost['live_vaccine']['industrial_cost_per_dose_usd']:.6f} |
| Adjuvant | {cost['live_vaccine']['adjuvant']} |
| Cold chain | {cost['live_vaccine']['cold_chain']} |

### Modalidade B: Proteina Secretada Purificada

| Item | Value |
|---|---|
| Growth + purification cost | ${cost['protein_vaccine']['growth_cost_per_l'] + cost['protein_vaccine']['purification_cost_per_run']:.2f}/L |
| Yield (purified) | {cost['protein_vaccine']['yield_purified_mg_per_l']:.1f} mg/L |
| Cost per mg | ${cost['protein_vaccine']['cost_per_mg_usd']:.2f} |
| Protein cost per dose | ${cost['protein_vaccine']['cost_per_dose_protein_usd']:.4f} |
| **Total cost per dose** | **${cost['protein_vaccine']['cost_per_dose_total_usd']:.4f}** |
| **Cost per animal (3 doses)** | **${cost['protein_vaccine']['cost_per_animal_usd']:.4f}** |
| Industrial total per dose | ${cost['protein_vaccine']['industrial_total_dose_usd']:.4f} |
| Adjuvant | {cost['protein_vaccine']['adjuvant']} |
| Cold chain | {cost['protein_vaccine']['cold_chain']} |

### Setup Costs (One-Time)

| Item | Cost |
|---|---|
| pLEXSY kit | ${cost['setup']['plexsy_kit_cost']:,.2f} |
| Amortized per batch | ${cost['setup']['plexsy_amortized_per_batch']:.2f} |
| BSL level | {cost['setup']['bsl_level']} |
| Growth temperature | {cost['setup']['growth_temperature']} |

### Platform Comparison

| Platform | Cost/dose | Cold chain | Adjuvant |
|---|---|---|---|
"""
    for pname, pdata in cost.get("platform_comparison", {}).items():
        dose_cost = pdata.get("cost_per_dose_usd", pdata.get("estimated_per_dose_usd", "N/A"))
        if isinstance(dose_cost, float):
            dose_cost = f"${dose_cost:.4f}"
        elif isinstance(dose_cost, str) and not dose_cost.startswith("$"):
            dose_cost = f"${dose_cost}"
        report += f"| {pname} | {dose_cost} | {pdata.get('cold_chain', 'N/A')} | {pdata.get('adjuvant', 'N/A')} |\n"

    # --- C5: Compatibilidade ASO ---
    report += f"""
---

## 5. ASO Compatibility Assessment (C5)

### Can MRL-ASO-001 and L. tarentolae vaccine be used together?

**Answer:** {compat['answer']}

### Recommended Protocol

{compat['protocol']}

### Risk Assessment

| Factor | Assessment |
|---|---|
| Thermodynamic selectivity | {compat['risk_assessment']['thermodynamic_selectivity']} |
| Mismatch protection | {compat['risk_assessment']['mismatch_protection']} |
| Temporal separation | {compat['risk_assessment']['temporal_separation']} |
| Overall risk | {compat['risk_assessment']['overall_risk']} |

### Clinical Implications

> {compat['clinical_implications']}

---

*Report generated by the Marley reverse vaccinology pipeline -- Platform C (L. tarentolae LEXSY)*
"""
    return report


# ---------------------------------------------------------------------------
# Avaliacao de compatibilidade ASO
# ---------------------------------------------------------------------------


def _assess_aso_compatibility(orth_result: C0.OrthogonalityResult) -> dict:
    """Avalia a compatibilidade entre ASO MRL-ASO-001 e vacina L. tarentolae.

    Baseado nos resultados da verificacao de ortogonalidade, determina:
    1. Se terapia combinada e viavel
    2. Protocolo recomendado
    3. Riscos e mitigacoes

    Args:
        orth_result: Resultado da verificacao de ortogonalidade.

    Returns:
        Dicionario com avaliacao completa.
    """
    if orth_result.passes_threshold:
        answer = (
            "SIM -- a terapia combinada e termodinamicamente viavel. "
            f"O ASO tem {orth_result.kd_ratio:.1f}x maior afinidade pelo "
            f"SL RNA de L. infantum vs L. tarentolae "
            f"(DeltaDeltaG = {orth_result.delta_dg:.2f} kcal/mol)."
        )
        protocol = (
            "**Protocolo sequencial recomendado:**\n\n"
            "1. **Fase 1 -- Tratamento ASO (dias 0-14):** Administrar MRL-ASO-001 "
            "para reduzir carga parasitaria de L. infantum. A seletividade "
            "termodinamica garante impacto minimo em L. tarentolae caso "
            "haja exposicao acidental.\n\n"
            "2. **Fase 2 -- Washout (dias 15-17):** Aguardar ~5 meias-vidas do "
            "ASO fosforotioado (~48h total) para clearance sistemico.\n\n"
            "3. **Fase 3 -- Vacinacao (dia 18+):** Administrar vacina "
            "L. tarentolae (viva ou proteina secretada). O ambiente "
            "pos-ASO, com carga parasitaria reduzida, pode potencializar "
            "a resposta vacinal.\n\n"
            "4. **Alternativa -- Administracao simultanea:** A diferenca de "
            f"afinidade ({orth_result.kd_ratio:.1f}x) pode permitir uso "
            "simultaneo, mas requer validacao experimental."
        )
        thermo = f"DeltaDeltaG = {orth_result.delta_dg:.2f} kcal/mol (acima do limiar de 2.0)"
        mismatch = (
            f"{orth_result.mismatches} mismatch(es) na regiao alvo reduz(em) "
            f"afinidade do ASO por L. tarentolae"
        )
        temporal = "Washout de 48h e suficiente para clearance do ASO"
        overall = "BAIXO -- terapia combinada e viavel com protocolo sequencial"
    else:
        answer = (
            "COM RESTRICOES -- a diferenca de afinidade e pequena "
            f"(DeltaDeltaG = {orth_result.delta_dg:.2f} kcal/mol, "
            f"abaixo do limiar de {orth_result.threshold_kcal} kcal/mol). "
            "O ASO pode afetar a viabilidade de L. tarentolae."
        )
        protocol = (
            "**Protocolo com precaucoes:**\n\n"
            "1. **NUNCA administrar simultaneamente.** O ASO pode matar "
            "L. tarentolae.\n\n"
            "2. **Opcao A -- Vacinar primeiro:** Administrar vacina L. tarentolae, "
            "aguardar resposta imune (2-4 semanas), depois iniciar ASO.\n\n"
            "3. **Opcao B -- Intervalo longo:** ASO primeiro, aguardar >1 semana "
            "de washout antes da vacinacao.\n\n"
            "4. **Opcao C -- Usar Platform B (E. coli):** Se a compatibilidade "
            "com ASO e prioritaria, a vacina recombinante em E. coli nao "
            "depende de SL RNA e e totalmente compativel.\n\n"
            "5. **Opcao D -- Modificar o ASO:** Introduzir modificacoes quimicas "
            "(LNA, 2'-OMe) nas posicoes de mismatch para aumentar seletividade "
            "por L. infantum."
        )
        thermo = (
            f"DeltaDeltaG = {orth_result.delta_dg:.2f} kcal/mol "
            f"(ABAIXO do limiar de {orth_result.threshold_kcal})"
        )
        mismatch = (
            f"Apenas {orth_result.mismatches} mismatch(es) -- protecao insuficiente"
        )
        temporal = "Requer intervalo prolongado (>1 semana) entre ASO e vacina"
        overall = "MODERADO a ALTO -- requer precaucoes ou plataforma alternativa"

    return {
        "answer": answer,
        "protocol": protocol,
        "risk_assessment": {
            "thermodynamic_selectivity": thermo,
            "mismatch_protection": mismatch,
            "temporal_separation": temporal,
            "overall_risk": overall,
        },
        "clinical_implications": orth_result.clinical_implications,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> dict:
    """Executa o pipeline completo da plataforma C (L. tarentolae).

    Etapas:
        C0. Verificacao de ortogonalidade do SL RNA (MANDATORIA)
        C1. Justificativa (incorporada no construto)
        C2. Redesenho do construto para L. tarentolae
        C3. Otimizacao de codons para trypanosomatideos
        C4. Modelo de custos
        C5. Avaliacao de compatibilidade ASO

    Returns:
        Dicionario completo com todos os resultados.
    """
    logger.info("=" * 60)
    logger.info("Platform C -- L. tarentolae Live Vaccine (LEXSY)")
    logger.info("=" * 60)

    # ===== C0: VERIFICACAO DE ORTOGONALIDADE (MANDATORIA) =====
    logger.info("--- C0: SL RNA Orthogonality Check (MANDATORY) ---")
    logger.info(
        "Comparing ASO MRL-ASO-001 binding to L. infantum vs L. tarentolae SL RNA..."
    )
    orth_result = C0.verify_sl_rna_orthogonality()

    if orth_result.passes_threshold:
        logger.info(
            "PASS: DeltaDeltaG = %.2f kcal/mol (>= %.1f threshold). "
            "ASO is selective for L. infantum.",
            orth_result.delta_dg,
            orth_result.threshold_kcal,
        )
    else:
        logger.warning(
            "WARNING: DeltaDeltaG = %.2f kcal/mol (< %.1f threshold). "
            "ASO may affect L. tarentolae viability!",
            orth_result.delta_dg,
            orth_result.threshold_kcal,
        )

    logger.info(
        "  L. infantum:   DG=%.2f kcal/mol, Tm=%.2f C",
        orth_result.dg_infantum,
        orth_result.tm_infantum,
    )
    logger.info(
        "  L. tarentolae: DG=%.2f kcal/mol, Tm=%.2f C",
        orth_result.dg_tarentolae,
        orth_result.tm_tarentolae,
    )
    logger.info(
        "  Kd ratio: %.1fx selectivity for L. infantum",
        orth_result.kd_ratio,
    )

    for w in orth_result.warnings:
        logger.warning("  %s", w)

    # ===== C2: REDESENHO DO CONSTRUTO =====
    logger.info("--- C2: Construct Redesign for L. tarentolae ---")
    construct = C2.build_tarentolae_construct()

    construct_warnings = C2.validate_construct(construct)
    for w in construct_warnings:
        logger.warning(w)

    logger.info(
        "Construct: %d aa, MW=%.1f Da, pI=%.2f, II=%.2f (%s)",
        construct.length_aa,
        construct.molecular_weight_da,
        construct.isoelectric_point,
        construct.instability_index,
        "stable" if construct.is_stable else "UNSTABLE",
    )

    # ===== C3: OTIMIZACAO DE CODONS =====
    logger.info("--- C3: Codon Optimization for L. tarentolae ---")
    codon_result = C3.run_codon_optimization(construct.protein_sequence)

    for w in codon_result.warnings:
        logger.warning(w)

    logger.info(
        "DNA: %d nt, GC=%.1f%% (%s), CAI=%.4f (%s), rare=%d",
        codon_result.length_nt,
        codon_result.gc_content * 100,
        "OK" if codon_result.gc_in_range else "OUT OF RANGE",
        codon_result.cai,
        "OK" if codon_result.cai_above_threshold else "LOW",
        len(codon_result.rare_codons),
    )

    # ===== C4: MODELO DE CUSTOS =====
    logger.info("--- C4: Cost Model ---")
    costs = C4.calculate_costs()

    logger.info(
        "Live vaccine: $%.6f/dose (%d doses/L)",
        costs.live_cost_per_dose,
        costs.live_doses_per_l,
    )
    logger.info(
        "Protein vaccine: $%.4f/dose (%.1f mg/L purified)",
        costs.protein_total_dose_cost,
        costs.protein_yield_mg_per_l,
    )

    # ===== C5: COMPATIBILIDADE ASO =====
    logger.info("--- C5: ASO Compatibility Assessment ---")
    compat = _assess_aso_compatibility(orth_result)
    logger.info("ASO compatibility: %s", compat["answer"][:100])

    # ===== MONTAR RESULTADO FINAL =====
    now = datetime.now(timezone.utc).isoformat()

    results = {
        "metadata": {
            "generated": now,
            "pipeline_version": "0.1.0",
            "platform": "C_tarentolae_lexsy",
            "orthogonality_status": "PASS" if orth_result.passes_threshold else "WARNING",
        },
        "orthogonality": C0.to_dict(orth_result),
        "construct": C2.to_dict(construct),
        "codon_optimization": C3.to_dict(codon_result),
        "cost": C4.to_dict(costs),
        "aso_compatibility": compat,
    }

    # ===== SALVAR RESULTADOS =====
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # JSON principal
    json_path = OUTPUT_DIR / "platform_c_results.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    logger.info("Results saved: %s", json_path)

    # Relatorio Markdown
    report_path = OUTPUT_DIR / "platform_c_report.md"
    report_md = _generate_markdown_report(results)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Report saved: %s", report_path)

    # FASTA do construto proteico
    fasta_path = OUTPUT_DIR / "tarentolae_construct.fasta"
    with open(fasta_path, "w", encoding="utf-8") as fh:
        fh.write(
            f">Marley_platform_C_tarentolae len={construct.length_aa} "
            f"MW={construct.molecular_weight_da:.0f}Da "
            f"pI={construct.isoelectric_point:.2f} "
            f"vector={construct.vector}\n"
        )
        for i in range(0, len(construct.protein_sequence), 70):
            fh.write(construct.protein_sequence[i:i + 70] + "\n")
    logger.info("FASTA saved: %s", fasta_path)

    # FASTA do inserto DNA (CDS otimizado)
    dna_path = OUTPUT_DIR / "tarentolae_insert.fasta"
    with open(dna_path, "w", encoding="utf-8") as fh:
        fh.write(
            f">Marley_platform_C_insert len={codon_result.length_nt}nt "
            f"GC={codon_result.gc_content:.1%} CAI={codon_result.cai:.4f} "
            f"organism=L.tarentolae\n"
        )
        for i in range(0, len(codon_result.dna_sequence), 70):
            fh.write(codon_result.dna_sequence[i:i + 70] + "\n")
    logger.info("DNA insert FASTA saved: %s", dna_path)

    # Inserto completo com flancos
    full_insert = C3.add_flanking_elements(codon_result.dna_sequence)
    full_path = OUTPUT_DIR / "tarentolae_full_insert.fasta"
    with open(full_path, "w", encoding="utf-8") as fh:
        fh.write(
            f">Marley_platform_C_full_insert len={len(full_insert)}nt "
            f"includes=SSU_flanks+UTRs+CDS+stop\n"
        )
        for i in range(0, len(full_insert), 70):
            fh.write(full_insert[i:i + 70] + "\n")
    logger.info("Full insert FASTA saved: %s", full_path)

    logger.info("=" * 60)
    logger.info(
        "Platform C pipeline complete. Orthogonality: %s",
        "PASS" if orth_result.passes_threshold else "WARNING",
    )
    logger.info("=" * 60)

    return results


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("Pipeline failed: %s", exc, exc_info=True)
        sys.exit(1)
