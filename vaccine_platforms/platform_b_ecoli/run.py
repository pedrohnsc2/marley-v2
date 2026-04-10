"""Plataforma B -- Vacina recombinante expressa em E. coli.

Ponto de entrada principal que orquestra os quatro sub-modulos:
    B1. Redesenho do constructo (construct.py)
    B2. Otimizacao de codons para E. coli K12 (codon_optimizer.py)
    B3. Protocolo de purificacao in silico (purification.py)
    B4. Modelo de custos (cost_model.py)

Le os dados dos epitopos do construct_card.json (plataforma mRNA) e
gera um constructo completo redesenhado para expressao procariotica.

Saidas:
    - results/platform_b_results.json -- dados estruturados completos
    - results/platform_b_report.md -- relatorio legivel em Markdown

Usage:
    python -m vaccine_platforms.platform_b_ecoli.run
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Final

from core.logger import get_logger

from . import construct as B1
from . import codon_optimizer as B2
from . import purification as B3
from . import cost_model as B4

logger = get_logger("platform_b_ecoli")

# ---------------------------------------------------------------------------
# Caminhos
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent.parent
CONSTRUCT_CARD_PATH: Final[Path] = PROJECT_ROOT / "results" / "construct" / "construct_card.json"
OUTPUT_DIR: Final[Path] = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_construct_card() -> dict:
    """Le o construct_card.json da plataforma mRNA.

    Este arquivo contem os epitopos selecionados pelo pipeline de
    vacinologia reversa -- sao compartilhados entre as plataformas.

    Returns:
        Dicionario com os dados do constructo mRNA.

    Raises:
        FileNotFoundError: Se o arquivo nao existe.
    """
    if not CONSTRUCT_CARD_PATH.exists():
        raise FileNotFoundError(
            f"construct_card.json nao encontrado em {CONSTRUCT_CARD_PATH}. "
            f"Execute o pipeline de construcao mRNA (pipeline/06_construct.py) primeiro."
        )

    with open(CONSTRUCT_CARD_PATH, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _generate_markdown_report(results: dict) -> str:
    """Gera o relatorio em Markdown para a plataforma B.

    Args:
        results: Dicionario completo de resultados.

    Returns:
        String com o relatorio formatado em Markdown.
    """
    c = results["construct"]
    co = results["codon_optimization"]
    p = results["purification"]
    cost = results["cost"]

    # Epitopos unicos
    epitopes = c.get("unique_epitopes", [])

    report = f"""# Platform B -- Recombinant E. coli Vaccine

**Generated:** {results['metadata']['generated']}
**Pipeline version:** {results['metadata']['pipeline_version']}

## 1. Construct Redesign (B1)

### Architecture

```
[NdeI/Met] - [His6] - [TEV site] - [L7/L12 adjuvant] - EAAAK - [Epitope cassette] - [Stop] - [XhoI]
```

| Property | Value |
|---|---|
| Total protein length | {c['length_aa']} aa |
| Molecular weight | {c['molecular_weight_da']:,.2f} Da ({c['molecular_weight_da']/1000:.1f} kDa) |
| Isoelectric point (pI) | {c['isoelectric_point']:.2f} |
| Instability index | {c['instability_index']:.2f} ({'stable' if c['is_stable'] else 'UNSTABLE'}) |
| Unique epitopes | {c['unique_epitope_count']} |
| Internal NdeI sites | {'None' if not c['has_internal_ndei'] else 'FOUND'} |
| Internal XhoI sites | {'None' if not c['has_internal_xhoi'] else 'FOUND'} |

### Modifications from mRNA Platform

| Component | mRNA (Platform A) | E. coli (Platform B) |
|---|---|---|
| Signal peptide | tPA (23 aa) | Removed |
| N-terminal tag | None | Met + His6 + TEV site |
| Codon optimization | *Canis lupus* | *E. coli* K12 |
| UTRs | 5' UTR + 3' UTR + poly(A) | None (bacterial) |
| Cloning sites | None | NdeI / XhoI (pET-28a) |
| Stop codon | TGA | TAA (preferred in E. coli) |

### Epitope Cassette

| # | Peptide | Linker |
|---|---|---|
"""
    for i, ep in enumerate(epitopes):
        linker = "AAY" if i < len(epitopes) - 1 else "(C-terminal)"
        report += f"| {i+1} | `{ep}` | {linker} |\n"

    report += f"""
### Protein Sequence

```
{_format_sequence(c['protein_sequence'], width=70)}
```

## 2. Codon Optimization (B2)

| Metric | Value | Target |
|---|---|---|
| DNA insert length | {co['length_nt']} nt | -- |
| GC content | {co['gc_content']:.1%} | 40-60% |
| CAI (Codon Adaptation Index) | {co['cai']:.4f} | > 0.80 |
| Rare codons (< 10%) | {len(co['rare_codons'])} | 0 |
| Internal NdeI sites | {'None' if not co['has_internal_ndei'] else 'FOUND'} | None |
| Internal XhoI sites | {'None' if not co['has_internal_xhoi'] else 'FOUND'} | None |
| Internal Shine-Dalgarno | {'None' if not co['has_shine_dalgarno'] else 'FOUND'} | None |
| Homopolymer runs > 8 nt | {len(co['homopolymer_runs'])} | 0 |
"""

    if co["restriction_sites_found"]:
        report += "\n**Restriction sites found:**\n"
        for site in co["restriction_sites_found"]:
            report += f"- {site['enzyme']} at position {site['position']}\n"

    if co["warnings"]:
        report += "\n**Warnings:**\n"
        for w in co["warnings"]:
            report += f"- {w}\n"

    report += f"""
### DNA Insert Sequence

```
{_format_sequence(co['dna_sequence'], width=70)}
```

## 3. Purification Protocol (B3)

**Expression system:** {p['expression_system']} / {p['expression_vector']}
**Induction:** {p['expression_conditions']['iptg_mm']} mM IPTG at OD600 = {p['expression_conditions']['od600_induction']}, {p['expression_conditions']['induction_temp_c']}C, {p['expression_conditions']['induction_time_h']}h

| Step | Method | Input (mg) | Output (mg) | Efficiency | Purity |
|---|---|---|---|---|---|
"""
    for step in p["steps"]:
        report += (
            f"| {step['step_number']}. {step['name']} | "
            f"{step['method'][:40]}{'...' if len(step['method']) > 40 else ''} | "
            f"{step['input_mg']:.1f} | {step['output_mg']:.1f} | "
            f"{step['efficiency_pct']:.0f}% | {step['purity_after_pct']:.0f}% |\n"
        )

    report += f"""
| Metric | Value |
|---|---|
| Culture volume | {p['culture_volume_l']} L |
| Initial yield (expression) | {p['initial_yield_mg']:.1f} mg |
| Final yield (purified) | {p['final_yield_mg']:.1f} mg |
| Final purity | {p['final_purity_pct']:.1f}% |
| Overall recovery | {p['overall_recovery_pct']:.1f}% |
| Expected yield range | {p['expected_yield_range_mg'][0]:.1f} - {p['expected_yield_range_mg'][1]:.1f} mg/L |

## 4. Cost Model (B4)

### Laboratory Scale (1L culture)

| Item | Cost (USD) |
|---|---|
| Fermentation | ${cost['lab_scale']['fermentation_cost_per_l']:.2f}/L |
| Purification | ${cost['lab_scale']['purification_cost_per_run']:.2f}/run |
| Total per L | ${cost['lab_scale']['total_cost_per_l']:.2f} |
| **Cost per mg** | **${cost['lab_scale']['cost_per_mg_usd']:.2f}** |
| **Cost per dose (50 ug)** | **${cost['per_dose_total']['total_lab_usd']:.4f}** |
| **Cost per animal (3 doses)** | **${cost['lab_scale']['cost_per_animal_usd']:.4f}** |

### Industrial Scale (1000L fermentor, {cost['industrial_scale']['scale_factor']}x factor)

| Item | Cost (USD) |
|---|---|
| **Cost per mg** | **${cost['industrial_scale']['cost_per_mg_usd']:.4f}** |
| **Cost per dose (50 ug)** | **${cost['per_dose_total']['total_industrial_usd']:.4f}** |
| **Cost per animal (3 doses)** | **${cost['industrial_scale']['cost_per_animal_usd']:.4f}** |

### Cold Chain

- **Storage:** {cost['cold_chain']}
- {cost['cold_chain_note']}

### Per-Dose Breakdown

| Component | Lab | Industrial |
|---|---|---|
| Protein | ${cost['per_dose_total']['protein_lab_usd']:.4f} | ${cost['per_dose_total']['protein_industrial_usd']:.4f} |
| QuilA adjuvant | ${cost['per_dose_total']['adjuvant_quila_usd']:.2f} | ${cost['per_dose_total']['adjuvant_quila_usd']:.2f} |
| Formulation | ${cost['per_dose_total']['formulation_usd']:.2f} | ${cost['per_dose_total']['formulation_usd']:.2f} |
| **Total** | **${cost['per_dose_total']['total_lab_usd']:.4f}** | **${cost['per_dose_total']['total_industrial_usd']:.4f}** |

---

*Report generated by the Marley reverse vaccinology pipeline -- Platform B (E. coli recombinant)*
"""
    return report


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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> dict:
    """Executa o pipeline completo da plataforma B (E. coli recombinante).

    Etapas:
        1. Le o construct_card.json com os epitopos da plataforma mRNA
        2. Redesenha o constructo para E. coli (B1)
        3. Otimiza codons para E. coli K12 (B2)
        4. Modela protocolo de purificacao (B3)
        5. Calcula modelo de custos (B4)
        6. Gera relatorio e salva resultados

    Returns:
        Dicionario completo com todos os resultados.
    """
    logger.info("=" * 60)
    logger.info("Platform B -- E. coli Recombinant Vaccine")
    logger.info("=" * 60)

    # ----- Carregar dados de entrada -----
    logger.info("Loading construct_card.json from mRNA platform...")
    card = _load_construct_card()
    epitopes_raw = card["epitopes"]
    logger.info(
        "Found %d epitope entries (%d unique peptides)",
        len(epitopes_raw),
        len(B1.extract_unique_epitopes(epitopes_raw)),
    )

    # ----- B1: Redesenho do constructo -----
    logger.info("--- B1: Construct Redesign ---")
    ecoli_construct = B1.build_ecoli_construct(epitopes_raw)

    warnings = B1.validate_construct(ecoli_construct)
    for w in warnings:
        logger.warning(w)

    logger.info(
        "Construct: %d aa, MW=%.1f Da, pI=%.2f, II=%.2f (%s)",
        ecoli_construct.length_aa,
        ecoli_construct.molecular_weight_da,
        ecoli_construct.isoelectric_point,
        ecoli_construct.instability_index,
        "stable" if ecoli_construct.is_stable else "UNSTABLE",
    )

    # ----- B2: Otimizacao de codons -----
    logger.info("--- B2: Codon Optimization for E. coli K12 ---")
    codon_result = B2.run_codon_optimization(ecoli_construct.protein_sequence)

    # Atualizar flags de sitios de restricao no constructo
    ecoli_construct.has_internal_ndei = codon_result.has_internal_ndei
    ecoli_construct.has_internal_xhoi = codon_result.has_internal_xhoi

    for w in codon_result.warnings:
        logger.warning(w)

    logger.info(
        "DNA: %d nt, GC=%.1f%%, CAI=%.4f, rare=%d",
        codon_result.length_nt,
        codon_result.gc_content * 100,
        codon_result.cai,
        len(codon_result.rare_codons),
    )

    # ----- B3: Protocolo de purificacao -----
    logger.info("--- B3: Purification Protocol Model ---")
    purification = B3.model_purification(
        molecular_weight_da=ecoli_construct.molecular_weight_da,
        culture_volume_l=1.0,
    )

    logger.info(
        "Purification: %.1f mg initial -> %.1f mg final (%.1f%% recovery), %.1f%% purity",
        purification.initial_yield_mg,
        purification.final_yield_mg,
        purification.overall_recovery_pct,
        purification.final_purity_pct,
    )

    # ----- B4: Modelo de custos -----
    logger.info("--- B4: Cost Model ---")
    costs = B4.calculate_costs(
        yield_mg_per_l=purification.final_yield_mg,
    )

    logger.info(
        "Lab cost: $%.2f/mg, $%.4f/dose; Industrial: $%.4f/dose",
        costs.lab_cost_per_mg,
        costs.total_dose_cost_lab,
        costs.total_dose_cost_industrial,
    )

    # ----- Montar resultado final -----
    now = datetime.now(timezone.utc).isoformat()

    results = {
        "metadata": {
            "generated": now,
            "pipeline_version": "0.1.0",
            "platform": "B_ecoli_recombinant",
            "source_construct": str(CONSTRUCT_CARD_PATH),
        },
        "construct": B1.to_dict(ecoli_construct),
        "codon_optimization": B2.to_dict(codon_result),
        "purification": B3.to_dict(purification),
        "cost": B4.to_dict(costs),
    }

    # ----- Salvar resultados -----
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    json_path = OUTPUT_DIR / "platform_b_results.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    logger.info("Results saved: %s", json_path)

    report_path = OUTPUT_DIR / "platform_b_report.md"
    report_md = _generate_markdown_report(results)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)
    logger.info("Report saved: %s", report_path)

    # Salvar FASTA do constructo proteico
    fasta_path = OUTPUT_DIR / "ecoli_construct.fasta"
    with open(fasta_path, "w", encoding="utf-8") as fh:
        fh.write(
            f">Marley_platform_B_ecoli len={ecoli_construct.length_aa} "
            f"MW={ecoli_construct.molecular_weight_da:.0f}Da "
            f"pI={ecoli_construct.isoelectric_point:.2f}\n"
        )
        for i in range(0, len(ecoli_construct.protein_sequence), 70):
            fh.write(ecoli_construct.protein_sequence[i:i + 70] + "\n")
    logger.info("FASTA saved: %s", fasta_path)

    # Salvar inserto DNA
    dna_path = OUTPUT_DIR / "ecoli_insert.fasta"
    with open(dna_path, "w", encoding="utf-8") as fh:
        fh.write(
            f">Marley_platform_B_insert len={codon_result.length_nt}nt "
            f"GC={codon_result.gc_content:.1%} CAI={codon_result.cai:.4f}\n"
        )
        for i in range(0, len(codon_result.dna_sequence), 70):
            fh.write(codon_result.dna_sequence[i:i + 70] + "\n")
    logger.info("DNA insert FASTA saved: %s", dna_path)

    logger.info("=" * 60)
    logger.info("Platform B pipeline complete.")
    logger.info("=" * 60)

    return results


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("Pipeline failed: %s", exc, exc_info=True)
        sys.exit(1)
