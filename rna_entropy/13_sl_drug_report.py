"""Relatorio comparativo: ASO vs. molecula pequena oral contra SL RNA.

Este modulo gera um relatorio final comparando as duas estrategias
terapeuticas contra o Spliced Leader RNA de L. infantum:

1. MRL-ASO-001 (antisense oligonucleotide) — injetavel, alta especificidade
2. Moleculas pequenas orais — identificadas por docking contra RNA 3D

O relatorio inclui:
- Ranking dos melhores compostos orais
- Comparacao ASO vs. small molecule (via oral, custo, estabilidade)
- Analise de Lipinski e propriedades ADMET estimadas
- Recomendacoes para validacao experimental
- Proximos passos (MD, ensaios in vitro)

Uso:
    python -m rna_entropy.13_sl_drug_report
    python -m rna_entropy.13_sl_drug_report --force
    python -m rna_entropy.13_sl_drug_report --dry-run
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from core.logger import get_logger

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

RESULTS_DIR: Final[Path] = Path("results/sl_drug")
DOCKING_JSON: Final[Path] = RESULTS_DIR / "sl_docking_results.json"
METADATA_JSON: Final[Path] = RESULTS_DIR / "sl_rna_3d_metadata.json"
ASO_JSON: Final[Path] = Path("results/aso/aso_candidates.json")

REPORT_MD: Final[Path] = RESULTS_DIR / "sl_drug_report.md"
COMPARISON_JSON: Final[Path] = RESULTS_DIR / "sl_drug_comparison.json"

logger = get_logger("sl_drug_report")


# ---------------------------------------------------------------------------
# Carregamento de dados
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> dict[str, Any] | list[Any] | None:
    """Carrega arquivo JSON se existir."""
    if not path.exists():
        logger.warning("Arquivo nao encontrado: %s", path)
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Geracao do relatorio
# ---------------------------------------------------------------------------


def _build_header() -> str:
    """Cabecalho do relatorio."""
    return (
        "# Marley - Relatorio de Drug Design contra SL RNA\n"
        "## Leishmania infantum - Leishmaniose Visceral Canina\n\n"
        f"**Data:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n\n"
        "---\n\n"
    )


def _build_target_section() -> str:
    """Secao sobre o alvo terapeutico."""
    return (
        "## 1. Alvo Terapeutico: Spliced Leader RNA\n\n"
        "O **Spliced Leader RNA** (SL RNA) e uma molecula de 39 nucleotideos "
        "que e adicionada ao 5' de **todos os mRNAs** de tripanossomatideos "
        "atraves de um processo chamado **trans-splicing**.\n\n"
        "### Caracteristicas do SL RNA de L. infantum\n\n"
        "| Propriedade | Valor |\n"
        "|------------|-------|\n"
        "| Sequencia | `AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG` |\n"
        "| Comprimento | 39 nt |\n"
        "| Conteudo GC | 28.2% (rico em AT: 72%) |\n"
        "| Presente em humanos | **NAO** (confirmado via BLAST) |\n"
        "| Funcao | Essencial para maturacao de todo mRNA |\n"
        "| Seletividade | Perfeita (exclusivo de tripanossomatideos) |\n\n"
        "### Por que o SL RNA e um alvo ideal?\n\n"
        "1. **Essencialidade**: Sem trans-splicing, nenhum mRNA e funcional. "
        "Bloquear o SL RNA mata o parasita.\n"
        "2. **Seletividade**: Humanos e caes NAO possuem trans-splicing. "
        "Toxicidade ao hospedeiro e minimizada.\n"
        "3. **Conservacao**: O SL RNA e identico em todas as especies de "
        "Leishmania, permitindo droga de amplo espectro.\n\n"
        "---\n\n"
    )


def _build_aso_section(aso_data: dict[str, Any] | None) -> str:
    """Secao sobre o ASO (MRL-ASO-001)."""
    text = (
        "## 2. Estrategia 1: Antisense Oligonucleotide (ASO)\n\n"
        "### MRL-ASO-001\n\n"
    )

    if aso_data and "candidates" in aso_data:
        top = aso_data["candidates"][0] if aso_data["candidates"] else None
        if top:
            text += (
                f"| Propriedade | Valor |\n"
                f"|------------|-------|\n"
                f"| ID | {top.get('aso_id', 'MRL-ASO-001')} |\n"
                f"| Sequencia | `{top.get('aso_sequence', 'N/A')}` |\n"
                f"| Comprimento | {top.get('length', 'N/A')} nt |\n"
                f"| Tm | {top.get('tm_celsius', 'N/A')} C |\n"
                f"| delta_G | {top.get('delta_g_kcal', 'N/A')} kcal/mol |\n"
                f"| GC | {top.get('gc_content', 'N/A')} |\n"
                f"| Score | {top.get('composite_score', 'N/A')} |\n\n"
            )
        text += (
            f"Total de candidatos ASO avaliados: "
            f"{aso_data.get('total_after_filter', 'N/A')}\n\n"
        )
    else:
        text += "*Dados do ASO nao disponiveis. Execute o modulo 08 primeiro.*\n\n"

    text += (
        "### Vantagens do ASO\n"
        "- Alta especificidade (Watson-Crick base pairing)\n"
        "- Mecanismo de acao bem compreendido\n"
        "- Possibilidade de modificacoes quimicas (LNA, PS backbone)\n\n"
        "### Limitacoes do ASO\n"
        "- **Requer injecao** (nao oral)\n"
        "- Custo elevado de producao\n"
        "- Meia-vida limitada in vivo sem modificacoes\n"
        "- Dificuldade de delivery ao parasita intracelular\n\n"
        "---\n\n"
    )

    return text


def _build_docking_section(
    docking_data: dict[str, Any] | None,
) -> str:
    """Secao sobre resultados de docking."""
    text = (
        "## 3. Estrategia 2: Molecula Pequena Oral (Small Molecule)\n\n"
        "### Abordagem\n\n"
        "Docking molecular de uma biblioteca de compostos RNA-targeting "
        "contra a estrutura 3D do SL RNA, buscando moleculas que:\n"
        "1. Ligam ao RNA com afinidade >= -6.0 kcal/mol\n"
        "2. Passam nas regras de Lipinski (biodisponibilidade oral)\n"
        "3. Bloqueiam trans-splicing por interferencia estereoquimica\n\n"
    )

    if docking_data is None:
        text += (
            "*Dados de docking nao disponiveis. "
            "Execute o modulo 12 primeiro.*\n\n"
        )
        return text + "---\n\n"

    # Informacoes gerais
    text += (
        "### Parametros do Docking\n\n"
        f"| Parametro | Valor |\n"
        f"|-----------|-------|\n"
        f"| Modo | {docking_data.get('docking_mode', 'N/A')} |\n"
        f"| Exhaustiveness | {docking_data.get('exhaustiveness', 'N/A')} |\n"
        f"| Biblioteca | {docking_data.get('library_size', 'N/A')} compostos |\n"
        f"| Lipinski pass | {docking_data.get('lipinski_passed', 'N/A')} |\n"
        f"| Pockets | {', '.join(docking_data.get('pockets_used', []))} |\n"
        f"| Total runs | {docking_data.get('total_docking_runs', 'N/A')} |\n\n"
    )

    # Top 10 resultados
    results = docking_data.get("results", [])
    if results:
        text += "### Top 10 Candidatos\n\n"
        text += (
            "| Rank | ID | Nome | Afinidade (kcal/mol) | Pocket | MW | "
            "Lipinski | Mecanismo |\n"
            "|------|-----|------|---------------------|--------|-----|"
            "---------|------------|\n"
        )

        for i, r in enumerate(results[:10], 1):
            lip = "SIM" if r.get("lipinski_pass") else "NAO"
            text += (
                f"| {i} | {r.get('compound_id', '')} | "
                f"{r.get('name', '')} | "
                f"{r.get('binding_affinity', 0):.2f} | "
                f"{r.get('pocket', '')} | "
                f"{r.get('mw', 0):.0f} | "
                f"{lip} | "
                f"{r.get('mechanism', '')[:30]} |\n"
            )
        text += "\n"

        # Melhores candidatos orais
        oral_hits = [
            r for r in results
            if r.get("lipinski_pass") and r.get("binding_affinity", 0) <= -6.0
        ]
        if oral_hits:
            text += (
                f"### Candidatos Orais Promissores "
                f"(Lipinski + afinidade <= -6.0 kcal/mol)\n\n"
                f"**{len(oral_hits)} compostos** passam ambos os criterios.\n\n"
            )
            for r in oral_hits[:5]:
                text += (
                    f"- **{r.get('compound_id', '')} ({r.get('name', '')})**\n"
                    f"  - Afinidade: {r.get('binding_affinity', 0):.2f} kcal/mol\n"
                    f"  - MW: {r.get('mw', 0):.0f} Da\n"
                    f"  - Mecanismo: {r.get('mechanism', '')}\n\n"
                )
        else:
            text += (
                "### Candidatos Orais Promissores\n\n"
                "Nenhum composto atinge ambos os criterios "
                "(Lipinski + afinidade <= -6.0 kcal/mol). "
                "Considerar:\n"
                "- Otimizacao de scaffolds com melhor afinidade\n"
                "- Relaxamento do cutoff para -5.5 kcal/mol\n"
                "- Moleculas maiores com profarmaco oral\n\n"
            )

    text += "---\n\n"
    return text


def _build_comparison_section(
    aso_data: dict[str, Any] | None,
    docking_data: dict[str, Any] | None,
) -> str:
    """Secao comparativa ASO vs. molecula pequena."""
    text = (
        "## 4. Comparacao: ASO vs. Molecula Pequena Oral\n\n"
        "| Criterio | ASO (MRL-ASO-001) | Molecula Pequena Oral |\n"
        "|----------|-------------------|----------------------|\n"
        "| Via de administracao | **Injecao** (SC/IV) | **Oral** |\n"
        "| Especificidade | Muito alta (Watson-Crick) | Moderada (shape-based) |\n"
        "| Custo de producao | Alto ($1000-5000/g) | Baixo ($10-100/g) |\n"
        "| Estabilidade | Requer modificacoes (LNA, PS) | Estavel como comprimido |\n"
        "| Delivery ao parasita | Dificil (intracelular) | Via absoricao GI |\n"
        "| Resistencia | Baixa (complementaridade exata) | Possivel (mutacoes no RNA) |\n"
        "| Tempo para clinica | 5-8 anos | 8-12 anos |\n"
        "| Precedente regulatorio | Antisense aprovados (Nusinersen) | RNA-targeting oral aprovado (Risdiplam) |\n"
        "| Adesao do paciente | Baixa (injecoes) | Alta (comprimido) |\n"
        "| Uso veterinario | Viavel (injecao e comum) | Preferivel (via oral) |\n\n"
    )

    # Recomendacao estrategica
    text += (
        "### Recomendacao Estrategica\n\n"
        "Para **leishmaniose visceral canina**, a estrategia ideal e uma "
        "**abordagem dual**:\n\n"
        "1. **Curto prazo (2-3 anos)**: Desenvolver o **ASO (MRL-ASO-001)** "
        "com modificacoes LNA para uso veterinario. Injecoes sao aceitas "
        "na clinica veterinaria e o mecanismo e mais previsivel.\n\n"
        "2. **Medio prazo (5-7 anos)**: Otimizar os melhores candidatos "
        "orais identificados pelo docking. Validar com ensaios de ligacao "
        "RNA (EMSA, SPR) e ensaios celulares contra promastigotas.\n\n"
        "3. **Longo prazo**: Se os compostos orais mostrarem eficacia, "
        "desenvolver para uso humano tambem (leishmaniose visceral humana "
        "tem as mesmas necessidades terapeuticas).\n\n"
        "---\n\n"
    )

    return text


def _build_validation_section() -> str:
    """Secao de validacao experimental proposta."""
    return (
        "## 5. Proximos Passos: Validacao Experimental\n\n"
        "### 5.1 Validacao In Silico (imediato)\n"
        "- [ ] Dinamica Molecular (MD) dos top 5 complexos RNA-ligando (50 ns)\n"
        "- [ ] Calculo de energia livre (MM-PBSA/MM-GBSA)\n"
        "- [ ] Re-docking com software especifico para RNA (rDock, DOCK6)\n"
        "- [ ] Predicao ADMET completa (pkCSM, SwissADME)\n\n"
        "### 5.2 Validacao In Vitro (3-6 meses)\n"
        "- [ ] Ensaio EMSA (Electrophoretic Mobility Shift Assay) "
        "para confirmar ligacao RNA\n"
        "- [ ] SPR (Surface Plasmon Resonance) para medir Kd\n"
        "- [ ] ITC (Isothermal Titration Calorimetry) para termodinamica\n"
        "- [ ] Ensaio de trans-splicing in vitro com SL RNA + compostos\n\n"
        "### 5.3 Validacao Celular (6-12 meses)\n"
        "- [ ] IC50 contra promastigotas de L. infantum\n"
        "- [ ] CC50 em celulas de mamifero (indice de seletividade)\n"
        "- [ ] Ensaio de amastigotas intracelulares em macrofagos\n"
        "- [ ] qRT-PCR para medir efeito no trans-splicing\n\n"
        "### 5.4 Validacao In Vivo (12-24 meses)\n"
        "- [ ] Farmacocinetica oral em camundongos\n"
        "- [ ] Modelo murino de leishmaniose visceral (BALB/c)\n"
        "- [ ] Modelo canino (se resultados murinos forem positivos)\n\n"
        "---\n\n"
    )


def _build_limitations_section() -> str:
    """Secao de limitacoes."""
    return (
        "## 6. Limitacoes\n\n"
        "### Docking contra RNA\n"
        "- AutoDock Vina foi parametrizado para proteinas. Scoring functions "
        "especificas para RNA (rDock, DOCK6, RLDock) sao mais apropriadas.\n"
        "- A estrutura 3D do SL RNA (39 nt) nao foi resolvida "
        "experimentalmente. Usamos modelo computacional.\n"
        "- RNA e flexivel; docking rigido subestima a adaptacao conformacional.\n\n"
        "### Biblioteca de Compostos\n"
        "- Biblioteca pequena (15 compostos curados). Um screening virtual "
        "completo usaria ~1M compostos (ZINC, Enamine REAL).\n"
        "- SMILES simplificados para alguns compostos (risdiplam, branaplam) "
        "nao representam a molecula completa.\n\n"
        "### Validacao\n"
        "- Nenhum resultado foi validado experimentalmente.\n"
        "- Afinidades de docking contra RNA tem correlacao limitada com Kd real.\n"
        "- Penetracao celular e acesso ao parasita intracelular nao foram "
        "modelados.\n\n"
        "---\n\n"
    )


def _build_references_section() -> str:
    """Secao de referencias."""
    return (
        "## 7. Referencias Tecnicas\n\n"
        "1. Ratmeyer LS et al. (1996) Sequence-specific thermodynamic and "
        "structural properties for DNA-RNA duplexes. *Biochemistry*.\n"
        "2. Warner KD et al. (2018) Principles for targeting RNA with "
        "drug-like small molecules. *Nat Rev Drug Discov* 17:547-558.\n"
        "3. Donlic A, Hargrove AE (2018) Targeting RNA in mammalian "
        "systems with small molecules. *WIREs RNA* 9:e1477.\n"
        "4. Childs-Disney JL et al. (2022) Targeting RNA structures with "
        "small molecules. *Nat Rev Drug Discov* 21:736-762.\n"
        "5. Liang XH et al. (2003) Trans and cis splicing in "
        "trypanosomatids: mechanism, factors, and regulation. "
        "*Eukaryotic Cell* 2:830-840.\n"
        "6. Trott O, Olson AJ (2010) AutoDock Vina: improving the speed "
        "and accuracy of docking. *J Comput Chem* 31:455-461.\n"
        "7. Ratni H et al. (2018) Discovery of Risdiplam, a Selective "
        "Survival of Motor Neuron-2 (SMN2) Gene Splicing Modifier. "
        "*J Med Chem* 61:6501-6517.\n\n"
        "---\n\n"
        "*Relatorio gerado automaticamente pelo pipeline Marley.*\n"
        "*Modulos: 11_sl_rna_3d, 12_sl_rna_docking, 13_sl_drug_report*\n"
    )


# ---------------------------------------------------------------------------
# API publica
# ---------------------------------------------------------------------------


def generate_sl_drug_report(
    force: bool = False,
    dry_run: bool = False,
) -> str:
    """Gera relatorio comparativo ASO vs. molecula pequena oral.

    Carrega resultados dos modulos 08 (ASO), 11 (3D), e 12 (docking)
    e gera um relatorio Markdown + JSON comparativo.

    Args:
        force: Regenerar mesmo se relatorio ja existe.
        dry_run: Modo teste (gera relatorio com dados disponiveis).

    Returns:
        Caminho do relatorio Markdown.
    """
    if REPORT_MD.exists() and not force:
        logger.info(
            "Relatorio ja existe em %s. Use --force para regenerar.",
            REPORT_MD,
        )
        return str(REPORT_MD)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Carregar dados
    aso_data = _load_json(ASO_JSON)
    docking_data = _load_json(DOCKING_JSON)
    metadata = _load_json(METADATA_JSON)

    if dry_run:
        logger.info("[DRY RUN] Gerando relatorio com dados disponiveis.")

    # Construir relatorio
    report_parts: list[str] = [
        _build_header(),
        _build_target_section(),
        _build_aso_section(aso_data),
        _build_docking_section(docking_data),
        _build_comparison_section(aso_data, docking_data),
        _build_validation_section(),
        _build_limitations_section(),
        _build_references_section(),
    ]

    report_text = "".join(report_parts)

    # Salvar Markdown
    with open(REPORT_MD, "w", encoding="utf-8") as fh:
        fh.write(report_text)
    logger.info("Relatorio Markdown salvo em %s", REPORT_MD)

    # Salvar JSON comparativo
    comparison: dict[str, Any] = {
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "strategies": {
            "aso": {
                "id": "MRL-ASO-001",
                "type": "Antisense Oligonucleotide",
                "route": "injection",
                "status": "designed",
                "top_candidate": (
                    aso_data["candidates"][0]
                    if aso_data and "candidates" in aso_data
                    and aso_data["candidates"]
                    else None
                ),
            },
            "small_molecule": {
                "type": "RNA-targeting small molecule",
                "route": "oral",
                "status": "virtual_screening",
                "top_candidates": (
                    docking_data.get("top_5", [])
                    if docking_data
                    else []
                ),
                "library_size": (
                    docking_data.get("library_size", 0)
                    if docking_data
                    else 0
                ),
            },
        },
        "recommendation": (
            "Abordagem dual: ASO para veterinaria (curto prazo) + "
            "molecula oral otimizada (medio prazo). "
            "Ambas as estrategias visam o mesmo alvo (SL RNA) "
            "por mecanismos complementares."
        ),
    }

    with open(COMPARISON_JSON, "w", encoding="utf-8") as fh:
        json.dump(comparison, fh, indent=2, ensure_ascii=False)
    logger.info("Comparacao JSON salva em %s", COMPARISON_JSON)

    # Resumo
    logger.info("=" * 60)
    logger.info("RELATORIO GERADO")
    logger.info("  Markdown: %s", REPORT_MD)
    logger.info("  JSON: %s", COMPARISON_JSON)
    if docking_data and "results" in docking_data:
        oral_hits = [
            r for r in docking_data["results"]
            if r.get("lipinski_pass") and r.get("binding_affinity", 0) <= -6.0
        ]
        logger.info("  Candidatos orais promissores: %d", len(oral_hits))
    logger.info("=" * 60)

    return str(REPORT_MD)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Gerar relatorio comparativo ASO vs. molecula pequena oral "
            "contra o SL RNA de L. infantum."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerar mesmo se relatorio ja existe.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Modo teste.",
    )
    args = parser.parse_args()

    result = generate_sl_drug_report(force=args.force, dry_run=args.dry_run)
    logger.info("Completo: %s", result)
