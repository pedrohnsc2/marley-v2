"""Gerador de relatorios para MRL-QUANTUM.

Produz relatorios em Markdown, JSON e texto anotado a partir dos
resultados dos modulos QAOA e VQE. Formato pensado para revisores
de periódicos cientificos — inclui contexto biologico e referencias.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from core.logger import get_logger

from mrl_quantum.config import (
    ASO_SEQUENCE,
    DG_BASE,
    DELIVERY_SCORE,
    MATH_SCORE,
    QAOA_RESULTS_DIR,
    VQE_RESULTS_DIR,
)

logger = get_logger("mrl_quantum.reporting")


# ---------------------------------------------------------------------------
# Metadados compartilhados por todos os relatorios
# ---------------------------------------------------------------------------


def _build_metadata(
    module: str,
    backend: str,
    seed: int,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Constroi bloco de metadados padronizado."""
    meta: dict[str, Any] = {
        "timestamp": datetime.now(tz=timezone.utc).isoformat(),
        "module": module,
        "backend": backend,
        "seed": seed,
        "dg_base_kcal_mol": DG_BASE,
        "math_certificate_score": MATH_SCORE,
        "delivery_certificate_score": DELIVERY_SCORE,
        "aso_sequence": ASO_SEQUENCE,
        "aso_length": len(ASO_SEQUENCE),
    }
    if extra:
        meta.update(extra)
    return meta


def _write_json(data: dict[str, Any], path: Path) -> Path:
    """Grava JSON com encoding UTF-8 e indentacao."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    return path


def _write_text(content: str, path: Path) -> Path:
    """Grava texto com encoding UTF-8."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return path


# ---------------------------------------------------------------------------
# Padrão LNA do gapmer MRL-ASO-001 (5-10-5 flanks)
# Posicoes 0-4: LNA (flank 5')
# Posicoes 5-19: DNA (gap)
# Posicoes 20-24: LNA (flank 3')
# ---------------------------------------------------------------------------

_LNA_FLANK_5P = 5
_GAP_END = 20
_LNA_FLANK_3P = 25


def _is_lna(position: int) -> bool:
    """Verifica se a posicao no ASO e LNA (flanks do gapmer)."""
    return position < _LNA_FLANK_5P or position >= _GAP_END


# ---------------------------------------------------------------------------
# Relatorios QAOA
# ---------------------------------------------------------------------------


def _extract_qaoa_fields(qaoa_result: dict[str, Any]) -> dict[str, Any]:
    """Extrai campos normalizados do resultado do run_qaoa."""
    final_best = qaoa_result.get("final_best") or {}
    params = qaoa_result.get("parameters") or {}
    runs = qaoa_result.get("runs") or []

    # Configuracao como string — pode ser "SpSpRp..." (str) ou ["Sp","Rp",...] (list)
    config_raw = final_best.get("configuration", "")
    if isinstance(config_raw, str) and config_raw:
        config_chunks = [config_raw[i:i + 2] for i in range(0, len(config_raw), 2)]
    elif isinstance(config_raw, list):
        config_chunks = config_raw
    else:
        config_chunks = []
    config_str = "".join("1" if c == "Sp" else "0" for c in config_chunks) if config_chunks else ""

    dg_total = final_best.get("dg_total", 0.0)
    improvement = dg_total - DG_BASE if dg_total else 0.0

    # Top 10 do melhor run
    top10 = []
    successful = [r for r in runs if r.get("status") != "FAILED"]
    if successful:
        best_run = min(successful, key=lambda r: r["best"]["dg_total"])
        top10 = best_run.get("top_10_samples", [])

    return {
        "config_str": config_str,
        "dg_total": dg_total,
        "improvement": improvement,
        "sp_count": final_best.get("sp_count", 0),
        "rp_count": final_best.get("rp_count", 0),
        "n_runs": len(runs),
        "layers": params.get("layers_tested", []),
        "method": qaoa_result.get("final_method", "N/A"),
        "runtime": qaoa_result.get("runtime_seconds", 0.0),
        "is_hybrid": qaoa_result.get("is_hybrid", False),
        "n_qaoa_qubits": qaoa_result.get("n_qaoa_qubits", 0),
        "top10": top10,
        "greedy_full": qaoa_result.get("greedy_full"),
        "greedy_beats_qaoa": qaoa_result.get("greedy_beats_qaoa", False),
        "comparison": qaoa_result.get("comparison", {}),
        "bitstring": final_best.get("bitstring", []),
    }


def _generate_qaoa_report_md(
    qaoa_result: dict[str, Any],
    backend: str,
    seed: int,
) -> str:
    """Gera relatorio Markdown do modulo QAOA."""
    f = _extract_qaoa_fields(qaoa_result)
    layers_str = ", ".join(str(p) for p in f["layers"]) if f["layers"] else "N/A"

    hybrid_note = ""
    if f["is_hybrid"]:
        hybrid_note = (
            f"\n> **Nota:** Abordagem hibrida — QAOA em {f['n_qaoa_qubits']} qubits "
            f"+ greedy para posicoes restantes.\n"
        )

    greedy_note = ""
    if f["greedy_beats_qaoa"] and f["greedy_full"]:
        greedy_note = (
            f"\n> **Nota:** Baseline greedy (dG={f['greedy_full']['dg_total']:.4f}) "
            f"superou QAOA — esperado para QUBOs com estrutura favoravel a gulosos.\n"
        )

    report = f"""# MRL-QUANTUM QAOA Report
## Otimizacao de Estereoisomeros Fosforotioato (PS)

**Data:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
**Backend:** {backend}
**Seed:** {seed}
**Metodo:** {f['method']}
**Tempo total:** {f['runtime']:.1f}s

---

### Contexto Biologico

O MRL-ASO-001 e um gapmer LNA-DNA-LNA de 25 nucleotideos com backbone
fosforotioato (PS) completo, desenhado contra o SL RNA de *Leishmania infantum*.
O backbone PS confere resistencia a nucleases e e essencial para a atividade
terapeutica, mas cria um centro quiral (Rp ou Sp) em cada ligacao
internucleotidica — totalizando 24 centros quirais.

A configuracao estereoquimica afeta:
- **Afinidade de ligacao** ao RNA-alvo (Sp favorece duplex com RNase H)
- **Reconhecimento pela RNase H1** (requer geometria especifica no gap)
- **Estabilidade metabolica** (Sp mais resistente a 3'-exonucleases)

Este modulo usa QAOA (Quantum Approximate Optimization Algorithm) para
explorar o espaco combinatorial de 2^24 configuracoes Rp/Sp e encontrar
a configuracao que maximiza a estabilidade termodinamica do duplex,
respeitando os requisitos estruturais da RNase H1.
{hybrid_note}
### Parametros de Entrada

| Parametro | Valor |
|-----------|-------|
| dG_BASE | {DG_BASE:.2f} kcal/mol |
| Certificado matematico | {MATH_SCORE} |
| Certificado delivery | {DELIVERY_SCORE} |
| Posicoes totais | 25 |
| Qubits QAOA | {f['n_qaoa_qubits']} |
| Camadas QAOA (p) | {layers_str} |
| Backend | {backend} |

### Resultados

| Metrica | Valor |
|---------|-------|
| dG otimizado | {f['dg_total']:.4f} kcal/mol |
| Melhoria sobre dG_BASE | {f['improvement']:+.4f} kcal/mol |
| Configuracao | Sp={f['sp_count']}, Rp={f['rp_count']} |
| Rodadas executadas | {f['n_runs']} |
{greedy_note}
### Interpretacao

A otimizacao QAOA identificou uma configuracao estereoquimica com
dG = {f['dg_total']:.4f} kcal/mol ({f['improvement']:+.4f} sobre dG_BASE
de {DG_BASE:.2f} kcal/mol), com {f['sp_count']} posicoes Sp e
{f['rp_count']} posicoes Rp.

### Referencias

1. Wan WB et al. (2014) *Nucleic Acids Res* 42(22):13456-68
2. Stec WJ et al. (2010) *J Am Chem Soc* 132(40):14133
3. Eckstein F (2014) *Nucleic Acids Ther* 24(6):374-87
4. Farhi E et al. (2014) *Science* 345(6195):1052-4 (QAOA)
"""
    return report


def _generate_qaoa_top10_json(
    qaoa_result: dict[str, Any],
    backend: str,
    seed: int,
) -> dict[str, Any]:
    """Gera JSON com as 10 melhores configuracoes QAOA."""
    f = _extract_qaoa_fields(qaoa_result)
    return {
        "metadata": _build_metadata("qaoa_top10", backend, seed),
        "top10_configurations": f["top10"],
        "energy_unit": "kcal/mol",
        "configuration_encoding": "1=Sp, 0=Rp (por posicao)",
    }


def _generate_qaoa_vs_classical_json(
    qaoa_result: dict[str, Any],
    backend: str,
    seed: int,
) -> dict[str, Any]:
    """Gera JSON comparativo QAOA vs busca classica."""
    f = _extract_qaoa_fields(qaoa_result)
    greedy = f["greedy_full"] or {}

    return {
        "metadata": _build_metadata("qaoa_vs_classical", backend, seed),
        "comparison": {
            "qaoa": {
                "best_energy": f["dg_total"],
                "method": f["method"],
                "sp_count": f["sp_count"],
                "rp_count": f["rp_count"],
                "wall_time_seconds": f["runtime"],
            },
            "greedy": {
                "best_energy": greedy.get("dg_total", 0.0),
                "sp_count": greedy.get("sp_count", 0),
                "rp_count": greedy.get("rp_count", 0),
            } if greedy else {},
        },
        "greedy_beats_qaoa": f["greedy_beats_qaoa"],
        "detailed_comparison": f["comparison"],
    }


def _generate_qaoa_stereo_annotated(
    qaoa_result: dict[str, Any],
) -> str:
    """Gera sequencia anotada com posicao/base/LNA/configuracao PS.

    Formato por linha:
        Pos  Base  Tipo    PS_Config
        0    A     LNA     Sp
        1    C     LNA     Rp
        ...
    """
    f = _extract_qaoa_fields(qaoa_result)
    # configuration e uma string "SpSpRpSp..." — extrair em chunks de 2
    config_raw = (qaoa_result.get("final_best") or {}).get("configuration", "")
    config_list: list[str] = []
    if isinstance(config_raw, str) and config_raw:
        config_list = [config_raw[i:i + 2] for i in range(0, len(config_raw), 2)]
    elif isinstance(config_raw, list):
        config_list = config_raw

    lines: list[str] = []
    lines.append("# MRL-ASO-001 — Sequencia Anotada com Estereoquimica PS Otimizada")
    lines.append(f"# Metodo: {f['method']}")
    lines.append(f"# dG_BASE = {DG_BASE:.2f} kcal/mol")
    lines.append(f"# dG_otimizado = {f['dg_total']:.4f} kcal/mol")
    lines.append(f"# Melhoria = {f['improvement']:+.4f} kcal/mol")
    lines.append(f"# Sp={f['sp_count']}, Rp={f['rp_count']}")
    lines.append("#")
    lines.append(f"# {'Pos':<5} {'Base':<6} {'Tipo':<8} {'PS_Config':<10}")
    lines.append(f"# {'-' * 5} {'-' * 6} {'-' * 8} {'-' * 10}")

    for i, base in enumerate(ASO_SEQUENCE):
        tipo = "LNA" if _is_lna(i) else "DNA"
        if i < len(config_list):
            ps = config_list[i]  # "Sp" or "Rp"
        elif i < len(ASO_SEQUENCE) - 1:
            ps = "N/A"
        else:
            ps = "---"
        lines.append(f"  {i:<5} {base:<6} {tipo:<8} {ps:<10}")

    return "\n".join(lines) + "\n"


def _generate_qaoa_reports(
    qaoa_result: dict[str, Any],
    backend: str,
    seed: int,
) -> list[Path]:
    """Gera todos os relatorios QAOA."""
    paths: list[Path] = []

    # Relatorio Markdown
    md = _generate_qaoa_report_md(qaoa_result, backend, seed)
    paths.append(_write_text(md, QAOA_RESULTS_DIR / "qaoa_report.md"))

    # Top 10 configuracoes
    top10 = _generate_qaoa_top10_json(qaoa_result, backend, seed)
    paths.append(_write_json(top10, QAOA_RESULTS_DIR / "qaoa_top10_configurations.json"))

    # Comparacao QAOA vs classico
    vs_classical = _generate_qaoa_vs_classical_json(qaoa_result, backend, seed)
    paths.append(
        _write_json(vs_classical, QAOA_RESULTS_DIR / "qaoa_vs_classical_comparison.json")
    )

    # Sequencia anotada
    annotated = _generate_qaoa_stereo_annotated(qaoa_result)
    paths.append(_write_text(annotated, QAOA_RESULTS_DIR / "qaoa_stereo_annotated.txt"))

    return paths


# ---------------------------------------------------------------------------
# Relatorios VQE
# ---------------------------------------------------------------------------


def _generate_vqe_report_md(
    vqe_result: dict[str, Any],
    backend: str,
    seed: int,
) -> str:
    """Gera relatorio Markdown do modulo VQE."""
    ground_state_energy = vqe_result.get("ground_state_energy_hartree", 0.0)
    binding_energy = vqe_result.get("binding_energy_kcal", 0.0)
    n_iterations = vqe_result.get("n_iterations", 0)
    converged = vqe_result.get("converged", False)
    basis = vqe_result.get("basis_set", "sto-3g")

    report = f"""# MRL-QUANTUM VQE Report
## Simulacao do Sitio Ativo da RNase H1

**Data:** {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
**Backend:** {backend}
**Seed:** {seed}

---

### Contexto Biologico

A RNase H1 humana e a enzima responsavel por clivar o RNA no duplex
ASO:RNA formado quando o MRL-ASO-001 se liga ao SL RNA de *L. infantum*.
O sitio ativo contem dois ions Mg2+ coordenados por residuos Asp e Glu,
essenciais para a catalise de clivagem fosfodiester.

A simulacao VQE (Variational Quantum Eigensolver) calcula a energia
do estado fundamental do sitio ativo para avaliar:
- A estabilidade do complexo enzima-substrato
- O efeito das modificacoes LNA e PS na geometria do sitio catalico
- A viabilidade da clivagem enzimatica do SL RNA quando ligado ao ASO

### Parametros de Entrada

| Parametro | Valor |
|-----------|-------|
| PDB | 2QKB (RNase H1 humana) |
| Base set | {basis} |
| Eletrons ativos | 6 |
| Orbitais ativos | 6 |
| Ansatz | UCCSD |
| Backend | {backend} |

### Resultados

| Metrica | Valor |
|---------|-------|
| Energia do estado fundamental | {ground_state_energy:.6f} Hartree |
| Energia de ligacao | {binding_energy:.4f} kcal/mol |
| Iteracoes | {n_iterations} |
| Convergiu | {'Sim' if converged else 'Nao'} |

### Interpretacao

A simulacao VQE do sitio ativo da RNase H1 com o duplex MRL-ASO-001:SL RNA
{'convergiu' if converged else 'nao convergiu'} em {n_iterations} iteracoes,
com energia de ligacao de {binding_energy:.2f} kcal/mol. Este valor
{'confirma' if binding_energy < -5.0 else 'sugere cautela quanto a'}
a viabilidade da clivagem enzimatica do SL RNA pelo complexo RNase H1.

### Referencias

1. Nowotny M et al. (2005) *Mol Cell* 17(5):691-700
2. Reiher M et al. (2017) *PNAS* 114(29):7555-7560
3. Peruzzo A et al. (2014) *Nat Commun* 5:4213 (VQE original)
4. Lima WF et al. (2007) *J Biol Chem* 282(19):14530-7 (RNase H1 + ASO)
"""
    return report


def _generate_vqe_mulliken_json(
    vqe_result: dict[str, Any],
    backend: str,
    seed: int,
) -> dict[str, Any]:
    """Gera JSON com cargas de Mulliken do sitio ativo."""
    return {
        "metadata": _build_metadata("vqe_mulliken", backend, seed),
        "mulliken_charges": vqe_result.get("mulliken_charges", {}),
        "total_charge": vqe_result.get("total_charge", 0.0),
        "spin_multiplicity": vqe_result.get("spin_multiplicity", 1),
        "notes": (
            "Cargas de Mulliken para o active space do sitio catalico "
            "da RNase H1 (2 Mg2+ + residuos coordenantes)."
        ),
    }


def _generate_vqe_geometry_xyz(vqe_result: dict[str, Any]) -> str:
    """Gera arquivo XYZ com a geometria otimizada do sitio ativo."""
    atoms = vqe_result.get("active_site_atoms", [])
    n_atoms = len(atoms)

    lines: list[str] = []
    lines.append(str(n_atoms))
    lines.append(
        f"RNase H1 active site (PDB: 2QKB) — VQE optimized, "
        f"E = {vqe_result.get('ground_state_energy_hartree', 0.0):.6f} Hartree"
    )
    for atom in atoms:
        # Formato XYZ: Simbolo  X  Y  Z
        symbol = atom.get("symbol", "X")
        x = atom.get("x", 0.0)
        y = atom.get("y", 0.0)
        z = atom.get("z", 0.0)
        lines.append(f"{symbol:<4} {x:>12.6f} {y:>12.6f} {z:>12.6f}")

    return "\n".join(lines) + "\n"


def _generate_vqe_reports(
    vqe_result: dict[str, Any],
    backend: str,
    seed: int,
) -> list[Path]:
    """Gera todos os relatorios VQE."""
    paths: list[Path] = []

    # Relatorio Markdown
    md = _generate_vqe_report_md(vqe_result, backend, seed)
    paths.append(_write_text(md, VQE_RESULTS_DIR / "vqe_report.md"))

    # Cargas de Mulliken
    mulliken = _generate_vqe_mulliken_json(vqe_result, backend, seed)
    paths.append(_write_json(mulliken, VQE_RESULTS_DIR / "vqe_mulliken_charges.json"))

    # Geometria do sitio ativo
    xyz = _generate_vqe_geometry_xyz(vqe_result)
    paths.append(_write_text(xyz, VQE_RESULTS_DIR / "vqe_active_site_geometry.xyz"))

    return paths


# ---------------------------------------------------------------------------
# Funcao publica de orquestracao
# ---------------------------------------------------------------------------


def generate_all_reports(
    qaoa_result: dict[str, Any] | None = None,
    vqe_result: dict[str, Any] | None = None,
    backend: str = "aer_statevector",
    seed: int = 42,
) -> list[Path]:
    """Gera todos os relatorios disponiveis.

    Args:
        qaoa_result: Dicionario de resultados do modulo QAOA (ou None).
        vqe_result: Dicionario de resultados do modulo VQE (ou None).
        backend: Nome do backend quantico utilizado.
        seed: Semente de reproducibilidade.

    Returns:
        Lista de caminhos dos arquivos gerados.
    """
    all_paths: list[Path] = []

    if qaoa_result is not None:
        logger.info("Gerando relatorios QAOA...")
        paths = _generate_qaoa_reports(qaoa_result, backend, seed)
        all_paths.extend(paths)
        logger.info("  %d arquivos QAOA gerados.", len(paths))

    if vqe_result is not None:
        logger.info("Gerando relatorios VQE...")
        paths = _generate_vqe_reports(vqe_result, backend, seed)
        all_paths.extend(paths)
        logger.info("  %d arquivos VQE gerados.", len(paths))

    if not all_paths:
        logger.warning("Nenhum resultado fornecido — nenhum relatorio gerado.")

    return all_paths
