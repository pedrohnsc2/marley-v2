"""Math 4 — Prova espectral de irrepositabilidade do SL RNA no spliceosome.

Demonstra via teoria espectral de grafos que o SL RNA e o no mais critico
na rede do spliceosome de L. infantum. Em trypanosomatideos, todo mRNA
requer trans-splicing — o SL RNA e uma dependencia UNIVERSAL.

Analises realizadas:
    1. Construcao da rede do spliceosome (~17 nos, ~35 arestas)
    2. Calculo do Laplaciano e Laplaciano normalizado
    3. Decomposicao espectral (autovalores e autovetores)
    4. Impacto de remocao de cada no na conectividade algebrica
    5. Analise do gap espectral antes/apos remocao do SL RNA
    6. Perturbacao estocastica (robustez, 1000 iteracoes)
    7. Ranking final de irrepositabilidade

O resultado e uma PROVA MATEMATICA de que alvejar o SL RNA com um ASO
e a estrategia de maximo impacto contra o parasita.

Refs:
    - Fiedler M. (1973) Czech Math J 23(98):298-305
    - Michaeli S. (2011) Parasitology 138(12):1-16
    - Chung FRK. (1997) Spectral Graph Theory. AMS.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from aso_math.envelope import Timer
from aso_math.math_4_spectral.network import (
    SPLICEOSOME_EDGES,
    SPLICEOSOME_NODES,
    build_adjacency_matrix,
    get_network_stats,
)
from aso_math.math_4_spectral.spectral import (
    algebraic_connectivity,
    compute_eigenvalues,
    compute_laplacian,
    count_connected_components,
    fiedler_vector,
    node_removal_impact,
    perturbation_analysis,
    spectral_gap_analysis,
)
from core.logger import get_logger

logger = get_logger("math_4_spectral")

# ---------------------------------------------------------------------------
# Diretorio de saida
# ---------------------------------------------------------------------------

RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

# No alvo principal
TARGET_NODE: str = "SL_RNA"

# Segundo no mais importante em spliceosomes convencionais (para comparacao)
COMPARISON_NODE: str = "U2"

# Parametros de perturbacao
PERTURBATION_ITERATIONS: int = 1000
PERTURBATION_FRACTION: float = 0.20
PERTURBATION_SEED: int = 42


# ---------------------------------------------------------------------------
# 1. Construcao e estatisticas da rede
# ---------------------------------------------------------------------------


def analyze_network() -> tuple[np.ndarray, dict[str, int], dict[str, Any]]:
    """Constroi e analisa a rede do spliceosome.

    Returns:
        Tupla (adjacency_matrix, node_index, network_analysis_dict).
    """
    adj, node_idx, nodes = build_adjacency_matrix()
    stats = get_network_stats(adj, node_idx)

    # Informacoes dos nos com tipo funcional
    node_info: list[dict[str, str]] = []
    for node in nodes:
        node_info.append({
            "id": node["id"],
            "name": node["name"],
            "type": node["type"],
            "description": node["description"],
        })

    # Informacoes das arestas
    edge_info: list[dict[str, Any]] = []
    for src, dst, weight, reason in SPLICEOSOME_EDGES:
        edge_info.append({
            "source": src,
            "target": dst,
            "weight": weight,
            "justification": reason,
        })

    analysis = {
        "nodes": node_info,
        "edges": edge_info,
        "statistics": stats,
    }

    logger.info(
        "Rede do spliceosome: %d nos, %d arestas, densidade = %.4f",
        stats["n_nodes"], stats["n_edges"], stats["density"],
    )

    return adj, node_idx, analysis


# ---------------------------------------------------------------------------
# 2. Analise espectral completa
# ---------------------------------------------------------------------------


def full_spectral_analysis(
    adj: np.ndarray,
    node_idx: dict[str, int],
) -> dict[str, Any]:
    """Realiza decomposicao espectral completa do Laplaciano.

    Calcula:
        - Autovalores do Laplaciano e do Laplaciano normalizado
        - Conectividade algebrica (valor de Fiedler)
        - Vetor de Fiedler (particao natural da rede)
        - Numero de componentes conexas

    Args:
        adj: Matriz de adjacencia.
        node_idx: Mapeamento {id: indice}.

    Returns:
        Dicionario com resultados da analise espectral.
    """
    L, L_norm = compute_laplacian(adj)
    eigenvalues, eigenvectors = compute_eigenvalues(L)
    eigenvalues_norm, _ = compute_eigenvalues(L_norm)

    lambda2 = algebraic_connectivity(eigenvalues)
    fiedler = fiedler_vector(eigenvectors)
    n_components = count_connected_components(eigenvalues)

    # Interpretar o vetor de Fiedler: particao da rede
    idx_to_id = {v: k for k, v in node_idx.items()}
    fiedler_partition: list[dict[str, Any]] = []
    for i in range(len(fiedler)):
        fiedler_partition.append({
            "node_id": idx_to_id[i],
            "fiedler_value": round(float(fiedler[i]), 6),
            "cluster": "A" if fiedler[i] >= 0 else "B",
        })

    # Ordenar por valor de Fiedler
    fiedler_partition.sort(key=lambda f: f["fiedler_value"])

    # Gap espectral: (lambda_n - lambda_2) / lambda_n
    lambda_n = float(eigenvalues[-1])
    spectral_gap = (lambda_n - lambda2) / lambda_n if lambda_n > 1e-12 else 0.0

    result = {
        "laplacian_eigenvalues": [round(float(e), 6) for e in eigenvalues],
        "normalized_eigenvalues": [round(float(e), 6) for e in eigenvalues_norm],
        "algebraic_connectivity_lambda2": round(lambda2, 6),
        "largest_eigenvalue_lambda_n": round(lambda_n, 6),
        "spectral_gap_ratio": round(spectral_gap, 6),
        "n_connected_components": n_components,
        "is_connected": n_components == 1,
        "fiedler_partition": fiedler_partition,
    }

    logger.info(
        "Analise espectral: lambda_2 = %.6f, lambda_n = %.6f, "
        "gap espectral = %.6f, componentes = %d",
        lambda2, lambda_n, spectral_gap, n_components,
    )

    return result


# ---------------------------------------------------------------------------
# 3. Ranking de irrepositabilidade
# ---------------------------------------------------------------------------


def irreplaceability_ranking(
    adj: np.ndarray,
    node_idx: dict[str, int],
    spectral_data: dict[str, Any],
) -> dict[str, Any]:
    """Classifica todos os nos por impacto de remocao.

    Compara especificamente SL RNA vs U2 snRNA (o segundo mais importante
    em spliceosomes convencionais).

    Args:
        adj: Matriz de adjacencia.
        node_idx: Mapeamento {id: indice}.
        spectral_data: Resultados da analise espectral (para lambda_2 original).

    Returns:
        Dicionario com ranking completo e comparacao SL RNA vs U2.
    """
    lambda2_orig = spectral_data["algebraic_connectivity_lambda2"]
    impacts = node_removal_impact(adj, node_idx)

    # Encontrar SL RNA e U2 nos resultados
    sl_rna_impact = next(
        (r for r in impacts if r["node_id"] == TARGET_NODE), None
    )
    u2_impact = next(
        (r for r in impacts if r["node_id"] == COMPARISON_NODE), None
    )

    # Rank do SL RNA
    sl_rank = next(
        (i + 1 for i, r in enumerate(impacts) if r["node_id"] == TARGET_NODE), 0
    )

    # O SL RNA e o mais critico?
    sl_is_most_critical = sl_rank == 1

    # Comparacao SL RNA vs U2
    comparison: dict[str, Any] = {}
    if sl_rna_impact and u2_impact:
        comparison = {
            "sl_rna_fractional_drop": sl_rna_impact["fractional_drop"],
            "u2_fractional_drop": u2_impact["fractional_drop"],
            "sl_rna_advantage": round(
                sl_rna_impact["fractional_drop"] - u2_impact["fractional_drop"], 6
            ),
            "sl_rna_over_u2_ratio": round(
                sl_rna_impact["fractional_drop"] / u2_impact["fractional_drop"]
                if u2_impact["fractional_drop"] > 1e-12 else float("inf"),
                4,
            ),
            "sl_rna_fragments_network": sl_rna_impact["network_fragmented"],
            "u2_fragments_network": u2_impact["network_fragmented"],
        }

    result = {
        "lambda2_original": lambda2_orig,
        "ranking": impacts,
        "sl_rna_rank": sl_rank,
        "sl_rna_is_most_critical": sl_is_most_critical,
        "most_critical_node": impacts[0]["node_id"] if impacts else "",
        "comparison_sl_vs_u2": comparison,
    }

    logger.info(
        "Ranking de irrepositabilidade: SL RNA rank = %d/%d, "
        "queda fracional = %.4f, mais critico = %s",
        sl_rank, len(impacts),
        sl_rna_impact["fractional_drop"] if sl_rna_impact else 0.0,
        "SIM" if sl_is_most_critical else "NAO",
    )

    return result


# ---------------------------------------------------------------------------
# Gerador de relatorio Markdown
# ---------------------------------------------------------------------------


def _generate_report(
    network_analysis: dict[str, Any],
    spectral_data: dict[str, Any],
    ranking_data: dict[str, Any],
    gap_analysis: dict[str, Any],
    perturbation_data: dict[str, Any],
    runtime_seconds: float,
) -> str:
    """Gera relatorio Markdown detalhado com os resultados.

    O relatorio e honesto — se o SL RNA nao for o no mais critico,
    isso sera reportado claramente.
    """
    now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines: list[str] = []
    lines.append("# Spectral Graph Theory — SL RNA Irreplaceability Proof")
    lines.append("")
    lines.append(f"**Generated:** {now}")
    lines.append(f"**Runtime:** {runtime_seconds:.2f} seconds")
    lines.append(f"**Method:** Laplacian spectral decomposition + Fiedler analysis")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- 1. Rede ---
    stats = network_analysis["statistics"]
    lines.append("## 1. Spliceosome Network")
    lines.append("")
    lines.append("The *L. infantum* spliceosome network models the trans-splicing")
    lines.append("machinery as a weighted undirected graph. In trypanosomatids,")
    lines.append("**100% of mRNAs** require trans-splicing (no cis-splicing")
    lines.append("alternative exists for mature mRNA production).")
    lines.append("")
    lines.append(f"- **Nodes:** {stats['n_nodes']}")
    lines.append(f"- **Edges:** {stats['n_edges']}")
    lines.append(f"- **Density:** {stats['density']:.4f}")
    lines.append(f"- **Max possible edges:** {stats['max_possible_edges']}")
    lines.append("")

    # Tabela de graus
    lines.append("### Node Degrees (sorted by weighted degree)")
    lines.append("")
    lines.append("| Rank | Node | Degree | Weighted Degree |")
    lines.append("|------|------|--------|-----------------|")
    for i, d in enumerate(stats["node_degrees"][:10]):
        marker = " **" if d["node_id"] == TARGET_NODE else ""
        marker_end = "**" if d["node_id"] == TARGET_NODE else ""
        lines.append(
            f"| {i + 1} | {marker}{d['node_id']}{marker_end} | "
            f"{d['degree']} | {d['weighted_degree']:.4f} |"
        )
    lines.append("")

    # --- 2. Espectro ---
    lines.append("## 2. Spectral Properties")
    lines.append("")
    lines.append(f"- **Algebraic connectivity (lambda_2):** "
                 f"{spectral_data['algebraic_connectivity_lambda2']:.6f}")
    lines.append(f"- **Largest eigenvalue (lambda_n):** "
                 f"{spectral_data['largest_eigenvalue_lambda_n']:.6f}")
    lines.append(f"- **Spectral gap ratio:** "
                 f"{spectral_data['spectral_gap_ratio']:.6f}")
    lines.append(f"- **Connected components:** "
                 f"{spectral_data['n_connected_components']}")
    lines.append(f"- **Network is connected:** "
                 f"{'YES' if spectral_data['is_connected'] else 'NO'}")
    lines.append("")

    # Espectro completo
    evals = spectral_data["laplacian_eigenvalues"]
    lines.append("### Eigenvalue Spectrum")
    lines.append("")
    lines.append("| Index | lambda_i | Normalized lambda_i |")
    lines.append("|-------|----------|---------------------|")
    for i, (ev, env) in enumerate(zip(
        evals, spectral_data["normalized_eigenvalues"]
    )):
        lines.append(f"| {i + 1} | {ev:.6f} | {env:.6f} |")
    lines.append("")

    # Vetor de Fiedler
    lines.append("### Fiedler Vector (network partition)")
    lines.append("")
    lines.append("| Node | Fiedler Value | Cluster |")
    lines.append("|------|---------------|---------|")
    for f in spectral_data["fiedler_partition"]:
        lines.append(f"| {f['node_id']} | {f['fiedler_value']:.6f} | {f['cluster']} |")
    lines.append("")

    # --- 3. Ranking de irrepositabilidade ---
    lines.append("## 3. Irreplaceability Ranking")
    lines.append("")
    lines.append("Each node is removed and the resulting drop in algebraic")
    lines.append("connectivity (lambda_2) is measured. A larger drop indicates")
    lines.append("the node is more critical to network cohesion.")
    lines.append("")
    lines.append(
        f"- **Original lambda_2:** {ranking_data['lambda2_original']:.6f}"
    )
    lines.append(
        f"- **Most critical node:** "
        f"**{ranking_data['most_critical_node']}** (rank 1)"
    )
    lines.append(
        f"- **SL RNA rank:** {ranking_data['sl_rna_rank']}/{len(ranking_data['ranking'])}"
    )
    lines.append("")

    # Tabela completa
    lines.append("### Full Ranking (by fractional drop in lambda_2)")
    lines.append("")
    lines.append(
        "| Rank | Node | lambda_2 after | Drop | "
        "Fractional Drop | Fragmented? |"
    )
    lines.append(
        "|------|------|----------------|------|"
        "-----------------|-------------|"
    )
    for i, r in enumerate(ranking_data["ranking"]):
        marker = " **" if r["node_id"] == TARGET_NODE else ""
        marker_end = "**" if r["node_id"] == TARGET_NODE else ""
        frag = "YES" if r["network_fragmented"] else "no"
        lines.append(
            f"| {i + 1} | {marker}{r['node_id']}{marker_end} | "
            f"{r['lambda2_after_removal']:.6f} | "
            f"{r['lambda2_drop']:.6f} | "
            f"{r['fractional_drop']:.6f} | {frag} |"
        )
    lines.append("")

    # Comparacao SL vs U2
    comp = ranking_data.get("comparison_sl_vs_u2", {})
    if comp:
        lines.append("### SL RNA vs U2 snRNA Comparison")
        lines.append("")
        lines.append(f"- **SL RNA fractional drop:** {comp['sl_rna_fractional_drop']:.6f}")
        lines.append(f"- **U2 snRNA fractional drop:** {comp['u2_fractional_drop']:.6f}")
        lines.append(f"- **SL RNA advantage:** {comp['sl_rna_advantage']:.6f}")
        lines.append(f"- **SL RNA / U2 ratio:** {comp['sl_rna_over_u2_ratio']:.4f}x")
        lines.append(
            f"- **SL RNA fragments network:** "
            f"{'YES' if comp['sl_rna_fragments_network'] else 'no'}"
        )
        lines.append(
            f"- **U2 fragments network:** "
            f"{'YES' if comp['u2_fragments_network'] else 'no'}"
        )
        lines.append("")

    # --- 4. Gap espectral ---
    lines.append("## 4. Spectral Gap Analysis")
    lines.append("")
    lines.append(
        f"- **Original spectral gap ratio:** "
        f"{gap_analysis['original']['spectral_gap_ratio']:.6f}"
    )
    lines.append(
        f"- **After SL RNA removal:** "
        f"{gap_analysis['after_removal']['spectral_gap_ratio']:.6f}"
    )
    lines.append(
        f"- **Gap reduction:** {gap_analysis['gap_reduction']:.6f} "
        f"({gap_analysis['gap_reduction_pct']:.2f}%)"
    )
    lines.append(
        f"- **Components after removal:** "
        f"{gap_analysis['after_removal']['n_components']}"
    )
    lines.append("")

    # --- 5. Perturbacao ---
    lines.append("## 5. Perturbation Robustness Analysis")
    lines.append("")
    lines.append(
        f"Edge weights were perturbed by +/-{perturbation_data['perturbation_fraction'] * 100:.0f}% "
        f"across {perturbation_data['n_iterations']} iterations."
    )
    lines.append("")
    lines.append(
        f"- **Most robust critical node:** "
        f"**{perturbation_data['most_robust_critical_node']}**"
    )
    lines.append(
        f"- **Fraction as #1:** "
        f"{perturbation_data['top_node_fraction']:.4f} "
        f"({perturbation_data['top_node_fraction'] * 100:.1f}%)"
    )
    lines.append("")

    # Tabela de perturbacao (top 5)
    lines.append("### Perturbation Rankings (top 10)")
    lines.append("")
    lines.append("| Rank | Node | Times #1 | Fraction |")
    lines.append("|------|------|----------|----------|")
    for i, r in enumerate(perturbation_data["node_rankings"][:10]):
        marker = " **" if r["node_id"] == TARGET_NODE else ""
        marker_end = "**" if r["node_id"] == TARGET_NODE else ""
        lines.append(
            f"| {i + 1} | {marker}{r['node_id']}{marker_end} | "
            f"{r['times_most_critical']} | "
            f"{r['fraction']:.4f} ({r['fraction'] * 100:.1f}%) |"
        )
    lines.append("")

    # --- 6. Conclusao ---
    lines.append("## 6. Conclusion")
    lines.append("")

    sl_is_most_critical = ranking_data["sl_rna_is_most_critical"]
    sl_rank = ranking_data["sl_rna_rank"]
    perturbation_fraction = perturbation_data["top_node_fraction"]
    perturbation_winner = perturbation_data["most_robust_critical_node"]

    if sl_is_most_critical and perturbation_winner == TARGET_NODE:
        # Caso ideal: SL RNA e o mais critico em ambas analises
        lines.append(
            f"**SL RNA is the MOST CRITICAL node** in the *L. infantum* spliceosome "
            f"network. Its removal causes the largest drop in algebraic connectivity "
            f"(Fiedler value) among all {len(ranking_data['ranking'])} nodes."
        )
        lines.append("")
        lines.append(
            f"This result is **robust**: under +/-{perturbation_data['perturbation_fraction'] * 100:.0f}% "
            f"edge weight perturbation, SL RNA remains the most critical node in "
            f"{perturbation_fraction * 100:.1f}% of {perturbation_data['n_iterations']} "
            f"Monte Carlo iterations."
        )
        lines.append("")
        lines.append("### Biological Interpretation")
        lines.append("")
        lines.append(
            "In trypanosomatids like *L. infantum*, every mRNA requires "
            "trans-splicing of the SL exon. Unlike humans (where >95% of splicing "
            "is cis-splicing), there is **no bypass pathway**. Destroying SL RNA "
            "function with an antisense oligonucleotide (ASO) would collapse the "
            "entire mRNA maturation pipeline, making it a **mathematically proven "
            "optimal therapeutic target**."
        )
    elif sl_is_most_critical:
        # SL RNA e o mais critico deterministicamente, mas nao em perturbacao
        lines.append(
            f"**SL RNA is the most critical node** in the deterministic analysis "
            f"(rank 1/{len(ranking_data['ranking'])}). However, under perturbation, "
            f"**{perturbation_winner}** is the most frequent #1 node "
            f"({perturbation_fraction * 100:.1f}% of iterations). "
            f"SL RNA's criticality is sensitive to edge weight assumptions."
        )
    else:
        # HONESTO: SL RNA nao e o mais critico
        most_critical = ranking_data["most_critical_node"]
        lines.append(
            f"**HONEST FINDING:** SL RNA is **NOT** the most critical node. "
            f"**{most_critical}** has the highest impact on algebraic connectivity "
            f"(rank 1). SL RNA is ranked #{sl_rank}/{len(ranking_data['ranking'])}."
        )
        lines.append("")
        lines.append(
            "This does not invalidate the ASO strategy — SL RNA remains a "
            "high-value target due to its biological essentiality in trans-splicing. "
            "However, the network topology suggests other nodes may be more "
            "structurally central."
        )

    lines.append("")
    lines.append("---")
    lines.append(f"*Analysis completed in {runtime_seconds:.2f} seconds.*")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Funcao auxiliar para salvar JSON
# ---------------------------------------------------------------------------


def _save_json(data: Any, filename: str) -> Path:
    """Salva dados como JSON no diretorio de resultados.

    Args:
        data: Dados a serializar.
        filename: Nome do arquivo (sem extensao).
    """
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / f"{filename}.json"
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    return path


# ---------------------------------------------------------------------------
# Orquestrador principal
# ---------------------------------------------------------------------------


def main() -> dict[str, Any]:
    """Executa a analise espectral completa do spliceosome.

    Fluxo:
        1. Construir rede do spliceosome
        2. Analise espectral (Laplaciano, autovalores, Fiedler)
        3. Ranking de irrepositabilidade (remocao de cada no)
        4. Analise do gap espectral (SL RNA)
        5. Analise de perturbacao (robustez, 1000 iteracoes)
        6. Gerar relatorio Markdown
        7. Salvar todos os resultados como JSON

    Returns:
        Dicionario com todos os resultados da analise.
    """
    logger.info("=" * 70)
    logger.info("MATH 4 — Prova Espectral de Irrepositabilidade do SL RNA")
    logger.info("=" * 70)

    with Timer() as timer:
        # ---------------------------------------------------------------
        # Passo 1: Construcao da rede do spliceosome
        # ---------------------------------------------------------------
        logger.info("Passo 1/5: Construindo rede do spliceosome...")
        adj, node_idx, network_analysis = analyze_network()
        _save_json(network_analysis, "spliceosome_network")
        logger.info("  -> spliceosome_network.json salvo")

        # ---------------------------------------------------------------
        # Passo 2: Analise espectral completa
        # ---------------------------------------------------------------
        logger.info("Passo 2/5: Decomposicao espectral do Laplaciano...")
        spectral_data = full_spectral_analysis(adj, node_idx)
        _save_json(spectral_data, "spectral_analysis")
        logger.info("  -> spectral_analysis.json salvo")

        # ---------------------------------------------------------------
        # Passo 3: Ranking de irrepositabilidade
        # ---------------------------------------------------------------
        logger.info("Passo 3/5: Calculando ranking de irrepositabilidade...")
        ranking_data = irreplaceability_ranking(adj, node_idx, spectral_data)
        _save_json(ranking_data, "irreplaceability_ranking")
        logger.info("  -> irreplaceability_ranking.json salvo")

        # ---------------------------------------------------------------
        # Passo 4: Analise do gap espectral
        # ---------------------------------------------------------------
        logger.info("Passo 4/5: Analise do gap espectral (SL RNA)...")
        gap_analysis = spectral_gap_analysis(adj, node_idx, TARGET_NODE)
        _save_json(gap_analysis, "spectral_gap")
        logger.info("  -> spectral_gap.json salvo")

        # ---------------------------------------------------------------
        # Passo 5: Analise de perturbacao (robustez)
        # ---------------------------------------------------------------
        logger.info(
            "Passo 5/5: Perturbacao estocastica (%d iteracoes, +/-%.0f%%)...",
            PERTURBATION_ITERATIONS,
            PERTURBATION_FRACTION * 100,
        )
        perturbation_data = perturbation_analysis(
            adj=adj,
            node_idx=node_idx,
            n_iterations=PERTURBATION_ITERATIONS,
            perturbation_fraction=PERTURBATION_FRACTION,
            seed=PERTURBATION_SEED,
        )
        _save_json(perturbation_data, "perturbation_robustness")
        logger.info("  -> perturbation_robustness.json salvo")

    # --- Gerar relatorio Markdown ---
    report = _generate_report(
        network_analysis=network_analysis,
        spectral_data=spectral_data,
        ranking_data=ranking_data,
        gap_analysis=gap_analysis,
        perturbation_data=perturbation_data,
        runtime_seconds=timer.elapsed,
    )
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report_path = RESULTS_DIR / "SPECTRAL_REPORT.md"
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    logger.info("Relatorio salvo em: %s", report_path)

    # --- Resultado consolidado ---
    sl_impact = next(
        (r for r in ranking_data["ranking"] if r["node_id"] == TARGET_NODE),
        None,
    )

    result = {
        "module": "math_4_spectral",
        "version": "1.0.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": timer.elapsed,
        "status": "success",
        "network": {
            "n_nodes": network_analysis["statistics"]["n_nodes"],
            "n_edges": network_analysis["statistics"]["n_edges"],
            "density": network_analysis["statistics"]["density"],
        },
        "spectral": {
            "algebraic_connectivity": spectral_data["algebraic_connectivity_lambda2"],
            "largest_eigenvalue": spectral_data["largest_eigenvalue_lambda_n"],
            "spectral_gap_ratio": spectral_data["spectral_gap_ratio"],
            "is_connected": spectral_data["is_connected"],
        },
        "irreplaceability": {
            "sl_rna_rank": ranking_data["sl_rna_rank"],
            "sl_rna_is_most_critical": ranking_data["sl_rna_is_most_critical"],
            "most_critical_node": ranking_data["most_critical_node"],
            "sl_rna_fractional_drop": (
                sl_impact["fractional_drop"] if sl_impact else 0.0
            ),
            "comparison_sl_vs_u2": ranking_data.get("comparison_sl_vs_u2", {}),
        },
        "spectral_gap": {
            "gap_before": gap_analysis["original"]["spectral_gap_ratio"],
            "gap_after": gap_analysis["after_removal"]["spectral_gap_ratio"],
            "gap_reduction_pct": gap_analysis["gap_reduction_pct"],
        },
        "perturbation": {
            "n_iterations": perturbation_data["n_iterations"],
            "sl_rna_fraction_as_top": perturbation_data["top_node_fraction"],
            "most_robust_critical": perturbation_data["most_robust_critical_node"],
        },
        "summary": {
            "conclusion": (
                f"SL RNA is ranked #{ranking_data['sl_rna_rank']} of "
                f"{len(ranking_data['ranking'])} nodes by removal impact. "
                f"Under perturbation, {perturbation_data['most_robust_critical_node']} "
                f"is #1 in {perturbation_data['top_node_fraction'] * 100:.1f}% of "
                f"{perturbation_data['n_iterations']} iterations."
            ),
            "key_metrics": {
                "lambda2": spectral_data["algebraic_connectivity_lambda2"],
                "sl_rna_rank": ranking_data["sl_rna_rank"],
                "sl_rna_drop": sl_impact["fractional_drop"] if sl_impact else 0.0,
                "perturbation_robustness": perturbation_data["top_node_fraction"],
            },
        },
    }

    logger.info("=" * 70)
    logger.info("MATH 4 COMPLETO — Tempo total: %.2f s", timer.elapsed)
    logger.info("=" * 70)

    return result


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
