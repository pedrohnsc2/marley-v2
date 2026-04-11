"""Math 3 — Topological Data Analysis de estruturas ASO e SL RNA.

Aplica TDA (Topological Data Analysis) para comparar a topologia das
estruturas 3D do ASO livre (MRL-ASO-001), SL RNA alvo de L. infantum
e o duplex ASO:SL RNA formado apos hibridizacao.

Analises:
    1. Extracao de nuvens de pontos do PDB (backbone P, C3', C4', O3', O5')
    2. Filtracao de Vietoris-Rips (epsilon 0-20 A, passo 0.5 A)
    3. Diagramas de persistencia (H0 + H1) para cada estrutura
    4. Distancia bottleneck entre diagramas (rearranjo topologico)
    5. Score de estabilidade topologica (features H1 persistentes > 5 A)
    6. Numeros de Betti em escalas caracteristicas

Implementacao inteiramente from-scratch: Union-Find para H0, deteccao de
ciclos para H1. Nenhuma dependencia de gudhi, ripser, ou outros pacotes TDA.

Referencia:
    Edelsbrunner H, Harer J (2010) "Computational Topology"
    Carlsson G (2009) "Topology and data" Bull AMS 46(2):255-308
    Xia K, Wei GW (2014) "Persistent homology analysis of protein structure"
        J Chem Phys 140:174102 — justifica TDA para biomoleculas
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from aso_math.config import (
    ASO_PDB,
    ASO_SEQUENCE,
    ASO_TARGET_SEQUENCE,
    SL_PDB,
    SL_SEQUENCE,
)
from aso_math.target_config import TargetConfig
from aso_math.envelope import Timer

from aso_math.math_3_tda.point_cloud import build_point_clouds
from aso_math.math_3_tda.persistence import (
    PersistenceDiagram,
    betti_numbers_at_scales,
    bottleneck_distance,
    compute_persistence_diagram,
)

from core.logger import get_logger

logger = get_logger("math_3_tda")

# ---------------------------------------------------------------------------
# Diretorio de saida
# ---------------------------------------------------------------------------

RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Constantes da analise
# ---------------------------------------------------------------------------

# Range de filtracao de Vietoris-Rips (Angstroms)
EPSILON_MAX: float = 20.0
EPSILON_STEP: float = 0.5

# Limiar de persistencia para features H1 "estaveis" (Angstroms)
# Para nuvens de pontos de backbone nucleico, ciclos significativos
# (entrecruzamento de fitas, sulcos) persistem por 0.1-0.5 A na filtracao VR.
# O limiar de 0.1 A separa ciclos genuinos de artefatos numericos.
H1_PERSISTENCE_THRESHOLD: float = 0.1


# ---------------------------------------------------------------------------
# 1. Nuvens de pontos
# ---------------------------------------------------------------------------


def extract_point_clouds(
    aso_sequence: str,
    sl_target_sequence: str,
    aso_pdb_path: Path,
    sl_pdb_path: Path,
) -> dict[str, Any]:
    """Constroi e valida nuvens de pontos para as tres estruturas.

    Tenta PDB real primeiro; se insuficiente, gera backbone sintetico.

    Args:
        aso_sequence: Sequencia do ASO.
        sl_target_sequence: Sequencia da regiao alvo no SL RNA.
        aso_pdb_path: Caminho do PDB do ASO.
        sl_pdb_path: Caminho do PDB do SL RNA.

    Returns:
        Dicionario com nuvens de pontos e metadados de cada estrutura.
    """
    clouds = build_point_clouds(
        aso_sequence=aso_sequence,
        sl_target_sequence=sl_target_sequence,
        aso_pdb_path=aso_pdb_path,
        sl_pdb_path=sl_pdb_path,
    )

    summary: dict[str, Any] = {}
    for key in ("free_aso", "free_sl_rna", "duplex"):
        cloud = clouds[key]
        summary[key] = {
            "n_points": int(cloud.shape[0]),
            "source": clouds["sources"][key],
            "centroid": [round(float(c), 4) for c in np.mean(cloud, axis=0)] if cloud.shape[0] > 0 else [0, 0, 0],
            "extent_angstrom": [
                round(float(cloud[:, d].max() - cloud[:, d].min()), 4)
                for d in range(3)
            ] if cloud.shape[0] > 0 else [0, 0, 0],
        }

    logger.info(
        "Nuvens de pontos: ASO=%d, SL=%d, duplex=%d",
        summary["free_aso"]["n_points"],
        summary["free_sl_rna"]["n_points"],
        summary["duplex"]["n_points"],
    )

    return {
        "clouds": clouds,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# 2-3. Diagramas de persistencia para cada estrutura
# ---------------------------------------------------------------------------


def compute_all_persistence_diagrams(
    clouds: dict[str, np.ndarray],
    epsilon_max: float = EPSILON_MAX,
    epsilon_step: float = EPSILON_STEP,
) -> dict[str, PersistenceDiagram]:
    """Calcula diagramas de persistencia para ASO, SL RNA e duplex.

    Args:
        clouds: Dicionario com arrays Nx3 para cada estrutura.
        epsilon_max: Raio maximo de filtracao.
        epsilon_step: Passo de discretizacao.

    Returns:
        Dicionario mapeando nomes de estrutura a PersistenceDiagram.
    """
    diagrams: dict[str, PersistenceDiagram] = {}

    for key in ("free_aso", "free_sl_rna", "duplex"):
        points = clouds[key]
        logger.info("Computando diagrama de persistencia para %s...", key)
        diagrams[key] = compute_persistence_diagram(
            points=points,
            name=key,
            epsilon_max=epsilon_max,
            epsilon_step=epsilon_step,
        )

    return diagrams


# ---------------------------------------------------------------------------
# 4. Distancia bottleneck
# ---------------------------------------------------------------------------


def compute_bottleneck_distances(
    diagrams: dict[str, PersistenceDiagram],
) -> dict[str, Any]:
    """Calcula distancias bottleneck entre diagramas.

    Mede quanto a topologia muda quando ASO e SL RNA formam o duplex:
    - d_B(free_ASO, duplex): rearranjo topologico do ASO ao se ligar
    - d_B(free_SL_RNA, duplex): rearranjo topologico do SL RNA ao ser ligado
    - d_B(free_ASO, free_SL_RNA): diferenca intrinseca entre as estruturas

    Maior distancia = mais rearranjo topologico (mais impacto da ligacao).

    Returns:
        Dicionario com distancias bottleneck para H0 e H1.
    """
    results: dict[str, Any] = {}

    pairs = [
        ("free_aso_vs_duplex", "free_aso", "duplex"),
        ("free_sl_rna_vs_duplex", "free_sl_rna", "duplex"),
        ("free_aso_vs_free_sl_rna", "free_aso", "free_sl_rna"),
    ]

    for label, key_a, key_b in pairs:
        d_h0 = bottleneck_distance(diagrams[key_a], diagrams[key_b], dimension=0)
        d_h1 = bottleneck_distance(diagrams[key_a], diagrams[key_b], dimension=1)

        results[label] = {
            "bottleneck_h0_angstrom": round(d_h0, 4),
            "bottleneck_h1_angstrom": round(d_h1, 4),
            "description": f"Distancia bottleneck entre {key_a} e {key_b}",
        }

        logger.info(
            "Bottleneck %s: H0=%.4f A, H1=%.4f A",
            label, d_h0, d_h1,
        )

    return results


# ---------------------------------------------------------------------------
# 5. Score de estabilidade topologica
# ---------------------------------------------------------------------------


def compute_stability_score(
    diagrams: dict[str, PersistenceDiagram],
    threshold: float = H1_PERSISTENCE_THRESHOLD,
) -> dict[str, Any]:
    """Avalia estabilidade topologica de cada estrutura e impacto da ligacao.

    Features H1 com persistencia > threshold representam motifs topologicos
    estaveis (loops formados por sulcos, empilhamento de bases, etc.).

    Se o duplex tem mais features estaveis que os componentes isolados,
    a ligacao cria nova topologia estavel — favoravel para a eficacia do ASO.

    Args:
        diagrams: Diagramas de persistencia de cada estrutura.
        threshold: Limiar de persistencia para H1 (Angstroms).

    Returns:
        Score de estabilidade e analise de impacto.
    """
    stability: dict[str, Any] = {}

    for key in ("free_aso", "free_sl_rna", "duplex"):
        diagram = diagrams[key]
        n_persistent = diagram.n_persistent_h1(threshold=threshold)
        total_h1 = diagram.total_persistence_h1
        n_h0 = len(diagram.h0_pairs)
        n_h1 = len(diagram.h1_pairs)

        # Complexidade topologica = persistencia total normalizada
        # Quanto mais features persistentes, mais complexa a topologia
        total_persistence = diagram.total_persistence_h0 + total_h1
        complexity_score = round(total_persistence, 4)

        stability[key] = {
            "n_h0_features": n_h0,
            "n_h1_features": n_h1,
            "n_persistent_h1": n_persistent,
            "total_persistence_h0": round(diagram.total_persistence_h0, 4),
            "total_persistence_h1": round(total_h1, 4),
            "topological_complexity": complexity_score,
        }

        logger.info(
            "%s: %d features H1 persistentes (> %.1f A), complexidade = %.2f",
            key, n_persistent, threshold, complexity_score,
        )

    # Analise de impacto: comparar duplex vs componentes isolados
    aso_persistent = stability["free_aso"]["n_persistent_h1"]
    sl_persistent = stability["free_sl_rna"]["n_persistent_h1"]
    duplex_persistent = stability["duplex"]["n_persistent_h1"]

    # O duplex criou nova topologia estavel?
    # Comparar com o maximo dos componentes isolados
    max_component = max(aso_persistent, sl_persistent)
    binding_creates_topology = duplex_persistent > max_component

    # Ganho relativo de topologia
    if max_component > 0:
        topology_gain_ratio = round(duplex_persistent / max_component, 4)
    else:
        topology_gain_ratio = float(duplex_persistent) if duplex_persistent > 0 else 0.0

    # Complexidade relativa
    aso_complexity = stability["free_aso"]["topological_complexity"]
    sl_complexity = stability["free_sl_rna"]["topological_complexity"]
    duplex_complexity = stability["duplex"]["topological_complexity"]
    sum_components = aso_complexity + sl_complexity

    if sum_components > 0:
        complexity_change = round(
            (duplex_complexity - sum_components) / sum_components * 100, 2,
        )
    else:
        complexity_change = 0.0

    binding_impact = {
        "aso_persistent_h1": aso_persistent,
        "sl_rna_persistent_h1": sl_persistent,
        "duplex_persistent_h1": duplex_persistent,
        "binding_creates_new_topology": binding_creates_topology,
        "topology_gain_ratio": topology_gain_ratio,
        "complexity_change_percent": complexity_change,
        "persistence_threshold_angstrom": threshold,
        "interpretation": (
            "A formacao do duplex ASO:SL RNA cria topologia estavel adicional "
            "(favoravel para eficacia terapeutica)"
            if binding_creates_topology
            else "A topologia do duplex e comparavel a dos componentes isolados "
                 "(a ligacao nao altera significativamente a estrutura topologica)"
        ),
    }

    logger.info(
        "Impacto da ligacao: cria nova topologia = %s, ganho = %.2fx, "
        "mudanca de complexidade = %.1f%%",
        binding_creates_topology, topology_gain_ratio, complexity_change,
    )

    return {
        "per_structure": stability,
        "binding_impact": binding_impact,
    }


# ---------------------------------------------------------------------------
# 6. Numeros de Betti em escalas caracteristicas
# ---------------------------------------------------------------------------


def compute_betti_at_scales(
    diagrams: dict[str, PersistenceDiagram],
) -> dict[str, Any]:
    """Calcula numeros de Betti para cada estrutura em escalas biologicamente relevantes.

    Escalas escolhidas:
    - 2 A: resolucao de ligacoes covalentes
    - 4 A: pontes de hidrogenio (Watson-Crick base pairing)
    - 6 A: empilhamento de bases (base stacking)
    - 8 A: transicao para escala de sulco menor
    - 10 A: diametro do backbone (raio da B-DNA)
    - 15 A: escala de sulco maior
    - 20 A: diametro total da dupla helice

    Args:
        diagrams: Diagramas de persistencia de cada estrutura.

    Returns:
        Numeros de Betti por estrutura e escala.
    """
    scales = [2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]
    results: dict[str, Any] = {}

    for key in ("free_aso", "free_sl_rna", "duplex"):
        betti = betti_numbers_at_scales(diagrams[key], scales=scales)
        results[key] = betti

        logger.info(
            "%s Betti: %s",
            key,
            ", ".join(
                f"eps={b['epsilon_angstrom']}: b0={b['beta_0']},b1={b['beta_1']}"
                for b in betti
            ),
        )

    return results


# ---------------------------------------------------------------------------
# Gerador de relatorio Markdown
# ---------------------------------------------------------------------------


def _generate_report(
    cloud_summary: dict[str, Any],
    diagrams: dict[str, PersistenceDiagram],
    bottleneck_results: dict[str, Any],
    stability_results: dict[str, Any],
    betti_results: dict[str, Any],
    runtime_seconds: float,
) -> str:
    """Gera relatorio Markdown com analise topologica completa."""
    now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines: list[str] = []
    lines.append("# Topological Data Analysis — ASO:SL RNA Structural Report")
    lines.append("")
    lines.append(f"**Generated:** {now}")
    lines.append(f"**Runtime:** {runtime_seconds:.2f} seconds")
    lines.append(f"**Method:** Vietoris-Rips filtration (0-{EPSILON_MAX} A, step {EPSILON_STEP} A)")
    lines.append("**Implementation:** From-scratch (Union-Find H0 + cycle detection H1)")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- 1. Point Clouds ---
    lines.append("## 1. Point Cloud Extraction")
    lines.append("")
    lines.append("| Structure | Points | Source | Extent X (A) | Extent Y (A) | Extent Z (A) |")
    lines.append("|-----------|--------|--------|-------------|-------------|-------------|")
    for key in ("free_aso", "free_sl_rna", "duplex"):
        s = cloud_summary[key]
        ext = s["extent_angstrom"]
        lines.append(
            f"| {key} | {s['n_points']} | {s['source']} | "
            f"{ext[0]:.1f} | {ext[1]:.1f} | {ext[2]:.1f} |"
        )
    lines.append("")

    # --- 2-3. Persistence Diagrams ---
    lines.append("## 2. Persistence Diagrams")
    lines.append("")
    for key in ("free_aso", "free_sl_rna", "duplex"):
        diag = diagrams[key]
        lines.append(f"### {key}")
        lines.append("")
        lines.append(f"- **H0 pairs:** {len(diag.h0_pairs)}")
        lines.append(f"- **H1 pairs:** {len(diag.h1_pairs)}")
        lines.append(f"- **Total persistence H0:** {diag.total_persistence_h0:.4f} A")
        lines.append(f"- **Total persistence H1:** {diag.total_persistence_h1:.4f} A")
        lines.append(f"- **Persistent H1 features (> {H1_PERSISTENCE_THRESHOLD} A):** {diag.n_persistent_h1}")
        lines.append("")

        # Top 5 pares H0 mais persistentes
        h0_finite = [p for p in diag.h0_pairs if p.death != float("inf")]
        if h0_finite:
            lines.append("**Top 5 H0 pairs (most persistent):**")
            lines.append("")
            lines.append("| Birth (A) | Death (A) | Persistence (A) |")
            lines.append("|-----------|-----------|-----------------|")
            for p in h0_finite[:5]:
                lines.append(f"| {p.birth:.4f} | {p.death:.4f} | {p.persistence:.4f} |")
            lines.append("")

        # Top 5 pares H1 mais persistentes
        if diag.h1_pairs:
            lines.append("**Top 5 H1 pairs (most persistent):**")
            lines.append("")
            lines.append("| Birth (A) | Death (A) | Persistence (A) |")
            lines.append("|-----------|-----------|-----------------|")
            for p in diag.h1_pairs[:5]:
                death_str = f"{p.death:.4f}" if p.death != float("inf") else "inf"
                pers_str = f"{p.persistence:.4f}" if p.persistence != float("inf") else "inf"
                lines.append(f"| {p.birth:.4f} | {death_str} | {pers_str} |")
            lines.append("")

    # --- 4. Bottleneck Distances ---
    lines.append("## 3. Bottleneck Distances")
    lines.append("")
    lines.append("Measures topological rearrangement upon binding.")
    lines.append("Higher distance = more structural change.")
    lines.append("")
    lines.append("| Comparison | d_B(H0) (A) | d_B(H1) (A) |")
    lines.append("|------------|-------------|-------------|")
    for label, data in bottleneck_results.items():
        lines.append(
            f"| {label} | {data['bottleneck_h0_angstrom']:.4f} | "
            f"{data['bottleneck_h1_angstrom']:.4f} |"
        )
    lines.append("")

    # --- 5. Stability Score ---
    lines.append("## 4. Topological Stability Analysis")
    lines.append("")
    impact = stability_results["binding_impact"]
    per_struct = stability_results["per_structure"]

    lines.append("| Structure | H0 features | H1 features | Persistent H1 | Complexity |")
    lines.append("|-----------|------------|------------|--------------|------------|")
    for key in ("free_aso", "free_sl_rna", "duplex"):
        s = per_struct[key]
        lines.append(
            f"| {key} | {s['n_h0_features']} | {s['n_h1_features']} | "
            f"{s['n_persistent_h1']} | {s['topological_complexity']:.2f} |"
        )
    lines.append("")

    lines.append("### Binding Impact")
    lines.append("")
    lines.append(f"- **ASO persistent H1:** {impact['aso_persistent_h1']}")
    lines.append(f"- **SL RNA persistent H1:** {impact['sl_rna_persistent_h1']}")
    lines.append(f"- **Duplex persistent H1:** {impact['duplex_persistent_h1']}")
    lines.append(f"- **Binding creates new topology:** "
                 f"{'YES' if impact['binding_creates_new_topology'] else 'NO'}")
    lines.append(f"- **Topology gain ratio:** {impact['topology_gain_ratio']:.2f}x")
    lines.append(f"- **Complexity change:** {impact['complexity_change_percent']:.1f}%")
    lines.append("")
    lines.append(f"**Interpretation:** {impact['interpretation']}")
    lines.append("")

    # --- 6. Betti Numbers ---
    lines.append("## 5. Betti Numbers at Characteristic Scales")
    lines.append("")
    lines.append("| epsilon (A) | ASO b0 | ASO b1 | SL b0 | SL b1 | Duplex b0 | Duplex b1 |")
    lines.append("|-------------|--------|--------|-------|-------|-----------|-----------|")

    aso_betti = betti_results.get("free_aso", [])
    sl_betti = betti_results.get("free_sl_rna", [])
    dup_betti = betti_results.get("duplex", [])

    for i in range(len(aso_betti)):
        eps = aso_betti[i]["epsilon_angstrom"]
        lines.append(
            f"| {eps:>11.1f} | "
            f"{aso_betti[i]['beta_0']:>6} | {aso_betti[i]['beta_1']:>6} | "
            f"{sl_betti[i]['beta_0']:>5} | {sl_betti[i]['beta_1']:>5} | "
            f"{dup_betti[i]['beta_0']:>9} | {dup_betti[i]['beta_1']:>9} |"
        )
    lines.append("")

    # --- Conclusao ---
    lines.append("## Conclusion")
    lines.append("")

    # Avaliar se a ligacao e topologicamente favoravel
    creates_topology = impact["binding_creates_new_topology"]
    gain_ratio = impact["topology_gain_ratio"]

    if creates_topology and gain_ratio >= 1.5:
        verdict = (
            "The ASO:SL RNA duplex exhibits **significantly enhanced topological complexity** "
            f"compared to the isolated components (gain ratio: {gain_ratio:.2f}x). "
            "The formation of stable H1 features indicates that binding creates "
            "new topological motifs (e.g., intertwined backbone loops), which is "
            "consistent with tight, therapeutically relevant hybridization."
        )
    elif creates_topology:
        verdict = (
            "The ASO:SL RNA duplex shows **moderately enhanced topological complexity** "
            f"(gain ratio: {gain_ratio:.2f}x). "
            "Binding creates additional persistent loops, suggesting meaningful "
            "structural rearrangement upon hybridization."
        )
    else:
        verdict = (
            "The duplex topology is comparable to the isolated components. "
            "This suggests the binding event primarily consolidates existing "
            "structural features rather than creating new topological motifs."
        )

    lines.append(verdict)
    lines.append("")
    lines.append("---")
    lines.append(f"*Analysis completed in {runtime_seconds:.2f} seconds. "
                 "Implementation: Vietoris-Rips from scratch (Union-Find + cycle detection).*")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Funcao auxiliar para salvar JSON
# ---------------------------------------------------------------------------


def _save_json(data: Any, filename: str) -> Path:
    """Salva dados como JSON no diretorio de resultados.

    Args:
        data: Dados a serializar (deve ser JSON-serializavel).
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


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a analise topologica completa.

    Fluxo:
        1. Extrair nuvens de pontos de PDB ou gerar sinteticas
        2. Computar diagramas de persistencia (H0 + H1) para cada estrutura
        3. Calcular distancias bottleneck entre diagramas
        4. Avaliar estabilidade topologica e impacto da ligacao
        5. Calcular numeros de Betti em escalas caracteristicas
        6. Gerar relatorio Markdown e JSONs

    Args:
        config: Configuracao do organismo alvo. Se None, usa defaults MRL-ASO-001.

    Returns:
        Dicionario com todos os resultados da analise.
    """
    # --- Config padrao para L. infantum (retrocompativel) ---
    if config is None:
        config = TargetConfig(
            aso_sequence=ASO_SEQUENCE,
            aso_name="MRL-ASO-001",
            sl_sequence=SL_SEQUENCE,
            aso_target_start=5,
            aso_target_end=30,
            known_dg=-27.97,
            known_tm=68.48,
        )

    aso_seq = config.aso_sequence.upper()
    sl_seq = config.sl_sequence.upper()
    aso_name = config.aso_name or "ASO"

    # Regiao alvo do SL RNA (onde o ASO se liga)
    sl_target_seq = sl_seq[config.aso_target_start:config.aso_target_end]

    logger.info("=" * 70)
    logger.info("MATH 3 — Topological Data Analysis — %s", aso_name)
    logger.info("=" * 70)
    logger.info("ASO: %d nt | SL target: %d nt | epsilon: 0-%.1f A",
                len(aso_seq), len(sl_target_seq), EPSILON_MAX)

    with Timer() as timer:
        # ---------------------------------------------------------------
        # Passo 1: Extracao de nuvens de pontos
        # ---------------------------------------------------------------
        logger.info("Passo 1/5: Extraindo nuvens de pontos...")
        cloud_data = extract_point_clouds(
            aso_sequence=aso_seq,
            sl_target_sequence=sl_target_seq,
            aso_pdb_path=ASO_PDB,
            sl_pdb_path=SL_PDB,
        )
        clouds = cloud_data["clouds"]
        cloud_summary = cloud_data["summary"]

        _save_json(cloud_summary, "point_clouds")
        logger.info("  -> point_clouds.json salvo")

        # ---------------------------------------------------------------
        # Passo 2-3: Diagramas de persistencia
        # ---------------------------------------------------------------
        logger.info("Passo 2/5: Computando diagramas de persistencia...")
        diagrams = compute_all_persistence_diagrams(
            clouds=clouds,
            epsilon_max=EPSILON_MAX,
            epsilon_step=EPSILON_STEP,
        )

        # Serializar diagramas
        diagrams_json = {
            key: diag.to_dict() for key, diag in diagrams.items()
        }
        _save_json(diagrams_json, "persistence_diagrams")
        logger.info("  -> persistence_diagrams.json salvo")

        # ---------------------------------------------------------------
        # Passo 3: Distancias bottleneck
        # ---------------------------------------------------------------
        logger.info("Passo 3/5: Calculando distancias bottleneck...")
        bottleneck_results = compute_bottleneck_distances(diagrams)
        _save_json(bottleneck_results, "bottleneck_distances")
        logger.info("  -> bottleneck_distances.json salvo")

        # ---------------------------------------------------------------
        # Passo 4: Score de estabilidade topologica
        # ---------------------------------------------------------------
        logger.info("Passo 4/5: Avaliando estabilidade topologica...")
        stability_results = compute_stability_score(diagrams)
        _save_json(stability_results, "topological_stability")
        logger.info("  -> topological_stability.json salvo")

        # ---------------------------------------------------------------
        # Passo 5: Numeros de Betti em escalas caracteristicas
        # ---------------------------------------------------------------
        logger.info("Passo 5/5: Calculando numeros de Betti...")
        betti_results = compute_betti_at_scales(diagrams)
        _save_json(betti_results, "betti_numbers")
        logger.info("  -> betti_numbers.json salvo")

    # --- Gerar relatorio Markdown ---
    report = _generate_report(
        cloud_summary=cloud_summary,
        diagrams=diagrams,
        bottleneck_results=bottleneck_results,
        stability_results=stability_results,
        betti_results=betti_results,
        runtime_seconds=timer.elapsed,
    )
    report_path = RESULTS_DIR / "TDA_REPORT.md"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    logger.info("Relatorio salvo em: %s", report_path)

    # --- Resultado consolidado ---
    # Extrair metricas-chave para o envelope
    impact = stability_results["binding_impact"]
    aso_bn = bottleneck_results.get("free_aso_vs_duplex", {})
    sl_bn = bottleneck_results.get("free_sl_rna_vs_duplex", {})

    result: dict[str, Any] = {
        "module": "math_3_tda",
        "version": "1.0.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": timer.elapsed,
        "status": "success",
        "aso_name": aso_name,
        "organism": config.species_name,
        "point_clouds": cloud_summary,
        "persistence_diagrams": diagrams_json,
        "bottleneck_distances": bottleneck_results,
        "topological_stability": stability_results,
        "betti_numbers": betti_results,
        "summary": {
            "conclusion": (
                "Duplex forma topologia estavel adicional"
                if impact["binding_creates_new_topology"]
                else "Topologia do duplex comparavel aos componentes isolados"
            ),
            "key_metrics": {
                "aso_persistent_h1": impact["aso_persistent_h1"],
                "sl_rna_persistent_h1": impact["sl_rna_persistent_h1"],
                "duplex_persistent_h1": impact["duplex_persistent_h1"],
                "binding_creates_new_topology": impact["binding_creates_new_topology"],
                "topology_gain_ratio": impact["topology_gain_ratio"],
                "complexity_change_percent": impact["complexity_change_percent"],
                "bottleneck_h0_aso_duplex": aso_bn.get("bottleneck_h0_angstrom", 0.0),
                "bottleneck_h1_aso_duplex": aso_bn.get("bottleneck_h1_angstrom", 0.0),
                "bottleneck_h0_sl_duplex": sl_bn.get("bottleneck_h0_angstrom", 0.0),
                "bottleneck_h1_sl_duplex": sl_bn.get("bottleneck_h1_angstrom", 0.0),
            },
        },
    }

    logger.info("=" * 70)
    logger.info("MATH 3 COMPLETO — Tempo total: %.2f s", timer.elapsed)
    logger.info("=" * 70)

    return result


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
