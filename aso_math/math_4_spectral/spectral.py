"""Analise espectral do grafo do spliceosome.

Implementa decomposicao espectral do Laplaciano para quantificar a
criticidade de cada no na rede. A conectividade algebrica (lambda_2,
valor de Fiedler) mede a coesao da rede — uma queda acentuada apos
remocao de um no indica que esse no e estruturalmente irrepositavel.

Teoria:
    - Laplaciano L = D - A, onde D = diag(graus), A = adjacencia
    - Laplaciano normalizado L_norm = D^{-1/2} L D^{-1/2}
    - Autovalores 0 = lambda_1 <= lambda_2 <= ... <= lambda_n
    - lambda_2 (valor de Fiedler) = conectividade algebrica
    - Vetor de Fiedler: autovetor de lambda_2 — indica particao natural

Refs:
    - Fiedler M. (1973) Czech Math J 23(98):298-305
    - Chung FRK. (1997) Spectral Graph Theory. AMS.
    - Estrada E. (2012) The Structure of Complex Networks. Oxford.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import linalg as la


# ---------------------------------------------------------------------------
# Laplaciano
# ---------------------------------------------------------------------------


def compute_laplacian(adj: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calcula o Laplaciano e o Laplaciano normalizado do grafo.

    Args:
        adj: Matriz de adjacencia simetrica ponderada.

    Returns:
        Tupla (L, L_norm) onde L = D - A e L_norm = D^{-1/2} L D^{-1/2}.
    """
    degrees = np.sum(adj, axis=1)
    D = np.diag(degrees)
    L = D - adj

    # Laplaciano normalizado
    # Para nos com grau 0, evitamos divisao por zero
    d_inv_sqrt = np.zeros_like(degrees)
    nonzero = degrees > 1e-12
    d_inv_sqrt[nonzero] = 1.0 / np.sqrt(degrees[nonzero])

    D_inv_sqrt = np.diag(d_inv_sqrt)
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    return L, L_norm


def compute_eigenvalues(
    L: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Calcula autovalores e autovetores do Laplaciano.

    Usa scipy.linalg.eigh para matrizes simetricas reais (mais eficiente
    e numericamente estavel que eig generico).

    Args:
        L: Matriz Laplaciana (simetrica, semi-definida positiva).

    Returns:
        Tupla (eigenvalues, eigenvectors) ordenados crescentemente.
    """
    eigenvalues, eigenvectors = la.eigh(L)

    # Corrigir autovalores muito pequenos para zero (artefato numerico)
    eigenvalues[np.abs(eigenvalues) < 1e-10] = 0.0

    return eigenvalues, eigenvectors


def algebraic_connectivity(eigenvalues: np.ndarray) -> float:
    """Retorna lambda_2 — a conectividade algebrica (valor de Fiedler).

    lambda_2 > 0 indica que o grafo e conexo.
    Quanto maior lambda_2, mais robusto o grafo a remocao de arestas.

    Args:
        eigenvalues: Autovalores ordenados crescentemente.

    Returns:
        Valor de lambda_2. Retorna 0.0 se o grafo tem menos de 2 nos.
    """
    if len(eigenvalues) < 2:
        return 0.0
    return float(eigenvalues[1])


def fiedler_vector(eigenvectors: np.ndarray) -> np.ndarray:
    """Retorna o vetor de Fiedler (autovetor de lambda_2).

    O sinal dos componentes indica a particao natural do grafo em
    dois clusters. Nos com mesmo sinal tendem a estar mais conectados
    entre si.

    Args:
        eigenvectors: Matriz de autovetores (colunas).

    Returns:
        Vetor de Fiedler (segunda coluna dos autovetores).
    """
    if eigenvectors.shape[1] < 2:
        return np.zeros(eigenvectors.shape[0])
    return eigenvectors[:, 1]


def count_connected_components(eigenvalues: np.ndarray) -> int:
    """Conta componentes conexas pela multiplicidade de lambda = 0.

    O numero de autovalores iguais a zero corresponde ao numero de
    componentes conexas do grafo.

    Args:
        eigenvalues: Autovalores do Laplaciano.

    Returns:
        Numero de componentes conexas.
    """
    return int(np.sum(np.abs(eigenvalues) < 1e-8))


# ---------------------------------------------------------------------------
# Analise de impacto de remocao de no
# ---------------------------------------------------------------------------


def node_removal_impact(
    adj: np.ndarray,
    node_idx: dict[str, int],
) -> list[dict[str, Any]]:
    """Calcula o impacto da remocao de cada no na conectividade algebrica.

    Para cada no:
        1. Remove o no da matriz de adjacencia
        2. Recalcula o Laplaciano
        3. Calcula lambda_2 da rede reduzida
        4. Mede a queda fracional em lambda_2

    Args:
        adj: Matriz de adjacencia original.
        node_idx: Mapeamento {id: indice}.

    Returns:
        Lista de dicionarios com impacto de cada no, ordenada por impacto
        decrescente (no mais critico primeiro).
    """
    # Importar aqui para evitar dependencia circular
    from aso_math.math_4_spectral.network import remove_node

    # Conectividade algebrica da rede original
    L_orig, _ = compute_laplacian(adj)
    evals_orig, _ = compute_eigenvalues(L_orig)
    lambda2_orig = algebraic_connectivity(evals_orig)
    components_orig = count_connected_components(evals_orig)

    idx_to_id = {v: k for k, v in node_idx.items()}
    results: list[dict[str, Any]] = []

    for node_id in node_idx:
        # Remover no
        adj_reduced, _ = remove_node(adj, node_idx, node_id)

        # Recalcular espectro
        L_red, _ = compute_laplacian(adj_reduced)
        evals_red, _ = compute_eigenvalues(L_red)
        lambda2_red = algebraic_connectivity(evals_red)
        components_red = count_connected_components(evals_red)

        # Queda fracional em lambda_2
        if lambda2_orig > 1e-12:
            fractional_drop = (lambda2_orig - lambda2_red) / lambda2_orig
        else:
            fractional_drop = 0.0

        # Aumento no numero de componentes
        component_increase = components_red - (components_orig - 1)
        # Subtraimos 1 porque remover um no de um grafo conexo com
        # N nos gera um grafo com N-1 nos. Se continua conexo,
        # components_red = 1 e components_orig - 1 = 0, logo aumento = 1.
        # Realmente o que queremos medir e: a rede fragmentou?
        # components_orig (rede inteira): 1
        # Ao remover um no: esperamos 1 componente se nao fragmentou
        # Se fragmentou, components_red > 1
        component_increase = components_red - 1  # quantas fragmentacoes

        results.append({
            "node_id": node_id,
            "lambda2_after_removal": round(lambda2_red, 6),
            "lambda2_drop": round(lambda2_orig - lambda2_red, 6),
            "fractional_drop": round(fractional_drop, 6),
            "components_after_removal": components_red,
            "network_fragmented": components_red > 1,
            "fragments_created": component_increase,
        })

    # Ordenar por queda fracional (maior impacto primeiro)
    results.sort(key=lambda r: r["fractional_drop"], reverse=True)

    return results


def spectral_gap_analysis(
    adj: np.ndarray,
    node_idx: dict[str, int],
    target_node: str = "SL_RNA",
) -> dict[str, Any]:
    """Analisa o gap espectral antes e depois da remocao do no alvo.

    O gap espectral (lambda_n - lambda_2) / lambda_n indica a robustez
    do grafo. Uma reducao significativa no gap sugere que o no removido
    era essencial para a coerencia espectral.

    Args:
        adj: Matriz de adjacencia.
        node_idx: Mapeamento {id: indice}.
        target_node: Id do no alvo (default: SL_RNA).

    Returns:
        Dicionario com analise do gap espectral.
    """
    from aso_math.math_4_spectral.network import remove_node

    # Espectro original
    L_orig, L_norm_orig = compute_laplacian(adj)
    evals_orig, evecs_orig = compute_eigenvalues(L_orig)
    evals_norm_orig, _ = compute_eigenvalues(L_norm_orig)

    lambda2_orig = algebraic_connectivity(evals_orig)
    lambda_n_orig = float(evals_orig[-1])

    # Gap espectral original
    if lambda_n_orig > 1e-12:
        gap_orig = (lambda_n_orig - lambda2_orig) / lambda_n_orig
    else:
        gap_orig = 0.0

    # Espectro apos remocao
    adj_red, _ = remove_node(adj, node_idx, target_node)
    L_red, L_norm_red = compute_laplacian(adj_red)
    evals_red, _ = compute_eigenvalues(L_red)
    evals_norm_red, _ = compute_eigenvalues(L_norm_red)

    lambda2_red = algebraic_connectivity(evals_red)
    lambda_n_red = float(evals_red[-1]) if len(evals_red) > 0 else 0.0

    if lambda_n_red > 1e-12:
        gap_red = (lambda_n_red - lambda2_red) / lambda_n_red
    else:
        gap_red = 0.0

    return {
        "target_node": target_node,
        "original": {
            "eigenvalues": [round(float(e), 6) for e in evals_orig],
            "eigenvalues_normalized": [round(float(e), 6) for e in evals_norm_orig],
            "lambda2": round(lambda2_orig, 6),
            "lambda_n": round(lambda_n_orig, 6),
            "spectral_gap_ratio": round(gap_orig, 6),
            "n_components": count_connected_components(evals_orig),
        },
        "after_removal": {
            "eigenvalues": [round(float(e), 6) for e in evals_red],
            "eigenvalues_normalized": [round(float(e), 6) for e in evals_norm_red],
            "lambda2": round(lambda2_red, 6),
            "lambda_n": round(lambda_n_red, 6),
            "spectral_gap_ratio": round(gap_red, 6),
            "n_components": count_connected_components(evals_red),
        },
        "gap_reduction": round(gap_orig - gap_red, 6),
        "gap_reduction_pct": round(
            (gap_orig - gap_red) / gap_orig * 100 if gap_orig > 1e-12 else 0.0, 2
        ),
    }


# ---------------------------------------------------------------------------
# Analise de perturbacao (robustez)
# ---------------------------------------------------------------------------


def perturbation_analysis(
    adj: np.ndarray,
    node_idx: dict[str, int],
    n_iterations: int = 1000,
    perturbation_fraction: float = 0.20,
    seed: int = 42,
) -> dict[str, Any]:
    """Analise de robustez por perturbacao estocastica dos pesos das arestas.

    Em cada iteracao:
        1. Perturba cada peso de aresta por um fator uniforme em
           [1 - perturbation_fraction, 1 + perturbation_fraction]
        2. Recalcula a analise de impacto de remocao de nos
        3. Identifica qual no tem o maior impacto

    Reporta a fracao de iteracoes em que cada no e o mais critico.

    Args:
        adj: Matriz de adjacencia original.
        node_idx: Mapeamento {id: indice}.
        n_iterations: Numero de iteracoes de perturbacao.
        perturbation_fraction: Fracao maxima de perturbacao (ex: 0.20 = +-20%).
        seed: Seed para reprodutibilidade.

    Returns:
        Dicionario com resultados da analise de perturbacao.
    """
    from aso_math.math_4_spectral.network import remove_node

    rng = np.random.default_rng(seed)
    n = adj.shape[0]

    # Contar quantas vezes cada no e o mais critico
    critical_counts: dict[str, int] = {nid: 0 for nid in node_idx}
    idx_to_id = {v: k for k, v in node_idx.items()}

    for _ in range(n_iterations):
        # Perturbar pesos
        noise = rng.uniform(
            1.0 - perturbation_fraction,
            1.0 + perturbation_fraction,
            size=(n, n),
        )
        # Simetrizar o ruido
        noise = (noise + noise.T) / 2.0

        adj_perturbed = adj * noise
        # Zerar a diagonal (sem self-loops)
        np.fill_diagonal(adj_perturbed, 0.0)
        # Manter zeros onde nao havia aresta
        adj_perturbed[adj == 0] = 0.0

        # Calcular lambda_2 da rede original perturbada
        L_pert, _ = compute_laplacian(adj_perturbed)
        evals_pert, _ = compute_eigenvalues(L_pert)
        lambda2_pert = algebraic_connectivity(evals_pert)

        # Para cada no, calcular impacto da remocao
        max_drop = -np.inf
        most_critical = ""

        for node_id in node_idx:
            adj_red, _ = remove_node(adj_perturbed, node_idx, node_id)
            L_red, _ = compute_laplacian(adj_red)
            evals_red, _ = compute_eigenvalues(L_red)
            lambda2_red = algebraic_connectivity(evals_red)

            drop = lambda2_pert - lambda2_red
            if drop > max_drop:
                max_drop = drop
                most_critical = node_id

        if most_critical:
            critical_counts[most_critical] += 1

    # Calcular fracoes
    results: list[dict[str, Any]] = []
    for node_id, count in critical_counts.items():
        results.append({
            "node_id": node_id,
            "times_most_critical": count,
            "fraction": round(count / n_iterations, 4),
        })

    # Ordenar por frequencia
    results.sort(key=lambda r: r["times_most_critical"], reverse=True)

    return {
        "n_iterations": n_iterations,
        "perturbation_fraction": perturbation_fraction,
        "seed": seed,
        "node_rankings": results,
        "most_robust_critical_node": results[0]["node_id"] if results else "",
        "top_node_fraction": results[0]["fraction"] if results else 0.0,
    }
