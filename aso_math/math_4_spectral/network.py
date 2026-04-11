"""Construcao do grafo do spliceosome de L. infantum.

Modela a rede de interacoes do maquinario de trans-splicing de
trypanosomatideos como um grafo ponderado nao-direcionado.

Notas biologicas:
    - Em trypanosomatideos, 100% dos mRNAs requerem trans-splicing
    - O SL RNA e universalmente necessario — hub central da rede
    - Os snRNAs (U1-U6) formam o core catalitico do spliceosome
    - As proteinas Sm formam o anel heptamerico que estabiliza os snRNPs
    - Prp8 e Snu114 sao componentes catalitcos do tri-snRNP U4/U6.U5

Refs:
    - Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
    - Mair G et al. (2000) RNA 6(2):163-169
    - Michaeli S. (2011) Parasitology 138(12):1-16
    - Will CL, Luhrmann R. (2011) Cold Spring Harb Perspect Biol 3(7):a003707
"""

from __future__ import annotations

from typing import Any

import numpy as np


# ---------------------------------------------------------------------------
# Definicao dos nos do spliceosome
# ---------------------------------------------------------------------------

# Cada no tem: id, nome, tipo funcional, descricao
# Tipos: "sl_rna", "snRNA", "sm_protein", "catalytic_protein", "sl_associated"

SPLICEOSOME_NODES: list[dict[str, str]] = [
    # SL RNA — hub universal do trans-splicing
    {"id": "SL_RNA", "name": "SL RNA", "type": "sl_rna",
     "description": "Spliced Leader RNA — doador universal do exon de 39 nt"},

    # snRNAs do spliceosome (core catalitico)
    {"id": "U1", "name": "U1 snRNA", "type": "snRNA",
     "description": "Reconhece o sitio 5' splice site (cis-splicing)"},
    {"id": "U2", "name": "U2 snRNA", "type": "snRNA",
     "description": "Liga-se ao branch point — essencial para ambos cis e trans"},
    {"id": "U4", "name": "U4 snRNA", "type": "snRNA",
     "description": "Base-pairs com U6, inibe atividade catalitica ate ativacao"},
    {"id": "U5", "name": "U5 snRNA", "type": "snRNA",
     "description": "Posiciona exons para ligacao — loop 1 conservado"},
    {"id": "U6", "name": "U6 snRNA", "type": "snRNA",
     "description": "Catalisa a transesterificacao — ribozima do spliceosome"},

    # Proteinas Sm — anel heptamerico que estabiliza snRNPs
    {"id": "SmB", "name": "SmB", "type": "sm_protein",
     "description": "Componente do anel Sm heptamerico"},
    {"id": "SmD1", "name": "SmD1", "type": "sm_protein",
     "description": "Componente do anel Sm — interacao direta com snRNAs"},
    {"id": "SmD2", "name": "SmD2", "type": "sm_protein",
     "description": "Componente do anel Sm — interacao direta com snRNAs"},
    {"id": "SmD3", "name": "SmD3", "type": "sm_protein",
     "description": "Componente do anel Sm heptamerico"},
    {"id": "SmE", "name": "SmE", "type": "sm_protein",
     "description": "Componente do anel Sm heptamerico"},
    {"id": "SmF", "name": "SmF", "type": "sm_protein",
     "description": "Componente do anel Sm heptamerico"},
    {"id": "SmG", "name": "SmG", "type": "sm_protein",
     "description": "Componente do anel Sm heptamerico"},

    # Proteinas cataliticas do spliceosome
    {"id": "Prp8", "name": "Prp8", "type": "catalytic_protein",
     "description": "Maior proteina do spliceosome — centro catalitico do U5 snRNP"},
    {"id": "Snu114", "name": "Snu114", "type": "catalytic_protein",
     "description": "GTPase do U5 snRNP — regula rearranjos conformacionais"},

    # Proteinas associadas ao SL RNA (especificas de trypanosomatideos)
    {"id": "SLA1", "name": "SLA1", "type": "sl_associated",
     "description": "SL RNA-associated protein 1 — biogenese do SL RNP"},
    {"id": "NHP2L1", "name": "NHP2L1", "type": "sl_associated",
     "description": "Proteina do SL RNP — estabiliza estrutura cap4 do SL RNA"},
]


# ---------------------------------------------------------------------------
# Definicao das arestas (interacoes biologicas)
# ---------------------------------------------------------------------------

# Cada aresta: (no_a, no_b, peso, justificativa)
# Peso reflete a forca/confiabilidade da interacao:
#   1.0 = interacao estrutural direta confirmada
#   0.8 = interacao funcional forte
#   0.6 = co-dependencia funcional moderada
#   0.4 = interacao indireta ou regulatoria
#   0.2 = co-localizacao ou associacao fraca

SPLICEOSOME_EDGES: list[tuple[str, str, float, str]] = [
    # ===================================================================
    # SL RNA — conexoes com TODOS os componentes do spliceosome
    # Justificativa: trans-splicing requer coordenacao do SL RNA com
    # cada componente. Em trypanosomatideos, todo evento de splicing
    # e trans-splicing (nao existe cis-splicing alternativo para mRNAs).
    # ===================================================================

    # SL RNA <-> snRNAs (interacoes funcionais diretas)
    ("SL_RNA", "U2", 1.0,
     "SL RNA e U2 cooperam no reconhecimento do branch point em trans-splicing"),
    ("SL_RNA", "U4", 0.8,
     "SL RNA depende do tri-snRNP U4/U6.U5 para catalisar trans-splicing"),
    ("SL_RNA", "U5", 0.9,
     "U5 loop 1 posiciona o exon do SL para ligacao ao pre-mRNA"),
    ("SL_RNA", "U6", 1.0,
     "U6 catalisa a transesterificacao; SL RNA fornece o exon doador"),
    ("SL_RNA", "U1", 0.6,
     "U1 auxilia no reconhecimento inicial do sitio de splice em trans"),

    # SL RNA <-> proteinas associadas (interacoes estruturais diretas)
    ("SL_RNA", "SLA1", 1.0,
     "SLA1 e componente direto do SL RNP — necessario para biogenese"),
    ("SL_RNA", "NHP2L1", 0.9,
     "NHP2L1 estabiliza a estrutura cap4 unica do SL RNA"),

    # SL RNA <-> proteinas Sm (o SL RNA possui sitio Sm)
    ("SL_RNA", "SmB", 0.8,
     "SL RNA tem sitio de ligacao Sm — interage com anel heptamerico"),
    ("SL_RNA", "SmD1", 0.8,
     "SmD1 liga diretamente ao sitio Sm do SL RNA"),

    # SL RNA <-> proteinas cataliticas
    ("SL_RNA", "Prp8", 0.7,
     "Prp8 no centro catalitico interage com substrato SL RNA"),
    ("SL_RNA", "Snu114", 0.5,
     "Snu114 regula conformacao durante trans-splicing com SL RNA"),

    # ===================================================================
    # Interacoes snRNA-snRNA (core catalitico do spliceosome)
    # ===================================================================

    ("U4", "U6", 1.0,
     "Base-pairing extensivo U4/U6 — di-snRNP estavel"),
    ("U4", "U5", 0.9,
     "Tri-snRNP U4/U6.U5 — complexo pre-catalitico"),
    ("U5", "U6", 0.9,
     "Tri-snRNP U4/U6.U5 — U5 e U6 cooperam na catalise"),
    ("U2", "U6", 0.8,
     "U2/U6 helix I e II — essenciais para centro catalitico"),
    ("U1", "U2", 0.6,
     "U1 recruta U2 ao branch point na montagem do spliceosome"),

    # ===================================================================
    # Interacoes proteinas Sm — anel heptamerico
    # As 7 proteinas Sm formam um anel que se liga aos snRNAs.
    # Modelamos as interacoes adjacentes no anel.
    # ===================================================================

    ("SmB", "SmD3", 0.8, "Adjacentes no anel Sm heptamerico"),
    ("SmD3", "SmD1", 0.8, "Adjacentes no anel Sm heptamerico"),
    ("SmD1", "SmD2", 0.8, "Adjacentes no anel Sm heptamerico"),
    ("SmD2", "SmF", 0.8, "Adjacentes no anel Sm heptamerico"),
    ("SmF", "SmE", 0.8, "Adjacentes no anel Sm heptamerico"),
    ("SmE", "SmG", 0.8, "Adjacentes no anel Sm heptamerico"),
    ("SmG", "SmB", 0.8, "Fecha o anel Sm heptamerico"),

    # ===================================================================
    # Proteinas Sm <-> snRNAs (o anel Sm estabiliza snRNPs)
    # O anel heptamerico completo envolve o sitio Sm de cada snRNA.
    # Diferentes subunidades fazem contatos diretos com diferentes snRNAs.
    # Ref: Kambach C et al. (1999) Cell 96(3):375-387
    # Ref: Pomeranz Krummel DA et al. (2009) Nature 458(7235):475-480
    # ===================================================================

    # U1 snRNP: SmD1 e SmD3 fazem contatos com o sitio Sm do U1
    ("SmD1", "U1", 0.5, "SmD1 contata sitio Sm do U1 snRNA"),
    ("SmD3", "U1", 0.5, "SmD3 contata sitio Sm do U1 snRNA"),

    # U2 snRNP: SmB e SmD1 fazem contatos com o sitio Sm do U2
    ("SmB", "U2", 0.5, "SmB contata sitio Sm do U2 snRNA"),
    ("SmD2", "U2", 0.5, "SmD2 contata sitio Sm do U2 snRNA"),

    # U4 snRNP: SmE e SmG fazem contatos com o sitio Sm do U4
    ("SmE", "U4", 0.5, "SmE contata sitio Sm do U4 snRNA"),
    ("SmG", "U4", 0.5, "SmG contata sitio Sm do U4 snRNA"),

    # U5 snRNP: SmF e SmD2 fazem contatos com o sitio Sm do U5
    ("SmF", "U5", 0.5, "SmF contata sitio Sm do U5 snRNA"),
    ("SmD1", "U5", 0.5, "SmD1 contata sitio Sm do U5 snRNA"),

    # ===================================================================
    # Proteinas cataliticas <-> snRNAs
    # ===================================================================

    ("Prp8", "U5", 1.0,
     "Prp8 e componente integral do U5 snRNP — maior proteina do spliceosome"),
    ("Prp8", "U6", 0.7,
     "Prp8 interage com U6 snRNA no centro catalitico"),
    ("Snu114", "U5", 0.9,
     "Snu114 e GTPase associada ao U5 snRNP"),
    ("Snu114", "Prp8", 0.8,
     "Snu114 e Prp8 cooperam no U5 snRNP"),

    # ===================================================================
    # Proteinas SL-associadas <-> Sm (montagem do SL RNP)
    # ===================================================================

    ("SLA1", "SmD1", 0.5,
     "SLA1 coopera com Sm na montagem do SL RNP"),
    ("NHP2L1", "SmB", 0.4,
     "NHP2L1 interage indiretamente com anel Sm no SL RNP"),
]


# ---------------------------------------------------------------------------
# Construcao da matriz de adjacencia
# ---------------------------------------------------------------------------


def build_node_index(nodes: list[dict[str, str]] | None = None) -> dict[str, int]:
    """Cria mapeamento de id do no para indice numerico.

    Args:
        nodes: Lista de nos. Se None, usa SPLICEOSOME_NODES.

    Returns:
        Dicionario {node_id: indice}.
    """
    if nodes is None:
        nodes = SPLICEOSOME_NODES
    return {node["id"]: i for i, node in enumerate(nodes)}


def build_adjacency_matrix(
    nodes: list[dict[str, str]] | None = None,
    edges: list[tuple[str, str, float, str]] | None = None,
) -> tuple[np.ndarray, dict[str, int], list[dict[str, str]]]:
    """Constroi a matriz de adjacencia ponderada do spliceosome.

    Retorna a matriz simetrica A[i,j] = peso da aresta entre i e j,
    o mapeamento de ids para indices, e a lista de nos utilizada.

    Args:
        nodes: Lista de nos. Se None, usa SPLICEOSOME_NODES.
        edges: Lista de arestas. Se None, usa SPLICEOSOME_EDGES.

    Returns:
        Tupla (adjacency_matrix, node_index, nodes_used).
    """
    if nodes is None:
        nodes = SPLICEOSOME_NODES
    if edges is None:
        edges = SPLICEOSOME_EDGES

    n = len(nodes)
    node_idx = build_node_index(nodes)
    adj = np.zeros((n, n), dtype=np.float64)

    for src, dst, weight, _reason in edges:
        i = node_idx[src]
        j = node_idx[dst]
        adj[i, j] = weight
        adj[j, i] = weight  # Grafo nao-direcionado

    return adj, node_idx, nodes


def remove_node(
    adj: np.ndarray,
    node_idx: dict[str, int],
    node_id: str,
) -> tuple[np.ndarray, dict[str, int]]:
    """Remove um no da matriz de adjacencia.

    Retorna a submatriz (n-1)x(n-1) e o mapeamento atualizado.

    Args:
        adj: Matriz de adjacencia original.
        node_idx: Mapeamento {id: indice} original.
        node_id: Id do no a remover.

    Returns:
        Tupla (nova_adjacency_matrix, novo_node_index).
    """
    idx_to_remove = node_idx[node_id]
    n = adj.shape[0]

    # Indices a manter
    keep = [i for i in range(n) if i != idx_to_remove]
    new_adj = adj[np.ix_(keep, keep)]

    # Reconstruir mapeamento
    id_to_old_idx = {v: k for k, v in node_idx.items()}
    new_node_idx: dict[str, int] = {}
    for new_i, old_i in enumerate(keep):
        old_id = id_to_old_idx[old_i]
        new_node_idx[old_id] = new_i

    return new_adj, new_node_idx


def get_network_stats(
    adj: np.ndarray,
    node_idx: dict[str, int],
) -> dict[str, Any]:
    """Calcula estatisticas basicas da rede.

    Args:
        adj: Matriz de adjacencia.
        node_idx: Mapeamento {id: indice}.

    Returns:
        Dicionario com numero de nos, arestas, densidade, graus.
    """
    n = adj.shape[0]

    # Numero de arestas (dividido por 2 porque a matriz e simetrica)
    n_edges = int(np.count_nonzero(adj) / 2)

    # Densidade = 2E / (N * (N-1))
    max_edges = n * (n - 1) / 2
    density = n_edges / max_edges if max_edges > 0 else 0.0

    # Grau ponderado de cada no
    degree_weighted = np.sum(adj, axis=1)

    # Grau nao-ponderado de cada no
    degree_unweighted = np.count_nonzero(adj, axis=1)

    # Mapa de graus por no
    idx_to_id = {v: k for k, v in node_idx.items()}
    degrees: list[dict[str, Any]] = []
    for i in range(n):
        degrees.append({
            "node_id": idx_to_id[i],
            "degree": int(degree_unweighted[i]),
            "weighted_degree": round(float(degree_weighted[i]), 4),
        })

    # Ordenar por grau ponderado (maior primeiro)
    degrees.sort(key=lambda d: d["weighted_degree"], reverse=True)

    return {
        "n_nodes": n,
        "n_edges": n_edges,
        "density": round(density, 4),
        "max_possible_edges": int(max_edges),
        "node_degrees": degrees,
    }
