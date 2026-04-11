"""Homologia persistente via filtracao de Vietoris-Rips (implementacao from-scratch).

Implementa TDA sem dependencias externas (sem gudhi/ripser), usando apenas numpy.

Algoritmos:
- H0 (componentes conexas): Union-Find com lista de arestas ordenada por distancia.
  A cada epsilon, arestas com peso <= epsilon sao adicionadas. O nascimento de
  cada componente e epsilon=0; a morte e o epsilon em que ela se funde com outra.
  A componente final (ultima a morrer) tem persistencia infinita.

- H1 (loops/ciclos): Deteccao de ciclos via arestas que nao alteram H0.
  Quando uma aresta conecta dois vertices ja na mesma componente, ela fecha
  um ciclo — registramos birth=epsilon_aresta. A morte e estimada pela
  proxima aresta que preenche o triangulo (kill edge).

Referencia:
    Edelsbrunner H, Harer J (2010) "Computational Topology: An Introduction"
    Zomorodian A (2005) "Topology for Computing"
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from core.logger import get_logger

logger = get_logger("math_3_tda.persistence")


# ---------------------------------------------------------------------------
# Union-Find (Disjoint Set Union) para H0
# ---------------------------------------------------------------------------


class UnionFind:
    """Estrutura de dados Union-Find com compressao de caminho e union by rank.

    Usada para rastrear componentes conexas durante a filtracao de Vietoris-Rips.
    Cada merge corresponde a morte de uma componente H0.
    """

    def __init__(self, n: int) -> None:
        """Inicializa n componentes disjuntas.

        Args:
            n: Numero de elementos (vertices do complexo simplicial).
        """
        self.parent = list(range(n))
        self.rank = [0] * n
        self.n_components = n

    def find(self, x: int) -> int:
        """Encontra o representante da componente de x (com compressao de caminho).

        Args:
            x: Indice do elemento.

        Returns:
            Indice do representante da componente.
        """
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x: int, y: int) -> bool:
        """Une as componentes de x e y, se distintas.

        Args:
            x: Indice do primeiro elemento.
            y: Indice do segundo elemento.

        Returns:
            True se a uniao ocorreu (componentes eram distintas), False caso contrario.
        """
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return False

        # Union by rank para manter arvores balanceadas
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

        self.n_components -= 1
        return True


# ---------------------------------------------------------------------------
# Estrutura para pares de persistencia
# ---------------------------------------------------------------------------


@dataclass
class PersistencePair:
    """Um par (birth, death) de um feature topologico.

    Attributes:
        birth: Valor de epsilon onde o feature nasce.
        death: Valor de epsilon onde o feature morre (inf se nunca morre).
        dimension: Dimensao homologica (0=componente, 1=loop, 2=vazio).
        persistence: Tempo de vida = death - birth.
    """
    birth: float
    death: float
    dimension: int

    @property
    def persistence(self) -> float:
        """Calcula a persistencia (tempo de vida) do feature."""
        if self.death == float("inf"):
            return float("inf")
        return self.death - self.birth

    def to_dict(self) -> dict[str, Any]:
        """Serializa para dicionario JSON-compativel."""
        return {
            "birth": round(self.birth, 4),
            "death": round(self.death, 4) if self.death != float("inf") else "inf",
            "persistence": round(self.persistence, 4) if self.persistence != float("inf") else "inf",
            "dimension": self.dimension,
        }


@dataclass
class PersistenceDiagram:
    """Diagrama de persistencia completo para uma estrutura.

    Contem todos os pares (birth, death) organizados por dimensao homologica.
    """
    name: str
    h0_pairs: list[PersistencePair] = field(default_factory=list)
    h1_pairs: list[PersistencePair] = field(default_factory=list)

    @property
    def total_persistence_h0(self) -> float:
        """Soma das persistencias finitas de H0."""
        return sum(
            p.persistence for p in self.h0_pairs
            if p.persistence != float("inf")
        )

    @property
    def total_persistence_h1(self) -> float:
        """Soma das persistencias de H1."""
        return sum(
            p.persistence for p in self.h1_pairs
            if p.persistence != float("inf")
        )

    def n_persistent_h1(self, threshold: float = 0.1) -> int:
        """Conta features H1 com persistencia acima de um limiar.

        O limiar padrao de 0.1 A e adequado para nuvens de pontos de backbone
        nucleico, onde ciclos genuinos (sulcos, entrecruzamento de fitas)
        persistem por 0.1-0.5 A na filtracao de Vietoris-Rips.

        Args:
            threshold: Limiar de persistencia minima (Angstroms).
        """
        return sum(
            1 for p in self.h1_pairs
            if p.persistence != float("inf") and p.persistence > threshold
        )

    def betti_at_epsilon(self, epsilon: float) -> dict[str, int]:
        """Calcula numeros de Betti em um dado epsilon.

        beta_0 = numero de componentes conexas vivas em epsilon
        beta_1 = numero de loops vivos em epsilon

        Args:
            epsilon: Valor do raio de filtracao.

        Returns:
            Dicionario com beta_0 e beta_1.
        """
        beta_0 = sum(
            1 for p in self.h0_pairs
            if p.birth <= epsilon and (p.death > epsilon or p.death == float("inf"))
        )
        beta_1 = sum(
            1 for p in self.h1_pairs
            if p.birth <= epsilon and (p.death > epsilon or p.death == float("inf"))
        )
        return {"beta_0": beta_0, "beta_1": beta_1}

    def to_dict(self) -> dict[str, Any]:
        """Serializa o diagrama completo para JSON."""
        return {
            "name": self.name,
            "h0_pairs": [p.to_dict() for p in self.h0_pairs],
            "h1_pairs": [p.to_dict() for p in self.h1_pairs],
            "total_persistence_h0": round(self.total_persistence_h0, 4),
            "total_persistence_h1": round(self.total_persistence_h1, 4),
            "n_persistent_h1_features": self.n_persistent_h1(),
        }


# ---------------------------------------------------------------------------
# Matriz de distancia
# ---------------------------------------------------------------------------


def compute_distance_matrix(points: np.ndarray) -> np.ndarray:
    """Calcula a matriz de distancias euclidianas entre todos os pares de pontos.

    Usa a expansao ||a - b||^2 = ||a||^2 + ||b||^2 - 2*a.b para eficiencia.

    Args:
        points: Array Nx3 de coordenadas.

    Returns:
        Matriz NxN de distancias simetricas.
    """
    # Produto escalar ponto a ponto
    sq_norms = np.sum(points ** 2, axis=1)
    # ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a.b
    dist_sq = sq_norms[:, np.newaxis] + sq_norms[np.newaxis, :] - 2.0 * (points @ points.T)
    # Corrigir possiveis negativos por erro numerico
    np.maximum(dist_sq, 0.0, out=dist_sq)
    return np.sqrt(dist_sq)


# ---------------------------------------------------------------------------
# Filtracao de Vietoris-Rips — H0 (Union-Find)
# ---------------------------------------------------------------------------


def compute_h0_persistence(dist_matrix: np.ndarray) -> list[PersistencePair]:
    """Calcula persistencia H0 (componentes conexas) via Union-Find.

    Algoritmo:
    1. Gera lista de todas as arestas (i, j, dist) com i < j
    2. Ordena por distancia crescente
    3. Percorre as arestas: se i e j estao em componentes diferentes,
       faz union — a componente menor "morre" em epsilon = dist
    4. A ultima componente (que nunca morre) tem death = inf

    Args:
        dist_matrix: Matriz NxN de distancias.

    Returns:
        Lista de PersistencePair para H0.
    """
    n = dist_matrix.shape[0]
    if n == 0:
        return []

    # Gerar lista de arestas (triangulo superior da matriz)
    edges: list[tuple[float, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            edges.append((dist_matrix[i, j], i, j))

    # Ordenar por distancia
    edges.sort(key=lambda e: e[0])

    uf = UnionFind(n)

    # Cada vertice nasce em epsilon=0. Rastrear morte de cada componente.
    # death_times[representante] = epsilon em que a componente foi absorvida
    death_times: dict[int, float] = {}

    for dist, i, j in edges:
        ri, rj = uf.find(i), uf.find(j)
        if ri != rj:
            # A componente menor morre — registramos a morte
            # Convencao: a componente com representante "mais novo" morre
            # Na pratica, o que importa e que uma das duas morre em epsilon=dist
            merged = uf.union(i, j)
            if merged:
                # Determinar qual representante sobreviveu
                new_root = uf.find(i)
                dead_root = rj if new_root == ri else ri
                death_times[dead_root] = dist

    # Montar pares de persistencia
    pairs: list[PersistencePair] = []

    # Todos os vertices originais nascem em epsilon=0
    # Componentes que morreram tem death registrado
    # A componente final tem death = inf
    all_original_roots = set(range(n))
    dead_roots = set(death_times.keys())
    alive_roots = all_original_roots - dead_roots

    for root in dead_roots:
        pairs.append(PersistencePair(
            birth=0.0,
            death=death_times[root],
            dimension=0,
        ))

    # Componentes que sobreviveram (tipicamente 1 se o grafo e conexo)
    for root in alive_roots:
        # Verificar se este root ainda e raiz no union-find final
        if uf.find(root) == root:
            pairs.append(PersistencePair(
                birth=0.0,
                death=float("inf"),
                dimension=0,
            ))

    # Ordenar por persistencia (decrescente) para facilitar analise
    pairs.sort(key=lambda p: p.persistence if p.persistence != float("inf") else 1e18, reverse=True)

    logger.info("H0: %d pares de persistencia computados", len(pairs))
    return pairs


# ---------------------------------------------------------------------------
# Filtracao de Vietoris-Rips — H1 (ciclos simplificados)
# ---------------------------------------------------------------------------


def compute_h1_persistence(
    dist_matrix: np.ndarray,
    epsilon_max: float = 20.0,
    epsilon_step: float = 0.5,
) -> list[PersistencePair]:
    """Calcula persistencia H1 (loops/ciclos) via deteccao de arestas ciclo.

    Abordagem simplificada:
    1. Percorre arestas em ordem crescente de distancia (como em H0)
    2. Uma aresta que conecta dois vertices JA na mesma componente fecha um ciclo
       -> birth do ciclo = distancia dessa aresta
    3. Para estimar a death: procura a menor aresta que forma um triangulo
       com a aresta do ciclo (preenchendo o buraco)

    Esta e uma aproximacao da homologia persistente real (que requer
    reducao de matriz de bordas), mas captura os ciclos mais significativos.

    Args:
        dist_matrix: Matriz NxN de distancias.
        epsilon_max: Raio maximo de filtracao.
        epsilon_step: Passo de discretizacao (nao usado diretamente, mas define
                      a resolucao de epsilon no contexto da filtracao).

    Returns:
        Lista de PersistencePair para H1.
    """
    n = dist_matrix.shape[0]
    if n < 3:
        return []

    # Gerar arestas ordenadas
    edges: list[tuple[float, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            d = dist_matrix[i, j]
            if d <= epsilon_max:
                edges.append((d, i, j))

    edges.sort(key=lambda e: e[0])

    uf = UnionFind(n)
    h1_pairs: list[PersistencePair] = []

    # Conjunto de arestas adicionadas (para busca de triangulos)
    added_edges: set[tuple[int, int]] = set()

    # Para cada aresta que fecha um ciclo, registrar e estimar death
    cycle_edges: list[tuple[float, int, int]] = []

    for dist, i, j in edges:
        if uf.find(i) == uf.find(j):
            # Aresta fecha um ciclo — H1 nasce aqui
            cycle_edges.append((dist, i, j))
        else:
            uf.union(i, j)

        added_edges.add((min(i, j), max(i, j)))

    # Estimar death de cada ciclo
    # O ciclo morre quando um 2-simplex (triangulo) preenche o buraco.
    # Para a aresta (i, j) que criou o ciclo em epsilon=birth, procurar
    # o vertice k que forma o triangulo com menor epsilon de aparicao.
    #
    # O triangulo (i, j, k) aparece em epsilon = max(dist(i,j), dist(i,k), dist(j,k)).
    # Para que o ciclo tenha persistencia nao-trivial, precisamos que o
    # triangulo apareca DEPOIS do ciclo: death > birth (estritamente).
    #
    # Se todos os triangulos aparecem exatamente em birth, o ciclo e
    # "instantaneo" — nao e um feature topologico significativo.
    # Nesse caso, registramos com persistencia zero (sera filtrado depois).

    for birth, i, j in cycle_edges:
        min_death = float("inf")

        for k in range(n):
            if k == i or k == j:
                continue

            d_ik = dist_matrix[i, k]
            d_jk = dist_matrix[j, k]

            # O triangulo (i, j, k) aparece quando todas as tres arestas existem
            # O epsilon necessario e o maximo das tres distancias
            triangle_epsilon = max(birth, d_ik, d_jk)

            if triangle_epsilon <= epsilon_max and triangle_epsilon < min_death:
                min_death = triangle_epsilon

        if min_death <= epsilon_max:
            h1_pairs.append(PersistencePair(
                birth=birth,
                death=min_death,
                dimension=1,
            ))
        else:
            # Ciclo que nunca morre dentro do range — persistencia truncada
            h1_pairs.append(PersistencePair(
                birth=birth,
                death=epsilon_max,
                dimension=1,
            ))

    # Filtrar ciclos triviais (persistencia = 0) — nao sao features significativos
    h1_pairs = [p for p in h1_pairs if p.persistence > 1e-10]

    # Ordenar por persistencia decrescente
    h1_pairs.sort(
        key=lambda p: p.persistence if p.persistence != float("inf") else 1e18,
        reverse=True,
    )

    # Limitar o numero de ciclos para evitar explosao computacional
    # Mantemos os 50 mais persistentes (os mais significativos topologicamente)
    max_cycles = 50
    if len(h1_pairs) > max_cycles:
        logger.info(
            "H1: truncando de %d para %d ciclos (mantendo os mais persistentes)",
            len(h1_pairs), max_cycles,
        )
        h1_pairs = h1_pairs[:max_cycles]

    logger.info("H1: %d pares de persistencia computados", len(h1_pairs))
    return h1_pairs


# ---------------------------------------------------------------------------
# Filtracao completa: combina H0 + H1
# ---------------------------------------------------------------------------


def compute_persistence_diagram(
    points: np.ndarray,
    name: str,
    epsilon_max: float = 20.0,
    epsilon_step: float = 0.5,
) -> PersistenceDiagram:
    """Calcula o diagrama de persistencia completo para uma nuvem de pontos.

    Executa:
    1. Amostragem de pontos se a nuvem e muito grande (> 200 pontos)
    2. Matriz de distancias
    3. Filtracao H0 (union-find)
    4. Filtracao H1 (ciclos simplificados)

    Args:
        points: Array Nx3 de coordenadas.
        name: Identificador da estrutura (para logging e JSON).
        epsilon_max: Raio maximo de filtracao em Angstroms.
        epsilon_step: Passo de discretizacao.

    Returns:
        PersistenceDiagram com todos os pares H0 e H1.
    """
    # Subamostrar se necessario para manter computacao tratavel
    # H1 tem complexidade O(n^3) no pior caso
    max_points = 200
    if points.shape[0] > max_points:
        logger.info(
            "%s: subsampling %d -> %d pontos para filtracao",
            name, points.shape[0], max_points,
        )
        indices = np.random.default_rng(seed=42).choice(
            points.shape[0], size=max_points, replace=False,
        )
        points = points[indices]

    logger.info("%s: calculando matriz de distancias (%d pontos)...", name, len(points))
    dist_matrix = compute_distance_matrix(points)

    logger.info("%s: computando H0 (componentes conexas)...", name)
    h0 = compute_h0_persistence(dist_matrix)

    logger.info("%s: computando H1 (loops/ciclos)...", name)
    h1 = compute_h1_persistence(dist_matrix, epsilon_max=epsilon_max, epsilon_step=epsilon_step)

    diagram = PersistenceDiagram(name=name, h0_pairs=h0, h1_pairs=h1)

    logger.info(
        "%s: diagrama completo — H0=%d pares, H1=%d pares, "
        "persistencia total H0=%.2f, H1=%.2f",
        name, len(h0), len(h1),
        diagram.total_persistence_h0, diagram.total_persistence_h1,
    )

    return diagram


# ---------------------------------------------------------------------------
# Distancia bottleneck entre diagramas de persistencia
# ---------------------------------------------------------------------------


def bottleneck_distance(
    diagram_a: PersistenceDiagram,
    diagram_b: PersistenceDiagram,
    dimension: int = 0,
) -> float:
    """Calcula a distancia bottleneck entre dois diagramas de persistencia.

    A distancia bottleneck e o infimo sobre todos os matchings perfeitos
    do supremo das distancias L_inf entre pares matched.

    Implementacao simplificada:
    1. Extrai pares finitos de cada diagrama na dimensao solicitada
    2. Projeta cada par nao-matched na diagonal (custo = persistencia / 2)
    3. Usa matching guloso (greedy) como aproximacao
       (O matching otimo requer algoritmo hungaro, mas greedy e suficiente
       para nuvens de pontos de tamanho moderado)

    Args:
        diagram_a: Primeiro diagrama.
        diagram_b: Segundo diagrama.
        dimension: Dimensao homologica (0 ou 1).

    Returns:
        Distancia bottleneck (L_inf entre diagramas).
    """
    # Extrair pares finitos na dimensao solicitada
    if dimension == 0:
        pairs_a = [(p.birth, p.death) for p in diagram_a.h0_pairs if p.death != float("inf")]
        pairs_b = [(p.birth, p.death) for p in diagram_b.h0_pairs if p.death != float("inf")]
    else:
        pairs_a = [(p.birth, p.death) for p in diagram_a.h1_pairs if p.death != float("inf")]
        pairs_b = [(p.birth, p.death) for p in diagram_b.h1_pairs if p.death != float("inf")]

    if not pairs_a and not pairs_b:
        return 0.0

    # Custo de projetar um ponto na diagonal (L_inf)
    def diagonal_cost(b: float, d: float) -> float:
        return (d - b) / 2.0

    # Custo L_inf entre dois pares
    def pair_cost(p1: tuple[float, float], p2: tuple[float, float]) -> float:
        return max(abs(p1[0] - p2[0]), abs(p1[1] - p2[1]))

    # Equalizar tamanhos adicionando projecoes na diagonal
    # Cada ponto nao-matched e projetado na diagonal com custo = persistencia/2
    n_a, n_b = len(pairs_a), len(pairs_b)

    # Se um diagrama e vazio, o bottleneck e a maxima projecao diagonal do outro
    if n_a == 0:
        return max(diagonal_cost(b, d) for b, d in pairs_b)
    if n_b == 0:
        return max(diagonal_cost(b, d) for b, d in pairs_a)

    # Matching greedy: para cada par em A, encontrar o par mais proximo em B
    # e vice-versa (matching bidirecional guloso)
    used_b: set[int] = set()
    max_cost = 0.0

    # Ordenar A por persistencia decrescente (priorizar features importantes)
    sorted_a = sorted(range(n_a), key=lambda i: pairs_a[i][1] - pairs_a[i][0], reverse=True)

    for idx_a in sorted_a:
        pa = pairs_a[idx_a]
        diag_cost_a = diagonal_cost(pa[0], pa[1])

        best_cost = diag_cost_a
        best_b: int | None = None

        for idx_b in range(n_b):
            if idx_b in used_b:
                continue
            cost = pair_cost(pa, pairs_b[idx_b])
            if cost < best_cost:
                best_cost = cost
                best_b = idx_b

        if best_b is not None:
            used_b.add(best_b)

        max_cost = max(max_cost, best_cost)

    # Pontos de B nao matched -> custo diagonal
    for idx_b in range(n_b):
        if idx_b not in used_b:
            pb = pairs_b[idx_b]
            cost = diagonal_cost(pb[0], pb[1])
            max_cost = max(max_cost, cost)

    return max_cost


# ---------------------------------------------------------------------------
# Numeros de Betti em escalas caracteristicas
# ---------------------------------------------------------------------------


def betti_numbers_at_scales(
    diagram: PersistenceDiagram,
    scales: list[float] | None = None,
) -> list[dict[str, Any]]:
    """Calcula numeros de Betti em escalas caracteristicas de epsilon.

    Escalas padrao: [2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0] Angstroms,
    correspondendo a:
    - 2-4 A: ligacoes covalentes e pontes de hidrogenio
    - 6-8 A: empilhamento de bases (stacking)
    - 10-15 A: diametro da helice B-DNA
    - 20 A: escala do sulco maior

    Args:
        diagram: Diagrama de persistencia.
        scales: Lista de valores epsilon (A). Se None, usa padrao.

    Returns:
        Lista de dicionarios com epsilon, beta_0, beta_1.
    """
    if scales is None:
        scales = [2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0]

    results: list[dict[str, Any]] = []
    for eps in scales:
        betti = diagram.betti_at_epsilon(eps)
        results.append({
            "epsilon_angstrom": eps,
            "beta_0": betti["beta_0"],
            "beta_1": betti["beta_1"],
        })

    return results
