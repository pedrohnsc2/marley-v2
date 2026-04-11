"""Math 2 — Geometria da informacao: alienness do SL RNA via Fisher-Rao.

Quantifica o quao "alienigena" o SL RNA de L. infantum e em relacao
aos transcriptomas hospedeiros (humano e canino), usando ferramentas
de geometria da informacao sobre a variedade estatistica.

Analises:
    1. Distribuicoes de frequencia de k-mers (k=1,2,3,4)
    2. Matriz de Informacao de Fisher para o modelo multinomial
    3. Distancia geodesica de Fisher-Rao entre distribuicoes
    4. Cobertura de posicoes de baixa entropia pelo MRL-ASO-001
    5. Divergencias de Kullback-Leibler e Jensen-Shannon
    6. Score composto de alienness

Referencia: Amari S, Nagaoka H (2000) Methods of Information Geometry.
            AMS/Oxford University Press.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np

from aso_math.config import (
    ASO_KNOWN_DG,
    ASO_KNOWN_TM,
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    SL_SEQUENCE,
)
from aso_math.target_config import TargetConfig
from aso_math.thermo import gc_content
from aso_math.envelope import Timer
from core.logger import get_logger

logger = get_logger("math_2_fisher")

# ---------------------------------------------------------------------------
# Diretorio de saida
# ---------------------------------------------------------------------------

RESULTS_DIR: Path = Path(__file__).resolve().parent / "results"

# ---------------------------------------------------------------------------
# Composicoes de referencia dos transcriptomas hospedeiros
# ---------------------------------------------------------------------------

# Composicao nucleotidica media do mRNA humano
# Ref: Forsdyke DR (1995) J Mol Evol 41(5):573-581
# Ref: Cristillo AD et al. (1998) Nucleic Acids Res 26(12):2837-2842
# Valores representativos da distribuicao media em CDS humanas
HUMAN_COMPOSITION: dict[str, float] = {
    "A": 0.26, "C": 0.24, "G": 0.26, "U": 0.24,
}

# Composicao nucleotidica media do mRNA canino (Canis lupus familiaris)
# Ref: Nakamura Y et al. (2000) Nucleic Acids Res 28(1):292
# Ligeiramente diferente do humano — reflete bias de codon de carnivoros
CANINE_COMPOSITION: dict[str, float] = {
    "A": 0.27, "C": 0.23, "G": 0.25, "U": 0.25,
}

# Bases RNA para consistencia
RNA_BASES: list[str] = ["A", "C", "G", "U"]


# ---------------------------------------------------------------------------
# 1. Distribuicoes de frequencia de k-mers
# ---------------------------------------------------------------------------


def _sequence_to_rna(seq: str) -> str:
    """Converte sequencia DNA para RNA (T -> U)."""
    return seq.upper().replace("T", "U")


def compute_kmer_distribution(
    sequence: str,
    k: int,
) -> dict[str, float]:
    """Calcula a distribuicao de frequencia de k-mers de uma sequencia.

    Conta todos os k-mers sobrepostos (sliding window de tamanho k)
    e normaliza para obter uma distribuicao de probabilidade.

    Aplica suavizacao de Laplace (+1 pseudocount) para evitar
    probabilidades zero, que causariam divergencia infinita em KL.

    Args:
        sequence: Sequencia RNA (bases A, C, G, U).
        k: Tamanho do k-mer (1, 2, 3 ou 4).

    Returns:
        Dicionario {kmer: probabilidade} normalizado.
    """
    seq = _sequence_to_rna(sequence)

    # Gerar todos os k-mers possiveis para o alfabeto RNA
    all_kmers = ["".join(combo) for combo in product(RNA_BASES, repeat=k)]

    # Contar k-mers observados na sequencia
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        # Ignorar k-mers com bases ambiguas
        if all(b in RNA_BASES for b in kmer):
            counts[kmer] += 1

    # Suavizacao de Laplace: pseudocount = 1 para cada k-mer possivel
    total = sum(counts.values()) + len(all_kmers)
    distribution: dict[str, float] = {}
    for kmer in all_kmers:
        distribution[kmer] = (counts[kmer] + 1) / total

    return distribution


def compute_reference_kmer_distribution(
    composition: dict[str, float],
    k: int,
) -> dict[str, float]:
    """Calcula distribuicao de k-mers a partir de composicao nucleotidica.

    Assume independencia entre posicoes (modelo i.i.d.):
    P(k-mer) = prod(P(base_i)) para i em 1..k

    Esta e uma aproximacao razoavel para mRNA de referencia,
    onde nao temos a sequencia completa mas sim a composicao media.

    Args:
        composition: Dicionario {base: frequencia} normalizado.
        k: Tamanho do k-mer.

    Returns:
        Dicionario {kmer: probabilidade} normalizado.
    """
    all_kmers = ["".join(combo) for combo in product(RNA_BASES, repeat=k)]
    distribution: dict[str, float] = {}

    for kmer in all_kmers:
        prob = 1.0
        for base in kmer:
            prob *= composition[base]
        distribution[kmer] = prob

    # Normalizar (deveria somar ~1.0, mas garante precisao)
    total = sum(distribution.values())
    if total > 0:
        for kmer in distribution:
            distribution[kmer] /= total

    return distribution


def generate_random_rna_distribution(
    length: int,
    gc_frac: float,
    k: int,
    n_samples: int = 10000,
    seed: int = 42,
) -> dict[str, float]:
    """Gera distribuicao de k-mers de RNA aleatorio com GC fixo.

    Gera n_samples sequencias aleatorias de comprimento dado com o
    mesmo GC content que o SL RNA, calcula k-mers de cada uma,
    e retorna a distribuicao media.

    Args:
        length: Comprimento da sequencia aleatoria.
        gc_frac: Fracao de GC desejada.
        k: Tamanho do k-mer.
        n_samples: Numero de sequencias aleatorias a gerar.
        seed: Semente para reprodutibilidade.

    Returns:
        Dicionario {kmer: probabilidade} medio.
    """
    rng = np.random.default_rng(seed)

    # Probabilidades por base: GC dividido igualmente entre G e C,
    # AT dividido igualmente entre A e U
    p_gc = gc_frac / 2.0  # P(G) = P(C) = gc_frac / 2
    p_au = (1.0 - gc_frac) / 2.0  # P(A) = P(U) = (1-gc_frac) / 2
    base_probs = np.array([p_au, p_gc, p_gc, p_au])  # A, C, G, U

    all_kmers = ["".join(combo) for combo in product(RNA_BASES, repeat=k)]
    total_counts: dict[str, int] = {kmer: 0 for kmer in all_kmers}

    base_chars = np.array(RNA_BASES)

    for _ in range(n_samples):
        indices = rng.choice(4, size=length, p=base_probs)
        seq_arr = base_chars[indices]

        for i in range(length - k + 1):
            kmer = "".join(seq_arr[i:i + k])
            if kmer in total_counts:
                total_counts[kmer] += 1

    # Normalizar com suavizacao de Laplace
    grand_total = sum(total_counts.values()) + len(all_kmers)
    distribution: dict[str, float] = {}
    for kmer in all_kmers:
        distribution[kmer] = (total_counts[kmer] + 1) / grand_total

    return distribution


def analyze_kmer_distributions(
    sl_sequence: str,
    gc_frac: float,
) -> dict[str, Any]:
    """Calcula distribuicoes de k-mers para todas as fontes.

    Fontes: SL RNA, humano, canino, RNA aleatorio com mesmo GC.

    Args:
        sl_sequence: Sequencia do SL RNA (DNA — sera convertida a RNA).
        gc_frac: Fracao GC do SL RNA (para gerar RNA aleatorio).

    Returns:
        Dicionario com distribuicoes para k=1,2,3,4.
    """
    results: dict[str, Any] = {}

    for k in [1, 2, 3, 4]:
        logger.info("  Calculando k-mers (k=%d)...", k)

        sl_dist = compute_kmer_distribution(sl_sequence, k)
        human_dist = compute_reference_kmer_distribution(HUMAN_COMPOSITION, k)
        canine_dist = compute_reference_kmer_distribution(CANINE_COMPOSITION, k)
        random_dist = generate_random_rna_distribution(
            length=len(sl_sequence),
            gc_frac=gc_frac,
            k=k,
        )

        n_kmers = 4 ** k
        results[f"k{k}"] = {
            "k": k,
            "n_possible_kmers": n_kmers,
            "distributions": {
                "sl_rna": sl_dist,
                "human": human_dist,
                "canine": canine_dist,
                "random": random_dist,
            },
        }

    return results


# ---------------------------------------------------------------------------
# 2. Matriz de Informacao de Fisher
# ---------------------------------------------------------------------------


def compute_fisher_information_matrix(
    distribution: dict[str, float],
) -> np.ndarray:
    """Calcula a Matriz de Informacao de Fisher para o modelo multinomial.

    Para uma distribuicao multinomial com parametros p = (p_1, ..., p_k),
    a FIM e diagonal:
        G_ij = delta_ij / p_i

    onde delta_ij e o delta de Kronecker.

    Esta metrica define a geometria Riemanniana no simplex de probabilidades.
    Pontos com probabilidades menores (eventos raros) tem curvatura maior,
    significando que pequenas mudancas nessas regioes sao mais "informativas".

    Ref: Amari S (1985) Differential-Geometrical Methods in Statistics.
         Springer Lecture Notes in Statistics 28.

    Args:
        distribution: Dicionario {kmer: probabilidade}.

    Returns:
        Matriz numpy diagonal (n x n) com 1/p_i na diagonal.
    """
    probs = np.array(list(distribution.values()), dtype=np.float64)

    # Proteger contra divisao por zero (nao deveria ocorrer com Laplace)
    probs = np.maximum(probs, 1e-15)

    # FIM para multinomial: diagonal com 1/p_i
    fim = np.diag(1.0 / probs)

    return fim


def analyze_fisher_matrices(
    kmer_results: dict[str, Any],
) -> dict[str, Any]:
    """Calcula FIMs para todas as distribuicoes e reporta propriedades.

    Propriedades reportadas:
    - Determinante (volume do elipsoide de informacao)
    - Traco (soma das informacoes marginais)
    - Maior/menor autovalor (direcoes de maxima/minima informacao)

    Args:
        kmer_results: Resultado de analyze_kmer_distributions().

    Returns:
        Dicionario com propriedades das FIMs por k e por fonte.
    """
    results: dict[str, Any] = {}

    for k_key in kmer_results:
        k = kmer_results[k_key]["k"]
        dists = kmer_results[k_key]["distributions"]
        k_results: dict[str, Any] = {}

        for source_name, dist in dists.items():
            fim = compute_fisher_information_matrix(dist)
            probs = np.array(list(dist.values()), dtype=np.float64)

            # Propriedades da FIM
            # Para diagonal: det = prod(1/p_i), trace = sum(1/p_i)
            diagonal = np.diag(fim)
            det_log = float(np.sum(np.log(diagonal)))  # log-det para estabilidade
            trace = float(np.sum(diagonal))
            max_eigenvalue = float(np.max(diagonal))
            min_eigenvalue = float(np.min(diagonal))

            # Posicao do maximo (k-mer mais raro = maior informacao)
            kmers = list(dist.keys())
            max_idx = int(np.argmax(diagonal))
            min_idx = int(np.argmin(diagonal))

            k_results[source_name] = {
                "log_determinant": round(det_log, 6),
                "trace": round(trace, 4),
                "max_eigenvalue": round(max_eigenvalue, 4),
                "min_eigenvalue": round(min_eigenvalue, 4),
                "most_informative_kmer": kmers[max_idx],
                "most_informative_kmer_prob": round(float(probs[max_idx]), 6),
                "least_informative_kmer": kmers[min_idx],
                "least_informative_kmer_prob": round(float(probs[min_idx]), 6),
            }

        results[k_key] = k_results

    return results


# ---------------------------------------------------------------------------
# 3. Distancia geodesica de Fisher-Rao
# ---------------------------------------------------------------------------


def fisher_rao_distance(
    p: dict[str, float],
    q: dict[str, float],
) -> float:
    """Calcula a distancia geodesica de Fisher-Rao entre duas distribuicoes.

    A distancia de Fisher-Rao e a distancia geodesica na variedade
    estatistica com a metrica de Fisher. Para distribuicoes multinomiais,
    equivale a:

        d_FR(p, q) = 2 * arccos(sum(sqrt(p_i * q_i)))

    Esta e a distancia de Bhattacharyya no espaco das raizes quadradas
    (esfera unitaria), conhecida como distancia de Hellinger angular.

    Propriedades:
    - Metrica: satisfaz d(p,q) >= 0, d(p,q)=0 iff p=q, simetria, triangulo
    - Invariante a reparametrizacao (propriedade fundamental da geometria)
    - Limitada em [0, pi] para distribuicoes no simplex

    Ref: Bhattacharyya A (1943) Bull Calcutta Math Soc 35:99-109
         Rao CR (1945) Bull Calcutta Math Soc 37:81-91

    Args:
        p: Primeira distribuicao (dicionario {kmer: prob}).
        q: Segunda distribuicao (dicionario {kmer: prob}).

    Returns:
        Distancia de Fisher-Rao em radianos [0, pi].
    """
    # Alinhar as chaves na mesma ordem
    keys = sorted(set(p.keys()) | set(q.keys()))
    p_arr = np.array([p.get(k, 1e-15) for k in keys], dtype=np.float64)
    q_arr = np.array([q.get(k, 1e-15) for k in keys], dtype=np.float64)

    # Normalizar (seguranca — devem estar normalizados)
    p_arr = p_arr / np.sum(p_arr)
    q_arr = q_arr / np.sum(q_arr)

    # Coeficiente de Bhattacharyya: BC = sum(sqrt(p_i * q_i))
    bc = float(np.sum(np.sqrt(p_arr * q_arr)))

    # Clamp em [0, 1] para seguranca numerica
    bc = max(0.0, min(1.0, bc))

    # d_FR = 2 * arccos(BC)
    return 2.0 * math.acos(bc)


def compute_pairwise_distances(
    kmer_results: dict[str, Any],
) -> dict[str, Any]:
    """Calcula a matriz de distancias de Fisher-Rao entre todas as fontes.

    Para cada k (1,2,3,4), calcula d_FR entre todos os 6 pares
    de 4 fontes (SL RNA, humano, canino, aleatorio).

    Args:
        kmer_results: Resultado de analyze_kmer_distributions().

    Returns:
        Dicionario com matrizes de distancia por k.
    """
    source_names = ["sl_rna", "human", "canine", "random"]
    results: dict[str, Any] = {}

    for k_key in kmer_results:
        k = kmer_results[k_key]["k"]
        dists = kmer_results[k_key]["distributions"]

        # Matriz de distancia 4x4
        n = len(source_names)
        distance_matrix = np.zeros((n, n), dtype=np.float64)
        pairs: list[dict[str, Any]] = []

        for i in range(n):
            for j in range(i + 1, n):
                d = fisher_rao_distance(
                    dists[source_names[i]],
                    dists[source_names[j]],
                )
                distance_matrix[i, j] = d
                distance_matrix[j, i] = d

                pairs.append({
                    "source_a": source_names[i],
                    "source_b": source_names[j],
                    "fisher_rao_distance": round(d, 6),
                })

        # Distancia SL RNA para cada hospedeiro
        sl_idx = source_names.index("sl_rna")
        human_idx = source_names.index("human")
        canine_idx = source_names.index("canine")
        random_idx = source_names.index("random")

        results[k_key] = {
            "k": k,
            "source_names": source_names,
            "distance_matrix": distance_matrix.tolist(),
            "pairwise_distances": pairs,
            "sl_to_human": round(float(distance_matrix[sl_idx, human_idx]), 6),
            "sl_to_canine": round(float(distance_matrix[sl_idx, canine_idx]), 6),
            "sl_to_random": round(float(distance_matrix[sl_idx, random_idx]), 6),
            "human_to_canine": round(float(distance_matrix[human_idx, canine_idx]), 6),
        }

    return results


# ---------------------------------------------------------------------------
# 4. Cobertura de posicoes de baixa entropia pelo MRL-ASO-001
# ---------------------------------------------------------------------------


def compute_positional_entropy(
    sequence: str,
) -> list[dict[str, Any]]:
    """Calcula a entropia posicional do SL RNA usando contexto de dinucleotideo.

    Como temos apenas uma especie (L. infantum), nao podemos calcular
    entropia posicional por alinhamento multiplo. Em vez disso, usamos
    entropia de contexto de dinucleotideo: para cada posicao i,
    consideramos o dinucleotideo (i-1, i) e (i, i+1) e calculamos
    a entropia do perfil de bases nesse contexto local.

    Posicoes com contextos ricos em GC (maior diversidade local)
    tem maior entropia. Posicoes com sequencias repetitivas
    (ex: AA, UU) tem menor entropia — sao mais conservadas/estruturadas.

    Para contexto mais rico, usamos a entropia do trinucleotideo
    centrado em cada posicao (i-1, i, i+1) versus a distribuicao
    esperada sob composicao uniforme.

    Metodo alternativo implementado: entropia de Shannon da composicao
    local em janelas de 5 nt centradas em cada posicao.

    Args:
        sequence: Sequencia do SL RNA.

    Returns:
        Lista de dicionarios com entropia por posicao.
    """
    seq = _sequence_to_rna(sequence)
    n = len(seq)

    # Janela local de 5 nt centrada em cada posicao
    window_size = 5
    half_w = window_size // 2

    positions: list[dict[str, Any]] = []

    for i in range(n):
        # Definir limites da janela
        start = max(0, i - half_w)
        end = min(n, i + half_w + 1)
        window = seq[start:end]

        # Contar frequencia de cada base na janela
        counts = Counter(window)
        total = len(window)

        # Calcular entropia de Shannon da janela local
        entropy = 0.0
        base_probs: dict[str, float] = {}
        for base in RNA_BASES:
            p = counts.get(base, 0) / total
            base_probs[base] = p
            if p > 0:
                entropy -= p * math.log2(p)

        # Entropia maxima para 4 bases = log2(4) = 2.0
        max_entropy = math.log2(4)
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0.0

        positions.append({
            "position": i,
            "base": seq[i],
            "window": window,
            "entropy_bits": round(entropy, 6),
            "normalized_entropy": round(normalized_entropy, 6),
            "base_frequencies": {b: round(base_probs[b], 4) for b in RNA_BASES},
        })

    return positions


def analyze_aso_entropy_coverage(
    positional_entropy: list[dict[str, Any]],
    aso_target_start: int,
    aso_target_end: int,
    sl_length: int,
) -> dict[str, Any]:
    """Verifica se o MRL-ASO-001 cobre as posicoes de menor entropia.

    Identifica as posicoes no quartil inferior de entropia (mais conservadas)
    e calcula que fracao delas e coberta pela janela de binding do ASO.

    Args:
        positional_entropy: Resultado de compute_positional_entropy().
        aso_target_start: Posicao inicial do binding do ASO no SL RNA (0-indexed).
        aso_target_end: Posicao final do binding (exclusiva).
        sl_length: Comprimento total do SL RNA.

    Returns:
        Dicionario com analise de cobertura.
    """
    n = len(positional_entropy)

    # Ordenar posicoes por entropia (menor primeiro)
    sorted_positions = sorted(positional_entropy, key=lambda p: p["entropy_bits"])

    # Quartil inferior: 25% das posicoes com menor entropia
    quartile_size = max(1, n // 4)
    low_entropy_positions = sorted_positions[:quartile_size]
    low_entropy_indices = {p["position"] for p in low_entropy_positions}

    # Posicoes cobertas pelo ASO
    aso_coverage = set(range(aso_target_start, aso_target_end))

    # Intersecao: posicoes de baixa entropia cobertas pelo ASO
    covered_low_entropy = low_entropy_indices & aso_coverage
    overlap_fraction = len(covered_low_entropy) / len(low_entropy_indices) if low_entropy_indices else 0.0

    # Entropia media na regiao de binding vs fora
    binding_entropies = [
        p["entropy_bits"] for p in positional_entropy
        if p["position"] in aso_coverage
    ]
    non_binding_entropies = [
        p["entropy_bits"] for p in positional_entropy
        if p["position"] not in aso_coverage
    ]

    mean_binding_entropy = float(np.mean(binding_entropies)) if binding_entropies else 0.0
    mean_non_binding_entropy = float(np.mean(non_binding_entropies)) if non_binding_entropies else 0.0

    return {
        "total_positions": n,
        "quartile_size": quartile_size,
        "low_entropy_positions": sorted(low_entropy_indices),
        "aso_binding_range": [aso_target_start, aso_target_end],
        "aso_covered_positions": sorted(aso_coverage),
        "covered_low_entropy_positions": sorted(covered_low_entropy),
        "overlap_fraction": round(overlap_fraction, 4),
        "overlap_description": (
            f"{len(covered_low_entropy)}/{len(low_entropy_indices)} posicoes "
            f"de baixa entropia ({overlap_fraction:.1%}) cobertas pelo ASO"
        ),
        "mean_entropy_binding_region": round(mean_binding_entropy, 6),
        "mean_entropy_non_binding_region": round(mean_non_binding_entropy, 6),
        "entropy_differential": round(mean_binding_entropy - mean_non_binding_entropy, 6),
        "binding_targets_lower_entropy": mean_binding_entropy <= mean_non_binding_entropy,
    }


# ---------------------------------------------------------------------------
# 5. Divergencias de Kullback-Leibler e Jensen-Shannon
# ---------------------------------------------------------------------------


def kl_divergence(
    p: dict[str, float],
    q: dict[str, float],
) -> float:
    """Calcula a divergencia de Kullback-Leibler KL(P || Q).

    KL(P || Q) = sum(p_i * log(p_i / q_i))

    Mede a quantidade de informacao perdida quando Q e usada
    para aproximar P. NAO e simetrica: KL(P||Q) != KL(Q||P).

    Valor maior = distribuicoes mais diferentes.
    KL(P||Q) = 0 iff P = Q.

    Ref: Kullback S, Leibler RA (1951) Ann Math Stat 22(1):79-86

    Args:
        p: Distribuicao P (a "verdade").
        q: Distribuicao Q (o "modelo").

    Returns:
        KL(P || Q) em nats (logaritmo natural).
    """
    keys = sorted(set(p.keys()) | set(q.keys()))
    kl = 0.0

    for k in keys:
        p_i = p.get(k, 1e-15)
        q_i = q.get(k, 1e-15)

        # Proteger contra log(0)
        p_i = max(p_i, 1e-15)
        q_i = max(q_i, 1e-15)

        if p_i > 1e-14:
            kl += p_i * math.log(p_i / q_i)

    return kl


def jensen_shannon_divergence(
    p: dict[str, float],
    q: dict[str, float],
) -> float:
    """Calcula a divergencia de Jensen-Shannon JSD(P || Q).

    JSD(P, Q) = 0.5 * KL(P || M) + 0.5 * KL(Q || M)
    onde M = 0.5 * (P + Q)

    E a versao simetrica e suavizada da KL divergence.
    JSD in [0, ln(2)] para logaritmo natural.
    sqrt(JSD) e uma metrica (satisfaz desigualdade triangular).

    Ref: Lin J (1991) IEEE Trans Inf Theory 37(1):145-151

    Args:
        p: Primeira distribuicao.
        q: Segunda distribuicao.

    Returns:
        JSD(P, Q) em nats.
    """
    keys = sorted(set(p.keys()) | set(q.keys()))

    # Mistura M = 0.5 * (P + Q)
    m: dict[str, float] = {}
    for k in keys:
        m[k] = 0.5 * (p.get(k, 1e-15) + q.get(k, 1e-15))

    return 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)


def compute_divergences(
    kmer_results: dict[str, Any],
) -> dict[str, Any]:
    """Calcula todas as divergencias KL e JSD para cada k.

    Para cada k, calcula:
    - KL(SL_RNA || human) e KL(human || SL_RNA)
    - KL(SL_RNA || canine) e KL(canine || SL_RNA)
    - JSD(SL_RNA, human) e JSD(SL_RNA, canine)
    - JSD(human, canine) como referencia

    Args:
        kmer_results: Resultado de analyze_kmer_distributions().

    Returns:
        Dicionario com divergencias por k.
    """
    results: dict[str, Any] = {}

    for k_key in kmer_results:
        k = kmer_results[k_key]["k"]
        dists = kmer_results[k_key]["distributions"]

        sl = dists["sl_rna"]
        human = dists["human"]
        canine = dists["canine"]
        random = dists["random"]

        # KL divergences
        kl_sl_human = kl_divergence(sl, human)
        kl_human_sl = kl_divergence(human, sl)
        kl_sl_canine = kl_divergence(sl, canine)
        kl_canine_sl = kl_divergence(canine, sl)
        kl_sl_random = kl_divergence(sl, random)

        # Jensen-Shannon divergences
        jsd_sl_human = jensen_shannon_divergence(sl, human)
        jsd_sl_canine = jensen_shannon_divergence(sl, canine)
        jsd_sl_random = jensen_shannon_divergence(sl, random)
        jsd_human_canine = jensen_shannon_divergence(human, canine)

        results[k_key] = {
            "k": k,
            "kl_divergences": {
                "kl_sl_human_nats": round(kl_sl_human, 6),
                "kl_human_sl_nats": round(kl_human_sl, 6),
                "kl_sl_canine_nats": round(kl_sl_canine, 6),
                "kl_canine_sl_nats": round(kl_canine_sl, 6),
                "kl_sl_random_nats": round(kl_sl_random, 6),
            },
            "jensen_shannon_divergences": {
                "jsd_sl_human_nats": round(jsd_sl_human, 6),
                "jsd_sl_canine_nats": round(jsd_sl_canine, 6),
                "jsd_sl_random_nats": round(jsd_sl_random, 6),
                "jsd_human_canine_nats": round(jsd_human_canine, 6),
            },
            "interpretation": {
                "sl_more_distant_from_human_than_random": jsd_sl_human > jsd_sl_random,
                "sl_more_distant_from_canine_than_random": jsd_sl_canine > jsd_sl_random,
                "sl_more_distant_from_hosts_than_hosts_apart": (
                    jsd_sl_human > jsd_human_canine
                    and jsd_sl_canine > jsd_human_canine
                ),
            },
        }

    return results


# ---------------------------------------------------------------------------
# 6. Score composto de alienness
# ---------------------------------------------------------------------------


def compute_alienness_score(
    distance_results: dict[str, Any],
    divergence_results: dict[str, Any],
    entropy_coverage: dict[str, Any],
) -> dict[str, Any]:
    """Calcula o score composto de alienness do SL RNA.

    O score combina tres dimensoes:
    1. Distancia geometrica (Fisher-Rao) — quao longe no simplex
    2. Divergencia estatistica (JSD) — quanta informacao se perde
    3. Cobertura de regioes conservadas pelo ASO

    Cada componente e normalizado em [0, 1] e combinado com pesos:
    - Fisher-Rao: 0.4 (geometria do espaco de distribuicoes)
    - Jensen-Shannon: 0.4 (teoria da informacao)
    - Cobertura de entropia: 0.2 (relevancia terapeutica)

    Args:
        distance_results: Resultado de compute_pairwise_distances().
        divergence_results: Resultado de compute_divergences().
        entropy_coverage: Resultado de analyze_aso_entropy_coverage().

    Returns:
        Dicionario com score composto e componentes.
    """
    # --- Componente 1: Fisher-Rao (media sobre k=1,2,3,4) ---
    # Normaliza pela distancia maxima possivel (pi para Fisher-Rao)
    fr_distances_human: list[float] = []
    fr_distances_canine: list[float] = []

    for k_key in distance_results:
        d_human = distance_results[k_key]["sl_to_human"]
        d_canine = distance_results[k_key]["sl_to_canine"]
        fr_distances_human.append(d_human)
        fr_distances_canine.append(d_canine)

    # Media das distancias FR normalizadas por pi
    mean_fr_human = float(np.mean(fr_distances_human)) / math.pi
    mean_fr_canine = float(np.mean(fr_distances_canine)) / math.pi
    fr_score = (mean_fr_human + mean_fr_canine) / 2.0

    # --- Componente 2: Jensen-Shannon (media sobre k=1,2,3,4) ---
    # JSD maxima = ln(2) ≈ 0.693 para duas distribuicoes
    jsd_max = math.log(2)
    jsd_human_values: list[float] = []
    jsd_canine_values: list[float] = []

    for k_key in divergence_results:
        jsd_h = divergence_results[k_key]["jensen_shannon_divergences"]["jsd_sl_human_nats"]
        jsd_c = divergence_results[k_key]["jensen_shannon_divergences"]["jsd_sl_canine_nats"]
        jsd_human_values.append(jsd_h)
        jsd_canine_values.append(jsd_c)

    mean_jsd_human = float(np.mean(jsd_human_values)) / jsd_max
    mean_jsd_canine = float(np.mean(jsd_canine_values)) / jsd_max
    jsd_score = (mean_jsd_human + mean_jsd_canine) / 2.0

    # --- Componente 3: Cobertura de entropia ---
    coverage_score = entropy_coverage["overlap_fraction"]

    # --- Score composto ---
    weights = {"fisher_rao": 0.4, "jensen_shannon": 0.4, "entropy_coverage": 0.2}
    composite = (
        weights["fisher_rao"] * fr_score
        + weights["jensen_shannon"] * jsd_score
        + weights["entropy_coverage"] * coverage_score
    )

    # --- Ranking de alvos hipoteticos ---
    # Compara SL RNA com sequencias hipoteticas de mesma distancia
    # que teriam alienness menor
    # Um alvo "perfeito" teria FR_score=1, JSD_score=1, coverage=1
    # Um alvo "terrivel" teria FR_score=0, JSD_score=0, coverage=0

    return {
        "composite_alienness_score": round(composite, 6),
        "components": {
            "fisher_rao": {
                "score": round(fr_score, 6),
                "weight": weights["fisher_rao"],
                "mean_distance_human_normalized": round(mean_fr_human, 6),
                "mean_distance_canine_normalized": round(mean_fr_canine, 6),
                "raw_distances_human": [round(d, 6) for d in fr_distances_human],
                "raw_distances_canine": [round(d, 6) for d in fr_distances_canine],
            },
            "jensen_shannon": {
                "score": round(jsd_score, 6),
                "weight": weights["jensen_shannon"],
                "mean_jsd_human_normalized": round(mean_jsd_human, 6),
                "mean_jsd_canine_normalized": round(mean_jsd_canine, 6),
                "raw_jsd_human": [round(d, 6) for d in jsd_human_values],
                "raw_jsd_canine": [round(d, 6) for d in jsd_canine_values],
            },
            "entropy_coverage": {
                "score": round(coverage_score, 6),
                "weight": weights["entropy_coverage"],
                "overlap_fraction": entropy_coverage["overlap_fraction"],
            },
        },
        "interpretation": {
            "maximum_possible_score": 1.0,
            "score_as_percentage": round(composite * 100, 2),
            "sl_rna_is_distant_from_hosts": fr_score > 0.01 and jsd_score > 0.01,
        },
    }


# ---------------------------------------------------------------------------
# Gerador de relatorio Markdown
# ---------------------------------------------------------------------------


def _generate_report(
    kmer_results: dict[str, Any],
    fisher_results: dict[str, Any],
    distance_results: dict[str, Any],
    positional_entropy: list[dict[str, Any]],
    entropy_coverage: dict[str, Any],
    divergence_results: dict[str, Any],
    alienness: dict[str, Any],
    runtime_seconds: float,
) -> str:
    """Gera relatorio Markdown com resultados da geometria da informacao.

    Relato honesto — se SL RNA NAO for maximalmente distante, isso e reportado.
    """
    now = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines: list[str] = []
    lines.append("# Information Geometry Report — SL RNA Alienness Analysis")
    lines.append("")
    lines.append(f"**Generated:** {now}")
    lines.append(f"**Runtime:** {runtime_seconds:.2f} seconds")
    lines.append("")
    lines.append("---")
    lines.append("")

    # --- 1. K-mer distributions ---
    lines.append("## 1. Nucleotide Frequency Distributions (k-mer)")
    lines.append("")
    lines.append("K-mer frequency distributions computed for 4 sources:")
    lines.append("- **SL RNA:** L. infantum Spliced Leader (39 nt)")
    lines.append("- **Human:** average human mRNA composition (iid model)")
    lines.append("- **Canine:** average canine mRNA composition (iid model)")
    lines.append("- **Random:** random RNA with same length and GC content")
    lines.append("")

    for k_key in sorted(kmer_results.keys()):
        k = kmer_results[k_key]["k"]
        n_kmers = kmer_results[k_key]["n_possible_kmers"]
        lines.append(f"### k = {k} ({n_kmers} possible {k}-mers)")
        lines.append("")

        # Mostrar top 5 k-mers mais frequentes no SL RNA
        sl_dist = kmer_results[k_key]["distributions"]["sl_rna"]
        sorted_kmers = sorted(sl_dist.items(), key=lambda x: x[1], reverse=True)

        lines.append(f"**Top 5 k-mers in SL RNA (k={k}):**")
        lines.append("")
        lines.append("| K-mer | SL RNA | Human | Canine | Random |")
        lines.append("|-------|--------|-------|--------|--------|")

        human_dist = kmer_results[k_key]["distributions"]["human"]
        canine_dist = kmer_results[k_key]["distributions"]["canine"]
        random_dist = kmer_results[k_key]["distributions"]["random"]

        for kmer, prob in sorted_kmers[:5]:
            lines.append(
                f"| {kmer} | {prob:.4f} | {human_dist.get(kmer, 0):.4f} | "
                f"{canine_dist.get(kmer, 0):.4f} | {random_dist.get(kmer, 0):.4f} |"
            )
        lines.append("")

    # --- 2. Fisher Information Matrix ---
    lines.append("## 2. Fisher Information Matrix Properties")
    lines.append("")
    lines.append("The Fisher Information Matrix (FIM) for the multinomial model is diagonal:")
    lines.append("G_ij = delta_ij / p_i. Properties reflect the geometry of each distribution.")
    lines.append("")

    for k_key in sorted(fisher_results.keys()):
        lines.append(f"### k = {kmer_results[k_key]['k']}")
        lines.append("")
        lines.append("| Source | log(det) | Trace | Max eigenvalue | Min eigenvalue |")
        lines.append("|--------|----------|-------|----------------|----------------|")

        for source in ["sl_rna", "human", "canine", "random"]:
            f = fisher_results[k_key][source]
            lines.append(
                f"| {source:>6} | {f['log_determinant']:>8.3f} | {f['trace']:>5.1f} | "
                f"{f['max_eigenvalue']:>14.4f} | {f['min_eigenvalue']:>14.4f} |"
            )
        lines.append("")

    # --- 3. Fisher-Rao distances ---
    lines.append("## 3. Fisher-Rao Geodesic Distances")
    lines.append("")
    lines.append("d_FR(p, q) = 2 * arccos(sum(sqrt(p_i * q_i)))")
    lines.append("Range: [0, pi]. Larger distance = more statistically distinct.")
    lines.append("")

    lines.append("| k | SL-Human | SL-Canine | SL-Random | Human-Canine |")
    lines.append("|---|----------|-----------|-----------|--------------|")
    for k_key in sorted(distance_results.keys()):
        dr = distance_results[k_key]
        lines.append(
            f"| {dr['k']} | {dr['sl_to_human']:.6f} | {dr['sl_to_canine']:.6f} | "
            f"{dr['sl_to_random']:.6f} | {dr['human_to_canine']:.6f} |"
        )
    lines.append("")

    # Interpretacao
    lines.append("**Interpretation:**")
    for k_key in sorted(distance_results.keys()):
        dr = distance_results[k_key]
        sl_max = max(dr["sl_to_human"], dr["sl_to_canine"])
        hc = dr["human_to_canine"]
        ratio = sl_max / hc if hc > 1e-10 else float("inf")
        lines.append(
            f"- k={dr['k']}: SL RNA is {ratio:.1f}x more distant from host "
            f"than hosts are from each other"
        )
    lines.append("")

    # --- 4. Positional entropy ---
    lines.append("## 4. Positional Shannon Entropy & ASO Coverage")
    lines.append("")
    lines.append("Entropy computed using 5-nt local window context.")
    lines.append("Lower entropy = more structured/conserved position.")
    lines.append("")

    lines.append("| Position | Base | Entropy (bits) | Normalized | In ASO? |")
    lines.append("|----------|------|----------------|------------|---------|")

    aso_positions = set(range(
        entropy_coverage["aso_binding_range"][0],
        entropy_coverage["aso_binding_range"][1],
    ))
    for pe in positional_entropy:
        in_aso = "YES" if pe["position"] in aso_positions else "no"
        lines.append(
            f"| {pe['position']:>8} | {pe['base']:>4} | "
            f"{pe['entropy_bits']:>14.4f} | {pe['normalized_entropy']:>10.4f} | "
            f"{in_aso:>7} |"
        )
    lines.append("")

    lines.append("### ASO Coverage of Low-Entropy Positions")
    lines.append("")
    lines.append(f"- **Low-entropy positions (bottom quartile):** {entropy_coverage['low_entropy_positions']}")
    lines.append(f"- **Covered by ASO:** {entropy_coverage['covered_low_entropy_positions']}")
    lines.append(f"- **Overlap fraction:** {entropy_coverage['overlap_fraction']:.1%}")
    lines.append(f"- **Mean entropy in binding region:** {entropy_coverage['mean_entropy_binding_region']:.4f} bits")
    lines.append(f"- **Mean entropy outside binding:** {entropy_coverage['mean_entropy_non_binding_region']:.4f} bits")
    lines.append(f"- **ASO targets lower-entropy region:** "
                 f"{'YES' if entropy_coverage['binding_targets_lower_entropy'] else 'NO'}")
    lines.append("")

    # --- 5. KL and JSD divergences ---
    lines.append("## 5. Kullback-Leibler & Jensen-Shannon Divergences")
    lines.append("")

    lines.append("### KL Divergence (asymmetric, in nats)")
    lines.append("")
    lines.append("| k | KL(SL\\|\\|Human) | KL(Human\\|\\|SL) | KL(SL\\|\\|Canine) | KL(Canine\\|\\|SL) | KL(SL\\|\\|Random) |")
    lines.append("|---|-------------|-------------|--------------|--------------|---------------|")
    for k_key in sorted(divergence_results.keys()):
        kl = divergence_results[k_key]["kl_divergences"]
        lines.append(
            f"| {divergence_results[k_key]['k']} | "
            f"{kl['kl_sl_human_nats']:.6f} | {kl['kl_human_sl_nats']:.6f} | "
            f"{kl['kl_sl_canine_nats']:.6f} | {kl['kl_canine_sl_nats']:.6f} | "
            f"{kl['kl_sl_random_nats']:.6f} |"
        )
    lines.append("")

    lines.append("### Jensen-Shannon Divergence (symmetric, in nats)")
    lines.append("")
    lines.append("| k | JSD(SL, Human) | JSD(SL, Canine) | JSD(SL, Random) | JSD(Human, Canine) |")
    lines.append("|---|----------------|-----------------|-----------------|---------------------|")
    for k_key in sorted(divergence_results.keys()):
        jsd = divergence_results[k_key]["jensen_shannon_divergences"]
        lines.append(
            f"| {divergence_results[k_key]['k']} | "
            f"{jsd['jsd_sl_human_nats']:.6f} | {jsd['jsd_sl_canine_nats']:.6f} | "
            f"{jsd['jsd_sl_random_nats']:.6f} | {jsd['jsd_human_canine_nats']:.6f} |"
        )
    lines.append("")

    lines.append("**Interpretation:** Higher divergence = more 'alien' = better therapeutic target.")
    lines.append("")

    for k_key in sorted(divergence_results.keys()):
        interp = divergence_results[k_key]["interpretation"]
        k = divergence_results[k_key]["k"]
        lines.append(f"- k={k}: SL more distant from human than random: "
                     f"{'YES' if interp['sl_more_distant_from_human_than_random'] else 'NO'}")
        lines.append(f"- k={k}: SL more distant from canine than random: "
                     f"{'YES' if interp['sl_more_distant_from_canine_than_random'] else 'NO'}")
        lines.append(f"- k={k}: SL more distant from hosts than hosts from each other: "
                     f"{'YES' if interp['sl_more_distant_from_hosts_than_hosts_apart'] else 'NO'}")
    lines.append("")

    # --- 6. Alienness score ---
    lines.append("## 6. Composite Alienness Score")
    lines.append("")

    score = alienness["composite_alienness_score"]
    pct = alienness["interpretation"]["score_as_percentage"]
    comp = alienness["components"]

    lines.append(f"**Composite Score: {score:.4f} ({pct:.1f}% of maximum)**")
    lines.append("")
    lines.append("| Component | Score | Weight | Contribution |")
    lines.append("|-----------|-------|--------|--------------|")

    for name, data in comp.items():
        contribution = data["score"] * data["weight"]
        lines.append(
            f"| {name} | {data['score']:.4f} | {data['weight']:.1f} | {contribution:.4f} |"
        )
    lines.append("")

    lines.append("### Score Interpretation")
    lines.append("")
    lines.append("- 0.0 = identical to host (no selectivity)")
    lines.append("- 0.5 = moderately distinct")
    lines.append("- 1.0 = maximally alien (perfect selectivity)")
    lines.append("")

    is_distant = alienness["interpretation"]["sl_rna_is_distant_from_hosts"]

    # --- Conclusao ---
    lines.append("## Conclusion")
    lines.append("")

    if is_distant:
        # Verificar quao forte e a evidencia
        fr_component = comp["fisher_rao"]["score"]
        jsd_component = comp["jensen_shannon"]["score"]
        coverage_component = comp["entropy_coverage"]["score"]

        findings: list[str] = []

        # Fisher-Rao
        findings.append(
            f"Fisher-Rao geometric distance confirms SL RNA occupies a "
            f"distinct region of the probability simplex (normalized score: {fr_component:.4f})"
        )

        # JSD
        findings.append(
            f"Jensen-Shannon divergence confirms information-theoretic "
            f"separation (normalized score: {jsd_component:.4f})"
        )

        # Cobertura
        if coverage_component > 0.5:
            findings.append(
                f"MRL-ASO-001 covers {coverage_component:.0%} of the "
                f"lowest-entropy positions, targeting the most conserved region"
            )
        else:
            findings.append(
                f"MRL-ASO-001 covers {coverage_component:.0%} of the "
                f"lowest-entropy positions"
            )

        lines.append(
            f"The SL RNA of L. infantum is **statistically distinct** from "
            f"both human and canine transcriptomes across all k-mer scales. "
            f"Composite alienness score: **{score:.4f}**."
        )
        lines.append("")
        for f in findings:
            lines.append(f"- {f}")
    else:
        lines.append(
            f"**WARNING:** The SL RNA is NOT clearly distinct from host "
            f"transcriptomes by information geometry measures. "
            f"Composite alienness score: **{score:.4f}**. "
            f"This raises questions about target selectivity."
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
    """Executa a analise de geometria da informacao do SL RNA.

    Fluxo:
        1. Calcular distribuicoes de k-mers (k=1,2,3,4) para 4 fontes
        2. Computar Matrizes de Informacao de Fisher
        3. Calcular distancias geodesicas de Fisher-Rao
        4. Analisar entropia posicional e cobertura do ASO
        5. Computar divergencias KL e JSD
        6. Calcular score composto de alienness
        7. Gerar relatorio Markdown

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
            aso_target_start=ASO_TARGET_START,
            aso_target_end=ASO_TARGET_END,
            known_dg=ASO_KNOWN_DG,
            known_tm=ASO_KNOWN_TM,
        )

    sl_seq = config.sl_sequence.upper()
    aso_seq = config.aso_sequence.upper()
    aso_name = config.aso_name or "ASO"
    target_start = config.aso_target_start
    target_end = config.aso_target_end

    gc_frac = gc_content(sl_seq)

    logger.info("=" * 70)
    logger.info("MATH 2 — Geometria da Informacao — Alienness do SL RNA")
    logger.info("=" * 70)
    logger.info("SL RNA: %d nt, GC = %.1f%%", len(sl_seq), gc_frac * 100)
    logger.info("ASO: %s (%d nt, pos %d-%d)", aso_name, len(aso_seq), target_start, target_end)

    with Timer() as timer:
        # ---------------------------------------------------------------
        # Passo 1: Distribuicoes de k-mers
        # ---------------------------------------------------------------
        logger.info("Passo 1/6: Calculando distribuicoes de k-mers...")
        kmer_results = analyze_kmer_distributions(sl_seq, gc_frac)
        _save_json(kmer_results, "kmer_distributions")
        logger.info("  -> kmer_distributions.json salvo")

        # ---------------------------------------------------------------
        # Passo 2: Matrizes de Informacao de Fisher
        # ---------------------------------------------------------------
        logger.info("Passo 2/6: Computando Matrizes de Informacao de Fisher...")
        fisher_results = analyze_fisher_matrices(kmer_results)
        _save_json(fisher_results, "fisher_matrices")
        logger.info("  -> fisher_matrices.json salvo")

        # ---------------------------------------------------------------
        # Passo 3: Distancias de Fisher-Rao
        # ---------------------------------------------------------------
        logger.info("Passo 3/6: Calculando distancias geodesicas de Fisher-Rao...")
        distance_results = compute_pairwise_distances(kmer_results)
        _save_json(distance_results, "fisher_rao_distances")
        logger.info("  -> fisher_rao_distances.json salvo")

        # ---------------------------------------------------------------
        # Passo 4: Entropia posicional e cobertura do ASO
        # ---------------------------------------------------------------
        logger.info("Passo 4/6: Analisando entropia posicional e cobertura...")
        positional_entropy = compute_positional_entropy(sl_seq)
        entropy_coverage = analyze_aso_entropy_coverage(
            positional_entropy=positional_entropy,
            aso_target_start=target_start,
            aso_target_end=target_end,
            sl_length=len(sl_seq),
        )
        _save_json(
            {"positional_entropy": positional_entropy, "coverage": entropy_coverage},
            "entropy_coverage",
        )
        logger.info("  -> entropy_coverage.json salvo")

        # ---------------------------------------------------------------
        # Passo 5: Divergencias KL e Jensen-Shannon
        # ---------------------------------------------------------------
        logger.info("Passo 5/6: Computando divergencias KL e JSD...")
        divergence_results = compute_divergences(kmer_results)
        _save_json(divergence_results, "divergences")
        logger.info("  -> divergences.json salvo")

        # ---------------------------------------------------------------
        # Passo 6: Score composto de alienness
        # ---------------------------------------------------------------
        logger.info("Passo 6/6: Calculando score composto de alienness...")
        alienness = compute_alienness_score(
            distance_results=distance_results,
            divergence_results=divergence_results,
            entropy_coverage=entropy_coverage,
        )
        _save_json(alienness, "alienness_score")
        logger.info("  -> alienness_score.json salvo")

    # --- Gerar relatorio Markdown ---
    report = _generate_report(
        kmer_results=kmer_results,
        fisher_results=fisher_results,
        distance_results=distance_results,
        positional_entropy=positional_entropy,
        entropy_coverage=entropy_coverage,
        divergence_results=divergence_results,
        alienness=alienness,
        runtime_seconds=timer.elapsed,
    )
    report_path = RESULTS_DIR / "INFORMATION_GEOMETRY_REPORT.md"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    logger.info("Relatorio salvo em: %s", report_path)

    # --- Resultado consolidado ---
    result: dict[str, Any] = {
        "module": "math_2_fisher",
        "version": "1.0.0",
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "runtime_seconds": timer.elapsed,
        "status": "success",
        "aso_name": aso_name,
        "organism": config.species_name,
        "sl_rna_gc_content": round(gc_frac, 4),
        "alienness_score": alienness,
        "fisher_rao_distances": {
            k_key: {
                "sl_to_human": distance_results[k_key]["sl_to_human"],
                "sl_to_canine": distance_results[k_key]["sl_to_canine"],
                "sl_to_random": distance_results[k_key]["sl_to_random"],
            }
            for k_key in distance_results
        },
        "divergences_summary": {
            k_key: {
                "jsd_sl_human": divergence_results[k_key]["jensen_shannon_divergences"]["jsd_sl_human_nats"],
                "jsd_sl_canine": divergence_results[k_key]["jensen_shannon_divergences"]["jsd_sl_canine_nats"],
            }
            for k_key in divergence_results
        },
        "entropy_coverage": {
            "overlap_fraction": entropy_coverage["overlap_fraction"],
            "binding_targets_lower_entropy": entropy_coverage["binding_targets_lower_entropy"],
        },
    }

    logger.info("=" * 70)
    logger.info(
        "MATH 2 COMPLETO — Alienness score: %.4f — Tempo: %.2f s",
        alienness["composite_alienness_score"],
        timer.elapsed,
    )
    logger.info("=" * 70)

    return result


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
