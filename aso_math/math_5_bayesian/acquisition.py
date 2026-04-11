"""Funcao de aquisicao Expected Improvement (EI) para otimizacao bayesiana.

A EI balanceia exploracao (regioes com alta incerteza) e exploitation
(regioes com media favoravel), guiando a busca para o proximo ponto
mais informativo no espaco de design.

Referencia: Jones DR, Schonlau M, Welch WJ (1998)
    "Efficient Global Optimization of Expensive Black-Box Functions"
    J Global Optimization 13(4):455-492
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from aso_math.math_5_bayesian.gaussian_process import GaussianProcess


# ---------------------------------------------------------------------------
# Expected Improvement
# ---------------------------------------------------------------------------


def expected_improvement(
    X_candidates: np.ndarray,
    gp: GaussianProcess,
    f_best: float,
    xi: float = 0.01,
) -> np.ndarray:
    """Calcula Expected Improvement para um conjunto de pontos candidatos.

    EI(x) = (f_best - mu(x)) * Phi(z) + sigma(x) * phi(z)

    onde:
        z = (f_best - mu(x)) / sigma(x)
        Phi = CDF da normal padrao
        phi = PDF da normal padrao

    Estamos MINIMIZANDO o objetivo composto, entao f_best e o menor
    valor observado e queremos mu(x) < f_best.

    O parametro xi controla exploracao vs exploitation:
        xi = 0: pura exploitation
        xi > 0: mais exploracao (evita convergencia prematura)

    Args:
        X_candidates: Pontos candidatos, shape (n, d).
        gp: Gaussian Process ajustado.
        f_best: Melhor valor observado ate agora (menor = melhor).
        xi: Parametro de exploracao (default 0.01).

    Returns:
        Vetor de EI, shape (n,). Valores maiores = mais promissores.
    """
    X_candidates = np.asarray(X_candidates, dtype=np.float64)
    if X_candidates.ndim == 1:
        X_candidates = X_candidates.reshape(1, -1)

    mu, sigma = gp.predict(X_candidates, return_std=True)

    # Evitar divisao por zero em regioes com variancia muito baixa
    # (pontos ja observados ou muito proximos dos dados de treino)
    mask = sigma > 1e-10
    ei = np.zeros_like(mu)

    if np.any(mask):
        improvement = f_best - mu[mask] - xi
        z = improvement / sigma[mask]
        ei[mask] = improvement * norm.cdf(z) + sigma[mask] * norm.pdf(z)

    # EI nao pode ser negativa (por construcao, mas arredondamento...)
    ei = np.maximum(ei, 0.0)

    return ei


# ---------------------------------------------------------------------------
# Selecao do proximo ponto via EI
# ---------------------------------------------------------------------------


def select_next_point(
    gp: GaussianProcess,
    f_best: float,
    bounds: np.ndarray,
    n_candidates: int = 5000,
    rng: np.random.Generator | None = None,
    xi: float = 0.01,
) -> np.ndarray:
    """Seleciona o proximo ponto a avaliar maximizando a EI.

    Usa uma estrategia de amostragem aleatoria: gera n_candidates pontos
    uniformes no espaco de busca e retorna o com maior EI. Para 4 dimensoes
    e 5000 candidatos, a cobertura e adequada sem precisar de otimizacao
    gradiente-based (que requer derivadas do GP).

    Args:
        gp: GP ajustado aos dados observados.
        f_best: Melhor valor objetivo observado.
        bounds: Limites do espaco, shape (d, 2) com [min, max] por dimensao.
        n_candidates: Numero de candidatos aleatorios a gerar.
        rng: Gerador de numeros aleatorios. Se None, cria um novo.
        xi: Parametro de exploracao para EI.

    Returns:
        Vetor de parametros do ponto selecionado, shape (d,).
    """
    if rng is None:
        rng = np.random.default_rng()

    d = bounds.shape[0]

    # Gerar candidatos uniformes dentro dos bounds
    candidates = np.zeros((n_candidates, d))
    for i in range(d):
        candidates[:, i] = rng.uniform(bounds[i, 0], bounds[i, 1], size=n_candidates)

    # Calcular EI para todos os candidatos
    ei_values = expected_improvement(candidates, gp, f_best, xi=xi)

    # Retornar o candidato com maior EI
    best_idx = np.argmax(ei_values)

    return candidates[best_idx]
