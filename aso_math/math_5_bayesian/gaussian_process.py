"""Gaussian Process (GP) implementado from scratch com kernel RBF.

Implementacao completa sem dependencias de ML externas. Usa apenas
numpy para algebra linear e scipy.linalg para decomposicao de Cholesky
(numericamente estavel). O GP modela a funcao objetivo como uma
distribuicao a posteriori sobre funcoes, permitindo quantificar
incerteza nas predicoes.

Referencia: Rasmussen & Williams (2006) "Gaussian Processes for Machine Learning"
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve


# ---------------------------------------------------------------------------
# Kernel RBF (Radial Basis Function / Squared Exponential)
# ---------------------------------------------------------------------------


def rbf_kernel(
    x1: np.ndarray,
    x2: np.ndarray,
    length_scale: float,
    signal_variance: float,
) -> np.ndarray:
    """Calcula a matriz de covariancia usando kernel RBF (squared exponential).

    k(x, x') = sigma_f^2 * exp(-||x - x'||^2 / (2 * l^2))

    O kernel RBF assume que funcoes proximas no espaco de entrada produzem
    valores proximos na saida — suavidade infinita. E o kernel padrao para
    GP quando nao ha informacao a priori sobre a estrutura da funcao.

    Args:
        x1: Matriz (n1, d) de pontos de entrada.
        x2: Matriz (n2, d) de pontos de entrada.
        length_scale: Escala de comprimento (l). Controla a suavidade —
            valores maiores = funcao mais suave.
        signal_variance: Variancia do sinal (sigma_f^2). Controla a
            amplitude das variacoes.

    Returns:
        Matriz de covariancia (n1, n2).
    """
    # Calcular distancia euclidiana ao quadrado entre todos os pares
    # ||x1_i - x2_j||^2 = ||x1_i||^2 + ||x2_j||^2 - 2 * x1_i . x2_j
    sq_dist = (
        np.sum(x1 ** 2, axis=1, keepdims=True)
        + np.sum(x2 ** 2, axis=1, keepdims=True).T
        - 2.0 * x1 @ x2.T
    )

    # Garantir que distancias sao nao-negativas (erros de ponto flutuante)
    sq_dist = np.maximum(sq_dist, 0.0)

    return signal_variance * np.exp(-sq_dist / (2.0 * length_scale ** 2))


# ---------------------------------------------------------------------------
# Gaussian Process regressor
# ---------------------------------------------------------------------------


class GaussianProcess:
    """Gaussian Process regressor com kernel RBF e ruido gaussiano.

    O modelo assume:
        y = f(x) + epsilon,  epsilon ~ N(0, noise_variance)

    A posteriori sobre f(x*) dado dados de treino (X, y):
        mu(x*) = k(x*, X) @ (K + sigma_n^2 I)^-1 @ y
        var(x*) = k(x*, x*) - k(x*, X) @ (K + sigma_n^2 I)^-1 @ k(X, x*)

    A inversao e feita via decomposicao de Cholesky para estabilidade numerica.

    Atributos:
        length_scale: Hiperparametro l do kernel RBF.
        signal_variance: Hiperparametro sigma_f^2 do kernel.
        noise_variance: Variancia do ruido de observacao sigma_n^2.
        X_train: Dados de treino (n, d).
        y_train: Targets de treino (n,).
        _L: Fator de Cholesky de (K + sigma_n^2 I).
        _alpha: Vetor auxiliar (K + sigma_n^2 I)^-1 @ y.
    """

    def __init__(
        self,
        length_scale: float = 1.0,
        signal_variance: float = 1.0,
        noise_variance: float = 1e-6,
    ) -> None:
        self.length_scale = length_scale
        self.signal_variance = signal_variance
        self.noise_variance = noise_variance

        self.X_train: np.ndarray | None = None
        self.y_train: np.ndarray | None = None
        self._cho_factor: tuple | None = None
        self._alpha: np.ndarray | None = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> None:
        """Ajusta o GP aos dados de treino.

        Computa a decomposicao de Cholesky de (K + sigma_n^2 I) e o vetor
        alpha = (K + sigma_n^2 I)^-1 @ y, que sao reutilizados em todas
        as predicoes subsequentes.

        Args:
            X: Dados de treino, shape (n, d).
            y: Targets de treino, shape (n,). Devem ser normalizados.
        """
        self.X_train = np.asarray(X, dtype=np.float64)
        self.y_train = np.asarray(y, dtype=np.float64)

        n = self.X_train.shape[0]

        # Matriz de covariancia do treino
        K = rbf_kernel(
            self.X_train, self.X_train,
            self.length_scale, self.signal_variance,
        )

        # Adicionar ruido na diagonal: K + sigma_n^2 * I
        # O jitter extra de 1e-8 garante estabilidade numerica
        K_noisy = K + (self.noise_variance + 1e-8) * np.eye(n)

        # Decomposicao de Cholesky: K_noisy = L @ L^T
        self._cho_factor = cho_factor(K_noisy, lower=True)

        # alpha = (K + sigma_n^2 I)^-1 @ y
        self._alpha = cho_solve(self._cho_factor, self.y_train)

    def predict(
        self,
        X_test: np.ndarray,
        return_std: bool = True,
    ) -> tuple[np.ndarray, np.ndarray] | np.ndarray:
        """Predicao a posteriori do GP em pontos de teste.

        Args:
            X_test: Pontos de teste, shape (m, d).
            return_std: Se True, retorna tambem o desvio padrao.

        Returns:
            Se return_std=True: (mu, std) ambos shape (m,).
            Se return_std=False: mu shape (m,).
        """
        if self.X_train is None or self._cho_factor is None:
            raise RuntimeError("GP nao foi ajustado — chamar fit() primeiro.")

        X_test = np.asarray(X_test, dtype=np.float64)
        if X_test.ndim == 1:
            X_test = X_test.reshape(1, -1)

        # Covariancia entre teste e treino
        K_star = rbf_kernel(
            X_test, self.X_train,
            self.length_scale, self.signal_variance,
        )

        # Media a posteriori: mu = K_star @ alpha
        mu = K_star @ self._alpha

        if not return_std:
            return mu

        # Variancia a posteriori
        # var = k(x*, x*) - K_star @ (K + sigma_n^2 I)^-1 @ K_star^T
        # Usando Cholesky: v = L^-1 @ K_star^T, var = k(x*,x*) - v^T @ v
        v = cho_solve(self._cho_factor, K_star.T)
        var = self.signal_variance - np.sum(K_star.T * v, axis=0)

        # Garantir variancia nao-negativa (erros numericos)
        var = np.maximum(var, 1e-12)
        std = np.sqrt(var)

        return mu, std

    def log_marginal_likelihood(self) -> float:
        """Calcula a log-verossimilhanca marginal (para diagnostico).

        log p(y|X, theta) = -0.5 * y^T @ alpha - sum(log(diag(L))) - n/2 * log(2*pi)

        Valores mais altos indicam melhor ajuste dos hiperparametros.
        """
        if self._cho_factor is None or self.y_train is None:
            return -np.inf

        n = len(self.y_train)
        L = self._cho_factor[0]

        # -0.5 * y^T @ alpha
        data_fit = -0.5 * self.y_train @ self._alpha

        # -sum(log(diag(L)))
        complexity = -np.sum(np.log(np.diag(L)))

        # -n/2 * log(2*pi)
        constant = -0.5 * n * np.log(2.0 * np.pi)

        return float(data_fit + complexity + constant)


# ---------------------------------------------------------------------------
# Otimizacao de hiperparametros (grid search simples)
# ---------------------------------------------------------------------------


def optimize_hyperparameters(
    X: np.ndarray,
    y: np.ndarray,
    length_scales: tuple[float, ...] = (0.1, 0.3, 0.5, 1.0, 2.0),
    signal_variances: tuple[float, ...] = (0.5, 1.0, 2.0),
    noise_variance: float = 1e-4,
) -> GaussianProcess:
    """Seleciona hiperparametros por maximizacao da log-verossimilhanca marginal.

    Grid search simples sobre length_scale e signal_variance. Para 250
    avaliacoes, a grid e pequena o suficiente para ser computacionalmente
    trivial, e a log-marginal likelihood e um criterio principled que
    balanceia automaticamente data fit vs. complexidade (Occam's razor).

    Args:
        X: Dados de treino (n, d).
        y: Targets de treino (n,).
        length_scales: Candidatos para length_scale.
        signal_variances: Candidatos para signal_variance.
        noise_variance: Variancia do ruido (fixa).

    Returns:
        GP ajustado com os melhores hiperparametros encontrados.
    """
    best_lml = -np.inf
    best_gp = None

    for ls in length_scales:
        for sv in signal_variances:
            gp = GaussianProcess(
                length_scale=ls,
                signal_variance=sv,
                noise_variance=noise_variance,
            )
            try:
                gp.fit(X, y)
                lml = gp.log_marginal_likelihood()
                if lml > best_lml:
                    best_lml = lml
                    best_gp = gp
            except np.linalg.LinAlgError:
                # Cholesky falhou — combinacao de hiperparametros invalida
                continue

    if best_gp is None:
        # Fallback: usar defaults
        best_gp = GaussianProcess(
            length_scale=1.0,
            signal_variance=1.0,
            noise_variance=noise_variance,
        )
        best_gp.fit(X, y)

    return best_gp
