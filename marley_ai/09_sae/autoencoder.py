"""Sparse Autoencoder implementado puramente em numpy.

Arquitetura: input -> encoder (linear + ReLU) -> bottleneck esparso
             -> decoder (linear) -> reconstrucao

A funcao de perda combina:
  - MSE: fidelidade da reconstrucao (o decoder deve reconstruir o input)
  - L1 penalty: esparsidade das ativacoes (poucos neuronios ativos por input)

Treinamento com SGD e proximal step para L1 (soft thresholding),
que produz ativacoes exatamente zero — esparsidade verdadeira,
nao apenas valores pequenos.

Ref: Cunningham H et al. (2023) "Sparse Autoencoders Find Highly
     Interpretable Features in Language Models" arXiv:2309.08600
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SAEWeights:
    """Pesos do Sparse Autoencoder.

    Armazena pesos, biases e historico de treinamento.
    W_enc e W_dec podem ser atados (tied) para reduzir parametros.
    """

    W_enc: np.ndarray          # (input_dim, hidden_dim)
    b_enc: np.ndarray          # (hidden_dim,)
    W_dec: np.ndarray          # (hidden_dim, input_dim)
    b_dec: np.ndarray          # (input_dim,)
    tied: bool = True          # Se True, W_dec = W_enc^T
    loss_history: list[float] = field(default_factory=list)
    mse_history: list[float] = field(default_factory=list)
    l1_history: list[float] = field(default_factory=list)
    sparsity_history: list[float] = field(default_factory=list)


def _relu(x: np.ndarray) -> np.ndarray:
    """ReLU: max(0, x). Ativacao nao-linear que promove esparsidade."""
    return np.maximum(0.0, x)


def _relu_grad(x: np.ndarray) -> np.ndarray:
    """Gradiente da ReLU: 1 onde x > 0, 0 caso contrario."""
    return (x > 0).astype(np.float64)


def _soft_threshold(x: np.ndarray, threshold: float) -> np.ndarray:
    """Proximal operator para L1: produz zeros exatos.

    S_lambda(x) = sign(x) * max(|x| - lambda, 0)

    Essencial para esparsidade verdadeira — SGD sozinho
    produz valores proximos de zero mas nunca exatamente zero.
    """
    return np.sign(x) * np.maximum(np.abs(x) - threshold, 0.0)


def initialize_weights(
    input_dim: int,
    hidden_dim: int,
    tied: bool = True,
    seed: int = 42,
) -> SAEWeights:
    """Inicializa pesos do SAE com Xavier/Glorot.

    Xavier initialization: W ~ U(-sqrt(6/(fan_in+fan_out)), sqrt(6/(fan_in+fan_out)))
    Adequada para redes com ReLU — mantém variancia das ativacoes
    estavel entre camadas, prevenindo gradientes explodindo/desaparecendo.

    Args:
        input_dim: dimensao do input (72 para nosso encoding)
        hidden_dim: dimensao do espaco esparso (512)
        tied: se True, W_dec sera transposta de W_enc
        seed: semente para reproducibilidade

    Returns:
        SAEWeights inicializados
    """
    rng = np.random.default_rng(seed)

    # Xavier/Glorot uniform initialization
    limit_enc = np.sqrt(6.0 / (input_dim + hidden_dim))
    W_enc = rng.uniform(-limit_enc, limit_enc, (input_dim, hidden_dim))

    if tied:
        W_dec = W_enc.T.copy()
    else:
        limit_dec = np.sqrt(6.0 / (hidden_dim + input_dim))
        W_dec = rng.uniform(-limit_dec, limit_dec, (hidden_dim, input_dim))

    # Bias do encoder inicia NEGATIVO para promover esparsidade inicial.
    # Com b_enc < 0, a pre-ativacao (x@W + b) e negativa para a maioria
    # dos neuronios, e ReLU os desliga. Apenas neuronios com forte
    # sinal positivo de x@W superam o bias negativo e ficam ativos.
    b_enc = np.full(hidden_dim, -0.5)
    b_dec = np.zeros(input_dim)

    return SAEWeights(
        W_enc=W_enc, b_enc=b_enc,
        W_dec=W_dec, b_dec=b_dec,
        tied=tied,
    )


def forward(
    X: np.ndarray,
    weights: SAEWeights,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Forward pass do SAE.

    x -> z_pre = x @ W_enc + b_enc
      -> z = ReLU(z_pre)           [ativacoes esparsas]
      -> x_hat = z @ W_dec + b_dec [reconstrucao]

    Args:
        X: input (N, input_dim)
        weights: pesos do SAE

    Returns:
        (x_hat, z, z_pre) — reconstrucao, ativacoes, pre-ativacoes
    """
    # Encoder
    z_pre = X @ weights.W_enc + weights.b_enc   # (N, hidden_dim)
    z = _relu(z_pre)                             # (N, hidden_dim)

    # Decoder (com tied weights, W_dec = W_enc^T)
    if weights.tied:
        W_dec = weights.W_enc.T
    else:
        W_dec = weights.W_dec
    x_hat = z @ W_dec + weights.b_dec            # (N, input_dim)

    return x_hat, z, z_pre


def compute_loss(
    X: np.ndarray,
    x_hat: np.ndarray,
    z: np.ndarray,
    sparsity_penalty: float,
) -> tuple[float, float, float]:
    """Calcula loss total = MSE + lambda * L1.

    MSE penaliza reconstrucao ruim — sem ela, o SAE ignora o input.
    L1 penaliza ativacoes densas — sem ela, o SAE usa todos os neuronios
    e nao descobre features interpretaveis.

    O balanco entre MSE e L1 e o hiperparametro mais critico:
      - lambda alto demais: ativacoes todas zero, reconstrucao ruim
      - lambda baixo demais: sem esparsidade, features nao interpretaveis

    Args:
        X: input original (N, D)
        x_hat: reconstrucao (N, D)
        z: ativacoes esparsas (N, H)
        sparsity_penalty: coeficiente lambda da L1

    Returns:
        (total_loss, mse_loss, l1_loss) — valores escalares
    """
    n = X.shape[0]
    mse = np.mean((X - x_hat) ** 2)
    l1 = np.mean(np.abs(z))
    total = mse + sparsity_penalty * l1
    return float(total), float(mse), float(l1)


def train(
    X: np.ndarray,
    weights: SAEWeights,
    sparsity_penalty: float = 0.05,
    learning_rate: float = 5e-3,
    n_epochs: int = 2000,
    target_sparsity: float = 0.15,
    verbose: bool = True,
    print_every: int = 100,
) -> SAEWeights:
    """Treina o SAE com SGD + ajuste adaptativo de bias para esparsidade.

    Algoritmo por epoca:
      1. Forward pass: calcula reconstrucao e ativacoes
      2. Backward pass: gradientes de MSE + L1 via backpropagation
      3. SGD step: atualiza pesos na direcao do gradiente
      4. Bias adaptativo: se esparsidade > alvo, empurra biases para baixo

    A esparsidade e controlada por dois mecanismos complementares:
      - L1 penalty no loss: pressao continua sobre ativacoes
      - Bias adaptativo: ajuste direto quando esparsidade excede alvo

    Args:
        X: dados de treinamento (N, input_dim), ja normalizados
        weights: pesos iniciais do SAE
        sparsity_penalty: coeficiente L1 (lambda)
        learning_rate: taxa de aprendizado para SGD
        n_epochs: numero de epocas de treinamento
        target_sparsity: fracao alvo de neuronios ativos (~0.10-0.20)
        verbose: se True, imprime metricas durante treinamento
        print_every: frequencia de impressao (a cada N epocas)

    Returns:
        SAEWeights treinados (mesmo objeto, modificado in-place)
    """
    n = X.shape[0]
    d = X.shape[1]

    for epoch in range(n_epochs):
        # --- Forward ---
        x_hat, z, z_pre = forward(X, weights)

        # --- Loss ---
        total_loss, mse_loss, l1_loss = compute_loss(
            X, x_hat, z, sparsity_penalty
        )

        # Fracao de neuronios ativos (esparsidade real)
        sparsity = float(np.mean(z > 0))

        # Registra historico
        weights.loss_history.append(total_loss)
        weights.mse_history.append(mse_loss)
        weights.l1_history.append(l1_loss)
        weights.sparsity_history.append(sparsity)

        # --- Backward (gradientes analiticos) ---
        # Gradiente do MSE: d(MSE)/d(x_hat) = 2*(x_hat - X) / (N*D)
        d_xhat = 2.0 * (x_hat - X) / (n * d)

        # Decoder weights para backprop
        if weights.tied:
            W_dec = weights.W_enc.T
        else:
            W_dec = weights.W_dec

        # Gradiente da reconstrucao em relacao a z
        d_z_mse = d_xhat @ W_dec.T

        # Gradiente L1 nas ativacoes: d(mean|z|)/d(z) = sign(z) / N
        # Escala compativel com MSE (ambos divididos por N)
        d_z_l1 = sparsity_penalty * np.sign(z) / n

        # Gradiente total em z, propagado pela ReLU
        d_zpre = (d_z_mse + d_z_l1) * _relu_grad(z_pre)

        # Gradientes dos pesos do encoder
        d_W_enc = X.T @ d_zpre / n
        d_b_enc = np.mean(d_zpre, axis=0)

        # Gradientes do decoder
        d_W_dec = z.T @ d_xhat / n
        d_b_dec = np.mean(d_xhat, axis=0)

        # --- SGD step ---
        weights.W_enc -= learning_rate * d_W_enc
        weights.b_enc -= learning_rate * d_b_enc

        if weights.tied:
            # Aplica gradiente do decoder no encoder (pesos atados)
            weights.W_enc -= learning_rate * d_W_dec.T
            weights.W_dec = weights.W_enc.T.copy()
        else:
            weights.W_dec -= learning_rate * d_W_dec

        weights.b_dec -= learning_rate * d_b_dec

        # --- Ajuste adaptativo de bias para esparsidade-alvo ---
        # Se esparsidade esta acima do alvo (muitos neuronios ativos),
        # empurra biases para baixo. Se abaixo, deixa o gradiente
        # do MSE recuperar neuronios naturalmente.
        if sparsity > target_sparsity:
            # Quantidade proporcional ao excesso de ativacao
            excess = sparsity - target_sparsity
            bias_push = learning_rate * excess * 2.0
            weights.b_enc -= bias_push

        # --- Log ---
        if verbose and (epoch % print_every == 0 or epoch == n_epochs - 1):
            print(
                f"  Epoca {epoch:4d}/{n_epochs} | "
                f"Loss: {total_loss:.6f} | "
                f"MSE: {mse_loss:.6f} | "
                f"L1: {l1_loss:.6f} | "
                f"Esparsidade: {sparsity:.1%}"
            )

    return weights


def get_activations(
    X: np.ndarray,
    weights: SAEWeights,
) -> np.ndarray:
    """Extrai ativacoes esparsas do encoder para cada input.

    Apos treinamento, as ativacoes sao o produto principal:
    cada dimensao corresponde a uma "feature" descoberta pelo SAE.

    Args:
        X: inputs (N, input_dim)
        weights: pesos treinados

    Returns:
        Ativacoes (N, hidden_dim) — maioria dos valores sera zero
    """
    _, z, _ = forward(X, weights)
    return z


def reconstruction_quality(
    X: np.ndarray,
    weights: SAEWeights,
) -> dict[str, float]:
    """Avalia qualidade da reconstrucao do SAE.

    Metricas:
      - MSE: erro medio quadratico
      - R2: coeficiente de determinacao (1.0 = perfeito)
      - Sparsity: fracao de neuronios ativos por amostra

    Args:
        X: inputs (N, input_dim)
        weights: pesos treinados

    Returns:
        Dict com metricas de qualidade
    """
    x_hat, z, _ = forward(X, weights)
    mse = float(np.mean((X - x_hat) ** 2))

    # R^2 — variancia explicada
    ss_res = np.sum((X - x_hat) ** 2)
    ss_tot = np.sum((X - np.mean(X, axis=0)) ** 2)
    r2 = 1.0 - ss_res / max(ss_tot, 1e-12)

    # Esparsidade: fracao media de neuronios ativos
    sparsity = float(np.mean(z > 0))

    # Neuronios nunca ativos (mortos) — indicam overcomplete demais
    dead_neurons = int(np.sum(np.all(z == 0, axis=0)))
    total_neurons = z.shape[1]

    return {
        "mse": round(mse, 6),
        "r2": round(float(r2), 4),
        "sparsity": round(sparsity, 4),
        "dead_neurons": dead_neurons,
        "alive_neurons": total_neurons - dead_neurons,
        "total_neurons": total_neurons,
    }
