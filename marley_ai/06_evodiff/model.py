"""Rede neural de denoising para difusao discreta.

Transformer encoder pequeno (2 camadas, 128 dim) que aprende a reconstruir
sequencias biologicas a partir de versoes corrompidas. Recebe a sequencia
corrompida x_t e o timestep t, e prediz logits sobre o alfabeto para
cada posicao.

Arquitetura:
    Input: one-hot(x_t) + timestep_embedding(t)
    -> Embedding layer (vocab_size -> d_model)
    -> Positional encoding (sinusoidal)
    -> Timestep embedding (MLP)
    -> 2x Transformer Encoder Layer
    -> Linear head (d_model -> vocab_size)
    -> Output: logits (B, L, vocab_size)

Dimensionado para treinar em MPS (Apple M3 Pro) em minutos com
sequencias curtas (9-27 posicoes).
"""

from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class SinusoidalPositionEncoding(nn.Module):
    """Codificacao posicional sinusoidal (Vaswani et al., 2017).

    Adiciona informacao de posicao a cada token da sequencia via funcoes
    seno/cosseno de diferentes frequencias. Permite ao transformer
    distinguir posicoes sem parametros aprendiveis.
    """

    def __init__(self, d_model: int, max_len: int = 512) -> None:
        super().__init__()
        # Pre-calcula a tabela de posicoes
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2, dtype=torch.float)
            * (-math.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        # Registra como buffer (nao e parametro treinavel)
        self.register_buffer("pe", pe.unsqueeze(0))  # (1, max_len, d_model)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Adiciona positional encoding ao input.

        Args:
            x: Tensor de shape (B, L, d_model).

        Returns:
            Tensor com positional encoding somado, mesma shape.
        """
        return x + self.pe[:, :x.size(1), :]


class TimestepEmbedding(nn.Module):
    """Embedding do timestep de difusao via sinusoidal + MLP.

    O timestep t e convertido em um vetor de dimensao d_model usando
    codificacao sinusoidal (mesma tecnica de positional encoding)
    seguida de uma MLP de 2 camadas para projecao nao-linear.

    Isso permite ao modelo condicionar sua predicao no nivel de
    corrupcao da sequencia.
    """

    def __init__(self, d_model: int) -> None:
        super().__init__()
        self.d_model = d_model
        # MLP: sinusoidal_dim -> 4*d_model -> d_model
        self.mlp = nn.Sequential(
            nn.Linear(d_model, d_model * 4),
            nn.GELU(),
            nn.Linear(d_model * 4, d_model),
        )

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        """Converte timestep inteiro em embedding continuo.

        Args:
            t: Timesteps inteiros, shape (B,).

        Returns:
            Embeddings do timestep, shape (B, d_model).
        """
        # Codificacao sinusoidal do timestep
        half_dim = self.d_model // 2
        emb_scale = math.log(10000.0) / (half_dim - 1)
        emb = torch.exp(
            torch.arange(half_dim, device=t.device, dtype=torch.float) * -emb_scale
        )
        emb = t.float().unsqueeze(1) * emb.unsqueeze(0)  # (B, half_dim)
        emb = torch.cat([torch.sin(emb), torch.cos(emb)], dim=-1)  # (B, d_model)

        return self.mlp(emb)


class SequenceDenoiser(nn.Module):
    """Transformer denoiser para difusao discreta em sequencias biologicas.

    Recebe uma sequencia corrompida (como indices inteiros) e um timestep,
    e prediz os logits do token original em cada posicao.

    A perda de treinamento e cross-entropy entre os logits preditos
    e os tokens originais (nao corrompidos).
    """

    def __init__(
        self,
        vocab_size: int,
        d_model: int = 128,
        n_heads: int = 4,
        n_layers: int = 2,
        dropout: float = 0.1,
        max_seq_len: int = 64,
    ) -> None:
        """Inicializa o denoiser.

        Args:
            vocab_size: Tamanho do vocabulario (alphabet + mask token).
            d_model: Dimensao do modelo (embedding e hidden).
            n_heads: Numero de cabecas de atencao.
            n_layers: Numero de camadas do transformer encoder.
            dropout: Taxa de dropout.
            max_seq_len: Comprimento maximo de sequencia suportado.
        """
        super().__init__()
        self.vocab_size = vocab_size
        self.d_model = d_model

        # Embedding de tokens (inclui token de mascara)
        self.token_embedding = nn.Embedding(vocab_size, d_model)

        # Positional encoding sinusoidal
        self.pos_encoding = SinusoidalPositionEncoding(d_model, max_len=max_seq_len)

        # Timestep embedding
        self.time_embedding = TimestepEmbedding(d_model)

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer,
            num_layers=n_layers,
        )

        # Cabeca de predicao: d_model -> vocab_size
        self.output_head = nn.Linear(d_model, vocab_size)

        # Inicializacao dos pesos (Xavier uniform)
        self._init_weights()

    def _init_weights(self) -> None:
        """Inicializa pesos com Xavier uniform para estabilidade."""
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(
        self,
        x_t: torch.Tensor,
        t: torch.Tensor,
    ) -> torch.Tensor:
        """Forward pass: prediz logits do token original a partir do corrompido.

        Args:
            x_t: Sequencia corrompida (indices), shape (B, L).
            t: Timestep de difusao, shape (B,).

        Returns:
            Logits sobre o vocabulario para cada posicao, shape (B, L, vocab_size).
        """
        # Embedding dos tokens: (B, L) -> (B, L, d_model)
        h = self.token_embedding(x_t)

        # Adiciona positional encoding
        h = self.pos_encoding(h)

        # Embedding do timestep: (B,) -> (B, d_model) -> (B, 1, d_model)
        t_emb = self.time_embedding(t).unsqueeze(1)

        # Soma o timestep embedding a todas as posicoes
        # (broadcast: (B, L, d_model) + (B, 1, d_model))
        h = h + t_emb

        # Transformer encoder
        h = self.transformer(h)  # (B, L, d_model)

        # Cabeca de predicao
        logits = self.output_head(h)  # (B, L, vocab_size)

        return logits


def compute_loss(
    model: SequenceDenoiser,
    x_0: torch.Tensor,
    alpha_bar: torch.Tensor,
    vocab_size: int,
) -> torch.Tensor:
    """Calcula a perda de treinamento (cross-entropy).

    1. Amostra timestep t uniforme para cada sequencia no batch
    2. Corrompe x_0 -> x_t usando forward diffusion
    3. Prediz logits com o modelo
    4. Cross-entropy entre logits e x_0 (tokens originais)

    Args:
        model: Rede de denoising.
        x_0: Sequencias originais (indices), shape (B, L).
        alpha_bar: Probabilidades cumulativas, shape (T,).
        vocab_size: Tamanho do vocabulario.

    Returns:
        Escalar: perda media cross-entropy.
    """
    # Importa aqui para evitar circular import
    from .diffusion import forward_diffusion

    batch_size = x_0.shape[0]
    n_steps = len(alpha_bar)
    device = x_0.device

    # Amostra timestep uniforme para cada elemento do batch
    t = torch.randint(0, n_steps, (batch_size,), device=device)

    # Corrompe a sequencia
    x_t = forward_diffusion(x_0, t, alpha_bar, vocab_size)

    # Predicao do modelo
    logits = model(x_t, t)  # (B, L, vocab_size)

    # Cross-entropy: compara predicao com sequencia original
    # Reshape para (B*L, vocab_size) e (B*L,)
    loss = F.cross_entropy(
        logits.reshape(-1, vocab_size),
        x_0.reshape(-1),
    )

    return loss
