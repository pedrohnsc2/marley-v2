"""Processo de difusao discreta para sequencias biologicas.

Implementa o forward process (corrupcao gradual) e o reverse process
(denoising iterativo) para gerar novas sequencias a partir de ruido.

O modelo segue a abordagem D3PM (Discrete Denoising Diffusion Probabilistic
Models) simplificada: a cada passo t, cada posicao da sequencia tem
probabilidade beta_t de ser substituida por um token aleatorio do alfabeto.

Ref: Austin J et al. (2021) Structured Denoising Diffusion Models — NeurIPS
"""

from __future__ import annotations

import torch
import torch.nn.functional as F


def linear_beta_schedule(
    n_steps: int,
    beta_start: float = 0.001,
    beta_end: float = 0.1,
) -> torch.Tensor:
    """Cria schedule linear de betas para o processo de difusao.

    beta_t cresce linearmente de beta_start a beta_end ao longo de T passos.
    Valores maiores de beta causam mais corrupcao por passo.

    Args:
        n_steps: Numero total de passos T.
        beta_start: Beta no passo t=0.
        beta_end: Beta no passo t=T-1.

    Returns:
        Tensor de shape (n_steps,) com os betas.
    """
    return torch.linspace(beta_start, beta_end, n_steps)


def compute_alpha_bar(betas: torch.Tensor) -> torch.Tensor:
    """Calcula alpha_bar_t = produto cumulativo de (1 - beta_i) para i=0..t.

    alpha_bar_t representa a probabilidade de um token NAO ter sido
    corrompido ate o passo t. Quanto menor, mais corrompida a sequencia.

    Args:
        betas: Schedule de betas, shape (T,).

    Returns:
        Tensor de shape (T,) com alpha_bar para cada passo.
    """
    alphas = 1.0 - betas
    return torch.cumprod(alphas, dim=0)


def forward_diffusion(
    x_0: torch.Tensor,
    t: torch.Tensor,
    alpha_bar: torch.Tensor,
    vocab_size: int,
) -> torch.Tensor:
    """Aplica corrupcao ao passo t (forward process).

    Para cada posicao na sequencia, com probabilidade (1 - alpha_bar_t)
    o token e substituido por um token aleatorio uniforme do alfabeto.

    Este e o processo de corrupcao "absorbing" simplificado: tokens
    corrompidos sao substituidos por tokens aleatorios (nao por um token
    especial de mascara). Isso forca o modelo a aprender a distribuicao
    completa do alfabeto.

    Args:
        x_0: Sequencias originais (indices), shape (B, L).
        t: Passo de difusao para cada amostra, shape (B,).
        alpha_bar: Probabilidades cumulativas, shape (T,).
        vocab_size: Tamanho do vocabulario (incluindo MASK se usado).

    Returns:
        Sequencias corrompidas x_t, shape (B, L).
    """
    batch_size, seq_len = x_0.shape
    device = x_0.device

    # alpha_bar para cada amostra do batch: shape (B, 1)
    ab_t = alpha_bar[t].unsqueeze(1)  # (B, 1)

    # Mascara: True onde o token sera mantido (nao corrompido)
    keep_mask = torch.rand(batch_size, seq_len, device=device) < ab_t

    # Tokens aleatorios para posicoes corrompidas
    random_tokens = torch.randint(0, vocab_size, (batch_size, seq_len), device=device)

    # Combina: mantém original onde keep_mask=True, aleatorio caso contrario
    x_t = torch.where(keep_mask, x_0, random_tokens)

    return x_t


@torch.no_grad()
def reverse_diffusion_sample(
    model: torch.nn.Module,
    seq_len: int,
    vocab_size: int,
    betas: torch.Tensor,
    alpha_bar: torch.Tensor,
    n_samples: int = 1,
    device: str = "cpu",
) -> torch.Tensor:
    """Gera sequencias via denoising iterativo (reverse process).

    Comeca com sequencias de tokens puramente aleatorios e, a cada passo
    (de T-1 ate 0), usa a rede neural para prever os tokens originais.
    A cada passo, a predicao e misturada com ruido residual proporcional
    ao nível de corrupcao restante.

    Estrategia: a cada passo t, o modelo prediz logits para cada posicao.
    Amostramos da distribuicao predita (com temperatura implícita via softmax)
    e mantemos ou re-ruidos com base em alpha_bar_{t-1}.

    Args:
        model: Rede de denoising (recebe x_t, t -> logits).
        seq_len: Comprimento da sequencia a gerar.
        vocab_size: Tamanho do vocabulario.
        betas: Schedule de betas, shape (T,).
        alpha_bar: alpha_bar cumulativo, shape (T,).
        n_samples: Numero de sequencias a gerar.
        device: Dispositivo de computacao.

    Returns:
        Sequencias geradas, shape (n_samples, seq_len) como indices inteiros.
    """
    model.eval()
    n_steps = len(betas)

    # Comeca com ruido puro: tokens aleatorios
    x_t = torch.randint(0, vocab_size, (n_samples, seq_len), device=device)

    for step in reversed(range(n_steps)):
        # Tensor de timestep para o batch inteiro
        t_tensor = torch.full((n_samples,), step, device=device, dtype=torch.long)

        # Predicao do modelo: logits para cada posicao
        logits = model(x_t, t_tensor)  # (B, L, vocab_size)

        # Amostra da distribuicao predita
        probs = F.softmax(logits, dim=-1)  # (B, L, vocab_size)
        flat_probs = probs.reshape(-1, vocab_size)
        predicted = torch.multinomial(flat_probs, num_samples=1).squeeze(-1)
        predicted = predicted.reshape(n_samples, seq_len)

        if step > 0:
            # Re-ruidar parcialmente: com probabilidade beta_{t-1} manter ruido
            # Isso evita colapso para sequencias deterministicas
            ab_prev = alpha_bar[step - 1]
            keep_predicted = torch.rand(n_samples, seq_len, device=device) < ab_prev
            random_tokens = torch.randint(0, vocab_size, (n_samples, seq_len), device=device)
            x_t = torch.where(keep_predicted, predicted, random_tokens)
        else:
            # Ultimo passo: aceita a predicao diretamente
            x_t = predicted

    return x_t


def tokens_to_sequences(
    token_indices: torch.Tensor,
    alphabet: str,
) -> list[str]:
    """Converte tensor de indices em lista de strings de sequencia.

    Indices fora do range do alfabeto (ex: token de mascara) sao
    mapeados para 'X' (desconhecido).

    Args:
        token_indices: Tensor de shape (B, L) com indices inteiros.
        alphabet: String do alfabeto (ex: "ACGT").

    Returns:
        Lista de B strings de sequencia.
    """
    sequences = []
    indices_np = token_indices.cpu().numpy()

    for row in indices_np:
        chars = []
        for idx in row:
            if 0 <= idx < len(alphabet):
                chars.append(alphabet[idx])
            else:
                chars.append("X")  # Token fora do alfabeto
        sequences.append("".join(chars))

    return sequences


def sequences_to_tokens(
    sequences: list[str],
    alphabet: str,
    device: str = "cpu",
) -> torch.Tensor:
    """Converte lista de sequencias em tensor de indices.

    Caracteres nao encontrados no alfabeto recebem indice len(alphabet),
    que corresponde ao token de mascara.

    Args:
        sequences: Lista de strings de sequencia.
        alphabet: String do alfabeto.
        device: Dispositivo de destino.

    Returns:
        Tensor de shape (B, max_len) com indices inteiros.
    """
    char_to_idx = {c: i for i, c in enumerate(alphabet)}
    mask_idx = len(alphabet)

    # Encontra comprimento maximo
    max_len = max(len(s) for s in sequences)

    tokens = []
    for seq in sequences:
        row = [char_to_idx.get(c.upper(), mask_idx) for c in seq]
        # Padding com token de mascara se necessario
        row.extend([mask_idx] * (max_len - len(row)))
        tokens.append(row)

    return torch.tensor(tokens, dtype=torch.long, device=device)
