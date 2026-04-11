"""Agente REINFORCE (policy gradient) para otimizacao de ASO.

Implementacao pura em numpy — sem dependencia de PyTorch ou
stable-baselines3. O agente usa uma politica linear com softmax
sobre o espaco de acoes discretas.

Algoritmo REINFORCE (Williams, 1992):
1. Coleta trajetoria completa (episodio) seguindo a politica atual
2. Calcula retornos descontados G_t = sum_{k=0}^{T-t} gamma^k * r_{t+k}
3. Subtrai baseline (media movel) para reducao de variancia
4. Atualiza pesos: theta += alpha * sum_t[ (G_t - b) * grad log pi(a_t|s_t) ]

A politica linear e suficiente para este problema porque:
- O espaco de acoes e moderado (~75 acoes para ASO de 25 nt)
- A funcao de recompensa e relativamente suave
- O objetivo e exploracao eficiente, nao performance de estado da arte

Ref: Williams RJ (1992) Simple statistical gradient-following algorithms
     for connectionist reinforcement learning. Machine Learning 8(3-4):229-256
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import RLConfig


# ---------------------------------------------------------------------------
# Buffer de trajetoria
# ---------------------------------------------------------------------------

@dataclass
class TrajectoryBuffer:
    """Armazena uma trajetoria completa (episodio) para calculo de gradiente.

    Cada entrada corresponde a um passo: (estado, acao, recompensa).
    Ao final do episodio, calcula retornos descontados e vantagens.
    """

    states: list[np.ndarray] = field(default_factory=list)
    actions: list[int] = field(default_factory=list)
    rewards: list[float] = field(default_factory=list)

    def store(self, state: np.ndarray, action: int, reward: float) -> None:
        """Armazena uma transicao (s, a, r)."""
        self.states.append(state.copy())
        self.actions.append(action)
        self.rewards.append(reward)

    def compute_returns(self, gamma: float) -> np.ndarray:
        """Calcula retornos descontados G_t para cada passo.

        G_t = r_t + gamma * r_{t+1} + gamma^2 * r_{t+2} + ...

        Calculado de tras para frente para eficiencia O(T).
        """
        n = len(self.rewards)
        returns = np.zeros(n, dtype=np.float64)
        running_return = 0.0

        for t in reversed(range(n)):
            running_return = self.rewards[t] + gamma * running_return
            returns[t] = running_return

        return returns

    def clear(self) -> None:
        """Limpa buffer para o proximo episodio."""
        self.states.clear()
        self.actions.clear()
        self.rewards.clear()

    def __len__(self) -> int:
        return len(self.rewards)


# ---------------------------------------------------------------------------
# Agente REINFORCE
# ---------------------------------------------------------------------------

class REINFORCEAgent:
    """Agente policy gradient com politica linear softmax.

    Politica: pi(a|s) = softmax(W @ s + b)[a]

    Onde:
        W: matriz de pesos (n_actions x obs_dim)
        b: vetor de bias (n_actions)
        s: vetor de observacao (obs_dim)

    O gradiente de log pi(a|s) para politica softmax e:
        grad_W log pi(a|s) = (e_a - pi(.|s)) @ s^T
        grad_b log pi(a|s) = (e_a - pi(.|s))

    Onde e_a e o vetor one-hot da acao selecionada.
    """

    def __init__(
        self,
        obs_dim: int,
        n_actions: int,
        config: RLConfig | None = None,
    ) -> None:
        self.obs_dim = obs_dim
        self.n_actions = n_actions
        self.config = config or RLConfig(
            module_slug="08_rl_ppo",
            module_name="ASO RL Optimization",
        )

        # Semente para reproducibilidade
        self.rng = np.random.default_rng(self.config.seed)

        # Pesos da politica linear — inicializacao pequena para explorar
        # uniformemente no inicio (softmax quase uniforme)
        scale = 0.01
        self.W = self.rng.normal(0, scale, size=(n_actions, obs_dim)).astype(np.float64)
        self.b = np.zeros(n_actions, dtype=np.float64)

        # Baseline de media movel (reducao de variancia)
        self.baseline: float = 0.0
        self.baseline_decay: float = self.config.baseline_decay

        # Buffer de trajetoria
        self.buffer = TrajectoryBuffer()

        # Historico de treinamento
        self.episode_rewards: list[float] = []
        self.episode_lengths: list[int] = []

    def select_action(self, obs: np.ndarray, greedy: bool = False) -> int:
        """Seleciona acao segundo a politica atual.

        Args:
            obs: Vetor de observacao do ambiente.
            greedy: Se True, seleciona a acao com maior probabilidade
                    (sem exploracao). Usado na avaliacao final.

        Returns:
            Indice da acao selecionada.
        """
        probs = self._compute_policy(obs)

        if greedy:
            return int(np.argmax(probs))

        # Amostragem estocastica da politica
        action = self.rng.choice(self.n_actions, p=probs)
        return int(action)

    def update(self) -> float:
        """Atualiza pesos da politica usando o buffer de trajetoria.

        Implementa REINFORCE com baseline:
            theta += alpha * sum_t[ (G_t - b) * grad log pi(a_t|s_t) ]

        A baseline e a media movel exponencial dos retornos,
        que reduz variancia sem introduzir vies.

        Returns:
            Perda media (para monitoramento).
        """
        if len(self.buffer) == 0:
            return 0.0

        # Calcula retornos descontados
        returns = self.buffer.compute_returns(self.config.gamma)

        # Normaliza retornos para estabilidade numerica
        # (tecnica padrao em policy gradient)
        returns_std = returns.std()
        if returns_std > 1e-8:
            returns_normalized = (returns - returns.mean()) / returns_std
        else:
            returns_normalized = returns - returns.mean()

        # Atualiza baseline (media movel exponencial do retorno do episodio)
        episode_return = returns[0]  # G_0 = retorno total do episodio
        self.baseline = (
            self.baseline_decay * self.baseline
            + (1 - self.baseline_decay) * episode_return
        )

        # Gradiente da politica
        lr = self.config.learning_rate
        total_loss = 0.0

        for t in range(len(self.buffer)):
            state = self.buffer.states[t]
            action = self.buffer.actions[t]
            advantage = returns_normalized[t]

            # Calcula probabilidades da politica
            probs = self._compute_policy(state)

            # Gradiente de log pi(a|s) para softmax linear:
            # grad_W = (e_a - pi) @ s^T
            # grad_b = (e_a - pi)
            one_hot_a = np.zeros(self.n_actions, dtype=np.float64)
            one_hot_a[action] = 1.0
            diff = one_hot_a - probs  # (n_actions,)

            # Atualiza pesos: theta += lr * advantage * grad log pi
            grad_W = np.outer(diff, state.astype(np.float64))  # (n_actions, obs_dim)
            grad_b = diff  # (n_actions,)

            self.W += lr * advantage * grad_W
            self.b += lr * advantage * grad_b

            # Perda para monitoramento (negativo da log-probabilidade ponderada)
            log_prob = np.log(probs[action] + 1e-10)
            total_loss -= advantage * log_prob

        avg_loss = total_loss / len(self.buffer)

        # Registra estatisticas do episodio
        self.episode_rewards.append(float(episode_return))
        self.episode_lengths.append(len(self.buffer))

        # Limpa buffer
        self.buffer.clear()

        return float(avg_loss)

    def _compute_policy(self, obs: np.ndarray) -> np.ndarray:
        """Calcula probabilidades da politica softmax.

        pi(a|s) = exp(z_a) / sum_a' exp(z_a')
        Onde z = W @ s + b (logits)

        Usa o truque de subtrair max(z) para estabilidade numerica.
        """
        obs_f64 = obs.astype(np.float64)
        logits = self.W @ obs_f64 + self.b

        # Estabilidade numerica: subtrair max para evitar overflow em exp
        logits_shifted = logits - logits.max()
        exp_logits = np.exp(logits_shifted)
        probs = exp_logits / exp_logits.sum()

        # Clampar para evitar probabilidades zero (causa log(0) = -inf)
        probs = np.clip(probs, 1e-8, 1.0)
        probs /= probs.sum()  # renormalizar apos clamp

        return probs

    def get_training_stats(self) -> dict[str, float]:
        """Retorna estatisticas de treinamento para diagnostico."""
        if not self.episode_rewards:
            return {
                "mean_reward": 0.0,
                "std_reward": 0.0,
                "max_reward": 0.0,
                "min_reward": 0.0,
                "baseline": self.baseline,
                "n_episodes": 0,
            }

        rewards = np.array(self.episode_rewards)
        last_50 = rewards[-50:] if len(rewards) >= 50 else rewards
        return {
            "mean_reward": float(last_50.mean()),
            "std_reward": float(last_50.std()),
            "max_reward": float(rewards.max()),
            "min_reward": float(rewards.min()),
            "baseline": self.baseline,
            "n_episodes": len(self.episode_rewards),
            "mean_last_50": float(last_50.mean()),
        }
