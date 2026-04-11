"""Ambiente de RL para otimizacao de modificacoes quimicas do ASO.

O ambiente modela o ASO MRL-ASO-001 como um vetor de posicoes, cada uma
com tres atributos: base (A/C/G/T), flag LNA (0/1), tipo de backbone
(PS=1/PO=0). O agente aplica modificacoes incrementais e recebe
recompensa baseada em metricas termodinamicas do aso_math.thermo.

Design gapmer: o ASO real usa flancos LNA (5+5) com gap DNA central (15),
backbone fosforotioato completo. O agente pode explorar variacoes dentro
dos limites biologicos definidos no RLConfig.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from aso_math.config import ASO_SEQUENCE, BASES
from aso_math.thermo import (
    compute_dg,
    compute_hairpin_dg,
    compute_self_dimer_dg,
    compute_tm,
    gc_content,
)
from .config import (
    BACKBONE_PO,
    BACKBONE_PS,
    MODIFICATION_TYPES,
    RLConfig,
)


# ---------------------------------------------------------------------------
# Representacao do estado do ASO
# ---------------------------------------------------------------------------

BASE_TO_IDX: dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}
IDX_TO_BASE: dict[int, str] = {v: k for k, v in BASE_TO_IDX.items()}
N_BASES: int = 4


@dataclass
class ASOState:
    """Estado mutavel do ASO durante um episodio de RL.

    Atributos:
        sequence: lista de caracteres da sequencia (len = n_positions)
        lna_mask: 1 se a posicao tem LNA, 0 se DNA (len = n_positions)
        backbone: 1 se PS, 0 se PO na posicao (len = n_positions)
    """

    sequence: list[str]
    lna_mask: list[int]
    backbone: list[int]

    @property
    def seq_str(self) -> str:
        """Sequencia como string para calculo termodinamico."""
        return "".join(self.sequence)

    @property
    def n_positions(self) -> int:
        return len(self.sequence)

    def copy(self) -> ASOState:
        """Retorna copia profunda do estado."""
        return ASOState(
            sequence=self.sequence.copy(),
            lna_mask=self.lna_mask.copy(),
            backbone=self.backbone.copy(),
        )


def create_initial_state(sequence: str = ASO_SEQUENCE) -> ASOState:
    """Cria estado inicial a partir da sequencia MRL-ASO-001.

    Design gapmer padrao: 5 LNA + 15 DNA + 5 LNA, backbone PS completo.
    Este e o ponto de partida do agente — o ASO referencia do pipeline.
    """
    n = len(sequence)
    seq_list = list(sequence.upper())

    # Gapmer 5-15-5: flancos LNA, centro DNA
    lna_mask = [0] * n
    for i in range(5):
        lna_mask[i] = 1          # flanco 5' LNA
        lna_mask[n - 1 - i] = 1  # flanco 3' LNA

    # Backbone PS completo (padrao para ASOs terapeuticos)
    backbone = [BACKBONE_PS] * n

    return ASOState(sequence=seq_list, lna_mask=lna_mask, backbone=backbone)


# ---------------------------------------------------------------------------
# Ambiente de RL
# ---------------------------------------------------------------------------


class ASODesignEnv:
    """Ambiente de aprendizado por reforco para otimizacao de ASO.

    Espaco de observacao:
        Vetor concatenado de:
        - One-hot da sequencia (n_positions * 4)
        - LNA mask (n_positions)
        - Backbone type (n_positions)
        Total: n_positions * 6

    Espaco de acoes:
        (posicao, tipo_modificacao) -> indice linear
        posicao: 0..n_positions-1
        tipo_modificacao:
            0 = change_base (cicla para proxima base)
            1 = toggle_lna
            2 = toggle_backbone
        Total: n_positions * 3 acoes discretas

    Recompensa:
        Combinacao ponderada de metricas termodinamicas.
        Valores normalizados para [0, 1] e somados com pesos do config.
    """

    def __init__(self, config: RLConfig | None = None) -> None:
        self.config = config or RLConfig(
            module_slug="08_rl_ppo",
            module_name="ASO RL Optimization",
        )
        self.initial_state = create_initial_state()
        self.state: ASOState = self.initial_state.copy()
        self.n_positions: int = self.state.n_positions
        self.n_action_types: int = len(MODIFICATION_TYPES)
        self.n_actions: int = self.n_positions * self.n_action_types
        self.obs_dim: int = self.n_positions * 6  # 4 one-hot + lna + backbone
        self.step_count: int = 0
        self.prev_reward: float = 0.0

        # Calcula recompensa do estado inicial (baseline)
        self.baseline_reward: float = self._compute_reward()

    def reset(self) -> np.ndarray:
        """Reinicia o ambiente para o estado inicial do MRL-ASO-001.

        Returns:
            Vetor de observacao do estado inicial.
        """
        self.state = self.initial_state.copy()
        self.step_count = 0
        self.prev_reward = self.baseline_reward
        return self.get_observation()

    def step(self, action: int) -> tuple[np.ndarray, float, bool, dict[str, Any]]:
        """Aplica uma acao ao ASO e retorna (obs, reward, done, info).

        Args:
            action: Indice linear da acao (0..n_actions-1).
                    Decodificado como (posicao, tipo_modificacao).

        Returns:
            Tupla (observacao, recompensa, episodio_terminado, info_dict).
        """
        # Decodificar acao
        position = action // self.n_action_types
        mod_type = action % self.n_action_types

        # Clampar posicao para seguranca
        position = min(position, self.n_positions - 1)

        # Aplicar modificacao
        info: dict[str, Any] = {"position": position, "modification": MODIFICATION_TYPES[mod_type]}
        self._apply_action(position, mod_type)

        # Calcular recompensa
        reward = self._compute_reward()

        # Recompensa incremental: diferenca em relacao ao passo anterior
        # incentiva melhorias graduais
        reward_delta = reward - self.prev_reward
        self.prev_reward = reward

        self.step_count += 1

        # Condicao de termino
        done = self.step_count >= self.config.max_steps_per_episode

        # Info adicional para diagnostico
        info["reward_absolute"] = reward
        info["reward_delta"] = reward_delta
        info["step"] = self.step_count
        info["sequence"] = self.state.seq_str

        return self.get_observation(), reward_delta, done, info

    def get_observation(self) -> np.ndarray:
        """Converte estado atual para vetor numerico de observacao.

        Formato: [one_hot_seq (n*4) | lna_mask (n) | backbone (n)]
        """
        obs = np.zeros(self.obs_dim, dtype=np.float32)

        for i, base in enumerate(self.state.sequence):
            base_idx = BASE_TO_IDX.get(base, 0)
            obs[i * N_BASES + base_idx] = 1.0

        offset_lna = self.n_positions * N_BASES
        for i, lna in enumerate(self.state.lna_mask):
            obs[offset_lna + i] = float(lna)

        offset_bb = offset_lna + self.n_positions
        for i, bb in enumerate(self.state.backbone):
            obs[offset_bb + i] = float(bb)

        return obs

    # -------------------------------------------------------------------
    # Acoes
    # -------------------------------------------------------------------

    def _apply_action(self, position: int, mod_type: int) -> None:
        """Aplica modificacao ao ASO na posicao dada.

        mod_type 0 (change_base): cicla para a proxima base (A->C->G->T->A)
        mod_type 1 (toggle_lna): inverte flag LNA na posicao
        mod_type 2 (toggle_backbone): inverte PS<->PO na posicao
        """
        if mod_type == 0:
            # Ciclar base: A -> C -> G -> T -> A
            current_idx = BASE_TO_IDX.get(self.state.sequence[position], 0)
            next_idx = (current_idx + 1) % N_BASES
            self.state.sequence[position] = IDX_TO_BASE[next_idx]

        elif mod_type == 1:
            # Toggle LNA
            self.state.lna_mask[position] = 1 - self.state.lna_mask[position]

        elif mod_type == 2:
            # Toggle backbone PS <-> PO
            self.state.backbone[position] = 1 - self.state.backbone[position]

    # -------------------------------------------------------------------
    # Recompensa
    # -------------------------------------------------------------------

    def _compute_reward(self) -> float:
        """Calcula recompensa multi-objetivo para o estado atual.

        Cada componente e normalizado para [0, 1] e ponderado pelos pesos
        definidos no config. A recompensa total esta em [0, 1].

        Componentes:
        1. dG binding: normalizado entre dg_worst e dg_best
        2. Tm: normalizado entre tm_worst e tm_best
        3. GC content: 1.0 se dentro de [0.40, 0.60], decai linearmente fora
        4. Hairpin: penaliza hairpins estaveis (dG muito negativo)
        5. Self-dimer: penaliza auto-dimeros estaveis
        6. Length penalty: 1.0 se dentro de [23, 27], 0.0 fora
        """
        seq = self.state.seq_str
        cfg = self.config
        weights = cfg.reward_weights

        # 1. delta_G de hibridizacao (mais negativo = melhor)
        dg = compute_dg(seq)
        r_dg = _normalize(dg, cfg.dg_worst, cfg.dg_best)

        # 2. Temperatura de melting (mais alto = melhor)
        tm = compute_tm(seq)
        r_tm = _normalize(tm, cfg.tm_worst, cfg.tm_best)

        # 3. GC content (40-60% e otimo)
        gc = gc_content(seq)
        r_gc = _gc_reward(gc, cfg.gc_optimal_low, cfg.gc_optimal_high)

        # 4. Hairpin stability (menos negativo = melhor, ou seja, mais proximo de 0)
        hairpin_dg = compute_hairpin_dg(seq)
        # hairpin_dg <= 0; quanto mais proximo de 0, melhor
        # Normalizamos: 0.0 -> reward 1.0, -10.0 -> reward 0.0
        r_hairpin = _normalize(hairpin_dg, -10.0, 0.0)

        # 5. Self-dimer (menos negativo = melhor)
        dimer_dg = compute_self_dimer_dg(seq)
        r_dimer = _normalize(dimer_dg, -15.0, 0.0)

        # 6. Penalidade de comprimento
        n = len(seq)
        r_length = 1.0 if cfg.length_min <= n <= cfg.length_max else 0.0

        # Recompensa total ponderada
        total = (
            weights["dg_binding"] * r_dg
            + weights["tm"] * r_tm
            + weights["gc_content"] * r_gc
            + weights["hairpin"] * r_hairpin
            + weights["self_dimer"] * r_dimer
            + weights["length_penalty"] * r_length
        )

        return total

    def get_metrics(self) -> dict[str, float]:
        """Retorna metricas termodinamicas detalhadas do estado atual.

        Util para diagnostico e comparacao de variantes.
        """
        seq = self.state.seq_str
        return {
            "sequence": seq,
            "length": len(seq),
            "dg_binding": compute_dg(seq),
            "tm": compute_tm(seq),
            "gc_content": round(gc_content(seq), 4),
            "hairpin_dg": compute_hairpin_dg(seq),
            "self_dimer_dg": compute_self_dimer_dg(seq),
            "n_lna": sum(self.state.lna_mask),
            "n_ps": sum(self.state.backbone),
            "reward": self._compute_reward(),
        }


# ---------------------------------------------------------------------------
# Funcoes auxiliares de normalizacao
# ---------------------------------------------------------------------------


def _normalize(value: float, worst: float, best: float) -> float:
    """Normaliza valor para [0, 1] entre worst e best.

    Clampado: valores alem de best retornam 1.0,
    valores alem de worst retornam 0.0.
    """
    if abs(best - worst) < 1e-12:
        return 0.5
    normalized = (value - worst) / (best - worst)
    return max(0.0, min(1.0, normalized))


def _gc_reward(gc: float, low: float, high: float) -> float:
    """Recompensa para GC content com zona otima.

    Retorna 1.0 se gc esta entre [low, high].
    Decai linearmente para 0.0 fora da zona otima.
    A tolerancia e de 0.20 (20%) para cada lado.
    """
    if low <= gc <= high:
        return 1.0
    tolerance = 0.20
    if gc < low:
        return max(0.0, 1.0 - (low - gc) / tolerance)
    # gc > high
    return max(0.0, 1.0 - (gc - high) / tolerance)
