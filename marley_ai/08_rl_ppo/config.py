"""Configuracao do modulo 08_rl_ppo — Otimizacao de ASO por REINFORCE.

Reproposto para otimizar modificacoes quimicas do ASO anti-SL RNA
usando aprendizado por reforco com policy gradient (REINFORCE).
O agente aprende a posicionar LNA, escolher backbone (PS/PO) e
ajustar a sequencia para maximizar eficacia terapeutica.

Funcao de recompensa: combinacao ponderada de metricas termodinamicas
calculadas por aso_math.thermo (compute_dg, compute_tm, gc_content,
compute_hairpin_dg, compute_self_dimer_dg).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Constantes do espaco de acoes
# ---------------------------------------------------------------------------

# Tipos de modificacao que o agente pode aplicar a cada posicao
MODIFICATION_TYPES: list[str] = [
    "change_base",       # 0: trocar base (A/C/G/T)
    "toggle_lna",        # 1: ativar/desativar LNA na posicao
    "toggle_backbone",   # 2: alternar entre PS e PO
]

# Bases possiveis para substituicao
BASES: list[str] = ["A", "C", "G", "T"]

# Tipos de backbone
BACKBONE_PS: int = 1   # fosforotioato (mais resistente a nucleases)
BACKBONE_PO: int = 0   # fosfodiester (nativo, menos resistente)


@dataclass(frozen=True)
class RLConfig(AIModuleConfig):
    """Configuracao para otimizacao de ASO via REINFORCE.

    O ambiente modela o ASO como um vetor de posicoes, cada uma com:
    - Base (A/C/G/T)
    - Flag LNA (0/1)
    - Tipo de backbone (PS=1 / PO=0)

    O agente seleciona acoes (posicao, tipo_modificacao) para maximizar
    a recompensa multi-objetivo baseada em termodinamica.
    """

    # --- REINFORCE ---
    learning_rate: float = 1e-3      # taxa de aprendizado para policy gradient
    gamma: float = 0.99              # fator de desconto temporal
    baseline_decay: float = 0.95     # decaimento da baseline de media movel

    # --- Ambiente ---
    max_steps_per_episode: int = 50  # acoes por episodio antes de terminar
    convergence_threshold: float = 0.01  # diferenca de recompensa para convergir

    # --- Treinamento ---
    n_episodes: int = 500            # episodios de treinamento
    eval_interval: int = 50          # intervalo para reportar progresso
    top_k_variants: int = 5          # numero de melhores variantes a reportar

    # --- Pesos da recompensa ---
    # Cada metrica contribui proporcionalmente para a recompensa total.
    # dG e Tm sao os mais importantes para eficacia terapeutica;
    # GC content garante estabilidade biofisica; hairpin e self-dimer
    # penalizam estruturas secundarias indesejaveis.
    reward_weights: dict = field(default_factory=lambda: {
        "dg_binding": 0.30,       # delta_G de hibridizacao (mais negativo = melhor)
        "tm": 0.20,               # temperatura de melting (mais alto = melhor)
        "gc_content": 0.15,       # fracao GC (40-60% = otimo)
        "hairpin": 0.15,          # estabilidade de hairpin (menos negativo = melhor)
        "self_dimer": 0.10,       # auto-dimerizacao (menos negativo = melhor)
        "length_penalty": 0.10,   # penalidade se comprimento fora de 23-27 nt
    })

    # --- Limites biologicos para normalizacao de recompensa ---
    dg_best: float = -35.0        # dG ideal (kcal/mol) — ligacao muito forte
    dg_worst: float = -10.0       # dG limiar inferior de funcionalidade
    tm_best: float = 85.0         # Tm ideal (Celsius)
    tm_worst: float = 40.0        # Tm minimo funcional
    gc_optimal_low: float = 0.40  # limite inferior GC otimo
    gc_optimal_high: float = 0.60 # limite superior GC otimo
    length_min: int = 23          # comprimento minimo do ASO
    length_max: int = 27          # comprimento maximo do ASO

    # --- Restricoes de design do gapmer ---
    # MRL-ASO-001 usa design 5-15-5 (LNA-DNA-LNA)
    # O agente pode variar, mas deve respeitar limites biologicos
    min_lna_flank: int = 3        # minimo de LNA em cada flanco
    max_lna_flank: int = 7        # maximo de LNA em cada flanco
    min_dna_gap: int = 8          # minimo de DNA no gap central (RNase H)

    # --- Caminhos ---
    checkpoints_dir: Path = field(
        default_factory=lambda: AI_ROOT / "08_rl_ppo" / "checkpoints"
    )
