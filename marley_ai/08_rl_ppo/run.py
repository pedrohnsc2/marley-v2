"""Execucao do modulo 08_rl_ppo — Otimizacao de ASO por REINFORCE.

Treina um agente REINFORCE para descobrir variantes otimizadas do
ASO MRL-ASO-001. O agente explora modificacoes quimicas (base, LNA,
backbone) e aprende a maximizar recompensa termodinamica multi-objetivo.

Uso:
    python -m marley_ai.08_rl_ppo.run
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from aso_math.config import ASO_SEQUENCE
from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

from .agent import REINFORCEAgent
from .config import RLConfig
from .environment import ASODesignEnv


# ---------------------------------------------------------------------------
# Registro do modulo no orquestrador
# ---------------------------------------------------------------------------

@register("08_rl_ppo")
class RLPPOModule:
    """Modulo de otimizacao de ASO por aprendizado por reforco."""

    def configure(self, config: Any) -> None:
        self.config = config

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica dependencias: aso_math.thermo deve estar disponivel."""
        try:
            from aso_math.thermo import compute_dg, compute_tm
            return {"valid": True, "missing": []}
        except ImportError:
            return {"valid": False, "missing": ["aso_math.thermo"]}

    def run(self) -> dict[str, Any]:
        return main(self.config)

    def get_dependencies(self) -> list[str]:
        return []  # Depende apenas de aso_math (nao de outros modulos AI)


# ---------------------------------------------------------------------------
# Treinamento
# ---------------------------------------------------------------------------

def _train_agent(
    env: ASODesignEnv,
    agent: REINFORCEAgent,
    config: RLConfig,
) -> dict[str, Any]:
    """Loop de treinamento REINFORCE.

    Executa n_episodes episodios, cada um com max_steps_per_episode passos.
    Coleta as melhores variantes descobertas durante o treinamento.

    Args:
        env: Ambiente de design de ASO.
        agent: Agente REINFORCE.
        config: Configuracao do treinamento.

    Returns:
        Dicionario com historico de treinamento e melhores variantes.
    """
    # Rastreia as melhores variantes encontradas
    best_variants: list[dict[str, Any]] = []
    reward_history: list[float] = []
    loss_history: list[float] = []

    print(f"[08_rl_ppo] Iniciando treinamento: {config.n_episodes} episodios")
    print(f"[08_rl_ppo] Sequencia base: {ASO_SEQUENCE} (len={len(ASO_SEQUENCE)})")
    print(f"[08_rl_ppo] Espaco de acoes: {env.n_actions} acoes discretas")
    print(f"[08_rl_ppo] Observacao: {env.obs_dim} dimensoes")
    print()

    # Metricas do ASO original para comparacao
    baseline_metrics = env.get_metrics()
    print(f"[08_rl_ppo] Baseline MRL-ASO-001:")
    print(f"  dG = {baseline_metrics['dg_binding']:.2f} kcal/mol")
    print(f"  Tm = {baseline_metrics['tm']:.2f} C")
    print(f"  GC = {baseline_metrics['gc_content']:.2%}")
    print(f"  Hairpin dG = {baseline_metrics['hairpin_dg']:.2f} kcal/mol")
    print(f"  Self-dimer dG = {baseline_metrics['self_dimer_dg']:.2f} kcal/mol")
    print(f"  Reward = {baseline_metrics['reward']:.4f}")
    print()

    for episode in range(config.n_episodes):
        obs = env.reset()
        episode_reward = 0.0

        for step in range(config.max_steps_per_episode):
            # Seleciona acao estocasticamente
            action = agent.select_action(obs)

            # Executa acao no ambiente
            next_obs, reward, done, info = env.step(action)

            # Armazena transicao no buffer
            agent.buffer.store(obs, action, reward)
            episode_reward += reward

            obs = next_obs
            if done:
                break

        # Atualiza politica ao final do episodio
        loss = agent.update()

        # Registra metricas do episodio
        reward_history.append(episode_reward)
        loss_history.append(loss)

        # Avalia variante final do episodio
        final_metrics = env.get_metrics()
        _update_best_variants(best_variants, final_metrics, config.top_k_variants)

        # Reporta progresso
        if (episode + 1) % config.eval_interval == 0 or episode == 0:
            stats = agent.get_training_stats()
            print(
                f"[08_rl_ppo] Episodio {episode + 1:4d}/{config.n_episodes} | "
                f"Reward ep: {episode_reward:+.4f} | "
                f"Mean(50): {stats['mean_last_50']:+.4f} | "
                f"Loss: {loss:.4f} | "
                f"Best seq: {final_metrics['sequence'][:15]}..."
            )

    return {
        "reward_history": reward_history,
        "loss_history": loss_history,
        "best_variants": best_variants,
        "baseline_metrics": baseline_metrics,
        "training_stats": agent.get_training_stats(),
    }


def _update_best_variants(
    variants: list[dict[str, Any]],
    metrics: dict[str, Any],
    top_k: int,
) -> None:
    """Mantem lista ordenada das top_k melhores variantes.

    Evita duplicatas (mesma sequencia). Ordena por reward decrescente.
    """
    # Verifica se ja existe variante com mesma sequencia
    seq = metrics["sequence"]
    for v in variants:
        if v["sequence"] == seq:
            # Atualiza se reward melhor
            if metrics["reward"] > v["reward"]:
                v.update(metrics)
            return

    # Adiciona nova variante
    variants.append(dict(metrics))

    # Ordena e mantem top_k
    variants.sort(key=lambda x: x["reward"], reverse=True)
    while len(variants) > top_k:
        variants.pop()


# ---------------------------------------------------------------------------
# Avaliacao final
# ---------------------------------------------------------------------------

def _evaluate_greedy(
    env: ASODesignEnv,
    agent: REINFORCEAgent,
    config: RLConfig,
) -> dict[str, Any]:
    """Executa um episodio greedy (sem exploracao) para avaliacao.

    Usa a politica treinada de forma deterministica para gerar a
    melhor variante que o agente consegue produzir.
    """
    obs = env.reset()
    total_reward = 0.0

    for step in range(config.max_steps_per_episode):
        action = agent.select_action(obs, greedy=True)
        obs, reward, done, info = env.step(action)
        total_reward += reward
        if done:
            break

    return {
        "greedy_reward": total_reward,
        "greedy_metrics": env.get_metrics(),
    }


# ---------------------------------------------------------------------------
# Ponto de entrada
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Otimiza ASO MRL-ASO-001 via REINFORCE.

    1. Cria ambiente com ASO referencia como ponto de partida
    2. Treina agente REINFORCE por N episodios
    3. Avalia politica treinada em modo greedy
    4. Reporta top 5 variantes vs baseline
    5. Salva resultados em JSON

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com resultados da otimizacao.
    """
    # Configuracao
    rl_config = RLConfig(
        module_slug="08_rl_ppo",
        module_name="ASO RL Optimization (REINFORCE)",
    )

    envelope = create_envelope("08_rl_ppo", version="0.2.0")

    with Timer() as timer:
        # Cria ambiente e agente
        env = ASODesignEnv(config=rl_config)
        agent = REINFORCEAgent(
            obs_dim=env.obs_dim,
            n_actions=env.n_actions,
            config=rl_config,
        )

        # Treinamento
        results = _train_agent(env, agent, rl_config)

        # Avaliacao greedy
        greedy_results = _evaluate_greedy(env, agent, rl_config)
        print()
        print("[08_rl_ppo] === Avaliacao Greedy (politica treinada) ===")
        gm = greedy_results["greedy_metrics"]
        print(f"  Sequencia: {gm['sequence']}")
        print(f"  dG = {gm['dg_binding']:.2f} kcal/mol")
        print(f"  Tm = {gm['tm']:.2f} C")
        print(f"  GC = {gm['gc_content']:.2%}")
        print(f"  Reward = {gm['reward']:.4f}")
        print()

        # Comparacao: top 5 variantes vs baseline
        baseline = results["baseline_metrics"]
        print("[08_rl_ppo] === Top 5 Variantes vs MRL-ASO-001 ===")
        print(f"  {'Rank':<5} {'Sequencia':<28} {'dG':>8} {'Tm':>8} {'GC':>7} {'Reward':>8} {'vs Base':>8}")
        print(f"  {'----':<5} {'--------':<28} {'------':>8} {'----':>8} {'----':>7} {'------':>8} {'-------':>8}")
        print(
            f"  {'BASE':<5} {baseline['sequence']:<28} "
            f"{baseline['dg_binding']:>8.2f} {baseline['tm']:>8.2f} "
            f"{baseline['gc_content']:>7.2%} {baseline['reward']:>8.4f} {'---':>8}"
        )
        for i, variant in enumerate(results["best_variants"]):
            delta = variant["reward"] - baseline["reward"]
            sign = "+" if delta >= 0 else ""
            print(
                f"  {i + 1:<5} {variant['sequence']:<28} "
                f"{variant['dg_binding']:>8.2f} {variant['tm']:>8.2f} "
                f"{variant['gc_content']:>7.2%} {variant['reward']:>8.4f} "
                f"{sign}{delta:>7.4f}"
            )
        print()

        # Monta envelope de resultados
        envelope["status"] = "completed"
        envelope["device"] = "cpu"
        envelope["dependencies"] = ["aso_math.thermo"]

        envelope["summary"]["conclusion"] = (
            f"Agente REINFORCE treinado por {rl_config.n_episodes} episodios. "
            f"Top variante: reward {results['best_variants'][0]['reward']:.4f} "
            f"(baseline: {baseline['reward']:.4f}). "
            f"Sequencia: {results['best_variants'][0]['sequence']}."
        )

        envelope["summary"]["key_metrics"] = {
            "n_episodes": rl_config.n_episodes,
            "baseline_reward": baseline["reward"],
            "best_reward": results["best_variants"][0]["reward"],
            "best_dg": results["best_variants"][0]["dg_binding"],
            "best_tm": results["best_variants"][0]["tm"],
            "improvement_pct": round(
                (results["best_variants"][0]["reward"] - baseline["reward"])
                / max(abs(baseline["reward"]), 1e-8) * 100, 2
            ),
        }

        envelope["metrics"] = results["training_stats"]

        envelope["data"] = {
            "baseline": baseline,
            "best_variants": results["best_variants"],
            "greedy_evaluation": greedy_results,
            "reward_history_summary": {
                "first_10_mean": float(np.mean(results["reward_history"][:10])),
                "last_10_mean": float(np.mean(results["reward_history"][-10:])),
                "overall_mean": float(np.mean(results["reward_history"])),
                "overall_std": float(np.std(results["reward_history"])),
                "max": float(np.max(results["reward_history"])),
                "min": float(np.min(results["reward_history"])),
            },
            "config": {
                "n_episodes": rl_config.n_episodes,
                "max_steps_per_episode": rl_config.max_steps_per_episode,
                "learning_rate": rl_config.learning_rate,
                "gamma": rl_config.gamma,
                "reward_weights": rl_config.reward_weights,
            },
        }

    envelope["runtime_seconds"] = timer.elapsed

    # Salva resultados
    output_dir = Path(__file__).resolve().parent / "results"
    output_path = write_result(envelope, output_dir=output_dir)
    print(f"[08_rl_ppo] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
