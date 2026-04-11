"""Execucao do modulo 07_contrastive — Aprendizado contrastivo epitopo-MHC.

Treina um modelo contrastivo (InfoNCE) para alinhar embeddings de
epitopos e alelos DLA no mesmo espaco vetorial. Avalia a correlacao
entre os scores aprendidos e os valores de IC50 experimentais.

Uso:
    python -m marley_ai.07_contrastive.run
"""

from __future__ import annotations

from typing import Any

import numpy as np

from vaccine_platforms.shared.epitopes import EPITOPES

from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

from .config import ContrastiveConfig
from .encoder import (
    DLA_ALLELE_PROPERTIES,
    build_kmer_vocab,
)
from .trainer import score_epitope_allele, train_contrastive


@register("07_contrastive")
class ContrastiveModule:
    """Modulo de aprendizado contrastivo epitopo-MHC canino.

    Aprende um espaco de embedding compartilhado onde epitopos e
    seus alelos DLA ligantes ficam proximos. Permite predizer
    afinidade de ligacao para pares nao vistos.
    """

    def __init__(self) -> None:
        self._config: ContrastiveConfig | None = None

    def configure(self, config: Any) -> None:
        """Recebe e armazena configuracao do modulo."""
        if isinstance(config, ContrastiveConfig):
            self._config = config
        else:
            self._config = ContrastiveConfig(
                module_slug="07_contrastive",
                module_name="Contrastive Epitope-MHC Learning",
            )

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica se os epitopos e numpy estao disponiveis."""
        missing = []
        if len(EPITOPES) == 0:
            missing.append("vaccine_platforms.shared.epitopes.EPITOPES")
        try:
            import numpy  # noqa: F401
        except ImportError:
            missing.append("numpy")
        return {"valid": len(missing) == 0, "missing": missing}

    def run(self) -> dict[str, Any]:
        """Executa o pipeline contrastivo completo."""
        if self._config is None:
            self.configure(None)
        assert self._config is not None
        return main(self._config)

    def get_dependencies(self) -> list[str]:
        """Sem dependencias de outros modulos AI."""
        return []


def _rankdata(x: np.ndarray) -> np.ndarray:
    """Calcula ranks com tratamento de empates (metodo average).

    Implementacao em numpy puro, equivalente a scipy.stats.rankdata.

    Args:
        x: array de valores para ranquear

    Returns:
        Array de ranks (1-based, empates recebem rank medio)
    """
    n = len(x)
    ranks = np.empty(n, dtype=np.float64)
    sorted_indices = np.argsort(x)
    i = 0
    while i < n:
        j = i
        # Encontrar grupo de empates
        while j < n - 1 and x[sorted_indices[j + 1]] == x[sorted_indices[j]]:
            j += 1
        # Rank medio para o grupo
        avg_rank = (i + j) / 2.0 + 1.0  # 1-based
        for k in range(i, j + 1):
            ranks[sorted_indices[k]] = avg_rank
        i = j + 1
    return ranks


def _spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """Correlacao de Spearman (rho) entre dois arrays.

    Spearman rho e a correlacao de Pearson aplicada aos ranks.
    Mede associacao monotonica (nao necessariamente linear).

    Args:
        x: primeiro array de valores
        y: segundo array de valores

    Returns:
        Coeficiente rho no intervalo [-1, 1]
    """
    rank_x = _rankdata(x)
    rank_y = _rankdata(y)
    n = len(x)

    # Pearson sobre os ranks
    mean_rx = np.mean(rank_x)
    mean_ry = np.mean(rank_y)
    dx = rank_x - mean_rx
    dy = rank_y - mean_ry

    numerator = np.sum(dx * dy)
    denominator = np.sqrt(np.sum(dx ** 2) * np.sum(dy ** 2))

    if denominator < 1e-10:
        return 0.0
    return float(numerator / denominator)


def _kendall_tau(x: np.ndarray, y: np.ndarray) -> float:
    """Correlacao de Kendall (tau-b) entre dois arrays.

    Conta pares concordantes e discordantes. Mais robusto que
    Spearman para amostras pequenas (como nossos 11 epitopos).

    Args:
        x: primeiro array de valores
        y: segundo array de valores

    Returns:
        Coeficiente tau no intervalo [-1, 1]
    """
    n = len(x)
    concordant = 0
    discordant = 0
    tied_x = 0
    tied_y = 0

    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]

            if dx == 0 and dy == 0:
                tied_x += 1
                tied_y += 1
            elif dx == 0:
                tied_x += 1
            elif dy == 0:
                tied_y += 1
            elif (dx > 0 and dy > 0) or (dx < 0 and dy < 0):
                concordant += 1
            else:
                discordant += 1

    n_pairs = n * (n - 1) / 2
    denom = np.sqrt((n_pairs - tied_x) * (n_pairs - tied_y))

    if denom < 1e-10:
        return 0.0
    return float((concordant - discordant) / denom)


def _compute_ranking_correlation(scores_by_epitope: dict[str, float]) -> dict[str, float]:
    """Calcula correlacao entre ranking por score e ranking por IC50.

    Se o modelo contrastivo aprendeu uma boa funcao de afinidade,
    epitopos com IC50 baixo (alta afinidade) devem ter scores altos.
    Portanto, esperamos correlacao negativa entre IC50 e score
    (ou correlacao positiva entre rank_IC50 e rank_score).

    Implementacao em numpy puro (sem scipy).

    Args:
        scores_by_epitope: mapeamento peptide -> score contrastivo

    Returns:
        Dicionario com spearman_rho e kendall_tau
    """
    # Alinhar epitopos que estao nos scores e nos EPITOPES
    peptides = []
    ic50_values = []
    score_values = []

    for ep in EPITOPES:
        if ep.peptide in scores_by_epitope:
            peptides.append(ep.peptide)
            ic50_values.append(ep.ic50)
            score_values.append(scores_by_epitope[ep.peptide])

    if len(peptides) < 3:
        return {"spearman_rho": 0.0, "kendall_tau": 0.0}

    ic50_arr = np.array(ic50_values)
    score_arr = np.array(score_values)

    # IC50 baixo = boa afinidade, score alto = boa afinidade
    # Esperamos correlacao negativa entre IC50 e score
    rho = _spearman_rho(ic50_arr, score_arr)
    tau = _kendall_tau(ic50_arr, score_arr)

    return {
        "spearman_rho": round(rho, 4),
        "kendall_tau": round(tau, 4),
    }


def main(config: ContrastiveConfig | None = None) -> dict[str, Any]:
    """Treina modelo contrastivo e avalia correlacao com IC50.

    Pipeline completo:
        1. Configurar parametros
        2. Gerar pares de treinamento (positivos dos 11 epitopos + negativos)
        3. Treinar MLP contrastivo com InfoNCE loss
        4. Scorer todos os 11 epitopos contra todos os 3 alelos DLA
        5. Avaliar: o ranking aprendido correlaciona com IC50?
        6. Salvar resultados

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com metricas de treinamento e scores.
    """
    if config is None:
        config = ContrastiveConfig(
            module_slug="07_contrastive",
            module_name="Contrastive Epitope-MHC Learning",
        )

    envelope = create_envelope("07_contrastive")
    envelope["device"] = "cpu"  # numpy puro, sem GPU

    print("[07_contrastive] Iniciando aprendizado contrastivo epitopo-MHC")
    print(f"  Embedding dim: {config.embedding_dim}, Hidden dim: {config.hidden_dim}")
    print(f"  Epocas: {config.n_epochs}, LR: {config.learning_rate}, Temperatura: {config.temperature}")
    print(f"  Negativos por positivo: {config.neg_ratio}, Hard neg fraction: {config.hard_neg_fraction}")
    print(f"  {len(EPITOPES)} epitopos canonicos, {len(DLA_ALLELE_PROPERTIES)} alelos DLA")
    print()

    with Timer() as timer:
        # --- Fase 1: Treinamento ---
        print("  [fase 1] Treinamento contrastivo...")
        result = train_contrastive(config)
        print()

        # --- Fase 2: Scoring de todos os pares epitopo-alelo ---
        print("  [fase 2] Scoring de todos os pares epitopo-DLA...")
        kmer_vocab = build_kmer_vocab(config.kmer_k)
        allele_names = sorted(DLA_ALLELE_PROPERTIES.keys())

        # Tabela de scores: epitopo x alelo
        score_table: list[dict[str, Any]] = []
        # Score do par original (epitopo com seu alelo ligante)
        binding_scores: dict[str, float] = {}

        for epitope in EPITOPES:
            row: dict[str, Any] = {
                "peptide": epitope.peptide,
                "gene_name": epitope.gene_name,
                "known_allele": epitope.allele,
                "ic50": epitope.ic50,
                "rank": epitope.rank,
                "scores": {},
            }

            best_allele = ""
            best_score = -2.0

            for allele in allele_names:
                score = score_epitope_allele(
                    epitope.peptide,
                    allele,
                    result.epitope_weights,
                    result.allele_weights,
                    config.kmer_k,
                    kmer_vocab,
                )
                row["scores"][allele] = round(score, 4)

                if score > best_score:
                    best_score = score
                    best_allele = allele

            row["predicted_best_allele"] = best_allele
            row["predicted_best_score"] = round(best_score, 4)
            row["correct_prediction"] = best_allele == epitope.allele

            # Score do par verdadeiro (para correlacao com IC50)
            binding_scores[epitope.peptide] = row["scores"][epitope.allele]

            score_table.append(row)

        # --- Fase 3: Avaliacao ---
        print("  [fase 3] Avaliacao de correlacao score vs IC50...")
        correlation = _compute_ranking_correlation(binding_scores)

        # Contar predicoes corretas de alelo
        n_correct = sum(1 for row in score_table if row["correct_prediction"])
        accuracy = n_correct / len(score_table) if score_table else 0.0

        # Imprimir tabela de resultados
        print()
        print("  " + "=" * 80)
        print(f"  {'Epitopo':<12s} {'IC50':>7s} {'Alelo Real':<15s} {'Score':>7s} {'Predito':<15s} {'OK':>3s}")
        print("  " + "-" * 80)
        for row in score_table:
            peptide = row["peptide"]
            ic50 = row["ic50"]
            known = row["known_allele"]
            score = row["scores"][known]
            predicted = row["predicted_best_allele"]
            correct = "SIM" if row["correct_prediction"] else "NAO"
            print(f"  {peptide:<12s} {ic50:>7.2f} {known:<15s} {score:>7.4f} {predicted:<15s} {correct:>3s}")
        print("  " + "=" * 80)
        print()

        print(f"  Acuracia de predicao de alelo: {n_correct}/{len(score_table)} ({accuracy:.1%})")
        print(f"  Spearman rho (IC50 vs score): {correlation['spearman_rho']:.4f}")
        print(f"  Kendall tau (IC50 vs score):  {correlation['kendall_tau']:.4f}")
        print()

        # Ranking dos epitopos por score do par verdadeiro (maior = melhor)
        ranking_by_score = sorted(
            [(row["peptide"], row["scores"][row["known_allele"]], row["ic50"])
             for row in score_table],
            key=lambda x: x[1],
            reverse=True,
        )
        print("  Ranking por score contrastivo (maior = melhor afinidade predita):")
        for i, (pep, scr, ic50) in enumerate(ranking_by_score, 1):
            print(f"    {i:>2d}. {pep:<12s}  score={scr:>7.4f}  IC50={ic50:>7.2f} nM")
        print()

        # --- Montar envelope ---
        envelope["status"] = "complete"
        envelope["summary"]["conclusion"] = (
            f"Modelo contrastivo treinado com {config.n_epochs} epocas. "
            f"Loss final: {result.final_loss:.4f}. "
            f"Acuracia de predicao de alelo: {accuracy:.1%}. "
            f"Correlacao Spearman IC50 vs score: {correlation['spearman_rho']:.4f}."
        )
        envelope["summary"]["key_metrics"] = {
            "embedding_dim": config.embedding_dim,
            "hidden_dim": config.hidden_dim,
            "n_epochs": config.n_epochs,
            "n_train_pairs": result.n_train_pairs,
            "final_loss": round(result.final_loss, 4),
            "allele_prediction_accuracy": round(accuracy, 4),
            "spearman_rho": correlation["spearman_rho"],
            "kendall_tau": correlation["kendall_tau"],
        }

        envelope["metrics"] = {
            "loss_history": [round(l, 4) for l in result.loss_history],
            "loss_first": round(result.loss_history[0], 4),
            "loss_final": round(result.final_loss, 4),
            "loss_reduction": round(result.loss_history[0] - result.final_loss, 4),
        }

        envelope["data"] = {
            "score_table": score_table,
            "ranking_by_score": [
                {"rank": i + 1, "peptide": pep, "score": round(scr, 4), "ic50": ic50}
                for i, (pep, scr, ic50) in enumerate(ranking_by_score)
            ],
            "correlation": correlation,
            "config": {
                "embedding_dim": config.embedding_dim,
                "hidden_dim": config.hidden_dim,
                "kmer_k": config.kmer_k,
                "temperature": config.temperature,
                "learning_rate": config.learning_rate,
                "momentum": config.momentum,
                "n_epochs": config.n_epochs,
                "neg_ratio": config.neg_ratio,
                "hard_neg_fraction": config.hard_neg_fraction,
                "mutation_positions": config.mutation_positions,
                "weight_decay": config.weight_decay,
                "seed": config.seed,
            },
        }

        envelope["dependencies"] = ["vaccine_platforms.shared.epitopes"]

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"[07_contrastive] Resultado salvo em {output_path}")
    print(f"[07_contrastive] Tempo total: {timer.elapsed:.2f}s")

    return envelope


if __name__ == "__main__":
    main()
