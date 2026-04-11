"""Execucao do modulo 06_evodiff — Geracao por difusao discreta.

Treina um modelo de difusao discreta do zero em PyTorch para gerar
variantes otimizadas de ASOs (nucleotideos) e epitopos (aminoacidos).

Modo 1 — ASO variant generation:
    - Treina em MRL-ASO-001 + biblioteca de mutantes (75 sequencias)
    - Gera 50 novas variantes de ASO
    - Avalia com funcoes termodinamicas (dG, Tm, GC%)
    - Reporta top 10 vs MRL-ASO-001

Modo 2 — Epitope variant generation:
    - Treina nos 11 epitopos canonicos + negativos embaralhados
    - Gera 50 variantes de epitopos
    - Avalia propriedades fisicoquimicas (MW, hidrofobicidade, carga)

Uso:
    python -m marley_ai.06_evodiff.run
"""

from __future__ import annotations

import importlib as _il
import random
from pathlib import Path
from typing import Any

import numpy as np
import torch

from aso_math.config import ASO_SEQUENCE, SL_SEQUENCE
from aso_math.thermo import compute_dg, compute_tm, gc_content
from marley_ai.config import AIModuleConfig, detect_device
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register
from vaccine_platforms.shared.epitopes import EPITOPES, get_epitope_sequences

# Importacao via importlib — Python nao permite "from marley_ai.06_evodiff..."
# porque "06" seria interpretado como literal numerico invalido.
_cfg = _il.import_module("marley_ai.06_evodiff.config")
_diff = _il.import_module("marley_ai.06_evodiff.diffusion")
_mdl = _il.import_module("marley_ai.06_evodiff.model")
_evl = _il.import_module("marley_ai.06_evodiff.evaluator")

EvoDiffConfig = _cfg.EvoDiffConfig
NUCLEOTIDE_ALPHABET = _cfg.NUCLEOTIDE_ALPHABET
AMINO_ACID_ALPHABET = _cfg.AMINO_ACID_ALPHABET

linear_beta_schedule = _diff.linear_beta_schedule
compute_alpha_bar = _diff.compute_alpha_bar
reverse_diffusion_sample = _diff.reverse_diffusion_sample
tokens_to_sequences = _diff.tokens_to_sequences
sequences_to_tokens = _diff.sequences_to_tokens

SequenceDenoiser = _mdl.SequenceDenoiser
compute_loss = _mdl.compute_loss

evaluate_batch = _evl.evaluate_batch
filter_and_rank = _evl.filter_and_rank


# ---------------------------------------------------------------------------
# Registro do modulo no orquestrador
# ---------------------------------------------------------------------------

@register("06_evodiff")
class EvoDiffModule:
    """Modulo de geracao de sequencias via difusao discreta."""

    def __init__(self) -> None:
        self._config: EvoDiffConfig | None = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao do orquestrador."""
        if isinstance(config, EvoDiffConfig):
            self._config = config
        else:
            self._config = EvoDiffConfig(
                module_slug="06_evodiff",
                module_name="Discrete Diffusion Generation",
            )

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica dependencias: torch + aso_math devem estar disponiveis."""
        missing = []
        try:
            import torch
        except ImportError:
            missing.append("torch")
        try:
            from aso_math.thermo import compute_dg
        except ImportError:
            missing.append("aso_math.thermo")
        return {"valid": len(missing) == 0, "missing": missing}

    def run(self) -> dict[str, Any]:
        """Executa pipeline completo de difusao."""
        return main(self._config)

    def get_dependencies(self) -> list[str]:
        """Depende de aso_math (termodinamica) e epitopos."""
        return []


# ---------------------------------------------------------------------------
# Construcao do dataset de treinamento
# ---------------------------------------------------------------------------

def _build_aso_training_set(seed: int = 42) -> list[str]:
    """Constroi biblioteca de treinamento para ASO a partir de MRL-ASO-001.

    Gera a sequencia wildtype + todas as mutacoes pontuais possiveis:
    para cada posicao (25 nt), substitui por cada uma das 3 bases
    alternativas, gerando 25 * 3 = 75 mutantes + 1 wildtype = 76 seqs.

    Args:
        seed: Semente para reproducibilidade.

    Returns:
        Lista de 76 sequencias de ASO (wildtype + 75 mutantes).
    """
    rng = random.Random(seed)
    bases = list(NUCLEOTIDE_ALPHABET)
    wt = ASO_SEQUENCE.upper()
    library = [wt]

    # Gera todos os mutantes pontuais (single-point mutations)
    for pos in range(len(wt)):
        for base in bases:
            if base != wt[pos]:
                mutant = wt[:pos] + base + wt[pos + 1:]
                library.append(mutant)

    # Embaralha para diversificar batches durante treinamento
    rng.shuffle(library)

    return library


def _build_epitope_training_set(seed: int = 42) -> list[str]:
    """Constroi dataset de treinamento para epitopos.

    Inclui os 11 epitopos canonicos + controles negativos (shuffled).
    Os negativos sao versoes embaralhadas dos epitopos originais,
    mantendo composicao aminoacidica mas destruindo a sequencia.

    Gera 3 shuffled por epitopo = 11 + 33 = 44 sequencias.

    Args:
        seed: Semente para reproducibilidade.

    Returns:
        Lista de sequencias de treinamento.
    """
    rng = random.Random(seed)
    epitope_seqs = get_epitope_sequences()
    training = list(epitope_seqs)

    # Gera controles embaralhados (3 por epitopo)
    for seq in epitope_seqs:
        for _ in range(3):
            shuffled = list(seq)
            rng.shuffle(shuffled)
            training.append("".join(shuffled))

    return training


# ---------------------------------------------------------------------------
# Treinamento
# ---------------------------------------------------------------------------

def _train_model(
    model: SequenceDenoiser,
    training_tokens: torch.Tensor,
    alpha_bar: torch.Tensor,
    config: EvoDiffConfig,
    device: str,
    vocab_size: int | None = None,
) -> list[float]:
    """Treina o modelo de denoising por N epocas.

    A cada epoca, amostra batches do dataset de treinamento,
    aplica corrupcao em timesteps aleatorios, e otimiza a cross-entropy
    entre a predicao do modelo e os tokens originais.

    Args:
        model: Rede de denoising.
        training_tokens: Tokens de treinamento, shape (N, L).
        alpha_bar: Schedule de difusao, shape (T,).
        config: Configuracao do modulo.
        device: Dispositivo de computacao.
        vocab_size: Tamanho do vocabulario. Se None, usa config.get_vocab_size().

    Returns:
        Lista de perdas por epoca.
    """
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)
    if vocab_size is None:
        vocab_size = config.get_vocab_size()
    n_seqs = training_tokens.shape[0]
    loss_history: list[float] = []

    model.train()

    for epoch in range(config.n_epochs):
        # Embaralha o dataset
        perm = torch.randperm(n_seqs, device=device)
        epoch_loss = 0.0
        n_batches = 0

        for start in range(0, n_seqs, config.batch_size):
            end = min(start + config.batch_size, n_seqs)
            batch_idx = perm[start:end]
            x_0 = training_tokens[batch_idx]

            loss = compute_loss(model, x_0, alpha_bar, vocab_size)

            optimizer.zero_grad()
            loss.backward()
            # Gradient clipping para estabilidade no MPS
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        avg_loss = epoch_loss / max(n_batches, 1)
        loss_history.append(avg_loss)

        # Reporta progresso a cada 20 epocas
        if (epoch + 1) % 20 == 0 or epoch == 0:
            print(f"  Epoca {epoch + 1:4d}/{config.n_epochs} | Loss: {avg_loss:.4f}")

    return loss_history


# ---------------------------------------------------------------------------
# Pipeline ASO
# ---------------------------------------------------------------------------

def _run_aso_pipeline(
    config: EvoDiffConfig,
    device: str,
) -> dict[str, Any]:
    """Pipeline completo de geracao de variantes de ASO.

    1. Constroi biblioteca de treinamento (MRL-ASO-001 + mutantes)
    2. Treina modelo de difusao em nucleotideos
    3. Gera 50 sequencias novas
    4. Avalia com funcoes termodinamicas
    5. Reporta top 10 vs baseline

    Args:
        config: Configuracao do modulo.
        device: Dispositivo de computacao.

    Returns:
        Dicionario com resultados do pipeline ASO.
    """
    print("\n[06_evodiff] === Pipeline ASO ===")
    alphabet = NUCLEOTIDE_ALPHABET
    vocab_size = config.get_vocab_size()
    seq_len = len(ASO_SEQUENCE)

    # 1. Dataset de treinamento
    print("[06_evodiff] Construindo biblioteca de treinamento...")
    training_seqs = _build_aso_training_set(seed=config.seed)
    print(f"  Sequencias de treinamento: {len(training_seqs)}")
    print(f"  Wildtype: {ASO_SEQUENCE} (len={seq_len})")

    training_tokens = sequences_to_tokens(training_seqs, alphabet, device=device)

    # 2. Modelo e schedule
    print(f"\n[06_evodiff] Inicializando modelo (d={config.d_model}, "
          f"h={config.n_heads}, L={config.n_layers})...")

    betas = linear_beta_schedule(
        config.n_diffusion_steps, config.beta_start, config.beta_end,
    ).to(device)
    alpha_bar = compute_alpha_bar(betas)

    model = SequenceDenoiser(
        vocab_size=vocab_size,
        d_model=config.d_model,
        n_heads=config.n_heads,
        n_layers=config.n_layers,
        dropout=config.dropout,
        max_seq_len=seq_len + 8,  # Margem para sequencias ligeiramente maiores
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Parametros: {n_params:,}")

    # 3. Treinamento
    print(f"\n[06_evodiff] Treinando por {config.n_epochs} epocas...")
    loss_history = _train_model(
        model, training_tokens, alpha_bar, config, device,
    )

    # 4. Geracao
    print(f"\n[06_evodiff] Gerando {config.n_samples} variantes de ASO...")
    generated_tokens = reverse_diffusion_sample(
        model=model,
        seq_len=seq_len,
        vocab_size=vocab_size,
        betas=betas,
        alpha_bar=alpha_bar,
        n_samples=config.n_samples,
        device=device,
    )
    generated_seqs = tokens_to_sequences(generated_tokens, alphabet)

    # Remove duplicatas mantendo ordem
    seen = set()
    unique_seqs = []
    for s in generated_seqs:
        if s not in seen:
            seen.add(s)
            unique_seqs.append(s)

    print(f"  Sequencias unicas: {len(unique_seqs)}/{config.n_samples}")

    # 5. Avaliacao
    print("\n[06_evodiff] Avaliando variantes...")
    evaluations = evaluate_batch(
        unique_seqs, training_seqs, mode="aso",
        gc_min=config.gc_min, gc_max=config.gc_max,
        homopolymer_max=config.homopolymer_max,
    )
    top_k, all_valid = filter_and_rank(evaluations, top_k=config.top_k)

    # Baseline MRL-ASO-001
    baseline = {
        "sequence": ASO_SEQUENCE,
        "dg": compute_dg(ASO_SEQUENCE),
        "tm": compute_tm(ASO_SEQUENCE),
        "gc": round(gc_content(ASO_SEQUENCE), 4),
    }

    # Relatorio
    print(f"\n[06_evodiff] Candidatos validos: {len(all_valid)}/{len(unique_seqs)}")
    print(f"\n[06_evodiff] === Top {min(len(top_k), config.top_k)} Variantes de ASO vs MRL-ASO-001 ===")
    print(f"  {'Rank':<5} {'Sequencia':<28} {'dG':>8} {'Tm':>8} {'GC':>7} {'Edit':>5} {'Score':>7}")
    print(f"  {'----':<5} {'--------':<28} {'------':>8} {'----':>8} {'----':>7} {'----':>5} {'-----':>7}")
    print(
        f"  {'BASE':<5} {baseline['sequence']:<28} "
        f"{baseline['dg']:>8.2f} {baseline['tm']:>8.2f} "
        f"{baseline['gc']:>7.2%} {'---':>5} {'---':>7}"
    )
    for i, v in enumerate(top_k):
        print(
            f"  {i + 1:<5} {v['sequence']:<28} "
            f"{v['dg']:>8.2f} {v['tm']:>8.2f} "
            f"{v['gc']:>7.2%} {v['edit_distance']:>5} {v['score']:>7.4f}"
        )

    return {
        "mode": "aso",
        "training_size": len(training_seqs),
        "generated": len(generated_seqs),
        "unique": len(unique_seqs),
        "valid": len(all_valid),
        "top_k": top_k,
        "all_valid": all_valid,
        "baseline": baseline,
        "loss_history": loss_history,
        "n_params": n_params,
    }


# ---------------------------------------------------------------------------
# Pipeline Epitopo
# ---------------------------------------------------------------------------

def _run_epitope_pipeline(
    config: EvoDiffConfig,
    device: str,
) -> dict[str, Any]:
    """Pipeline completo de geracao de variantes de epitopos.

    1. Constroi dataset (11 epitopos + negativos embaralhados)
    2. Treina modelo de difusao em aminoacidos
    3. Gera 50 peptideos candidatos
    4. Avalia propriedades fisicoquimicas
    5. Reporta top 10

    Args:
        config: Configuracao do modulo (sobrescreve sequence_type).
        device: Dispositivo de computacao.

    Returns:
        Dicionario com resultados do pipeline epitopo.
    """
    print("\n[06_evodiff] === Pipeline Epitopo ===")
    alphabet = AMINO_ACID_ALPHABET
    # Vocabulario para aminoacidos: 20 AAs + 1 MASK
    vocab_size = len(alphabet) + 1
    # Epitopos sao 9-mers (padrao MHC-I)
    seq_len = 9

    # 1. Dataset de treinamento
    print("[06_evodiff] Construindo dataset de treinamento...")
    training_seqs = _build_epitope_training_set(seed=config.seed)
    print(f"  Sequencias de treinamento: {len(training_seqs)}")
    print(f"  Epitopos originais: {len(get_epitope_sequences())}")

    training_tokens = sequences_to_tokens(training_seqs, alphabet, device=device)

    # 2. Modelo e schedule
    print(f"\n[06_evodiff] Inicializando modelo epitopo (d={config.d_model}, "
          f"h={config.n_heads}, L={config.n_layers})...")

    betas = linear_beta_schedule(
        config.n_diffusion_steps, config.beta_start, config.beta_end,
    ).to(device)
    alpha_bar = compute_alpha_bar(betas)

    model = SequenceDenoiser(
        vocab_size=vocab_size,
        d_model=config.d_model,
        n_heads=config.n_heads,
        n_layers=config.n_layers,
        dropout=config.dropout,
        max_seq_len=seq_len + 8,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Parametros: {n_params:,}")

    # 3. Treinamento (mais epocas para dataset menor)
    n_epochs_epitope = config.n_epochs
    print(f"\n[06_evodiff] Treinando por {n_epochs_epitope} epocas...")

    # Passa vocab_size explicitamente (20 AAs + 1 MASK, nao o default nucleotideo)
    loss_history = _train_model(
        model, training_tokens, alpha_bar, config, device,
        vocab_size=vocab_size,
    )

    # 4. Geracao
    print(f"\n[06_evodiff] Gerando {config.n_samples} variantes de epitopos...")
    generated_tokens = reverse_diffusion_sample(
        model=model,
        seq_len=seq_len,
        vocab_size=vocab_size,
        betas=betas,
        alpha_bar=alpha_bar,
        n_samples=config.n_samples,
        device=device,
    )
    generated_seqs = tokens_to_sequences(generated_tokens, alphabet)

    # Remove duplicatas
    seen = set()
    unique_seqs = []
    for s in generated_seqs:
        if s not in seen:
            seen.add(s)
            unique_seqs.append(s)

    print(f"  Sequencias unicas: {len(unique_seqs)}/{config.n_samples}")

    # 5. Avaliacao
    print("\n[06_evodiff] Avaliando variantes de epitopos...")
    epitope_originals = get_epitope_sequences()
    evaluations = evaluate_batch(
        unique_seqs, epitope_originals, mode="epitope",
        min_hydro=config.min_hydrophobicity,
        max_hydro=config.max_hydrophobicity,
    )
    top_k, all_valid = filter_and_rank(evaluations, top_k=config.top_k)

    # Relatorio
    print(f"\n[06_evodiff] Candidatos validos: {len(all_valid)}/{len(unique_seqs)}")
    print(f"\n[06_evodiff] === Top {min(len(top_k), config.top_k)} Variantes de Epitopos ===")
    print(f"  {'Rank':<5} {'Sequencia':<12} {'MW':>8} {'Hydro':>7} {'Charge':>7} {'Edit':>5} {'Score':>7}")
    print(f"  {'----':<5} {'--------':<12} {'----':>8} {'-----':>7} {'------':>7} {'----':>5} {'-----':>7}")
    for i, v in enumerate(top_k):
        print(
            f"  {i + 1:<5} {v['sequence']:<12} "
            f"{v['mw']:>8.1f} {v['hydrophobicity']:>7.2f} "
            f"{v['charge']:>7.1f} {v['edit_distance']:>5} {v['score']:>7.4f}"
        )

    return {
        "mode": "epitope",
        "training_size": len(training_seqs),
        "generated": len(generated_seqs),
        "unique": len(unique_seqs),
        "valid": len(all_valid),
        "top_k": top_k,
        "all_valid": all_valid,
        "epitope_originals": epitope_originals,
        "loss_history": loss_history,
        "n_params": n_params,
    }


# ---------------------------------------------------------------------------
# Ponto de entrada
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Pipeline completo de geracao por difusao discreta.

    Executa ambos os pipelines (ASO e epitopos), combina resultados
    em um envelope unificado, e salva como JSON.

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com candidatos gerados.
    """
    if config is None:
        config = EvoDiffConfig(
            module_slug="06_evodiff",
            module_name="Discrete Diffusion Generation",
        )

    # Detecta dispositivo disponivel
    device = detect_device()
    print(f"[06_evodiff] Dispositivo: {device}")

    # Semente para reproducibilidade
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    random.seed(config.seed)

    envelope = create_envelope("06_evodiff", version="0.2.0")
    envelope["device"] = device

    with Timer() as timer:

        # --- Pipeline 1: ASO variants ---
        aso_results = _run_aso_pipeline(config, device)

        # --- Pipeline 2: Epitope variants ---
        epitope_results = _run_epitope_pipeline(config, device)

        # --- Monta envelope de resultados ---
        envelope["status"] = "completed"
        envelope["dependencies"] = [
            "aso_math.thermo",
            "aso_math.config",
            "vaccine_platforms.shared.epitopes",
        ]

        # Metricas de treinamento
        envelope["metrics"] = {
            "aso_model": {
                "n_params": aso_results["n_params"],
                "final_loss": aso_results["loss_history"][-1] if aso_results["loss_history"] else 0.0,
                "training_size": aso_results["training_size"],
                "n_generated": aso_results["generated"],
                "n_unique": aso_results["unique"],
                "n_valid": aso_results["valid"],
            },
            "epitope_model": {
                "n_params": epitope_results["n_params"],
                "final_loss": epitope_results["loss_history"][-1] if epitope_results["loss_history"] else 0.0,
                "training_size": epitope_results["training_size"],
                "n_generated": epitope_results["generated"],
                "n_unique": epitope_results["unique"],
                "n_valid": epitope_results["valid"],
            },
        }

        # Dados — top-K candidatos de cada pipeline
        envelope["data"] = {
            "aso": {
                "baseline": aso_results["baseline"],
                "top_candidates": aso_results["top_k"],
                "all_valid_count": aso_results["valid"],
                "loss_history_summary": {
                    "first": aso_results["loss_history"][0] if aso_results["loss_history"] else 0.0,
                    "last": aso_results["loss_history"][-1] if aso_results["loss_history"] else 0.0,
                    "min": min(aso_results["loss_history"]) if aso_results["loss_history"] else 0.0,
                },
            },
            "epitope": {
                "original_epitopes": epitope_results.get("epitope_originals", []),
                "top_candidates": epitope_results["top_k"],
                "all_valid_count": epitope_results["valid"],
                "loss_history_summary": {
                    "first": epitope_results["loss_history"][0] if epitope_results["loss_history"] else 0.0,
                    "last": epitope_results["loss_history"][-1] if epitope_results["loss_history"] else 0.0,
                    "min": min(epitope_results["loss_history"]) if epitope_results["loss_history"] else 0.0,
                },
            },
            "config": {
                "n_diffusion_steps": config.n_diffusion_steps,
                "d_model": config.d_model,
                "n_heads": config.n_heads,
                "n_layers": config.n_layers,
                "n_epochs": config.n_epochs,
                "learning_rate": config.learning_rate,
                "n_samples": config.n_samples,
                "seed": config.seed,
                "device": device,
            },
        }

        # Resumo
        n_aso_valid = aso_results["valid"]
        n_epi_valid = epitope_results["valid"]
        best_aso_score = aso_results["top_k"][0]["score"] if aso_results["top_k"] else 0.0
        best_aso_seq = aso_results["top_k"][0]["sequence"] if aso_results["top_k"] else "N/A"
        best_epi_score = epitope_results["top_k"][0]["score"] if epitope_results["top_k"] else 0.0

        envelope["summary"]["conclusion"] = (
            f"Difusao discreta gerou {n_aso_valid} variantes validas de ASO "
            f"e {n_epi_valid} variantes de epitopos. "
            f"Melhor ASO: {best_aso_seq} (score {best_aso_score:.4f}). "
            f"Melhor epitopo score: {best_epi_score:.4f}."
        )
        envelope["summary"]["key_metrics"] = {
            "aso_valid_candidates": n_aso_valid,
            "epitope_valid_candidates": n_epi_valid,
            "best_aso_score": best_aso_score,
            "best_epitope_score": best_epi_score,
            "diffusion_steps": config.n_diffusion_steps,
            "model_dim": config.d_model,
        }

    envelope["runtime_seconds"] = timer.elapsed

    # Salva resultados
    output_dir = Path(__file__).resolve().parent / "results"
    output_path = write_result(envelope, output_dir=output_dir)
    print(f"\n[06_evodiff] Resultado salvo em {output_path}")
    print(f"[06_evodiff] Tempo total: {timer.elapsed:.1f}s")

    return envelope


if __name__ == "__main__":
    main()
