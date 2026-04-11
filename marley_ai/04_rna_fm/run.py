"""Execucao do modulo 04_rna_fm — Analise de RNA com deep learning.

Pipeline completo:
    1. Treina encoder customizado de RNA (Masked Nucleotide Modeling)
    2. Prediz estrutura secundaria do SL RNA (livre e com ASO)
    3. Analisa conservacao cross-species e unicidade do SL RNA
    4. Gera embeddings e compara SL RNA vs controles humanos/caninos

Uso:
    python -m marley_ai.04_rna_fm.run
"""

from __future__ import annotations

import random
from typing import Any

import numpy as np
import torch

from marley_ai.config import AIModuleConfig, detect_device
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

from .config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    CANINE_CONTROL_RNA,
    ENCODER_PARAMS,
    HUMAN_SNRNA,
    RNAFMConfig,
    SL_SEQUENCES,
)
from .model import RNAEncoder, train_rna_encoder
from .structure import (
    analyze_aso_structural_impact,
    predict_structure,
)
from .sl_analysis import run_sl_analysis


# ---------------------------------------------------------------------------
# Registro do modulo no orquestrador
# ---------------------------------------------------------------------------

@register("04_rna_fm")
class RNAFM:
    """Modulo de analise de RNA com encoder customizado e Nussinov."""

    def __init__(self) -> None:
        self.config: RNAFMConfig | None = None
        self.model: RNAEncoder | None = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao do modulo."""
        self.config = config

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica dependencias — este modulo e auto-contido."""
        return {"valid": True, "missing": []}

    def run(self) -> dict[str, Any]:
        """Executa pipeline completo."""
        return main(self.config)

    def get_dependencies(self) -> list[str]:
        """Nenhuma dependencia de outros modulos."""
        return []


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Executa pipeline completo de analise de RNA para o SL RNA de Leishmania.

    Etapas:
        1. Detecta device e configura reproducibilidade
        2. Treina encoder customizado de RNA (MNM — ~100 epochs)
        3. Prediz estrutura secundaria (Nussinov) — livre e com ASO
        4. Analisa conservacao, embeddings e unicidade posicional
        5. Grava envelope JSON com todos os resultados

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com resultados completos da analise de RNA.
    """
    envelope = create_envelope("04_rna_fm")

    # --- Detecta device ---
    device = detect_device()
    envelope["device"] = device
    print(f"[04_rna_fm] Device: {device}")

    # --- Reproducibilidade ---
    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if device == "mps":
        torch.mps.manual_seed(seed)

    with Timer() as timer:
        try:
            # ==============================================================
            # ETAPA 1: Treinar encoder customizado de RNA
            # ==============================================================
            print("\n[04_rna_fm] === ETAPA 1: Treinamento do encoder de RNA ===")

            # Combina todas as sequencias para treinamento
            training_sequences: dict[str, str] = {}
            training_sequences.update(SL_SEQUENCES)
            training_sequences.update(HUMAN_SNRNA)
            training_sequences.update(CANINE_CONTROL_RNA)

            model, losses = train_rna_encoder(
                sequences=training_sequences,
                device=device,
                num_epochs=ENCODER_PARAMS["num_epochs"],
                learning_rate=ENCODER_PARAMS["learning_rate"],
                verbose=True,
            )

            envelope["metrics"]["training"] = {
                "n_sequences_original": len(training_sequences),
                "n_epochs": len(losses),
                "initial_loss": round(losses[0], 4),
                "final_loss": round(losses[-1], 4),
                "loss_reduction_pct": round(
                    (1 - losses[-1] / losses[0]) * 100 if losses[0] > 0 else 0, 2
                ),
                "device": device,
                "embed_dim": ENCODER_PARAMS["embed_dim"],
                "num_layers": ENCODER_PARAMS["num_layers"],
                "num_heads": ENCODER_PARAMS["num_heads"],
            }

            print(f"\n  Treinamento concluido: loss {losses[0]:.4f} -> {losses[-1]:.4f}")

            # ==============================================================
            # ETAPA 2: Predicao de estrutura secundaria
            # ==============================================================
            print("\n[04_rna_fm] === ETAPA 2: Estrutura secundaria (Nussinov) ===")

            structure_results = analyze_aso_structural_impact()

            print(f"  SL RNA livre:")
            print(f"    Sequencia:   {structure_results['sl_sequence']}")
            print(f"    Estrutura:   {structure_results['free_structure']['dot_bracket']}")
            print(f"    Pares:       {structure_results['free_structure']['n_pairs']}")
            print(f"  SL RNA + ASO (pos {ASO_TARGET_START}-{ASO_TARGET_END - 1} bloqueadas):")
            print(f"    Estrutura:   {structure_results['bound_structure']['dot_bracket']}")
            print(f"    Pares:       {structure_results['bound_structure']['n_pairs']}")
            print(f"  Impacto:")
            impact = structure_results["structural_impact"]
            print(f"    Pares disruptados:   {impact['pairs_disrupted']}")
            print(f"    Rearranjos:          {impact['pairs_new_rearrangement']}")
            print(f"    Score change:        {impact['score_change']}")

            envelope["data"]["structure"] = _serialize_structure(structure_results)

            # ==============================================================
            # ETAPA 3: Analise completa do SL RNA
            # ==============================================================
            print("\n[04_rna_fm] === ETAPA 3: Analise de unicidade do SL RNA ===")

            sl_results = run_sl_analysis(model, device=device)

            # Conservacao
            cons = sl_results["conservation"]
            print(f"\n  Conservacao cross-species:")
            print(f"    Identidade media: {cons['mean_identity'] * 100:.1f}%")
            print(f"    Posicoes variantes: {cons['n_variant_positions']}/{len(SL_SEQUENCES['L_infantum'])}")

            # Embeddings
            emb = sl_results["embedding_distances"]
            print(f"\n  Distancias no espaco de embeddings:")
            print(f"    Separacao inter/intra: {emb['separation_ratio']:.2f}x")
            print(f"    SL -> humano: {emb['mean_sl_to_human_distance']:.4f}")
            print(f"    SL -> canino: {emb['mean_sl_to_canine_distance']:.4f}")

            # Unicidade
            uniq = sl_results["positional_uniqueness"]
            print(f"\n  Unicidade posicional:")
            print(f"    Regiao ASO: {uniq['mean_aso_region_uniqueness']:.4f}")
            print(f"    Fora ASO:   {uniq['mean_non_aso_uniqueness']:.4f}")
            print(f"    Enrichment: {uniq['aso_uniqueness_enrichment']:.2f}x")

            envelope["data"]["conservation"] = cons
            envelope["data"]["embedding_distances"] = _serialize_embedding_data(emb)
            envelope["data"]["positional_uniqueness"] = _serialize_uniqueness(uniq)

            # ==============================================================
            # ETAPA 4: Resumo e conclusao
            # ==============================================================
            envelope["status"] = "complete"
            envelope["summary"]["conclusion"] = _build_conclusion(
                structure_results, cons, emb, uniq
            )
            envelope["summary"]["key_metrics"] = {
                "sl_length": len(SL_SEQUENCES["L_infantum"]),
                "model": "custom_rna_encoder",
                "embed_dim": ENCODER_PARAMS["embed_dim"],
                "training_loss_final": round(losses[-1], 4),
                "free_structure_pairs": structure_results["free_structure"]["n_pairs"],
                "bound_structure_pairs": structure_results["bound_structure"]["n_pairs"],
                "pairs_disrupted_by_aso": impact["pairs_disrupted"],
                "cross_species_identity": cons["mean_identity"],
                "embedding_separation_ratio": emb["separation_ratio"],
                "aso_region_uniqueness": uniq["mean_aso_region_uniqueness"],
            }
            envelope["warnings"] = []

        except Exception as exc:
            envelope["status"] = "failed"
            envelope["summary"]["conclusion"] = f"Falha na execucao: {exc}"
            envelope["warnings"].append(str(exc))
            raise

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"\n[04_rna_fm] Resultado salvo em {output_path}")
    print(f"[04_rna_fm] Runtime: {timer.elapsed:.2f}s")

    return envelope


# ---------------------------------------------------------------------------
# Funcoes auxiliares de serializacao (tensores -> listas para JSON)
# ---------------------------------------------------------------------------

def _serialize_structure(results: dict) -> dict:
    """Converte resultados de estrutura para formato serializavel (JSON).

    Substitui tuplas e objetos numpy por listas e floats nativos.
    """
    serialized = {}
    for key, value in results.items():
        if key in ("free_structure", "bound_structure"):
            sub = dict(value)
            if "pairs" in sub:
                sub["pairs"] = [list(p) for p in sub["pairs"]]
            serialized[key] = sub
        elif key == "structural_impact":
            serialized[key] = {k: _to_native(v) for k, v in value.items()}
        else:
            serialized[key] = _to_native(value)
    return serialized


def _serialize_embedding_data(data: dict) -> dict:
    """Converte dados de embedding para formato JSON-serializavel."""
    serialized = {}
    for key, value in data.items():
        if key == "selected_pairs":
            serialized[key] = [
                {k: _to_native(v) for k, v in pair.items()}
                for pair in value
            ]
        else:
            serialized[key] = _to_native(value)
    return serialized


def _serialize_uniqueness(data: dict) -> dict:
    """Converte dados de unicidade para formato JSON-serializavel."""
    serialized = {}
    for key, value in data.items():
        if key == "top_unique_positions_in_aso":
            serialized[key] = [[pos, _to_native(score)] for pos, score in value]
        else:
            serialized[key] = _to_native(value)
    return serialized


def _to_native(value: Any) -> Any:
    """Converte valores numpy/torch para tipos nativos Python (JSON-safe)."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, torch.Tensor):
        return value.tolist()
    if isinstance(value, dict):
        return {k: _to_native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_native(v) for v in value]
    return value


# ---------------------------------------------------------------------------
# Conclusao sintetica
# ---------------------------------------------------------------------------

def _build_conclusion(
    structure: dict,
    conservation: dict,
    embeddings: dict,
    uniqueness: dict,
) -> str:
    """Gera conclusao sintetica integrando todos os resultados.

    Args:
        structure: Resultados de estrutura secundaria.
        conservation: Resultados de conservacao cross-species.
        embeddings: Resultados de distancias no espaco de embeddings.
        uniqueness: Resultados de unicidade posicional.

    Returns:
        String com conclusao integrativa.
    """
    impact = structure["structural_impact"]
    n_disrupted = impact["pairs_disrupted"]
    n_free = structure["free_structure"]["n_pairs"]
    identity = conservation["mean_identity"]
    sep_ratio = embeddings["separation_ratio"]
    aso_uniq = uniqueness["mean_aso_region_uniqueness"]

    return (
        f"Analise de RNA concluida para o SL RNA de L. infantum (39 nt). "
        f"ESTRUTURA: O ASO MRL-ASO-001 disrupta {n_disrupted} de {n_free} pares "
        f"de bases previstos pelo algoritmo de Nussinov, confirmando competicao "
        f"com a estrutura intramolecular. "
        f"CONSERVACAO: O SL RNA e {identity * 100:.1f}% conservado entre "
        f"{conservation['n_species']} especies de Leishmania, com "
        f"{conservation['n_variant_positions']} posicoes variantes. "
        f"UNICIDADE: No espaco de embeddings do encoder customizado, "
        f"SL RNAs se separam dos controles humanos/caninos com ratio "
        f"{sep_ratio:.2f}x. A regiao de ligacao do ASO tem score de "
        f"unicidade {aso_uniq:.2f}, confirmando a seletividade do alvo."
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
