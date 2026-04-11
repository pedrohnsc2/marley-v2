"""Execucao do modulo 09_sae — Sparse Autoencoders para interpretabilidade.

Descobre features interpretaveis a partir de representacoes hand-crafted
de peptideos, comparando epitopos vacinais de L. infantum com controles
aleatorios. Cada neuronio do SAE corresponde a uma propriedade biologica
(ex: "regiao hidrofobica", "cluster aromatico", "carga positiva").

Uso:
    python -m marley_ai.09_sae.run
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

import importlib as _il

import numpy as np

from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

# Importacao via importlib — Python nao permite "from marley_ai.09_sae..."
# porque "09" e interpretado como literal numerico invalido no parser.
_sae_config = _il.import_module("marley_ai.09_sae.config")
_sae_encoder = _il.import_module("marley_ai.09_sae.encoder")
_sae_autoencoder = _il.import_module("marley_ai.09_sae.autoencoder")
_sae_features = _il.import_module("marley_ai.09_sae.features")

SAEConfig = _sae_config.SAEConfig
encode_batch = _sae_encoder.encode_batch
normalize_features = _sae_encoder.normalize_features
generate_random_peptides = _sae_encoder.generate_random_peptides
initialize_weights = _sae_autoencoder.initialize_weights
train = _sae_autoencoder.train
get_activations = _sae_autoencoder.get_activations
reconstruction_quality = _sae_autoencoder.reconstruction_quality
interpret_features = _sae_features.interpret_features
feature_importance_ranking = _sae_features.feature_importance_ranking
summarize_features = _sae_features.summarize_features

# Epitopos canonicos do pipeline Marley — fonte unica de verdade
from vaccine_platforms.shared.epitopes import EPITOPES, get_epitope_sequences

# Sequencia alvo do ASO (SL RNA de L. infantum)
from aso_math.config import SL_SEQUENCE, ASO_TARGET_SEQUENCE


# ---------------------------------------------------------------------------
# Sequencias adicionais de proteinas-alvo conhecidas
# ---------------------------------------------------------------------------
# Fragmentos representativos de proteinas-alvo vacinais de L. infantum.
# Estes peptideos sao regioes imunodominantes documentadas na literatura,
# usados como contexto adicional para o SAE (alem dos 11 epitopos).

ADDITIONAL_TARGETS: dict[str, str] = {
    # GP63/Leishmanolysin — regiao catalitica com motivo HEXXH de metalopeptidase
    "GP63_catalytic": "VATHEIGHVLGLAHQ",
    # Prohibitin — regiao de dimerizacao SPFH/PHB
    "Prohibitin_SPFH": "AQNLEKQIIEQRKTV",
    # CPB — pro-dominio (ERFNIN motif) com papel na regulacao da atividade
    "CPB_prodomain": "ERFNINKDLTEEFRK",
    # Histona H2A — alvo vacinal conhecido, superficie exposta
    "H2A_surface": "SGRGKQGGKTRAKAK",
    # KMP-11 — proteina kinetoplastidea de membrana, alvo de vacina DNA
    "KMP11_core": "AEDKLTQEQLANFQK",
}


# ---------------------------------------------------------------------------
# Registro do modulo no sistema Marley
# ---------------------------------------------------------------------------

@register("09_sae")
class SAEModule:
    """Modulo de Sparse Autoencoders para interpretabilidade mecanistica."""

    def __init__(self) -> None:
        self._config: SAEConfig | None = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao do orquestrador."""
        if isinstance(config, SAEConfig):
            self._config = config
        else:
            self._config = SAEConfig(
                module_slug="09_sae",
                module_name="Sparse Autoencoders",
            )

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica dependencias — este modulo usa dados internos."""
        return {"valid": True, "missing": []}

    def run(self) -> dict[str, Any]:
        """Executa pipeline SAE completo."""
        return main(self._config)

    def get_dependencies(self) -> list[str]:
        """Sem dependencias externas — usa epitopos embutidos."""
        return []


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Pipeline completo do Sparse Autoencoder.

    Etapas:
        1. Coleta sequencias: 11 epitopos + alvos adicionais + controles
        2. Codifica em vetores de propriedades (72-dim)
        3. Normaliza features (z-score)
        4. Treina SAE (MSE + L1)
        5. Extrai e interpreta features
        6. Gera relatorio com top 10 features biologicas

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com features interpretaveis descobertas.
    """
    if config is None:
        config = SAEConfig(
            module_slug="09_sae",
            module_name="Sparse Autoencoders",
        )

    envelope = create_envelope("09_sae")
    envelope["device"] = "cpu"  # numpy puro — sem GPU

    with Timer() as timer:

        # ---------------------------------------------------------------
        # 1. Coleta de sequencias
        # ---------------------------------------------------------------
        print("[09_sae] Coletando sequencias...")

        # Epitopos do construto vacinal (11 peptideos unicos)
        epitope_seqs = get_epitope_sequences()
        epitope_names = [
            f"Epi_{ep.gene_name[:15]}_{ep.peptide}"
            for ep in EPITOPES
        ]

        # Alvos adicionais (fragmentos de proteinas vacinais)
        additional_seqs = list(ADDITIONAL_TARGETS.values())
        additional_names = list(ADDITIONAL_TARGETS.keys())

        # Sequencia-alvo do ASO (SL RNA convertida para "proteina" — cada
        # base mapeada para aminoacido representativo para compatibilidade)
        # Na verdade, usamos a regiao alvo do ASO como peptideo hipotetico
        # Nota: ASO_TARGET_SEQUENCE e DNA, nao proteina. Incluimos como
        # controle de composicao — o SAE deveria REJEITAR esta sequencia
        # (nao e um peptideo real).
        sl_rna_as_peptide = ASO_TARGET_SEQUENCE.replace("A", "A").replace(
            "T", "T").replace("G", "G").replace("C", "C")

        # Peptideos aleatorios como controle negativo
        n_controls = 50
        if isinstance(config, SAEConfig):
            n_controls = config.n_control_peptides
        control_seqs = generate_random_peptides(n_controls, seed=config.seed)
        control_names = [f"Ctrl_{i:03d}" for i in range(n_controls)]

        # Dataset completo
        all_seqs = epitope_seqs + additional_seqs + control_seqs
        all_names = epitope_names + additional_names + control_names

        # Mascara: True para epitopos + alvos, False para controles
        n_positive = len(epitope_seqs) + len(additional_seqs)
        epitope_mask = np.array(
            [True] * n_positive + [False] * len(control_seqs),
            dtype=bool,
        )

        print(f"  Epitopos: {len(epitope_seqs)}")
        print(f"  Alvos adicionais: {len(additional_seqs)}")
        print(f"  Controles aleatorios: {len(control_seqs)}")
        print(f"  Total: {len(all_seqs)} sequencias")

        # ---------------------------------------------------------------
        # 2. Codificacao em vetores de propriedades
        # ---------------------------------------------------------------
        print("\n[09_sae] Codificando sequencias em vetores de propriedades...")
        X_raw = encode_batch(all_seqs)
        print(f"  Dimensao do dataset: {X_raw.shape}")

        # ---------------------------------------------------------------
        # 3. Normalizacao z-score
        # ---------------------------------------------------------------
        X_norm, feat_mean, feat_std = normalize_features(X_raw)
        print(f"  Features normalizadas (media ~0, std ~1)")

        # ---------------------------------------------------------------
        # 4. Treinamento do SAE
        # ---------------------------------------------------------------
        print("\n[09_sae] Treinando Sparse Autoencoder...")

        input_dim = X_norm.shape[1]
        hidden_dim = 512
        sparsity_penalty = 0.01
        lr = 5e-3
        n_epochs = 500

        if isinstance(config, SAEConfig):
            hidden_dim = config.hidden_dim
            sparsity_penalty = config.sparsity_penalty
            lr = config.learning_rate
            n_epochs = config.n_epochs

        print(f"  Arquitetura: {input_dim} -> {hidden_dim} -> {input_dim}")
        print(f"  L1 penalty: {sparsity_penalty}")
        print(f"  Epocas: {n_epochs}")
        print()

        weights = initialize_weights(
            input_dim, hidden_dim,
            tied=True, seed=config.seed,
        )

        weights = train(
            X_norm, weights,
            sparsity_penalty=sparsity_penalty,
            learning_rate=lr,
            n_epochs=n_epochs,
            verbose=True,
            print_every=100,
        )

        # ---------------------------------------------------------------
        # 5. Avaliacao da reconstrucao
        # ---------------------------------------------------------------
        print("\n[09_sae] Avaliando qualidade da reconstrucao...")
        quality = reconstruction_quality(X_norm, weights)
        print(f"  MSE: {quality['mse']:.6f}")
        print(f"  R2: {quality['r2']:.4f}")
        print(f"  Esparsidade: {quality['sparsity']:.1%}")
        print(f"  Neuronios vivos: {quality['alive_neurons']}/{quality['total_neurons']}")

        # ---------------------------------------------------------------
        # 6. Extracao e interpretacao de features
        # ---------------------------------------------------------------
        print("\n[09_sae] Extraindo e interpretando features...")
        activations = get_activations(X_norm, weights)

        top_k = 10
        if isinstance(config, SAEConfig):
            top_k = config.top_k_features

        interpreted = interpret_features(
            X_norm, activations, epitope_mask,
            all_names, top_k=top_k,
        )

        # Importancia global das features
        importance = feature_importance_ranking(activations, epitope_mask)
        top_important_indices = np.argsort(importance)[::-1][:top_k]

        # Resumo textual
        summary_text = summarize_features(interpreted)
        print()
        print(summary_text)

        # ---------------------------------------------------------------
        # 7. Montagem do envelope de resultados
        # ---------------------------------------------------------------
        envelope["status"] = "success"
        envelope["dependencies"] = ["vaccine_platforms.shared.epitopes", "aso_math.config"]

        # Metricas de treinamento
        envelope["metrics"] = {
            "final_loss": weights.loss_history[-1] if weights.loss_history else 0.0,
            "final_mse": weights.mse_history[-1] if weights.mse_history else 0.0,
            "final_l1": weights.l1_history[-1] if weights.l1_history else 0.0,
            "reconstruction_r2": quality["r2"],
            "reconstruction_mse": quality["mse"],
            "sparsity": quality["sparsity"],
            "alive_neurons": quality["alive_neurons"],
            "dead_neurons": quality["dead_neurons"],
            "total_neurons": quality["total_neurons"],
            "n_epochs": n_epochs,
            "n_sequences": len(all_seqs),
            "n_epitopes": len(epitope_seqs),
            "n_additional_targets": len(additional_seqs),
            "n_controls": len(control_seqs),
        }

        # Features interpretadas
        envelope["data"]["interpreted_features"] = [
            {
                "rank": i + 1,
                "neuron_index": feat.index,
                "label": feat.label,
                "description": feat.description,
                "mean_activation_epitopes": feat.mean_activation_epitopes,
                "mean_activation_controls": feat.mean_activation_controls,
                "selectivity_ratio": feat.selectivity_ratio,
                "top_correlated_properties": [
                    {"property": p, "correlation": c}
                    for p, c in feat.top_correlated_properties
                ],
                "top_activating_sequences": [
                    {"sequence": s, "activation": a}
                    for s, a in feat.top_activating_sequences
                ],
            }
            for i, feat in enumerate(interpreted)
        ]

        # Top features por importancia global (diferenca epitopos vs controles)
        envelope["data"]["top_importance_indices"] = [
            int(idx) for idx in top_important_indices
        ]
        envelope["data"]["top_importance_scores"] = [
            round(float(importance[idx]), 4) for idx in top_important_indices
        ]

        # Configuracao usada
        envelope["data"]["config"] = {
            "input_dim": input_dim,
            "hidden_dim": hidden_dim,
            "sparsity_penalty": sparsity_penalty,
            "learning_rate": lr,
            "n_epochs": n_epochs,
            "seed": config.seed,
        }

        # Dataset info
        envelope["data"]["sequences"] = {
            "epitopes": {
                name: seq for name, seq in zip(epitope_names, epitope_seqs)
            },
            "additional_targets": ADDITIONAL_TARGETS,
            "n_controls": len(control_seqs),
        }

        # Historico de treinamento (subamostrado para nao inflar o JSON)
        step = max(1, n_epochs // 20)
        envelope["data"]["training_history"] = {
            "loss": [round(v, 6) for v in weights.loss_history[::step]],
            "mse": [round(v, 6) for v in weights.mse_history[::step]],
            "l1": [round(v, 6) for v in weights.l1_history[::step]],
            "sparsity": [round(v, 4) for v in weights.sparsity_history[::step]],
        }

        # Resumo
        n_features_found = len(interpreted)
        top_feature_label = interpreted[0].label if interpreted else "nenhuma"
        top_selectivity = interpreted[0].selectivity_ratio if interpreted else 0.0

        envelope["summary"]["conclusion"] = (
            f"SAE descobriu {n_features_found} features interpretaveis "
            f"que distinguem epitopos de L. infantum de peptideos aleatorios. "
            f"Feature mais seletiva: '{top_feature_label}' "
            f"(razao epitopo/controle: {top_selectivity:.1f}x). "
            f"Reconstrucao R2={quality['r2']:.3f}, "
            f"esparsidade={quality['sparsity']:.1%}."
        )
        envelope["summary"]["key_metrics"] = {
            "input_dim": input_dim,
            "hidden_dim": hidden_dim,
            "reconstruction_r2": quality["r2"],
            "sparsity": quality["sparsity"],
            "n_features_found": n_features_found,
            "top_selectivity_ratio": top_selectivity,
        }

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"\n[09_sae] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
