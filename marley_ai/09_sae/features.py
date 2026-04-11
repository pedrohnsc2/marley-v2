"""Interpretacao das features descobertas pelo Sparse Autoencoder.

Apos treinamento, cada neuronio do hidden layer representa uma "feature"
aprendida. Este modulo identifica o que cada feature significa biologicamente:

  1. Correlacao com propriedades fisico-quimicas conhecidas
  2. Quais sequencias ativam cada feature mais fortemente
  3. O que distingue epitopos de peptideos aleatorios
  4. Descricoes em linguagem natural para cada feature
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import importlib as _il

# Importacao via importlib — prefixo numerico "09" impede import direto
_encoder = _il.import_module("marley_ai.09_sae.encoder")
get_property_names = _encoder.get_property_names


# ---------------------------------------------------------------------------
# Estrutura de uma feature interpretada
# ---------------------------------------------------------------------------

@dataclass
class InterpretedFeature:
    """Uma feature descoberta pelo SAE com interpretacao biologica.

    Campos:
        index: indice do neuronio no hidden layer
        label: descricao curta da feature (ex: "hydrophobic_patch")
        description: explicacao biologica em portugues
        mean_activation_epitopes: ativacao media nos epitopos
        mean_activation_controls: ativacao media nos controles
        selectivity_ratio: razao epitopo/controle (>1 = preferencia por epitopos)
        top_correlated_properties: propriedades do input mais correlacionadas
        top_activating_sequences: sequencias que mais ativam esta feature
    """

    index: int
    label: str
    description: str
    mean_activation_epitopes: float
    mean_activation_controls: float
    selectivity_ratio: float
    top_correlated_properties: list[tuple[str, float]]
    top_activating_sequences: list[tuple[str, float]]


def correlate_with_properties(
    X: np.ndarray,
    activations: np.ndarray,
) -> np.ndarray:
    """Calcula correlacao de Pearson entre features do SAE e propriedades do input.

    Para cada neuronio do hidden layer, mede quanto sua ativacao
    esta correlacionada com cada propriedade fisico-quimica do input.
    Isto revela O QUE cada feature "aprendeu a detectar".

    Args:
        X: inputs normalizados (N, input_dim)
        activations: ativacoes do SAE (N, hidden_dim)

    Returns:
        Matriz de correlacao (hidden_dim, input_dim) — valores em [-1, 1]
    """
    hidden_dim = activations.shape[1]
    input_dim = X.shape[1]
    corr = np.zeros((hidden_dim, input_dim))

    for h in range(hidden_dim):
        z_h = activations[:, h]
        # Neuronio morto — correlacao indefinida
        if np.std(z_h) < 1e-12:
            continue
        for d in range(input_dim):
            x_d = X[:, d]
            if np.std(x_d) < 1e-12:
                continue
            # Pearson = cov(z,x) / (std(z) * std(x))
            corr[h, d] = np.corrcoef(z_h, x_d)[0, 1]

    return corr


def identify_selective_features(
    activations: np.ndarray,
    epitope_mask: np.ndarray,
    min_activation: float = 0.01,
) -> list[dict]:
    """Identifica features que distinguem epitopos de controles.

    Para cada neuronio, compara sua ativacao media em epitopos vs
    peptideos aleatorios. Features com alta selectivity_ratio sao
    candidatas a representar propriedades biologicamente relevantes
    para imunogenicidade.

    Args:
        activations: (N, hidden_dim) ativacoes do SAE
        epitope_mask: (N,) booleano — True para epitopos, False para controles
        min_activation: ativacao minima para considerar o neuronio "ativo"

    Returns:
        Lista de dicts com indice, ativacao media por grupo, e razao de seletividade.
        Ordenada por selectivity_ratio decrescente.
    """
    n_epitopes = np.sum(epitope_mask)
    n_controls = np.sum(~epitope_mask)

    if n_epitopes == 0 or n_controls == 0:
        return []

    mean_epi = np.mean(activations[epitope_mask], axis=0)
    mean_ctrl = np.mean(activations[~epitope_mask], axis=0)

    results = []
    for h in range(activations.shape[1]):
        # Ignora neuronios mortos
        if mean_epi[h] < min_activation and mean_ctrl[h] < min_activation:
            continue

        # Razao de seletividade — evita divisao por zero
        denominator = max(mean_ctrl[h], 1e-8)
        ratio = mean_epi[h] / denominator

        results.append({
            "index": h,
            "mean_epitope": round(float(mean_epi[h]), 6),
            "mean_control": round(float(mean_ctrl[h]), 6),
            "selectivity_ratio": round(float(ratio), 4),
        })

    # Ordena por seletividade (preferencia por epitopos)
    results.sort(key=lambda r: r["selectivity_ratio"], reverse=True)
    return results


def _generate_label(
    corr_row: np.ndarray,
    property_names: list[str],
    top_n: int = 3,
) -> tuple[str, str, list[tuple[str, float]]]:
    """Gera rotulo e descricao para uma feature com base em correlacoes.

    A feature e rotulada pela propriedade com maior correlacao absoluta.
    A descricao explica o significado biologico da combinacao de
    propriedades correlacionadas.

    Mapeamento propriedade -> significado biologico:
      - hydrophobicity: regioes transmembrana, nucleos proteicos, interfaces
      - charge: interacoes eletrostaticas, ligacao ao DNA/RNA
      - mol_weight: tamanho dos residuos, acessibilidade
      - volume: empacotamento, cavidades cataliticas
      - flexibility: loops, desordem, sitios de modificacao
      - aromaticity: stacking pi-pi, ligacao a ligantes aromaticos
    """
    # Mapeamento de propriedades para descricoes biologicas
    bio_meaning = {
        "hydrophobicity": "hidrofobicidade (dominios transmembrana, nucleos)",
        "charge": "carga eletrica (interacoes ionicas, ligacao acidos nucleicos)",
        "mol_weight": "peso molecular dos residuos (tamanho/acessibilidade)",
        "volume": "volume dos residuos (empacotamento, cavidades cataliticas)",
        "flexibility": "flexibilidade da cadeia (loops, desordem intrinseca)",
        "aromaticity": "aromaticidade (stacking pi-pi, ligacao a ligantes)",
    }

    # Indices das propriedades com maior correlacao absoluta
    abs_corr = np.abs(corr_row)
    top_indices = np.argsort(abs_corr)[::-1][:top_n]

    top_props = [
        (property_names[i], round(float(corr_row[i]), 4))
        for i in top_indices
    ]

    # Rotulo baseado na propriedade dominante
    dominant_prop = property_names[top_indices[0]]
    dominant_corr = corr_row[top_indices[0]]

    # Extrai o componente biologico do nome da propriedade
    # formato: "window_stat_property" -> extrai "property"
    parts = dominant_prop.split("_")
    if len(parts) >= 3:
        window = parts[0]
        stat = parts[1]
        bio_prop = "_".join(parts[2:])
    else:
        window, stat, bio_prop = "full", "mean", dominant_prop

    # Direcao da correlacao
    direction = "alta" if dominant_corr > 0 else "baixa"

    # Rotulo compacto
    label = f"{bio_prop}_{window}_{stat}_{'pos' if dominant_corr > 0 else 'neg'}"

    # Descricao em portugues
    bio_desc = bio_meaning.get(bio_prop, bio_prop)
    description = (
        f"Detecta {direction} {bio_desc} na regiao {window}. "
        f"Correlacao dominante: {dominant_prop} (r={dominant_corr:.3f})."
    )

    return label, description, top_props


def interpret_features(
    X: np.ndarray,
    activations: np.ndarray,
    epitope_mask: np.ndarray,
    sequence_names: list[str],
    top_k: int = 10,
) -> list[InterpretedFeature]:
    """Pipeline completo de interpretacao das features do SAE.

    Para cada feature (neuronio) ativo:
      1. Calcula correlacao com propriedades fisico-quimicas
      2. Identifica sequencias que mais a ativam
      3. Compara ativacao em epitopos vs controles
      4. Gera rotulo e descricao em linguagem natural

    Args:
        X: inputs normalizados (N, input_dim)
        activations: ativacoes do SAE (N, hidden_dim)
        epitope_mask: (N,) booleano — True para epitopos
        sequence_names: nomes das sequencias (para rotular top activators)
        top_k: numero de features mais seletivas a interpretar

    Returns:
        Lista de InterpretedFeature, ordenada por seletividade
    """
    property_names = get_property_names()

    # Correlacao features x propriedades
    corr_matrix = correlate_with_properties(X, activations)

    # Features seletivas para epitopos
    selective = identify_selective_features(activations, epitope_mask)

    # Interpretar as top_k mais seletivas
    interpreted: list[InterpretedFeature] = []

    for feat_info in selective[:top_k]:
        h = feat_info["index"]

        # Rotulo e descricao baseados em correlacoes
        label, description, top_props = _generate_label(
            corr_matrix[h], property_names
        )

        # Top sequencias que ativam esta feature
        activation_col = activations[:, h]
        top_seq_indices = np.argsort(activation_col)[::-1][:5]
        top_seqs = [
            (sequence_names[i], round(float(activation_col[i]), 4))
            for i in top_seq_indices
        ]

        interpreted.append(InterpretedFeature(
            index=h,
            label=label,
            description=description,
            mean_activation_epitopes=feat_info["mean_epitope"],
            mean_activation_controls=feat_info["mean_control"],
            selectivity_ratio=feat_info["selectivity_ratio"],
            top_correlated_properties=top_props,
            top_activating_sequences=top_seqs,
        ))

    return interpreted


def feature_importance_ranking(
    activations: np.ndarray,
    epitope_mask: np.ndarray,
) -> np.ndarray:
    """Calcula importancia de cada feature por diferenca de ativacao.

    Importancia = |mean_epitope - mean_control| / (std_pooled + eps)

    Analoga a estatistica t, mas simplificada. Features com alta
    importancia sao candidatas a biomarcadores de imunogenicidade.

    Args:
        activations: (N, hidden_dim)
        epitope_mask: (N,) booleano

    Returns:
        Vetor (hidden_dim,) de importancias — maior = mais discriminativo
    """
    epi = activations[epitope_mask]
    ctrl = activations[~epitope_mask]

    mean_diff = np.abs(np.mean(epi, axis=0) - np.mean(ctrl, axis=0))
    # Desvio padrao combinado (pooled)
    std_pooled = np.sqrt(
        (np.var(epi, axis=0) + np.var(ctrl, axis=0)) / 2.0 + 1e-12
    )

    return mean_diff / std_pooled


def summarize_features(
    features: list[InterpretedFeature],
) -> str:
    """Gera resumo textual das features descobertas.

    Formato legivel para inclusao no relatorio/envelope JSON.

    Args:
        features: lista de features interpretadas

    Returns:
        String formatada com o resumo
    """
    lines = [
        "=" * 65,
        "  FEATURES DESCOBERTAS PELO SPARSE AUTOENCODER",
        "=" * 65,
        "",
    ]

    for i, feat in enumerate(features, 1):
        lines.append(f"  Feature #{i} (neuronio {feat.index})")
        lines.append(f"    Rotulo: {feat.label}")
        lines.append(f"    {feat.description}")
        lines.append(
            f"    Ativacao media — epitopos: {feat.mean_activation_epitopes:.4f} | "
            f"controles: {feat.mean_activation_controls:.4f}"
        )
        lines.append(f"    Razao de seletividade: {feat.selectivity_ratio:.2f}x")
        lines.append(f"    Top propriedades correlacionadas:")
        for prop_name, corr_val in feat.top_correlated_properties[:3]:
            lines.append(f"      - {prop_name}: r={corr_val:.3f}")
        lines.append(f"    Top sequencias ativadoras:")
        for seq_name, act_val in feat.top_activating_sequences[:3]:
            lines.append(f"      - {seq_name}: ativacao={act_val:.4f}")
        lines.append("")

    lines.append("=" * 65)
    return "\n".join(lines)
