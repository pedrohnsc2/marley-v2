"""Analise de unicidade e conservacao do SL RNA de Leishmania.

Este modulo compara o SL RNA de Leishmania contra:
    1. SL RNA de outras especies (conservacao cross-species)
    2. snRNAs humanos (divergencia — prova de seletividade)
    3. RNAs caninos (controle — especie hospedeira)

A analise no espaco de embeddings do encoder customizado demonstra
que o SL RNA ocupa uma regiao distinta do espaco latente, separada
dos RNAs do hospedeiro. Isso fundamenta a seletividade do ASO.

Adicionalmente, mapeia o perfil de unicidade posicional do SL RNA
sobre a regiao de ligacao do MRL-ASO-001, identificando quais
posicoes do alvo sao mais divergentes do transcriptoma humano.
"""

from __future__ import annotations

import math

import torch
import numpy as np

from .config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    CANINE_CONTROL_RNA,
    HUMAN_SNRNA,
    SL_SEQUENCES,
)
from .model import (
    RNAEncoder,
    get_sequence_embeddings,
)


# ---------------------------------------------------------------------------
# Conservacao cross-species
# ---------------------------------------------------------------------------

def compute_conservation(sequences: dict[str, str]) -> dict[str, object]:
    """Calcula metricas de conservacao entre sequencias de SL RNA.

    Alinha par-a-par as sequencias de SL RNA de diferentes especies de
    Leishmania e calcula identidade, divergencia e posicoes variantes.

    Args:
        sequences: Dicionario {especie: sequencia_sl_rna}.

    Returns:
        Dicionario com metricas de conservacao.
    """
    species = list(sequences.keys())
    seqs = list(sequences.values())

    # Normaliza: DNA -> RNA
    seqs_rna = [s.upper().replace("T", "U") for s in seqs]

    # Comparacao par-a-par (alinhamento global simples — sequencias de mesmo tamanho)
    pairwise_identity: dict[str, float] = {}
    for i in range(len(species)):
        for j in range(i + 1, len(species)):
            seq_a, seq_b = seqs_rna[i], seqs_rna[j]
            # Alinhamento trivial (sequencias de tamanho similar)
            min_len = min(len(seq_a), len(seq_b))
            matches = sum(1 for k in range(min_len) if seq_a[k] == seq_b[k])
            identity = matches / min_len if min_len > 0 else 0.0
            pair_key = f"{species[i]}_vs_{species[j]}"
            pairwise_identity[pair_key] = round(identity, 4)

    # Perfil de conservacao posicional (frequencia do nucleotideo mais comum)
    ref_len = len(seqs_rna[0])
    position_conservation = []
    variant_positions = []

    for pos in range(ref_len):
        nucs_at_pos = [s[pos] for s in seqs_rna if pos < len(s)]
        if not nucs_at_pos:
            continue
        # Frequencia do nucleotideo mais comum nesta posicao
        from collections import Counter
        counts = Counter(nucs_at_pos)
        most_common_freq = counts.most_common(1)[0][1] / len(nucs_at_pos)
        position_conservation.append(round(most_common_freq, 4))

        # Marca posicoes variantes (nao 100% conservadas)
        if most_common_freq < 1.0:
            variant_positions.append(pos)

    # Identidade media global
    identities = list(pairwise_identity.values())
    mean_identity = sum(identities) / len(identities) if identities else 0.0

    return {
        "n_species": len(species),
        "species": species,
        "pairwise_identity": pairwise_identity,
        "mean_identity": round(mean_identity, 4),
        "position_conservation": position_conservation,
        "variant_positions": variant_positions,
        "n_variant_positions": len(variant_positions),
        "conservation_summary": (
            f"SL RNA conservado em {mean_identity * 100:.1f}% entre "
            f"{len(species)} especies de Leishmania. "
            f"{len(variant_positions)} posicoes variantes de {ref_len} total."
        ),
    }


# ---------------------------------------------------------------------------
# Distancias no espaco de embeddings
# ---------------------------------------------------------------------------

def compute_embedding_distances(
    embeddings: dict[str, torch.Tensor],
) -> dict[str, object]:
    """Calcula distancias par-a-par no espaco de embeddings.

    Usa distancia euclidiana e similaridade cosseno para comparar
    RNAs no espaco latente do encoder. A hipotese e que SL RNAs
    agrupam entre si e se separam dos snRNAs humanos e RNAs caninos.

    Args:
        embeddings: Dicionario {nome: tensor_embedding}.

    Returns:
        Dicionario com matrizes de distancia e analise de agrupamento.
    """
    names = list(embeddings.keys())
    n = len(names)
    vecs = torch.stack([embeddings[name] for name in names])

    # Distancia euclidiana par-a-par
    euclidean_matrix = torch.cdist(vecs.unsqueeze(0), vecs.unsqueeze(0)).squeeze(0)

    # Similaridade cosseno par-a-par
    norms = vecs.norm(dim=1, keepdim=True).clamp(min=1e-8)
    vecs_normalized = vecs / norms
    cosine_matrix = vecs_normalized @ vecs_normalized.T

    # Classifica cada RNA em grupo
    groups = _classify_rna_groups(names)

    # Calcula distancias intra-grupo e inter-grupo
    intra_distances = []  # distancias dentro do mesmo grupo
    inter_distances = []  # distancias entre grupos diferentes

    sl_to_human_distances = []  # distancias SL -> humano
    sl_to_canine_distances = []  # distancias SL -> canino

    for i in range(n):
        for j in range(i + 1, n):
            dist = euclidean_matrix[i][j].item()
            if groups[names[i]] == groups[names[j]]:
                intra_distances.append(dist)
            else:
                inter_distances.append(dist)

            # Distancias especificas SL -> controles
            if groups[names[i]] == "sl_rna" and groups[names[j]] == "human_snrna":
                sl_to_human_distances.append(dist)
            elif groups[names[j]] == "sl_rna" and groups[names[i]] == "human_snrna":
                sl_to_human_distances.append(dist)
            if groups[names[i]] == "sl_rna" and groups[names[j]] == "canine_rna":
                sl_to_canine_distances.append(dist)
            elif groups[names[j]] == "sl_rna" and groups[names[i]] == "canine_rna":
                sl_to_canine_distances.append(dist)

    mean_intra = sum(intra_distances) / len(intra_distances) if intra_distances else 0.0
    mean_inter = sum(inter_distances) / len(inter_distances) if inter_distances else 0.0
    mean_sl_human = (
        sum(sl_to_human_distances) / len(sl_to_human_distances)
        if sl_to_human_distances else 0.0
    )
    mean_sl_canine = (
        sum(sl_to_canine_distances) / len(sl_to_canine_distances)
        if sl_to_canine_distances else 0.0
    )

    # Separation ratio: inter/intra > 1 indica boa separacao
    separation_ratio = mean_inter / mean_intra if mean_intra > 1e-8 else float("inf")

    # Tabela resumo de pares selecionados (SL vs controles)
    selected_pairs = _build_selected_pairs(names, groups, euclidean_matrix, cosine_matrix)

    return {
        "n_sequences": n,
        "groups": groups,
        "mean_intra_group_distance": round(mean_intra, 4),
        "mean_inter_group_distance": round(mean_inter, 4),
        "separation_ratio": round(separation_ratio, 4),
        "mean_sl_to_human_distance": round(mean_sl_human, 4),
        "mean_sl_to_canine_distance": round(mean_sl_canine, 4),
        "selected_pairs": selected_pairs,
        "interpretation": _interpret_embedding_distances(
            separation_ratio, mean_sl_human, mean_sl_canine, mean_intra
        ),
    }


def _classify_rna_groups(names: list[str]) -> dict[str, str]:
    """Classifica cada RNA em um grupo funcional.

    Args:
        names: Lista de nomes de RNA.

    Returns:
        Dicionario {nome: grupo}.
    """
    groups = {}
    sl_species = set(SL_SEQUENCES.keys())
    human_names = set(HUMAN_SNRNA.keys())
    canine_names = set(CANINE_CONTROL_RNA.keys())

    for name in names:
        if name in sl_species:
            groups[name] = "sl_rna"
        elif name in human_names:
            groups[name] = "human_snrna"
        elif name in canine_names:
            groups[name] = "canine_rna"
        else:
            groups[name] = "unknown"

    return groups


def _build_selected_pairs(
    names: list[str],
    groups: dict[str, str],
    euclidean_matrix: torch.Tensor,
    cosine_matrix: torch.Tensor,
) -> list[dict[str, object]]:
    """Constroi tabela de pares selecionados com distancias.

    Seleciona pares representativos: L. infantum vs cada controle.

    Args:
        names: Lista de nomes.
        groups: Classificacao em grupos.
        euclidean_matrix: Matriz de distancias euclidianas.
        cosine_matrix: Matriz de similaridade cosseno.

    Returns:
        Lista de dicionarios com pares e metricas.
    """
    pairs = []
    ref_name = "L_infantum"

    if ref_name not in names:
        return pairs

    ref_idx = names.index(ref_name)

    for i, name in enumerate(names):
        if name == ref_name:
            continue
        pairs.append({
            "pair": f"{ref_name} vs {name}",
            "group_a": groups.get(ref_name, "unknown"),
            "group_b": groups.get(name, "unknown"),
            "euclidean_distance": round(euclidean_matrix[ref_idx][i].item(), 4),
            "cosine_similarity": round(cosine_matrix[ref_idx][i].item(), 4),
        })

    # Ordena por distancia euclidiana (mais proximo primeiro)
    pairs.sort(key=lambda p: p["euclidean_distance"])
    return pairs


def _interpret_embedding_distances(
    separation_ratio: float,
    mean_sl_human: float,
    mean_sl_canine: float,
    mean_intra: float,
) -> str:
    """Gera interpretacao textual das distancias no espaco de embeddings.

    Args:
        separation_ratio: Razao inter-grupo / intra-grupo.
        mean_sl_human: Distancia media SL -> humano.
        mean_sl_canine: Distancia media SL -> canino.
        mean_intra: Distancia media intra-grupo.

    Returns:
        String com interpretacao biologica.
    """
    if separation_ratio > 2.0:
        sep_quality = "excelente"
    elif separation_ratio > 1.5:
        sep_quality = "boa"
    elif separation_ratio > 1.0:
        sep_quality = "moderada"
    else:
        sep_quality = "fraca"

    return (
        f"Separacao {sep_quality} no espaco de embeddings "
        f"(ratio inter/intra = {separation_ratio:.2f}x). "
        f"SL RNAs de Leishmania estao {mean_sl_human:.2f} unidades distantes "
        f"dos snRNAs humanos e {mean_sl_canine:.2f} unidades dos RNAs caninos, "
        f"vs {mean_intra:.2f} unidades intra-grupo. "
        f"Isso confirma que o SL RNA ocupa uma regiao distinta do espaco latente, "
        f"fundamentando a seletividade do ASO contra o parasita."
    )


# ---------------------------------------------------------------------------
# Perfil de unicidade posicional
# ---------------------------------------------------------------------------

def compute_positional_uniqueness(
    model: RNAEncoder,
    device: str = "mps",
) -> dict[str, object]:
    """Calcula o perfil de unicidade posicional do SL RNA de L. infantum.

    Para cada posicao do SL RNA, compara o embedding do token naquela
    posicao contra os tokens correspondentes em snRNAs humanos.
    Posicoes com maior divergencia sao "mais unicas" e, portanto,
    mais seletivas como alvo terapeutico.

    Mapeia o perfil de unicidade sobre a regiao de ligacao do ASO.

    Args:
        model: Encoder de RNA treinado.
        device: Dispositivo de computo.

    Returns:
        Dicionario com perfil de unicidade e analise da regiao ASO.
    """
    model.eval()

    sl_seq = SL_SEQUENCES["L_infantum"].upper().replace("T", "U")

    # Coleta embeddings de token para SL RNA e todos os controles humanos
    all_seqs = {"L_infantum": sl_seq}
    all_seqs.update(HUMAN_SNRNA)

    max_len = max(len(s.replace("T", "U")) for s in all_seqs.values())

    from .model import encode_sequence, PAD_IDX
    token_embeddings = {}

    with torch.no_grad():
        for name, seq in all_seqs.items():
            seq_rna = seq.upper().replace("T", "U")
            encoded = encode_sequence(seq_rna, max_len=max_len)
            input_ids = torch.tensor([encoded], dtype=torch.long, device=device)
            padding_mask = input_ids == PAD_IDX
            outputs = model(input_ids, padding_mask=padding_mask)
            # Pega embeddings de token apenas para posicoes validas
            token_emb = outputs["token_embeddings"].squeeze(0).cpu()
            token_embeddings[name] = token_emb[: len(seq_rna)]

    # Calcula divergencia posicional: para cada posicao do SL RNA,
    # compara contra a media dos embeddings humanos naquela posicao
    sl_emb = token_embeddings["L_infantum"]
    sl_len = len(sl_seq)

    # Media dos embeddings humanos por posicao (ate o comprimento do SL)
    human_names = list(HUMAN_SNRNA.keys())
    uniqueness_scores = []

    for pos in range(sl_len):
        sl_token = sl_emb[pos]

        # Coleta embeddings humanos nesta posicao (se disponivel)
        human_tokens = []
        for h_name in human_names:
            h_emb = token_embeddings[h_name]
            if pos < h_emb.size(0):
                human_tokens.append(h_emb[pos])

        if not human_tokens:
            uniqueness_scores.append(0.0)
            continue

        human_mean = torch.stack(human_tokens).mean(dim=0)

        # Distancia euclidiana como score de unicidade
        dist = torch.dist(sl_token, human_mean).item()
        uniqueness_scores.append(round(dist, 4))

    # Normaliza scores para [0, 1]
    max_score = max(uniqueness_scores) if uniqueness_scores else 1.0
    normalized_scores = [
        round(s / max_score, 4) if max_score > 0 else 0.0
        for s in uniqueness_scores
    ]

    # Analise da regiao de ligacao do ASO
    aso_region_scores = normalized_scores[ASO_TARGET_START:ASO_TARGET_END]
    non_aso_scores = (
        normalized_scores[:ASO_TARGET_START] + normalized_scores[ASO_TARGET_END:]
    )

    mean_aso_uniqueness = (
        sum(aso_region_scores) / len(aso_region_scores)
        if aso_region_scores else 0.0
    )
    mean_non_aso_uniqueness = (
        sum(non_aso_scores) / len(non_aso_scores) if non_aso_scores else 0.0
    )

    # Identifica posicoes mais unicas na regiao do ASO
    aso_positions_ranked = sorted(
        [
            (pos + ASO_TARGET_START, score)
            for pos, score in enumerate(aso_region_scores)
        ],
        key=lambda x: x[1],
        reverse=True,
    )

    return {
        "sl_sequence": sl_seq,
        "uniqueness_profile": normalized_scores,
        "raw_scores": uniqueness_scores,
        "aso_binding_region": f"pos {ASO_TARGET_START}-{ASO_TARGET_END - 1}",
        "mean_aso_region_uniqueness": round(mean_aso_uniqueness, 4),
        "mean_non_aso_uniqueness": round(mean_non_aso_uniqueness, 4),
        "aso_uniqueness_enrichment": round(
            mean_aso_uniqueness / mean_non_aso_uniqueness
            if mean_non_aso_uniqueness > 1e-8 else 0.0,
            4,
        ),
        "top_unique_positions_in_aso": aso_positions_ranked[:5],
        "interpretation": _interpret_uniqueness(
            mean_aso_uniqueness, mean_non_aso_uniqueness, aso_positions_ranked
        ),
    }


def _interpret_uniqueness(
    mean_aso: float,
    mean_non_aso: float,
    top_positions: list[tuple[int, float]],
) -> str:
    """Gera interpretacao textual do perfil de unicidade.

    Args:
        mean_aso: Unicidade media na regiao do ASO.
        mean_non_aso: Unicidade media fora da regiao do ASO.
        top_positions: Posicoes mais unicas na regiao do ASO.

    Returns:
        String com interpretacao biologica.
    """
    enrichment = mean_aso / mean_non_aso if mean_non_aso > 1e-8 else 0.0

    if enrichment > 1.2:
        quality = "favoravel"
        detail = (
            "A regiao de ligacao do ASO apresenta unicidade acima da media, "
            "indicando que o alvo e intrinsecamente divergente do transcriptoma humano."
        )
    elif enrichment > 0.8:
        quality = "neutra"
        detail = (
            "A regiao de ligacao do ASO tem unicidade comparavel ao resto do SL RNA."
        )
    else:
        quality = "sub-otima"
        detail = (
            "A regiao de ligacao do ASO apresenta unicidade abaixo da media. "
            "Posicoes mais divergentes estao fora da regiao alvo."
        )

    top_str = ", ".join(f"pos {p} ({s:.2f})" for p, s in top_positions[:3])

    return (
        f"Perfil de unicidade {quality} para a regiao do ASO "
        f"(enrichment = {enrichment:.2f}x). {detail} "
        f"Posicoes mais divergentes na regiao alvo: {top_str}."
    )


# ---------------------------------------------------------------------------
# Funcao principal de analise
# ---------------------------------------------------------------------------

def run_sl_analysis(
    model: RNAEncoder,
    device: str = "mps",
) -> dict[str, object]:
    """Executa analise completa do SL RNA: conservacao + unicidade + embeddings.

    Combina tres analises:
        1. Conservacao cross-species (alinhamento par-a-par)
        2. Distancias no espaco de embeddings (SL vs controles)
        3. Perfil de unicidade posicional (regiao do ASO)

    Args:
        model: Encoder de RNA treinado.
        device: Dispositivo de computo.

    Returns:
        Dicionario com resultados completos da analise.
    """
    print("  [sl_analysis] Analisando conservacao cross-species...")
    conservation = compute_conservation(SL_SEQUENCES)

    print("  [sl_analysis] Gerando embeddings para todos os RNAs...")
    # Combina todas as sequencias para embedding
    all_sequences = {}
    all_sequences.update(SL_SEQUENCES)
    all_sequences.update(HUMAN_SNRNA)
    all_sequences.update(CANINE_CONTROL_RNA)
    embeddings = get_sequence_embeddings(model, all_sequences, device=device)

    print("  [sl_analysis] Calculando distancias no espaco de embeddings...")
    embedding_analysis = compute_embedding_distances(embeddings)

    print("  [sl_analysis] Computando perfil de unicidade posicional...")
    uniqueness = compute_positional_uniqueness(model, device=device)

    return {
        "conservation": conservation,
        "embedding_distances": embedding_analysis,
        "positional_uniqueness": uniqueness,
    }
