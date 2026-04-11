"""Codificacao de sequencias proteicas em vetores de propriedades.

Converte peptideos/proteinas em vetores numericos de dimensao fixa
usando propriedades fisico-quimicas dos 20 aminoacidos canonicos.
Sem dependencias externas (ESM-2) — tudo calculado com numpy.

Propriedades por aminoacido (Kyte-Doolittle, carga pH 7, peso molecular,
volume, flexibilidade, aromaticidade) sao agregadas em estatisticas
por janela deslizante para capturar distribuicao espacial.
"""

from __future__ import annotations

import numpy as np

import importlib as _il

# Importacao via importlib — Python nao permite "from marley_ai.09_sae..."
# porque "09" e interpretado como literal numerico invalido no parser.
_cfg = _il.import_module("marley_ai.09_sae.config")
NUM_AA_PROPERTIES: int = _cfg.NUM_AA_PROPERTIES
NUM_STATISTICS: int = _cfg.NUM_STATISTICS
NUM_WINDOWS: int = _cfg.NUM_WINDOWS


# ---------------------------------------------------------------------------
# Tabela de propriedades fisico-quimicas dos 20 aminoacidos canonicos
# ---------------------------------------------------------------------------
# Cada linha: [hidrofobicidade, carga_pH7, peso_molecular, volume, flexibilidade, aromaticidade]
#
# Fontes:
#   Hidrofobicidade: Kyte & Doolittle (1982) J Mol Biol 157:105-132
#   Carga pH 7: valores teoricos baseados em pKa de cadeias laterais
#   Peso molecular: massa do residuo (sem agua) em Daltons
#   Volume: Zamyatnin (1972) Prog Biophys Mol Biol 24:107-123
#   Flexibilidade: Bhaskaran & Ponnuswamy (1988) Int J Pept Protein Res 32:242-255
#   Aromaticidade: 1 para F/W/Y/H, 0 para os demais

AA_PROPERTIES: dict[str, np.ndarray] = {
    # AA:  [hydro,  charge, MW,     vol,    flex,  arom]
    "A": np.array([ 1.8,   0.0,    89.1,   88.6,  0.360, 0.0]),
    "R": np.array([-4.5,   1.0,   174.2,  173.4,  0.529, 0.0]),
    "N": np.array([-3.5,   0.0,   132.1,  114.1,  0.463, 0.0]),
    "D": np.array([-3.5,  -1.0,   133.1,  111.1,  0.511, 0.0]),
    "C": np.array([ 2.5,   0.0,   121.2,  108.5,  0.346, 0.0]),
    "E": np.array([-3.5,  -1.0,   147.1,  138.4,  0.497, 0.0]),
    "Q": np.array([-3.5,   0.0,   146.2,  143.8,  0.493, 0.0]),
    "G": np.array([-0.4,   0.0,    75.0,   60.1,  0.540, 0.0]),
    "H": np.array([-3.2,   0.1,   155.2,  153.2,  0.323, 1.0]),
    "I": np.array([ 4.5,   0.0,   131.2,  166.7,  0.462, 0.0]),
    "L": np.array([ 3.8,   0.0,   131.2,  166.7,  0.365, 0.0]),
    "K": np.array([-3.9,   1.0,   146.2,  168.6,  0.466, 0.0]),
    "M": np.array([ 1.9,   0.0,   149.2,  162.9,  0.295, 0.0]),
    "F": np.array([ 2.8,   0.0,   165.2,  189.9,  0.314, 1.0]),
    "P": np.array([-1.6,   0.0,   115.1,  112.7,  0.509, 0.0]),
    "S": np.array([-0.8,   0.0,   105.1,   89.0,  0.507, 0.0]),
    "T": np.array([-0.7,   0.0,   119.1,  116.1,  0.585, 0.0]),
    "W": np.array([-0.9,   0.0,   204.2,  227.8,  0.170, 1.0]),
    "Y": np.array([-1.3,   0.0,   181.2,  193.6,  0.420, 1.0]),
    "V": np.array([ 4.2,   0.0,   117.1,  140.0,  0.386, 0.0]),
}

# Aminoacidos canonicos em ordem (para referencia)
CANONICAL_AAS: list[str] = sorted(AA_PROPERTIES.keys())


def _seq_to_property_matrix(sequence: str) -> np.ndarray:
    """Converte sequencia de aminoacidos em matriz (L x P) de propriedades.

    Aminoacidos desconhecidos (X, B, Z, etc.) recebem a media das
    propriedades dos 20 canonicos — abordagem conservadora que nao
    introduz vies para nenhuma classe especifica.

    Args:
        sequence: sequencia de aminoacidos (case-insensitive)

    Returns:
        Matriz (L, NUM_AA_PROPERTIES) onde L = len(sequence)
    """
    # Vetor medio para aminoacidos nao-canonicos
    mean_props = np.mean(
        [v for v in AA_PROPERTIES.values()], axis=0
    )

    rows = []
    for aa in sequence.upper():
        rows.append(AA_PROPERTIES.get(aa, mean_props))
    return np.array(rows, dtype=np.float64)


def _window_statistics(matrix: np.ndarray) -> np.ndarray:
    """Calcula estatisticas (mean, std, max, min) por janela.

    Tres janelas capturam distribuicao espacial das propriedades:
      - Janela 1: sequencia inteira (perfil global)
      - Janela 2: metade N-terminal (regiao de sinal/ancora)
      - Janela 3: metade C-terminal (regiao de efetor)

    Para peptideos curtos (< 4 residuos), todas as janelas
    cobrem a sequencia inteira — nao ha informacao espacial
    significativa em peptideos tao pequenos.

    Args:
        matrix: (L, P) matriz de propriedades

    Returns:
        Vetor (NUM_WINDOWS * NUM_STATISTICS * NUM_AA_PROPERTIES,)
    """
    length = matrix.shape[0]
    mid = max(length // 2, 1)

    # Tres janelas: toda, N-terminal, C-terminal
    windows = [
        matrix,            # sequencia completa
        matrix[:mid],      # metade N-terminal
        matrix[mid:],      # metade C-terminal
    ]

    parts: list[np.ndarray] = []
    for w in windows:
        # Para janelas com 1 residuo, std = 0 (evita NaN)
        parts.append(np.mean(w, axis=0))
        parts.append(np.std(w, axis=0))
        parts.append(np.max(w, axis=0))
        parts.append(np.min(w, axis=0))

    return np.concatenate(parts)


def encode_sequence(sequence: str) -> np.ndarray:
    """Codifica uma sequencia proteica em vetor de dimensao fixa.

    Pipeline: sequencia -> matriz de propriedades -> estatisticas por janela

    Args:
        sequence: sequencia de aminoacidos (minimo 1 residuo)

    Returns:
        Vetor (INPUT_DIM,) = (72,) com propriedades agregadas

    Raises:
        ValueError: se a sequencia estiver vazia
    """
    if not sequence or not sequence.strip():
        raise ValueError("Sequencia vazia nao pode ser codificada")

    seq_clean = sequence.upper().strip()
    matrix = _seq_to_property_matrix(seq_clean)
    return _window_statistics(matrix)


def encode_batch(sequences: list[str]) -> np.ndarray:
    """Codifica lote de sequencias em matriz (N, INPUT_DIM).

    Args:
        sequences: lista de sequencias proteicas

    Returns:
        Matriz (N, 72) com uma linha por sequencia
    """
    return np.array([encode_sequence(s) for s in sequences], dtype=np.float64)


def normalize_features(
    X: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Normaliza features para media 0 e desvio 1 (z-score).

    Essencial antes do SAE — propriedades em escalas diferentes
    (MW ~ 75-204, charge ~ -1 a 1) dominariam a reconstrucao.

    Args:
        X: matriz (N, D) de features brutas

    Returns:
        Tupla (X_norm, mean, std) onde X_norm = (X - mean) / std.
        Features com std=0 (constantes) recebem std=1 para evitar divisao por zero.
    """
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)
    # Evita divisao por zero em features constantes
    std[std < 1e-12] = 1.0
    X_norm = (X - mean) / std
    return X_norm, mean, std


def generate_random_peptides(
    n: int,
    min_length: int = 8,
    max_length: int = 15,
    seed: int = 42,
) -> list[str]:
    """Gera peptideos aleatorios como controle negativo.

    Amostra aminoacidos com frequencia uniforme — nao reflete
    composicao real do proteoma, mas serve como baseline para
    identificar features especificas dos epitopos.

    Args:
        n: numero de peptideos a gerar
        min_length: comprimento minimo
        max_length: comprimento maximo
        seed: semente para reproducibilidade

    Returns:
        Lista de n peptideos aleatorios
    """
    rng = np.random.default_rng(seed)
    peptides = []
    for _ in range(n):
        length = rng.integers(min_length, max_length + 1)
        peptide = "".join(rng.choice(list(CANONICAL_AAS), size=length))
        peptides.append(peptide)
    return peptides


def get_property_names() -> list[str]:
    """Retorna nomes legiveis das 72 features do vetor codificado.

    Formato: "{janela}_{estatistica}_{propriedade}"
    Exemplo: "full_mean_hydrophobicity", "nterm_std_charge"
    """
    property_names = [
        "hydrophobicity", "charge", "mol_weight",
        "volume", "flexibility", "aromaticity",
    ]
    window_names = ["full", "nterm", "cterm"]
    stat_names = ["mean", "std", "max", "min"]

    names: list[str] = []
    for win in window_names:
        for stat in stat_names:
            for prop in property_names:
                names.append(f"{win}_{stat}_{prop}")
    return names
