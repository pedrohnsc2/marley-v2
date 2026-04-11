"""Espaco de design para ASOs tipo gapmer LNA-DNA-LNA.

Define a parametrizacao do espaco de busca e as funcoes para converter
vetores de parametros em sequencias ASO concretas. Garante que todas
as restricoes de design gapmer sao satisfeitas (LNA apenas nas flanking
regions, gap central de DNA puro para ativacao de RNase H).
"""

from __future__ import annotations

import numpy as np

from aso_math.config import SL_SEQUENCE
from aso_math.thermo import reverse_complement


# ---------------------------------------------------------------------------
# Limites do espaco de design
# ---------------------------------------------------------------------------

# Comprimento do ASO: 20 a 30 nt
ASO_LENGTH_MIN: int = 20
ASO_LENGTH_MAX: int = 30

# Gap central de DNA: 10 a 20 posicoes (requisito para RNase H)
GAP_SIZE_MIN: int = 10
GAP_SIZE_MAX: int = 20

# Flanco LNA minimo: 2 posicoes de cada lado (para estabilidade)
LNA_FLANK_MIN: int = 2

# Posicao de inicio no SL RNA (0-indexed)
START_POS_MIN: int = 0
START_POS_MAX: int = 14  # para 25-mer em SL de 39 nt


# ---------------------------------------------------------------------------
# Representacao do design
# ---------------------------------------------------------------------------


class ASODesign:
    """Representacao completa de um design ASO gapmer.

    Atributos:
        length: Comprimento total do ASO (20-30 nt).
        gap_size: Tamanho do gap de DNA central.
        start_pos: Posicao de inicio no SL RNA alvo (0-indexed).
        lna_5prime: Numero de posicoes LNA no flanco 5'.
        lna_3prime: Numero de posicoes LNA no flanco 3'.
        sequence: Sequencia DNA do ASO (complemento reverso do alvo).
        lna_mask: Vetor booleano indicando posicoes com LNA.
    """

    __slots__ = (
        "length", "gap_size", "start_pos",
        "lna_5prime", "lna_3prime",
        "sequence", "lna_mask",
    )

    def __init__(
        self,
        length: int,
        gap_size: int,
        start_pos: int,
        lna_5prime: int,
        lna_3prime: int,
        sl_sequence: str = SL_SEQUENCE,
    ) -> None:
        self.length = length
        self.gap_size = gap_size
        self.start_pos = start_pos
        self.lna_5prime = lna_5prime
        self.lna_3prime = lna_3prime

        # Gerar a sequencia ASO a partir do alvo no SL RNA
        target_region = sl_sequence[start_pos:start_pos + length].upper()
        self.sequence = reverse_complement(target_region)

        # Construir mascara LNA: True nas flanking regions, False no gap
        self.lna_mask = np.zeros(length, dtype=bool)
        self.lna_mask[:lna_5prime] = True
        self.lna_mask[length - lna_3prime:] = True

    @property
    def total_lna(self) -> int:
        """Numero total de posicoes com modificacao LNA."""
        return int(np.sum(self.lna_mask))

    @property
    def architecture_str(self) -> str:
        """Representacao textual da arquitetura gapmer (ex: '4-17-4')."""
        return f"{self.lna_5prime}-{self.gap_size}-{self.lna_3prime}"

    def to_dict(self) -> dict:
        """Serializa o design para JSON."""
        return {
            "sequence": self.sequence,
            "length": self.length,
            "gap_size": self.gap_size,
            "start_pos": self.start_pos,
            "lna_5prime": self.lna_5prime,
            "lna_3prime": self.lna_3prime,
            "total_lna": self.total_lna,
            "architecture": self.architecture_str,
            "lna_mask": self.lna_mask.tolist(),
        }


# ---------------------------------------------------------------------------
# Parametrizacao e amostragem
# ---------------------------------------------------------------------------


def params_to_design(params: np.ndarray, sl_sequence: str = SL_SEQUENCE) -> ASODesign | None:
    """Converte um vetor de 4 parametros continuos em um ASODesign valido.

    Parametros (todos em [0, 1], escalonados internamente):
        params[0]: comprimento do ASO (escalonado para 20-30)
        params[1]: fracao do gap em relacao ao comprimento disponivel
        params[2]: posicao de inicio no SL RNA (escalonado para range valido)
        params[3]: assimetria do flanco LNA (0 = simetrico, 1 = max 5')

    Retorna None se as restricoes nao puderem ser satisfeitas.
    """
    if len(params) != 4:
        return None

    # Descretizar comprimento
    length = int(round(params[0] * (ASO_LENGTH_MAX - ASO_LENGTH_MIN) + ASO_LENGTH_MIN))
    length = max(ASO_LENGTH_MIN, min(ASO_LENGTH_MAX, length))

    # Posicao de inicio: depende do comprimento e do SL RNA
    sl_len = len(sl_sequence)
    max_start = sl_len - length
    if max_start < 0:
        return None
    start_pos = int(round(params[2] * max_start))
    start_pos = max(0, min(max_start, start_pos))

    # Gap size: entre GAP_SIZE_MIN e min(GAP_SIZE_MAX, length - 2*LNA_FLANK_MIN)
    max_gap = length - 2 * LNA_FLANK_MIN
    if max_gap < GAP_SIZE_MIN:
        return None
    effective_gap_max = min(GAP_SIZE_MAX, max_gap)
    gap_size = int(round(params[1] * (effective_gap_max - GAP_SIZE_MIN) + GAP_SIZE_MIN))
    gap_size = max(GAP_SIZE_MIN, min(effective_gap_max, gap_size))

    # Distribuir LNA nas flanking regions
    total_flank = length - gap_size
    if total_flank < 2 * LNA_FLANK_MIN:
        return None

    # Assimetria: como distribuir os flancos entre 5' e 3'
    # params[3] controla a fracao alocada ao 5'
    lna_5prime = int(round(params[3] * (total_flank - 2 * LNA_FLANK_MIN)) + LNA_FLANK_MIN)
    lna_5prime = max(LNA_FLANK_MIN, min(total_flank - LNA_FLANK_MIN, lna_5prime))
    lna_3prime = total_flank - lna_5prime

    # Validacao final
    if lna_5prime < LNA_FLANK_MIN or lna_3prime < LNA_FLANK_MIN:
        return None
    if gap_size < GAP_SIZE_MIN:
        return None
    if lna_5prime + gap_size + lna_3prime != length:
        return None

    return ASODesign(
        length=length,
        gap_size=gap_size,
        start_pos=start_pos,
        lna_5prime=lna_5prime,
        lna_3prime=lna_3prime,
        sl_sequence=sl_sequence,
    )


def sample_random_designs(
    n_samples: int,
    rng: np.random.Generator,
    sl_sequence: str = SL_SEQUENCE,
) -> list[tuple[np.ndarray, ASODesign]]:
    """Amostra designs aleatorios validos no espaco de busca.

    Gera parametros uniformes em [0, 1]^4 e filtra designs invalidos.
    Continua amostrando ate obter n_samples designs validos.

    Args:
        n_samples: Numero de designs validos desejados.
        rng: Gerador de numeros aleatorios numpy.
        sl_sequence: Sequencia do SL RNA alvo.

    Returns:
        Lista de tuplas (vetor_parametros, ASODesign).
    """
    designs: list[tuple[np.ndarray, ASODesign]] = []
    max_attempts = n_samples * 20  # limite de seguranca contra loop infinito
    attempts = 0

    while len(designs) < n_samples and attempts < max_attempts:
        params = rng.uniform(0.0, 1.0, size=4)
        design = params_to_design(params, sl_sequence)
        if design is not None:
            designs.append((params, design))
        attempts += 1

    return designs


def mrl_aso_001_params(sl_sequence: str = SL_SEQUENCE) -> tuple[np.ndarray, ASODesign]:
    """Retorna os parametros normalizados correspondentes ao MRL-ASO-001.

    MRL-ASO-001: 25 nt, gap 17 (4-17-4), start pos 5 no SL RNA.
    """
    length = 25
    gap_size = 17
    start_pos = 5
    lna_5prime = 4
    lna_3prime = 4

    # Converter para espaco [0, 1] (inverso de params_to_design)
    p0 = (length - ASO_LENGTH_MIN) / (ASO_LENGTH_MAX - ASO_LENGTH_MIN)

    sl_len = len(sl_sequence)
    max_start = sl_len - length
    p2 = start_pos / max_start if max_start > 0 else 0.0

    max_gap = length - 2 * LNA_FLANK_MIN
    effective_gap_max = min(GAP_SIZE_MAX, max_gap)
    p1 = (gap_size - GAP_SIZE_MIN) / (effective_gap_max - GAP_SIZE_MIN)

    total_flank = length - gap_size
    p3 = (lna_5prime - LNA_FLANK_MIN) / (total_flank - 2 * LNA_FLANK_MIN) if total_flank > 2 * LNA_FLANK_MIN else 0.5

    params = np.array([p0, p1, p2, p3])

    design = ASODesign(
        length=length,
        gap_size=gap_size,
        start_pos=start_pos,
        lna_5prime=lna_5prime,
        lna_3prime=lna_3prime,
        sl_sequence=sl_sequence,
    )

    return params, design
