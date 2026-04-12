"""Enumerador de configuracoes LNA-DNA-LNA gapmer validas para ASOs de 25 nt.

Enumera todas as combinacoes de flancos LNA 5' e 3' com gap de DNA central
que satisfazem as restricoes biologicas para ativacao de RNase H.

Restricoes:
    - LNA 5': minimo 2, maximo 5 posicoes
    - LNA 3': minimo 2, maximo 5 posicoes
    - DNA gap: minimo 8 nt (requisito RNase H, Crooke et al. 2017)
    - Total LNA: maximo 10 posicoes (limiar de toxicidade)
    - Total: lna_5p + dna_gap + lna_3p == comprimento do ASO

Referencia:
    McTigue PM et al. (2004) Biochemistry 43(18):5388-5405
    Crooke ST et al. (2017) Nat Rev Drug Discov
"""

from __future__ import annotations

from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Constantes de design
# ---------------------------------------------------------------------------

LNA_FLANK_MIN: int = 2
LNA_FLANK_MAX: int = 5
DNA_GAP_MIN: int = 8
TOTAL_LNA_MAX: int = 10


# ---------------------------------------------------------------------------
# Dataclass de configuracao
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GapWindowConfig:
    """Configuracao de gapmer LNA-DNA-LNA.

    Attributes:
        lna_5p: Numero de posicoes LNA no flanco 5'.
        dna_gap: Numero de posicoes de DNA no gap central.
        lna_3p: Numero de posicoes LNA no flanco 3'.
    """

    lna_5p: int
    dna_gap: int
    lna_3p: int

    @property
    def total_length(self) -> int:
        """Comprimento total do ASO (LNA 5' + DNA gap + LNA 3')."""
        return self.lna_5p + self.dna_gap + self.lna_3p

    @property
    def total_lna(self) -> int:
        """Numero total de posicoes LNA."""
        return self.lna_5p + self.lna_3p

    def notation(self) -> str:
        """Retorna notacao legivel: ex. '4-17-4'."""
        return f"{self.lna_5p}-{self.dna_gap}-{self.lna_3p}"


# ---------------------------------------------------------------------------
# Enumeracao
# ---------------------------------------------------------------------------


def enumerate_gap_configs(
    aso_length: int = 25,
    lna_min: int = LNA_FLANK_MIN,
    lna_max: int = LNA_FLANK_MAX,
    gap_min: int = DNA_GAP_MIN,
    total_lna_max: int = TOTAL_LNA_MAX,
) -> list[GapWindowConfig]:
    """Enumera todas as configuracoes LNA-DNA-LNA gapmer validas.

    Para cada combinacao de (lna_5p, lna_3p) em [lna_min, lna_max] x [lna_min, lna_max],
    calcula o gap de DNA resultante e verifica as restricoes.

    Args:
        aso_length: Comprimento total do ASO. Default: 25 nt.
        lna_min: Minimo de posicoes LNA por flanco. Default: 2.
        lna_max: Maximo de posicoes LNA por flanco. Default: 5.
        gap_min: Minimo de DNA gap para RNase H. Default: 8.
        total_lna_max: Maximo total de LNA. Default: 10.

    Returns:
        Lista de GapWindowConfig validas, ordenada por gap decrescente.
    """
    configs: list[GapWindowConfig] = []

    for lna_5p in range(lna_min, lna_max + 1):
        for lna_3p in range(lna_min, lna_max + 1):
            dna_gap = aso_length - lna_5p - lna_3p

            # Restricao: gap minimo para ativacao de RNase H
            if dna_gap < gap_min:
                continue

            # Restricao: total de LNA nao pode exceder limiar de toxicidade
            if lna_5p + lna_3p > total_lna_max:
                continue

            # Restricao: gap deve ser positivo (redundante mas explicito)
            if dna_gap <= 0:
                continue

            configs.append(GapWindowConfig(
                lna_5p=lna_5p,
                dna_gap=dna_gap,
                lna_3p=lna_3p,
            ))

    # Ordenar por gap decrescente (maior gap = maior eficiencia RNase H potencial)
    configs.sort(key=lambda c: (-c.dna_gap, c.lna_5p, c.lna_3p))

    return configs
