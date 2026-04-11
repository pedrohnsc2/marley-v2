"""Modelo de liberacao pH-responsiva de LNP com lipidio ionizavel.

O mecanismo central das LNPs modernas e a liberacao pH-responsiva:
- Em pH 7.4 (sangue): lipidio ionizavel NEUTRO -> LNP estavel, stealth
- Em pH <6.5 (endossomo): lipidio PROTONADO -> desestabiliza membrana -> libera carga

Para MRL-ASO-001, isto e particularmente vantajoso porque:
1. O ASO precisa chegar ao fagolisossomo (pH 4.5)
2. A LNP e endocitada pelo macrofago infectado
3. No endossomo tardio/fagolisossomo, pH cai abaixo do pKa do MC3
4. MC3 protona -> interage com fosfolipidios endossomais -> poro hexagonal
5. ASO liberado DIRETAMENTE no compartimento onde o SL RNA esta

Modelo sigmoidal de liberacao:
    f_released(pH) = 1 / (1 + exp(k * (pH - pKa)))

onde k controla a inclinacao da curva (Hill-like) e pKa e o ponto
de inflexao (50% liberacao). Para DLin-MC3-DMA: pKa = 6.44.

Referencias:
- Jayaraman M et al. (2012) Angew Chem 51(34):8529-8533 — pKa de MC3
- Kulkarni JA et al. (2018) Nanoscale 10(10):4567-4573 — mecanismo de release
- Sahay G et al. (2013) Nat Biotechnol 31(7):653-658 — escape endossomal
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Constantes do modelo de liberacao
# ---------------------------------------------------------------------------

# pKa do lipidio ionizavel DLin-MC3-DMA
# Ref: Jayaraman M et al. (2012) — valor determinado por TNS assay
PKA_MC3: Final[float] = 6.44

# Inclinacao da curva sigmoidal (Hill coefficient analog)
# Valor mais alto = transicao mais abrupta
# Para LNPs com MC3: k ~3-5 (transicao cooperativa)
# Ref: Kulkarni JA et al. (2018) — ajuste a dados experimentais de TNS
HILL_SLOPE: Final[float] = 4.0

# Compartimentos intracelulares e seus pH
# Nota: o fagolisossomo e o destino final para macrofagos infectados
COMPARTMENT_PH: Final[dict[str, float]] = {
    "blood_plasma": 7.4,
    "early_endosome": 6.5,
    "late_endosome": 5.0,
    "phagolysosome": 4.5,
}

# Tempo de transito estimado entre compartimentos (minutos)
# Ref: Sahay G et al. (2013) — cinetica de trafego intracelular
TRANSIT_TIMES_MIN: Final[dict[str, float]] = {
    "blood_plasma": 0.0,        # tempo zero (ponto de partida)
    "early_endosome": 5.0,      # 5 min apos endocitose
    "late_endosome": 15.0,      # 15 min apos endocitose
    "phagolysosome": 30.0,      # 30 min (fusao com fagossomo)
}


# ---------------------------------------------------------------------------
# Dataclasses de resultado
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReleaseAtPH:
    """Fracao de ASO liberada em um dado pH/compartimento.

    O modelo sigmoidal captura a cooperatividade da transicao:
    - pH > pKa + 1: <10% liberacao (MC3 neutro -> LNP estavel)
    - pH ~ pKa: 50% liberacao (ponto de inflexao)
    - pH < pKa - 1: >90% liberacao (MC3 protonado -> desestabilizacao)
    """

    compartment: str
    ph: float
    fraction_protonated: float    # fracao de MC3 protonada
    fraction_released: float      # fracao do ASO liberada
    transit_time_min: float       # tempo para chegar a este compartimento
    release_status: str           # "minimal" / "partial" / "complete"


@dataclass(frozen=True)
class ReleaseProfile:
    """Perfil completo de liberacao ao longo do pathway intracelular.

    Resume a cinetica de liberacao desde o sangue ate o fagolisossomo,
    destacando onde e quando a carga e disponibilizada.
    """

    pka_ionizable: float
    hill_slope: float
    compartment_releases: dict[str, ReleaseAtPH]
    release_at_target: float      # fracao liberada no fagolisossomo
    target_compartment: str
    advantage_for_macrophage: str  # explicacao da vantagem


# ---------------------------------------------------------------------------
# Funcoes de modelo
# ---------------------------------------------------------------------------


def compute_fraction_protonated(ph: float, pka: float = PKA_MC3) -> float:
    """Calcula fracao de lipidio ionizavel protonado via Henderson-Hasselbalch.

    f = 1 / (1 + 10^(pH - pKa))

    A pH 7.4: f ~ 0.01 (1% protonado -> neutro -> stealth)
    A pH 4.5: f ~ 0.999 (99.9% protonado -> cationico -> libera)

    Args:
        ph: pH do meio.
        pka: pKa do lipidio ionizavel (default: 6.44 para MC3).

    Returns:
        Fracao protonada (0 a 1).
    """
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def compute_fraction_released(
    ph: float,
    pka: float = PKA_MC3,
    k: float = HILL_SLOPE,
) -> float:
    """Calcula fracao de ASO liberada pela LNP em funcao do pH.

    Modelo sigmoidal:
        f_released = 1 / (1 + exp(k * (pH - pKa)))

    A curva e mais abrupta que Henderson-Hasselbalch puro porque
    a desestabilizacao da LNP envolve cooperatividade:
    1. MC3 protona
    2. MC3+ interage com aniônios da membrana endossomal
    3. Fase hexagonal invertida (H_II) forma-se
    4. LNP desintegra e libera conteudo

    Args:
        ph: pH do meio.
        pka: pKa do lipidio ionizavel.
        k: Inclinacao da curva (Hill slope analog).

    Returns:
        Fracao liberada (0 a 1).
    """
    exponent = k * (ph - pka)
    # Protecao contra overflow numerico
    if exponent > 500:
        return 0.0
    if exponent < -500:
        return 1.0
    return 1.0 / (1.0 + math.exp(exponent))


def classify_release(fraction: float) -> str:
    """Classifica o nivel de liberacao.

    Args:
        fraction: Fracao liberada (0 a 1).

    Returns:
        "minimal" (<10%), "partial" (10-80%), ou "complete" (>80%).
    """
    if fraction < 0.10:
        return "minimal"
    if fraction < 0.80:
        return "partial"
    return "complete"


def compute_release_profile() -> ReleaseProfile:
    """Calcula perfil completo de liberacao ao longo da jornada intracelular.

    Para cada compartimento (sangue -> endossomo precoce -> endossomo tardio
    -> fagolisossomo), calcula a fracao de MC3 protonada e a fracao de ASO
    liberada.

    O insight chave: a liberacao MAXIMA ocorre no fagolisossomo (pH 4.5),
    que e EXATAMENTE onde o ASO precisa atuar contra o SL RNA de
    L. infantum. A LNP age como um cavalo de Troia biologico.

    Returns:
        ReleaseProfile com dados completos.
    """
    releases: dict[str, ReleaseAtPH] = {}

    for compartment, ph in COMPARTMENT_PH.items():
        f_prot = compute_fraction_protonated(ph)
        f_rel = compute_fraction_released(ph)
        status = classify_release(f_rel)
        transit = TRANSIT_TIMES_MIN[compartment]

        releases[compartment] = ReleaseAtPH(
            compartment=compartment,
            ph=ph,
            fraction_protonated=round(f_prot, 6),
            fraction_released=round(f_rel, 6),
            transit_time_min=transit,
            release_status=status,
        )

    # Liberacao no alvo (fagolisossomo)
    target = releases["phagolysosome"]

    advantage = (
        "LNP pH-responsive release is ADVANTAGEOUS for macrophage delivery: "
        f"at pH {target.ph} (phagolysosome), {target.fraction_released*100:.1f}% "
        f"of the ASO cargo is released. Since the target RNA (SL RNA) is located "
        f"in the SAME compartment where L. infantum amastigotes reside, "
        f"endosomal escape is NOT required — the LNP releases MRL-ASO-001 "
        f"directly where it needs to act."
    )

    return ReleaseProfile(
        pka_ionizable=PKA_MC3,
        hill_slope=HILL_SLOPE,
        compartment_releases=releases,
        release_at_target=target.fraction_released,
        target_compartment="phagolysosome",
        advantage_for_macrophage=advantage,
    )
