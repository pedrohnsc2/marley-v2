"""Comparacao de MRL-ASO-001 com ASOs aprovados pelo FDA.

Dados de referencia para contextualizar a estabilidade de MRL-ASO-001
no ambiente acido do fagolisossomo frente a ASOs clinicamente validados.

Os tres ASOs comparadores foram escolhidos por razoes especificas:
- Nusinersen (Spinraza): referencia de eficacia de gapmer 2'-MOE,
  administracao intratecal em pH 7.3 (baseline favoravel)
- Mipomersen (Kynamro): melhor comparador — subcutaneo, enfrenta pH 4.5
  em lisossomos hepaticos, mesma quimica de base
- Inotersen (Tegsedi): subcutaneo, mesmo desafio de delivery lisossomal

Nota: MRL-ASO-001 enfrenta desafio ADICIONAL vs Kynamro/Tegsedi porque
o alvo e INTRAparasitario (dentro do fagolisossomo), nao citoplasmatico
do hepatocito. Porem, o SL RNA do parasita esta no MESMO compartimento
acido, eliminando a necessidade de escape endossomal.

Referencias:
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
- Bennett CF (2019) Annu Rev Pharmacol Toxicol 59:447-464
- Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True)
class FDAApprovedASO:
    """Dados de um ASO aprovado pelo FDA para comparacao.

    Apenas propriedades publicadas na literatura aberta.
    Nenhum dado proprietario ou confidencial.
    """

    generic_name: str
    brand_name: str
    fda_approval_year: int
    indication: str
    target_gene: str
    length_nt: int
    chemistry: str               # tipo de modificacao quimica
    backbone: str                # PO, PS, ou misto
    sugar_modification: str      # 2'-MOE, LNA, etc.
    gapmer_design: str           # ex: "5-10-5" para flancos de 5 nt
    route_of_administration: str
    target_tissue_ph: float      # pH do compartimento alvo
    reported_half_life_tissue_hours: float  # meia-vida em tecido
    mechanism: str               # RNase H, steric block, etc.
    key_reference: str           # referencia primaria


# ---------------------------------------------------------------------------
# Dados dos ASOs aprovados
# Fonte: labels FDA, revisoes publicadas, dados publicados
# ---------------------------------------------------------------------------

NUSINERSEN: Final[FDAApprovedASO] = FDAApprovedASO(
    generic_name="nusinersen",
    brand_name="Spinraza",
    fda_approval_year=2016,
    indication="Spinal muscular atrophy (SMA)",
    target_gene="SMN2 pre-mRNA (exon 7 inclusion)",
    length_nt=18,
    chemistry="2'-O-methoxyethyl (2'-MOE) phosphorothioate",
    backbone="full PS",
    sugar_modification="2'-MOE (uniform)",
    gapmer_design="N/A (uniform 2'-MOE, splice-switching)",
    route_of_administration="intrathecal",
    target_tissue_ph=7.3,
    # Meia-vida no CNS: 135-177 dias (Geary RS et al. 2015)
    reported_half_life_tissue_hours=135.0 * 24.0,
    mechanism="splice-switching (steric block)",
    key_reference="Finkel RS et al. (2017) NEJM 377(18):1723-1732",
)

MIPOMERSEN: Final[FDAApprovedASO] = FDAApprovedASO(
    generic_name="mipomersen",
    brand_name="Kynamro",
    fda_approval_year=2013,
    indication="Homozygous familial hypercholesterolemia (HoFH)",
    target_gene="APOB-100 mRNA",
    length_nt=20,
    chemistry="2'-O-methoxyethyl (2'-MOE) gapmer phosphorothioate",
    backbone="full PS",
    sugar_modification="2'-MOE flanks, DNA gap",
    gapmer_design="5-10-5",
    route_of_administration="subcutaneous",
    target_tissue_ph=4.5,
    # Meia-vida em tecido hepatico: ~30 dias (Geary RS et al. 2015)
    reported_half_life_tissue_hours=30.0 * 24.0,
    mechanism="RNase H1-mediated mRNA degradation",
    key_reference="Raal FJ et al. (2010) Lancet 375(9719):998-1006",
)

INOTERSEN: Final[FDAApprovedASO] = FDAApprovedASO(
    generic_name="inotersen",
    brand_name="Tegsedi",
    fda_approval_year=2018,
    indication="Hereditary transthyretin amyloidosis (hATTR)",
    target_gene="TTR mRNA",
    length_nt=20,
    chemistry="2'-O-methoxyethyl (2'-MOE) gapmer phosphorothioate",
    backbone="full PS",
    sugar_modification="2'-MOE flanks, DNA gap",
    gapmer_design="5-10-5",
    route_of_administration="subcutaneous",
    target_tissue_ph=4.5,
    # Meia-vida em tecido hepatico: ~23-29 dias (Benson MD et al. 2018)
    reported_half_life_tissue_hours=26.0 * 24.0,
    mechanism="RNase H1-mediated mRNA degradation",
    key_reference="Benson MD et al. (2018) NEJM 379(1):22-31",
)


# ---------------------------------------------------------------------------
# Dados de MRL-ASO-001 para comparacao
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MRLASOComparison:
    """Dados do MRL-ASO-001 formatados para comparacao com FDA ASOs.

    Nota: a meia-vida aqui e a estimativa do modelo cinetico
    (module_a_stability/models.py), nao um dado experimental.
    """

    name: str
    length_nt: int
    chemistry: str
    backbone: str
    sugar_modification: str
    gapmer_design: str
    route_of_administration: str
    target_tissue_ph: float
    estimated_half_life_hours: float
    mechanism: str
    unique_advantage: str


def build_mrl_comparison(estimated_half_life_hours: float) -> MRLASOComparison:
    """Constroi objeto de comparacao para MRL-ASO-001.

    A meia-vida estimada vem do modelo cinetico em models.py.
    O campo unique_advantage destaca o diferencial biologico.

    Args:
        estimated_half_life_hours: Meia-vida estimada pelo modelo cinetico.

    Returns:
        Dados do MRL-ASO-001 formatados para comparacao.
    """
    return MRLASOComparison(
        name="MRL-ASO-001",
        length_nt=25,
        chemistry="LNA-DNA-LNA gapmer phosphorothioate",
        backbone="full PS",
        sugar_modification="LNA flanks (5+5), DNA gap (15)",
        gapmer_design="5-15-5",
        route_of_administration="to be determined (parenteral expected)",
        target_tissue_ph=4.5,
        estimated_half_life_hours=estimated_half_life_hours,
        mechanism="RNase H1-mediated SL RNA degradation",
        unique_advantage=(
            "Target (SL RNA) and ASO colocalize in phagolysosome — "
            "no endosomal escape required. LNA flanks provide superior "
            "nuclease resistance vs 2'-MOE at pH 4.5."
        ),
    )


def build_comparison_table(
    mrl_half_life_hours: float,
) -> dict[str, dict[str, str | float | int]]:
    """Constroi tabela comparativa completa entre MRL-ASO-001 e ASOs FDA.

    Retorna dicionario com dados normalizados para comparacao direta.
    Cada entrada tem os mesmos campos para facilitar tabulacao.

    A comparacao mais relevante e com mipomersen, que enfrenta o mesmo
    desafio de pH 4.5 em lisossomos e usa mecanismo RNase H identico.

    Args:
        mrl_half_life_hours: Meia-vida estimada para MRL-ASO-001.

    Returns:
        Dicionario com dados comparativos de todos os ASOs.
    """
    mrl = build_mrl_comparison(mrl_half_life_hours)

    table: dict[str, dict[str, str | float | int]] = {}

    # Nusinersen
    table["nusinersen"] = {
        "brand_name": NUSINERSEN.brand_name,
        "approval_year": NUSINERSEN.fda_approval_year,
        "indication": NUSINERSEN.indication,
        "length_nt": NUSINERSEN.length_nt,
        "chemistry": NUSINERSEN.chemistry,
        "backbone": NUSINERSEN.backbone,
        "gapmer_design": NUSINERSEN.gapmer_design,
        "route": NUSINERSEN.route_of_administration,
        "target_ph": NUSINERSEN.target_tissue_ph,
        "half_life_hours": NUSINERSEN.reported_half_life_tissue_hours,
        "mechanism": NUSINERSEN.mechanism,
        "reference": NUSINERSEN.key_reference,
    }

    # Mipomersen
    table["mipomersen"] = {
        "brand_name": MIPOMERSEN.brand_name,
        "approval_year": MIPOMERSEN.fda_approval_year,
        "indication": MIPOMERSEN.indication,
        "length_nt": MIPOMERSEN.length_nt,
        "chemistry": MIPOMERSEN.chemistry,
        "backbone": MIPOMERSEN.backbone,
        "gapmer_design": MIPOMERSEN.gapmer_design,
        "route": MIPOMERSEN.route_of_administration,
        "target_ph": MIPOMERSEN.target_tissue_ph,
        "half_life_hours": MIPOMERSEN.reported_half_life_tissue_hours,
        "mechanism": MIPOMERSEN.mechanism,
        "reference": MIPOMERSEN.key_reference,
    }

    # Inotersen
    table["inotersen"] = {
        "brand_name": INOTERSEN.brand_name,
        "approval_year": INOTERSEN.fda_approval_year,
        "indication": INOTERSEN.indication,
        "length_nt": INOTERSEN.length_nt,
        "chemistry": INOTERSEN.chemistry,
        "backbone": INOTERSEN.backbone,
        "gapmer_design": INOTERSEN.gapmer_design,
        "route": INOTERSEN.route_of_administration,
        "target_ph": INOTERSEN.target_tissue_ph,
        "half_life_hours": INOTERSEN.reported_half_life_tissue_hours,
        "mechanism": INOTERSEN.mechanism,
        "reference": INOTERSEN.key_reference,
    }

    # MRL-ASO-001
    table["mrl_aso_001"] = {
        "brand_name": "N/A (preclinical)",
        "approval_year": 0,
        "indication": "Canine visceral leishmaniasis",
        "length_nt": mrl.length_nt,
        "chemistry": mrl.chemistry,
        "backbone": mrl.backbone,
        "gapmer_design": mrl.gapmer_design,
        "route": mrl.route_of_administration,
        "target_ph": mrl.target_tissue_ph,
        "half_life_hours": mrl.estimated_half_life_hours,
        "mechanism": mrl.mechanism,
        "reference": "This work (computational)",
        "unique_advantage": mrl.unique_advantage,
    }

    return table
