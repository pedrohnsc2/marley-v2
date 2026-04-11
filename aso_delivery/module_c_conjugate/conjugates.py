"""Modelos de conjugados para entrega de ASO a macrofagos.

Define propriedades fisico-quimicas, receptores-alvo e scores
para cada estrategia de conjugacao avaliada no Modulo C.

A logica biologica: macrofagos infectados por L. infantum expressam
receptores de superficie especificos que podem ser explorados para
entrega seletiva do ASO ao fagolisossomo. O receptor manose (MRC1/CD206)
e particularmente atrativo porque:
1. Alta expressao em macrofagos (especialmente ativados)
2. Via de internalizacao termina no lisossomo/fagolisossomo
3. Carreamento de manose e sinteticamente acessivel
4. Ausente em maioria das celulas nao-fagociticas

Cada conjugado e modelado como dataclass imutavel com propriedades
calculadas a partir de constantes publicadas na literatura.

Referencias:
- East L, Isacke CM (2002) Biochim Biophys Acta 1572(2-3):364-386 — MRC1
- Nair JK et al. (2014) JACS 136(49):16958-16961 — GalNAc-ASO
- Wolfrum C et al. (2007) Nat Biotechnol 25(10):1149-1157 — colesterol-ASO
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO geral
- Sehgal A et al. (2013) J Hepatol 59(6):1354-1359 — GalNAc hepatocitos
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes de receptores em macrofagos
# Expressao relativa: 1.0 = expressao basal em macrofago nao-ativado
# Valores derivados de dados de expressao publicados (FPKM / TPM normalizados)
# ---------------------------------------------------------------------------

# MRC1 (CD206) — receptor manose
# Alta expressao em macrofagos M2/alternativamente ativados
# Ref: Martinez-Pomares L (2012) J Leukoc Biol 92(6):1177-1186
MRC1_EXPRESSION_MACROPHAGE: Final[float] = 8.5      # alta
MRC1_EXPRESSION_HEPATOCYTE: Final[float] = 0.1      # negligivel
MRC1_EXPRESSION_FIBROBLAST: Final[float] = 0.2      # negligivel

# LDLR — receptor LDL (via colesterol)
# Expressao moderada em macrofagos, aumenta durante infeccao
# Ref: Goldstein JL, Brown MS (2009) Arterioscler Thromb Vasc Biol 29(4):431-438
LDLR_EXPRESSION_MACROPHAGE: Final[float] = 3.2      # moderada
LDLR_EXPRESSION_HEPATOCYTE: Final[float] = 9.0      # muito alta
LDLR_EXPRESSION_FIBROBLAST: Final[float] = 1.5      # basal

# ASGPR — receptor asialoglicoproteina (via GalNAc)
# Expresso EXCLUSIVAMENTE em hepatocitos — ausente em macrofagos
# Ref: Spiess M (1990) Biochemistry 29(43):10009-10018
ASGPR_EXPRESSION_MACROPHAGE: Final[float] = 0.0     # ZERO
ASGPR_EXPRESSION_HEPATOCYTE: Final[float] = 9.5     # altissima
ASGPR_EXPRESSION_FIBROBLAST: Final[float] = 0.0     # ZERO

# Scavenger receptor classe A (SR-A / CD204) — via palmitato
# Alta expressao em macrofagos, reconhece lipideos modificados
# Ref: Peiser L et al. (2002) Curr Opin Immunol 14(1):123-128
SRA_EXPRESSION_MACROPHAGE: Final[float] = 7.0       # alta
SRA_EXPRESSION_HEPATOCYTE: Final[float] = 0.5       # baixa
SRA_EXPRESSION_FIBROBLAST: Final[float] = 0.3       # baixa


# ---------------------------------------------------------------------------
# Constantes fisico-quimicas dos conjugados
# ---------------------------------------------------------------------------

# Peso molecular do ASO nu (fosforotioato, 25-mer)
# MW_PS ~= 330 Da/nt para PS backbone
# Ref: Crooke ST et al. (2017) — MW medio de PS-ASO
MW_PER_NT_PS: Final[float] = 330.0       # Da/nt
ASO_LENGTH_NT: Final[int] = 25
MW_ASO_NAKED: Final[float] = MW_PER_NT_PS * ASO_LENGTH_NT  # ~8250 Da

# logP do ASO nu (PS backbone, altamente hidrofilico e anionico)
# Oligonucleotideos fosforotioato sao muito polares, logP negativo
# Ref: Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51
LOGP_ASO_NAKED: Final[float] = -3.5

# Constantes de ligantes/conjugados
MW_MANNOSE: Final[float] = 180.16       # D-manose, Da
MW_TRIMANNOSE: Final[float] = 504.44    # trimanose cluster, Da
MW_C6_LINKER: Final[float] = 114.14     # acido 6-aminohexanoico, Da
MW_CHOLESTEROL: Final[float] = 386.65   # colesterol, Da
MW_PALMITATE: Final[float] = 256.42     # acido palmitico, Da
MW_GALNAC: Final[float] = 221.21       # N-acetilgalactosamina, Da
MW_GALNAC_TRIMER: Final[float] = 663.63  # cluster trivalente, Da

# logP de fragmentos conjugados
LOGP_MANNOSE: Final[float] = -2.7       # muito hidrofílico
LOGP_TRIMANNOSE: Final[float] = -4.5    # mais hidrofilico ainda
LOGP_CHOLESTEROL: Final[float] = 7.0    # altamente lipofilico
LOGP_PALMITATE: Final[float] = 6.4      # altamente lipofilico
LOGP_GALNAC: Final[float] = -3.0        # hidrofilico
LOGP_C6_LINKER: Final[float] = 0.5      # levemente lipofilico


# ---------------------------------------------------------------------------
# Kd aparente para receptores (valores simplificados)
# Kd menor = maior afinidade
# ---------------------------------------------------------------------------

# Kd da manose para MRC1: ~2-10 microM para monovalente
# Cluster trivalente melhora para ~50-200 nM (efeito avidez)
# Ref: Taylor ME et al. (1992) J Biol Chem 267(3):1719-1726
KD_MANNOSE_MRC1_UM: Final[float] = 5.0       # microM (monovalente)
KD_TRIMANNOSE_MRC1_UM: Final[float] = 0.1    # microM (trivalente, avidez)

# Kd do colesterol para LDLR: mediada por LDL, ~10-50 nM
# Mas conjugado direto: ~1-5 microM (sem contexto lipoproteico)
# Ref: Wolfrum C et al. (2007) Nat Biotechnol 25(10):1149-1157
KD_CHOLESTEROL_LDLR_UM: Final[float] = 2.0   # microM

# Kd do palmitato para SR-A: ~1-10 microM
# Ref: Peiser L et al. (2002) Curr Opin Immunol 14(1):123-128
KD_PALMITATE_SRA_UM: Final[float] = 5.0      # microM

# Kd do GalNAc trimer para ASGPR: ~2-5 nM (altissima afinidade)
# Ref: Nair JK et al. (2014) JACS 136(49):16958-16961
KD_GALNAC_ASGPR_UM: Final[float] = 0.003     # microM (3 nM)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReceptorProfile:
    """Perfil de expressao de um receptor em diferentes tipos celulares.

    Permite calcular a seletividade de um conjugado para macrofagos
    vs outros tipos celulares (especialmente hepatocitos, onde a
    maioria dos ASOs se acumula por efeito de primeira passagem).
    """

    receptor_name: str
    gene_symbol: str
    expression_macrophage: float     # expressao relativa em macrofagos
    expression_hepatocyte: float     # expressao relativa em hepatocitos
    expression_fibroblast: float     # expressao relativa em fibroblastos
    internalization_pathway: str     # via de internalizacao
    endpoint_compartment: str        # compartimento final apos internalizacao

    @property
    def macrophage_selectivity_vs_hepatocyte(self) -> float:
        """Razao de expressao macrofago/hepatocito.

        Valores > 1 indicam seletividade para macrofagos.
        Valores < 1 indicam acumulo preferencial no figado.
        """
        if self.expression_hepatocyte == 0:
            return float("inf")
        return self.expression_macrophage / self.expression_hepatocyte


@dataclass(frozen=True)
class ConjugateProperties:
    """Propriedades fisico-quimicas de um conjugado ASO-ligante.

    Calcula MW total, logP estimado e constante de associacao
    a membrana a partir das propriedades dos componentes.
    """

    name: str
    ligand_name: str
    linker: str                      # tipo de linker (ex: "C6 amino")
    attachment_site: str             # posicao de conjugacao (3', 5', interna)
    mw_ligand_da: float              # MW do ligante em Daltons
    mw_linker_da: float              # MW do linker em Daltons
    logp_ligand: float               # logP do fragmento ligante
    target_receptor: str             # receptor alvo
    kd_receptor_um: float            # Kd para o receptor em microM
    receptor_expression: float       # expressao do receptor em macrofagos (relativa)
    synthesis_complexity: int        # escala 1-5 (1 = simples)
    estimated_cost_per_mg_usd: float # custo estimado por mg
    literature_precedent: str        # referencia a precedentes clinicos/pre-clinicos

    @property
    def mw_total_da(self) -> float:
        """MW total do conjugado ASO + linker + ligante."""
        return MW_ASO_NAKED + self.mw_linker_da + self.mw_ligand_da

    @property
    def logp_estimated(self) -> float:
        """logP estimado do conjugado completo.

        Aproximacao simplificada: logP(conjugado) ~ logP(ASO) + contribuicao
        do ligante ponderada pela fracao de massa. Nao e aditivo simples
        porque o ASO domina a hidrofilicidade, mas o conjugado lipofilico
        aumenta a associacao a membrana.
        """
        # Fracao de massa do ligante no conjugado
        f_ligand = (self.mw_ligand_da + self.mw_linker_da) / self.mw_total_da
        # Contribuicao ponderada
        return LOGP_ASO_NAKED * (1 - f_ligand) + self.logp_ligand * f_ligand

    @property
    def membrane_association_constant(self) -> float:
        """Constante de associacao a membrana (adimensional, relativa ao ASO nu).

        Modelo simplificado: Ka_membrane ~ 10^(delta_logP / 2)
        onde delta_logP = logP(conjugado) - logP(ASO nu).
        Conjugados lipofilicos aumentam a interacao com membranas.
        """
        delta_logp = self.logp_estimated - LOGP_ASO_NAKED
        return 10.0 ** (delta_logp / 2.0)


@dataclass(frozen=True)
class ConjugateScore:
    """Score final de um conjugado para entrega a macrofagos.

    Combina afinidade pelo receptor, expressao do receptor em macrofagos,
    seletividade vs hepatocitos, e viabilidade sintetica em um score
    composto que permite comparacao direta entre estrategias.
    """

    conjugate_name: str
    receptor_score: float            # log10(1/Kd) normalizado [0-1]
    expression_score: float          # expressao normalizada [0-1]
    selectivity_score: float         # seletividade mac/hep [0-1]
    pathway_score: float             # 1.0 se termina no fagolisossomo, 0.5 se endossomo
    synthesis_feasibility: float     # (5 - complexity) / 4.0, [0-1]
    cost_score: float                # normalizado inversamente ao custo [0-1]
    recommendation_score: float      # score composto final [0-1]
    suitable_for_macrophage: bool    # recomendacao binaria
    rationale: str                   # justificativa textual


# ---------------------------------------------------------------------------
# Funcoes de calculo
# ---------------------------------------------------------------------------


def build_receptor_profiles() -> dict[str, ReceptorProfile]:
    """Constroi perfis de todos os receptores relevantes.

    Cada receptor e avaliado quanto a expressao em macrofagos,
    hepatocitos e fibroblastos, alem da via de internalizacao.

    Returns:
        Dicionario receptor_name -> ReceptorProfile.
    """
    return {
        "MRC1": ReceptorProfile(
            receptor_name="Mannose receptor C-type 1",
            gene_symbol="MRC1 (CD206)",
            expression_macrophage=MRC1_EXPRESSION_MACROPHAGE,
            expression_hepatocyte=MRC1_EXPRESSION_HEPATOCYTE,
            expression_fibroblast=MRC1_EXPRESSION_FIBROBLAST,
            internalization_pathway="Clathrin-mediated endocytosis -> early endosome -> lysosome",
            endpoint_compartment="lysosome/phagolysosome",
        ),
        "LDLR": ReceptorProfile(
            receptor_name="Low-density lipoprotein receptor",
            gene_symbol="LDLR",
            expression_macrophage=LDLR_EXPRESSION_MACROPHAGE,
            expression_hepatocyte=LDLR_EXPRESSION_HEPATOCYTE,
            expression_fibroblast=LDLR_EXPRESSION_FIBROBLAST,
            internalization_pathway="Clathrin-mediated endocytosis -> endosome -> lysosome",
            endpoint_compartment="lysosome",
        ),
        "ASGPR": ReceptorProfile(
            receptor_name="Asialoglycoprotein receptor",
            gene_symbol="ASGPR1/ASGPR2",
            expression_macrophage=ASGPR_EXPRESSION_MACROPHAGE,
            expression_hepatocyte=ASGPR_EXPRESSION_HEPATOCYTE,
            expression_fibroblast=ASGPR_EXPRESSION_FIBROBLAST,
            internalization_pathway="Clathrin-mediated endocytosis -> endosome -> lysosome",
            endpoint_compartment="hepatocyte lysosome (liver-specific)",
        ),
        "SR-A": ReceptorProfile(
            receptor_name="Scavenger receptor class A",
            gene_symbol="MSR1 (CD204)",
            expression_macrophage=SRA_EXPRESSION_MACROPHAGE,
            expression_hepatocyte=SRA_EXPRESSION_HEPATOCYTE,
            expression_fibroblast=SRA_EXPRESSION_FIBROBLAST,
            internalization_pathway="Non-clathrin endocytosis -> endosome -> lysosome",
            endpoint_compartment="lysosome",
        ),
    }


def build_conjugate_properties() -> dict[str, ConjugateProperties]:
    """Constroi propriedades fisico-quimicas de todos os conjugados avaliados.

    Cinco estrategias: manose (mono), trimanose (cluster), colesterol,
    palmitato e GalNAc. O ASO nu tambem e incluido como controle.

    Returns:
        Dicionario conjugate_name -> ConjugateProperties.
    """
    return {
        "mannose_mono": ConjugateProperties(
            name="Mannose-ASO (monovalent)",
            ligand_name="D-mannose",
            linker="C6 aminohexyl",
            attachment_site="3' terminal",
            mw_ligand_da=MW_MANNOSE,
            mw_linker_da=MW_C6_LINKER,
            logp_ligand=LOGP_MANNOSE,
            target_receptor="MRC1 (CD206)",
            kd_receptor_um=KD_MANNOSE_MRC1_UM,
            receptor_expression=MRC1_EXPRESSION_MACROPHAGE,
            synthesis_complexity=2,
            estimated_cost_per_mg_usd=150.0,
            literature_precedent=(
                "Kawakami S et al. (2008) J Control Release 131(3):153-158 — "
                "mannosylated liposomes for macrophage delivery in leishmaniasis"
            ),
        ),
        "trimannose_cluster": ConjugateProperties(
            name="Trimannose-ASO (trivalent cluster)",
            ligand_name="alpha-1,3/alpha-1,6-trimannose",
            linker="C6 aminohexyl + branching scaffold",
            attachment_site="3' terminal",
            mw_ligand_da=MW_TRIMANNOSE,
            mw_linker_da=MW_C6_LINKER * 1.5,  # scaffold adicional
            logp_ligand=LOGP_TRIMANNOSE,
            target_receptor="MRC1 (CD206)",
            kd_receptor_um=KD_TRIMANNOSE_MRC1_UM,
            receptor_expression=MRC1_EXPRESSION_MACROPHAGE,
            synthesis_complexity=4,
            estimated_cost_per_mg_usd=450.0,
            literature_precedent=(
                "Taylor ME et al. (1992) J Biol Chem 267(3):1719-1726 — "
                "multivalent mannose dramatically improves MRC1 binding; "
                "Irache JM et al. (2008) J Control Release 128(1):15-25 — "
                "mannosylated nanoparticles for macrophage targeting"
            ),
        ),
        "cholesterol": ConjugateProperties(
            name="Cholesterol-ASO",
            ligand_name="Cholesterol",
            linker="C6 aminohexyl",
            attachment_site="3' terminal",
            mw_ligand_da=MW_CHOLESTEROL,
            mw_linker_da=MW_C6_LINKER,
            logp_ligand=LOGP_CHOLESTEROL,
            target_receptor="LDLR",
            kd_receptor_um=KD_CHOLESTEROL_LDLR_UM,
            receptor_expression=LDLR_EXPRESSION_MACROPHAGE,
            synthesis_complexity=2,
            estimated_cost_per_mg_usd=120.0,
            literature_precedent=(
                "Wolfrum C et al. (2007) Nat Biotechnol 25(10):1149-1157 — "
                "cholesterol-conjugated siRNA silences gene in vivo; "
                "Soutschek J et al. (2004) Nature 432(7014):173-178 — "
                "cholesterol-siRNA systemic delivery"
            ),
        ),
        "palmitate": ConjugateProperties(
            name="Palmitate-ASO",
            ligand_name="Palmitic acid (C16:0)",
            linker="Amide bond (direct)",
            attachment_site="3' terminal",
            mw_ligand_da=MW_PALMITATE,
            mw_linker_da=18.0,  # perda de H2O na amidacao
            logp_ligand=LOGP_PALMITATE,
            target_receptor="SR-A (CD204)",
            kd_receptor_um=KD_PALMITATE_SRA_UM,
            receptor_expression=SRA_EXPRESSION_MACROPHAGE,
            synthesis_complexity=1,
            estimated_cost_per_mg_usd=80.0,
            literature_precedent=(
                "Toulme JJ et al. (1994) Biochimie 76(3-4):153-154 — "
                "palmitate conjugation for macrophage-targeted antisense; "
                "Bijsterbosch MK et al. (2000) Nucleic Acids Res 28(14):2717-2725 — "
                "lipid-conjugated ODN macrophage uptake"
            ),
        ),
        "galnac": ConjugateProperties(
            name="GalNAc-ASO (trivalent)",
            ligand_name="N-acetylgalactosamine trimer",
            linker="Tris-linker + PEG spacer",
            attachment_site="3' terminal",
            mw_ligand_da=MW_GALNAC_TRIMER,
            mw_linker_da=MW_C6_LINKER * 2.0,  # tris scaffold + PEG
            logp_ligand=LOGP_GALNAC,
            target_receptor="ASGPR",
            kd_receptor_um=KD_GALNAC_ASGPR_UM,
            receptor_expression=ASGPR_EXPRESSION_MACROPHAGE,  # ZERO
            synthesis_complexity=5,
            estimated_cost_per_mg_usd=800.0,
            literature_precedent=(
                "Nair JK et al. (2014) JACS 136(49):16958-16961 — "
                "GalNAc-siRNA platform (Alnylam); "
                "Prakash TP et al. (2014) Nucleic Acids Res 42(13):8796-8807 — "
                "GalNAc-ASO hepatocyte delivery. "
                "NOTE: ASGPR is liver-specific, NOT expressed on macrophages"
            ),
        ),
    }


def compute_uptake_fold_vs_naked(
    kd_um: float,
    receptor_expression: float,
    membrane_ka: float,
) -> float:
    """Calcula captacao esperada relativa ao ASO nu (fold change).

    Modelo simplificado que combina tres contribuicoes:
    1. Afinidade pelo receptor (menor Kd = maior contribuicao)
    2. Expressao do receptor na celula alvo
    3. Associacao a membrana do conjugado

    A captacao do ASO nu por macrofagos ocorre via endocitose
    nao-especifica e interacao com scavenger receptors (PS backbone
    e polianionico). Estimamos captacao basal = 1.0x.

    Conjugados melhoram a captacao por adicionar mecanismo receptor-mediado.

    Args:
        kd_um: Constante de dissociacao em microM.
        receptor_expression: Expressao relativa do receptor (0-10).
        membrane_ka: Constante de associacao a membrana (relativa ao ASO nu).

    Returns:
        Fold change de captacao vs ASO nu.
    """
    # Se receptor nao e expresso, nao ha ganho receptor-mediado
    if receptor_expression == 0:
        # Apenas efeito de membrana (muito pequeno para conjugados hidrofilicos)
        return max(membrane_ka * 0.1, 0.1)

    # Score de afinidade: inversamente proporcional ao Kd
    # Escala logaritmica porque Kd varia em ordens de magnitude
    # Kd = 0.003 uM (GalNAc) -> affinity_score = 2.52
    # Kd = 5.0 uM (mannose mono) -> affinity_score = -0.70
    # Normalizado para que Kd = 1 uM -> affinity_score = 0
    affinity_score = math.log10(1.0 / kd_um)

    # Contribuicao receptor-mediada: expressao * afinidade
    # Normalizado para que receptor_expression=10, Kd=0.001 -> ~50x
    receptor_contribution = receptor_expression * (1 + affinity_score) * 0.5

    # Contribuicao de membrana (menor efeito)
    membrane_contribution = membrane_ka * 0.2

    # Captacao total = 1 (basal) + contribuicoes
    fold = 1.0 + max(receptor_contribution, 0) + membrane_contribution

    return round(fold, 2)


def compute_conjugate_score(
    conjugate: ConjugateProperties,
    receptor: ReceptorProfile,
    uptake_fold: float,
    max_uptake_fold: float,
) -> ConjugateScore:
    """Calcula score composto para um conjugado.

    O score final e uma media ponderada de seis componentes:
    - receptor_score (25%): afinidade pelo receptor alvo
    - expression_score (25%): expressao do receptor em macrofagos
    - selectivity_score (20%): seletividade macrofago vs hepatocito
    - pathway_score (15%): conjugado entrega ao compartimento correto?
    - synthesis_feasibility (10%): viabilidade sintetica
    - cost_score (5%): custo relativo

    Args:
        conjugate: Propriedades do conjugado.
        receptor: Perfil do receptor alvo.
        uptake_fold: Captacao estimada (fold vs naked).
        max_uptake_fold: Maior fold entre todos os conjugados (para normalizacao).

    Returns:
        ConjugateScore com score detalhado e recomendacao.
    """
    # 1. Receptor score: log10(1/Kd) normalizado para [0,1]
    # Kd em uM: menor = melhor. log10(1/0.003) = 2.52, log10(1/5) = -0.70
    # Normalizo para range esperado [-1, 3]
    raw_affinity = math.log10(1.0 / conjugate.kd_receptor_um)
    receptor_score = min(max((raw_affinity + 1.0) / 4.0, 0.0), 1.0)

    # 2. Expression score: expressao em macrofago / 10 (escala 0-10)
    expression_score = min(conjugate.receptor_expression / 10.0, 1.0)

    # 3. Selectivity score: razao macrofago/hepatocito
    selectivity_ratio = receptor.macrophage_selectivity_vs_hepatocyte
    if selectivity_ratio == float("inf"):
        # Receptor so existe em macrofagos (caso SR-A) -> score maximo
        # Mas tambem pode ser ASGPR com expressao 0 em macrofago -> tratar
        if receptor.expression_macrophage > 0:
            selectivity_score = 1.0
        else:
            selectivity_score = 0.0
    else:
        # Log scale: ratio=1 -> 0.5, ratio=10 -> ~0.75, ratio=100 -> 1.0
        selectivity_score = min(math.log10(max(selectivity_ratio, 0.01) + 1) / 2.0, 1.0)

    # 4. Pathway score: destino final e fagolisossomo?
    endpoint = receptor.endpoint_compartment.lower()
    if "phagolysosome" in endpoint or "lysosome" in endpoint:
        pathway_score = 1.0
    elif "endosome" in endpoint:
        pathway_score = 0.5
    else:
        pathway_score = 0.2

    # Caso especial: ASGPR entrega ao HEPATOCITO, nao ao macrofago
    if receptor.gene_symbol.startswith("ASGPR") and receptor.expression_macrophage == 0:
        pathway_score = 0.0

    # 5. Synthesis feasibility: (5 - complexity) / 4
    synthesis_feasibility = (5 - conjugate.synthesis_complexity) / 4.0

    # 6. Cost score: inversamente proporcional ao custo
    # $80/mg -> 1.0, $800/mg -> 0.1
    cost_score = min(80.0 / conjugate.estimated_cost_per_mg_usd, 1.0)

    # Score composto (ponderado)
    recommendation_score = (
        0.25 * receptor_score
        + 0.25 * expression_score
        + 0.20 * selectivity_score
        + 0.15 * pathway_score
        + 0.10 * synthesis_feasibility
        + 0.05 * cost_score
    )

    # Criterios de adequacao para macrofago
    # FALHA se: expressao do receptor = 0, OU selectivity < 0.3, OU pathway = 0
    suitable = (
        expression_score > 0.0
        and selectivity_score > 0.1
        and pathway_score > 0.0
        and uptake_fold > 1.5
    )

    # Gerar justificativa textual
    rationale = _generate_rationale(
        conjugate, receptor, uptake_fold, suitable,
        receptor_score, expression_score, selectivity_score,
    )

    return ConjugateScore(
        conjugate_name=conjugate.name,
        receptor_score=round(receptor_score, 4),
        expression_score=round(expression_score, 4),
        selectivity_score=round(selectivity_score, 4),
        pathway_score=round(pathway_score, 4),
        synthesis_feasibility=round(synthesis_feasibility, 4),
        cost_score=round(cost_score, 4),
        recommendation_score=round(recommendation_score, 4),
        suitable_for_macrophage=suitable,
        rationale=rationale,
    )


def _generate_rationale(
    conjugate: ConjugateProperties,
    receptor: ReceptorProfile,
    uptake_fold: float,
    suitable: bool,
    receptor_score: float,
    expression_score: float,
    selectivity_score: float,
) -> str:
    """Gera justificativa textual para a recomendacao de um conjugado.

    Explica POR QUE o conjugado e adequado ou inadequado para
    entrega a macrofagos, baseado nos scores individuais.
    """
    parts: list[str] = []

    if not suitable:
        # Identificar razao principal de exclusao
        if expression_score == 0:
            parts.append(
                f"UNSUITABLE: {receptor.gene_symbol} is not expressed on macrophages "
                f"(expression = {receptor.expression_macrophage:.1f}/10). "
                f"This receptor is {receptor.endpoint_compartment}-specific."
            )
        elif selectivity_score < 0.1:
            parts.append(
                f"UNSUITABLE: {receptor.gene_symbol} has poor macrophage selectivity "
                f"(macrophage/hepatocyte ratio = "
                f"{receptor.macrophage_selectivity_vs_hepatocyte:.2f})."
            )
        elif uptake_fold < 1.5:
            parts.append(
                f"MARGINAL: {conjugate.name} shows minimal uptake improvement "
                f"({uptake_fold:.1f}x vs naked ASO)."
            )
        return " ".join(parts)

    # Conjugado adequado — detalhar pontos fortes
    parts.append(
        f"SUITABLE: {receptor.gene_symbol} is expressed on macrophages "
        f"({receptor.expression_macrophage:.1f}/10) with "
        f"{'high' if selectivity_score > 0.6 else 'moderate'} selectivity "
        f"vs hepatocytes."
    )

    if conjugate.kd_receptor_um < 1.0:
        parts.append(
            f"Strong receptor affinity (Kd = {conjugate.kd_receptor_um:.3f} uM)."
        )
    else:
        parts.append(
            f"Moderate receptor affinity (Kd = {conjugate.kd_receptor_um:.1f} uM)."
        )

    parts.append(
        f"Expected uptake: {uptake_fold:.1f}x over naked ASO. "
        f"Internalization via {receptor.internalization_pathway}."
    )

    return " ".join(parts)


def compute_naked_aso_score() -> ConjugateScore:
    """Score de referencia para o ASO nu (sem conjugado).

    O ASO nu com backbone PS entra em macrofagos via:
    1. Endocitose nao-especifica (interacao com glicocalice)
    2. Scavenger receptors (PS e polianionico)

    A captacao e real mas ineficiente — define o baseline de 1.0x.

    Returns:
        ConjugateScore para ASO nu (referencia).
    """
    return ConjugateScore(
        conjugate_name="Naked ASO (PS backbone)",
        receptor_score=0.20,
        expression_score=0.50,
        selectivity_score=0.30,
        pathway_score=0.50,
        synthesis_feasibility=1.0,
        cost_score=1.0,
        recommendation_score=0.42,
        suitable_for_macrophage=True,
        rationale=(
            "BASELINE: Naked PS-ASO enters macrophages via non-specific endocytosis "
            "and scavenger receptor interactions. Uptake is real but inefficient "
            "(1.0x reference). The polyanionic PS backbone provides some inherent "
            "macrophage tropism via Class A scavenger receptors."
        ),
    )
