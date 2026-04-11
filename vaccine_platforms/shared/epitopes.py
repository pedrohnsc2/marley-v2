"""Fonte unica de verdade (SSOT) para os 11 epitopos canonicos.

Os epitopos foram selecionados pelo pipeline de vacinologia reversa do Marley
a partir do proteoma de Leishmania infantum JPCM5. Cada peptideo foi predito
como ligante forte de alelos DLA (Dog Leukocyte Antigen) usando NetMHCpan,
com IC50 < 500 nM e rank percentual < 2% na maioria dos casos.

A arquitetura do construto multi-epitopo segue o padrao classico:
    [Peptideo Sinal]-[Adjuvante]-EAAAK-[Epi1]-AAY-[Epi2]-...-AAY-[EpiN]

O linker AAY entre epitopos CTL e padrao na literatura de vacinas
multi-epitopo -- a clivagem pelo proteassoma ocorre preferencialmente
apos tirosina (Y), liberando os epitopos individuais para apresentacao
via MHC-I.

IMPORTANTE: Estes dados sao IMUTAVEIS. Qualquer modulo que precise dos
epitopos deve importar deste arquivo. Nunca copie as sequencias para
outro local.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True, slots=True)
class Epitope:
    """Epitopo CTL predito para Canis lupus familiaris (DLA).

    Campos correspondem diretamente ao output do NetMHCpan filtrado
    pelo pipeline Marley (modulo 04 - construct design).

    Atributos:
        peptide: sequencia aminoacidica do epitopo (9-mer tipico)
        gene_id: identificador do gene no genoma de L. infantum JPCM5
        gene_name: produto genico anotado no TriTrypDB
        allele: alelo DLA que apresenta este epitopo (MHC-I canino)
        ic50: afinidade de ligacao predita em nanomolar (quanto menor, melhor)
        rank: rank percentual no NetMHCpan (< 0.5% = forte, < 2% = fraco)
        position: posicao inicial do epitopo na proteina de origem
    """

    peptide: str
    gene_id: str
    gene_name: str
    allele: str
    ic50: float
    rank: float
    position: int


# ---------------------------------------------------------------------------
# Os 11 epitopos unicos selecionados pelo pipeline Marley.
#
# Nota: o construto original possui 15 entradas no construct_card.json,
# porem LLTANVCYK aparece 5 vezes (de 5 copias do gene CPB em cromossomos
# diferentes). Aqui mantemos apenas a entrada representativa (gene LINF_080014800)
# para evitar redundancia. O construto vacinal usa todas as 15 copias,
# pois a repeticao de LLTANVCYK e intencional para aumentar a imunogenicidade.
# ---------------------------------------------------------------------------

EPITOPES: Final[tuple[Epitope, ...]] = (
    # --- A ordem segue o construto vacinal original (vaccine_construct.fasta) ---
    # A disposicao foi otimizada pelo modulo 04 para minimizar energia livre
    # de juncao entre epitopos adjacentes (heuristica greedy de ordenacao).

    # 1. emp24/gp25L/p24 family/GOLD (cromossomo 24)
    # Proteina da familia p24/GOLD, envolvida no transporte vesicular do
    # parasita. Afinidade muito alta para DLA-88 (IC50 = 11.16 nM).
    Epitope(
        peptide="RMMRSLTPF",
        gene_id="LINF_240013900",
        gene_name="emp24/gp25L/p24 family/GOLD",
        allele="DLA-8803401",
        ic50=11.16,
        rank=0.01,
        position=56,
    ),

    # 2. emp24/gp25L/p24 family/GOLD (cromossomo 24)
    Epitope(
        peptide="YIYETFASM",
        gene_id="LINF_240013900",
        gene_name="emp24/gp25L/p24 family/GOLD",
        allele="DLA-8850101",
        ic50=18.78,
        rank=0.02,
        position=160,
    ),

    # 3. Proteina hipotetica conservada (cromossomo 23)
    # Alvo interessante: proteina conservada sem funcao anotada, mas com
    # multiplos epitopos fortes -- sugere importancia funcional no parasita.
    Epitope(
        peptide="SLMCVFYFK",
        gene_id="LINF_230010300",
        gene_name="hypothetical protein - conserved",
        allele="DLA-8850801",
        ic50=33.78,
        rank=0.02,
        position=173,
    ),

    # 4. Proteina hipotetica conservada (cromossomo 23)
    Epitope(
        peptide="YLAALVPAL",
        gene_id="LINF_230010300",
        gene_name="hypothetical protein - conserved",
        allele="DLA-8850101",
        ic50=38.17,
        rank=0.02,
        position=200,
    ),

    # 5. Proibina (cromossomo 16) -- chaperone mitocondrial essencial
    # Altamente conservada entre especies de Leishmania.
    Epitope(
        peptide="LIIEDLSLV",
        gene_id="LINF_160022200",
        gene_name="prohibitin",
        allele="DLA-8850101",
        ic50=64.12,
        rank=0.04,
        position=158,
    ),

    # 6. Proteina hipotetica conservada (cromossomo 23)
    Epitope(
        peptide="FAFSVSARR",
        gene_id="LINF_230010300",
        gene_name="hypothetical protein - conserved",
        allele="DLA-8850801",
        ic50=64.65,
        rank=0.04,
        position=252,
    ),

    # 7. emp24/gp25L/p24 family/GOLD (cromossomo 24)
    Epitope(
        peptide="MILGTFVRL",
        gene_id="LINF_240013900",
        gene_name="emp24/gp25L/p24 family/GOLD",
        allele="DLA-8850101",
        ic50=72.87,
        rank=0.05,
        position=207,
    ),

    # 8. Proibina (cromossomo 16)
    Epitope(
        peptide="MQNVTFVPK",
        gene_id="LINF_160022200",
        gene_name="prohibitin",
        allele="DLA-8850801",
        ic50=75.37,
        rank=0.04,
        position=249,
    ),

    # 9. GP63 / Leishmanolysin (cromossomo 10)
    # Metalopeptidase de superficie mais abundante do parasita. Alvo classico
    # da vacinologia anti-Leishmania, com papel na evasao imune e invasao
    # de macrofagos.
    Epitope(
        peptide="RILESISNV",
        gene_id="LINF_100010400",
        gene_name="GP63 - leishmanolysin",
        allele="DLA-8850101",
        ic50=97.44,
        rank=0.06,
        position=277,
    ),

    # 10. Proibina (cromossomo 16)
    Epitope(
        peptide="ILYNKISGL",
        gene_id="LINF_160022200",
        gene_name="prohibitin",
        allele="DLA-8803401",
        ic50=97.63,
        rank=0.02,
        position=32,
    ),

    # 11. Cisteina peptidase B / CPB (cromossomo 08)
    # Fator de virulencia bem caracterizado, envolvido na degradacao de
    # proteinas do hospedeiro e modulacao da resposta imune. O gene possui
    # multiplas copias em tandem (5 copias em L. infantum JPCM5), todas
    # codificando o mesmo epitopo nesta regiao.
    Epitope(
        peptide="LLTANVCYK",
        gene_id="LINF_080014800",
        gene_name="cysteine peptidase B (CPB)",
        allele="DLA-8850801",
        ic50=117.96,
        rank=0.09,
        position=379,
    ),
)

# ---------------------------------------------------------------------------
# Constantes do construto vacinal
# ---------------------------------------------------------------------------

# Linkers proteoliticos para construtos multi-epitopo.
# AAY: clivado pelo proteassoma, libera epitopos CTL (MHC classe I)
# GPGPG: clivado por catepsinas, libera epitopos HTL (MHC classe II)
LINKERS: Final[dict[str, str]] = {
    "CTL": "AAY",
    "HTL": "GPGPG",
}

# Adjuvante L7/L12 (proteina ribosomal bacteriana).
# Potente ativador de TLR4, induz resposta Th1 -- essencial contra
# Leishmania, onde imunidade celular (nao humoral) e protetora.
ADJUVANT: Final[str] = "L7/L12"

# Sequencia completa do adjuvante L7/L12 ribosomal, extraida do construto
# vacinal validado (residuos 24-153 do construto final).
ADJUVANT_SEQUENCE: Final[str] = (
    "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGA"
    "APAGAAVEAAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAK"
    "DLVDGAPKPLLEKVAKEAADEAKAKLEAAGATVTVK"
)

# Linker entre adjuvante e cassete de epitopos.
# EAAAK forma alfa-helice rigida, separando dominio do adjuvante
# dos epitopos para evitar interferencia estereoquimica.
ADJUVANT_EPITOPE_LINKER: Final[str] = "EAAAK"

# Peptideo sinal tPA (tissue Plasminogen Activator) humano.
# Direciona a proteina para a via secretoria, aumentando a
# apresentacao via MHC-I por cross-priming e a secrecao do
# construto pela celula transfectada.
SIGNAL_PEPTIDE_TPA: Final[str] = "MDAMKRGLCCVLLLCGAVFVSAS"


# ---------------------------------------------------------------------------
# Funcoes auxiliares
# ---------------------------------------------------------------------------

def get_epitope_sequences() -> list[str]:
    """Retorna lista de sequencias peptidicas dos 11 epitopos unicos.

    A ordem segue a mesma do construto vacinal original, que foi
    otimizada para minimizar energia livre de juncao entre epitopos
    adjacentes (heuristica greedy do modulo 04).
    """
    return [ep.peptide for ep in EPITOPES]


def get_unique_source_genes() -> dict[str, list[Epitope]]:
    """Agrupa epitopos por gene de origem.

    Util para analise de cobertura antigenica -- quantos alvos
    proteicos distintos do parasita estao representados no construto.
    """
    genes: dict[str, list[Epitope]] = {}
    for ep in EPITOPES:
        genes.setdefault(ep.gene_name, []).append(ep)
    return genes


def get_construct_sequence(
    *,
    include_signal: bool = True,
    cpb_repeats: int = 5,
) -> str:
    """Monta a sequencia completa do construto multi-epitopo.

    Arquitetura:
        [tPA signal]-[L7/L12 adjuvante]-EAAAK-[Epi1]-AAY-[Epi2]-...-AAY-[EpiN]

    O parametro cpb_repeats controla quantas copias do epitopo LLTANVCYK
    sao incluidas (default=5, como no construto original). Biologicamente,
    a repeticao amplifica a resposta CTL contra CPB, que e fator de
    virulencia critico de Leishmania.

    Args:
        include_signal: se True, inclui o peptideo sinal tPA no N-terminal
        cpb_repeats: numero de repeticoes do epitopo LLTANVCYK de CPB (1-5)

    Returns:
        Sequencia aminoacidica completa do construto vacinal
    """
    if not 1 <= cpb_repeats <= 5:
        raise ValueError(
            f"cpb_repeats deve ser entre 1 e 5, recebeu {cpb_repeats}"
        )

    parts: list[str] = []

    # Peptideo sinal (opcional -- removido no construto maduro)
    if include_signal:
        parts.append(SIGNAL_PEPTIDE_TPA)

    # Adjuvante L7/L12
    parts.append(ADJUVANT_SEQUENCE)

    # Linker rigido entre adjuvante e cassete de epitopos
    parts.append(ADJUVANT_EPITOPE_LINKER)

    # Cassete de epitopos: os primeiros 10 epitopos unicos (exceto CPB)
    # seguidos por N copias do epitopo CPB, todos separados por AAY
    linker = LINKERS["CTL"]
    non_cpb = [ep.peptide for ep in EPITOPES if ep.peptide != "LLTANVCYK"]
    cpb_peptide = "LLTANVCYK"

    epitope_list = non_cpb + [cpb_peptide] * cpb_repeats
    parts.append(linker.join(epitope_list))

    return "".join(parts)


def get_epitope_by_gene(gene_id: str) -> list[Epitope]:
    """Busca epitopos por ID do gene de origem.

    Args:
        gene_id: identificador do gene (ex: 'LINF_240013900')

    Returns:
        Lista de epitopos derivados deste gene (pode ser vazia)
    """
    return [ep for ep in EPITOPES if ep.gene_id == gene_id]


def get_epitope_by_allele(allele: str) -> list[Epitope]:
    """Busca epitopos por alelo DLA apresentador.

    Util para avaliar cobertura allelica -- idealmente, o construto
    deve ter epitopos para multiplos alelos DLA para cobrir a
    diversidade genetica da populacao canina.

    Args:
        allele: nome do alelo DLA (ex: 'DLA-8850101')

    Returns:
        Lista de epitopos apresentados por este alelo
    """
    return [ep for ep in EPITOPES if ep.allele == allele]
