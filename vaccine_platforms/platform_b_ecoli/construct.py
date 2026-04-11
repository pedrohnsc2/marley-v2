"""B1 -- Redesenho do constructo vacinal para expressao em E. coli.

Remove o peptideo sinal tPA (especifico de celulas de mamifero) e adiciona
componentes necessarios para expressao e purificacao bacteriana: His6-tag
no N-terminal, sitio de clivagem TEV, e sitios de restricao NdeI/XhoI para
clonagem em pET-28a(+).

A sequencia do adjuvante L7/L12 e todos os 11 epitopos unicos sao mantidos
intactos.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes biologicas
# ---------------------------------------------------------------------------

# Adjuvante L7/L12 de Mycobacterium tuberculosis -- resposta Th1 robusta
# contra patogenos intracelulares como Leishmania
L7L12_ADJUVANT: Final[str] = (
    "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAEEQ"
    "SEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAKAKL"
    "EAAGATVTVK"
)

# Linkers -- AAY para epitopos CTL (proteassomal), GPGPG para HTL
LINKER_ADJUVANT: Final[str] = "EAAAK"
LINKER_CTL: Final[str] = "AAY"
LINKER_HTL: Final[str] = "GPGPG"

# Componentes de expressao bacteriana
HIS6_TAG: Final[str] = "HHHHHH"
TEV_SITE: Final[str] = "ENLYFQS"  # A Ser que resta apos clivagem pela TEV

# Sitios de restricao para clonagem
NDEI_SITE: Final[str] = "CATATG"  # Contem ATG do inicio da traducao
XHOI_SITE: Final[str] = "CTCGAG"

# Peptideo sinal tPA da plataforma mRNA -- sera REMOVIDO
TPA_SIGNAL: Final[str] = "MDAMKRGLCCVLLLCGAVFVSAS"

# Pesos moleculares medios dos aminoacidos (Da)
# Fonte: Gasteiger et al. 2005, ExPASy ProtParam
AA_MW: Final[dict[str, float]] = {
    "A": 89.094,  "R": 174.203, "N": 132.119, "D": 133.104,
    "C": 121.154, "Q": 146.146, "E": 147.130, "G": 75.067,
    "H": 155.156, "I": 131.175, "L": 131.175, "K": 146.189,
    "M": 149.208, "F": 165.192, "P": 115.132, "S": 105.093,
    "T": 119.119, "W": 204.228, "Y": 181.191, "V": 117.148,
}
WATER_MW: Final[float] = 18.015  # Perda de agua por ligacao peptidica

# pKa dos grupos ionizaveis para calculo do ponto isoeletrico
# Valores de Lehninger (classicos, usados pelo ExPASy ProtParam)
PK_N_TERM: Final[float] = 9.69
PK_C_TERM: Final[float] = 2.34
PK_SIDE: Final[dict[str, float]] = {
    "D": 3.65,   # Asp
    "E": 4.25,   # Glu
    "C": 8.18,   # Cys
    "Y": 10.07,  # Tyr
    "H": 6.00,   # His
    "K": 10.53,  # Lys
    "R": 12.48,  # Arg
}

# Indice de instabilidade DIWV (Guruprasad et al. 1990)
# Tabela 10x10 de pesos dipeptidicos
DIWV: Final[dict[str, float]] = {
    "WW": 1.0, "WC": 1.0, "WM": 24.68, "WH": 24.68, "WY": 1.0,
    "WF": 1.0, "WQ": 1.0, "WN": 13.34, "WI": 1.0, "WR": 1.0,
    "WD": 1.0, "WP": 1.0, "WT": -14.03, "WK": 1.0, "WE": 1.0,
    "WV": -7.49, "WS": 1.0, "WG": -9.37, "WA": -14.03, "WL": 13.34,
    "CW": 24.68, "CC": 1.0, "CM": 33.6, "CH": 33.6, "CY": 1.0,
    "CF": 1.0, "CQ": -6.54, "CN": 1.0, "CI": 1.0, "CR": 1.0,
    "CD": 20.26, "CP": 20.26, "CT": 33.6, "CK": 1.0, "CE": 1.0,
    "CV": -6.54, "CS": 1.0, "CG": 1.0, "CA": 1.0, "CL": 20.26,
    "MW": 1.0, "MC": 1.0, "MM": -1.88, "MH": 58.28, "MY": 24.68,
    "MF": 1.0, "MQ": -6.54, "MN": 1.0, "MI": 1.0, "MR": -6.54,
    "MD": 1.0, "MP": 44.94, "MT": -1.88, "MK": 1.0, "ME": 1.0,
    "MV": 1.0, "MS": 44.94, "MG": 1.0, "MA": 13.34, "ML": 1.0,
    "HW": -1.88, "HC": 1.0, "HM": 1.0, "HH": 1.0, "HY": 44.94,
    "HF": -9.37, "HQ": 1.0, "HN": 24.68, "HI": 44.94, "HR": -1.88,
    "HD": 1.0, "HP": -1.88, "HT": -6.54, "HK": 24.68, "HE": 1.0,
    "HV": 1.0, "HS": 1.0, "HG": -9.37, "HA": 1.0, "HL": 1.0,
    "YW": -9.37, "YC": 1.0, "YM": 44.94, "YH": 13.34, "YY": 13.34,
    "YF": 1.0, "YQ": 1.0, "YN": 1.0, "YI": 1.0, "YR": -15.91,
    "YD": 24.68, "YP": 13.34, "YT": -7.49, "YK": 1.0, "YE": -6.54,
    "YV": 1.0, "YS": 1.0, "YG": -7.49, "YA": 24.68, "YL": 1.0,
    "FW": 1.0, "FC": 1.0, "FM": 1.0, "FH": 1.0, "FY": 33.6,
    "FF": 1.0, "FQ": 1.0, "FN": 1.0, "FI": 1.0, "FR": 1.0,
    "FD": 13.34, "FP": 20.26, "FT": 1.0, "FK": -14.03, "FE": 1.0,
    "FV": 1.0, "FS": 1.0, "FG": 1.0, "FA": 1.0, "FL": 1.0,
    "QW": 1.0, "QC": -6.54, "QM": 1.0, "QH": 1.0, "QY": -6.54,
    "QF": -6.54, "QQ": 20.26, "QN": 1.0, "QI": 1.0, "QR": 1.0,
    "QD": 20.26, "QP": 20.26, "QT": 1.0, "QK": 1.0, "QE": 20.26,
    "QV": -6.54, "QS": 44.94, "QG": 1.0, "QA": 1.0, "QL": 1.0,
    "NW": -9.37, "NC": -1.88, "NM": 1.0, "NH": 1.0, "NY": 1.0,
    "NF": -14.03, "NQ": -6.54, "NN": 1.0, "NI": 44.94, "NR": 1.0,
    "ND": 1.0, "NP": -1.88, "NT": -7.49, "NK": 24.68, "NE": 1.0,
    "NV": 1.0, "NS": 1.0, "NG": -14.03, "NA": 1.0, "NL": 1.0,
    "IW": 1.0, "IC": 1.0, "IM": 1.0, "IH": 13.34, "IY": 1.0,
    "IF": 1.0, "IQ": 1.0, "IN": 1.0, "II": 1.0, "IR": 1.0,
    "ID": 1.0, "IP": -1.88, "IT": 1.0, "IK": -7.49, "IE": 44.94,
    "IV": -7.49, "IS": 1.0, "IG": 1.0, "IA": 1.0, "IL": 20.26,
    "RW": 58.28, "RC": 1.0, "RM": 1.0, "RH": 20.26, "RY": -6.54,
    "RF": 1.0, "RQ": 20.26, "RN": 13.34, "RI": 1.0, "RR": 58.28,
    "RD": 1.0, "RP": 20.26, "RT": 1.0, "RK": 1.0, "RE": 1.0,
    "RV": 1.0, "RS": 44.94, "RG": -7.49, "RA": 1.0, "RL": 1.0,
    "DW": 1.0, "DC": 1.0, "DM": 1.0, "DH": 1.0, "DY": 1.0,
    "DF": -6.54, "DQ": 1.0, "DN": 1.0, "DI": 1.0, "DR": -6.54,
    "DD": 1.0, "DP": 1.0, "DT": -14.03, "DK": -7.49, "DE": 1.0,
    "DV": 1.0, "DS": 20.26, "DG": 1.0, "DA": 1.0, "DL": 1.0,
    "PW": -1.88, "PC": -6.54, "PM": -6.54, "PH": 1.0, "PY": 1.0,
    "PF": 20.26, "PQ": 20.26, "PN": 1.0, "PI": 1.0, "PR": -6.54,
    "PD": -6.54, "PP": 20.26, "PT": 1.0, "PK": 1.0, "PE": 18.38,
    "PV": 20.26, "PS": 20.26, "PG": 1.0, "PA": 20.26, "PL": 1.0,
    "TW": -14.03, "TC": 1.0, "TM": 1.0, "TH": 1.0, "TY": 1.0,
    "TF": 13.34, "TQ": -6.54, "TN": -14.03, "TI": 1.0, "TR": 1.0,
    "TD": 1.0, "TP": 1.0, "TT": 1.0, "TK": 1.0, "TE": 20.26,
    "TV": 1.0, "TS": 1.0, "TG": -7.49, "TA": 1.0, "TL": 1.0,
    "KW": 1.0, "KC": 1.0, "KM": 33.6, "KH": 1.0, "KY": 1.0,
    "KF": 1.0, "KQ": 24.68, "KN": 1.0, "KI": -7.49, "KR": 33.6,
    "KD": 1.0, "KP": -6.54, "KT": 1.0, "KK": 1.0, "KE": 1.0,
    "KV": -7.49, "KS": 1.0, "KG": -7.49, "KA": 1.0, "KL": -7.49,
    "EW": -14.03, "EC": 44.94, "EM": 1.0, "EH": -6.54, "EY": 1.0,
    "EF": 1.0, "EQ": 20.26, "EN": 1.0, "EI": 20.26, "ER": 1.0,
    "ED": 20.26, "EP": 20.26, "ET": 1.0, "EK": 1.0, "EE": 33.6,
    "EV": 1.0, "ES": 20.26, "EG": 1.0, "EA": 11.0, "EL": 1.0,
    "VW": 1.0, "VC": 1.0, "VM": 1.0, "VH": 1.0, "VY": -6.54,
    "VF": 1.0, "VQ": 1.0, "VN": 1.0, "VI": 1.0, "VR": 1.0,
    "VD": -14.03, "VP": 20.26, "VT": -7.49, "VK": -1.88, "VE": 1.0,
    "VV": 1.0, "VS": 1.0, "VG": -7.49, "VA": 1.0, "VL": 1.0,
    "SW": 1.0, "SC": 1.0, "SM": 1.0, "SH": 1.0, "SY": 1.0,
    "SF": 1.0, "SQ": 20.26, "SN": 1.0, "SI": 1.0, "SR": 20.26,
    "SD": 1.0, "SP": 44.94, "ST": 1.0, "SK": 1.0, "SE": 20.26,
    "SV": 1.0, "SS": 20.26, "SG": 1.0, "SA": 1.0, "SL": 1.0,
    "GW": 13.34, "GC": 1.0, "GM": 1.0, "GH": 1.0, "GY": -7.49,
    "GF": 1.0, "GQ": 1.0, "GN": -7.49, "GI": -7.49, "GR": 1.0,
    "GD": 1.0, "GP": 1.0, "GT": -7.49, "GK": -7.49, "GE": -6.54,
    "GV": 1.0, "GS": 1.0, "GG": 13.34, "GA": -7.49, "GL": 1.0,
    "AW": 1.0, "AC": 44.94, "AM": 1.0, "AH": -7.49, "AY": 1.0,
    "AF": 1.0, "AQ": 1.0, "AN": 1.0, "AI": 1.0, "AR": 1.0,
    "AD": -7.49, "AP": 20.26, "AT": 1.0, "AK": 1.0, "AE": 1.0,
    "AV": 1.0, "AS": 1.0, "AG": 1.0, "AA": 1.0, "AL": 1.0,
    "LW": 24.68, "LC": 1.0, "LM": 1.0, "LH": 1.0, "LY": 1.0,
    "LF": 1.0, "LQ": 33.6, "LN": 1.0, "LI": 1.0, "LR": 20.26,
    "LD": 1.0, "LP": 20.26, "LT": 1.0, "LK": -7.49, "LE": 1.0,
    "LV": 1.0, "LS": 1.0, "LG": 1.0, "LA": 1.0, "LL": 1.0,
}


# ---------------------------------------------------------------------------
# Dataclasses de resultado
# ---------------------------------------------------------------------------


@dataclass
class EcoliConstruct:
    """Constructo redesenhado para expressao em E. coli."""

    protein_sequence: str
    length_aa: int
    molecular_weight_da: float
    isoelectric_point: float
    instability_index: float
    is_stable: bool
    has_internal_ndei: bool
    has_internal_xhoi: bool
    components: dict[str, str] = field(default_factory=dict)

    # Detalhamento das regioes do constructo
    his_tag: str = ""
    tev_site: str = ""
    adjuvant: str = ""
    epitope_cassette: str = ""
    unique_epitopes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Funcoes de calculo fisico-quimico
# ---------------------------------------------------------------------------


def calculate_molecular_weight(sequence: str) -> float:
    """Calcula o peso molecular de uma proteina pela soma dos residuos.

    Usa os pesos moleculares medios dos aminoacidos e subtrai a agua
    liberada em cada ligacao peptidica (n-1 moleculas de agua para n
    residuos).

    Args:
        sequence: Sequencia proteica em letra unica maiuscula.

    Returns:
        Peso molecular em Daltons.
    """
    total = sum(AA_MW.get(aa, 0.0) for aa in sequence)
    # Cada ligacao peptidica libera uma molecula de agua
    water_loss = (len(sequence) - 1) * WATER_MW
    return round(total - water_loss, 2)


def calculate_isoelectric_point(sequence: str) -> float:
    """Calcula o ponto isoeletrico (pI) pelo metodo de biseccao.

    O pI e o pH no qual a carga liquida da proteina e zero.
    Utiliza os valores de pKa de Lehninger, consistente com o ExPASy
    ProtParam.

    Args:
        sequence: Sequencia proteica em letra unica maiuscula.

    Returns:
        pI estimado com precisao de 0.01 unidades de pH.
    """
    # Contar residuos ionizaveis
    aa_count: dict[str, int] = {}
    for aa in sequence:
        if aa in PK_SIDE:
            aa_count[aa] = aa_count.get(aa, 0) + 1

    def _charge_at_ph(ph: float) -> float:
        """Carga liquida no pH dado usando a equacao de Henderson-Hasselbalch."""
        # N-terminal positivo
        charge = 1.0 / (1.0 + 10 ** (ph - PK_N_TERM))
        # C-terminal negativo
        charge -= 1.0 / (1.0 + 10 ** (PK_C_TERM - ph))
        # Cadeias laterais
        for aa, pka in PK_SIDE.items():
            n = aa_count.get(aa, 0)
            if n == 0:
                continue
            if aa in ("D", "E", "C", "Y"):
                # Acidos -- carregados negativamente quando desprotonados
                charge -= n / (1.0 + 10 ** (pka - ph))
            else:
                # Basicos (H, K, R) -- carregados positivamente quando protonados
                charge += n / (1.0 + 10 ** (ph - pka))
        return charge

    # Biseccao entre pH 0 e pH 14
    low, high = 0.0, 14.0
    while (high - low) > 0.001:
        mid = (low + high) / 2.0
        if _charge_at_ph(mid) > 0:
            low = mid
        else:
            high = mid

    return round((low + high) / 2.0, 2)


def calculate_instability_index(sequence: str) -> float:
    """Calcula o indice de instabilidade de Guruprasad (1990).

    Proteinas com II < 40 sao classificadas como estaveis in vivo.
    O indice e calculado a partir dos pesos estatisticos de dipeptideos
    (tabela DIWV) observados em proteinas com meia-vida conhecida.

    Args:
        sequence: Sequencia proteica em letra unica maiuscula.

    Returns:
        Indice de instabilidade (adimensional).
    """
    if len(sequence) < 2:
        return 0.0

    total = 0.0
    for i in range(len(sequence) - 1):
        dipeptide = sequence[i] + sequence[i + 1]
        total += DIWV.get(dipeptide, 1.0)

    return round((10.0 / len(sequence)) * total, 2)


# ---------------------------------------------------------------------------
# Montagem do constructo
# ---------------------------------------------------------------------------


def extract_unique_epitopes(epitopes_raw: list[dict]) -> list[str]:
    """Extrai epitopos unicos na ordem de aparicao original.

    O construct_card.json pode ter epitopos repetidos (ex: LLTANVCYK aparece
    5 vezes porque foi encontrado em 5 genes CPB paralogos). Para o constructo
    precisamos de cada peptideo unico uma unica vez.

    Args:
        epitopes_raw: Lista de dicts de epitopos do construct_card.json.

    Returns:
        Lista de sequencias peptidicas unicas preservando a ordem.
    """
    seen: set[str] = set()
    unique: list[str] = []
    for ep in epitopes_raw:
        peptide = ep["peptide"]
        if peptide not in seen:
            seen.add(peptide)
            unique.append(peptide)
    return unique


def build_ecoli_construct(
    epitopes_raw: list[dict],
    adjuvant: str = L7L12_ADJUVANT,
    linker_adj: str = LINKER_ADJUVANT,
    linker_ctl: str = LINKER_CTL,
) -> EcoliConstruct:
    """Monta o constructo proteico para expressao em E. coli.

    Arquitetura:
        Met - His6 - TEV - L7/L12 - EAAAK - [Epi1-AAY-Epi2-AAY-...-EpiN]

    O Met inicial ja esta embutido no sitio NdeI (CATATG) do vetor.
    A His6-tag permite purificacao por IMAC (Ni-NTA).
    O sitio TEV (ENLYFQS) permite remocao da tag, deixando Ser no N-terminal.

    Args:
        epitopes_raw: Lista de dicts de epitopos do construct_card.json.
        adjuvant: Sequencia do adjuvante (padrao: L7/L12).
        linker_adj: Linker entre adjuvante e cassete de epitopos.
        linker_ctl: Linker entre epitopos CTL consecutivos.

    Returns:
        EcoliConstruct com todos os campos calculados.
    """
    unique_epitopes = extract_unique_epitopes(epitopes_raw)

    # Montar cassete de epitopos com linkers AAY entre eles
    epitope_cassette = linker_ctl.join(unique_epitopes)

    # Sequencia proteica completa:
    # M + His6 + TEV + adjuvante + EAAAK + cassete de epitopos
    protein_parts = [
        "M",           # Met do inicio (codificado pelo ATG do NdeI)
        HIS6_TAG,      # HHHHHH -- tag de afinidade
        TEV_SITE,      # ENLYFQS -- sitio de clivagem TEV
        adjuvant,      # L7/L12 completo
        linker_adj,    # EAAAK -- linker flexivel
        epitope_cassette,
    ]
    protein_sequence = "".join(protein_parts)

    # Calcular propriedades fisico-quimicas
    mw = calculate_molecular_weight(protein_sequence)
    pi = calculate_isoelectric_point(protein_sequence)
    ii = calculate_instability_index(protein_sequence)

    # Mapa de componentes para rastreabilidade
    components = {
        "met_start": "M",
        "his6_tag": HIS6_TAG,
        "tev_site": TEV_SITE,
        "adjuvant_l7l12": adjuvant,
        "linker_adjuvant": linker_adj,
        "epitope_cassette": epitope_cassette,
        "linker_ctl": linker_ctl,
    }

    return EcoliConstruct(
        protein_sequence=protein_sequence,
        length_aa=len(protein_sequence),
        molecular_weight_da=mw,
        isoelectric_point=pi,
        instability_index=ii,
        is_stable=ii < 40.0,
        has_internal_ndei=False,  # Sera verificado apos otimizacao de codons
        has_internal_xhoi=False,  # Sera verificado apos otimizacao de codons
        components=components,
        his_tag=HIS6_TAG,
        tev_site=TEV_SITE,
        adjuvant=adjuvant,
        epitope_cassette=epitope_cassette,
        unique_epitopes=unique_epitopes,
    )


def validate_construct(construct: EcoliConstruct) -> list[str]:
    """Executa validacoes sobre o constructo montado.

    Verifica:
    - Inicia com Met
    - Contem His6 e sitio TEV
    - Indice de instabilidade < 40 (estavel)
    - pI nao esta proximo de 7 (precipitacao no armazenamento)
    - Peso molecular dentro do esperado para expressao em E. coli

    Args:
        construct: O constructo montado.

    Returns:
        Lista de avisos/erros encontrados (vazia se tudo ok).
    """
    warnings: list[str] = []

    if not construct.protein_sequence.startswith("M"):
        warnings.append("ERRO: Sequencia nao inicia com Met")

    if HIS6_TAG not in construct.protein_sequence:
        warnings.append("ERRO: His6-tag nao encontrada na sequencia")

    if TEV_SITE not in construct.protein_sequence:
        warnings.append("ERRO: Sitio TEV nao encontrado na sequencia")

    if not construct.is_stable:
        warnings.append(
            f"AVISO: Indice de instabilidade = {construct.instability_index:.2f} "
            f"(>40 = instavel)"
        )

    if 6.5 <= construct.isoelectric_point <= 7.5:
        warnings.append(
            f"AVISO: pI = {construct.isoelectric_point:.2f} esta proximo de "
            f"pH neutro — risco de precipitacao em tampao fisiologico"
        )

    # Proteinas >100 kDa podem ter problemas de expressao em E. coli
    if construct.molecular_weight_da > 100_000:
        warnings.append(
            f"AVISO: MW = {construct.molecular_weight_da:.0f} Da (>100 kDa) — "
            f"pode ter expressao reduzida em E. coli"
        )

    return warnings


def to_dict(construct: EcoliConstruct) -> dict:
    """Serializa o constructo para formato JSON compativel com o pipeline.

    Args:
        construct: O constructo montado.

    Returns:
        Dicionario com todos os campos relevantes.
    """
    return {
        "protein_sequence": construct.protein_sequence,
        "length_aa": construct.length_aa,
        "molecular_weight_da": construct.molecular_weight_da,
        "isoelectric_point": construct.isoelectric_point,
        "instability_index": construct.instability_index,
        "is_stable": construct.is_stable,
        "has_internal_ndei": construct.has_internal_ndei,
        "has_internal_xhoi": construct.has_internal_xhoi,
        "unique_epitope_count": len(construct.unique_epitopes),
        "unique_epitopes": construct.unique_epitopes,
        "components": construct.components,
    }
