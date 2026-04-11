"""C2 -- Redesenho do construto vacinal para expressao em L. tarentolae.

Diferentemente da plataforma E. coli (B), L. tarentolae e um eucarioto
que realiza trans-splicing e adicao de poli-A, possui maquinaria de
glicosilacao e ancoras GPI, e pode secretar proteinas pela via
secretoria convencional.

Modificacoes em relacao ao construto mRNA (Plataforma A):
    - REMOVER: peptideo sinal tPA (sinal mamifero, nao reconhecido)
    - ADICIONAR: SAP1 de L. mexicana (sinal de secrecao para Leishmania)
    - ADICIONAR: 5'UTR e 3'UTR de L. tarentolae (para trans-splicing correto)
    - ADICIONAR: sequencias flanqueadoras do SSU rRNA (integracao cromossomal)
    - ADICIONAR: His6-Strep double tag (purificacao por afinidade tandem)
    - MANTER: adjuvante L7/L12, todos os 11 epitopos, linkers AAY/EAAAK
    - VETOR: pLEXSY-sat2 (resistencia a nourseotricina)

Referencia: Breitling et al. (2002) Protein Expr Purif 25(2):209-218
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

from vaccine_platforms.shared.epitopes import (
    ADJUVANT_EPITOPE_LINKER,
    ADJUVANT_SEQUENCE,
    EPITOPES,
    LINKERS,
)
from vaccine_platforms.shared.validators import (
    calculate_molecular_weight,
    calculate_isoelectric_point,
    calculate_instability_index,
    calculate_gravy,
    calculate_aromaticity,
)

# ---------------------------------------------------------------------------
# Constantes biologicas especificas de L. tarentolae / LEXSY
# ---------------------------------------------------------------------------

# Peptideo sinal SAP1 de L. mexicana (Secreted Acid Phosphatase 1)
# Funciona como sinal de secrecao no RE de trypanosomatideos.
# Ref: Carrion et al. (2009) Mol Biochem Parasitol 168(2):196-202
SAP1_SIGNAL: Final[str] = "MKFFVFALFAVALCSAEA"

# Adjuvante L7/L12 (mesmo das plataformas A e B)
L7L12_ADJUVANT: Final[str] = ADJUVANT_SEQUENCE

# Linkers
LINKER_ADJUVANT: Final[str] = ADJUVANT_EPITOPE_LINKER  # EAAAK
LINKER_CTL: Final[str] = LINKERS["CTL"]                # AAY

# Tags de purificacao
HIS6_TAG: Final[str] = "HHHHHH"
STREP_TAG: Final[str] = "WSHPQFEK"  # Strep-tag II (IBA/Sigma)

# Linker entre His6 e Strep (flexivel, sem estrutura secundaria)
TAG_LINKER: Final[str] = "GSGSGS"

# Peptideo sinal tPA -- sera REMOVIDO (referencia)
TPA_SIGNAL: Final[str] = "MDAMKRGLCCVLLLCGAVFVSAS"

# Sequencia do 5'UTR de L. tarentolae para trans-splicing correto.
# O mini-exon do SL (39 nt) e adicionado por trans-splicing ao 5'UTR
# do mRNA. No vetor pLEXSY, a regiao intergênica fornece o sinal
# de aceitacao para o splicing.
# Ref: Manual pLEXSY, Jena Bioscience (2023)
UTR5_TARENTOLAE: Final[str] = (
    "CTGGATCCCCAAGCTTGGGTCCCCCTGCAATTTGAGCTGTTTTTGCAGAGGGCCAGCTCA"
    "GCGGCCG"
)

# 3'UTR de L. tarentolae (regiao de poliadenilacao)
# Contem o sinal de poliadenilacao para processamento correto do mRNA
# em trypanosomatideos.
UTR3_TARENTOLAE: Final[str] = (
    "TGAGTTTTTGTGTTAGAAAGAGCAGGTTTTAATGTATCGTTCGCAGATCAGAA"
    "CACTTTCCAATAG"
)

# Sequencias flanqueadoras do sitio de integracao SSU rRNA.
# A integracao cromossomal no locus do 18S rRNA e o metodo padrao
# do sistema LEXSY. O construto e ladeado por ~1 kb de homologia
# ao SSU rRNA para recombinacao homologa.
# Aqui representamos apenas os primeiros 60 nt de cada flanco.
SSU_5_FLANK: Final[str] = (
    "GATCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGT"
)
SSU_3_FLANK: Final[str] = (
    "CCTGCGGAAGGATCATTACTGGCCGGAAATCCTTTGGGTTCCTATAGGTCGAAACGACTT"
)

# Vetor e selecao
VECTOR_NAME: Final[str] = "pLEXSY-sat2"
SELECTION_MARKER: Final[str] = "nourseothricin"
SELECTION_CONCENTRATION: Final[str] = "100 ug/mL"

# Justificativa da plataforma L. tarentolae
PLATFORM_RATIONALE: Final[dict[str, str]] = {
    "biosafety": "BSL-1 (nao patogenico para mamiferos; infecta lagartixas)",
    "growth_conditions": "26C, meio BHI suplementado, sem CO2, tempo de duplicacao ~8h",
    "ptm_capability": (
        "N-glicosilacao, ancoras GPI, acetilacao N-terminal -- "
        "equivalentes a L. infantum patogenica"
    ),
    "intrinsic_adjuvant": (
        "Moleculas do parasita (LPG, GP63, GIPLs) estimulam TLR2/TLR4 "
        "do sistema imune inato, eliminando necessidade de adjuvante exogeno "
        "para a forma viva"
    ),
    "literature_precedent": (
        "HPV (Breitling 2002), HCV (Niimi 2018), SARS-CoV-2 (Pirdel 2021) -- "
        "proteinas vacinais expressas com sucesso em L. tarentolae"
    ),
    "commercial_system": (
        "LEXSY (Jena Bioscience): kit comercial maduro com vetores, cepas "
        "e protocolos otimizados para expressao constitutiva e induzivel"
    ),
    "vector": VECTOR_NAME,
    "selection": f"{SELECTION_MARKER} ({SELECTION_CONCENTRATION})",
}


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass
class TarentolaeConstruct:
    """Construto redesenhado para expressao em L. tarentolae."""

    # Sequencia proteica completa (SAP1 + adjuvante + epitopos + tags)
    protein_sequence: str
    length_aa: int
    molecular_weight_da: float
    isoelectric_point: float
    instability_index: float
    is_stable: bool
    gravy: float
    aromaticity: float

    # Componentes individuais
    signal_peptide: str
    adjuvant: str
    epitope_cassette: str
    his6_tag: str
    strep_tag: str
    tag_linker: str
    unique_epitopes: list[str] = field(default_factory=list)

    # Sequencias flanqueadoras (DNA, nao proteina)
    utr5: str = ""
    utr3: str = ""
    ssu_5_flank: str = ""
    ssu_3_flank: str = ""

    # Metadados
    vector: str = VECTOR_NAME
    selection: str = SELECTION_MARKER
    components: dict[str, str] = field(default_factory=dict)
    rationale: dict[str, str] = field(default_factory=dict)
    modifications_from_mrna: list[dict[str, str]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Montagem do construto
# ---------------------------------------------------------------------------


def extract_unique_epitopes_from_ssot() -> list[str]:
    """Extrai epitopos unicos do SSOT (shared/epitopes.py).

    Usa a fonte unica de verdade em vez do construct_card.json
    para garantir consistencia entre plataformas.

    Returns:
        Lista de sequencias peptidicas dos 11 epitopos unicos.
    """
    return [ep.peptide for ep in EPITOPES]


def build_tarentolae_construct() -> TarentolaeConstruct:
    """Monta o construto proteico para expressao em L. tarentolae.

    Arquitetura (N -> C-terminal):
        SAP1 - L7/L12 - EAAAK - [Epi1-AAY-Epi2-...-EpiN] - GSGSGS - His6 - Strep

    O peptideo sinal SAP1 e clivado no RE, resultando na proteina
    madura secretada com a estrutura:
        L7/L12 - EAAAK - [cassete de epitopos] - GSGSGS - His6 - Strep

    A dupla tag His6-Strep permite purificacao tandem de alta pureza:
        1. IMAC (Ni-NTA): captura grosseira pela His6
        2. Strep-Tactin: polimento pela Strep-tag
        3. Pureza esperada: >95%

    Returns:
        TarentolaeConstruct com todos os campos calculados.
    """
    unique_epitopes = extract_unique_epitopes_from_ssot()

    # Montar cassete de epitopos com linkers AAY
    epitope_cassette = LINKER_CTL.join(unique_epitopes)

    # Sequencia proteica completa (inclui peptideo sinal)
    protein_parts = [
        SAP1_SIGNAL,       # Peptideo sinal SAP1 de L. mexicana
        L7L12_ADJUVANT,    # Adjuvante L7/L12 ribosomal
        LINKER_ADJUVANT,   # EAAAK -- linker rigido
        epitope_cassette,  # 11 epitopos separados por AAY
        TAG_LINKER,        # GSGSGS -- linker flexivel
        HIS6_TAG,          # HHHHHH -- tag de afinidade metal
        STREP_TAG,         # WSHPQFEK -- Strep-tag II
    ]
    protein_sequence = "".join(protein_parts)

    # Calcular propriedades fisico-quimicas
    mw = calculate_molecular_weight(protein_sequence)
    pi = calculate_isoelectric_point(protein_sequence)
    ii = calculate_instability_index(protein_sequence)
    gravy = calculate_gravy(protein_sequence)
    arom = calculate_aromaticity(protein_sequence)

    # Mapa de componentes para rastreabilidade
    components = {
        "signal_peptide_sap1": SAP1_SIGNAL,
        "adjuvant_l7l12": L7L12_ADJUVANT,
        "linker_adjuvant": LINKER_ADJUVANT,
        "epitope_cassette": epitope_cassette,
        "linker_ctl": LINKER_CTL,
        "tag_linker": TAG_LINKER,
        "his6_tag": HIS6_TAG,
        "strep_tag": STREP_TAG,
    }

    # Documentar modificacoes em relacao ao construto mRNA
    modifications = [
        {
            "component": "signal_peptide",
            "mrna_platform": f"tPA ({TPA_SIGNAL})",
            "tarentolae_platform": f"SAP1 ({SAP1_SIGNAL})",
            "reason": (
                "tPA e sinal de secrecao mamifero, nao reconhecido pela "
                "maquinaria do RE de trypanosomatideos. SAP1 de L. mexicana "
                "e sinal nativo de Leishmania."
            ),
        },
        {
            "component": "purification_tag",
            "mrna_platform": "None (mRNA traduzido in vivo)",
            "tarentolae_platform": f"His6-GSGSGS-Strep ({HIS6_TAG}{TAG_LINKER}{STREP_TAG})",
            "reason": (
                "Proteina recombinante secretada precisa de purificacao. "
                "Tag dupla permite purificacao tandem com pureza >95%."
            ),
        },
        {
            "component": "utrs",
            "mrna_platform": "5'UTR humano + 3'UTR com poli(A)",
            "tarentolae_platform": "Regioes intergenicas de L. tarentolae para trans-splicing",
            "reason": (
                "mRNAs de trypanosomatideos requerem trans-splicing do SL RNA "
                "no 5' e poliadenilacao no 3'. UTRs de mamiferos nao funcionam."
            ),
        },
        {
            "component": "integration",
            "mrna_platform": "None (mRNA transiente)",
            "tarentolae_platform": "SSU rRNA locus (recombinacao homologa)",
            "reason": (
                "Integracao cromossomal no locus 18S rRNA garante expressao "
                "estavel e constitutiva sem necessidade de pressao seletiva "
                "continua (embora nourseotricina seja usada para selecao inicial)."
            ),
        },
        {
            "component": "codon_usage",
            "mrna_platform": "Canis lupus familiaris (otimizado para traducao in vivo)",
            "tarentolae_platform": "L. tarentolae (GC-rich, trypanosomatid bias)",
            "reason": (
                "Trypanosomatideos tem forte vies de codons GC-rich (GC ~60%). "
                "Codons otimizados para mamiferos tem expressao reduzida."
            ),
        },
        {
            "component": "vector",
            "mrna_platform": "None (mRNA sintetico)",
            "tarentolae_platform": f"{VECTOR_NAME} ({SELECTION_MARKER})",
            "reason": (
                "pLEXSY-sat2 contem as regioes flanqueadoras do SSU rRNA, "
                "o marcador de selecao sat2 (nourseotricina acetiltransferase), "
                "e o promotor constitutivo para expressao em L. tarentolae."
            ),
        },
    ]

    # Validacoes
    warnings: list[str] = []

    if not protein_sequence.startswith("M"):
        warnings.append("ERRO: Sequencia nao inicia com Met")

    if ii >= 40.0:
        warnings.append(
            f"AVISO: Indice de instabilidade = {ii:.2f} (>40 = instavel)"
        )

    if 6.5 <= pi <= 7.5:
        warnings.append(
            f"AVISO: pI = {pi:.2f} proximo de pH neutro — risco de precipitacao"
        )

    if mw > 80_000:
        warnings.append(
            f"AVISO: MW = {mw:.0f} Da (>80 kDa) — pode ter expressao reduzida "
            f"em L. tarentolae; considerar fragmentacao do construto"
        )

    return TarentolaeConstruct(
        protein_sequence=protein_sequence,
        length_aa=len(protein_sequence),
        molecular_weight_da=mw,
        isoelectric_point=pi,
        instability_index=ii,
        is_stable=ii < 40.0,
        gravy=gravy,
        aromaticity=arom,
        signal_peptide=SAP1_SIGNAL,
        adjuvant=L7L12_ADJUVANT,
        epitope_cassette=epitope_cassette,
        his6_tag=HIS6_TAG,
        strep_tag=STREP_TAG,
        tag_linker=TAG_LINKER,
        unique_epitopes=unique_epitopes,
        utr5=UTR5_TARENTOLAE,
        utr3=UTR3_TARENTOLAE,
        ssu_5_flank=SSU_5_FLANK,
        ssu_3_flank=SSU_3_FLANK,
        vector=VECTOR_NAME,
        selection=SELECTION_MARKER,
        components=components,
        rationale=PLATFORM_RATIONALE,
        modifications_from_mrna=modifications,
        warnings=warnings,
    )


def validate_construct(construct: TarentolaeConstruct) -> list[str]:
    """Executa validacoes adicionais sobre o construto montado.

    Verifica:
    - Peptideo sinal SAP1 esta presente
    - Adjuvante L7/L12 intacto
    - Todos os 11 epitopos presentes
    - Tags de purificacao no C-terminal
    - Propriedades fisico-quimicas dentro dos limites

    Args:
        construct: O construto montado.

    Returns:
        Lista de avisos/erros encontrados.
    """
    warnings: list[str] = list(construct.warnings)

    if SAP1_SIGNAL not in construct.protein_sequence:
        warnings.append("ERRO: Peptideo sinal SAP1 nao encontrado")

    if L7L12_ADJUVANT not in construct.protein_sequence:
        warnings.append("ERRO: Adjuvante L7/L12 nao encontrado intacto")

    if HIS6_TAG not in construct.protein_sequence:
        warnings.append("ERRO: His6-tag nao encontrada")

    if STREP_TAG not in construct.protein_sequence:
        warnings.append("ERRO: Strep-tag nao encontrada")

    # Verificar todos os 11 epitopos
    expected_epitopes = [ep.peptide for ep in EPITOPES]
    for ep_seq in expected_epitopes:
        if ep_seq not in construct.protein_sequence:
            warnings.append(f"ERRO: Epitopo {ep_seq} nao encontrado no construto")

    # Verificar que nao contem tPA
    if TPA_SIGNAL in construct.protein_sequence:
        warnings.append(
            "ERRO: Peptideo sinal tPA encontrado -- deve ser substituido por SAP1"
        )

    return warnings


def to_dict(construct: TarentolaeConstruct) -> dict:
    """Serializa o construto para formato JSON.

    Args:
        construct: O construto montado.

    Returns:
        Dicionario serializavel para JSON.
    """
    return {
        "protein_sequence": construct.protein_sequence,
        "length_aa": construct.length_aa,
        "molecular_weight_da": construct.molecular_weight_da,
        "isoelectric_point": construct.isoelectric_point,
        "instability_index": construct.instability_index,
        "is_stable": construct.is_stable,
        "gravy": construct.gravy,
        "aromaticity": construct.aromaticity,
        "unique_epitope_count": len(construct.unique_epitopes),
        "unique_epitopes": construct.unique_epitopes,
        "components": construct.components,
        "dna_elements": {
            "utr5_length_nt": len(construct.utr5),
            "utr3_length_nt": len(construct.utr3),
            "ssu_5_flank_length_nt": len(construct.ssu_5_flank),
            "ssu_3_flank_length_nt": len(construct.ssu_3_flank),
        },
        "vector": construct.vector,
        "selection": construct.selection,
        "rationale": construct.rationale,
        "modifications_from_mrna": construct.modifications_from_mrna,
        "architecture": (
            f"[SAP1 signal ({len(construct.signal_peptide)} aa)] - "
            f"[L7/L12 adjuvant ({len(construct.adjuvant)} aa)] - "
            f"[EAAAK] - "
            f"[{len(construct.unique_epitopes)} epitopes + AAY linkers "
            f"({len(construct.epitope_cassette)} aa)] - "
            f"[{construct.tag_linker}] - [His6] - [Strep-tag II]"
        ),
    }
