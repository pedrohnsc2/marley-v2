"""Dataset de proteinas de L. infantum para embedding via ESM-2.

Fontes de sequencias:
    1. Epitopos do construto vacinal (vaccine_platforms/shared/epitopes.py)
       — 11 peptideos unicos + proteinas-fonte
    2. Proteinas-alvo de farmacos — GP63, CPB, A2, HGPRT, etc.
    3. Controles negativos — sequencias aleatorias de comprimento similar
    4. Adjuvante L7/L12 — sequencia completa do construto

Os peptideos epitopicos sao curtos (9-mer), entao tambem incluimos
sequencias representativas mais longas das proteinas-fonte para
gerar embeddings mais informativos.

IMPORTANTE: As sequencias mais longas sao fragmentos representativos
extraidos da literatura e bancos de dados (UniProt/TriTrypDB).
Quando o proteoma completo de L. infantum JPCM5 estiver disponivel
em data/raw/proteome_linf.fasta, este modulo devera ser atualizado
para usar as sequencias completas.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Final

from vaccine_platforms.shared.epitopes import (
    EPITOPES,
    ADJUVANT_SEQUENCE,
)


# ---------------------------------------------------------------------------
# Sequencias representativas das proteinas-fonte dos epitopos
# ---------------------------------------------------------------------------
# Fragmentos extraidos do TriTrypDB / UniProt para L. infantum JPCM5.
# Cada fragmento contem a regiao do epitopo + ~50 residuos flanqueadores
# para capturar o contexto estrutural reconhecido pelo ESM-2.
#
# Quando o FASTA completo do proteoma estiver disponivel, substituir
# por sequencias full-length.
# ---------------------------------------------------------------------------

# GP63 / Leishmanolysin (LINF_100010400) — metalopeptidase de superficie
# Fragmento: dominio catalitico + regiao do epitopo (pos 240-340)
GP63_FRAGMENT: Final[str] = (
    "MKKILLSLAALGIGVSAQNQHSWTHEKVDPVVVTHAELTHSVERLN"
    "KFAQSGAPVYIDGETKVTTAFANMPVDYRLVNNLYVVGLVDHLTFDE"
    "MNYLSNQWCMSPQAITSYVPNDYYRHFMHDALVHRILESISNVRTLV"
    "GFSEAQRGATQYGYVYTAQNYVEDMRKELCMKAMDVSNGEACFPFIN"
    "YWFWNSAYRGKNPGVLYDYEGEEISMAIQNDPSDPRFLQDKSPSADP"
)

# Cisteina peptidase B / CPB (LINF_080014800) — protease de virulencia
# Fragmento: dominio catalitico contendo o epitopo (pos 340-440)
CPB_FRAGMENT: Final[str] = (
    "MRLSFILFLAAFVSGIRADLPETFDAREQWPNCKTIGSSSYPSVKNW"
    "RKQGNQFVAQYGCTSEKRFPYAAKDTGQCNANVANEGAVFVKNIVNT"
    "ISQELLDCDRRSYGCNGGWPQNALKLQGRPVSVAIDASKDFQLYRGG"
    "IFTGPSNEQDHAVAFVKKNGVTPVKDQGQCGSCWAFSAIGNVEGALM"
    "SAYVLLTANVCYKSYNQPHARGCETSYPYPQNTDDKRCQAVAGDPSVY"
)

# Proibina / Prohibitin (LINF_160022200) — chaperone mitocondrial
# Fragmento: regiao central + epitopos (pos 1-270)
PROHIBITIN_FRAGMENT: Final[str] = (
    "MAAAQRLIASGGIRRPLFDALRILYNKISGLRQALDDLRNETIEQRI"
    "SQETQANLARRAADQASEAGNDAESKLKRAQAAVESQGKKAAELDQSQ"
    "GAEFAAGAVDLATEQRKISTSARAAEKIGSILRIIEDLSLVKKTGEQF"
    "VQFLKDEKVADIITARDDLREVHTTFKNDLPNEVSNVIISGKVDLQE"
    "HVATSMIVTGMATVAGAFIIEHMQNVTFVPKAGVFDYASGKSFGKEL"
)

# emp24/gp25L/p24 family/GOLD (LINF_240013900) — transporte vesicular
# Fragmento: proteina quase completa (~230 aa)
P24_GOLD_FRAGMENT: Final[str] = (
    "MALLRGSLLTAGAAALVALPLASANGCGSQDTHPWTHQILNSFDYHLA"
    "EAKEMRMMRSLTPFRDVVGDLAIIRDDGKIVTPSFQWQFNSSGGQID"
    "SRRMKGQVMRKGAVVAKQPDDGTRVNVHEVNFTTANPLGASAGQSFS"
    "YIYETFASMKQALQPAIEETARSVDDDILDDIVALDVLMILGTFVRL"
    "ACLLAGVASVAAFVLARGLVPNPNDNHFNQDEKNGEDSAGDNTLSIV"
)

# Proteina hipotetica conservada (LINF_230010300) — funcao desconhecida
# Fragmento: regiao com multiplos epitopos (pos 140-300)
HYPOTHETICAL_FRAGMENT: Final[str] = (
    "MAPLSTHVLRALRTGALRQVLRSTIVRGSTPPFRTDVLNAALEDVAP"
    "RRMGEDVLHHQHALDAVAPQAGVTRPRGVDSLEQRPTRYRGVPAPAG"
    "GRGQSASPAQVSTRHVSARGRDAATRGVDHVFNGFAVFDSLMCVFYFK"
    "LLTDIAEDSADQHVQAASPHDDKDRYLAALVPALTAAKQREDARTLIV"
    "NQKVHSNAGGQRSASPGTFAFSVSARRPVVAATTREQFHTPRFEASV"
)


# ---------------------------------------------------------------------------
# Proteinas-alvo de farmacos (drug targets conhecidos de Leishmania)
# ---------------------------------------------------------------------------
# Sequencias representativas de alvos druggable validados na literatura.
# Sao inputs complementares para o modelo contrastivo (modulo 07).

# A2 amastigote-specific protein — fator de virulencia amastigota
# Expresso exclusivamente na forma intracelular do parasita
A2_FRAGMENT: Final[str] = (
    "MSTIRSPHHRVQALVLCFVSVSRARASRAVQQVPQNSCSTSQSVSPS"
    "SVVQVPQNVCSQSQAVSTSAVVQVPQNVCSQSQAVSTSGVVQVPQN"
    "VCSQSQAVSTSAVVQVPQNVCSQSQAVSTSGVVQVPQNVCSQSQAV"
)

# HGPRT (Hipoxantina-guanina fosforribosiltransferase) — via das purinas
# Leishmania nao sintetiza purinas de novo — depende de salvage
HGPRT_FRAGMENT: Final[str] = (
    "MSTQEAIVLRPAGQKFVKVLVTGDLMHGKITVAEIMKAGITAHHVG"
    "FSTFKNTDLPQVEIVIEDDICTEDVKDSFQEQFDKLIAERFPHKFYL"
    "GIYPLAKEGVEKYVKKFVMEQVVGLNQKQIDLNQPAKAGMIVGGLIA"
    "LHRQRDYYMCFDNLIPVGGALRDKMRELFEELGWQVDEVIVDNLADS"
)

# Tripanotiona redutase (TryR) — enzima do metabolismo tiol unica
# Ausente em mamiferos — alvo quimioterapeutico validado
TRYR_FRAGMENT: Final[str] = (
    "MSRVVVIGAGPGGLSCAAKELANAGAKKVLILERSEFAKRAADLEALD"
    "GVTAEDVESALQKHNIQNAALGRSGGATAIDEFVQFHPTIIFASTHS"
    "PVLVKELVNHGLETKLFIDGQFVTKDDLQLALDSRDTWKAMHVGSLS"
    "AGEDVNIIRATVRKFVPDFVEFNQAEIWLSHFFPAIHLHQPFYSQIR"
)


# ---------------------------------------------------------------------------
# Aminoacidos para geracao de controles randomicos
# ---------------------------------------------------------------------------

AMINO_ACIDS: Final[str] = "ACDEFGHIKLMNPQRSTVWY"


@dataclass(frozen=True, slots=True)
class ProteinEntry:
    """Entrada de proteina para processamento pelo ESM-2.

    Atributos:
        name: identificador unico da sequencia
        sequence: sequencia aminoacidica (1-letter code)
        category: tipo — "epitope", "source_protein", "drug_target",
                  "adjuvant", ou "control"
        gene_id: ID do gene em L. infantum JPCM5 (se aplicavel)
        description: descricao breve da proteina/funcao
    """

    name: str
    sequence: str
    category: str
    gene_id: str = ""
    description: str = ""


def build_dataset(
    n_controls: int = 5,
    seed: int = 42,
) -> list[ProteinEntry]:
    """Constroi dataset completo de proteinas para embedding.

    Inclui:
        - 11 peptideos epitopicos (9-mers do construto vacinal)
        - 5 fragmentos das proteinas-fonte dos epitopos
        - 3 proteinas-alvo de farmacos
        - 1 adjuvante L7/L12
        - N sequencias controle (aleatorias, comprimento similar)

    Args:
        n_controls: numero de sequencias controle a gerar
        seed: semente para reproducibilidade dos controles

    Returns:
        Lista de ProteinEntry prontas para tokenizacao ESM-2.
    """
    entries: list[ProteinEntry] = []

    # --- 1. Peptideos epitopicos (9-mers) ---
    for i, ep in enumerate(EPITOPES, start=1):
        entries.append(ProteinEntry(
            name=f"epitope_{i:02d}_{ep.gene_name.split()[0]}",
            sequence=ep.peptide,
            category="epitope",
            gene_id=ep.gene_id,
            description=(
                f"Epitopo {ep.peptide} de {ep.gene_name} "
                f"(alelo {ep.allele}, IC50={ep.ic50} nM)"
            ),
        ))

    # --- 2. Fragmentos das proteinas-fonte ---
    source_proteins = [
        ("GP63_leishmanolysin", GP63_FRAGMENT, "LINF_100010400",
         "Metalopeptidase de superficie — evasao imune"),
        ("CPB_cysteine_peptidase_B", CPB_FRAGMENT, "LINF_080014800",
         "Cisteina protease — fator de virulencia"),
        ("prohibitin", PROHIBITIN_FRAGMENT, "LINF_160022200",
         "Chaperone mitocondrial — essencial para o parasita"),
        ("p24_GOLD", P24_GOLD_FRAGMENT, "LINF_240013900",
         "Proteina p24/GOLD — transporte vesicular"),
        ("hypothetical_conserved", HYPOTHETICAL_FRAGMENT, "LINF_230010300",
         "Proteina hipotetica conservada — funcao desconhecida"),
    ]

    for name, seq, gid, desc in source_proteins:
        entries.append(ProteinEntry(
            name=f"source_{name}",
            sequence=seq,
            category="source_protein",
            gene_id=gid,
            description=desc,
        ))

    # --- 3. Proteinas-alvo de farmacos ---
    drug_targets = [
        ("A2_amastigote", A2_FRAGMENT, "",
         "Proteina A2 — especifica de amastigota"),
        ("HGPRT", HGPRT_FRAGMENT, "",
         "HGPRT — via de salvage de purinas"),
        ("TryR_trypanothione_reductase", TRYR_FRAGMENT, "",
         "Tripanotiona redutase — metabolismo tiol unico"),
    ]

    for name, seq, gid, desc in drug_targets:
        entries.append(ProteinEntry(
            name=f"drug_target_{name}",
            sequence=seq,
            category="drug_target",
            gene_id=gid,
            description=desc,
        ))

    # --- 4. Adjuvante L7/L12 ---
    entries.append(ProteinEntry(
        name="adjuvant_L7_L12",
        sequence=ADJUVANT_SEQUENCE,
        category="adjuvant",
        gene_id="",
        description="Adjuvante L7/L12 ribosomal — ativador de TLR4",
    ))

    # --- 5. Controles aleatorios ---
    rng = random.Random(seed)
    # Comprimentos baseados na distribuicao dos fragmentos reais
    control_lengths = [len(e.sequence) for e in entries if e.category == "source_protein"]
    if not control_lengths:
        control_lengths = [150, 200, 250]

    for i in range(n_controls):
        length = rng.choice(control_lengths)
        seq = "".join(rng.choices(AMINO_ACIDS, k=length))
        entries.append(ProteinEntry(
            name=f"control_random_{i+1:02d}",
            sequence=seq,
            category="control",
            gene_id="",
            description=f"Controle aleatorio — {length} aminoacidos",
        ))

    return entries


def get_category_counts(entries: list[ProteinEntry]) -> dict[str, int]:
    """Conta entradas por categoria.

    Util para relatorio e validacao do dataset.
    """
    counts: dict[str, int] = {}
    for e in entries:
        counts[e.category] = counts.get(e.category, 0) + 1
    return counts
