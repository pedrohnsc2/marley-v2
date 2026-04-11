"""Configuracao do modulo 04_rna_fm — Analise de RNA com deep learning.

Define sequencias biologicas de referencia, controles negativos e
hiperparametros do encoder customizado de RNA. As sequencias SL RNA
de multiplas especies de Leishmania permitem analise de conservacao
cross-species, enquanto snRNAs/rRNAs/tRNAs humanos e caninos servem
como controles negativos para prova de unicidade.

Refs:
    - Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
    - Michaeli S (2011) Parasitology 138(12):1473-1487
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from marley_ai.config import AIModuleConfig, AI_ROOT


# ---------------------------------------------------------------------------
# Sequencias SL RNA de especies de Leishmania
# Conservado ha ~500 milhoes de anos em trypanosomatideos.
# O miniexon (SL RNA) e adicionado a todo mRNA via trans-splicing.
# Ref: Michaeli S (2011) Parasitology 138(12):1473-1487
# ---------------------------------------------------------------------------

SL_SEQUENCES: Final[dict[str, str]] = {
    # L. infantum — alvo principal do projeto Marley
    "L_infantum": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    # L. major — modelo experimental classico, genoma de referencia
    "L_major": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG",
    # L. donovani — agente da leishmaniose visceral (kala-azar)
    "L_donovani": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    # L. braziliensis — agente da leishmaniose mucocutanea
    "L_braziliensis": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG",
}

# ---------------------------------------------------------------------------
# Sequencias de controle negativo — snRNAs humanos
# Small nuclear RNAs envolvidos no spliceosome de mamiferos.
# Mecanismo de splicing fundamentalmente diferente (cis-splicing).
# Ref: Will CL & Luhrmann R (2011) Cold Spring Harb Perspect Biol
# ---------------------------------------------------------------------------

HUMAN_SNRNA: Final[dict[str, str]] = {
    # U1 snRNA — reconhece sitio 5' splice via pareamento de bases
    "U1_snRNA": (
        "AUACUUACCUGGCAGGGGAGAUACCAUGAUCACGAAGGUG"
        "GUUUUCCCAGGGCGAGGCUUAUCCAUUGCACUCCGGAUGU"
        "GAUGACUUCCAAUU"
    ),
    # U2 snRNA — reconhece branch point, essencial para splicing
    "U2_snRNA": (
        "AUCCUUUUGCUUUGGCUAAAGAUCAAGUGUAGUAUCUGUUC"
        "CUUUAUAUAUAUUAAAUGGAUUUUUGGGAUCUUUAGCAAAG"
    ),
    # U4 snRNA — forma di-snRNP com U6
    "U4_snRNA": (
        "AGCUUUGCGCAGUGGCAGUGUGCAAUAUGCUUCUUGGAGCAA"
        "GGCUCAUAUGCUAGAUUGGCAGAAUUCGGCCGGGUUGCG"
    ),
    # U5 snRNA — interage com ambos exons no spliceosome
    "U5_snRNA": (
        "AUACUAUUUCCAAUUUUUUGUCCAAUCCUUGGUACAAAGAAU"
        "UUUUCGCAAGCCUCCUUUUCUACUGGAUCGGAUUUC"
    ),
    # U6 snRNA — atividade catalitica no splicing (ribozima)
    "U6_snRNA": (
        "GUGCUCCCCGCGCGACUGCUAAUGAAAGCGUUACGGCACG"
        "GUAUACUGCCCAGGGAAAACCGCUGAAGGCUUUUCACGGCA"
    ),
}

# ---------------------------------------------------------------------------
# Controles negativos adicionais — RNAs caninos (hospedeiro)
# tRNAs e rRNAs de Canis lupus familiaris.
# Ref: Lindblad-Toh K et al. (2005) Nature 438(7069):803-819
# ---------------------------------------------------------------------------

CANINE_CONTROL_RNA: Final[dict[str, str]] = {
    # tRNA-Ala (anticodon AGC) — molecula adaptadora da traducao
    "tRNA_Ala_canine": (
        "GGGGAAUUAGCUCAAGUGGUAGAGCGCUUGCUUAGCAUGC"
        "AAGAGGUCCUGGGUUCAAUCCCCAGAUUCUCCACC"
    ),
    # tRNA-Gly (anticodon GCC) — alta expressao constitutiva
    "tRNA_Gly_canine": (
        "GCAUUGGUGGUUCAGUGGUAGAAUUCUCGCCUGCCACGCGGG"
        "AGGCCCGGGUUCAAAUUCCCGGCCAAUGCACC"
    ),
    # 5.8S rRNA — componente do ribossomo 60S
    "rRNA_5_8S_canine": (
        "GACUACUUGGAUCGAUGAUCGAUGAAAUGCGAUAAUUAGUGUG"
        "AAUUGCAGAAUUCAGUGAAUCAUCGAAUCUCUUGAACGCAA"
    ),
    # 5S rRNA — componente do ribossomo 60S, gene multigene
    "rRNA_5S_canine": (
        "GCCUACGGCCAUACCAGCCUGAAAGCACCGCAUUCCCCGUCCG"
        "AUCUGACUCUUCAGAUCCUUCCCUCCACCCCACUCU"
    ),
}

# ---------------------------------------------------------------------------
# MRL-ASO-001 — melhor candidato do pipeline Marley
# Ref: aso_math/config.py (rank 1, score 0.7647)
# ---------------------------------------------------------------------------

ASO_SEQUENCE: Final[str] = "ACAGAAACTGATACTTATATAGCGT"
ASO_TARGET_START: Final[int] = 5   # posicao 0-indexed no SL RNA
ASO_TARGET_END: Final[int] = 30

# ---------------------------------------------------------------------------
# Hiperparametros do encoder customizado de RNA
# Dimensionados para treinar rapidamente em Apple M3 Pro (36GB RAM).
# Arquitetura: Transformer encoder com predicao de nucleotideo mascarado
# (inspirado em BERT/RNA-FM, mas treinavel from scratch no dataset pequeno).
# ---------------------------------------------------------------------------

ENCODER_PARAMS: Final[dict[str, int | float]] = {
    "vocab_size": 5,           # A, U, G, C + [MASK]
    "max_seq_len": 256,        # comprimento maximo de sequencia
    "embed_dim": 128,          # dimensao dos embeddings
    "num_heads": 4,            # atencao multi-cabeca
    "num_layers": 3,           # profundidade do transformer
    "feedforward_dim": 256,    # camada feedforward interna
    "dropout": 0.1,            # regularizacao
    "mask_fraction": 0.15,     # fracao de posicoes mascaradas (BERT-style)
    "learning_rate": 1e-3,     # taxa de aprendizado (Adam)
    "num_epochs": 100,         # epocas de treinamento
    "batch_size": 16,          # tamanho do mini-batch
}

# ---------------------------------------------------------------------------
# Regras de pareamento de bases para algoritmo de Nussinov
# Ref: Nussinov R et al. (1978) SIAM J Appl Math 35(1):68-82
# ---------------------------------------------------------------------------

BASE_PAIR_RULES: Final[dict[tuple[str, str], float]] = {
    ("A", "U"): 1.0,   # Watson-Crick: 2 pontes de hidrogenio
    ("U", "A"): 1.0,
    ("G", "C"): 1.0,   # Watson-Crick: 3 pontes de hidrogenio (mais forte)
    ("C", "G"): 1.0,
    ("G", "U"): 0.5,   # Wobble: 2 pontes de hidrogenio (mais fraco)
    ("U", "G"): 0.5,
}

# Distancia minima no loop para pareamento (evita hairpins muito curtos)
MIN_LOOP_LENGTH: Final[int] = 4


# ---------------------------------------------------------------------------
# Configuracao frozen dataclass — segue padrao AIModuleConfig
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RNAFMConfig(AIModuleConfig):
    """Configuracao especifica para o modulo de analise de RNA.

    Combina encoder customizado de RNA (transformer treinavel) com
    predicao de estrutura secundaria (Nussinov) e analise de unicidade
    do SL RNA de Leishmania.
    """

    # --- Modelo ---
    rna_model: str = "custom_rna_encoder"
    embedding_dim: int = ENCODER_PARAMS["embed_dim"]

    # --- Sequencia alvo ---
    sl_sequence: str = SL_SEQUENCES["L_infantum"]
    predict_structure: bool = True
    predict_accessibility: bool = True

    # --- Caminhos ---
    rna_embeddings_dir: Path = field(
        default_factory=lambda: AI_ROOT / "04_rna_fm" / "embeddings"
    )
    model_checkpoint_dir: Path = field(
        default_factory=lambda: AI_ROOT / "04_rna_fm" / "checkpoints"
    )
