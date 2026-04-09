"""Configuracao centralizada do suite de validacao matematica aso_math.

Todas as constantes biologicas, termodinamicas e de caminho sao definidas
aqui para evitar duplicacao e garantir consistencia entre os 5 modulos.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

# ---------------------------------------------------------------------------
# Caminhos do projeto
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
ASO_MATH_ROOT: Final[Path] = Path(__file__).resolve().parent
RESULTS_DIR: Final[Path] = ASO_MATH_ROOT / "results"
CANDIDATES_JSON: Final[Path] = PROJECT_ROOT / "results" / "aso" / "aso_candidates.json"
ASO_PDB: Final[Path] = PROJECT_ROOT / "data" / "structures" / "MRL_ASO_001.pdb"
SL_PDB: Final[Path] = PROJECT_ROOT / "data" / "sl_rna_3d" / "sl_rna.pdb"
TRANSCRIPTOME_HUMAN: Final[Path] = PROJECT_ROOT / "data" / "raw" / "transcriptome_human.fasta"
TRANSCRIPTOME_LINF: Final[Path] = PROJECT_ROOT / "data" / "raw" / "transcriptome_linf.fasta"

# ---------------------------------------------------------------------------
# Sequencias biologicas
# ---------------------------------------------------------------------------

# Spliced Leader RNA de L. infantum (39 nt, GC = 28.2%)
# Adicionado a todo mRNA do parasita via trans-splicing.
# Conservado ha ~500 milhoes de anos em trypanosomatideos.
# Ausente em mamiferos — alvo seletivo por construcao.
SL_SEQUENCE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"
SL_LENGTH: Final[int] = len(SL_SEQUENCE)

# MRL-ASO-001: melhor candidato do pipeline (rank 1, score 0.7647)
ASO_SEQUENCE: Final[str] = "ACAGAAACTGATACTTATATAGCGT"
ASO_LENGTH: Final[int] = len(ASO_SEQUENCE)
ASO_TARGET_START: Final[int] = 5   # posicao 0-indexed no SL RNA
ASO_TARGET_END: Final[int] = 30
ASO_TARGET_SEQUENCE: Final[str] = SL_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END]

# Valores conhecidos do pipeline original (usados como validacao)
ASO_KNOWN_TM: Final[float] = 68.48       # Celsius
ASO_KNOWN_DG: Final[float] = -27.97      # kcal/mol
ASO_KNOWN_GC: Final[float] = 0.32

# ---------------------------------------------------------------------------
# Parametros termodinamicos — SantaLucia 1998
# Referencia: SantaLucia J Jr. (1998) PNAS 95(4):1460-1465
# DNA/DNA duplex parameters.
# Chave: dinucleotideo XY/X'Y' (X'Y' = complemento)
# Valores: (delta_H em kcal/mol, delta_S em cal/(mol*K))
# ---------------------------------------------------------------------------

NN_PARAMS: Final[dict[str, tuple[float, float]]] = {
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}

# Parametros de iniciacao (adicionados uma vez por duplex)
NN_INIT: Final[tuple[float, float]] = (0.1, -2.8)

# Mapeamento de complemento DNA (inclui U para input RNA)
COMPLEMENT: Final[dict[str, str]] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A",
}

# Bases possiveis
BASES: Final[str] = "ACGT"

# ---------------------------------------------------------------------------
# Constantes termodinamicas
# ---------------------------------------------------------------------------

R_GAS: Final[float] = 1.987               # cal / (mol * K)
ASO_CONCENTRATION: Final[float] = 250e-9  # 250 nM padrao
T_PHYSIOLOGICAL: Final[float] = 310.15    # 37 C em Kelvin

# ---------------------------------------------------------------------------
# Parametros de Leishmania
# ---------------------------------------------------------------------------

# Taxa de mutacao por nucleotideo por replicacao
# Ref: Rogers MB et al. (2011) PLoS Genetics 7(8):e1002237
# MA-line estimate: ~1-5 x 10^-9 por base por divisao mitotica
MUTATION_RATE: Final[float] = 2.0e-9

# Tempo de geracao (divisao de amastigota em macrofago)
GENERATION_TIME_HOURS: Final[float] = 12.0

# Numero de copias do SL RNA em tandem array
# Ref: Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
SL_RNA_COPY_NUMBER: Final[int] = 150  # ~100-200 copias

# Limiar de dG acima do qual o ASO nao funciona mais
# Ref: Crooke et al. (2017) Nucleic Acids Res — ASOs precisam de ~14-16 bp
# para ativacao de RNase H, correspondendo a dG ~ -15 a -20 kcal/mol
DG_FUNCTIONAL_THRESHOLD: Final[float] = -15.0  # kcal/mol

# ---------------------------------------------------------------------------
# Faixa de comprimento para varredura
# ---------------------------------------------------------------------------

LENGTH_SCAN_MIN: Final[int] = 18
LENGTH_SCAN_MAX: Final[int] = 27

# ---------------------------------------------------------------------------
# Identificadores de modulos (usados no JSON envelope)
# ---------------------------------------------------------------------------

MODULES: Final[dict[str, str]] = {
    "01": "01_thermodynamic_landscape",
    "02": "02_selectivity_proof",
    "03": "03_evolutionary_conservation",
    "04": "04_exhaustive_optimization",
    "05": "05_resistance_model",
}
