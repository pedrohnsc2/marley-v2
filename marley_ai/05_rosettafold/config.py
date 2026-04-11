"""Configuracao do modulo 05_rosettafold — Analise estrutural do complexo ASO:SL RNA.

Define parametros geometricos para acidos nucleicos: forma-A (RNA:DNA hibrido),
modificacoes LNA (C3'-endo constrainado), backbone fosforotioato (PS), e
constantes termodinamicas para decomposicao de energia de ligacao.

Ref: Saenger W (1984) Principles of Nucleic Acid Structure. Springer.
     Egli M et al. (2005) Chem Biol 12(6):669-675 (LNA geometry).
     Eckstein F (2014) NAR 42(22):13762-13774 (PS backbone).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from marley_ai.config import AIModuleConfig, AI_ROOT, PROJECT_ROOT


# ---------------------------------------------------------------------------
# Sequencias biologicas (espelhadas de aso_math/config.py para independencia)
# ---------------------------------------------------------------------------

# SL RNA de L. infantum — alvo do ASO (39 nt, DNA template)
SL_SEQUENCE: Final[str] = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG"
SL_LENGTH: Final[int] = len(SL_SEQUENCE)

# SL RNA como RNA (U no lugar de T) — forma biologica in vivo
SL_RNA_SEQUENCE: Final[str] = SL_SEQUENCE.replace("T", "U")

# MRL-ASO-001: melhor candidato do pipeline (rank 1)
ASO_SEQUENCE: Final[str] = "ACAGAAACTGATACTTATATAGCGT"
ASO_LENGTH: Final[int] = len(ASO_SEQUENCE)

# Regiao-alvo no SL RNA (posicoes 5-30, 0-indexed)
ASO_TARGET_START: Final[int] = 5
ASO_TARGET_END: Final[int] = 30
TARGET_SEQUENCE: Final[str] = SL_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END]

# Sequencia complementar do alvo no SL RNA (fita sense do RNA)
TARGET_RNA_SEQUENCE: Final[str] = TARGET_SEQUENCE.replace("T", "U")

# ---------------------------------------------------------------------------
# Mapeamento de complementaridade
# ---------------------------------------------------------------------------

# Complemento DNA padrao (inclui U para input RNA)
COMPLEMENT: Final[dict[str, str]] = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A",
}

# Complemento RNA (retorna ribonucleotideo)
RNA_COMPLEMENT: Final[dict[str, str]] = {
    "A": "U", "U": "A", "G": "C", "C": "G", "T": "A",
}

# ---------------------------------------------------------------------------
# Parametros geometricos — forma-A (RNA:DNA hibrido)
# Ref: Saenger (1984), Tabela 9.2; Egli & Manoharan (2023) Acc Chem Res
# ---------------------------------------------------------------------------

# Ascensao por par de bases ao longo do eixo helicoidal (Angstroms)
A_FORM_RISE: Final[float] = 2.81

# Angulo de rotacao por par de bases (graus)
A_FORM_TWIST: Final[float] = 32.7

# Inclinacao das bases em relacao ao eixo (graus) — forma-A tem inclinacao
A_FORM_INCLINATION: Final[float] = 19.0

# Deslocamento-x das bases em relacao ao eixo helicoidal (Angstroms)
A_FORM_X_DISPLACEMENT: Final[float] = -4.4

# Numero de pares de bases por volta completa
A_FORM_BP_PER_TURN: Final[float] = 360.0 / A_FORM_TWIST  # ~11.0

# Diametro do duplex (Angstroms)
A_FORM_DIAMETER: Final[float] = 23.0

# Largura dos sulcos (Angstroms)
A_FORM_MAJOR_GROOVE_WIDTH: Final[float] = 2.7   # sulco maior — estreito na forma-A
A_FORM_MINOR_GROOVE_WIDTH: Final[float] = 11.0  # sulco menor — largo na forma-A

# Profundidade dos sulcos (Angstroms)
A_FORM_MAJOR_GROOVE_DEPTH: Final[float] = 13.5
A_FORM_MINOR_GROOVE_DEPTH: Final[float] = 2.8

# ---------------------------------------------------------------------------
# Parametros geometricos — forma-B (DNA:DNA, referencia de comparacao)
# ---------------------------------------------------------------------------

B_FORM_RISE: Final[float] = 3.38
B_FORM_TWIST: Final[float] = 36.0
B_FORM_INCLINATION: Final[float] = -1.2
B_FORM_X_DISPLACEMENT: Final[float] = 0.7
B_FORM_BP_PER_TURN: Final[float] = 360.0 / B_FORM_TWIST  # 10.0

# ---------------------------------------------------------------------------
# Parametros de modificacao LNA
# Ref: Egli M et al. (2005) Chem Biol 12:669-675
#      Koshkin AA et al. (1998) Tetrahedron 54:3607-3630
# ---------------------------------------------------------------------------

# Posicoes LNA no MRL-ASO-001 (convencionalmente, 3 na ponta 5' e 3 na 3')
# Posicoes 0-indexed no ASO de 24 nt
LNA_POSITIONS_5PRIME: Final[list[int]] = [0, 1, 2]
LNA_POSITIONS_3PRIME: Final[list[int]] = [21, 22, 23]
LNA_POSITIONS: Final[list[int]] = LNA_POSITIONS_5PRIME + LNA_POSITIONS_3PRIME

# Conformacao do acucar LNA — constrainado C3'-endo (Norte)
# Pseudorotation phase angle P (graus) — C3'-endo ideal
LNA_SUGAR_PUCKER_PHASE: Final[float] = 18.0  # graus — tipo Norte puro
LNA_SUGAR_PUCKER_AMPLITUDE: Final[float] = 38.0  # graus — amplitude maxima

# DNA normal — mistura C2'-endo (Sul) / C3'-endo
DNA_SUGAR_PUCKER_PHASE: Final[float] = 162.0  # graus — tipo Sul (C2'-endo)

# RNA — C3'-endo (Norte), similar ao LNA
RNA_SUGAR_PUCKER_PHASE: Final[float] = 18.0  # graus

# Efeito do LNA na Tm: +3 a +5 C por modificacao LNA (vs DNA)
LNA_TM_INCREMENT_PER_MOD: Final[float] = 4.0  # graus Celsius

# LNA forca geometria A-form no vizinhanca — raio de influencia (pb)
LNA_INFLUENCE_RADIUS: Final[int] = 2  # pares de bases vizinhos afetados

# ---------------------------------------------------------------------------
# Parametros de backbone fosforotioato (PS)
# Ref: Eckstein F (2014) NAR 42:13762-13774
# ---------------------------------------------------------------------------

# Backbone PS: S substitui um O nao-ponte no fosfato
# Angulo de torsao alpha ligeiramente diferente vs PO (fosfodiester)
PS_BOND_ANGLE_SHIFT: Final[float] = 3.5  # graus — desvio do PO ideal

# Comprimento de ligacao P-S vs P-O (Angstroms)
PS_BOND_LENGTH: Final[float] = 2.05   # P-S
PO_BOND_LENGTH: Final[float] = 1.62   # P-O (referencia)

# Raio de van der Waals do S vs O (Angstroms)
VDW_SULFUR: Final[float] = 1.80
VDW_OXYGEN: Final[float] = 1.52

# Efeito na Tm: PS reduz ~0.5 C por modificacao vs PO
PS_TM_DECREMENT_PER_MOD: Final[float] = -0.5  # graus Celsius

# ---------------------------------------------------------------------------
# Parametros de atomo para construcao PDB
# ---------------------------------------------------------------------------

# Raios atomicos para calculo de SASA (Angstroms)
ATOMIC_RADII: Final[dict[str, float]] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "P": 1.80,
    "S": 1.80,
    "H": 1.20,
}

# Raio da sonda de solvente para SASA (Angstroms)
SOLVENT_PROBE_RADIUS: Final[float] = 1.4

# ---------------------------------------------------------------------------
# Parametros termodinamicos para decomposicao de energia
# ---------------------------------------------------------------------------

# Parametros nearest-neighbor RNA:DNA hibrido (Sugimoto et al. 1995)
# Ref: Sugimoto N et al. (1995) Biochemistry 34:11211-11216
# Chave: rXdY — r=RNA, d=DNA, XY lidos na direcao 5'->3' da fita RNA
# Valores: (dH kcal/mol, dS cal/(mol*K))
NN_RNA_DNA_HYBRID: Final[dict[str, tuple[float, float]]] = {
    "rAA/dTT": (-7.8, -21.9),
    "rAC/dTG": (-5.9, -12.3),
    "rAG/dTC": (-9.1, -23.5),
    "rAU/dTA": (-8.3, -23.9),
    "rCA/dGT": (-9.0, -26.1),
    "rCC/dGG": (-9.3, -23.2),
    "rCG/dGC": (-16.3, -47.1),
    "rCU/dGA": (-7.0, -19.7),
    "rGA/dCT": (-5.5, -13.5),
    "rGC/dCG": (-8.0, -17.1),
    "rGG/dCC": (-12.8, -31.9),
    "rGU/dCA": (-7.8, -21.6),
    "rUA/dAT": (-7.8, -23.2),
    "rUC/dAG": (-8.6, -22.9),
    "rUG/dAC": (-10.4, -28.4),
    "rUU/dAA": (-11.5, -36.4),
}

# Iniciacao de duplex RNA:DNA (Sugimoto 1995)
NN_HYBRID_INIT: Final[tuple[float, float]] = (1.9, -3.9)

# Energia de ligacao de hidrogenio por tipo de par (kcal/mol)
# Ref: Sponer J et al. (2001) J Phys Chem A 105:10759
HBOND_ENERGY: Final[dict[str, float]] = {
    "GC": -6.0,   # 3 pontes H (G-C Watson-Crick)
    "CG": -6.0,
    "AT": -4.0,   # 2 pontes H (A-T Watson-Crick)
    "TA": -4.0,
    "AU": -4.0,   # 2 pontes H (A-U Watson-Crick, equivalente RNA)
    "UA": -4.0,
    "GU": -3.5,   # wobble pair (1 H-bond forte + 1 fraca)
    "UG": -3.5,
    "GT": -3.5,
    "TG": -3.5,
}

# Energia de stacking por dinucleotideo (kcal/mol, simplificado)
# Ref: Turner DH & Mathews DH (2010) NAR 38:D209-D215
STACKING_ENERGY: Final[dict[str, float]] = {
    "AA": -1.0, "AC": -1.5, "AG": -1.3, "AT": -0.9, "AU": -0.9,
    "CA": -1.6, "CC": -2.1, "CG": -2.0, "CT": -1.4, "CU": -1.4,
    "GA": -1.3, "GC": -2.3, "GG": -1.8, "GT": -1.2, "GU": -1.2,
    "TA": -0.6, "TC": -1.0, "TG": -0.8, "TT": -0.5,
    "UA": -0.6, "UC": -1.0, "UG": -0.8, "UU": -0.5, "UT": -0.5,
}

# Constantes fisicas
R_GAS: Final[float] = 1.987              # cal/(mol*K)
T_PHYSIOLOGICAL: Final[float] = 310.15   # 37 C em Kelvin
KBOLTZMANN: Final[float] = 0.001987      # kcal/(mol*K)

# Constante dieletrica para calculos eletrostaticos
DIELECTRIC_WATER: Final[float] = 80.0
DIELECTRIC_PROTEIN: Final[float] = 4.0

# Carga do fosfato por residuo (unidades de carga elementar)
PHOSPHATE_CHARGE: Final[float] = -1.0

# Penalidade de dessolvatacao por par de bases enterrado (kcal/mol)
DESOLVATION_PENALTY_PER_BP: Final[float] = 0.2

# ---------------------------------------------------------------------------
# Caminhos
# ---------------------------------------------------------------------------

# Estruturas existentes
STRUCTURES_DIR: Final[Path] = PROJECT_ROOT / "data" / "structures"
SL_RNA_PDB: Final[Path] = PROJECT_ROOT / "data" / "sl_rna_3d" / "sl_rna.pdb"
ASO_PDB: Final[Path] = STRUCTURES_DIR / "MRL_ASO_001.pdb"

# Saida do modulo
MODULE_DIR: Final[Path] = AI_ROOT / "05_rosettafold"
OUTPUT_STRUCTURES_DIR: Final[Path] = MODULE_DIR / "structures"
OUTPUT_RESULTS_DIR: Final[Path] = MODULE_DIR / "results"


# ---------------------------------------------------------------------------
# Dataclass de configuracao do modulo
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StructuralAnalysisConfig(AIModuleConfig):
    """Configuracao para analise estrutural do complexo ASO:SL RNA.

    Substitui RoseTTAFoldConfig — nao depende de RoseTTAFold2NA,
    implementa analise geometrica computacional direta.
    """

    # --- Identidade ---
    module_slug: str = "05_rosettafold"
    module_name: str = "Structural Analysis - ASO:SL RNA Complex"

    # --- Geometria ---
    helix_form: str = "A"           # A-form para RNA:DNA hibrido
    rise_per_bp: float = A_FORM_RISE
    twist_per_bp: float = A_FORM_TWIST
    inclination: float = A_FORM_INCLINATION
    x_displacement: float = A_FORM_X_DISPLACEMENT

    # --- Modificacoes ---
    lna_positions: tuple[int, ...] = tuple(LNA_POSITIONS)
    use_ps_backbone: bool = True    # fosforotioato no gap central

    # --- Caminhos de entrada ---
    sl_rna_pdb: Path = field(default_factory=lambda: SL_RNA_PDB)
    aso_pdb: Path = field(default_factory=lambda: ASO_PDB)
    structures_dir: Path = field(default_factory=lambda: STRUCTURES_DIR)

    # --- Caminhos de saida ---
    output_structures: Path = field(default_factory=lambda: OUTPUT_STRUCTURES_DIR)
    output_results: Path = field(default_factory=lambda: OUTPUT_RESULTS_DIR)
