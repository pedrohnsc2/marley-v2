"""Decomposicao de energia de ligacao do complexo ASO:SL RNA.

Calcula a energia livre de ligacao (deltaG) decomposta em componentes:
- Base stacking (nearest-neighbor)
- Pontes de hidrogenio (Watson-Crick + wobble)
- Interacoes eletrostaticas (backbone fosforotioato)
- Penalidade de dessolvatacao
- Contribuicao LNA (estabilizacao termica adicional)

O objetivo e validar computacionalmente que MRL-ASO-001 tem
energia de ligacao suficiente para ativacao de RNase H (< -15 kcal/mol).

Ref: SantaLucia J Jr (1998) PNAS 95:1460-1465.
     Sugimoto N et al. (1995) Biochemistry 34:11211-11216.
     Kurreck J et al. (2002) NAR 30:1911-1918 (LNA-ASO).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import importlib as _il

import numpy as np

# Importacao via importlib — "05" nao e valido como identificador Python
_cfg = _il.import_module("marley_ai.05_rosettafold.config")

ASO_SEQUENCE = _cfg.ASO_SEQUENCE
ASO_TARGET_END = _cfg.ASO_TARGET_END
ASO_TARGET_START = _cfg.ASO_TARGET_START
COMPLEMENT = _cfg.COMPLEMENT
DESOLVATION_PENALTY_PER_BP = _cfg.DESOLVATION_PENALTY_PER_BP
DIELECTRIC_WATER = _cfg.DIELECTRIC_WATER
HBOND_ENERGY = _cfg.HBOND_ENERGY
KBOLTZMANN = _cfg.KBOLTZMANN
LNA_POSITIONS = _cfg.LNA_POSITIONS
LNA_TM_INCREMENT_PER_MOD = _cfg.LNA_TM_INCREMENT_PER_MOD
NN_HYBRID_INIT = _cfg.NN_HYBRID_INIT
NN_RNA_DNA_HYBRID = _cfg.NN_RNA_DNA_HYBRID
PHOSPHATE_CHARGE = _cfg.PHOSPHATE_CHARGE
PS_TM_DECREMENT_PER_MOD = _cfg.PS_TM_DECREMENT_PER_MOD
R_GAS = _cfg.R_GAS
RNA_COMPLEMENT = _cfg.RNA_COMPLEMENT
SL_RNA_SEQUENCE = _cfg.SL_RNA_SEQUENCE
STACKING_ENERGY = _cfg.STACKING_ENERGY
T_PHYSIOLOGICAL = _cfg.T_PHYSIOLOGICAL
StructuralAnalysisConfig = _cfg.StructuralAnalysisConfig


# ---------------------------------------------------------------------------
# Resultado da decomposicao de energia
# ---------------------------------------------------------------------------

@dataclass
class EnergyComponent:
    """Uma componente individual da energia de ligacao."""
    name: str
    value_kcal: float     # kcal/mol
    description: str = ""
    n_contributions: int = 0


@dataclass
class BindingEnergyResult:
    """Resultado completo da decomposicao de energia de ligacao."""
    # Componentes individuais
    stacking_energy: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("base_stacking", 0.0)
    )
    hbond_energy: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("hydrogen_bonding", 0.0)
    )
    electrostatic_energy: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("electrostatic", 0.0)
    )
    solvation_penalty: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("desolvation_penalty", 0.0)
    )
    lna_contribution: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("lna_stabilization", 0.0)
    )
    ps_contribution: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("ps_backbone", 0.0)
    )
    nn_thermodynamic: EnergyComponent = field(
        default_factory=lambda: EnergyComponent("nearest_neighbor", 0.0)
    )

    # Totais
    total_dg: float = 0.0            # kcal/mol — energia livre total
    total_dh: float = 0.0            # kcal/mol — entalpia total
    total_ds: float = 0.0            # cal/(mol*K) — entropia total
    predicted_tm: float = 0.0        # Celsius — temperatura de melting predita

    # Comparacao com valores experimentais
    experimental_dg: float = -27.97  # kcal/mol (do pipeline aso_math)
    experimental_tm: float = 68.48   # Celsius (do pipeline aso_math)
    dg_deviation: float = 0.0        # % desvio do experimental
    tm_deviation: float = 0.0        # Celsius desvio do experimental

    # Classificacao
    is_functional: bool = False      # dG < -15 kcal/mol (limiar para RNase H)
    binding_class: str = ""          # "strong", "moderate", "weak"

    # Detalhamento por posicao
    per_position_stacking: list[float] = field(default_factory=list)
    per_position_hbond: list[float] = field(default_factory=list)
    lna_vs_dna_comparison: dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Calculo nearest-neighbor (RNA:DNA hibrido)
# ---------------------------------------------------------------------------

def compute_nn_thermodynamics(
    rna_sequence: str,
    dna_sequence: str,
) -> tuple[float, float, float, float]:
    """Calcula dH, dS, dG e Tm usando parametros nearest-neighbor para hibrido RNA:DNA.

    Usa os parametros de Sugimoto et al. (1995) para hibridos RNA:DNA,
    que sao mais precisos que os parametros DNA:DNA de SantaLucia para ASOs.

    O ASO (DNA) hibrida com o SL RNA (RNA). A fita RNA e lida 5'->3' e
    o DNA e a fita complementar antiparalela.

    Args:
        rna_sequence: Sequencia RNA 5'->3' (letras ACGU).
        dna_sequence: Sequencia DNA complementar 5'->3' (letras ACGT).

    Returns:
        Tupla (dH kcal/mol, dS cal/(mol*K), dG kcal/mol, Tm Celsius).
    """
    # Garantir que as sequencias sao complementares e antiparalelas
    # RNA: 5' r1 r2 r3 ... rN 3'
    # DNA: 3' d1 d2 d3 ... dN 5'  (ou seja, DNA lido 5'->3' e a reversa)
    n = min(len(rna_sequence), len(dna_sequence))

    # Reverter o DNA para alinhar com o RNA (5'->3' vs 3'->5')
    dna_rev = dna_sequence[::-1]

    # Inicializacao
    dh_init, ds_init = NN_HYBRID_INIT
    total_dh = dh_init
    total_ds = ds_init

    # Iterar sobre dinucleotideos
    for i in range(n - 1):
        # Par de dinucleotideo: rXY/dX'Y'
        r1 = rna_sequence[i]
        r2 = rna_sequence[i + 1]
        d1 = dna_rev[i]
        d2 = dna_rev[i + 1]

        key = f"r{r1}{r2}/d{d1}{d2}"

        if key in NN_RNA_DNA_HYBRID:
            dh, ds = NN_RNA_DNA_HYBRID[key]
        else:
            # Fallback: usar media dos parametros disponiveis
            avg_dh = np.mean([v[0] for v in NN_RNA_DNA_HYBRID.values()])
            avg_ds = np.mean([v[1] for v in NN_RNA_DNA_HYBRID.values()])
            dh, ds = float(avg_dh), float(avg_ds)

        total_dh += dh
        total_ds += ds

    # dG a 37 C
    dg_37 = total_dh - (T_PHYSIOLOGICAL * total_ds / 1000.0)

    # Tm (assumindo concentracoes iguais de ASO e alvo, modelo two-state)
    # Tm = dH / (dS + R * ln(Ct/4))
    # Ct = concentracao total = 250 nM (padrao)
    ct = 250e-9  # mol/L
    tm_kelvin = (total_dh * 1000.0) / (total_ds + R_GAS * math.log(ct / 4.0))
    tm_celsius = tm_kelvin - 273.15

    return total_dh, total_ds, dg_37, tm_celsius


# ---------------------------------------------------------------------------
# Energia de stacking por posicao
# ---------------------------------------------------------------------------

def compute_stacking_energy(rna_sequence: str) -> tuple[float, list[float]]:
    """Calcula energia de stacking por posicao usando parametros simplificados.

    Stacking e a principal forca estabilizante em acidos nucleicos,
    contribuindo com ~60-70% da energia total de ligacao.

    Args:
        rna_sequence: Sequencia RNA 5'->3'.

    Returns:
        Tupla (energia total kcal/mol, lista de energias por posicao).
    """
    per_position: list[float] = []
    total = 0.0

    for i in range(len(rna_sequence) - 1):
        dinuc = rna_sequence[i] + rna_sequence[i + 1]
        energy = STACKING_ENERGY.get(dinuc, -1.0)  # default conservador
        per_position.append(energy)
        total += energy

    return total, per_position


# ---------------------------------------------------------------------------
# Energia de pontes de hidrogenio
# ---------------------------------------------------------------------------

def compute_hbond_energy(
    rna_sequence: str,
    dna_sequence: str,
) -> tuple[float, list[float]]:
    """Calcula energia de pontes de hidrogenio por par de bases.

    Cada par Watson-Crick contribui com energia proporcional ao
    numero de pontes H: G-C (~-6 kcal/mol, 3 H-bonds),
    A-T/U (~-4 kcal/mol, 2 H-bonds).

    Nota: estes valores representam a contribuicao entalpica das H-bonds,
    que e parcialmente compensada pela dessolvatacao. O efeito liquido
    e menor que os valores brutos.

    Args:
        rna_sequence: Sequencia RNA 5'->3'.
        dna_sequence: Sequencia DNA complementar 5'->3'.

    Returns:
        Tupla (energia total kcal/mol, lista por posicao).
    """
    n = min(len(rna_sequence), len(dna_sequence))
    dna_rev = dna_sequence[::-1]  # alinhar com RNA

    per_position: list[float] = []
    total = 0.0

    for i in range(n):
        rna_base = rna_sequence[i]
        dna_base = dna_rev[i]

        pair = rna_base + dna_base
        energy = HBOND_ENERGY.get(pair, -2.0)  # default para pares nao-canonicos
        per_position.append(energy)
        total += energy

    return total, per_position


# ---------------------------------------------------------------------------
# Contribuicao eletrostatica
# ---------------------------------------------------------------------------

def compute_electrostatic_energy(
    n_bp: int,
    use_ps: bool = True,
) -> float:
    """Estima contribuicao eletrostatica backbone-backbone.

    Os backbones de fosfato de ambas as fitas carregam carga -1 por residuo.
    A repulsao e parcialmente blindada por contra-ions (Na+, Mg2+) e
    pela alta constante dieletrica da agua.

    O backbone PS (fosforotioato) tem carga ligeiramente diferente do PO
    devido a substituicao S->O.

    Modelo simplificado: Debye-Huckel com blindagem ionica.

    Args:
        n_bp: Numero de pares de bases no duplex.
        use_ps: Se True, considera backbone PS.

    Returns:
        Energia eletrostatica em kcal/mol (positiva = repulsiva).
    """
    # Constantes do modelo Debye-Huckel
    # Comprimento de Debye em condicoes fisiologicas (~150 mM NaCl, 37 C)
    debye_length = 7.8  # Angstroms

    # Distancia media entre fosfatos na mesma fita (rise)
    p_p_distance = 5.9  # Angstroms (distancia P-P na forma-A)

    # Distancia entre fosfatos das fitas opostas
    cross_strand_distance = 17.0  # Angstroms (forma-A)

    # Constante de Coulomb (kcal*A/(mol*e^2)) em agua
    coulomb_constant = 332.0 / DIELECTRIC_WATER

    # Contribuicao intra-fita (repulsao dentro de cada fita)
    e_intra = 0.0
    for i in range(n_bp - 1):
        dist = p_p_distance * (i + 1)
        # Potencial de Debye-Huckel
        e_intra += coulomb_constant * math.exp(-dist / debye_length) / dist

    # Duas fitas contribuem
    e_intra *= 2.0

    # Contribuicao inter-fita (repulsao entre fitas)
    e_inter = 0.0
    for i in range(n_bp):
        for j in range(n_bp):
            dx = abs(i - j) * A_FORM_RISE
            dist = math.sqrt(cross_strand_distance ** 2 + dx ** 2)
            e_inter += coulomb_constant * math.exp(-dist / debye_length) / dist

    # Fator PS: fosforotioato reduz a carga efetiva em ~10%
    ps_factor = 0.90 if use_ps else 1.0

    total = (e_intra + e_inter) * ps_factor

    return total


# A_FORM_RISE ja importado via _cfg no topo do arquivo
A_FORM_RISE = _cfg.A_FORM_RISE


# ---------------------------------------------------------------------------
# Penalidade de dessolvatacao
# ---------------------------------------------------------------------------

def compute_desolvation_penalty(n_bp: int) -> float:
    """Estima penalidade de dessolvatacao na formacao do duplex.

    Quando o ASO hibrida com o SL RNA, moleculas de agua na interface
    sao deslocadas. Isso tem custo entropico (~0.2 kcal/mol por bp).

    Na pratica, a dessolvatacao e parcialmente compensada pelo ganho
    de entropia da agua liberada (efeito hidrofobico).

    Args:
        n_bp: Numero de pares de bases.

    Returns:
        Penalidade de dessolvatacao em kcal/mol (positiva = desfavoravel).
    """
    return n_bp * DESOLVATION_PENALTY_PER_BP


# ---------------------------------------------------------------------------
# Contribuicao LNA
# ---------------------------------------------------------------------------

def compute_lna_contribution(
    lna_positions: list[int],
    n_bp: int,
) -> tuple[float, dict[str, float]]:
    """Calcula a contribuicao das modificacoes LNA para a estabilidade.

    LNA (Locked Nucleic Acid) estabiliza o duplex por:
    1. Pre-organizacao conformacional: acucar constrainado C3'-endo
       reduz a penalidade entropica de ligacao
    2. Efeito de vizinhanca: LNA induz conformacao A-form nos vizinhos
    3. Incremento de Tm: +3 a +5 C por modificacao

    O efeito na dG pode ser estimado como:
    dG_LNA = -R * T * ln(K_LNA / K_DNA)
    Onde K_LNA/K_DNA ~ exp(dTm * dH / (R * Tm^2))

    Ref: Koshkin AA et al. (1998) Tetrahedron 54:3607-3630.
         Vester B & Wengel J (2004) Biochemistry 43:13233-13241.

    Args:
        lna_positions: Posicoes LNA (0-indexed).
        n_bp: Numero total de pares de bases.

    Returns:
        Tupla (dG_LNA kcal/mol, dicionario comparativo LNA vs DNA).
    """
    n_lna = len(lna_positions)

    if n_lna == 0:
        return 0.0, {"n_lna": 0, "dg_lna": 0.0, "dtm_lna": 0.0}

    # Incremento de Tm por LNA
    total_dtm = n_lna * LNA_TM_INCREMENT_PER_MOD  # graus Celsius

    # Converter dTm em contribuicao para dG
    # dG = dH * (1 - T / Tm), entao ddG ~ dH * T * dTm / Tm^2
    # Usando dH tipico de ~-8 kcal/mol por bp e Tm ~ 341 K (68 C)
    dh_per_bp = -8.0  # kcal/mol (valor medio)
    tm_ref = 273.15 + 68.48  # Kelvin (Tm experimental)

    # Contribuicao ao dG por estabilizacao LNA
    dg_lna = n_lna * dh_per_bp * T_PHYSIOLOGICAL * LNA_TM_INCREMENT_PER_MOD / (
        tm_ref ** 2
    )

    # Pre-organizacao conformacional: reducao da penalidade entropica
    # ~0.5 kcal/mol por LNA devido a pre-organizacao do acucar
    dg_preorg = -0.5 * n_lna

    dg_total_lna = dg_lna + dg_preorg

    comparison = {
        "n_lna": float(n_lna),
        "n_dna_gap": float(n_bp - n_lna),
        "dg_lna_total": dg_total_lna,
        "dg_preorganization": dg_preorg,
        "dg_tm_effect": dg_lna,
        "dtm_total": total_dtm,
        "fraction_lna": n_lna / n_bp if n_bp > 0 else 0.0,
    }

    return dg_total_lna, comparison


# ---------------------------------------------------------------------------
# Contribuicao PS (fosforotioato)
# ---------------------------------------------------------------------------

def compute_ps_contribution(
    n_ps: int,
    n_bp: int,
) -> float:
    """Calcula o efeito do backbone fosforotioato na estabilidade.

    PS reduz levemente a Tm (~-0.5 C por modificacao) mas aumenta
    drasticamente a resistencia a nucleases. O trade-off e favoravel
    porque LNA compensa a perda termica.

    Args:
        n_ps: Numero de ligacoes PS no gap central.
        n_bp: Numero total de pares de bases.

    Returns:
        Contribuicao ao dG em kcal/mol (positiva = desestabilizante).
    """
    if n_ps == 0:
        return 0.0

    # Efeito na Tm
    total_dtm = n_ps * PS_TM_DECREMENT_PER_MOD  # negativo

    # Converter dTm em dG (similar ao LNA mas com sinal oposto)
    dh_per_bp = -8.0
    tm_ref = 273.15 + 68.48

    dg_ps = n_ps * dh_per_bp * T_PHYSIOLOGICAL * PS_TM_DECREMENT_PER_MOD / (
        tm_ref ** 2
    )

    return dg_ps


# ---------------------------------------------------------------------------
# Decomposicao completa
# ---------------------------------------------------------------------------

def decompose_binding_energy(
    config: StructuralAnalysisConfig | None = None,
) -> BindingEnergyResult:
    """Executa decomposicao completa da energia de ligacao ASO:SL RNA.

    Calcula cada componente separadamente e compara o total com
    o valor experimental do pipeline aso_math (-27.97 kcal/mol).

    O design gapmer do MRL-ASO-001:
    - Flancos LNA (pos 0-2, 21-23): alta afinidade, resistencia a nucleases
    - Gap DNA (pos 3-20): permite recrutamento de RNase H
    - Backbone PS: protecao contra exonucleases

    Args:
        config: Configuracao do modulo.

    Returns:
        BindingEnergyResult com decomposicao completa.
    """
    if config is None:
        config = StructuralAnalysisConfig()

    result = BindingEnergyResult()

    # Sequencia RNA alvo (SL RNA na regiao de hibridizacao)
    target_rna = SL_RNA_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END]
    n_bp = min(len(target_rna), len(ASO_SEQUENCE))

    # --- 1. Nearest-neighbor thermodynamics (referencia principal) ---
    dh, ds, dg_nn, tm_nn = compute_nn_thermodynamics(target_rna, ASO_SEQUENCE)
    result.nn_thermodynamic = EnergyComponent(
        name="nearest_neighbor",
        value_kcal=dg_nn,
        description=(
            f"Parametros nearest-neighbor RNA:DNA hibrido (Sugimoto 1995). "
            f"dH={dh:.1f} kcal/mol, dS={ds:.1f} cal/(mol*K)"
        ),
        n_contributions=n_bp - 1,
    )
    result.total_dh = dh
    result.total_ds = ds

    # --- 2. Base stacking ---
    stacking_total, stacking_per_pos = compute_stacking_energy(target_rna)
    result.stacking_energy = EnergyComponent(
        name="base_stacking",
        value_kcal=stacking_total,
        description=(
            "Interacoes de stacking pi-pi entre bases consecutivas. "
            "Contribuicao dominante (~60-70% da estabilidade)."
        ),
        n_contributions=len(stacking_per_pos),
    )
    result.per_position_stacking = stacking_per_pos

    # --- 3. Pontes de hidrogenio ---
    hbond_total, hbond_per_pos = compute_hbond_energy(target_rna, ASO_SEQUENCE)
    result.hbond_energy = EnergyComponent(
        name="hydrogen_bonding",
        value_kcal=hbond_total,
        description=(
            "Pontes de hidrogenio Watson-Crick entre bases complementares. "
            f"G-C: 3 H-bonds (-6 kcal/mol), A-T/U: 2 H-bonds (-4 kcal/mol)."
        ),
        n_contributions=len(hbond_per_pos),
    )
    result.per_position_hbond = hbond_per_pos

    # --- 4. Eletrostatica ---
    e_elec = compute_electrostatic_energy(n_bp, use_ps=config.use_ps_backbone)
    result.electrostatic_energy = EnergyComponent(
        name="electrostatic",
        value_kcal=e_elec,
        description=(
            "Repulsao eletrostatica entre backbones de fosfato. "
            "Blindada por contra-ions em condicoes fisiologicas."
        ),
        n_contributions=n_bp,
    )

    # --- 5. Dessolvatacao ---
    e_desolv = compute_desolvation_penalty(n_bp)
    result.solvation_penalty = EnergyComponent(
        name="desolvation_penalty",
        value_kcal=e_desolv,
        description=(
            "Penalidade entropica por deslocamento de agua na interface. "
            f"~{DESOLVATION_PENALTY_PER_BP} kcal/mol por par de bases."
        ),
        n_contributions=n_bp,
    )

    # --- 6. Contribuicao LNA ---
    dg_lna, lna_comparison = compute_lna_contribution(
        list(config.lna_positions), n_bp
    )
    result.lna_contribution = EnergyComponent(
        name="lna_stabilization",
        value_kcal=dg_lna,
        description=(
            "Estabilizacao por modificacoes LNA nas extremidades do gapmer. "
            f"{len(config.lna_positions)} posicoes LNA "
            f"(+{LNA_TM_INCREMENT_PER_MOD} C/mod)."
        ),
        n_contributions=len(config.lna_positions),
    )
    result.lna_vs_dna_comparison = lna_comparison

    # --- 7. Contribuicao PS ---
    # Gap central: posicoes que NAO sao LNA
    n_ps = n_bp - len(config.lna_positions)
    dg_ps = compute_ps_contribution(n_ps, n_bp)
    result.ps_contribution = EnergyComponent(
        name="ps_backbone",
        value_kcal=dg_ps,
        description=(
            f"Efeito do backbone fosforotioato ({n_ps} ligacoes PS). "
            f"Leve desestabilizacao ({PS_TM_DECREMENT_PER_MOD} C/mod) "
            "compensada por resistencia a nucleases."
        ),
        n_contributions=n_ps,
    )

    # --- Total: usar nearest-neighbor como base, ajustar com LNA e PS ---
    result.total_dg = dg_nn + dg_lna + dg_ps

    # Tm ajustada por LNA e PS
    n_lna = len(config.lna_positions)
    tm_lna_adjustment = n_lna * LNA_TM_INCREMENT_PER_MOD
    tm_ps_adjustment = n_ps * PS_TM_DECREMENT_PER_MOD
    result.predicted_tm = tm_nn + tm_lna_adjustment + tm_ps_adjustment

    # --- Comparacao com experimental ---
    result.experimental_dg = -27.97
    result.experimental_tm = 68.48

    if abs(result.experimental_dg) > 0:
        result.dg_deviation = abs(
            (result.total_dg - result.experimental_dg) / result.experimental_dg * 100
        )

    result.tm_deviation = result.predicted_tm - result.experimental_tm

    # --- Classificacao ---
    result.is_functional = result.total_dg < -15.0
    if result.total_dg < -25.0:
        result.binding_class = "strong"
    elif result.total_dg < -15.0:
        result.binding_class = "moderate"
    else:
        result.binding_class = "weak"

    return result
