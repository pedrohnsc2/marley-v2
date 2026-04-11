"""Modelos de interacao ASO-membrana e propriedades fisico-quimicas.

Implementa tres componentes:
1. Propriedades fisico-quimicas do ASO (MW, logP, logD, carga, raio)
2. Modelo de interacao com bicamada lipidica (eletrostatico, hidrofobico, Born)
3. Coeficiente de particao membrana/agua

A motivacao biologica: oligonucleotideos fosforotioato (PS) sao
moleculas grandes (~8.5 kDa), altamente carregadas negativamente
(carga ~-24 a pH 7.4) e hidrofílicas. A membrana celular e uma
barreira fundamentalmente hostil para essas moleculas.

POREM, o backbone PS confere propriedades que facilitam a captacao:
- Aumento de hidrofobicidade vs PO (enxofre e menos polar que oxigenio)
- Ligacao a proteinas plasmaticas (albumina) que mediam entrega ao tecido
- Reconhecimento por receptores scavenger em macrofagos (tropismo natural)

Referencias:
- Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51 — farmacocinetica de PS-ASOs
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — mecanismos de uptake
- Liang XH et al. (2015) Nucleic Acids Res 43(5):2927-2945 — gymnosis e uptake
- Parsegian VA (1969) Nature 221(5183):844-846 — energia de Born em membranas
- McLaughlin S (1989) Annu Rev Biophys Biophys Chem 18:113-136 — eletrostatica
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Constantes fisico-quimicas
# ---------------------------------------------------------------------------

# Constante de Boltzmann em kcal/(mol*K)
KB_KCAL: Final[float] = 1.987e-3

# Temperatura fisiologica (37 C)
T_PHYSIOLOGICAL_K: Final[float] = 310.15

# Constante de Avogadro
AVOGADRO: Final[float] = 6.022e23

# Carga elementar em Coulombs
ELEMENTARY_CHARGE_C: Final[float] = 1.602e-19

# Permissividade do vacuo (C^2 / (N * m^2))
EPSILON_0: Final[float] = 8.854e-12

# Constante dieletrica da agua a 37 C
EPSILON_WATER: Final[float] = 74.0

# Constante dieletrica do interior da membrana lipidica
EPSILON_MEMBRANE: Final[float] = 2.0

# Espessura da bicamada lipidica (metros)
MEMBRANE_THICKNESS_M: Final[float] = 4.0e-9

# Potencial de superficie da membrana (mV, negativo por fosfolipideos)
# Ref: McLaughlin S (1989) — membranas biologicas tipicas: -20 a -40 mV
MEMBRANE_SURFACE_POTENTIAL_MV: Final[float] = -30.0

# Comprimento de Debye a forca ionica fisiologica (0.15 M NaCl, 37 C)
# kappa^-1 ~ 0.78 nm
DEBYE_LENGTH_M: Final[float] = 7.8e-10

# Fator de conversao: 1 kcal/mol = 4184 J/mol
KCAL_TO_J: Final[float] = 4184.0

# ---------------------------------------------------------------------------
# Pesos moleculares de componentes do ASO
# ---------------------------------------------------------------------------

# Peso molecular medio de nucleotideo com backbone fosfodiester (PO)
# Base media (~330 Da) + fosfato + acucar
MW_NUCLEOTIDE_PO: Final[float] = 330.0

# Adicao por ligacao PS: substituicao O -> S acrescenta ~16 Da por ligacao
# (S = 32, O = 16, net +16 Da)
MW_PS_ADDITION_PER_LINK: Final[float] = 16.0

# Adicao por modificacao LNA: ponte metileno C2'-C4' acrescenta ~14 Da
# (CH2 bridge)
MW_LNA_ADDITION_PER_RESIDUE: Final[float] = 14.0

# Contribuicao de logP por tipo de backbone (por ligacao)
# Ref: Geary RS et al. (2015) — PS e mais hidrofobico que PO
# PO: ~-1.0 por ligacao (muito polar, carregado)
# PS: ~-0.75 por ligacao (enxofre e menos polar que oxigenio)
LOGP_PER_LINK_PO: Final[float] = -1.0
LOGP_PER_LINK_PS: Final[float] = -0.75

# Contribuicao de LNA ao logP (restricao covalente reduz exposicao ao solvente)
LOGP_LNA_CORRECTION_PER_RESIDUE: Final[float] = 0.02

# pKa das ligacoes fosfato/fosforotioato para calculo de carga
# Ref: Eckstein F (2014) — pKa do fosfotioato ~1.5 (diester)
PKA_PS_LINKAGE: Final[float] = 1.5

# Raio de van der Waals medio por nucleotideo (Angstroms)
# Ref: baseado em modelos de acidos nucleicos — ~3.4 A entre bases empilhadas
VDW_RADIUS_PER_NT_A: Final[float] = 3.4

# Comprimento de persistencia de ssDNA/ssRNA (nanometros)
# Ref: Murphy MC et al. (2004) Biophys J 86(4):2530-2537
PERSISTENCE_LENGTH_NM: Final[float] = 1.5


# ---------------------------------------------------------------------------
# 1. Propriedades fisico-quimicas do ASO
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PhysicochemicalProperties:
    """Propriedades fisico-quimicas calculadas para um oligonucleotideo.

    Inclui peso molecular, logP/logD, carga liquida e raio hidrodinamico.
    Esses parametros determinam se a molecula pode interagir com e
    atravessar a membrana celular.

    A intuicao: moleculas com MW > 500 Da geralmente nao cruzam membranas
    por difusao passiva (regra de Lipinski). ASOs tem MW ~8500 Da e
    carga -24, tornando difusao passiva essencialmente impossivel.
    O uptake depende INTEIRAMENTE de mecanismos ativos (endocitose).
    """

    name: str
    length_nt: int
    n_ps_linkages: int
    n_lna_residues: int
    molecular_weight_da: float
    logp: float                    # coeficiente de particao octanol/agua
    logd_ph74: float               # logD a pH 7.4 (forma ionizada)
    logd_ph45: float               # logD a pH 4.5 (fagolisossomo)
    net_charge_ph74: float         # carga liquida a pH 7.4
    net_charge_ph45: float         # carga liquida a pH 4.5
    hydrodynamic_radius_nm: float  # raio hidrodinamico estimado
    lipinski_violations: int       # numero de violacoes da regra de 5


def compute_molecular_weight(
    length_nt: int,
    n_ps_linkages: int,
    n_lna_residues: int,
) -> float:
    """Calcula peso molecular de um oligonucleotideo com modificacoes PS e LNA.

    MW = (n_nt * MW_nucleotideo_base) + (n_PS * addicao_PS) + (n_LNA * addicao_LNA)

    Para MRL-ASO-001 (25 nt, 24 PS, 10 LNA):
    MW ~= 25 * 330 + 24 * 16 + 10 * 14 = 8250 + 384 + 140 = 8774 Da

    Nota: este e um calculo aproximado. O MW exato depende da composicao
    de bases especifica (A/T/G/C tem pesos diferentes). A aproximacao
    e suficiente para estimar propriedades de membrana.

    Args:
        length_nt: Comprimento em nucleotideos.
        n_ps_linkages: Numero de ligacoes fosforotioato.
        n_lna_residues: Numero de residuos LNA.

    Returns:
        Peso molecular em Daltons.
    """
    mw_base = length_nt * MW_NUCLEOTIDE_PO
    mw_ps = n_ps_linkages * MW_PS_ADDITION_PER_LINK
    mw_lna = n_lna_residues * MW_LNA_ADDITION_PER_RESIDUE
    return round(mw_base + mw_ps + mw_lna, 1)


def compute_logp(
    n_linkages: int,
    n_lna_residues: int,
    backbone_type: str = "PS",
) -> float:
    """Calcula logP (coeficiente de particao octanol/agua) para o ASO.

    Para oligonucleotideos, logP e extremamente negativo porque a
    molecula e altamente carregada e hidrofilica. O backbone PS
    e ligeiramente menos polar que PO (enxofre vs oxigenio).

    logP = sum(contribuicao_por_ligacao) + correcao_LNA

    Ref: Geary RS et al. (2015) — PS-ASOs tem logP ~ -15 a -25
    (dependendo do comprimento e composicao)

    Args:
        n_linkages: Numero de ligacoes internucleotidicas.
        n_lna_residues: Numero de residuos LNA.
        backbone_type: "PO" ou "PS".

    Returns:
        logP estimado (adimensional, muito negativo para ASOs).
    """
    logp_per_link = LOGP_PER_LINK_PS if backbone_type == "PS" else LOGP_PER_LINK_PO
    logp_backbone = n_linkages * logp_per_link
    logp_lna = n_lna_residues * LOGP_LNA_CORRECTION_PER_RESIDUE
    return round(logp_backbone + logp_lna, 2)


def compute_logd(logp: float, net_charge: float) -> float:
    """Calcula logD (coeficiente de distribuicao) considerando ionizacao.

    logD = logP - log10(1 + 10^|charge|)

    Para moleculas altamente carregadas (ASOs), logD << logP porque
    as cargas negativas favorecem fortemente a fase aquosa.

    A intuicao: em pH fisiologico, cada carga negativa adicional
    reduz logD em ~1 unidade (fator de 10x mais hidrofílico).

    Nota: para oligonucleotideos com ~24 cargas, logD pode chegar
    a -40 ou menos. Isso confirma que difusao passiva e impossivel.

    Args:
        logp: logP da molecula neutra.
        net_charge: Carga liquida no pH de interesse.

    Returns:
        logD no pH especificado.
    """
    # Para moleculas polianionicas, a correcao e dominada pela carga
    # Cada carga desfavorece a fase organica (octanol) por ~0.5 unidades de log
    # Ref: abordagem simplificada para polianions
    charge_correction = abs(net_charge) * 0.5
    return round(logp - charge_correction, 2)


def compute_net_charge(ph: float, n_ps_linkages: int) -> float:
    """Calcula carga liquida do ASO no pH especificado.

    Cada ligacao PS tem pKa ~1.5 e carrega carga -1 quando
    desprotonada (pH > pKa). Em pH fisiologico (7.4) e no
    fagolisossomo (4.5), praticamente todas as ligacoes estao
    desprotonadas porque pH >> pKa.

    carga_por_ligacao = -1 * fracao_desprotonada(pH, pKa)
    fracao_desprotonada = 1 / (1 + 10^(pKa - pH))

    Para pH 7.4: f_desprot = 1 / (1 + 10^(1.5-7.4)) = 0.999998 -> -24.0
    Para pH 4.5: f_desprot = 1 / (1 + 10^(1.5-4.5)) = 0.999 -> -23.98

    Conclusao: a carga e essencialmente -24 em qualquer pH biologico.

    Args:
        ph: pH do ambiente.
        n_ps_linkages: Numero de ligacoes PS.

    Returns:
        Carga liquida (negativa para ASOs).
    """
    f_deprotonated = 1.0 / (1.0 + 10.0 ** (PKA_PS_LINKAGE - ph))
    return round(-n_ps_linkages * f_deprotonated, 2)


def compute_hydrodynamic_radius(length_nt: int) -> float:
    """Estima raio hidrodinamico de um oligonucleotideo em solucao.

    Modelo de cadeia semiflexivel (worm-like chain):
    R_h = sqrt(L_p * L_c / 3)

    Onde:
    - L_p = comprimento de persistencia (~1.5 nm para ssDNA)
    - L_c = comprimento de contorno = n_nt * 0.59 nm (distancia entre fosfatos)

    Para MRL-ASO-001 (25 nt):
    L_c = 25 * 0.59 = 14.75 nm
    R_h = sqrt(1.5 * 14.75 / 3) = sqrt(7.375) = 2.72 nm

    Ref: Murphy MC et al. (2004) Biophys J 86(4):2530-2537

    Args:
        length_nt: Comprimento em nucleotideos.

    Returns:
        Raio hidrodinamico em nanometros.
    """
    # Distancia entre fosfatos consecutivos em ssDNA: ~0.59 nm
    phosphate_spacing_nm = 0.59
    contour_length = length_nt * phosphate_spacing_nm
    rh = math.sqrt(PERSISTENCE_LENGTH_NM * contour_length / 3.0)
    return round(rh, 2)


def compute_lipinski_violations(
    mw: float,
    logp: float,
    net_charge: float,
) -> int:
    """Conta violacoes da regra dos 5 de Lipinski (adaptada para ASOs).

    Regra original (moleculas pequenas):
    1. MW <= 500 Da
    2. logP <= 5
    3. H-bond donors <= 5
    4. H-bond acceptors <= 10

    ASOs violam TODAS essas regras porque sao macromoleculas biologicas,
    nao moleculas pequenas. O ponto e quantificar o quanto o ASO
    esta fora do espaco de drogas convencionais para justificar
    a necessidade de mecanismos de uptake ativos.

    Adaptacao para ASOs:
    1. MW > 500 Da -> violacao
    2. |logP| > 5 -> violacao (extremamente polar OU apolar)
    3. |carga| > 3 -> violacao (altamente carregado)
    4. MW > 1000 Da -> violacao adicional (macromolecula)

    Args:
        mw: Peso molecular em Daltons.
        logp: logP calculado.
        net_charge: Carga liquida.

    Returns:
        Numero de violacoes (0-4).
    """
    violations = 0
    if mw > 500:
        violations += 1
    if abs(logp) > 5:
        violations += 1
    if abs(net_charge) > 3:
        violations += 1
    if mw > 1000:
        violations += 1
    return violations


def compute_physicochemical_properties(
    name: str,
    length_nt: int,
    n_ps_linkages: int,
    n_lna_residues: int,
    backbone_type: str = "PS",
) -> PhysicochemicalProperties:
    """Calcula todas as propriedades fisico-quimicas de um ASO.

    Agrega todos os calculos individuais em um unico perfil completo.
    Este e o ponto de entrada principal para a secao 1 do modulo.

    Args:
        name: Identificador do ASO.
        length_nt: Comprimento em nucleotideos.
        n_ps_linkages: Numero de ligacoes PS.
        n_lna_residues: Numero de residuos LNA.
        backbone_type: "PO" ou "PS".

    Returns:
        Perfil completo de propriedades fisico-quimicas.
    """
    mw = compute_molecular_weight(length_nt, n_ps_linkages, n_lna_residues)
    logp = compute_logp(length_nt - 1, n_lna_residues, backbone_type)
    charge_74 = compute_net_charge(7.4, n_ps_linkages)
    charge_45 = compute_net_charge(4.5, n_ps_linkages)
    logd_74 = compute_logd(logp, charge_74)
    logd_45 = compute_logd(logp, charge_45)
    rh = compute_hydrodynamic_radius(length_nt)
    violations = compute_lipinski_violations(mw, logp, charge_74)

    return PhysicochemicalProperties(
        name=name,
        length_nt=length_nt,
        n_ps_linkages=n_ps_linkages,
        n_lna_residues=n_lna_residues,
        molecular_weight_da=mw,
        logp=logp,
        logd_ph74=logd_74,
        logd_ph45=logd_45,
        net_charge_ph74=charge_74,
        net_charge_ph45=charge_45,
        hydrodynamic_radius_nm=rh,
        lipinski_violations=violations,
    )


# ---------------------------------------------------------------------------
# 2. Modelo de interacao com bicamada lipidica
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MembraneInteractionEnergy:
    """Energias de interacao entre ASO e bicamada lipidica.

    Tres termos de energia sao calculados:
    1. Eletrostatico: repulsao entre ASO(-24) e membrana(-30 mV)
    2. Hidrofobico: contribuicao favoravel do backbone PS (parcialmente apolar)
    3. Born (solvatacao): custo de transferir cargas de agua para lipideo

    A soma desses termos determina a barreira energetica total para
    insercao do ASO na membrana. Uma barreira alta (>> kT) significa
    que difusao passiva e cineticamente impossivel.

    Ref: Parsegian VA (1969) Nature 221:844-846 — modelo de Born
         McLaughlin S (1989) Annu Rev Biophys Biophys Chem 18:113-136
    """

    name: str
    backbone_type: str
    electrostatic_kcal: float     # repulsao eletrostatica (positivo)
    hydrophobic_kcal: float       # insercao hidrofobica (negativo = favoravel)
    born_solvation_kcal: float    # custo de desolvatacao (positivo)
    total_barrier_kcal: float     # barreira total (positivo = desfavoravel)
    barrier_in_kt: float          # barreira em unidades de kT
    passive_diffusion_feasible: bool  # True se barreira < 20 kT


def compute_electrostatic_repulsion(
    net_charge: float,
    membrane_potential_mv: float = MEMBRANE_SURFACE_POTENTIAL_MV,
) -> float:
    """Calcula energia de repulsao eletrostatica ASO-membrana.

    O ASO (carga ~-24) se aproxima de uma membrana com potencial
    de superficie negativo (~-30 mV). A energia de interacao e:

    E_elec = z * e * psi_0 * N_A / 1000

    Onde:
    - z = carga do ASO (em unidades elementares)
    - e = carga elementar (1.602e-19 C)
    - psi_0 = potencial de superficie (V)
    - N_A = numero de Avogadro
    - Dividido por 1000 para converter J -> kJ, depois por 4.184 para kcal

    A repulsao e FORTE porque ambos sao negativos: molecula negativa
    se aproximando de superficie negativa. Porem, o screening ionico
    (ions Na+, K+ no meio) reduz essa repulsao exponencialmente com
    a distancia (comprimento de Debye).

    Aplicamos fator de screening de Debye para distancia de contato:
    E_screened = E_bare * exp(-d / lambda_D)

    onde d = raio do ASO (~2.7 nm) e lambda_D = 0.78 nm.

    Ref: McLaughlin S (1989) — eletrostatica de membranas biologicas

    Args:
        net_charge: Carga liquida do ASO (negativa).
        membrane_potential_mv: Potencial de superficie em mV.

    Returns:
        Energia de repulsao eletrostatica em kcal/mol (positivo).
    """
    psi_v = membrane_potential_mv / 1000.0  # mV -> V

    # Energia de interacao carga-superficie (sem screening)
    # E = z * e * psi * N_A (em Joules/mol)
    energy_j_per_mol = net_charge * ELEMENTARY_CHARGE_C * psi_v * AVOGADRO
    energy_kcal = energy_j_per_mol / KCAL_TO_J

    # Screening de Debye: a forca ionica fisiologica reduz a interacao
    # Para distancia de contato ASO-membrana (~2.7 nm), o fator e:
    contact_distance_m = 2.7e-9
    screening_factor = math.exp(-contact_distance_m / DEBYE_LENGTH_M)

    # A interacao e repulsiva (carga negativa * potencial negativo = positivo)
    # mas a formula da energy_j_per_mol com sinais negativos da valor positivo
    # ((-24) * 1.602e-19 * (-0.03) * 6.022e23 = +694 J/mol)
    return round(abs(energy_kcal) * screening_factor, 4)


def compute_hydrophobic_insertion(
    n_linkages: int,
    backbone_type: str = "PS",
) -> float:
    """Calcula energia de insercao hidrofobica do backbone na membrana.

    O backbone PS tem carater hidrofobico PARCIAL devido ao enxofre:
    - O atomo S e maior e mais polarizavel que O
    - A distribuicao de carga e P-S e mais difusa que em P-O
    - Isso confere afinidade parcial pelo interior apolar da membrana

    Modelo simplificado:
    E_hydrophobic = n_ligacoes * delta_G_transferencia

    Para PO: delta_G ~= +0.3 kcal/mol por ligacao (desfavoravel, muito polar)
    Para PS: delta_G ~= -0.1 kcal/mol por ligacao (levemente favoravel)

    Nota: a contribuicao hidrofobica do PS e PEQUENA comparada a barreira
    eletrostatica e de Born, mas e biologicamente relevante porque contribui
    para a ligacao a proteinas plasmaticas e receptores scavenger.

    Ref: Geary RS et al. (2015) — PS confere propriedades protein-binding
         Liang XH et al. (2015) — interacao PS com proteinas celulares

    Args:
        n_linkages: Numero de ligacoes internucleotidicas.
        backbone_type: "PO" ou "PS".

    Returns:
        Energia de insercao em kcal/mol (negativo = favoravel).
    """
    # Delta G de transferencia agua -> lipideo por ligacao
    if backbone_type == "PS":
        dg_per_link = -0.1  # levemente favoravel (enxofre e parcialmente apolar)
    else:
        dg_per_link = 0.3   # desfavoravel (fosfato e muito polar)

    return round(n_linkages * dg_per_link, 4)


def compute_born_solvation_energy(
    net_charge: float,
    radius_nm: float,
) -> float:
    """Calcula energia de solvatacao de Born para transferencia agua -> membrana.

    O modelo de Born descreve o custo energetico de transferir uma
    carga de um meio com alta constante dieletrica (agua, eps=74)
    para um meio com baixa constante dieletrica (membrana, eps=2).

    E_Born = (z^2 * e^2 * N_A) / (8 * pi * eps_0 * r) * (1/eps_m - 1/eps_w)

    Para um ASO com carga -24 e raio ~2.7 nm:
    E_Born >> 100 kcal/mol

    Este e o termo DOMINANTE que impede difusao passiva de ASOs.
    A membrana lipidica simplesmente nao pode acomodar 24 cargas
    negativas em seu interior apolar (eps=2).

    Nota: o calculo assume transferencia completa para o interior
    da membrana. Na realidade, o ASO interage com a superficie
    (headgroups), nao penetra completamente. Portanto, usamos
    um fator de correcao (fracao efetiva de inserção ~10%) para
    obter uma estimativa mais realista da barreira.

    Ref: Parsegian VA (1969) Nature 221:844-846

    Args:
        net_charge: Carga liquida do ASO.
        radius_nm: Raio efetivo em nanometros.

    Returns:
        Energia de Born em kcal/mol (positivo = custo energetico).
    """
    radius_m = radius_nm * 1e-9

    # Fator de Born: (z^2 * e^2 * N_A) / (8 * pi * eps_0 * r)
    numerator = net_charge**2 * ELEMENTARY_CHARGE_C**2 * AVOGADRO
    denominator = 8.0 * math.pi * EPSILON_0 * radius_m

    # Diferenca de dieletrico: (1/eps_membrana - 1/eps_agua)
    dielectric_diff = (1.0 / EPSILON_MEMBRANE) - (1.0 / EPSILON_WATER)

    # Energia total em J/mol -> kcal/mol
    energy_j = (numerator / denominator) * dielectric_diff
    energy_kcal = energy_j / KCAL_TO_J

    # Fator de correcao: apenas ~10% da molecula interage efetivamente
    # com o interior hidrofobico (a maioria fica na superficie aquosa)
    insertion_fraction = 0.10

    return round(energy_kcal * insertion_fraction, 2)


def compute_membrane_interaction(
    name: str,
    net_charge: float,
    n_linkages: int,
    radius_nm: float,
    backbone_type: str = "PS",
) -> MembraneInteractionEnergy:
    """Calcula perfil completo de interacao ASO-membrana.

    Combina os tres termos de energia (eletrostatico, hidrofobico, Born)
    para determinar a barreira total de insercao na membrana.

    Uma barreira > 20 kT (~12 kcal/mol) e considerada intransponivel
    por difusao termica, exigindo mecanismos ativos de uptake.

    Args:
        name: Identificador do ASO/quimica.
        net_charge: Carga liquida no pH de interesse.
        n_linkages: Numero de ligacoes internucleotidicas.
        radius_nm: Raio hidrodinamico em nm.
        backbone_type: "PO" ou "PS".

    Returns:
        Perfil de energias de interacao com a membrana.
    """
    e_elec = compute_electrostatic_repulsion(net_charge)
    e_hydro = compute_hydrophobic_insertion(n_linkages, backbone_type)
    e_born = compute_born_solvation_energy(net_charge, radius_nm)

    total = e_elec + e_hydro + e_born

    # Converter para unidades kT (a 37 C)
    kt_kcal = KB_KCAL * T_PHYSIOLOGICAL_K
    barrier_kt = total / kt_kcal

    # Difusao passiva e viavel se barreira < 20 kT
    passive_ok = barrier_kt < 20.0

    return MembraneInteractionEnergy(
        name=name,
        backbone_type=backbone_type,
        electrostatic_kcal=e_elec,
        hydrophobic_kcal=e_hydro,
        born_solvation_kcal=e_born,
        total_barrier_kcal=round(total, 2),
        barrier_in_kt=round(barrier_kt, 1),
        passive_diffusion_feasible=passive_ok,
    )


# ---------------------------------------------------------------------------
# 3. Coeficiente de particao membrana/agua
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PartitionCoefficient:
    """Coeficiente de particao membrana/agua para diferentes quimicas.

    K_partition = exp(-dG_barrier / RT)

    Valores tipicos:
    - Moleculas pequenas hidrofobicas: K ~ 10^1 a 10^3
    - Peptideos: K ~ 10^-2 a 10^1
    - ASOs PO: K ~ 10^-50 (essencialmente zero)
    - ASOs PS: K ~ 10^-40 (marginalmente melhor, mas ainda zero)

    A conclusao e clara: particao direta e impossivel.
    O uptake ocorre EXCLUSIVAMENTE por endocitose.
    """

    name: str
    backbone_type: str
    barrier_kcal: float
    log_k_partition: float
    k_partition: float       # se > 1e-300, senao 0.0 (underflow)
    passive_permeability_cm_s: float  # estimativa de permeabilidade


def compute_partition_coefficient(
    name: str,
    backbone_type: str,
    barrier_kcal: float,
) -> PartitionCoefficient:
    """Calcula coeficiente de particao membrana/agua.

    K = exp(-dG / RT)

    Para barreira de 50 kcal/mol a 37 C:
    K = exp(-50 / 0.616) = exp(-81.2) = 2.4e-36

    Isso significa que para cada 10^36 moleculas na agua,
    UMA estaria na membrana no equilibrio. Em termos praticos,
    nenhuma molecula de ASO cruza a membrana por difusao.

    A permeabilidade e estimada pelo modelo de solubilidade-difusao:
    P = K * D / delta

    Onde:
    - K = coeficiente de particao
    - D = coeficiente de difusao na membrana (~1e-8 cm^2/s para macromoleculas)
    - delta = espessura da membrana (4 nm = 4e-7 cm)

    Ref: Stein WD (1986) Transport and Diffusion Across Cell Membranes

    Args:
        name: Identificador.
        backbone_type: Tipo de backbone.
        barrier_kcal: Barreira energetica total em kcal/mol.

    Returns:
        Coeficiente de particao e permeabilidade estimada.
    """
    rt = KB_KCAL * T_PHYSIOLOGICAL_K  # ~0.616 kcal/mol

    # log10(K) = -dG / (2.303 * RT)
    log_k = -barrier_kcal / (2.303 * rt)

    # K numerico (protecao contra underflow)
    exponent = -barrier_kcal / rt
    if exponent < -700:
        k_value = 0.0
    else:
        k_value = math.exp(exponent)

    # Permeabilidade: P = K * D / delta
    # D para macromoleculas na membrana ~1e-8 cm^2/s
    # delta = 4 nm = 4e-7 cm
    d_membrane = 1.0e-8   # cm^2/s
    delta_cm = MEMBRANE_THICKNESS_M * 100.0  # m -> cm
    if k_value > 0:
        perm = k_value * d_membrane / delta_cm
    else:
        perm = 0.0

    return PartitionCoefficient(
        name=name,
        backbone_type=backbone_type,
        barrier_kcal=barrier_kcal,
        log_k_partition=round(log_k, 1),
        k_partition=k_value,
        passive_permeability_cm_s=perm,
    )
