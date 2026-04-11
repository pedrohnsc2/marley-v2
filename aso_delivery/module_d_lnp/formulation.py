"""Modelos de formulacao LNP para encapsulacao de ASOs.

Implementa quatro modelos independentes:
1. Composicao lipidica padrao (4 componentes, razoes molares)
2. Otimizacao de razao N/P (eficiencia de encapsulacao)
3. Modificacoes para targeting de macrofagos (manose-PEG)
4. Estabilidade e armazenamento (taxa de agregacao termica)

A motivacao biologica: LNPs sao o vetor de entrega mais eficaz
para acidos nucleicos. Para MRL-ASO-001, a vantagem e dupla:
(a) protecao contra degradacao no soro, e (b) entrega seletiva
a macrofagos via receptores de manose. Como o alvo (SL RNA)
reside no fagolisossomo, a liberacao pH-responsiva do lipidio
ionizavel e vantajosa — a LNP libera a carga EXATAMENTE no
compartimento onde o ASO precisa atuar.

Referencias:
- Cullis PR, Hope MJ (2017) Mol Ther 25(7):1467-1475 — LNP design
- Jayaraman M et al. (2012) Angew Chem 51(34):8529-8533 — DLin-MC3-DMA
- Semple SC et al. (2010) Nat Biotechnol 28(2):172-176 — LNP optimization
- Patel S et al. (2020) J Control Release 327:146-160 — macrophage targeting
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes dos componentes lipidicos
# ---------------------------------------------------------------------------

# Pesos moleculares dos componentes (g/mol)
# Ref: PubChem CIDs / fabricantes
MW_DLIN_MC3_DMA: Final[float] = 642.09   # DLin-MC3-DMA (lipidio ionizavel)
MW_DSPC: Final[float] = 790.15           # 1,2-distearoil-sn-glicero-3-fosfocolina
MW_CHOLESTEROL: Final[float] = 386.65    # colesterol
MW_PEG_LIPID: Final[float] = 2693.0     # DMG-PEG2000 (PEG-lipidio padrao)
MW_MANNOSE_PEG: Final[float] = 2855.0   # manose-PEG2000-DSPE (com ligante mannose)

# Razoes molares padrao (Onpattro-like, Jayaraman 2012)
MOLAR_RATIO_IONIZABLE: Final[float] = 0.50
MOLAR_RATIO_DSPC: Final[float] = 0.10
MOLAR_RATIO_CHOLESTEROL: Final[float] = 0.385
MOLAR_RATIO_PEG: Final[float] = 0.015

# pKa do lipidio ionizavel — parametro critico para liberacao pH-responsiva
# Ref: Jayaraman M et al. (2012) — pKa otimo entre 6.2-6.5 para
# eficacia in vivo maxima; DLin-MC3-DMA tem pKa 6.44
PKA_DLIN_MC3: Final[float] = 6.44

# Numero de aminas ionizaveis por molecula de DLin-MC3-DMA
AMINES_PER_MC3: Final[int] = 1

# Peso molecular medio de um nucleotideo PS (com backbone fosforotioato)
# Fosforo + base + acucar + enxofre ~ 345 g/mol medio
MW_NT_PS: Final[float] = 345.0

# Constante de encapsulacao — calibrada para PS-ASOs
# Para mRNA: k ~0.5-0.8; para PS-ASOs: ligacao mais forte ao lipidio
# ionizavel devido a carga negativa permanente do backbone PS
# Ref: Semple et al. (2010) Nat Biotechnol 28(2):172-176
K_ENCAPSULATION: Final[float] = 0.55

# Parametros de tamanho de particula (nm)
# Diametro base: tamanho minimo da LNP vazia
# Crescimento por N/P: particulas maiores com excesso de lipidio
BASE_DIAMETER_NM: Final[float] = 65.0
DIAMETER_GROWTH_PER_NP: Final[float] = 2.5   # nm por unidade de N/P

# Limite de tamanho para captacao fagocitica
MIN_PHAGOCYTIC_NM: Final[float] = 50.0
MAX_PHAGOCYTIC_NM: Final[float] = 200.0

# Fator de aumento de uptake por manose (receptor CD206 de macrofagos)
# Ref: Patel S et al. (2020) — 2-5x para macrofagos M2 (peritoneais)
# Macrofagos caninos infectados por Leishmania expressam CD206 elevado
MANNOSE_UPTAKE_FACTOR_MIN: Final[float] = 2.0
MANNOSE_UPTAKE_FACTOR_MAX: Final[float] = 5.0
MANNOSE_UPTAKE_FACTOR_MEAN: Final[float] = 3.5

# Parametros de estabilidade / armazenamento
# Taxa de agregacao (k_agg em 1/dia) a 25 C como referencia
# Ref: Schoenmaker L et al. (2021) Int J Pharm 601:120586
K_AGG_25C: Final[float] = 0.015     # 1/dia a 25 C
EA_AGGREGATION: Final[float] = 60.0  # kJ/mol — energia de ativacao (Arrhenius)
R_GAS_KJ: Final[float] = 8.314e-3    # kJ/(mol*K) — constante dos gases

# Temperatura de referencia (25 C em Kelvin)
T_REF_K: Final[float] = 298.15


# ---------------------------------------------------------------------------
# Dataclasses de resultado
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LNPComposition:
    """Composicao de uma formulacao LNP com 4 componentes.

    Armazena razoes molares, pesos moleculares efetivos,
    e propriedades fisico-quimicas calculadas.
    """

    # Componentes
    ionizable_lipid: str
    helper_lipid: str
    sterol: str
    peg_lipid: str

    # Razoes molares (somam 1.0)
    mol_fraction_ionizable: float
    mol_fraction_helper: float
    mol_fraction_sterol: float
    mol_fraction_peg: float

    # Peso molecular medio ponderado (g/mol)
    weighted_average_mw: float

    # Cargas a diferentes pH
    charge_ph_7_4: float   # carga liquida a pH 7.4 (sangue)
    charge_ph_4_5: float   # carga liquida a pH 4.5 (fagolisossomo)

    # pKa do lipidio ionizavel
    pka_ionizable: float


@dataclass(frozen=True)
class NPRatioResult:
    """Resultado da analise de uma razao N/P especifica.

    N/P = razao entre aminas (N) cationicas do lipidio ionizavel e
    fosfatos (P) aniônicos do ASO. Parametro critico para:
    - Eficiencia de encapsulacao (muito baixo -> pouco complexo)
    - Tamanho da particula (muito alto -> agregacao)
    - Toxicidade (excesso de lipidio cationico e citotoxico)
    """

    np_ratio: float
    n_amines: float             # moles de amina
    n_phosphates: int           # moles de fosfato (= comprimento ASO - 1)
    encapsulation_efficiency: float   # 0-1
    particle_diameter_nm: float
    zeta_potential_mv: float    # mV (potencial zeta)
    suitable_for_phagocytosis: bool   # 50-200 nm?
    optimal: bool               # na faixa otima para PS-ASOs?


@dataclass(frozen=True)
class MacrophageTargeting:
    """Resultado da analise de targeting com manose-PEG.

    Macrofagos infectados por L. infantum expressam receptor de manose
    (CD206) em nivel elevado. Substituir PEG-lipidio convencional por
    manose-PEG-DSPE permite reconhecimento ativo pelo receptor.
    """

    mannose_peg_fraction: float       # fracao de PEG-lipidio que e manose-PEG
    uptake_fold_increase: float       # fator de aumento de captacao
    receptor_target: str              # ex: "CD206 (mannose receptor)"
    particle_diameter_nm: float       # verificar se ainda e fagocitavel
    suitable_for_phagocytosis: bool


@dataclass(frozen=True)
class StorageStability:
    """Resultado da analise de estabilidade em armazenamento.

    Modela agregacao por Arrhenius: k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T))
    Considera tres temperaturas clinicamente relevantes:
    - -20 C: armazenamento longo prazo
    -   4 C: refrigerado (cadeia fria)
    -  25 C: temperatura ambiente (relevante para Brasil)
    """

    temperature_celsius: float
    temperature_kelvin: float
    k_aggregation_per_day: float
    half_life_days: float
    fraction_intact_30_days: float
    fraction_intact_90_days: float
    fraction_intact_365_days: float
    lyophilization_feasible: bool


# ---------------------------------------------------------------------------
# 1. Composicao lipidica
# ---------------------------------------------------------------------------


def compute_lnp_composition() -> LNPComposition:
    """Calcula composicao e propriedades da formulacao LNP padrao.

    Formulacao baseada em Onpattro (patisiran), adaptada para ASO:
    - DLin-MC3-DMA (50%) — lipidio ionizavel (neutro pH 7.4, cationico pH <6.5)
    - DSPC (10%) — lipidio auxiliar estrutural (bicamada rigida)
    - Colesterol (38.5%) — estabilidade de membrana e fusogenicidade
    - DMG-PEG2000 (1.5%) — stealth (evita opsonizacao e clearance pelo RES)

    A carga liquida muda com pH porque DLin-MC3-DMA protona abaixo
    do pKa (6.44). A pH 7.4: neutro. A pH 4.5: cationico.
    Isto permite complexacao com ASO aniônico em pH acido durante
    formulacao, e liberacao no endossomo/fagolisossomo.

    Returns:
        LNPComposition com todas as propriedades calculadas.
    """
    # Peso molecular medio ponderado pela fracao molar
    avg_mw = (
        MOLAR_RATIO_IONIZABLE * MW_DLIN_MC3_DMA
        + MOLAR_RATIO_DSPC * MW_DSPC
        + MOLAR_RATIO_CHOLESTEROL * MW_CHOLESTEROL
        + MOLAR_RATIO_PEG * MW_PEG_LIPID
    )

    # Carga a pH 7.4: fracao protonada do lipidio ionizavel via Henderson-Hasselbalch
    # f_protonated = 1 / (1 + 10^(pH - pKa))
    f_prot_7_4 = 1.0 / (1.0 + 10.0 ** (7.4 - PKA_DLIN_MC3))
    charge_7_4 = MOLAR_RATIO_IONIZABLE * AMINES_PER_MC3 * f_prot_7_4

    # Carga a pH 4.5: quase totalmente protonado
    f_prot_4_5 = 1.0 / (1.0 + 10.0 ** (4.5 - PKA_DLIN_MC3))
    charge_4_5 = MOLAR_RATIO_IONIZABLE * AMINES_PER_MC3 * f_prot_4_5

    return LNPComposition(
        ionizable_lipid="DLin-MC3-DMA",
        helper_lipid="DSPC",
        sterol="Cholesterol",
        peg_lipid="DMG-PEG2000",
        mol_fraction_ionizable=MOLAR_RATIO_IONIZABLE,
        mol_fraction_helper=MOLAR_RATIO_DSPC,
        mol_fraction_sterol=MOLAR_RATIO_CHOLESTEROL,
        mol_fraction_peg=MOLAR_RATIO_PEG,
        weighted_average_mw=round(avg_mw, 2),
        charge_ph_7_4=round(charge_7_4, 6),
        charge_ph_4_5=round(charge_4_5, 6),
        pka_ionizable=PKA_DLIN_MC3,
    )


# ---------------------------------------------------------------------------
# 2. Otimizacao de razao N/P
# ---------------------------------------------------------------------------


def compute_np_ratio(
    np_ratio: float,
    aso_length: int,
) -> NPRatioResult:
    """Calcula propriedades da LNP para uma dada razao N/P.

    N/P = (moles de amina ionizavel) / (moles de fosfato no ASO)

    Para ASOs com backbone PS, cada ligacao internucleotidica carrega
    uma carga negativa (-1). Um ASO de 25 nt tem 24 ligacoes PS,
    portanto 24 cargas negativas.

    Modelo de encapsulacao:
        EE = 1 - exp(-k * N/P)
    onde k = 0.55 para PS-ASOs (maior que mRNA porque PS tem
    maior afinidade eletrostatica com lipidio cationico).

    Modelo de tamanho:
        d = d_base + growth * N/P
    Particulas maiores com mais lipidio. Limite de 50-200 nm para
    captacao fagocitica por macrofagos.

    Modelo de potencial zeta (simplificado):
        zeta = -30 + 10 * (N/P - 1) para N/P < 4
        zeta satura em ~+10 mV para N/P alto
    Em pH 7.4, zeta deve ser proximo de neutro (stealth).

    Args:
        np_ratio: Razao N/P a avaliar.
        aso_length: Comprimento do ASO em nucleotideos.

    Returns:
        NPRatioResult com metricas calculadas.
    """
    n_phosphates = aso_length - 1  # ligacoes internucleotidicas

    # Aminas totais para esta razao N/P
    n_amines = np_ratio * n_phosphates

    # Eficiencia de encapsulacao: modelo exponencial saturante
    ee = 1.0 - math.exp(-K_ENCAPSULATION * np_ratio)
    ee = min(ee, 0.99)  # teto fisico: nunca 100% perfeito

    # Tamanho da particula (nm)
    diameter = BASE_DIAMETER_NM + DIAMETER_GROWTH_PER_NP * np_ratio

    # Potencial zeta (mV) — modelo sigmoidal simplificado
    # A pH 7.4 (onde medimos zeta), MC3 esta maioritariamente neutro
    # mas excesso de lipidio ionizavel gera carga residual positiva
    zeta = -30.0 + 40.0 / (1.0 + math.exp(-0.8 * (np_ratio - 3.0)))

    # Critérios de avaliacao
    suitable = MIN_PHAGOCYTIC_NM <= diameter <= MAX_PHAGOCYTIC_NM
    # Faixa otima para PS-ASOs: N/P 4-8 (diferente de mRNA que e 6-12)
    optimal = 4.0 <= np_ratio <= 8.0

    return NPRatioResult(
        np_ratio=np_ratio,
        n_amines=round(n_amines, 1),
        n_phosphates=n_phosphates,
        encapsulation_efficiency=round(ee, 4),
        particle_diameter_nm=round(diameter, 1),
        zeta_potential_mv=round(zeta, 2),
        suitable_for_phagocytosis=suitable,
        optimal=optimal,
    )


def scan_np_ratios(
    aso_length: int,
    np_min: float = 1.0,
    np_max: float = 20.0,
    step: float = 1.0,
) -> list[NPRatioResult]:
    """Varre razoes N/P e retorna resultados para cada uma.

    Permite identificar o ponto otimo de formulacao:
    - N/P muito baixo -> encapsulacao pobre
    - N/P muito alto -> particulas grandes, toxicas, custo alto
    - Otimo para PS-ASOs: tipicamente N/P 4-8

    Args:
        aso_length: Comprimento do ASO (nucleotideos).
        np_min: N/P minimo a testar.
        np_max: N/P maximo a testar.
        step: Incremento de N/P.

    Returns:
        Lista de NPRatioResult ordenada por N/P.
    """
    results: list[NPRatioResult] = []
    np = np_min
    while np <= np_max + 1e-9:
        results.append(compute_np_ratio(np, aso_length))
        np += step
    return results


def find_optimal_np(results: list[NPRatioResult]) -> NPRatioResult:
    """Seleciona a razao N/P otima da varredura.

    Criterio: maior encapsulacao dentro da faixa otima (N/P 4-8)
    e tamanho compativel com fagocitose.

    Args:
        results: Lista de NPRatioResult da varredura.

    Returns:
        NPRatioResult com o melhor N/P.
    """
    optimal_candidates = [
        r for r in results
        if r.optimal and r.suitable_for_phagocytosis
    ]
    if not optimal_candidates:
        # Fallback: melhor encapsulacao geral com tamanho aceitavel
        optimal_candidates = [r for r in results if r.suitable_for_phagocytosis]
    if not optimal_candidates:
        # Ultimo recurso: melhor encapsulacao absoluta
        optimal_candidates = results

    return max(optimal_candidates, key=lambda r: r.encapsulation_efficiency)


# ---------------------------------------------------------------------------
# 3. Targeting de macrofagos
# ---------------------------------------------------------------------------


def compute_macrophage_targeting(
    base_diameter_nm: float,
    mannose_peg_fraction: float = 0.5,
) -> MacrophageTargeting:
    """Modela o efeito de manose-PEG-DSPE no targeting de macrofagos.

    Macrofagos infectados por L. infantum apresentam expressao elevada
    do receptor de manose (CD206/MRC1). Este receptor medeia endocitose
    de particulas manose-decoradas, direcionando-as ao sistema
    endossomo-lisossomo-fagolisossomo.

    A substituicao parcial de DMG-PEG2000 por manose-PEG-DSPE
    funcionaliza a superficie da LNP sem alterar significativamente
    o tamanho ou a estabilidade coloidal.

    O fator de aumento de uptake e modelado como funcao linear da
    fracao de manose-PEG, interpolando entre os limites de literatura
    (2-5x para macrofagos M2 peritoneais).

    Args:
        base_diameter_nm: Diametro da LNP sem modificacao.
        mannose_peg_fraction: Fracao do PEG-lipidio substituida por
            manose-PEG (0 a 1). Default 0.5 (50%).

    Returns:
        MacrophageTargeting com metricas calculadas.
    """
    # Fator de uptake: interpolacao linear entre min e max
    uptake = (
        MANNOSE_UPTAKE_FACTOR_MIN
        + mannose_peg_fraction
        * (MANNOSE_UPTAKE_FACTOR_MAX - MANNOSE_UPTAKE_FACTOR_MIN)
    )

    # Manose-PEG-DSPE e ligeiramente maior que DMG-PEG2000
    # Efeito no diametro: aumento modesto (~2-5 nm) proporcional a fracao
    diameter_increase = mannose_peg_fraction * 4.0  # nm
    final_diameter = base_diameter_nm + diameter_increase

    suitable = MIN_PHAGOCYTIC_NM <= final_diameter <= MAX_PHAGOCYTIC_NM

    return MacrophageTargeting(
        mannose_peg_fraction=mannose_peg_fraction,
        uptake_fold_increase=round(uptake, 2),
        receptor_target="CD206 (mannose receptor / MRC1)",
        particle_diameter_nm=round(final_diameter, 1),
        suitable_for_phagocytosis=suitable,
    )


# ---------------------------------------------------------------------------
# 4. Estabilidade em armazenamento
# ---------------------------------------------------------------------------


def compute_storage_stability(temperature_celsius: float) -> StorageStability:
    """Modela estabilidade da LNP em armazenamento por Arrhenius.

    A taxa de agregacao/fusao de LNPs segue cinetica de Arrhenius:
        k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T))

    onde k_ref e a taxa a 25 C e Ea e a energia de ativacao.
    A meia-vida e: t1/2 = ln(2) / k(T)

    Considercao para Brasil: temperaturas ambiente de 25-35 C sao
    comuns em regioes endemicas de leishmaniose visceral (Nordeste).
    Liofilizacao elimina a dependencia de cadeia fria e e essencial
    para distribuicao em areas rurais.

    Liofilizacao e viavel para LNPs quando:
    - Crioprotetores (sacarose/trealose) sao adicionados
    - O diametro apos reconstituicao nao aumenta >20%
    - A encapsulacao permanece >80% apos reconstituicao
    Ref: Ball RL et al. (2017) Drug Deliv Transl Res 7(1):89-103

    Args:
        temperature_celsius: Temperatura de armazenamento (C).

    Returns:
        StorageStability com metricas calculadas.
    """
    t_kelvin = temperature_celsius + 273.15

    # Taxa de agregacao por Arrhenius
    k_agg = K_AGG_25C * math.exp(
        EA_AGGREGATION / R_GAS_KJ * (1.0 / T_REF_K - 1.0 / t_kelvin)
    )

    # Meia-vida (dias)
    half_life = math.log(2) / k_agg if k_agg > 0 else float("inf")

    # Fracao intacta em varios tempos (cinetica de primeira ordem)
    f_30 = math.exp(-k_agg * 30)
    f_90 = math.exp(-k_agg * 90)
    f_365 = math.exp(-k_agg * 365)

    # Liofilizacao: viavel para qualquer temperatura de armazenamento pos-liofilizacao
    # mas so faz sentido economicamente se armazenamento a frio nao e possivel
    # Para Brasil: ESSENCIAL (cadeia fria fragil em areas endemicas rurais)
    lyophilization_feasible = True  # LNPs podem ser liofilizadas com crioprotetores

    return StorageStability(
        temperature_celsius=temperature_celsius,
        temperature_kelvin=round(t_kelvin, 2),
        k_aggregation_per_day=round(k_agg, 8),
        half_life_days=round(half_life, 1),
        fraction_intact_30_days=round(f_30, 4),
        fraction_intact_90_days=round(f_90, 4),
        fraction_intact_365_days=round(f_365, 4),
        lyophilization_feasible=lyophilization_feasible,
    )


def compute_storage_profiles() -> list[StorageStability]:
    """Calcula estabilidade para tres temperaturas clinicamente relevantes.

    Returns:
        Lista com StorageStability para -20 C, 4 C, e 25 C.
    """
    temperatures = [-20.0, 4.0, 25.0]
    return [compute_storage_stability(t) for t in temperatures]
