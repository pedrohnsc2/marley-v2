"""A3 -- Formulacao de LNP (Lipid Nanoparticle) veterinaria para mRNA.

Modela a formulacao de nanoparticulas lipidicas para entrega do mRNA
vacinal em caes (Canis lupus familiaris), com foco em:

    - Dose: calculo por peso corporal (ug mRNA/kg) para cao de 15 kg
    - Composicao lipidica: ionizable lipid/DSPC/cholesterol/PEG-lipid
    - Razao N/P: otimizacao para encapsulamento de mRNA
    - Estabilidade: -70°C, -20°C, 2-8°C (liofilizado)
    - Custo: escala laboratorial vs industrial
    - Comparacao com vacinas veterinarias de mRNA existentes

As premissas de formulacao sao baseadas nas vacinas de mRNA-LNP
aprovadas para humanos (BNT162b2, mRNA-1273) adaptadas para uso
veterinario. Nao existem vacinas de mRNA aprovadas para caes ate
2025, mas ensaios pre-clinicos em veterinaria estao em andamento.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes de formulacao
# ---------------------------------------------------------------------------

# Peso corporal alvo (cao de porte medio, referencia para dose)
TARGET_BODY_WEIGHT_KG: Final[float] = 15.0

# Faixa de dose de mRNA por kg (baseada em vacinas COVID-19 humanas)
# BNT162b2: 30 ug/dose (~0.4 ug/kg para adulto de 70 kg)
# mRNA-1273: 100 ug/dose (~1.4 ug/kg)
# Para veterinaria, doses maiores sao tipicas: 10-100 ug/kg
DOSE_RANGE_UG_PER_KG: Final[tuple[float, float]] = (10.0, 100.0)
RECOMMENDED_DOSE_UG_PER_KG: Final[float] = 30.0  # ug/kg, valor recomendado

# Composicao lipidica padrao (molar %)
# Baseada na formulacao da Moderna (SM-102:DSPC:Chol:PEG-DMG)
LNP_COMPOSITION: Final[dict[str, float]] = {
    "ionizable_lipid_SM102": 50.0,   # mol%
    "DSPC": 10.0,                     # mol%
    "cholesterol": 38.5,              # mol%
    "PEG_DMG_2000": 1.5,             # mol%
}

# Pesos moleculares dos lipideos (g/mol)
LIPID_MW: Final[dict[str, float]] = {
    "ionizable_lipid_SM102": 586.8,   # SM-102
    "DSPC": 790.2,                     # 1,2-distearoyl-sn-glycero-3-phosphocholine
    "cholesterol": 386.7,              # Colesterol
    "PEG_DMG_2000": 2509.2,           # PEG-DMG 2000
}

# Faixa de razao N/P (nitrogeno do lipidio ionizavel : fosfato do RNA)
NP_RATIO_RANGE: Final[tuple[float, float]] = (6.0, 12.0)
OPTIMAL_NP_RATIO: Final[float] = 6.0  # BNT162b2 usa N/P ~ 6

# Massa molecular media de um nucleotideo de mRNA (com cap e modificacoes)
NT_AVG_MW: Final[float] = 330.0  # g/mol (media para A, G, C, psiU com backbone)

# Numero de cargas negativas por nucleotideo (um fosfato = uma carga)
PHOSPHATE_PER_NT: Final[int] = 1

# Numero de nitrogenos ionizaveis por molecula de SM-102
N_PER_IONIZABLE_LIPID: Final[int] = 1

# Encapsulamento tipico (% de mRNA encapsulado nas LNPs)
ENCAPSULATION_EFFICIENCY: Final[float] = 0.95  # 95% e tipico para mRNA-LNP

# Regime vacinal
DOSES_PER_ANIMAL: Final[int] = 2  # Primo + 1 reforco (tipico para mRNA)

# ---------------------------------------------------------------------------
# Custos de lipideos (USD por grama, escala laboratorial)
# Fonte: catologos AventiPharma, BroadPharm, MedChemExpress (2024-2025)
# ---------------------------------------------------------------------------

LIPID_COST_PER_G_LAB: Final[dict[str, float]] = {
    "ionizable_lipid_SM102": 2500.00,  # ~$2500/g laboratorial
    "DSPC": 150.00,                     # ~$150/g
    "cholesterol": 5.00,                # ~$5/g (commodity)
    "PEG_DMG_2000": 800.00,            # ~$800/g
}

# Fator de reducao para escala industrial (contrato com fabricante)
INDUSTRIAL_COST_FACTOR: Final[float] = 0.05  # 5% do custo lab

# Custo de IVT (transcricao in vitro) por mg de mRNA
IVT_COST_PER_MG_LAB: Final[float] = 500.00     # USD/mg escala lab
IVT_COST_PER_MG_INDUSTRIAL: Final[float] = 5.00  # USD/mg escala industrial

# Custos adicionais por dose
FORMULATION_COST_PER_DOSE_LAB: Final[float] = 15.00
FORMULATION_COST_PER_DOSE_INDUSTRIAL: Final[float] = 1.50
QC_COST_PER_DOSE_LAB: Final[float] = 20.00
QC_COST_PER_DOSE_INDUSTRIAL: Final[float] = 2.00

# ---------------------------------------------------------------------------
# Condicoes de armazenamento
# ---------------------------------------------------------------------------

STORAGE_CONDITIONS: Final[list[dict[str, str | float]]] = [
    {
        "label": "Ultra-frozen (-70 C)",
        "temperature_c": -70.0,
        "shelf_life_months": 6,
        "cold_chain": "Ultra-cold chain required (-60 to -80 C)",
        "logistics_note": (
            "Requer freezer -80 C; cadeia fria cara e complexa para "
            "areas rurais endemicas no Brasil (Nordeste, Norte)"
        ),
        "feasibility": "Laboratorial e centros urbanos",
    },
    {
        "label": "Frozen (-20 C)",
        "temperature_c": -20.0,
        "shelf_life_months": 3,
        "cold_chain": "Standard freezer (-15 to -25 C)",
        "logistics_note": (
            "Freezer domestico/comercial; viavel para clinicas "
            "veterinarias urbanas e hospitais veterinarios"
        ),
        "feasibility": "Ampla distribuicao urbana",
    },
    {
        "label": "Refrigerated (2-8 C, lyophilized)",
        "temperature_c": 5.0,
        "shelf_life_months": 12,
        "cold_chain": "Standard refrigerator (2-8 C)",
        "logistics_note": (
            "Liofilizacao permite armazenamento em geladeira padrao; "
            "reconstituicao antes do uso; melhor opcao para campo "
            "e campanhas de vacinacao em areas endemicas"
        ),
        "feasibility": "Distribuicao ampla incluindo areas rurais",
    },
]


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class DoseCalculation:
    """Calculo de dose para um animal alvo."""

    body_weight_kg: float
    dose_ug_per_kg: float
    total_dose_ug: float
    total_dose_nmol: float  # nanomoles de mRNA
    mrna_length_nt: int
    doses_per_animal: int
    total_mrna_per_animal_ug: float


@dataclass(frozen=True, slots=True)
class LNPComposition:
    """Composicao detalhada da LNP para uma dose."""

    molar_ratios: dict[str, float]
    np_ratio: float
    ionizable_lipid_nmol: float
    dspc_nmol: float
    cholesterol_nmol: float
    peg_lipid_nmol: float
    total_lipid_mass_ug: float
    lipid_to_mrna_ratio: float
    encapsulation_efficiency: float


@dataclass(frozen=True, slots=True)
class NPRatioPoint:
    """Ponto de avaliacao para uma razao N/P especifica."""

    np_ratio: float
    ionizable_lipid_nmol: float
    total_lipid_mass_ug: float
    estimated_encapsulation_pct: float
    estimated_transfection_score: float
    is_optimal: bool


@dataclass(frozen=True, slots=True)
class CostEstimate:
    """Estimativa de custo por dose."""

    scale: str  # "laboratory" ou "industrial"
    ivt_cost_usd: float
    lipid_cost_usd: float
    formulation_cost_usd: float
    qc_cost_usd: float
    total_per_dose_usd: float
    total_per_animal_usd: float
    lipid_breakdown: dict[str, float]


@dataclass
class LNPFormulation:
    """Resultado completo da formulacao de LNP."""

    dose: DoseCalculation
    composition: LNPComposition
    np_ratio_scan: list[NPRatioPoint]
    storage_options: list[dict]
    cost_lab: CostEstimate
    cost_industrial: CostEstimate
    veterinary_comparison: dict
    recommendation: str


# ---------------------------------------------------------------------------
# Funcoes de calculo
# ---------------------------------------------------------------------------


def _calculate_dose(mrna_length_nt: int) -> DoseCalculation:
    """Calcula a dose de mRNA para um cao de 15 kg.

    A dose e baseada em ug de mRNA por kg de peso corporal.
    O valor recomendado de 30 ug/kg esta entre BNT162b2 (~0.4 ug/kg)
    e mRNA-1273 (~1.4 ug/kg) escalonado para veterinaria, onde
    doses maiores sao comuns pela menor necessidade de minimizar
    efeitos adversos e pelo custo relativo menor por animal.

    Args:
        mrna_length_nt: comprimento do mRNA em nucleotideos

    Returns:
        DoseCalculation com todos os parametros
    """
    total_dose_ug = RECOMMENDED_DOSE_UG_PER_KG * TARGET_BODY_WEIGHT_KG

    # Converter ug para nmol usando peso molecular do mRNA
    mrna_mw = mrna_length_nt * NT_AVG_MW  # g/mol
    total_dose_g = total_dose_ug * 1e-6
    total_dose_mol = total_dose_g / mrna_mw
    total_dose_nmol = total_dose_mol * 1e9

    total_per_animal = total_dose_ug * DOSES_PER_ANIMAL

    return DoseCalculation(
        body_weight_kg=TARGET_BODY_WEIGHT_KG,
        dose_ug_per_kg=RECOMMENDED_DOSE_UG_PER_KG,
        total_dose_ug=total_dose_ug,
        total_dose_nmol=round(total_dose_nmol, 4),
        mrna_length_nt=mrna_length_nt,
        doses_per_animal=DOSES_PER_ANIMAL,
        total_mrna_per_animal_ug=total_per_animal,
    )


def _calculate_lnp_composition(dose: DoseCalculation) -> LNPComposition:
    """Calcula a composicao lipidica da LNP para a dose especificada.

    A razao N/P (nitrogeno ionizavel : fosfato do RNA) determina
    a quantidade de lipidio ionizavel necessaria. Os demais lipideos
    sao calculados pelas razoes molares fixas (50:10:38.5:1.5).

    Args:
        dose: calculo de dose com quantidade de mRNA

    Returns:
        LNPComposition detalhada
    """
    # Numero de fosfatos no mRNA (1 por nucleotideo)
    n_phosphates = dose.mrna_length_nt  # por molecula de mRNA

    # Nmol de fosfato total na dose
    phosphate_nmol = dose.total_dose_nmol * n_phosphates

    # Nmol de lipidio ionizavel necessario (N/P ratio)
    ionizable_nmol = phosphate_nmol * OPTIMAL_NP_RATIO / N_PER_IONIZABLE_LIPID

    # Calcular demais lipideos pelas razoes molares
    # SM-102 = 50 mol%, entao 1 mol% = ionizable_nmol / 50
    one_mol_pct = ionizable_nmol / LNP_COMPOSITION["ionizable_lipid_SM102"]
    dspc_nmol = one_mol_pct * LNP_COMPOSITION["DSPC"]
    chol_nmol = one_mol_pct * LNP_COMPOSITION["cholesterol"]
    peg_nmol = one_mol_pct * LNP_COMPOSITION["PEG_DMG_2000"]

    # Massa total de lipideos (em ug)
    total_lipid_ug = (
        ionizable_nmol * LIPID_MW["ionizable_lipid_SM102"] * 1e-3
        + dspc_nmol * LIPID_MW["DSPC"] * 1e-3
        + chol_nmol * LIPID_MW["cholesterol"] * 1e-3
        + peg_nmol * LIPID_MW["PEG_DMG_2000"] * 1e-3
    )

    lipid_to_mrna = total_lipid_ug / dose.total_dose_ug if dose.total_dose_ug > 0 else 0.0

    return LNPComposition(
        molar_ratios=dict(LNP_COMPOSITION),
        np_ratio=OPTIMAL_NP_RATIO,
        ionizable_lipid_nmol=round(ionizable_nmol, 2),
        dspc_nmol=round(dspc_nmol, 2),
        cholesterol_nmol=round(chol_nmol, 2),
        peg_lipid_nmol=round(peg_nmol, 2),
        total_lipid_mass_ug=round(total_lipid_ug, 2),
        lipid_to_mrna_ratio=round(lipid_to_mrna, 2),
        encapsulation_efficiency=ENCAPSULATION_EFFICIENCY,
    )


def _scan_np_ratios(dose: DoseCalculation) -> list[NPRatioPoint]:
    """Varre razoes N/P de 6 a 12 e estima encapsulamento e transfeccao.

    Para cada N/P ratio, calcula a quantidade de lipidio ionizavel
    e estima eficiencia de encapsulamento e score de transfeccao
    baseado em dados publicados de otimizacao de LNP.

    Modelo empirico:
    - Encapsulamento: maximo em N/P 6-8, decai levemente acima
    - Transfeccao: maximo em N/P 6-10, citotoxicidade acima de 10

    Args:
        dose: calculo de dose

    Returns:
        Lista de NPRatioPoint para cada N/P testado
    """
    points: list[NPRatioPoint] = []
    n_phosphates = dose.mrna_length_nt
    phosphate_nmol = dose.total_dose_nmol * n_phosphates

    for np_int in range(6, 13):
        np_ratio = float(np_int)
        ion_nmol = phosphate_nmol * np_ratio / N_PER_IONIZABLE_LIPID

        # Massa total de lipideos para esta razao N/P
        one_mol_pct = ion_nmol / LNP_COMPOSITION["ionizable_lipid_SM102"]
        total_lipid_ug = (
            ion_nmol * LIPID_MW["ionizable_lipid_SM102"] * 1e-3
            + one_mol_pct * LNP_COMPOSITION["DSPC"] * LIPID_MW["DSPC"] * 1e-3
            + one_mol_pct * LNP_COMPOSITION["cholesterol"] * LIPID_MW["cholesterol"] * 1e-3
            + one_mol_pct * LNP_COMPOSITION["PEG_DMG_2000"] * LIPID_MW["PEG_DMG_2000"] * 1e-3
        )

        # Modelos empiricos simplificados
        # Encapsulamento: ~95% em N/P 6, ~97% em N/P 8, ~93% em N/P 12
        encap = 0.95 + 0.02 * (1.0 - abs(np_ratio - 8.0) / 6.0)
        encap = min(0.98, max(0.88, encap))

        # Score de transfeccao: pico em N/P ~8, decai nas extremidades
        # e penalizado por citotoxicidade em N/P > 10
        base_transfection = 1.0 - 0.05 * (np_ratio - 8.0) ** 2 / 10.0
        cytotox_penalty = max(0.0, (np_ratio - 10.0) * 0.1)
        transfection = max(0.0, min(1.0, base_transfection - cytotox_penalty))

        is_optimal = np_ratio == OPTIMAL_NP_RATIO

        points.append(NPRatioPoint(
            np_ratio=np_ratio,
            ionizable_lipid_nmol=round(ion_nmol, 2),
            total_lipid_mass_ug=round(total_lipid_ug, 2),
            estimated_encapsulation_pct=round(encap * 100, 1),
            estimated_transfection_score=round(transfection, 3),
            is_optimal=is_optimal,
        ))

    return points


def _calculate_costs(dose: DoseCalculation) -> tuple[CostEstimate, CostEstimate]:
    """Calcula custos de producao em escala laboratorial e industrial.

    Componentes de custo:
    1. IVT (transcricao in vitro): producao do mRNA
    2. Lipideos: SM-102, DSPC, colesterol, PEG-DMG
    3. Formulacao: mistura, extrusion, purificacao
    4. QC: controle de qualidade (tamanho, encapsulamento, esterilidade)

    Args:
        dose: calculo de dose

    Returns:
        Tupla (custo_lab, custo_industrial)
    """
    dose_mg = dose.total_dose_ug / 1000.0

    # --- Custo de lipideos por dose ---
    # Calcular massa de cada lipideo necessaria (em gramas)
    composition = _calculate_lnp_composition(dose)

    lipid_masses_g = {
        "ionizable_lipid_SM102": composition.ionizable_lipid_nmol * LIPID_MW["ionizable_lipid_SM102"] * 1e-12,
        "DSPC": composition.dspc_nmol * LIPID_MW["DSPC"] * 1e-12,
        "cholesterol": composition.cholesterol_nmol * LIPID_MW["cholesterol"] * 1e-12,
        "PEG_DMG_2000": composition.peg_lipid_nmol * LIPID_MW["PEG_DMG_2000"] * 1e-12,
    }

    # Lab scale
    lipid_cost_lab = sum(
        mass * LIPID_COST_PER_G_LAB[lipid]
        for lipid, mass in lipid_masses_g.items()
    )
    lipid_breakdown_lab = {
        lipid: round(mass * LIPID_COST_PER_G_LAB[lipid], 4)
        for lipid, mass in lipid_masses_g.items()
    }

    ivt_lab = dose_mg * IVT_COST_PER_MG_LAB
    total_lab = ivt_lab + lipid_cost_lab + FORMULATION_COST_PER_DOSE_LAB + QC_COST_PER_DOSE_LAB

    cost_lab = CostEstimate(
        scale="laboratory",
        ivt_cost_usd=round(ivt_lab, 2),
        lipid_cost_usd=round(lipid_cost_lab, 4),
        formulation_cost_usd=FORMULATION_COST_PER_DOSE_LAB,
        qc_cost_usd=QC_COST_PER_DOSE_LAB,
        total_per_dose_usd=round(total_lab, 2),
        total_per_animal_usd=round(total_lab * DOSES_PER_ANIMAL, 2),
        lipid_breakdown=lipid_breakdown_lab,
    )

    # Industrial scale
    lipid_cost_ind = sum(
        mass * LIPID_COST_PER_G_LAB[lipid] * INDUSTRIAL_COST_FACTOR
        for lipid, mass in lipid_masses_g.items()
    )
    lipid_breakdown_ind = {
        lipid: round(mass * LIPID_COST_PER_G_LAB[lipid] * INDUSTRIAL_COST_FACTOR, 6)
        for lipid, mass in lipid_masses_g.items()
    }

    ivt_ind = dose_mg * IVT_COST_PER_MG_INDUSTRIAL
    total_ind = ivt_ind + lipid_cost_ind + FORMULATION_COST_PER_DOSE_INDUSTRIAL + QC_COST_PER_DOSE_INDUSTRIAL

    cost_ind = CostEstimate(
        scale="industrial",
        ivt_cost_usd=round(ivt_ind, 4),
        lipid_cost_usd=round(lipid_cost_ind, 6),
        formulation_cost_usd=FORMULATION_COST_PER_DOSE_INDUSTRIAL,
        qc_cost_usd=QC_COST_PER_DOSE_INDUSTRIAL,
        total_per_dose_usd=round(total_ind, 4),
        total_per_animal_usd=round(total_ind * DOSES_PER_ANIMAL, 4),
        lipid_breakdown=lipid_breakdown_ind,
    )

    return cost_lab, cost_ind


def _veterinary_mrna_comparison() -> dict:
    """Compara com vacinas veterinarias de mRNA existentes ou em desenvolvimento.

    Ate 2025, nao ha vacinas de mRNA aprovadas para caes. Existem
    vacinas de mRNA para suinos (PRRS) e equinos em fase experimental.

    Returns:
        Dicionario com comparacao estruturada
    """
    return {
        "approved_veterinary_mrna_vaccines": [],
        "note": (
            "Ate 2025, nenhuma vacina de mRNA foi aprovada para uso veterinario "
            "em caes. As vacinas de mRNA aprovadas sao exclusivamente para humanos "
            "(COVID-19: BNT162b2/Comirnaty, mRNA-1273/Spikevax)."
        ),
        "experimental_veterinary_mrna": [
            {
                "target": "PRRS virus (suinos)",
                "species": "Sus scrofa (suino)",
                "status": "Pre-clinico",
                "developer": "Varios grupos academicos",
                "reference": "Lutz et al., npj Vaccines, 2023",
            },
            {
                "target": "Influenza equina",
                "species": "Equus caballus (equino)",
                "status": "Fase I/II",
                "developer": "CureVac / colaboracao academica",
                "reference": "Schnee et al., PLoS Pathog, 2016",
            },
            {
                "target": "Raiva (prova de conceito)",
                "species": "Multiplas especies",
                "status": "Pre-clinico",
                "developer": "CureVac",
                "reference": "Alberer et al., Lancet, 2017",
            },
        ],
        "regulatory_pathway": (
            "No Brasil: MAPA (Ministerio da Agricultura) regulamenta biologicos "
            "veterinarios. Nao existe via regulatoria especifica para mRNA "
            "veterinario -- seria necessario processo de novo. "
            "Nos EUA: USDA-APHIS Center for Veterinary Biologics (CVB) aceita "
            "novas tecnologias via 9 CFR Part 113, com requerimento de "
            "demonstracao de seguranca, pureza e potencia."
        ),
        "advantage_over_protein": (
            "mRNA-LNP elimina necessidade de expressao recombinante, purificacao "
            "e refolding. Producao mais rapida e escalavel. Contudo, cold chain "
            "mais exigente e custo de NTPs modificados sao desvantagens para "
            "mercado veterinario (sensivel a preco)."
        ),
    }


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def run_lnp_formulation(mrna_length_nt: int) -> LNPFormulation:
    """Executa a analise completa de formulacao de LNP veterinaria.

    Etapas:
        1. Calcula dose para cao de 15 kg
        2. Formula composicao lipidica otimizada
        3. Varre razoes N/P de 6 a 12
        4. Avalia opcoes de armazenamento
        5. Estima custos em duas escalas
        6. Compara com vacinas veterinarias existentes

    Args:
        mrna_length_nt: comprimento do mRNA em nucleotideos

    Returns:
        LNPFormulation com resultados completos
    """
    # Passo 1: Dose
    dose = _calculate_dose(mrna_length_nt)

    # Passo 2: Composicao LNP
    composition = _calculate_lnp_composition(dose)

    # Passo 3: Scan de N/P ratios
    np_scan = _scan_np_ratios(dose)

    # Passo 4: Custos
    cost_lab, cost_ind = _calculate_costs(dose)

    # Passo 5: Comparacao veterinaria
    vet_comparison = _veterinary_mrna_comparison()

    # Passo 6: Recomendacao
    recommendation = (
        f"Para um cao de {TARGET_BODY_WEIGHT_KG:.0f} kg, recomenda-se "
        f"{dose.total_dose_ug:.0f} ug de mRNA por dose "
        f"({RECOMMENDED_DOSE_UG_PER_KG:.0f} ug/kg), formulado em LNP com "
        f"composicao SM-102:DSPC:Chol:PEG-DMG (50:10:38.5:1.5 mol%) "
        f"e razao N/P = {OPTIMAL_NP_RATIO:.0f}. "
        f"Regime de {DOSES_PER_ANIMAL} doses (primo + reforco). "
        f"Para distribuicao em areas endemicas brasileiras, a formulacao "
        f"liofilizada (2-8 C, 12 meses de validade) e a mais adequada, "
        f"apesar do custo adicional de liofilizacao. "
        f"Custo estimado por animal: ${cost_lab.total_per_animal_usd:.2f} (lab) / "
        f"${cost_ind.total_per_animal_usd:.2f} (industrial)."
    )

    return LNPFormulation(
        dose=dose,
        composition=composition,
        np_ratio_scan=np_scan,
        storage_options=list(STORAGE_CONDITIONS),
        cost_lab=cost_lab,
        cost_industrial=cost_ind,
        veterinary_comparison=vet_comparison,
        recommendation=recommendation,
    )


def to_dict(formulation: LNPFormulation) -> dict:
    """Serializa a formulacao de LNP para JSON.

    Args:
        formulation: resultado da formulacao

    Returns:
        Dicionario serializavel para JSON
    """
    d = formulation.dose
    c = formulation.composition

    return {
        "dose": {
            "target_species": "Canis lupus familiaris",
            "body_weight_kg": d.body_weight_kg,
            "dose_ug_per_kg": d.dose_ug_per_kg,
            "dose_range_ug_per_kg": list(DOSE_RANGE_UG_PER_KG),
            "total_dose_ug": d.total_dose_ug,
            "total_dose_nmol": d.total_dose_nmol,
            "mrna_length_nt": d.mrna_length_nt,
            "doses_per_animal": d.doses_per_animal,
            "total_mrna_per_animal_ug": d.total_mrna_per_animal_ug,
        },
        "lnp_composition": {
            "molar_ratios_pct": c.molar_ratios,
            "np_ratio": c.np_ratio,
            "ionizable_lipid_nmol": c.ionizable_lipid_nmol,
            "dspc_nmol": c.dspc_nmol,
            "cholesterol_nmol": c.cholesterol_nmol,
            "peg_lipid_nmol": c.peg_lipid_nmol,
            "total_lipid_mass_ug": c.total_lipid_mass_ug,
            "lipid_to_mrna_ratio_w_w": c.lipid_to_mrna_ratio,
            "encapsulation_efficiency_pct": c.encapsulation_efficiency * 100,
        },
        "np_ratio_optimization": [
            {
                "np_ratio": p.np_ratio,
                "ionizable_lipid_nmol": p.ionizable_lipid_nmol,
                "total_lipid_mass_ug": p.total_lipid_mass_ug,
                "estimated_encapsulation_pct": p.estimated_encapsulation_pct,
                "estimated_transfection_score": p.estimated_transfection_score,
                "is_optimal": p.is_optimal,
            }
            for p in formulation.np_ratio_scan
        ],
        "storage_options": formulation.storage_options,
        "cost_laboratory": {
            "ivt_usd": formulation.cost_lab.ivt_cost_usd,
            "lipids_usd": formulation.cost_lab.lipid_cost_usd,
            "formulation_usd": formulation.cost_lab.formulation_cost_usd,
            "qc_usd": formulation.cost_lab.qc_cost_usd,
            "total_per_dose_usd": formulation.cost_lab.total_per_dose_usd,
            "total_per_animal_usd": formulation.cost_lab.total_per_animal_usd,
            "lipid_breakdown": formulation.cost_lab.lipid_breakdown,
        },
        "cost_industrial": {
            "ivt_usd": formulation.cost_industrial.ivt_cost_usd,
            "lipids_usd": formulation.cost_industrial.lipid_cost_usd,
            "formulation_usd": formulation.cost_industrial.formulation_cost_usd,
            "qc_usd": formulation.cost_industrial.qc_cost_usd,
            "total_per_dose_usd": formulation.cost_industrial.total_per_dose_usd,
            "total_per_animal_usd": formulation.cost_industrial.total_per_animal_usd,
            "lipid_breakdown": formulation.cost_industrial.lipid_breakdown,
        },
        "veterinary_mrna_comparison": formulation.veterinary_comparison,
        "recommendation": formulation.recommendation,
    }
