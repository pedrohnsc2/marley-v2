"""C4 -- Modelo de custos para producao da vacina em L. tarentolae.

Estima custos de producao para uma vacina baseada em L. tarentolae
recombinante, tanto na forma viva (organismo inteiro como adjuvante
intrinseco) quanto na forma de proteina secretada purificada.

Vantagens de custo sobre outras plataformas:
    - Crescimento a 26C (sem incubadora de CO2)
    - Meio BHI relativamente barato (~$8/L)
    - BSL-1 (sem infraestrutura de biosseguranca especial)
    - Adjuvante intrinseco (moleculas do parasita ativam TLR2/4)
    - Integracao cromossomal = linhagem estavel (um unico evento de
      transfeccao gera cultura perpetua)

Duas modalidades de uso:
    A. Vacina viva (organismo inteiro) -- mais barato, com adjuvante intrinseco
    B. Proteina secretada purificada -- mais controle de qualidade, precisa
       de adjuvante exogeno (QuilA)

Custos baseados em catalogos Sigma-Aldrich, Jena Bioscience e
fornecedores brasileiros (2024-2025).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes de custo (USD)
# ---------------------------------------------------------------------------

# Meio de crescimento (BHI suplementado)
BHI_MEDIUM_PER_L: Final[float] = 8.00        # Brain Heart Infusion
HEMIN_PER_L: Final[float] = 2.00            # Suplemento de hemina
PENICILLIN_PER_L: Final[float] = 1.50       # Antibiotico padrao
NOURSEOTHRICIN_PER_L: Final[float] = 3.00   # Selecao sat2
CONSUMABLES_PER_L: Final[float] = 3.50      # Frascos, agitacao, etc
GROWTH_COST_PER_L: Final[float] = (
    BHI_MEDIUM_PER_L + HEMIN_PER_L + PENICILLIN_PER_L
    + NOURSEOTHRICIN_PER_L + CONSUMABLES_PER_L
)  # ~$18/L

# Equipamento (sem necessidade de CO2 ou temperatura especial)
# Incubadora padrao a 26C -- custo eletrico negligivel vs 37C com CO2
EQUIPMENT_NOTE: Final[str] = (
    "Crescimento a 26C em frascos de cultura padrao (Erlenmeyer ou T-flask). "
    "Sem necessidade de incubadora de CO2 ou equipamento especializado. "
    "Agitacao orbital a 120 rpm e suficiente."
)

# pLEXSY kit (custo unico de setup)
PLEXSY_KIT_COST: Final[float] = 1500.00     # Kit completo Jena Bioscience
PLEXSY_KIT_NOTE: Final[str] = (
    "Kit pLEXSY-sat2 (Jena Bioscience) inclui vetor, cepa host, "
    "reagentes de transfeccao e nourseotricina para selecao. "
    "Custo unico que gera linhagem estavel para producao perpetua."
)
# Numero de batches estimados a partir de um kit (linhagem e permanente)
PLEXSY_BATCHES_FROM_KIT: Final[int] = 100
PLEXSY_AMORTIZED_PER_BATCH: Final[float] = PLEXSY_KIT_COST / PLEXSY_BATCHES_FROM_KIT

# Purificacao (para modalidade de proteina secretada)
# His6 + Strep-Tactin tandem
NI_NTA_RESIN_COST: Final[float] = 200.00     # 25 mL resina
NI_NTA_REUSES: Final[int] = 5
NI_NTA_COST_PER_USE: Final[float] = NI_NTA_RESIN_COST / NI_NTA_REUSES

STREP_TACTIN_RESIN_COST: Final[float] = 350.00  # 5 mL resina Strep-Tactin XT
STREP_TACTIN_REUSES: Final[int] = 10
STREP_TACTIN_COST_PER_USE: Final[float] = STREP_TACTIN_RESIN_COST / STREP_TACTIN_REUSES

PURIFICATION_BUFFERS: Final[float] = 10.00
PURIFICATION_CONSUMABLES: Final[float] = 10.00

# Rendimento estimado de proteina secretada por litro de cultura
# L. tarentolae secreta proteinas recombinantes a ~1-10 mg/L (tipico)
# Com otimizacao: 5-20 mg/L (secrecao SAP1 otimizada)
ESTIMATED_YIELD_LOW_MG_L: Final[float] = 5.0
ESTIMATED_YIELD_HIGH_MG_L: Final[float] = 20.0
ESTIMATED_YIELD_TYPICAL_MG_L: Final[float] = 10.0
PURIFICATION_RECOVERY: Final[float] = 0.60   # 60% recuperacao tandem

# Adjuvante (so necessario para proteina purificada; vivo tem intrinseco)
QUILA_COST_PER_DOSE: Final[float] = 1.50    # QuilA saponina por dose
FORMULATION_COST_PER_DOSE: Final[float] = 0.50

# Dose
DOSE_UG_PROTEIN: Final[float] = 50.0        # ug de proteina por dose
DOSE_LIVE_CFU: Final[float] = 1e7           # 10^7 parasitas por dose (viva)
DOSES_PER_ANIMAL: Final[int] = 3            # Primo + 2 reforcos

# Cold chain
COLD_CHAIN_LIVE: Final[str] = "2-8C (cultura viva) ou liofilizado (temp ambiente)"
COLD_CHAIN_PROTEIN: Final[str] = "2-8C (proteina purificada)"

# Fator de escala industrial
INDUSTRIAL_SCALE_FACTOR: Final[float] = 0.12  # Um pouco maior que E. coli (parasita mais complexo)


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass
class TarentolaCostModel:
    """Modelo de custos completo para a vacina baseada em L. tarentolae."""

    # Custos comuns
    growth_cost_per_l: float
    plexsy_kit_cost: float
    plexsy_amortized_per_batch: float

    # Modalidade A: vacina viva (organismo inteiro)
    live_cost_per_l: float
    live_doses_per_l: float
    live_cost_per_dose: float
    live_cost_per_animal: float

    # Modalidade B: proteina secretada purificada
    purification_cost_per_run: float
    protein_yield_mg_per_l: float
    protein_cost_per_mg: float
    protein_cost_per_dose: float
    protein_total_dose_cost: float    # Inclui adjuvante
    protein_cost_per_animal: float

    # Industrial
    live_industrial_cost_per_dose: float
    protein_industrial_cost_per_dose: float
    protein_industrial_total_dose: float

    # Metadados
    dose_ug_protein: float
    dose_live_cfu: float
    doses_per_animal: int
    cold_chain_live: str
    cold_chain_protein: str
    industrial_scale_factor: float

    # Detalhamento
    cost_breakdown: dict
    comparison_vs_platforms: dict


# ---------------------------------------------------------------------------
# Calculo
# ---------------------------------------------------------------------------


def calculate_costs(
    protein_yield_mg_per_l: float = ESTIMATED_YIELD_TYPICAL_MG_L,
) -> TarentolaCostModel:
    """Calcula os custos de producao para ambas as modalidades.

    Modalidade A (vacina viva):
        Custo = meio + selecao + consumiveis por litro de cultura
        1L a ~5 x 10^8 celulas/mL = 5 x 10^11 celulas total
        Dose = 10^7 celulas -> ~50,000 doses por litro

    Modalidade B (proteina secretada):
        Custo = crescimento + purificacao (His6/Strep tandem)
        Rendimento tipico: 10 mg/L (apos purificacao: 6 mg/L)
        Dose = 50 ug -> ~120 doses por litro (purificado)

    Args:
        protein_yield_mg_per_l: Rendimento bruto de proteina (mg/L).

    Returns:
        TarentolaCostModel com todos os custos.
    """
    # --- Modalidade A: vacina viva ---
    # Custo por litro de cultura
    live_cost_per_l = GROWTH_COST_PER_L + PLEXSY_AMORTIZED_PER_BATCH

    # Doses por litro: ~5e8 celulas/mL * 1000 mL / 1e7 celulas/dose
    cells_per_l = 5e8 * 1000  # 5 x 10^11 celulas por litro
    live_doses_per_l = cells_per_l / DOSE_LIVE_CFU  # ~50,000 doses

    live_cost_per_dose = live_cost_per_l / live_doses_per_l
    live_cost_per_animal = live_cost_per_dose * DOSES_PER_ANIMAL

    # --- Modalidade B: proteina secretada ---
    purification_cost = (
        NI_NTA_COST_PER_USE
        + STREP_TACTIN_COST_PER_USE
        + PURIFICATION_BUFFERS
        + PURIFICATION_CONSUMABLES
    )

    # Rendimento apos purificacao
    purified_yield_mg_l = protein_yield_mg_per_l * PURIFICATION_RECOVERY
    protein_total_cost_per_l = GROWTH_COST_PER_L + PLEXSY_AMORTIZED_PER_BATCH + purification_cost

    protein_cost_per_mg = (
        protein_total_cost_per_l / purified_yield_mg_l
        if purified_yield_mg_l > 0 else 0.0
    )
    dose_mg = DOSE_UG_PROTEIN / 1000.0
    protein_cost_per_dose = protein_cost_per_mg * dose_mg
    protein_total_dose_cost = protein_cost_per_dose + QUILA_COST_PER_DOSE + FORMULATION_COST_PER_DOSE
    protein_cost_per_animal = protein_total_dose_cost * DOSES_PER_ANIMAL

    # --- Escala industrial ---
    live_industrial = live_cost_per_dose * INDUSTRIAL_SCALE_FACTOR
    protein_industrial_dose = protein_cost_per_dose * INDUSTRIAL_SCALE_FACTOR
    protein_industrial_total = (
        protein_industrial_dose + QUILA_COST_PER_DOSE + FORMULATION_COST_PER_DOSE
    )

    # --- Comparacao com plataformas A e B ---
    comparison = {
        "platform_a_mrna": {
            "note": "mRNA sintetico: alto custo de sintese e cold chain (-20 a -80C)",
            "estimated_per_dose_usd": "5.00-15.00",
            "cold_chain": "-20C a -80C",
            "adjuvant": "LNP (nanoparticula lipidica) incluida",
        },
        "platform_b_ecoli": {
            "note": "Recombinante E. coli: bom custo, sem PTMs eucarioticos",
            "estimated_per_dose_usd": "2.00-4.00",
            "cold_chain": "2-8C",
            "adjuvant": "QuilA ($1.50/dose)",
        },
        "platform_c_tarentolae_live": {
            "note": "L. tarentolae viva: custo minimo, adjuvante intrinseco",
            "cost_per_dose_usd": round(live_cost_per_dose, 6),
            "cold_chain": COLD_CHAIN_LIVE,
            "adjuvant": "Intrinseco (moleculas do parasita)",
        },
        "platform_c_tarentolae_protein": {
            "note": "L. tarentolae secretada: PTMs corretas, custo intermediario",
            "cost_per_dose_usd": round(protein_total_dose_cost, 4),
            "cold_chain": COLD_CHAIN_PROTEIN,
            "adjuvant": "QuilA ($1.50/dose)",
        },
    }

    # --- Detalhamento ---
    breakdown = {
        "growth_medium": {
            "bhi_per_l": BHI_MEDIUM_PER_L,
            "hemin_per_l": HEMIN_PER_L,
            "penicillin_per_l": PENICILLIN_PER_L,
            "nourseothricin_per_l": NOURSEOTHRICIN_PER_L,
            "consumables_per_l": CONSUMABLES_PER_L,
            "subtotal_per_l": GROWTH_COST_PER_L,
        },
        "plexsy_kit": {
            "one_time_cost": PLEXSY_KIT_COST,
            "estimated_batches": PLEXSY_BATCHES_FROM_KIT,
            "amortized_per_batch": round(PLEXSY_AMORTIZED_PER_BATCH, 2),
            "note": PLEXSY_KIT_NOTE,
        },
        "purification_tandem": {
            "ni_nta_per_use": round(NI_NTA_COST_PER_USE, 2),
            "strep_tactin_per_use": round(STREP_TACTIN_COST_PER_USE, 2),
            "buffers": PURIFICATION_BUFFERS,
            "consumables": PURIFICATION_CONSUMABLES,
            "subtotal": round(purification_cost, 2),
            "recovery_pct": PURIFICATION_RECOVERY * 100,
        },
        "live_vaccine": {
            "cost_per_l": round(live_cost_per_l, 2),
            "cells_per_l": f"{cells_per_l:.0e}",
            "dose_cfu": f"{DOSE_LIVE_CFU:.0e}",
            "doses_per_l": round(live_doses_per_l, 0),
            "cost_per_dose": round(live_cost_per_dose, 6),
            "no_adjuvant_needed": True,
        },
        "protein_vaccine": {
            "yield_crude_mg_per_l": protein_yield_mg_per_l,
            "yield_purified_mg_per_l": round(purified_yield_mg_l, 2),
            "dose_ug": DOSE_UG_PROTEIN,
            "doses_per_l": round((purified_yield_mg_l * 1000) / DOSE_UG_PROTEIN, 0),
            "protein_cost_per_dose": round(protein_cost_per_dose, 4),
            "adjuvant_per_dose": QUILA_COST_PER_DOSE,
            "formulation_per_dose": FORMULATION_COST_PER_DOSE,
            "total_cost_per_dose": round(protein_total_dose_cost, 4),
        },
        "equipment": EQUIPMENT_NOTE,
        "bsl_level": "BSL-1 (custo de infraestrutura minimo)",
    }

    return TarentolaCostModel(
        growth_cost_per_l=GROWTH_COST_PER_L,
        plexsy_kit_cost=PLEXSY_KIT_COST,
        plexsy_amortized_per_batch=round(PLEXSY_AMORTIZED_PER_BATCH, 2),
        live_cost_per_l=round(live_cost_per_l, 2),
        live_doses_per_l=round(live_doses_per_l, 0),
        live_cost_per_dose=round(live_cost_per_dose, 6),
        live_cost_per_animal=round(live_cost_per_animal, 6),
        purification_cost_per_run=round(purification_cost, 2),
        protein_yield_mg_per_l=round(purified_yield_mg_l, 2),
        protein_cost_per_mg=round(protein_cost_per_mg, 4),
        protein_cost_per_dose=round(protein_cost_per_dose, 4),
        protein_total_dose_cost=round(protein_total_dose_cost, 4),
        protein_cost_per_animal=round(protein_cost_per_animal, 4),
        live_industrial_cost_per_dose=round(live_industrial, 6),
        protein_industrial_cost_per_dose=round(protein_industrial_dose, 4),
        protein_industrial_total_dose=round(protein_industrial_total, 4),
        dose_ug_protein=DOSE_UG_PROTEIN,
        dose_live_cfu=DOSE_LIVE_CFU,
        doses_per_animal=DOSES_PER_ANIMAL,
        cold_chain_live=COLD_CHAIN_LIVE,
        cold_chain_protein=COLD_CHAIN_PROTEIN,
        industrial_scale_factor=INDUSTRIAL_SCALE_FACTOR,
        cost_breakdown=breakdown,
        comparison_vs_platforms=comparison,
    )


def to_dict(model: TarentolaCostModel) -> dict:
    """Serializa o modelo de custos para JSON.

    Args:
        model: Modelo de custos calculado.

    Returns:
        Dicionario serializavel para JSON.
    """
    return {
        "live_vaccine": {
            "growth_cost_per_l": model.growth_cost_per_l,
            "doses_per_l": model.live_doses_per_l,
            "cost_per_dose_usd": model.live_cost_per_dose,
            "cost_per_animal_usd": model.live_cost_per_animal,
            "industrial_cost_per_dose_usd": model.live_industrial_cost_per_dose,
            "cold_chain": model.cold_chain_live,
            "adjuvant": "intrinsic (parasite molecules)",
        },
        "protein_vaccine": {
            "growth_cost_per_l": model.growth_cost_per_l,
            "purification_cost_per_run": model.purification_cost_per_run,
            "yield_purified_mg_per_l": model.protein_yield_mg_per_l,
            "cost_per_mg_usd": model.protein_cost_per_mg,
            "cost_per_dose_protein_usd": model.protein_cost_per_dose,
            "cost_per_dose_total_usd": model.protein_total_dose_cost,
            "cost_per_animal_usd": model.protein_cost_per_animal,
            "industrial_cost_per_dose_usd": model.protein_industrial_cost_per_dose,
            "industrial_total_dose_usd": model.protein_industrial_total_dose,
            "cold_chain": model.cold_chain_protein,
            "adjuvant": f"QuilA (${QUILA_COST_PER_DOSE:.2f}/dose)",
        },
        "setup": {
            "plexsy_kit_cost": model.plexsy_kit_cost,
            "plexsy_amortized_per_batch": model.plexsy_amortized_per_batch,
            "bsl_level": "BSL-1",
            "growth_temperature": "26C",
        },
        "dosing": {
            "dose_ug_protein": model.dose_ug_protein,
            "dose_live_cfu": f"{model.dose_live_cfu:.0e}",
            "doses_per_animal": model.doses_per_animal,
        },
        "breakdown": model.cost_breakdown,
        "platform_comparison": model.comparison_vs_platforms,
    }
