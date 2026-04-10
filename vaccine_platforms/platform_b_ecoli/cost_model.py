"""B4 -- Modelo de custos para producao da vacina recombinante em E. coli.

Estima custos de producao em escala laboratorial e industrial para uma
vacina proteica recombinante expressa em E. coli BL21(DE3).

Premissas:
- Dose: 50 ug de proteina purificada + adjuvante QuilA
- Regime: 3 doses por animal (primo-vacinacao + 2 reforcos)
- Cold chain: 2-8°C (refrigerador padrao, sem ultracongelamento)
- Escala industrial reduz custos em ~10x (economia de escala)

Custos baseados em catalogos Sigma-Aldrich, GE Healthcare e fornecedores
brasileiros (2024-2025).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

# ---------------------------------------------------------------------------
# Constantes de custo (USD)
# ---------------------------------------------------------------------------

# Fermentacao (custo por litro de cultura)
LB_MEDIUM_PER_L: Final[float] = 5.00       # Meio LB autoclavado
ANTIBIOTICS_PER_L: Final[float] = 2.00     # Kanamicina (pET-28a)
IPTG_PER_L: Final[float] = 3.00            # 0.5 mM final
CONSUMABLES_PER_L: Final[float] = 5.00     # Frascos, agitacao, etc
FERMENTATION_COST_PER_L: Final[float] = (
    LB_MEDIUM_PER_L + ANTIBIOTICS_PER_L + IPTG_PER_L + CONSUMABLES_PER_L
)  # ~$15/L

# Purificacao
NINTA_RESIN_COST: Final[float] = 200.00    # 25 mL resina Ni-NTA
NINTA_RESIN_VOLUME_ML: Final[float] = 25.0
NINTA_REUSES: Final[int] = 5               # Reutilizavel ate 5 vezes
NINTA_COST_PER_USE: Final[float] = NINTA_RESIN_COST / NINTA_REUSES

TEV_PROTEASE_COST: Final[float] = 50.00    # $50 por 1000 unidades
TEV_UNITS_PER_VIAL: Final[float] = 1000.0
TEV_MG_PER_VIAL: Final[float] = 1.0        # 1 vial cliva ~1 mg substrato
# Nota: TEV comercial e caro (~$50/mg substrato). Labs academicos
# frequentemente produzem TEV in-house a custo negligivel.
# Aqui usamos o custo comercial como cenario conservador.

SEC_COLUMN_COST: Final[float] = 2000.00    # Superdex 75 10/300 GL
SEC_COLUMN_USES: Final[int] = 100
SEC_COST_PER_USE: Final[float] = SEC_COLUMN_COST / SEC_COLUMN_USES

# Buffers e consumiveis de purificacao
PURIFICATION_BUFFERS_PER_RUN: Final[float] = 10.00
PURIFICATION_CONSUMABLES_PER_RUN: Final[float] = 15.00

# Adjuvante e formulacao
QUILA_COST_PER_DOSE: Final[float] = 1.50   # QuilA saponina por dose
FORMULATION_COST_PER_DOSE: Final[float] = 0.50  # Excipientes, pH, filtracao

# Dose por animal
DOSE_UG: Final[float] = 50.0               # ug de proteina por dose
DOSES_PER_ANIMAL: Final[int] = 3           # Primo + 2 reforcos

# Fator de escala industrial
# Fermentadores de 10L -> 100L -> 1000L reduzem custo significativamente
# A economia de escala em bioprocesso tipicamente e 0.1x a 0.15x do lab
INDUSTRIAL_SCALE_FACTOR: Final[float] = 0.10

# Cold chain
COLD_CHAIN_TYPE: Final[str] = "2-8°C"
COLD_CHAIN_NOTE: Final[str] = (
    "Proteina recombinante liofilizada ou em solucao; "
    "refrigerador padrao, sem necessidade de ultracongelamento (-80°C); "
    "vantagem significativa sobre vacinas de mRNA (que requerem -20°C a -80°C)"
)


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass
class CostModel:
    """Modelo de custos completo para a vacina recombinante."""

    # Custos laboratoriais
    fermentation_cost_per_l: float
    purification_cost_per_run: float
    total_lab_cost_per_l: float
    yield_mg_per_l: float
    lab_cost_per_mg: float
    lab_cost_per_dose: float
    lab_cost_per_animal: float

    # Custos industriais
    industrial_cost_per_mg: float
    industrial_cost_per_dose: float
    industrial_cost_per_animal: float

    # Adjuvante e formulacao
    adjuvant_cost_per_dose: float
    formulation_cost_per_dose: float
    total_dose_cost_lab: float
    total_dose_cost_industrial: float

    # Metadados
    dose_ug: float
    doses_per_animal: int
    cold_chain: str
    cold_chain_note: str
    scale_factor: float

    # Detalhamento
    cost_breakdown: dict


# ---------------------------------------------------------------------------
# Calculo
# ---------------------------------------------------------------------------


def calculate_costs(
    yield_mg_per_l: float,
    dose_ug: float = DOSE_UG,
) -> CostModel:
    """Calcula os custos de producao em escala laboratorial e industrial.

    O modelo considera:
    - Custo de fermentacao por litro de cultura
    - Custo de purificacao por batch (IMAC + TEV + reverso + SEC)
    - Rendimento final em mg/L apos purificacao
    - Dose por animal (ug)
    - Fator de escala industrial (0.1x)

    Args:
        yield_mg_per_l: Rendimento final apos purificacao (mg/L cultura).
        dose_ug: Microgramas de proteina por dose.

    Returns:
        CostModel com todos os custos calculados.
    """
    # Custo de purificacao por batch (1L)
    # TEV: custo proporcional a quantidade de proteina produzida
    tev_cost = (yield_mg_per_l / TEV_MG_PER_VIAL) * TEV_PROTEASE_COST

    purification_cost = (
        NINTA_COST_PER_USE
        + tev_cost
        + SEC_COST_PER_USE
        + PURIFICATION_BUFFERS_PER_RUN
        + PURIFICATION_CONSUMABLES_PER_RUN
    )

    # Custo total lab por litro de cultura
    total_lab_per_l = FERMENTATION_COST_PER_L + purification_cost

    # Custo por mg (escala lab)
    cost_per_mg_lab = total_lab_per_l / yield_mg_per_l if yield_mg_per_l > 0 else 0.0

    # Custo por dose (escala lab) -- dose em ug, converter para mg
    dose_mg = dose_ug / 1000.0
    cost_per_dose_lab = cost_per_mg_lab * dose_mg

    # Custo por animal (3 doses)
    cost_per_animal_lab = cost_per_dose_lab * DOSES_PER_ANIMAL

    # Custos industriais (fator de escala)
    cost_per_mg_industrial = cost_per_mg_lab * INDUSTRIAL_SCALE_FACTOR
    cost_per_dose_industrial = cost_per_dose_lab * INDUSTRIAL_SCALE_FACTOR
    cost_per_animal_industrial = cost_per_animal_lab * INDUSTRIAL_SCALE_FACTOR

    # Custo total por dose incluindo adjuvante e formulacao
    total_dose_lab = cost_per_dose_lab + QUILA_COST_PER_DOSE + FORMULATION_COST_PER_DOSE
    total_dose_industrial = cost_per_dose_industrial + QUILA_COST_PER_DOSE + FORMULATION_COST_PER_DOSE

    # Detalhamento
    breakdown = {
        "fermentation": {
            "lb_medium": LB_MEDIUM_PER_L,
            "antibiotics_kanamicin": ANTIBIOTICS_PER_L,
            "iptg_05mm": IPTG_PER_L,
            "consumables": CONSUMABLES_PER_L,
            "subtotal_per_l": FERMENTATION_COST_PER_L,
        },
        "purification": {
            "ni_nta_resin_per_use": round(NINTA_COST_PER_USE, 2),
            "tev_protease": round(tev_cost, 2),
            "sec_column_per_use": SEC_COST_PER_USE,
            "buffers": PURIFICATION_BUFFERS_PER_RUN,
            "consumables": PURIFICATION_CONSUMABLES_PER_RUN,
            "subtotal_per_run": round(purification_cost, 2),
        },
        "adjuvant_and_formulation": {
            "quila_saponin_per_dose": QUILA_COST_PER_DOSE,
            "formulation_per_dose": FORMULATION_COST_PER_DOSE,
            "subtotal_per_dose": QUILA_COST_PER_DOSE + FORMULATION_COST_PER_DOSE,
        },
        "scale_comparison": {
            "lab_1l_culture": {
                "total_cost_usd": round(total_lab_per_l, 2),
                "yield_mg": round(yield_mg_per_l, 2),
                "doses_produced": round((yield_mg_per_l * 1000) / dose_ug, 0),
                "cost_per_dose_usd": round(total_dose_lab, 4),
            },
            "industrial_1000l_culture": {
                "estimated_total_cost_usd": round(
                    total_lab_per_l * 1000 * INDUSTRIAL_SCALE_FACTOR, 2
                ),
                "estimated_yield_mg": round(yield_mg_per_l * 1000, 0),
                "doses_produced": round(
                    (yield_mg_per_l * 1000 * 1000) / dose_ug, 0
                ),
                "cost_per_dose_usd": round(total_dose_industrial, 4),
            },
        },
    }

    return CostModel(
        fermentation_cost_per_l=FERMENTATION_COST_PER_L,
        purification_cost_per_run=round(purification_cost, 2),
        total_lab_cost_per_l=round(total_lab_per_l, 2),
        yield_mg_per_l=round(yield_mg_per_l, 2),
        lab_cost_per_mg=round(cost_per_mg_lab, 4),
        lab_cost_per_dose=round(cost_per_dose_lab, 4),
        lab_cost_per_animal=round(cost_per_animal_lab, 4),
        industrial_cost_per_mg=round(cost_per_mg_industrial, 4),
        industrial_cost_per_dose=round(cost_per_dose_industrial, 4),
        industrial_cost_per_animal=round(cost_per_animal_industrial, 4),
        adjuvant_cost_per_dose=QUILA_COST_PER_DOSE,
        formulation_cost_per_dose=FORMULATION_COST_PER_DOSE,
        total_dose_cost_lab=round(total_dose_lab, 4),
        total_dose_cost_industrial=round(total_dose_industrial, 4),
        dose_ug=dose_ug,
        doses_per_animal=DOSES_PER_ANIMAL,
        cold_chain=COLD_CHAIN_TYPE,
        cold_chain_note=COLD_CHAIN_NOTE,
        scale_factor=INDUSTRIAL_SCALE_FACTOR,
        cost_breakdown=breakdown,
    )


def to_dict(model: CostModel) -> dict:
    """Serializa o modelo de custos para JSON.

    Args:
        model: Modelo de custos calculado.

    Returns:
        Dicionario serializavel para JSON.
    """
    return {
        "lab_scale": {
            "fermentation_cost_per_l": model.fermentation_cost_per_l,
            "purification_cost_per_run": model.purification_cost_per_run,
            "total_cost_per_l": model.total_lab_cost_per_l,
            "yield_mg_per_l": model.yield_mg_per_l,
            "cost_per_mg_usd": model.lab_cost_per_mg,
            "cost_per_dose_usd": model.lab_cost_per_dose,
            "cost_per_animal_usd": model.lab_cost_per_animal,
        },
        "industrial_scale": {
            "scale_factor": model.scale_factor,
            "cost_per_mg_usd": model.industrial_cost_per_mg,
            "cost_per_dose_usd": model.industrial_cost_per_dose,
            "cost_per_animal_usd": model.industrial_cost_per_animal,
        },
        "per_dose_total": {
            "protein_lab_usd": model.lab_cost_per_dose,
            "protein_industrial_usd": model.industrial_cost_per_dose,
            "adjuvant_quila_usd": model.adjuvant_cost_per_dose,
            "formulation_usd": model.formulation_cost_per_dose,
            "total_lab_usd": model.total_dose_cost_lab,
            "total_industrial_usd": model.total_dose_cost_industrial,
        },
        "dosing": {
            "dose_ug": model.dose_ug,
            "doses_per_animal": model.doses_per_animal,
        },
        "cold_chain": model.cold_chain,
        "cold_chain_note": model.cold_chain_note,
        "cost_breakdown": model.cost_breakdown,
    }
