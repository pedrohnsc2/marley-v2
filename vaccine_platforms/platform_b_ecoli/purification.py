"""B3 -- Modelo de protocolo de purificacao in silico.

Simula o pipeline de purificacao cromatografica de uma proteina
recombinante His6-tagged expressa em E. coli BL21(DE3) / pET-28a(+).

Pipeline:
    1. IMAC (Ni-NTA) -- captura por afinidade ao His6-tag
    2. Clivagem TEV -- remocao do His6-tag (gera Ser N-terminal)
    3. IMAC reverso -- proteina clivada no flow-through
    4. SEC (Size Exclusion Chromatography) -- polimento final

Os rendimentos e purezas sao estimativas baseadas em literatura e
experiencia pratica com proteinas recombinantes de 25-50 kDa em E. coli.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Constantes do modelo
# ---------------------------------------------------------------------------

# Rendimento tipico de expressao em pET-28a / BL21(DE3) com IPTG
# Para proteinas soluveis de 25-50 kDa, 20-50 mg/L e conservador
EXPRESSION_YIELD_LOW: Final[float] = 20.0   # mg/L
EXPRESSION_YIELD_HIGH: Final[float] = 50.0  # mg/L
EXPRESSION_YIELD_TYPICAL: Final[float] = 35.0  # mg/L (ponto medio)

# Condicoes de expressao recomendadas
EXPRESSION_STRAIN: Final[str] = "BL21(DE3)"
EXPRESSION_VECTOR: Final[str] = "pET-28a(+)"
IPTG_CONCENTRATION_MM: Final[float] = 0.5   # mM
INDUCTION_TEMP_C: Final[float] = 18.0       # °C (baixa para solubilidade)
INDUCTION_TIME_H: Final[float] = 16.0       # horas (overnight a 18°C)
OD600_INDUCTION: Final[float] = 0.6         # Fase log media

# Eficiencias de cada etapa (fracao do material que entra)
# Baseado em: Structural Genomics Consortium benchmarks
IMAC_CAPTURE_EFFICIENCY: Final[float] = 0.90   # 90% de captura
IMAC_CAPTURE_PURITY: Final[float] = 0.80       # 80% puro apos IMAC
TEV_CLEAVAGE_EFFICIENCY: Final[float] = 0.95   # 95% de clivagem
TEV_RECOVERY: Final[float] = 0.90              # 90% recuperado (perda em diluicao)
REVERSE_IMAC_EFFICIENCY: Final[float] = 0.85   # 85% no flow-through
REVERSE_IMAC_PURITY: Final[float] = 0.95       # 95% puro apos reverso
SEC_RECOVERY: Final[float] = 0.80              # 80% recuperado
SEC_PURITY: Final[float] = 0.98                # >95%, tipicamente 98%

# TEV protease
TEV_RECOGNITION: Final[str] = "ENLYFQ/S"  # Cliva entre Q e S
TEV_RATIO: Final[str] = "1:50"            # w/w TEV:substrato
TEV_INCUBATION_TEMP: Final[float] = 4.0   # °C (overnight a 4°C)
TEV_INCUBATION_TIME_H: Final[float] = 16.0

# SEC
SEC_COLUMN: Final[str] = "Superdex 75 10/300 GL"
SEC_BUFFER: Final[str] = "20 mM Tris-HCl pH 8.0, 150 mM NaCl"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class PurificationStep:
    """Uma etapa individual do protocolo de purificacao."""

    step_number: int
    name: str
    method: str
    input_mg: float
    output_mg: float
    efficiency: float
    purity_after: float
    buffer: str
    conditions: str
    notes: str = ""


@dataclass
class PurificationProtocol:
    """Protocolo completo de purificacao com rendimento e custo."""

    steps: list[PurificationStep]
    expression_system: str
    expression_vector: str
    culture_volume_l: float
    initial_yield_mg: float
    final_yield_mg: float
    final_purity_pct: float
    overall_recovery_pct: float
    expected_yield_range: tuple[float, float]


# ---------------------------------------------------------------------------
# Modelo de purificacao
# ---------------------------------------------------------------------------


def model_purification(
    molecular_weight_da: float,
    culture_volume_l: float = 1.0,
) -> PurificationProtocol:
    """Modela o protocolo de purificacao para a proteina recombinante.

    O modelo assume expressao soluvel em E. coli BL21(DE3) com inducao
    a baixa temperatura (18°C overnight). A proteina e purificada por:

    1. Lise celular (sonicacao ou prensa francesa)
    2. IMAC Ni-NTA (His6-tag captura)
    3. Clivagem TEV (remocao do tag)
    4. IMAC reverso (separacao tag/proteina)
    5. SEC (polimento e troca de tampao)

    Args:
        molecular_weight_da: Peso molecular da proteina em Daltons.
        culture_volume_l: Volume de cultura em litros.

    Returns:
        PurificationProtocol com todas as etapas e metricas.
    """
    # Producao inicial estimada
    initial_yield = EXPRESSION_YIELD_TYPICAL * culture_volume_l
    yield_low = EXPRESSION_YIELD_LOW * culture_volume_l
    yield_high = EXPRESSION_YIELD_HIGH * culture_volume_l

    steps: list[PurificationStep] = []
    current_amount = initial_yield

    # ----- Etapa 1: Lise celular -----
    # Nao e uma etapa cromatografica mas e necessaria
    lysis_recovery = 0.95  # 95% de recuperacao na lise
    amount_after_lysis = current_amount * lysis_recovery

    steps.append(PurificationStep(
        step_number=0,
        name="Lise Celular",
        method="Sonicacao (6x 30s pulsos, 50% amplitude, gelo)",
        input_mg=current_amount,
        output_mg=amount_after_lysis,
        efficiency=lysis_recovery,
        purity_after=0.10,  # ~10% da proteina total e o recombinante
        buffer="50 mM Tris-HCl pH 8.0, 300 mM NaCl, 10 mM imidazol, 1 mM PMSF",
        conditions="4°C, centrifugacao 20.000g 30 min",
        notes="Lisozima 1 mg/mL pre-incubada 30 min; DNase I 10 ug/mL",
    ))
    current_amount = amount_after_lysis

    # ----- Etapa 2: IMAC Ni-NTA (captura) -----
    amount_after_imac = current_amount * IMAC_CAPTURE_EFFICIENCY

    steps.append(PurificationStep(
        step_number=1,
        name="IMAC Ni-NTA (Captura)",
        method="Cromatografia de afinidade a metal imobilizado",
        input_mg=current_amount,
        output_mg=amount_after_imac,
        efficiency=IMAC_CAPTURE_EFFICIENCY,
        purity_after=IMAC_CAPTURE_PURITY,
        buffer=(
            "Lavagem: 50 mM Tris pH 8.0, 300 mM NaCl, 20 mM imidazol; "
            "Eluicao: gradiente 50-250 mM imidazol"
        ),
        conditions="4°C, fluxo 1 mL/min, coluna HisTrap FF 5mL (GE Healthcare)",
        notes=(
            "His6-tag se coordena com Ni2+ via ligacao de coordenacao; "
            "eluicao competitiva com imidazol; dialise overnight contra "
            "tampao TEV"
        ),
    ))
    current_amount = amount_after_imac

    # ----- Etapa 3: Clivagem TEV -----
    amount_after_tev = current_amount * TEV_CLEAVAGE_EFFICIENCY * TEV_RECOVERY

    steps.append(PurificationStep(
        step_number=2,
        name="Clivagem TEV",
        method="Protease TEV (Tobacco Etch Virus NIa)",
        input_mg=current_amount,
        output_mg=amount_after_tev,
        efficiency=TEV_CLEAVAGE_EFFICIENCY * TEV_RECOVERY,
        purity_after=IMAC_CAPTURE_PURITY,  # Purity unchanged by cleavage
        buffer="50 mM Tris pH 8.0, 150 mM NaCl, 1 mM DTT, 0.5 mM EDTA",
        conditions=(
            f"{TEV_INCUBATION_TEMP}°C, {TEV_INCUBATION_TIME_H}h; "
            f"razao TEV:substrato {TEV_RATIO} (w/w)"
        ),
        notes=(
            f"Reconhecimento: {TEV_RECOGNITION}; apos clivagem resta "
            "Ser no N-terminal da proteina alvo; TEV e His6-tagged "
            "(sera removida no passo seguinte)"
        ),
    ))
    current_amount = amount_after_tev

    # ----- Etapa 4: IMAC reverso -----
    # Proteina clivada (sem His6) passa direto; TEV e fragmentos His6 ficam
    amount_after_reverse = current_amount * REVERSE_IMAC_EFFICIENCY

    steps.append(PurificationStep(
        step_number=3,
        name="IMAC Reverso (Subtractivo)",
        method="Re-passagem em Ni-NTA",
        input_mg=current_amount,
        output_mg=amount_after_reverse,
        efficiency=REVERSE_IMAC_EFFICIENCY,
        purity_after=REVERSE_IMAC_PURITY,
        buffer="50 mM Tris pH 8.0, 300 mM NaCl, 20 mM imidazol",
        conditions="4°C, coleta do flow-through",
        notes=(
            "Proteina clivada nao tem His6-tag, entao nao liga na resina; "
            "TEV protease (His-tagged) e proteina nao-clivada ficam retidas"
        ),
    ))
    current_amount = amount_after_reverse

    # ----- Etapa 5: SEC (polimento) -----
    amount_after_sec = current_amount * SEC_RECOVERY

    steps.append(PurificationStep(
        step_number=4,
        name="SEC (Polimento)",
        method="Cromatografia de Exclusao Molecular",
        input_mg=current_amount,
        output_mg=amount_after_sec,
        efficiency=SEC_RECOVERY,
        purity_after=SEC_PURITY,
        buffer=SEC_BUFFER,
        conditions=f"4°C, coluna {SEC_COLUMN}, fluxo 0.5 mL/min",
        notes=(
            f"MW esperado: {molecular_weight_da / 1000:.1f} kDa (monomero); "
            "verificar por SDS-PAGE se ha oligomeros; "
            "fracoes centrais do pico simetrico"
        ),
    ))
    current_amount = amount_after_sec

    # Calcular metricas globais
    overall_recovery = current_amount / initial_yield if initial_yield > 0 else 0.0

    return PurificationProtocol(
        steps=steps,
        expression_system=EXPRESSION_STRAIN,
        expression_vector=EXPRESSION_VECTOR,
        culture_volume_l=culture_volume_l,
        initial_yield_mg=round(initial_yield, 2),
        final_yield_mg=round(current_amount, 2),
        final_purity_pct=round(SEC_PURITY * 100, 1),
        overall_recovery_pct=round(overall_recovery * 100, 1),
        expected_yield_range=(
            round(yield_low * overall_recovery, 2),
            round(yield_high * overall_recovery, 2),
        ),
    )


def to_dict(protocol: PurificationProtocol) -> dict:
    """Serializa o protocolo de purificacao para JSON.

    Args:
        protocol: Protocolo modelado.

    Returns:
        Dicionario serializavel para JSON.
    """
    steps_list = []
    for step in protocol.steps:
        steps_list.append({
            "step_number": step.step_number,
            "name": step.name,
            "method": step.method,
            "input_mg": round(step.input_mg, 2),
            "output_mg": round(step.output_mg, 2),
            "efficiency_pct": round(step.efficiency * 100, 1),
            "purity_after_pct": round(step.purity_after * 100, 1),
            "buffer": step.buffer,
            "conditions": step.conditions,
            "notes": step.notes,
        })

    return {
        "expression_system": protocol.expression_system,
        "expression_vector": protocol.expression_vector,
        "expression_conditions": {
            "strain": EXPRESSION_STRAIN,
            "vector": EXPRESSION_VECTOR,
            "iptg_mm": IPTG_CONCENTRATION_MM,
            "induction_temp_c": INDUCTION_TEMP_C,
            "induction_time_h": INDUCTION_TIME_H,
            "od600_induction": OD600_INDUCTION,
        },
        "culture_volume_l": protocol.culture_volume_l,
        "initial_yield_mg": protocol.initial_yield_mg,
        "final_yield_mg": protocol.final_yield_mg,
        "final_purity_pct": protocol.final_purity_pct,
        "overall_recovery_pct": protocol.overall_recovery_pct,
        "expected_yield_range_mg": list(protocol.expected_yield_range),
        "steps": steps_list,
    }
