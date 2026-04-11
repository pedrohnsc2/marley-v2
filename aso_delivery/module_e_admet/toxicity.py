"""Avaliacao de toxicidade de PS-ASOs em Canis lupus familiaris.

Modela os efeitos adversos de classe conhecidos de ASOs com backbone
fosforotioato e modificacoes LNA, com foco em seguranca veterinaria.

NOTA SOBRE HONESTIDADE: este modulo reporta TODOS os riscos conhecidos
de PS-ASOs sem minimiza-los. PS-ASOs tem um perfil de seguranca bem
documentado em humanos (>10 aprovacoes FDA), mas efeitos adversos
reais existem e devem ser monitorados em uso veterinario.

Efeitos de classe de PS-ASOs (literatura consolidada):
1. Reacoes no local de injecao (inflamacao SC local)
2. Trombocitopenia (ativacao de complemento por PS)
3. Hepatotoxicidade (acumulacao hepatica de PS-ASOs)
4. Nefrotoxicidade (acumulacao tubular proximal)
5. Efeitos pro-inflamatorios (CpG -> TLR9 em caes)
6. Coagulopatia (interacao PS com fatores de coagulacao)

CASO ESPECIAL - TLR9: a ativacao de TLR9 por motivos CpG e normalmente
listada como efeito adverso, mas para MRL-ASO-001 e um MECANISMO
TERAPEUTICO INTENCIONAL (funcao dual). O desafio e controlar a dose
para manter o efeito imunosstimulatório na faixa terapeutica.

Referencias:
- Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563
- Henry SP et al. (2008) Toxicology 252(1-3):97-105
- Frazier KS (2015) Toxicol Pathol 43(1):78-89
- Inotersen label (FDA 2018) — thrombocytopenia warning
- Volpi S et al. (2012) J Immunol 188(12):5890-5897 — TLR9 in dogs
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final


# ---------------------------------------------------------------------------
# Constantes de seguranca
# ---------------------------------------------------------------------------

# NOAEL (No Observed Adverse Effect Level) em caes para PS-ASOs
# Ref: Henry SP et al. (2008) Toxicology 252(1-3):97-105
# Estudos de 13 semanas em caes com mipomersen-like PS-ASOs
NOAEL_MG_KG_WEEK: Final[float] = 40.0

# Dose maxima tolerada (MTD) em caes
# Doses acima de ~100 mg/kg/semana produzem toxicidade significativa
MTD_MG_KG_WEEK: Final[float] = 100.0

# Dose terapeutica planejada
THERAPEUTIC_DOSE_MG_KG_WEEK: Final[float] = 5.0


# ---------------------------------------------------------------------------
# Dataclasses de resultado
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ToxicityEndpoint:
    """Avaliacao de um endpoint toxico especifico.

    Cada endpoint tem:
    - Mecanismo pelo qual ocorre
    - Incidencia esperada na dose terapeutica
    - Severidade estimada
    - Estrategia de monitoramento
    - Se e reversivel apos descontinuacao
    """

    name: str
    mechanism: str
    incidence_at_therapeutic_dose: str
    severity: str                    # mild, moderate, severe
    dose_dependent: bool
    reversible: bool
    monitoring_strategy: str
    clinical_significance: str
    risk_score: float                # 0-10 (0 = sem risco, 10 = risco grave)


@dataclass(frozen=True)
class TLR9Assessment:
    """Avaliacao especifica do eixo TLR9 — dual function.

    Para MRL-ASO-001, a ativacao de TLR9 e simultaneamente:
    - Efeito adverso potencial (inflamacao sistemica se excessiva)
    - Mecanismo terapeutico desejado (imunoestimulacao anti-Leishmania)

    O desafio e manter a ativacao na faixa terapeutica.
    """

    cpg_motifs_in_sequence: int
    tlr9_activation_expected: bool
    therapeutic_benefit: str
    toxicity_risk: str
    dose_window_description: str
    monitoring_markers: list[str]
    dual_function_assessment: str


@dataclass(frozen=True)
class OrganToxicityProfile:
    """Perfil de toxicidade por orgao.

    PS-ASOs acumulam em orgaos especificos (figado > rim > baco),
    portanto o perfil toxico e orgao-dependente.
    """

    organ: str
    accumulation_factor: float       # Kp (tecido:plasma)
    primary_cell_type_affected: str
    mechanism_of_toxicity: str
    biomarkers: list[str]
    threshold_dose_mg_kg: float      # dose onde toxicidade comeca
    expected_at_5_mg_kg: str
    reversibility: str


@dataclass(frozen=True)
class SafetyProfile:
    """Perfil de seguranca completo do ASO.

    Integra todos os endpoints toxicos, perfis por orgao,
    avaliacao de TLR9, e recomendacoes de monitoramento.
    """

    therapeutic_index: float
    noael_mg_kg_week: float
    mtd_mg_kg_week: float
    dose_mg_kg_week: float
    safety_margin: float
    endpoints: list[ToxicityEndpoint]
    organ_profiles: list[OrganToxicityProfile]
    tlr9_assessment: TLR9Assessment
    overall_risk_classification: str
    monitoring_schedule: dict[str, str]
    contraindications: list[str]
    drug_interactions: list[str]


# ===========================================================================
# Funcoes de avaliacao de toxicidade
# ===========================================================================


def count_cpg_motifs(sequence: str) -> int:
    """Conta motivos CpG na sequencia do ASO.

    Motivos CpG (dinucleotideo 5'-CG-3') sao reconhecidos por TLR9
    no endossomo de celulas imunes. Em caes, TLR9 e especialmente
    responsivo a CpG nao-metilados (como em ASOs sinteticos).

    Para MRL-ASO-001: CpG e uma feature, nao um bug — mas precisa
    ser quantificado para prever o grau de ativacao imune.

    Args:
        sequence: Sequencia do ASO (DNA).

    Returns:
        Numero de dinucleotideos CpG na sequencia.
    """
    seq = sequence.upper()
    count = 0
    for i in range(len(seq) - 1):
        if seq[i] == "C" and seq[i + 1] == "G":
            count += 1
    return count


def assess_tlr9(aso_sequence: str) -> TLR9Assessment:
    """Avalia o impacto da ativacao de TLR9 pelo ASO.

    TLR9 reconhece DNA nao-metilado com motivos CpG — exatamente
    o que PS-ASOs sinteticos apresentam. Em caes:

    - TLR9 esta em celulas B, DCs plasmacitoides, macrofagos
    - Ativacao produz IFN-alpha, IL-6, TNF-alpha
    - Em doses baixas: imunoestimulacao benefica (anti-Leishmania)
    - Em doses altas: inflamacao sistemica (citokine storm)

    Para MRL-ASO-001, a ativacao de TLR9 complementa o knockdown
    do SL RNA: enquanto o ASO degrada o RNA parasitario via RNase H,
    o backbone PS + CpG ativa a imunidade inata do hospedeiro.

    Args:
        aso_sequence: Sequencia do ASO.

    Returns:
        Avaliacao completa do eixo TLR9.
    """
    n_cpg = count_cpg_motifs(aso_sequence)

    return TLR9Assessment(
        cpg_motifs_in_sequence=n_cpg,
        tlr9_activation_expected=n_cpg > 0,
        therapeutic_benefit=(
            "TLR9 activation by PS-ASO CpG motifs induces IFN-alpha and "
            "pro-inflammatory cytokines from pDCs and macrophages. This enhances "
            "anti-Leishmania immunity by: (1) activating macrophage killing of "
            "intracellular amastigotes (iNOS upregulation, NO production), "
            "(2) promoting Th1 polarization (critical for leishmaniasis resolution), "
            "(3) enhancing antigen presentation. This dual function (ASO + adjuvant) "
            "is UNIQUE to MRL-ASO-001 among proposed leishmaniasis treatments."
        ),
        toxicity_risk=(
            f"With {n_cpg} CpG motif(s), moderate TLR9 activation is expected. "
            "At therapeutic dose (5 mg/kg/week), activation should remain in the "
            "immunostimulatory range without reaching cytokine storm levels. "
            "RISK: at loading dose (10 mg/kg 2x/week), transient fever, "
            "lymphadenopathy, and elevated acute phase proteins are possible. "
            "These are EXPECTED pharmacological effects, not idiosyncratic toxicity."
        ),
        dose_window_description=(
            "Therapeutic window for TLR9 activation:\n"
            "  - Below 1 mg/kg/week: minimal immune activation (subtherapeutic)\n"
            "  - 2-10 mg/kg/week: optimal immunostimulation (target range)\n"
            "  - 10-40 mg/kg/week: elevated cytokines, manageable inflammation\n"
            "  - Above 40 mg/kg/week: risk of systemic inflammatory response\n"
            "Proposed dose (5 mg/kg/week) is in the middle of the optimal range."
        ),
        monitoring_markers=[
            "Serum IFN-alpha (TLR9 activation marker)",
            "IL-6 and TNF-alpha (inflammatory cytokines)",
            "C-reactive protein (acute phase response)",
            "Body temperature (transient post-injection fever)",
            "Lymph node size (reactive lymphadenopathy)",
        ],
        dual_function_assessment=(
            f"MRL-ASO-001 has {n_cpg} CpG motif(s) in a 25-nt PS backbone. "
            "The combination of RNase H-mediated SL RNA degradation AND "
            "TLR9-mediated innate immune activation creates a dual mechanism "
            "of action: (1) direct parasite RNA knockdown, (2) host immune "
            "potentiation against Leishmania. No other proposed leishmaniasis "
            "therapy combines target-specific RNA degradation with built-in "
            "immune adjuvant activity. The PS backbone, normally considered a "
            "liability for TLR9 activation, becomes a therapeutic ASSET in this "
            "infectious disease context."
        ),
    )


def assess_organ_toxicity() -> list[OrganToxicityProfile]:
    """Avalia toxicidade orgao-especifica do PS-ASO.

    PS-ASOs acumulam em orgaos com endotelio fenestrado:
    - Figado: Kp ~30x (celulas de Kupffer, hepatocitos)
    - Rim: Kp ~25x (tubulo proximal — reabsorcao)
    - Baco: Kp ~15x (zona marginal)

    Para cada orgao, avaliamos o mecanismo de toxicidade,
    biomarkers de monitoramento, e limiar de dose.

    Returns:
        Lista de perfis toxicologicos por orgao.
    """
    profiles: list[OrganToxicityProfile] = []

    # Figado — maior acumulacao, maior preocupacao
    profiles.append(OrganToxicityProfile(
        organ="liver",
        accumulation_factor=30.0,
        primary_cell_type_affected="hepatocytes and Kupffer cells",
        mechanism_of_toxicity=(
            "PS-ASO accumulation in hepatocytes via receptor-mediated endocytosis "
            "(stabilin-2, ASGPR). At high concentrations, ASO aggregates in "
            "lysosomes can cause basophilic granulation and mild transaminase "
            "elevation. Kupffer cell activation may contribute to inflammatory "
            "signaling. NOTE: hepatic accumulation is THERAPEUTICALLY BENEFICIAL "
            "for leishmaniasis (liver is a primary site of L. infantum infection)."
        ),
        biomarkers=["ALT", "AST", "ALP", "GGT", "total bilirubin"],
        threshold_dose_mg_kg=20.0,
        expected_at_5_mg_kg=(
            "Minimal hepatotoxicity expected. Possible mild ALT elevation "
            "(1.5-2x ULN) during loading phase, self-resolving. Hepatic "
            "accumulation at this dose is within the range seen in clinical "
            "PS-ASOs (mipomersen, inotersen) without dose-limiting toxicity."
        ),
        reversibility="Fully reversible upon dose reduction or discontinuation",
    ))

    # Rim — segundo orgao de acumulacao
    profiles.append(OrganToxicityProfile(
        organ="kidney",
        accumulation_factor=25.0,
        primary_cell_type_affected="proximal tubule epithelial cells",
        mechanism_of_toxicity=(
            "PS-ASOs are filtered glomerularly (small metabolites) and reabsorbed "
            "by proximal tubule cells via megalin/cubilin receptors. Intracellular "
            "accumulation can cause tubular basophilic granulation. At high doses "
            "(>20 mg/kg/week), tubular degeneration and proteinuria may occur. "
            "Inotersen (FDA label) carries a black box warning for glomerulonephritis, "
            "though this is antibody-mediated and rare."
        ),
        biomarkers=[
            "serum creatinine", "BUN", "urine protein:creatinine ratio",
            "urinalysis (sediment)", "cystatin C",
        ],
        threshold_dose_mg_kg=25.0,
        expected_at_5_mg_kg=(
            "Low nephrotoxicity risk at therapeutic dose. Possible minimal "
            "tubular basophilic granulation (histological finding without "
            "clinical significance). No clinically relevant creatinine "
            "elevation expected. Monitor urine protein:creatinine ratio "
            "monthly during treatment."
        ),
        reversibility="Reversible upon discontinuation; tubular regeneration in 2-4 weeks",
    ))

    # Baco
    profiles.append(OrganToxicityProfile(
        organ="spleen",
        accumulation_factor=15.0,
        primary_cell_type_affected="marginal zone macrophages and red pulp",
        mechanism_of_toxicity=(
            "Splenic accumulation of PS-ASOs activates marginal zone macrophages "
            "and can cause mild splenomegaly at high doses. TLR9 activation in "
            "splenic B cells contributes to inflammatory signaling. "
            "NOTE: for leishmaniasis, splenic accumulation is THERAPEUTICALLY "
            "BENEFICIAL — the spleen is a major reservoir of L. infantum "
            "amastigotes in canine visceral leishmaniasis."
        ),
        biomarkers=["spleen size (ultrasound)", "platelet count", "lymphocyte subsets"],
        threshold_dose_mg_kg=30.0,
        expected_at_5_mg_kg=(
            "Minimal splenic toxicity. Therapeutic accumulation in splenic "
            "macrophages is expected and desired for anti-leishmanial activity. "
            "Monitor spleen size by ultrasound at baseline and monthly."
        ),
        reversibility="Fully reversible upon discontinuation",
    ))

    # Medula ossea
    profiles.append(OrganToxicityProfile(
        organ="bone_marrow",
        accumulation_factor=5.0,
        primary_cell_type_affected="megakaryocytes and myeloid progenitors",
        mechanism_of_toxicity=(
            "PS-ASOs at high doses inhibit megakaryopoiesis, reducing platelet "
            "production. Mechanism: PS backbone interacts with growth factor "
            "receptors on megakaryocyte progenitors. Inotersen (FDA) carries "
            "a black box warning for thrombocytopenia (platelet count <100k in "
            "~24% of patients at 300 mg/week). Mipomersen: ~5% incidence. "
            "HONEST ASSESSMENT: thrombocytopenia is the most clinically "
            "significant class effect of PS-ASOs."
        ),
        biomarkers=[
            "platelet count (CBC)", "mean platelet volume",
            "immature platelet fraction", "bone marrow biopsy (if severe)",
        ],
        threshold_dose_mg_kg=10.0,
        expected_at_5_mg_kg=(
            "Mild thrombocytopenia possible (platelet count 100-150k, vs "
            "normal canine range 200-500k). Clinically significant "
            "thrombocytopenia (<100k) unlikely at 5 mg/kg/week but cannot "
            "be excluded. MANDATORY: monitor platelet count weekly during "
            "loading phase, biweekly during maintenance. Dose reduction "
            "or hold if platelets <100k."
        ),
        reversibility=(
            "Reversible upon dose reduction in 1-2 weeks. If severe (<50k), "
            "discontinue and monitor daily. Recovery expected in 5-7 days."
        ),
    ))

    return profiles


def assess_class_effects(aso_sequence: str) -> list[ToxicityEndpoint]:
    """Avalia efeitos adversos de classe de PS-ASOs.

    Lista completa e HONESTA dos efeitos conhecidos:
    - Nao minimiza riscos
    - Inclui dados de labels FDA (mipomersen, inotersen)
    - Reporta incidencias reais da literatura

    Args:
        aso_sequence: Sequencia do ASO (para contar CpG).

    Returns:
        Lista de endpoints toxicos avaliados.
    """
    n_cpg = count_cpg_motifs(aso_sequence)
    endpoints: list[ToxicityEndpoint] = []

    # 1. Reacoes no local de injecao
    endpoints.append(ToxicityEndpoint(
        name="Injection site reactions (ISR)",
        mechanism=(
            "Local inflammatory response to SC depot of PS-ASO. PS backbone "
            "activates innate immune cells at injection site. Erythema, "
            "induration, pruritus, and pain are common."
        ),
        incidence_at_therapeutic_dose="70-90% of patients (mipomersen: 84%, inotersen: 72%)",
        severity="mild",
        dose_dependent=True,
        reversible=True,
        monitoring_strategy=(
            "Visual inspection of injection site. Rotate injection sites "
            "(left/right lateral thorax, dorsal neck). Apply cold compress "
            "post-injection if needed."
        ),
        clinical_significance=(
            "Most common adverse effect but rarely dose-limiting. Self-resolving "
            "in 2-5 days. Tolerance develops with repeated dosing in most subjects."
        ),
        risk_score=3.0,
    ))

    # 2. Trombocitopenia
    endpoints.append(ToxicityEndpoint(
        name="Thrombocytopenia",
        mechanism=(
            "PS backbone interaction with megakaryocyte growth factor receptors "
            "reduces platelet production. May also involve complement activation "
            "(C3a, C5a) leading to platelet aggregation and consumption. "
            "Inotersen: black box warning for thrombocytopenia."
        ),
        incidence_at_therapeutic_dose=(
            "5-25% depending on dose. Mipomersen: ~5% at 200 mg/week. "
            "Inotersen: 24% had platelet <100k at 300 mg/week. "
            "At canine dose of 5 mg/kg/week (~75 mg/week for 15 kg dog): "
            "estimated 5-10% risk of clinically significant thrombocytopenia."
        ),
        severity="moderate",
        dose_dependent=True,
        reversible=True,
        monitoring_strategy=(
            "MANDATORY: CBC with platelet count weekly during loading phase, "
            "biweekly during maintenance. Hold dose if platelets <100k/uL. "
            "Discontinue if platelets <50k/uL. Do not restart until >150k."
        ),
        clinical_significance=(
            "Most clinically significant class effect. Can be dose-limiting. "
            "In the context of canine leishmaniasis (where thrombocytopenia is "
            "already common due to disease), careful baseline assessment and "
            "monitoring is ESSENTIAL. Pre-existing thrombocytopenia from "
            "leishmaniasis may increase risk."
        ),
        risk_score=6.0,
    ))

    # 3. Hepatotoxicidade
    endpoints.append(ToxicityEndpoint(
        name="Hepatotoxicity",
        mechanism=(
            "Hepatocyte accumulation of PS-ASO causes lysosomal overload "
            "and basophilic granulation. At high concentrations, inflammatory "
            "signaling (NF-kB activation) leads to transaminase elevation. "
            "Steatosis reported with mipomersen (related to ApoB knockdown, "
            "not applicable to MRL-ASO-001)."
        ),
        incidence_at_therapeutic_dose=(
            "Transaminase elevation (>3x ULN): 10-15% (mipomersen data). "
            "At lower canine dose (5 mg/kg/week): estimated <5% risk of "
            "clinically significant transaminase elevation."
        ),
        severity="mild to moderate",
        dose_dependent=True,
        reversible=True,
        monitoring_strategy=(
            "Serum ALT, AST at baseline, weekly during loading, monthly "
            "during maintenance. Dose reduction if ALT >5x ULN. "
            "Discontinue if ALT >10x ULN or if accompanied by clinical signs."
        ),
        clinical_significance=(
            "Manageable with monitoring. Hepatic accumulation is actually "
            "desired for anti-leishmanial activity (liver is primary site of "
            "L. infantum infection). The trade-off is favorable: hepatic ASO "
            "exposure treats the infection but requires liver enzyme monitoring."
        ),
        risk_score=4.0,
    ))

    # 4. Nefrotoxicidade
    endpoints.append(ToxicityEndpoint(
        name="Nephrotoxicity",
        mechanism=(
            "Proximal tubule accumulation via megalin-mediated reabsorption "
            "of filtered ASO metabolites. Causes tubular basophilic granulation "
            "and, at high doses, tubular degeneration. Glomerulonephritis "
            "(immune-complex mediated) reported rarely with inotersen."
        ),
        incidence_at_therapeutic_dose=(
            "Tubular basophilic granulation: common (histological finding, "
            "usually without clinical significance). Clinically significant "
            "nephrotoxicity: <5% at therapeutic doses. Glomerulonephritis: "
            "rare (<1% with inotersen at higher human doses)."
        ),
        severity="mild",
        dose_dependent=True,
        reversible=True,
        monitoring_strategy=(
            "Serum creatinine, BUN, urine protein:creatinine ratio at baseline "
            "and monthly. Urinalysis for sediment. If proteinuria develops, "
            "reduce dose and monitor weekly."
        ),
        clinical_significance=(
            "Low risk at proposed dose. Renal accumulation is less therapeutically "
            "relevant than hepatic/splenic for leishmaniasis. Monitor as "
            "standard precaution for PS-ASO class."
        ),
        risk_score=3.0,
    ))

    # 5. Efeitos pro-inflamatorios (CpG/TLR9)
    endpoints.append(ToxicityEndpoint(
        name="Pro-inflammatory effects (TLR9/CpG)",
        mechanism=(
            f"PS-ASO with {n_cpg} CpG motif(s) activates TLR9 in plasmacytoid "
            "DCs, B cells, and macrophages. Produces IFN-alpha, IL-6, TNF-alpha. "
            "In dogs, TLR9 is particularly responsive to unmethylated CpG DNA. "
            "This is BOTH a therapeutic mechanism (anti-Leishmania immunity) "
            "and a potential toxicity (systemic inflammation)."
        ),
        incidence_at_therapeutic_dose=(
            "Transient post-injection cytokine elevation: >50% expected. "
            "Clinically significant inflammation (fever, malaise): 10-20%. "
            "Self-limiting flu-like symptoms: 15-30% (extrapolated from "
            "CpG-ODN clinical trials in dogs)."
        ),
        severity="mild to moderate",
        dose_dependent=True,
        reversible=True,
        monitoring_strategy=(
            "Body temperature 6-24h post-injection. Serum CRP at each visit. "
            "Cytokine panel (IFN-alpha, IL-6) if clinical signs develop. "
            "Dose reduction if persistent fever >40C or CRP >5x baseline."
        ),
        clinical_significance=(
            "DUAL FUNCTION: the same TLR9 activation that causes transient "
            "inflammation is THERAPEUTICALLY BENEFICIAL for leishmaniasis. "
            "The clinical challenge is dose optimization to maximize anti-parasitic "
            "immune activation while minimizing systemic inflammatory symptoms. "
            "This is a FEATURE of MRL-ASO-001, not merely a side effect."
        ),
        risk_score=4.0,
    ))

    # 6. Coagulopatia
    endpoints.append(ToxicityEndpoint(
        name="Coagulopathy (complement activation)",
        mechanism=(
            "PS backbone activates the alternative complement pathway at high "
            "concentrations (>10 ug/mL plasma). C3a and C5a anaphylatoxins "
            "can cause transient hypotension, bronchospasm (rare), and "
            "complement-mediated platelet activation. Primarily a concern "
            "with IV bolus dosing; less relevant for SC administration."
        ),
        incidence_at_therapeutic_dose=(
            "Clinically significant complement activation: <2% with SC dosing "
            "(Cmax after SC is 10-100x lower than IV bolus). Subclinical "
            "complement activation markers may be elevated in ~10%."
        ),
        severity="mild (SC route)",
        dose_dependent=True,
        reversible=True,
        monitoring_strategy=(
            "aPTT, PT/INR at baseline and monthly. Complement C3, C4 if "
            "clinical signs of complement activation (urticaria, hypotension). "
            "SC route inherently limits peak plasma concentrations."
        ),
        clinical_significance=(
            "Low risk with SC administration. The slow absorption from SC depot "
            "(Tmax 2-4h) prevents the high peak plasma concentrations that "
            "trigger complement activation. This is a key advantage of SC "
            "over IV dosing for PS-ASOs."
        ),
        risk_score=2.0,
    ))

    return endpoints


def compute_safety_profile(aso_sequence: str) -> SafetyProfile:
    """Calcula perfil de seguranca completo do ASO.

    Integra todas as avaliacoes:
    - Endpoints toxicos de classe
    - Perfis por orgao
    - Avaliacao de TLR9
    - Indice terapeutico
    - Recomendacoes de monitoramento

    A classificacao de risco global considera:
    - Indice terapeutico (dose NOAEL / dose terapeutica)
    - Score medio dos endpoints
    - Presenca de efeitos severos
    - Contexto da doenca (leishmaniose visceral e FATAL sem tratamento)

    Args:
        aso_sequence: Sequencia do ASO.

    Returns:
        Perfil de seguranca integrado.
    """
    endpoints = assess_class_effects(aso_sequence)
    organ_profiles = assess_organ_toxicity()
    tlr9 = assess_tlr9(aso_sequence)

    # Indice terapeutico
    ti = NOAEL_MG_KG_WEEK / THERAPEUTIC_DOSE_MG_KG_WEEK

    # Margem de seguranca (MTD / dose)
    safety_margin = MTD_MG_KG_WEEK / THERAPEUTIC_DOSE_MG_KG_WEEK

    # Score medio de risco
    mean_risk = sum(e.risk_score for e in endpoints) / len(endpoints) if endpoints else 0

    # Classificacao global
    # Contexto: leishmaniose visceral canina tem mortalidade >90% sem tratamento
    # Um risco moderado e aceitavel para uma doenca letal
    if ti >= 8.0 and mean_risk < 4.0:
        risk_class = "ACCEPTABLE — favorable therapeutic index for a lethal disease"
    elif ti >= 4.0 and mean_risk < 6.0:
        risk_class = "MANAGEABLE — requires monitoring but justified for visceral leishmaniasis"
    else:
        risk_class = "CAUTION — close monitoring required, consider dose adjustment"

    # Cronograma de monitoramento
    monitoring = {
        "baseline": (
            "CBC with platelet count, serum chemistry (ALT, AST, ALP, GGT, "
            "creatinine, BUN), urine protein:creatinine ratio, coagulation "
            "panel (aPTT, PT/INR), abdominal ultrasound (spleen, liver), "
            "Leishmania parasite load (qPCR), cytokine panel (IFN-alpha, IL-6)"
        ),
        "loading_phase_weekly": (
            "CBC with platelet count, body temperature post-injection, "
            "injection site assessment, ALT/AST"
        ),
        "maintenance_biweekly": (
            "CBC with platelet count, serum chemistry panel, "
            "urine protein:creatinine ratio"
        ),
        "maintenance_monthly": (
            "Full serum chemistry, abdominal ultrasound, Leishmania qPCR "
            "(treatment response), cytokine panel"
        ),
        "end_of_treatment": (
            "Full baseline panel repeated. Leishmania qPCR for clearance. "
            "Follow-up at 1, 3, 6 months post-treatment for relapse."
        ),
    }

    # Contraindicacoes
    contraindications = [
        "Severe thrombocytopenia (<50k/uL) — pre-existing from leishmaniasis or other cause",
        "Severe hepatic failure (ALT >10x ULN) — cannot tolerate additional hepatic ASO load",
        "Severe renal failure (GFR <30 mL/min) — impaired excretion of ASO metabolites",
        "Active hemorrhage or coagulopathy — PS-ASOs may exacerbate bleeding risk",
        "Pregnancy — teratogenicity not evaluated; do not use in pregnant animals",
        "Known hypersensitivity to phosphorothioate oligonucleotides",
    ]

    # Interacoes medicamentosas
    drug_interactions = [
        "Anticoagulants (heparin, warfarin): theoretical additive bleeding risk — monitor closely",
        "NSAIDs: additive platelet dysfunction risk — avoid concurrent use if possible",
        "Allopurinol: no pharmacokinetic interaction expected (different metabolism pathways)",
        "Miltefosine: no pharmacokinetic interaction; potential additive nephrotoxicity — monitor renal function",
        "Antimonials (meglumine antimoniate): potential additive cardiotoxicity — ECG monitoring recommended",
        "Immunosuppressants: may blunt TLR9-mediated therapeutic immune activation — use with caution",
    ]

    return SafetyProfile(
        therapeutic_index=round(ti, 1),
        noael_mg_kg_week=NOAEL_MG_KG_WEEK,
        mtd_mg_kg_week=MTD_MG_KG_WEEK,
        dose_mg_kg_week=THERAPEUTIC_DOSE_MG_KG_WEEK,
        safety_margin=round(safety_margin, 1),
        endpoints=endpoints,
        organ_profiles=organ_profiles,
        tlr9_assessment=tlr9,
        overall_risk_classification=risk_class,
        monitoring_schedule=monitoring,
        contraindications=contraindications,
        drug_interactions=drug_interactions,
    )
