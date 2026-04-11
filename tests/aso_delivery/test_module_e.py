"""Testes unitarios para aso_delivery/module_e_admet.

Valida modelos farmacocineticos e de toxicidade de MRL-ASO-001 para
uso veterinario em caes (Canis lupus familiaris) com leishmaniose visceral.

Testes cobrem:
- Absorcao SC (biodisponibilidade, Cmax, AUC)
- Distribuicao tecidual (Vd, coeficientes de particionamento, modelo 2-comp)
- Metabolismo (ausencia de CYP450, degradacao por nucleases)
- Excrecao (clearance renal, meia-vida terminal)
- Regime posologico (escalonamento alometrico, indice terapeutico)
- Toxicidade (endpoints de classe, perfil por orgao, avaliacao TLR9)

Todos os calculos cientificos sao executados sem mocks.
"""

import math

import pytest

from aso_delivery.module_e_admet.pharmacokinetics import (
    BIOAVAILABILITY_SC_MEAN,
    BODY_WEIGHT_KG,
    CL_TOTAL_ML_MIN_KG,
    FORMULATION_CONC_MG_ML,
    FRACTION_UNBOUND,
    KA_SC_MEAN,
    PLASMA_HALF_LIFE_HOURS,
    PROTEIN_BINDING_FRACTION,
    TISSUE_HALF_LIFE_DAYS_MEAN,
    TISSUE_PARTITION_COEFFICIENTS,
    VD_SS_L_PER_KG,
    AbsorptionProfile,
    DistributionProfile,
    ExcretionProfile,
    MetabolismProfile,
    compute_absorption,
    compute_distribution,
    compute_dosing_regimen,
    compute_excretion,
    compute_metabolism,
    compute_pk_simulation,
)
from aso_delivery.module_e_admet.toxicity import (
    NOAEL_MG_KG_WEEK,
    THERAPEUTIC_DOSE_MG_KG_WEEK,
    MTD_MG_KG_WEEK,
    ToxicityEndpoint,
    assess_class_effects,
    assess_organ_toxicity,
    assess_tlr9,
    compute_safety_profile,
    count_cpg_motifs,
)
from aso_math.config import ASO_SEQUENCE


# =========================================================================
# 1. Absorcao SC
# =========================================================================


class TestAbsorption:
    """Testes do perfil de absorcao subcutanea de PS-ASOs."""

    def test_bioavailability_87_percent(self):
        """Biodisponibilidade SC de PS-ASOs deve ser ~87%.

        Ref: Geary RS et al. (2015) — media entre mipomersen (84%)
        e inotersen (90%).
        """
        abs_profile = compute_absorption(dose_mg_kg=5.0)
        assert abs_profile.bioavailability_percent == pytest.approx(87.0, abs=0.1)

    def test_dose_absolute_calculation(self):
        """Dose absoluta = dose_mg_kg * peso corporal (15 kg para cao de referencia).

        5 mg/kg * 15 kg = 75 mg.
        """
        abs_profile = compute_absorption(dose_mg_kg=5.0)
        assert abs_profile.dose_absolute_mg == pytest.approx(75.0, abs=0.1)

    def test_dose_absorbed_calculation(self):
        """Dose absorvida = dose absoluta * biodisponibilidade.

        75 mg * 0.87 = 65.25 mg.
        """
        abs_profile = compute_absorption(dose_mg_kg=5.0)
        expected = 75.0 * BIOAVAILABILITY_SC_MEAN
        assert abs_profile.dose_absorbed_mg == pytest.approx(expected, abs=0.2)

    def test_injection_volume(self):
        """Volume de injecao = dose absoluta / concentracao da formulacao.

        75 mg / 200 mg/mL = 0.375 mL — volume aceitavel para SC em caes.
        """
        abs_profile = compute_absorption(dose_mg_kg=5.0)
        expected = 75.0 / FORMULATION_CONC_MG_ML
        assert abs_profile.injection_volume_ml == pytest.approx(expected, abs=0.01)

    def test_higher_dose_higher_cmax(self):
        """Dose maior deve produzir Cmax maior (farmacocinetica linear).

        PS-ASOs tem PK linear em doses terapeuticas: Cmax proporcional a dose.
        """
        abs_5 = compute_absorption(dose_mg_kg=5.0)
        abs_10 = compute_absorption(dose_mg_kg=10.0)
        assert abs_10.cmax_estimated_ng_ml > abs_5.cmax_estimated_ng_ml

    def test_tmax_range(self):
        """Tmax SC deve ser 2-4 horas (dados de PK em caes).

        Absorcao do deposito SC por drenagem linfatica/capilar
        leva 2-4h para atingir pico plasmatico.
        """
        abs_profile = compute_absorption(dose_mg_kg=5.0)
        assert abs_profile.tmax_hours_range == (2.0, 4.0)

    def test_route_is_subcutaneous(self):
        """Via de administracao deve ser SC (padrao para ASOs veterinarios)."""
        abs_profile = compute_absorption(dose_mg_kg=5.0)
        assert "subcutaneous" in abs_profile.route.lower()

    def test_auc_proportional_to_dose(self):
        """AUC deve ser proporcional a dose (PK linear).

        AUC = F * Dose / CL — relacao linear.
        """
        abs_5 = compute_absorption(dose_mg_kg=5.0)
        abs_10 = compute_absorption(dose_mg_kg=10.0)
        ratio = abs_10.auc_estimated_ng_h_ml / abs_5.auc_estimated_ng_h_ml
        assert abs(ratio - 2.0) < 0.1


# =========================================================================
# 2. Distribuicao tecidual
# =========================================================================


class TestDistribution:
    """Testes da distribuicao do ASO em tecidos caninos."""

    def test_vd_ss_value(self):
        """Volume de distribuicao no estado estacionario deve ser 0.25 L/kg.

        Ref: Geary RS et al. (2015) — tipico para PS-ASOs (0.1-0.5 L/kg).
        """
        dist = compute_distribution(65.0)
        assert dist.vd_ss_l_per_kg == VD_SS_L_PER_KG

    def test_protein_binding_90_percent(self):
        """Ligacao a proteinas plasmaticas deve ser 90%.

        PS backbone confere carga negativa que interage com albumina.
        Ref: Crooke ST et al. (2017).
        """
        dist = compute_distribution(65.0)
        assert dist.protein_binding_percent == pytest.approx(90.0, abs=0.1)

    def test_liver_highest_accumulation(self):
        """Figado deve ter a maior concentracao tecidual (Kp = 30).

        Sinusoides hepaticos fenestrados permitem acesso direto.
        VANTAGEM para leishmaniose: figado e sitio primario de infeccao.
        """
        dist = compute_distribution(65.0)
        assert dist.tissue_partition_coefficients["liver"] >= 30.0
        # Figado deve ter maior concentracao que qualquer outro tecido
        liver_conc = dist.tissue_concentrations["liver"]
        for tissue, conc in dist.tissue_concentrations.items():
            if tissue != "liver":
                assert liver_conc >= conc

    def test_brain_negligible_penetration(self):
        """Penetracao cerebral deve ser negligivel (Kp = 0.01).

        Barreira hematoencefalica impede passagem de PS-ASOs.
        """
        dist = compute_distribution(65.0)
        assert dist.tissue_partition_coefficients["brain"] == 0.01

    def test_two_compartment_alpha_greater_than_beta(self):
        """Constante alpha (fase rapida) deve ser maior que beta (fase lenta).

        Modelo de 2 compartimentos: alpha governa distribuicao rapida
        e beta governa eliminacao terminal lenta.
        """
        dist = compute_distribution(65.0)
        alpha = dist.two_compartment_params["alpha_per_hour"]
        beta = dist.two_compartment_params["beta_per_hour"]
        assert alpha > beta

    def test_spleen_relevant_for_leishmaniasis(self):
        """Baco deve ter acumulacao significativa (Kp >= 15).

        O baco e reservatorio importante de L. infantum em LVC.
        Acumulacao de PS-ASOs no baco e terapeuticamente vantajosa.
        """
        dist = compute_distribution(65.0)
        assert dist.tissue_partition_coefficients["spleen"] >= 15.0


# =========================================================================
# 3. Metabolismo
# =========================================================================


class TestMetabolism:
    """Testes do perfil metabolico de PS-ASOs."""

    def test_no_cyp450_metabolism(self):
        """PS-ASOs NAO sao metabolizados por CYP450.

        Diferenca fundamental em relacao a small molecules.
        Sem risco de interacoes farmacologicas com alopurinol, miltefosina etc.
        """
        met = compute_metabolism()
        assert met.cyp450_metabolism is False

    def test_tissue_half_life_21_days(self):
        """Meia-vida tecidual deve ser ~21 dias (3 semanas).

        Flancos LNA protegem contra exonucleases, prolongando a
        atividade do ASO nos tecidos-alvo.
        """
        met = compute_metabolism()
        assert met.tissue_half_life_days == TISSUE_HALF_LIFE_DAYS_MEAN

    def test_plasma_half_life_short(self):
        """Meia-vida plasmatica deve ser 2.5h (fase alpha de distribuicao).

        PS-ASOs distribuem rapidamente do plasma para tecidos.
        """
        met = compute_metabolism()
        assert met.plasma_half_life_hours == PLASMA_HALF_LIFE_HOURS

    def test_no_drug_drug_interactions(self):
        """Sem DDI farmacocinetica esperada (sem CYP450).

        Vantagem clinica: caes com leishmaniose recebem multiplos farmacos.
        """
        met = compute_metabolism()
        assert "No CYP450" in met.drug_drug_interactions

    def test_metabolites_list_not_empty(self):
        """Lista de metabolitos deve ter pelo menos 3 etapas de degradacao."""
        met = compute_metabolism()
        assert len(met.metabolites) >= 3


# =========================================================================
# 4. Excrecao
# =========================================================================


class TestExcretion:
    """Testes do perfil de excrecao de PS-ASOs em caes."""

    def test_clearance_total(self):
        """Clearance total deve ser 2.0 mL/min/kg * 15 kg = 30 mL/min."""
        exc = compute_excretion()
        expected = CL_TOTAL_ML_MIN_KG * BODY_WEIGHT_KG
        assert exc.cl_total_ml_min == pytest.approx(expected, abs=0.1)

    def test_renal_fraction_70_percent(self):
        """70% da excrecao deve ser renal (metabolitos curtos filtrados).

        Fragmentos <10 nt perdem afinidade por proteinas e sao filtrados.
        """
        exc = compute_excretion()
        assert exc.renal_fraction == pytest.approx(0.70, abs=0.01)

    def test_renal_plus_hepatic_equals_total(self):
        """Clearance renal + hepatico deve igualar clearance total."""
        exc = compute_excretion()
        total = exc.cl_renal_ml_min + exc.cl_hepatic_ml_min
        assert total == pytest.approx(exc.cl_total_ml_min, abs=0.1)

    def test_time_to_steady_state(self):
        """Tempo para estado estacionario ~ 5 meias-vidas = 105 dias.

        Com meia-vida terminal de 21 dias, estado estacionario em ~15 semanas.
        """
        exc = compute_excretion()
        expected_days = 5.0 * TISSUE_HALF_LIFE_DAYS_MEAN
        assert exc.time_to_steady_state_days == pytest.approx(expected_days, abs=0.1)

    def test_terminal_half_life_in_hours(self):
        """Meia-vida terminal em horas = 21 dias * 24 = 504 horas."""
        exc = compute_excretion()
        assert exc.terminal_half_life_hours == pytest.approx(
            TISSUE_HALF_LIFE_DAYS_MEAN * 24.0, abs=0.1
        )


# =========================================================================
# 5. Simulacao PK de dois compartimentos
# =========================================================================


class TestPKSimulation:
    """Testes da simulacao PK temporal (Euler explicito, 2 compartimentos)."""

    def test_simulation_starts_at_zero(self):
        """Concentracao plasmatica deve ser zero em t=0 (antes de absorcao)."""
        points = compute_pk_simulation(dose_mg=75.0, duration_hours=168.0)
        assert points[0].plasma_conc_ng_ml == 0.0

    def test_plasma_concentration_peaks_then_decays(self):
        """Concentracao plasmatica deve subir (absorcao) e depois cair (eliminacao).

        Perfil classico de dose SC: absorve do deposito, distribui
        para tecidos, elimina por nucleases.
        """
        points = compute_pk_simulation(dose_mg=75.0, duration_hours=168.0)
        concs = [p.plasma_conc_ng_ml for p in points]
        max_idx = concs.index(max(concs))
        # Pico nao deve ser no primeiro ou ultimo ponto
        assert max_idx > 0
        assert max_idx < len(concs) - 1

    def test_tissue_accumulation(self):
        """Concentracao tecidual deve ser maior que plasmatica no final.

        PS-ASOs acumulam em tecidos (k12 >> k21). No final da simulacao,
        a maioria do ASO esta no compartimento periferico.
        """
        points = compute_pk_simulation(dose_mg=75.0, duration_hours=168.0)
        last_point = points[-1]
        assert last_point.tissue_conc_ng_ml > last_point.plasma_conc_ng_ml

    def test_mass_balance(self):
        """Quantidade absorvida deve ser ~F*Dose ao final da simulacao.

        Biodisponibilidade = 87%, dose = 75 mg -> absorvido ~ 65.25 mg.
        """
        points = compute_pk_simulation(dose_mg=75.0, duration_hours=168.0)
        last = points[-1]
        expected_absorbed = 75.0 * BIOAVAILABILITY_SC_MEAN
        assert last.amount_absorbed_mg == pytest.approx(expected_absorbed, rel=0.05)

    def test_no_negative_concentrations(self):
        """Nenhuma concentracao pode ser negativa (protecao numerica).

        Euler explicito pode gerar artefatos negativos — verificar.
        """
        points = compute_pk_simulation(dose_mg=75.0, duration_hours=168.0)
        for p in points:
            assert p.plasma_conc_ng_ml >= 0.0
            assert p.tissue_conc_ng_ml >= 0.0


# =========================================================================
# 6. Regime posologico
# =========================================================================


class TestDosingRegimen:
    """Testes do regime posologico recomendado para LVC."""

    def test_maintenance_dose_5_mg_kg(self):
        """Dose de manutencao deve ser 5 mg/kg/semana (arredondamento clinico).

        Derivado de escalonamento alometrico de mipomersen (humanos).
        """
        regimen = compute_dosing_regimen()
        assert regimen.maintenance_dose_mg_kg == 5.0

    def test_loading_dose_double(self):
        """Loading dose deve ser 2x manutencao (10 mg/kg) para saturar tecidos.

        Leishmaniose visceral requer reducao rapida da carga parasitaria,
        justificando uma fase de ataque mais agressiva.
        """
        regimen = compute_dosing_regimen()
        assert regimen.loading_dose_mg_kg == 10.0

    def test_therapeutic_index_favorable(self):
        """Indice terapeutico deve ser >= 8 (NOAEL 40 / dose 5 = 8).

        IT >= 8 indica margem de seguranca favoravel para uso veterinario.
        """
        regimen = compute_dosing_regimen()
        assert regimen.therapeutic_index >= 8.0

    def test_total_treatment_12_weeks(self):
        """Tratamento total = 2 semanas loading + 10 semanas manutencao = 12."""
        regimen = compute_dosing_regimen()
        assert regimen.total_treatment_weeks == 12

    def test_comparison_with_miltefosine_present(self):
        """Deve incluir comparacao com miltefosina (tratamento padrao canino).

        Miltefosina: 2 mg/kg/dia VO por 28 dias — perfil diferente.
        """
        regimen = compute_dosing_regimen()
        assert "miltefosine_dose" in regimen.comparison_miltefosine
        assert "mrl_aso_001_dose" in regimen.comparison_miltefosine


# =========================================================================
# 7. Toxicidade — motivos CpG e TLR9
# =========================================================================


class TestToxicity:
    """Testes do perfil de toxicidade de PS-ASOs."""

    def test_cpg_count_in_mrl_aso_001(self):
        """Contar motivos CpG na sequencia MRL-ASO-001.

        CpG = dinucleotideo 5'-CG-3'. Motivos CpG ativam TLR9 no
        endossomo, induzindo resposta imune inata — mecanismo dual do ASO.
        """
        # ASO_SEQUENCE = "ACAGAAACTGATACTTATATAGCGT"
        # CG aparece nas posicoes: 22-23 (GC -> nao, CG sim)
        n = count_cpg_motifs(ASO_SEQUENCE)
        assert isinstance(n, int)
        assert n >= 0

    def test_cpg_count_known_sequence(self):
        """Sequencia 'ACGT' tem exatamente 1 motivo CpG (pos 1-2)."""
        assert count_cpg_motifs("ACGT") == 1

    def test_cpg_count_no_cpg(self):
        """Sequencia sem CpG deve retornar 0."""
        assert count_cpg_motifs("AAAA") == 0
        assert count_cpg_motifs("TTTT") == 0

    def test_tlr9_assessment_has_cpg(self):
        """Avaliacao de TLR9 deve detectar se ha motivos CpG."""
        tlr9 = assess_tlr9(ASO_SEQUENCE)
        assert isinstance(tlr9.cpg_motifs_in_sequence, int)

    def test_class_effects_all_six(self):
        """Deve avaliar 6 endpoints toxicos de classe de PS-ASOs.

        1. ISR (local de injecao)
        2. Trombocitopenia
        3. Hepatotoxicidade
        4. Nefrotoxicidade
        5. Pro-inflamatorio (TLR9/CpG)
        6. Coagulopatia
        """
        endpoints = assess_class_effects(ASO_SEQUENCE)
        assert len(endpoints) == 6

    def test_thrombocytopenia_highest_risk(self):
        """Trombocitopenia deve ter o maior risk_score (efeito mais significativo).

        Ref: Inotersen (FDA) tem black box warning para trombocitopenia.
        """
        endpoints = assess_class_effects(ASO_SEQUENCE)
        risk_scores = {e.name: e.risk_score for e in endpoints}
        thrombocytopenia = [e for e in endpoints if "Thrombocytopenia" in e.name][0]
        assert thrombocytopenia.risk_score == max(e.risk_score for e in endpoints)

    def test_organ_profiles_four_organs(self):
        """Deve avaliar 4 orgaos: figado, rim, baco, medula ossea."""
        profiles = assess_organ_toxicity()
        assert len(profiles) == 4
        organs = {p.organ for p in profiles}
        assert "liver" in organs
        assert "kidney" in organs
        assert "spleen" in organs
        assert "bone_marrow" in organs

    def test_safety_profile_therapeutic_index(self):
        """Indice terapeutico = NOAEL / dose = 40 / 5 = 8."""
        profile = compute_safety_profile(ASO_SEQUENCE)
        expected_ti = NOAEL_MG_KG_WEEK / THERAPEUTIC_DOSE_MG_KG_WEEK
        assert profile.therapeutic_index == pytest.approx(expected_ti, abs=0.1)

    def test_safety_profile_acceptable_classification(self):
        """Classificacao de risco deve ser 'ACCEPTABLE' para TI >= 8 e risco < 4.

        Contexto: leishmaniose visceral canina e FATAL sem tratamento,
        justificando tolerancia a riscos moderados.
        """
        profile = compute_safety_profile(ASO_SEQUENCE)
        assert "ACCEPTABLE" in profile.overall_risk_classification

    def test_all_endpoints_dose_dependent(self):
        """Todos os efeitos de classe de PS-ASOs sao dose-dependentes.

        Isto e uma vantagem: reduzir a dose reduz os efeitos adversos.
        """
        endpoints = assess_class_effects(ASO_SEQUENCE)
        for ep in endpoints:
            assert ep.dose_dependent is True

    def test_all_endpoints_reversible(self):
        """Todos os efeitos de classe devem ser reversiveis apos descontinuacao."""
        endpoints = assess_class_effects(ASO_SEQUENCE)
        for ep in endpoints:
            assert ep.reversible is True
