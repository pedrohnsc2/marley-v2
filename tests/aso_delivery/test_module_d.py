"""Testes unitarios para aso_delivery/module_d_lnp.

Valida os modelos de formulacao LNP e liberacao pH-responsiva para
encapsulacao de MRL-ASO-001 visando macrofagos infectados por L. infantum.

Testes cobrem:
- Composicao lipidica (razoes molares, pesos moleculares)
- Razao N/P e eficiencia de encapsulacao
- Cinetica de liberacao pH-responsiva (Henderson-Hasselbalch + sigmoidal)
- Targeting de macrofagos via manose-PEG
- Estabilidade em armazenamento (cinetica de Arrhenius)

Todos os calculos cientificos sao executados sem mocks — os valores
sao verificados contra propriedades fisico-quimicas conhecidas.
"""

import math

import pytest

from aso_delivery.module_d_lnp.formulation import (
    BASE_DIAMETER_NM,
    DIAMETER_GROWTH_PER_NP,
    EA_AGGREGATION,
    K_AGG_25C,
    K_ENCAPSULATION,
    MANNOSE_UPTAKE_FACTOR_MAX,
    MANNOSE_UPTAKE_FACTOR_MEAN,
    MANNOSE_UPTAKE_FACTOR_MIN,
    MAX_PHAGOCYTIC_NM,
    MIN_PHAGOCYTIC_NM,
    MOLAR_RATIO_CHOLESTEROL,
    MOLAR_RATIO_DSPC,
    MOLAR_RATIO_IONIZABLE,
    MOLAR_RATIO_PEG,
    MW_CHOLESTEROL,
    MW_DLIN_MC3_DMA,
    MW_DSPC,
    MW_PEG_LIPID,
    PKA_DLIN_MC3,
    R_GAS_KJ,
    T_REF_K,
    LNPComposition,
    NPRatioResult,
    compute_lnp_composition,
    compute_macrophage_targeting,
    compute_np_ratio,
    compute_storage_profiles,
    compute_storage_stability,
    find_optimal_np,
    scan_np_ratios,
)
from aso_delivery.module_d_lnp.release import (
    COMPARTMENT_PH,
    HILL_SLOPE,
    PKA_MC3,
    TRANSIT_TIMES_MIN,
    classify_release,
    compute_fraction_protonated,
    compute_fraction_released,
    compute_release_profile,
)


# =========================================================================
# 1. Composicao lipidica
# =========================================================================


class TestLNPComposition:
    """Testes da composicao lipidica da formulacao LNP (Onpattro-like)."""

    def test_molar_fractions_sum_to_one(self):
        """Razoes molares dos 4 componentes devem somar 1.0.

        Principio basico de formulacao: as fracoes molares de todos
        os componentes em uma mistura devem somar exatamente 1.
        """
        total = (
            MOLAR_RATIO_IONIZABLE
            + MOLAR_RATIO_DSPC
            + MOLAR_RATIO_CHOLESTEROL
            + MOLAR_RATIO_PEG
        )
        assert abs(total - 1.0) < 1e-9

    def test_composition_components(self):
        """Formulacao deve ter os 4 componentes corretos (DLin-MC3-DMA, DSPC,
        colesterol, DMG-PEG2000)."""
        comp = compute_lnp_composition()
        assert comp.ionizable_lipid == "DLin-MC3-DMA"
        assert comp.helper_lipid == "DSPC"
        assert comp.sterol == "Cholesterol"
        assert comp.peg_lipid == "DMG-PEG2000"

    def test_weighted_average_mw(self):
        """Peso molecular medio ponderado deve corresponder ao calculo manual.

        MW_avg = sum(x_i * MW_i) para os 4 componentes.
        """
        comp = compute_lnp_composition()
        expected_mw = (
            MOLAR_RATIO_IONIZABLE * MW_DLIN_MC3_DMA
            + MOLAR_RATIO_DSPC * MW_DSPC
            + MOLAR_RATIO_CHOLESTEROL * MW_CHOLESTEROL
            + MOLAR_RATIO_PEG * MW_PEG_LIPID
        )
        assert abs(comp.weighted_average_mw - round(expected_mw, 2)) < 0.01

    def test_charge_neutral_at_physiological_ph(self):
        """A pH 7.4 (sangue), a LNP deve ter carga baixa (stealth).

        DLin-MC3-DMA tem pKa 6.44. A pH 7.4, a maior parte do lipidio
        ionizavel esta desprotonado, conferindo carga liquida baixa.
        Henderson-Hasselbalch: f = 1/(1+10^(7.4-6.44)) ~ 0.099
        Carga liquida = 0.50 * 0.099 ~ 0.05 (stealth suficiente).
        """
        comp = compute_lnp_composition()
        # Carga baixa: <0.10 (vs ~0.49 a pH 4.5)
        assert comp.charge_ph_7_4 < 0.10, (
            "Carga a pH 7.4 deve ser baixa para propriedade stealth"
        )

    def test_charge_positive_at_phagolysosome_ph(self):
        """A pH 4.5 (fagolisossomo), o lipidio ionizavel deve estar protonado.

        No fagolisossomo onde L. infantum reside, o pH acido protona
        MC3 quase completamente, gerando carga positiva que desestabiliza
        a LNP e libera o ASO.
        """
        comp = compute_lnp_composition()
        # Fracao protonada a pH 4.5: 1/(1+10^(4.5-6.44)) = ~0.989
        # Carga liquida = 0.50 * 1 * 0.989 = ~0.495
        assert comp.charge_ph_4_5 > 0.45

    def test_pka_of_ionizable_lipid(self):
        """pKa do DLin-MC3-DMA deve ser 6.44 (Jayaraman 2012)."""
        comp = compute_lnp_composition()
        assert comp.pka_ionizable == 6.44


# =========================================================================
# 2. Razao N/P e eficiencia de encapsulacao
# =========================================================================


class TestNPRatio:
    """Testes de otimizacao da razao N/P (amina:fosfato)."""

    def test_encapsulation_increases_with_np_ratio(self):
        """Eficiencia de encapsulacao deve aumentar com a razao N/P.

        Modelo exponencial saturante: EE = 1 - exp(-k * N/P).
        Mais lipidio cationico -> maior complexacao com o ASO aniônico.
        """
        low = compute_np_ratio(2.0, 25)
        high = compute_np_ratio(8.0, 25)
        assert high.encapsulation_efficiency > low.encapsulation_efficiency

    def test_encapsulation_never_exceeds_99_percent(self):
        """EE tem teto de 99% — encapsulacao perfeita e fisicamente impossivel."""
        extreme = compute_np_ratio(50.0, 25)
        assert extreme.encapsulation_efficiency <= 0.99

    def test_phosphate_count_for_25mer(self):
        """ASO de 25 nt com backbone PS tem 24 ligacoes fosforotioato.

        Cada ligacao internucleotidica carrega uma carga negativa (-1),
        total de (comprimento - 1) cargas para o ASO.
        """
        result = compute_np_ratio(4.0, 25)
        assert result.n_phosphates == 24

    def test_particle_diameter_increases_with_np(self):
        """Diametro da particula cresce linearmente com N/P.

        Mais lipidio = particulas maiores. Modelo: d = d_base + growth * N/P.
        """
        r1 = compute_np_ratio(2.0, 25)
        r2 = compute_np_ratio(10.0, 25)
        assert r2.particle_diameter_nm > r1.particle_diameter_nm
        # Verificar modelo linear
        expected_diff = DIAMETER_GROWTH_PER_NP * (10.0 - 2.0)
        actual_diff = r2.particle_diameter_nm - r1.particle_diameter_nm
        assert abs(actual_diff - expected_diff) < 0.2

    def test_optimal_range_is_4_to_8(self):
        """Faixa otima de N/P para PS-ASOs e 4-8 (diferente de mRNA: 6-12).

        PS-ASOs tem carga negativa permanente (fosforotioato), necessitando
        de menos lipidio ionizavel por carga do que mRNA.
        """
        in_range = compute_np_ratio(6.0, 25)
        below_range = compute_np_ratio(2.0, 25)
        above_range = compute_np_ratio(12.0, 25)
        assert in_range.optimal is True
        assert below_range.optimal is False
        assert above_range.optimal is False

    def test_scan_covers_full_range(self):
        """Varredura N/P de 1 a 20 (step 1) deve gerar 20 resultados."""
        results = scan_np_ratios(25, np_min=1.0, np_max=20.0, step=1.0)
        assert len(results) == 20

    def test_find_optimal_selects_from_optimal_range(self):
        """find_optimal_np deve retornar resultado com N/P entre 4-8.

        O criterio e: maior encapsulacao dentro da faixa otima E com
        tamanho compativel para fagocitose (50-200 nm).
        """
        results = scan_np_ratios(25)
        best = find_optimal_np(results)
        assert best.optimal is True
        assert best.suitable_for_phagocytosis is True

    def test_zeta_potential_trend(self):
        """Potencial zeta deve ser mais positivo com maior N/P.

        Modelo sigmoidal: zeta = -30 + 40/(1+exp(-0.8*(N/P - 3)))
        Excesso de lipidio cationico torna a superficie mais positiva.
        """
        low = compute_np_ratio(1.0, 25)
        high = compute_np_ratio(10.0, 25)
        assert high.zeta_potential_mv > low.zeta_potential_mv


# =========================================================================
# 3. Liberacao pH-responsiva
# =========================================================================


class TestPHRelease:
    """Testes do modelo de liberacao pH-responsiva da LNP."""

    def test_fraction_protonated_at_blood_ph(self):
        """A pH 7.4, MC3 esta ~10% protonado (Henderson-Hasselbalch).

        f = 1/(1+10^(pH-pKa)) = 1/(1+10^(7.4-6.44)) = 1/(1+10^0.96) ~ 0.099
        Maioria desprotonada confere propriedade stealth a LNP.
        """
        f = compute_fraction_protonated(7.4)
        assert f < 0.15

    def test_fraction_protonated_at_phagolysosome_ph(self):
        """A pH 4.5, MC3 esta >98% protonado.

        f = 1/(1+10^(4.5-6.44)) = ~0.989
        Protonacao massiva desestabiliza a LNP e libera a carga.
        """
        f = compute_fraction_protonated(4.5)
        assert f > 0.98

    def test_fraction_protonated_at_pka(self):
        """No pKa (6.44), exatamente 50% esta protonado.

        Henderson-Hasselbalch: f = 1/(1+10^0) = 0.5
        """
        f = compute_fraction_protonated(6.44)
        assert abs(f - 0.5) < 0.001

    def test_release_minimal_at_blood_ph(self):
        """Liberacao deve ser <10% a pH 7.4 (classificacao 'minimal').

        LNP deve permanecer intacta na circulacao para proteger o ASO.
        """
        f = compute_fraction_released(7.4)
        assert f < 0.10
        assert classify_release(f) == "minimal"

    def test_release_complete_at_phagolysosome(self):
        """Liberacao deve ser >80% a pH 4.5 (classificacao 'complete').

        No fagolisossomo, onde L. infantum reside, a LNP deve liberar
        todo o ASO para acao no local do parasita.
        """
        f = compute_fraction_released(4.5)
        assert f > 0.80
        assert classify_release(f) == "complete"

    def test_release_monotonically_increases_with_decreasing_ph(self):
        """Liberacao deve aumentar conforme pH diminui (mais acido = mais liberacao).

        Isto garante que a LNP funciona como 'cavalo de Troia': estavel
        fora da celula, libera dentro do compartimento acido.
        """
        phs = [7.4, 6.5, 5.0, 4.5]
        releases = [compute_fraction_released(ph) for ph in phs]
        for i in range(len(releases) - 1):
            assert releases[i + 1] > releases[i]

    def test_release_profile_all_compartments(self):
        """Perfil completo deve ter dados para 4 compartimentos.

        Jornada intracelular: sangue -> endossomo precoce ->
        endossomo tardio -> fagolisossomo.
        """
        profile = compute_release_profile()
        assert len(profile.compartment_releases) == 4
        assert "blood_plasma" in profile.compartment_releases
        assert "phagolysosome" in profile.compartment_releases

    def test_release_profile_target_is_phagolysosome(self):
        """Compartimento alvo deve ser o fagolisossomo.

        O fagolisossomo e onde amastigotas de L. infantum residem
        e onde o SL RNA alvo esta localizado.
        """
        profile = compute_release_profile()
        assert profile.target_compartment == "phagolysosome"
        assert profile.release_at_target > 0.80

    def test_classify_release_boundaries(self):
        """Classificacao deve respeitar limiares: <10% minimal, <80% partial, >=80% complete."""
        assert classify_release(0.0) == "minimal"
        assert classify_release(0.09) == "minimal"
        assert classify_release(0.10) == "partial"
        assert classify_release(0.79) == "partial"
        assert classify_release(0.80) == "complete"
        assert classify_release(1.0) == "complete"


# =========================================================================
# 4. Targeting de macrofagos
# =========================================================================


class TestMacrophageTargeting:
    """Testes de targeting via receptor de manose (CD206)."""

    def test_uptake_factor_within_literature_range(self):
        """Fator de uptake com manose-PEG deve estar entre 2-5x (Patel 2020).

        Macrofagos infectados por L. infantum expressam CD206 elevado,
        permitindo entrega seletiva da LNP.
        """
        targeting = compute_macrophage_targeting(80.0, mannose_peg_fraction=0.5)
        assert MANNOSE_UPTAKE_FACTOR_MIN <= targeting.uptake_fold_increase <= MANNOSE_UPTAKE_FACTOR_MAX

    def test_uptake_increases_with_mannose_fraction(self):
        """Mais manose-PEG na superficie -> maior uptake por macrofagos.

        Interpolacao linear: uptake = min + fraction * (max - min).
        """
        low = compute_macrophage_targeting(80.0, mannose_peg_fraction=0.2)
        high = compute_macrophage_targeting(80.0, mannose_peg_fraction=0.8)
        assert high.uptake_fold_increase > low.uptake_fold_increase

    def test_diameter_increase_is_modest(self):
        """Manose-PEG aumenta o diametro em ate ~4 nm (100% substituicao).

        A modificacao de superficie nao deve alterar dramaticamente
        o tamanho da particula, mantendo-a na faixa fagocitica.
        """
        targeting = compute_macrophage_targeting(80.0, mannose_peg_fraction=1.0)
        diameter_increase = targeting.particle_diameter_nm - 80.0
        assert 0 < diameter_increase <= 5.0

    def test_phagocytosis_suitability(self):
        """Particula com manose-PEG deve permanecer na faixa 50-200 nm.

        Macrofagos fagocitam particulas otimamente entre 50-200 nm.
        """
        targeting = compute_macrophage_targeting(80.0, mannose_peg_fraction=0.5)
        assert targeting.suitable_for_phagocytosis is True

    def test_receptor_target_is_cd206(self):
        """Receptor alvo deve ser CD206 (receptor de manose / MRC1)."""
        targeting = compute_macrophage_targeting(80.0)
        assert "CD206" in targeting.receptor_target


# =========================================================================
# 5. Estabilidade em armazenamento
# =========================================================================


class TestStorageStability:
    """Testes de estabilidade da LNP em diferentes temperaturas."""

    def test_arrhenius_lower_temp_slower_aggregation(self):
        """Temperatura mais baixa -> menor taxa de agregacao (Arrhenius).

        k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T))
        Armazenamento refrigerado preserva melhor a formulacao.
        """
        cold = compute_storage_stability(-20.0)
        warm = compute_storage_stability(25.0)
        assert cold.k_aggregation_per_day < warm.k_aggregation_per_day

    def test_half_life_longer_at_cold(self):
        """Meia-vida de armazenamento deve ser muito maior a -20C que a 25C.

        Essencial para distribuicao em areas endemicas de leishmaniose
        no Brasil, onde cadeia fria pode ser fragil.
        """
        cold = compute_storage_stability(-20.0)
        warm = compute_storage_stability(25.0)
        assert cold.half_life_days > warm.half_life_days * 5

    def test_fraction_intact_decreases_over_time(self):
        """Fracao intacta deve diminuir: 30 dias > 90 dias > 365 dias.

        Cinetica de primeira ordem: f(t) = exp(-k*t).
        """
        s = compute_storage_stability(4.0)
        assert s.fraction_intact_30_days > s.fraction_intact_90_days
        assert s.fraction_intact_90_days > s.fraction_intact_365_days

    def test_room_temperature_reference_rate(self):
        """A 25C, a taxa de agregacao deve ser k_ref = 0.015/dia.

        Esta e a temperatura de referencia no modelo de Arrhenius.
        """
        s = compute_storage_stability(25.0)
        assert abs(s.k_aggregation_per_day - K_AGG_25C) < 1e-6

    def test_storage_profiles_three_temperatures(self):
        """Deve calcular perfis para -20C, 4C, e 25C."""
        profiles = compute_storage_profiles()
        assert len(profiles) == 3
        temps = [p.temperature_celsius for p in profiles]
        assert -20.0 in temps
        assert 4.0 in temps
        assert 25.0 in temps

    def test_lyophilization_always_feasible(self):
        """Liofilizacao deve ser viavel para qualquer temperatura.

        LNPs com crioprotetores (sacarose/trealose) podem ser liofilizadas,
        eliminando dependencia de cadeia fria — essencial para areas rurais
        endemicas do Nordeste brasileiro.
        """
        for temp in [-20.0, 4.0, 25.0]:
            s = compute_storage_stability(temp)
            assert s.lyophilization_feasible is True

    def test_kelvin_conversion(self):
        """Conversao Celsius -> Kelvin deve ser correta (T_K = T_C + 273.15)."""
        s = compute_storage_stability(25.0)
        assert abs(s.temperature_kelvin - 298.15) < 0.01
