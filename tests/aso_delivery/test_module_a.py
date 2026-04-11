"""Testes unitarios para aso_delivery/module_a_stability.

Valida os modelos de estabilidade do ASO MRL-ASO-001 em pH acido:
1. Correcao termodinamica por protonacao (Henderson-Hasselbalch)
2. Cinetica de degradacao por nucleases (primeira ordem)
3. Estabilidade conformacional de nucleotideos LNA

Os testes usam valores reais — nenhum calculo cientifico e mockado.
Cada teste inclui comentario explicando a intuicao biologica.
"""

import math

import pytest

from aso_delivery.module_a_stability.models import (
    PKA_CYTOSINE_N3,
    PKA_ADENINE_N1,
    DDG_PROTONATION_C,
    DDG_PROTONATION_A,
    R_KCAL,
    T_PHYSIOLOGICAL_K,
    K_PO_PH45,
    K_PS_PH45,
    K_LNA_FACTOR,
    DNA_C3_ENDO_NEUTRAL,
    DNA_C3_ENDO_PH_SHIFT_PER_UNIT,
    fraction_protonated,
    compute_protonation_correction,
    compute_dg_at_ph,
    compute_tm_at_ph,
    compute_fraction_bound,
    compute_half_life,
    compute_fraction_remaining,
    compute_nuclease_resistance,
    compute_gapmer_degradation,
    compute_lna_conformation,
)


# ---------------------------------------------------------------------------
# 1. Henderson-Hasselbalch — fracao protonada
# ---------------------------------------------------------------------------


class TestFractionProtonated:
    """Testes para a funcao fraction_protonated (Henderson-Hasselbalch)."""

    def test_ph_equal_pka_gives_half(self):
        """Quando pH = pKa, exatamente metade das moleculas esta protonada.
        Principio fundamental de Henderson-Hasselbalch: ponto de inflexao.
        """
        assert fraction_protonated(4.2, 4.2) == pytest.approx(0.5, abs=1e-6)

    def test_very_low_ph_fully_protonated(self):
        """Em pH muito abaixo do pKa, quase todas as bases estao protonadas.
        No estomago (pH ~2), citosinas estariam 99%+ protonadas —
        desestabilizacao maxima do duplex.
        """
        f = fraction_protonated(1.0, PKA_CYTOSINE_N3)
        assert f > 0.99

    def test_physiological_ph_minimal_protonation(self):
        """Em pH fisiologico (7.4), muito acima do pKa da citosina (4.2),
        a protonacao e negligivel. O duplex ASO:alvo e estavel no sangue.
        """
        f = fraction_protonated(7.4, PKA_CYTOSINE_N3)
        assert f < 0.001

    def test_phagolysosome_ph_cytosine(self):
        """No fagolisossomo (pH 4.5), a fracao de citosinas protonadas
        e significativa mas nao dominante (~33%). Isso e o cenario real
        para MRL-ASO-001 dentro do macrofago infectado.
        """
        f = fraction_protonated(4.5, PKA_CYTOSINE_N3)
        # pH 4.5, pKa 4.2: f = 1/(1 + 10^0.3) = 1/(1+1.995) = 0.334
        assert f == pytest.approx(1.0 / (1.0 + 10.0 ** 0.3), abs=1e-4)

    def test_adenine_less_protonated_than_cytosine_at_same_ph(self):
        """Adenina tem pKa (3.5) menor que citosina (4.2), portanto
        e MENOS protonada no mesmo pH. A desestabilizacao por protonacao
        de adenina e menor — citosina e o elo mais fraco.
        """
        f_c = fraction_protonated(4.5, PKA_CYTOSINE_N3)
        f_a = fraction_protonated(4.5, PKA_ADENINE_N1)
        assert f_a < f_c

    def test_monotonically_decreasing_with_ph(self):
        """A fracao protonada deve diminuir monotonicamente com o aumento
        do pH. Quanto mais basico, menos protons disponiveis.
        """
        phs = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        fractions = [fraction_protonated(ph, PKA_CYTOSINE_N3) for ph in phs]
        for i in range(len(fractions) - 1):
            assert fractions[i] > fractions[i + 1]


# ---------------------------------------------------------------------------
# 2. Correcao de protonacao — dG, Tm e fracao ligada
# ---------------------------------------------------------------------------


class TestProtonationCorrection:
    """Testes para compute_protonation_correction e funcoes derivadas."""

    def test_neutral_ph_no_correction(self):
        """Em pH neutro (7.4), a correcao de protonacao e desprezivel.
        O ASO funciona normalmente no sangue.
        """
        corr = compute_protonation_correction(
            ph=7.4, n_cytosines=2, n_adenines=8,
        )
        # ddg_total deve ser negligivel (< 0.01 kcal/mol)
        assert corr.ddg_total < 0.01

    def test_acid_ph_positive_correction(self):
        """Em pH acido (4.5), ddg_total e positivo (desestabiliza o duplex).
        A protonacao de bases enfraquece pontes de hidrogenio Watson-Crick.
        """
        corr = compute_protonation_correction(
            ph=4.5, n_cytosines=2, n_adenines=8,
        )
        assert corr.ddg_total > 0

    def test_more_cytosines_more_penalty(self):
        """Mais citosinas = maior penalidade de protonacao.
        Citosina tem pKa mais alto (4.2) que adenina (3.5),
        portanto e mais protonada no fagolisossomo.
        """
        corr_2c = compute_protonation_correction(ph=4.5, n_cytosines=2, n_adenines=0)
        corr_5c = compute_protonation_correction(ph=4.5, n_cytosines=5, n_adenines=0)
        assert corr_5c.ddg_total > corr_2c.ddg_total

    def test_zero_bases_no_penalty(self):
        """Sem citosinas nem adeninas, nao ha penalidade de protonacao.
        (caso limite — sequencia hipotetica so de G e T)
        """
        corr = compute_protonation_correction(ph=4.5, n_cytosines=0, n_adenines=0)
        assert corr.ddg_total == 0.0

    def test_dg_at_ph_less_negative_in_acid(self):
        """dG corrigido para pH acido e MENOS negativo que em pH neutro.
        Menos negativo = ligacao mais fraca. O duplex perde estabilidade
        termodinamica no fagolisossomo.
        """
        dg_neutral = -28.0
        dg_acid = compute_dg_at_ph(
            dg_neutral=dg_neutral, ph=4.5, n_cytosines=2, n_adenines=8,
        )
        assert dg_acid > dg_neutral  # menos negativo

    def test_dg_still_functional_at_phagolysosome(self):
        """Para MRL-ASO-001 (2 C, 8 A), dG no fagolisossomo (pH 4.5)
        deve permanecer abaixo do limiar funcional (-15 kcal/mol).
        O ASO continua funcional mesmo no ambiente mais acido.
        """
        dg_acid = compute_dg_at_ph(
            dg_neutral=-28.0, ph=4.5, n_cytosines=2, n_adenines=8,
        )
        assert dg_acid < -15.0

    def test_tm_decreases_in_acid(self):
        """Tm corrigido para pH acido e MENOR que Tm neutro.
        Desestabilizacao termodinamica -> menor temperatura de fusao.
        """
        dg_neutral = -28.0
        dg_acid = compute_dg_at_ph(
            dg_neutral=dg_neutral, ph=4.5, n_cytosines=2, n_adenines=8,
        )
        tm_neutral = 108.5
        tm_acid = compute_tm_at_ph(
            tm_neutral=tm_neutral,
            dg_neutral=dg_neutral,
            dg_at_ph=dg_acid,
            n_bp=25,
        )
        assert tm_acid < tm_neutral

    def test_fraction_bound_strongly_negative_dg(self):
        """Com dG muito negativo (-28 kcal/mol), quase 100% do ASO esta ligado.
        A constante de equilibrio e tao favoravel que virtualmente toda
        molecula esta no complexo ASO:alvo.
        """
        f = compute_fraction_bound(-28.0)
        assert f > 0.99

    def test_fraction_bound_zero_dg(self):
        """Com dG = 0, exatamente metade esta ligada (equilibrio).
        Ponto de transicao: ligacao nao e favoravel nem desfavoravel.
        """
        f = compute_fraction_bound(0.0)
        assert f == pytest.approx(0.5, abs=0.01)

    def test_fraction_bound_positive_dg(self):
        """Com dG positivo, menos da metade esta ligada.
        Ligacao termodinamicamente desfavoravel — ASO nao funcionaria.
        """
        f = compute_fraction_bound(5.0)
        assert f < 0.5

    def test_fraction_bound_overflow_protection(self):
        """Valores extremos de dG nao devem causar overflow numerico.
        O modelo precisa ser robusto para dG -> -infinito e +infinito.
        """
        assert compute_fraction_bound(-1000.0) == 1.0
        assert compute_fraction_bound(1000.0) == 0.0


# ---------------------------------------------------------------------------
# 3. Cinetica de degradacao por nucleases
# ---------------------------------------------------------------------------


class TestNucleaseDegradation:
    """Testes para compute_half_life, compute_fraction_remaining e
    compute_nuclease_resistance."""

    def test_half_life_po_very_short(self):
        """Backbone PO nativo tem meia-vida MUITO curta em pH acido.
        Nucleases lisossomais degradam PO em horas — inutil para terapia.
        """
        t_half = compute_half_life(K_PO_PH45)
        # k=0.5/h -> t_half = ln(2)/0.5 = 1.39 h
        assert t_half == pytest.approx(math.log(2) / 0.5, abs=0.01)
        assert t_half < 2.0  # menos de 2 horas

    def test_half_life_ps_much_longer(self):
        """Backbone PS tem meia-vida 500x maior que PO.
        O enxofre impede coordenacao do metal catalitico da nuclease.
        """
        t_half_po = compute_half_life(K_PO_PH45)
        t_half_ps = compute_half_life(K_PS_PH45)
        assert t_half_ps > t_half_po * 400

    def test_half_life_zero_k_infinite(self):
        """Constante de degradacao zero = molecula infinitamente estavel.
        Caso limite para backbone teoricamente indestrutivel.
        """
        assert compute_half_life(0.0) == float("inf")

    def test_fraction_remaining_at_t0_is_one(self):
        """No tempo zero, 100% do ASO esta intacto (inicio do experimento)."""
        assert compute_fraction_remaining(K_PS_PH45, 0.0) == 1.0

    def test_fraction_remaining_decays_exponentially(self):
        """A fracao remanescente deve decair exponencialmente.
        Cinetica de primeira ordem: cada meia-vida reduz pela metade.
        """
        t_half = compute_half_life(K_PS_PH45)
        f_at_half = compute_fraction_remaining(K_PS_PH45, t_half)
        assert f_at_half == pytest.approx(0.5, abs=0.01)

    def test_nuclease_resistance_po_fails_therapeutic(self):
        """PO nao atinge a janela terapeutica (meia-vida < 24h).
        Backbone nativo e inadequado para entrega ao fagolisossomo.
        """
        profile = compute_nuclease_resistance("PO", K_PO_PH45)
        assert profile.therapeutic_window_met is False

    def test_nuclease_resistance_ps_meets_therapeutic(self):
        """PS atinge a janela terapeutica (meia-vida > 24h).
        O backbone fosforotioato e essencial para funcao no fagolisossomo.
        """
        profile = compute_nuclease_resistance("PS", K_PS_PH45)
        assert profile.therapeutic_window_met is True

    def test_gapmer_degradation_model(self):
        """O gapmer 5-15-5 tem constante efetiva entre PS e LNA pura.
        Flancos LNA protegem extremidades; gap DNA central e o elo limitante.
        """
        gapmer = compute_gapmer_degradation(
            lna_5prime=5, dna_gap=15, lna_3prime=5, k_ps=K_PS_PH45,
        )
        k_lna_expected = K_PS_PH45 * K_LNA_FACTOR
        assert gapmer.k_lna_flank == pytest.approx(k_lna_expected, rel=1e-4)
        # k_effective deve estar entre k_lna e k_ps (media ponderada)
        assert gapmer.k_effective < K_PS_PH45
        assert gapmer.k_effective > k_lna_expected
        # Meia-vida do gapmer > PS puro
        assert gapmer.half_life_hours > compute_half_life(K_PS_PH45)


# ---------------------------------------------------------------------------
# 4. Estabilidade conformacional LNA
# ---------------------------------------------------------------------------


class TestLNAConformation:
    """Testes para compute_lna_conformation."""

    def test_lna_c3_endo_maintained_at_any_ph(self):
        """LNA mantem conformacao C3'-endo independente do pH.
        A ponte metileno covalente nao e afetada por protonacao —
        vantagem estrutural fundamental do design gapmer.
        """
        profile = compute_lna_conformation(ph=4.5, n_lna=10, n_dna=15)
        assert profile.c3_endo_fraction_lna >= 0.99

    def test_dna_c3_endo_decreases_in_acid(self):
        """DNA livre perde conformacao C3'-endo em pH acido.
        Protonacao perturba o equilibrio N/S do acucar — DNA e
        mais vulneravel que LNA neste aspecto.
        """
        profile_neutral = compute_lna_conformation(ph=7.4, n_lna=0, n_dna=15)
        profile_acid = compute_lna_conformation(ph=4.5, n_lna=0, n_dna=15)
        assert profile_acid.c3_endo_fraction_dna < profile_neutral.c3_endo_fraction_dna

    def test_geometry_maintained_flag(self):
        """A flag geometry_maintained deve ser True quando >= 95% dos LNA
        estao em C3'-endo. Para MRL-ASO-001, isso e sempre verdade.
        """
        profile = compute_lna_conformation(ph=4.5, n_lna=10, n_dna=15)
        assert profile.geometry_maintained is True

    def test_dna_c3_endo_has_minimum_floor(self):
        """Mesmo em pH extremamente acido, a fracao C3'-endo do DNA
        nao cai abaixo de 10% (limite fisico do modelo).
        """
        profile = compute_lna_conformation(ph=1.0, n_lna=0, n_dna=15)
        assert profile.c3_endo_fraction_dna >= 0.10

    def test_neutral_ph_dna_c3_endo_baseline(self):
        """Em pH neutro (7.4), a fracao C3'-endo do DNA deve ser ~36%.
        Valor de referencia: Altona & Sundaralingam (1972).
        """
        profile = compute_lna_conformation(ph=7.4, n_lna=0, n_dna=15)
        assert profile.c3_endo_fraction_dna == pytest.approx(
            DNA_C3_ENDO_NEUTRAL, abs=0.01,
        )

    def test_penalty_increases_with_acidity(self):
        """A penalidade energetica aumenta conforme pH diminui.
        Mais protonacao -> mais perda de pre-organizacao -> maior custo.
        """
        profile_65 = compute_lna_conformation(ph=6.5, n_lna=10, n_dna=15)
        profile_45 = compute_lna_conformation(ph=4.5, n_lna=10, n_dna=15)
        assert profile_45.free_energy_penalty >= profile_65.free_energy_penalty
