"""Testes unitarios para aso_delivery/module_c_conjugate.

Valida os modelos de conjugados para entrega seletiva a macrofagos:
1. Perfis de receptores (MRC1, LDLR, ASGPR, SR-A)
2. Propriedades fisico-quimicas dos conjugados
3. Captacao estimada (fold vs ASO nu)
4. Score composto de recomendacao

A conclusao biologica: trimanose e o melhor conjugado para entrega
a macrofagos infectados por L. infantum (alta afinidade MRC1,
internalizacao para fagolisossomo). GalNAc e um ANTI-RESULTADO —
ASGPR e hepatocito-especifico, inutil para macrofagos.
"""

import math

import pytest

from aso_delivery.module_c_conjugate.conjugates import (
    MW_ASO_NAKED,
    LOGP_ASO_NAKED,
    MW_MANNOSE,
    MW_TRIMANNOSE,
    MW_C6_LINKER,
    MW_CHOLESTEROL,
    MW_PALMITATE,
    MW_GALNAC_TRIMER,
    KD_MANNOSE_MRC1_UM,
    KD_TRIMANNOSE_MRC1_UM,
    KD_GALNAC_ASGPR_UM,
    ASGPR_EXPRESSION_MACROPHAGE,
    MRC1_EXPRESSION_MACROPHAGE,
    SRA_EXPRESSION_MACROPHAGE,
    build_receptor_profiles,
    build_conjugate_properties,
    compute_uptake_fold_vs_naked,
    compute_conjugate_score,
    compute_naked_aso_score,
    ConjugateProperties,
    ReceptorProfile,
)


# ---------------------------------------------------------------------------
# 1. Perfis de receptores
# ---------------------------------------------------------------------------


class TestReceptorProfiles:
    """Testes para build_receptor_profiles e propriedades dos receptores."""

    def test_four_receptors_defined(self):
        """Deve haver exatamente 4 receptores no mapa: MRC1, LDLR, ASGPR, SR-A."""
        profiles = build_receptor_profiles()
        assert len(profiles) == 4
        assert set(profiles.keys()) == {"MRC1", "LDLR", "ASGPR", "SR-A"}

    def test_mrc1_high_macrophage_expression(self):
        """MRC1 (CD206) deve ter ALTA expressao em macrofagos.
        O receptor manose e marcador classico de macrofagos M2,
        que sao o tipo predominante na leishmaniose visceral.
        """
        profiles = build_receptor_profiles()
        assert profiles["MRC1"].expression_macrophage > 7.0

    def test_asgpr_zero_macrophage_expression(self):
        """ASGPR NAO e expresso em macrofagos (expressao = 0).
        Este receptor e EXCLUSIVO de hepatocitos — conjugar GalNAc
        ao ASO seria inutil para entrega a macrofagos. Este e o
        resultado negativo mais importante do Modulo C.
        """
        profiles = build_receptor_profiles()
        assert profiles["ASGPR"].expression_macrophage == 0.0

    def test_mrc1_selectivity_vs_hepatocyte(self):
        """MRC1 deve ter alta seletividade macrofago/hepatocito.
        Razao > 10 significa que para cada molecula captada por
        hepatocito, ~85 sao captadas por macrofago.
        """
        profiles = build_receptor_profiles()
        ratio = profiles["MRC1"].macrophage_selectivity_vs_hepatocyte
        assert ratio > 10.0

    def test_asgpr_selectivity_zero_or_inverse(self):
        """ASGPR com expressao 0 em macrofago: seletividade = 0 (nao macrofago).
        Divisao por hepatocito (9.5) com numerador 0 = 0.
        """
        profiles = build_receptor_profiles()
        ratio = profiles["ASGPR"].macrophage_selectivity_vs_hepatocyte
        assert ratio == 0.0

    def test_mrc1_endpoint_is_lysosome(self):
        """MRC1 internaliza para o lisossomo/fagolisossomo — exatamente
        onde L. infantum reside. O ASO seria entregue diretamente ao
        compartimento do parasita.
        """
        profiles = build_receptor_profiles()
        endpoint = profiles["MRC1"].endpoint_compartment.lower()
        assert "lysosome" in endpoint or "phagolysosome" in endpoint


# ---------------------------------------------------------------------------
# 2. Propriedades fisico-quimicas dos conjugados
# ---------------------------------------------------------------------------


class TestConjugateProperties:
    """Testes para build_conjugate_properties e computed properties."""

    def test_five_conjugates_defined(self):
        """Deve haver exatamente 5 conjugados avaliados."""
        conjugates = build_conjugate_properties()
        assert len(conjugates) == 5

    def test_mw_total_greater_than_naked(self):
        """Todos os conjugados devem ter MW > ASO nu.
        Conjugacao sempre adiciona massa (ligante + linker).
        """
        conjugates = build_conjugate_properties()
        for name, conj in conjugates.items():
            assert conj.mw_total_da > MW_ASO_NAKED, (
                f"{name}: MW total ({conj.mw_total_da}) <= ASO nu ({MW_ASO_NAKED})"
            )

    def test_cholesterol_increases_logp(self):
        """Colesterol e altamente lipofilico (logP ~7).
        O conjugado colesterol-ASO deve ter logP mais alto (menos negativo)
        que o ASO nu, aumentando associacao a membrana.
        """
        conjugates = build_conjugate_properties()
        chol = conjugates["cholesterol"]
        assert chol.logp_estimated > LOGP_ASO_NAKED

    def test_mannose_keeps_logp_low(self):
        """Manose e muito hidrofilica (logP ~-2.7).
        Conjugado manose-ASO deve ter logP similar ao ASO nu —
        a seletividade vem da via receptora, nao da hidrofobicidade.
        """
        conjugates = build_conjugate_properties()
        mannose = conjugates["mannose_mono"]
        # logP nao deve mudar drasticamente
        assert abs(mannose.logp_estimated - LOGP_ASO_NAKED) < 2.0

    def test_membrane_association_cholesterol_highest(self):
        """Colesterol deve ter a MAIOR constante de associacao a membrana.
        Colesterol e componente natural da bicamada lipidica — insercao
        espontanea. Isso confere tropismo amplo (nao seletivo).
        """
        conjugates = build_conjugate_properties()
        ka_chol = conjugates["cholesterol"].membrane_association_constant
        ka_mannose = conjugates["mannose_mono"].membrane_association_constant
        assert ka_chol > ka_mannose

    def test_mannose_mw_total_correct(self):
        """MW do conjugado manose-ASO deve ser ASO + manose + linker C6."""
        conjugates = build_conjugate_properties()
        mannose = conjugates["mannose_mono"]
        expected = MW_ASO_NAKED + MW_MANNOSE + MW_C6_LINKER
        assert mannose.mw_total_da == pytest.approx(expected, abs=0.1)

    def test_trimannose_kd_much_lower_than_mono(self):
        """Trimanose tem Kd 50x menor que manose monovalente (efeito avidez).
        Multivalencia multiplica a afinidade geometricamente.
        """
        conjugates = build_conjugate_properties()
        kd_mono = conjugates["mannose_mono"].kd_receptor_um
        kd_tri = conjugates["trimannose_cluster"].kd_receptor_um
        assert kd_tri < kd_mono / 10.0


# ---------------------------------------------------------------------------
# 3. Captacao estimada (fold vs naked)
# ---------------------------------------------------------------------------


class TestUptakeFold:
    """Testes para compute_uptake_fold_vs_naked."""

    def test_no_receptor_minimal_uptake(self):
        """Quando receptor nao e expresso (expressao = 0), nao ha ganho
        receptor-mediado. A captacao e residual (~0.1x).
        Exemplo: GalNAc em macrofago (ASGPR ausente).
        """
        fold = compute_uptake_fold_vs_naked(
            kd_um=0.003,  # GalNAc tem afinidade altissima...
            receptor_expression=0.0,  # ...mas receptor nao existe no macrofago
            membrane_ka=1.0,
        )
        # Sem receptor, captacao e muito baixa
        assert fold < 1.0

    def test_high_affinity_high_expression_good_uptake(self):
        """Alta afinidade + alta expressao = alta captacao.
        Trimanose + MRC1 em macrofagos: o melhor cenario.
        """
        fold = compute_uptake_fold_vs_naked(
            kd_um=KD_TRIMANNOSE_MRC1_UM,  # 0.1 uM
            receptor_expression=MRC1_EXPRESSION_MACROPHAGE,  # 8.5
            membrane_ka=1.0,
        )
        assert fold > 5.0  # melhoria significativa

    def test_lower_kd_better_uptake(self):
        """Menor Kd = maior afinidade = maior captacao.
        Comparando trimanose (0.1 uM) vs manose mono (5 uM).
        """
        fold_mono = compute_uptake_fold_vs_naked(
            kd_um=KD_MANNOSE_MRC1_UM,
            receptor_expression=MRC1_EXPRESSION_MACROPHAGE,
            membrane_ka=1.0,
        )
        fold_tri = compute_uptake_fold_vs_naked(
            kd_um=KD_TRIMANNOSE_MRC1_UM,
            receptor_expression=MRC1_EXPRESSION_MACROPHAGE,
            membrane_ka=1.0,
        )
        assert fold_tri > fold_mono

    def test_higher_expression_better_uptake(self):
        """Mais receptores na superficie = mais captacao.
        Macrofago vs fibroblasto para o mesmo conjugado.
        """
        fold_high = compute_uptake_fold_vs_naked(
            kd_um=5.0, receptor_expression=8.5, membrane_ka=1.0,
        )
        fold_low = compute_uptake_fold_vs_naked(
            kd_um=5.0, receptor_expression=1.0, membrane_ka=1.0,
        )
        assert fold_high > fold_low

    def test_membrane_ka_contributes(self):
        """Maior associacao a membrana melhora a captacao (efeito menor).
        Conjugados lipofilicos (colesterol) tem Ka de membrana elevado.
        """
        fold_low_ka = compute_uptake_fold_vs_naked(
            kd_um=2.0, receptor_expression=3.0, membrane_ka=0.5,
        )
        fold_high_ka = compute_uptake_fold_vs_naked(
            kd_um=2.0, receptor_expression=3.0, membrane_ka=5.0,
        )
        assert fold_high_ka > fold_low_ka


# ---------------------------------------------------------------------------
# 4. Score composto de recomendacao
# ---------------------------------------------------------------------------


class TestConjugateScore:
    """Testes para compute_conjugate_score e compute_naked_aso_score."""

    def _get_receptor(self, name: str) -> ReceptorProfile:
        """Helper: obtem perfil de um receptor pelo nome."""
        return build_receptor_profiles()[name]

    def _get_conjugate(self, name: str) -> ConjugateProperties:
        """Helper: obtem propriedades de um conjugado pelo nome."""
        return build_conjugate_properties()[name]

    def test_naked_aso_score_baseline(self):
        """O ASO nu deve ter score de referencia (recommendation ~0.42).
        Define o baseline para comparacao com conjugados.
        """
        score = compute_naked_aso_score()
        assert score.suitable_for_macrophage is True
        assert 0.3 < score.recommendation_score < 0.6

    def test_galnac_unsuitable_for_macrophage(self):
        """GalNAc NAO deve ser recomendado para macrofagos.
        ASGPR nao e expresso em macrofagos — o conjugado mais eficiente
        para hepatocitos e completamente inutil para o nosso alvo.
        """
        conj = self._get_conjugate("galnac")
        receptor = self._get_receptor("ASGPR")
        fold = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        score = compute_conjugate_score(conj, receptor, fold, max_uptake_fold=20.0)
        assert score.suitable_for_macrophage is False
        assert score.expression_score == 0.0

    def test_trimannose_suitable_for_macrophage(self):
        """Trimanose DEVE ser recomendado para macrofagos.
        Alta expressao de MRC1 + alta afinidade + entrega ao fagolisossomo.
        """
        conj = self._get_conjugate("trimannose_cluster")
        receptor = self._get_receptor("MRC1")
        fold = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        score = compute_conjugate_score(conj, receptor, fold, max_uptake_fold=20.0)
        assert score.suitable_for_macrophage is True
        assert score.pathway_score == 1.0  # fagolisossomo

    def test_trimannose_scores_higher_than_naked(self):
        """Trimanose deve ter score total maior que ASO nu.
        O conjugado melhora a entrega em todos os aspectos relevantes.
        """
        conj = self._get_conjugate("trimannose_cluster")
        receptor = self._get_receptor("MRC1")
        fold = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        score_tri = compute_conjugate_score(conj, receptor, fold, max_uptake_fold=20.0)
        score_naked = compute_naked_aso_score()
        assert score_tri.recommendation_score > score_naked.recommendation_score

    def test_score_components_bounded_zero_to_one(self):
        """Todos os componentes do score devem estar em [0, 1].
        Normalizacao garante comparabilidade entre conjugados.
        """
        conj = self._get_conjugate("mannose_mono")
        receptor = self._get_receptor("MRC1")
        fold = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        score = compute_conjugate_score(conj, receptor, fold, max_uptake_fold=20.0)
        assert 0.0 <= score.receptor_score <= 1.0
        assert 0.0 <= score.expression_score <= 1.0
        assert 0.0 <= score.selectivity_score <= 1.0
        assert 0.0 <= score.pathway_score <= 1.0
        assert 0.0 <= score.synthesis_feasibility <= 1.0
        assert 0.0 <= score.cost_score <= 1.0
        assert 0.0 <= score.recommendation_score <= 1.0

    def test_galnac_pathway_score_zero(self):
        """GalNAc entrega ao hepatocito, nao ao macrofago.
        O pathway_score deve ser 0 porque ASGPR e hepatocito-especifico.
        """
        conj = self._get_conjugate("galnac")
        receptor = self._get_receptor("ASGPR")
        fold = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        score = compute_conjugate_score(conj, receptor, fold, max_uptake_fold=20.0)
        assert score.pathway_score == 0.0

    def test_palmitate_suitable_via_sra(self):
        """Palmitato deve ser recomendado para macrofagos via SR-A.
        Scavenger receptor classe A tem alta expressao em macrofagos.
        Precedente historico: Toulme et al. (1994).
        """
        conj = self._get_conjugate("palmitate")
        receptor = self._get_receptor("SR-A")
        fold = compute_uptake_fold_vs_naked(
            kd_um=conj.kd_receptor_um,
            receptor_expression=conj.receptor_expression,
            membrane_ka=conj.membrane_association_constant,
        )
        score = compute_conjugate_score(conj, receptor, fold, max_uptake_fold=20.0)
        assert score.suitable_for_macrophage is True
