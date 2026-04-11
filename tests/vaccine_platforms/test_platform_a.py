"""Testes para a Plataforma A -- Vacina de mRNA otimizada.

Testa os sub-modulos A1 (analise de epitopos), A2 (pseudouridina),
A3 (formulacao LNP) e A4 (analise estrategica) a nivel unitario,
sem executar o pipeline completo (que depende de construct_card.json).
"""

import pytest

from vaccine_platforms.shared.epitopes import (
    ADJUVANT_SEQUENCE,
    EPITOPES,
    get_construct_sequence,
    get_epitope_sequences,
    get_unique_source_genes,
)
from vaccine_platforms.shared.validators import (
    calculate_cai,
    calculate_gc_content,
    calculate_instability_index,
    calculate_molecular_weight,
)
from vaccine_platforms.platform_a_mrna import epitope_analysis as A1


# ===========================================================================
# A1: Analise de redundancia de epitopos
# ===========================================================================


class TestEpitopeRedundancy:
    """Testa analise de redundancia (A1) usando os epitopos da SSOT."""

    def test_run_redundancy_analysis_retorna_resultado(self):
        """A analise de redundancia deve completar sem erro.

        Executa a comparacao par-a-par dos 11 epitopos para identificar
        pares redundantes (>70% similaridade + mesmo alelo DLA).
        """
        result = A1.run_redundancy_analysis()
        assert result is not None
        assert result.total_epitopes == 11

    def test_contagem_pares(self):
        """Deve gerar C(11,2) = 55 comparacoes par-a-par.

        Cada par de epitopos e comparado exatamente uma vez (combinacao,
        nao permutacao).
        """
        result = A1.run_redundancy_analysis()
        assert len(result.pairwise_comparisons) == 55

    def test_genes_unicos(self):
        """Deve identificar pelo menos 4 genes fonte distintos.

        Os epitopos vem de: emp24/GOLD, hipotetica conservada, proibina,
        GP63/leishmanolysin, e CPB. Cobertura antigenica robusta.
        """
        result = A1.run_redundancy_analysis()
        assert result.unique_genes >= 4

    def test_alelos_unicos(self):
        """Deve identificar pelo menos 2 alelos DLA distintos.

        Cobertura alelica e importante: caes tem diversidade de MHC
        significativa. Quanto mais alelos cobertos, maior a fracao da
        populacao canina protegida.
        """
        result = A1.run_redundancy_analysis()
        assert result.unique_alleles >= 2

    def test_minimum_epitopes_para_cobertura(self):
        """Numero minimo de epitopos para 95% de cobertura deve ser <= 11.

        Se 11 epitopos cobrem 100%, o minimo para 95% sera <= 11 por definicao.
        """
        result = A1.run_redundancy_analysis()
        assert result.minimum_epitopes_95pct <= 11
        assert result.minimum_epitopes_95pct >= 1


# ===========================================================================
# Propriedades biofisicas do construto mRNA
# ===========================================================================


class TestConstructBiophysics:
    """Testa propriedades biofisicas do construto vacinal completo.

    O construto multi-epitopo tem ~335 aa e deve ter propriedades
    compativeis com expressao, estabilidade e imunogenicidade.
    """

    def test_construto_comprimento_razoavel(self):
        """Construto deve ter entre 250 e 500 aa.

        Muito curto: imunogenicidade insuficiente.
        Muito longo: dificuldade de producao e folding.
        """
        construct = get_construct_sequence(include_signal=True)
        length = len(construct)
        assert 250 <= length <= 500, f"Comprimento inesperado: {length} aa"

    def test_peso_molecular_razoavel(self):
        """Peso molecular deve estar entre 25 e 60 kDa.

        Proteinas vacinais tipicamente ficam nesta faixa para boa
        expressao em sistemas de producao recombinante.
        """
        construct = get_construct_sequence(include_signal=True)
        mw = calculate_molecular_weight(construct)
        mw_kda = mw / 1000.0
        assert 25 <= mw_kda <= 60, f"MW inesperado: {mw_kda:.1f} kDa"

    def test_instabilidade_calculavel(self):
        """Indice de instabilidade deve ser calculavel e numerico.

        Nao validamos se e estavel ou instavel aqui -- apenas se o calculo
        funciona corretamente com a sequencia completa.
        """
        construct = get_construct_sequence(include_signal=True)
        ii = calculate_instability_index(construct)
        assert isinstance(ii, float)
        # Indice tipico varia de -40 a +80
        assert -100 < ii < 100

    def test_adjuvante_na_posicao_correta(self):
        """O adjuvante L7/L12 deve vir apos o peptideo sinal no construto.

        Arquitetura esperada: [tPA signal]-[L7/L12]-EAAAK-[epitopos].
        O adjuvante nao deve estar no inicio absoluto (sinal vem antes).
        """
        construct = get_construct_sequence(include_signal=True)
        adj_pos = construct.find(ADJUVANT_SEQUENCE)
        assert adj_pos > 0, "Adjuvante nao encontrado ou na posicao 0"

    def test_todos_epitopos_presentes(self):
        """Todos os 11 epitopos unicos devem estar presentes no construto.

        A integridade do cassete de epitopos e fundamental -- perder
        epitopos significaria reduzir a cobertura antigenica.
        """
        construct = get_construct_sequence(include_signal=True)
        for ep in EPITOPES:
            assert ep.peptide in construct, (
                f"Epitopo {ep.peptide} ausente do construto"
            )
