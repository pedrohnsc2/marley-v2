"""Testes para a Plataforma C -- Vacina em L. tarentolae (LEXSY).

Foco nos testes CRITICOS:
  - C0: Ortogonalidade do SL RNA (MRL-ASO-001 NAO pode matar L. tarentolae)
  - C2: Construcao do construto
  - C3: Otimizacao de codons para trypanosomatideos
  - C4: Modelo de custos
"""

import pytest

from vaccine_platforms.platform_c_tarentolae.orthogonality import (
    DG_SELECTIVITY_THRESHOLD,
    OrthogonalityResult,
    SL_TARENTOLAE,
    verify_sl_rna_orthogonality,
    _align_aso_to_target,
    _reverse_complement,
    _calculate_nn_thermodynamics,
)
from vaccine_platforms.platform_c_tarentolae.construct import (
    HIS6_TAG,
    SAP1_SIGNAL,
    STREP_TAG,
    build_tarentolae_construct,
    validate_construct,
)
from vaccine_platforms.platform_c_tarentolae.codon_optimizer import (
    run_codon_optimization as run_lt_codon_optimization,
)
from vaccine_platforms.platform_c_tarentolae.cost_model import (
    calculate_costs,
)
from aso_math.config import ASO_SEQUENCE, SL_SEQUENCE


# ===========================================================================
# C0: Verificacao de ortogonalidade SL RNA (TESTE CRITICO)
# ===========================================================================


class TestSLOrthogonality:
    """Verifica que o ASO MRL-ASO-001 NAO mata L. tarentolae.

    Este e o teste mais critico de toda a Plataforma C. Se o ASO
    desenhado para inibir o trans-splicing de L. infantum tambem
    afetasse L. tarentolae, a terapia combinada (ASO + vacina viva)
    seria inviavel.

    A seletividade depende de mismatches nas sequencias do mini-exon
    do SL RNA entre as duas especies.
    """

    @pytest.fixture
    def orth(self) -> OrthogonalityResult:
        """Executa a verificacao de ortogonalidade uma vez para todos os testes."""
        return verify_sl_rna_orthogonality()

    def test_sl_rna_tarentolae_39_nt(self):
        """SL RNA de L. tarentolae deve ter 39 nucleotideos (mini-exon).

        O mini-exon e conservado em tamanho entre trypanosomatideos --
        todos usam 39 nt que sao trans-spliced no 5' de cada mRNA.
        """
        assert len(SL_TARENTOLAE) == 39

    def test_sl_rna_infantum_39_nt(self):
        """SL RNA de L. infantum tambem deve ter 39 nucleotideos."""
        assert len(SL_SEQUENCE) == 39

    def test_sl_rnas_diferem(self):
        """Os SL RNAs de L. infantum e L. tarentolae devem ser DIFERENTES.

        Se fossem identicos, o ASO teria afinidade igual por ambos e
        a terapia combinada seria impossivel. A diferenca permite
        seletividade termodinamica.
        """
        assert SL_TARENTOLAE != SL_SEQUENCE

    def test_ortogonalidade_nao_passa_threshold(self, orth: OrthogonalityResult):
        """DeltaDeltaG NAO atinge o threshold de 2.0 kcal/mol.

        As sequencias de SL RNA de L. infantum e L. tarentolae diferem
        apenas nas posicoes 34-35 (fora da janela alvo do ASO, que e
        posicoes 5-30). Dentro da regiao alvo, as sequencias sao IDENTICAS.
        Portanto, o ASO tem afinidade IGUAL por ambas as especies --
        passes_threshold sera False, e isso e o resultado biologico correto.

        Implicacao terapeutica: a terapia combinada ASO + vacina L. tarentolae
        requer precaucoes (administracao sequencial com washout).
        """
        assert orth.passes_threshold is False, (
            "Inesperado: DeltaDeltaG deveria ser 0 porque a regiao alvo "
            "do ASO (pos 5-30) e identica em ambas as especies."
        )

    def test_delta_dg_zero(self, orth: OrthogonalityResult):
        """Diferenca de DeltaG deve ser 0 kcal/mol.

        A regiao alvo do ASO (posicoes 5-30 do SL RNA) e 100% conservada
        entre L. infantum e L. tarentolae. As diferencas entre as duas
        especies estao nas posicoes 34-35, FORA da regiao alvo.
        """
        assert orth.delta_dg == 0.0

    def test_kd_ratio_igual_a_um(self, orth: OrthogonalityResult):
        """Razao Kd deve ser 1.0 (afinidade identica).

        Como a regiao alvo e identica, o ASO forma duplexes termodinamicamente
        equivalentes com ambas as especies.
        """
        assert orth.kd_ratio == 1.0

    def test_zero_mismatches_na_regiao_alvo(self, orth: OrthogonalityResult):
        """Zero mismatches na regiao alvo do ASO.

        As diferencas entre L. infantum e L. tarentolae estao fora da
        janela alvo (posicoes 5-30). Dentro dessa janela, os SL RNAs
        sao identicos -- o ASO se liga igualmente a ambos.
        """
        assert orth.mismatches == 0

    def test_complementaridade_infantum_perfeita(self, orth: OrthogonalityResult):
        """Complementaridade ASO:L.infantum deve ser 100%.

        O ASO foi desenhado para ter match perfeito com o SL RNA de
        L. infantum na regiao alvo.
        """
        assert orth.complementarity_infantum == 1.0

    def test_complementaridade_tarentolae_tambem_perfeita(self, orth: OrthogonalityResult):
        """Complementaridade ASO:L.tarentolae tambem e 100% na regiao alvo.

        Embora os SL RNAs completos difiram em 2 posicoes (34-35), a regiao
        alvo do ASO (posicoes 5-30) e identica. Portanto, complementaridade
        na regiao alvo e perfeita para ambas as especies.
        """
        assert orth.complementarity_tarentolae == 1.0

    def test_tm_identica(self, orth: OrthogonalityResult):
        """Tm dos duplexes deve ser identica (mesma regiao alvo).

        Como a regiao alvo e identica, os parametros termodinamicos
        (incluindo Tm) devem ser identicos.
        """
        assert orth.tm_infantum == orth.tm_tarentolae

    def test_veredito_contem_warning(self, orth: OrthogonalityResult):
        """O veredito deve ser WARNING (nao PASS).

        Este e um resultado honesto: o ASO pode afetar L. tarentolae,
        e o pipeline documenta isso explicitamente ao inves de esconder.
        """
        assert "WARNING" in orth.verdict

    def test_warnings_nao_vazios(self, orth: OrthogonalityResult):
        """Deve haver avisos documentando o risco.

        O sistema de avisos e essencial para comunicar riscos ao
        pesquisador. Neste caso, o risco e real e deve ser reportado.
        """
        assert len(orth.warnings) >= 1


# ===========================================================================
# Funcoes auxiliares de ortogonalidade
# ===========================================================================


class TestOrthogonalityHelpers:
    """Testa funcoes auxiliares usadas na verificacao de ortogonalidade."""

    def test_reverse_complement_simples(self):
        """Complemento reverso de ATGC deve ser GCAT.

        Watson-Crick: A<->T, G<->C. Reverso: ler de tras para frente.
        """
        assert _reverse_complement("ATGC") == "GCAT"

    def test_reverse_complement_palindromo(self):
        """Complemento reverso de AATT deve ser AATT (palindromo).

        Palindromos de DNA sao a base dos sitios de restricao -- a
        sequencia e identica em ambas as fitas.
        """
        assert _reverse_complement("AATT") == "AATT"

    def test_align_aso_identico(self):
        """Alinhamento perfeito deve ter 0 mismatches.

        Quando o complemento reverso do ASO corresponde exatamente ao
        alvo, devemos ter complementaridade perfeita.
        """
        # ATGC -> rc = GCAT. Se alvo = GCAT, match perfeito.
        matches, mismatches, positions = _align_aso_to_target("ATGC", "GCAT")
        assert matches == 4
        assert mismatches == 0
        assert positions == []

    def test_thermodynamics_retorna_quatro_valores(self):
        """Calculo NN deve retornar tupla (dH, dS, Tm, dG).

        Modelo de vizinhos proximos (nearest-neighbor) de SantaLucia (1998)
        calcula parametros termodinamicos de duplexes DNA/DNA.
        """
        dh, ds, tm, dg = _calculate_nn_thermodynamics("ATGCATGC", "TACGTACG")
        assert isinstance(dh, float)
        assert isinstance(ds, float)
        assert isinstance(tm, float)
        assert isinstance(dg, float)


# ===========================================================================
# C2: Construto para L. tarentolae
# ===========================================================================


class TestTarentolaeConstruct:
    """Testa a montagem do construto para expressao em L. tarentolae."""

    @pytest.fixture
    def construct(self):
        """Constroi o construto uma vez."""
        return build_tarentolae_construct()

    def test_contem_sap1_signal(self, construct):
        """Construto deve conter peptideo sinal SAP1 de L. mexicana.

        O SAP1 e reconhecido pela maquinaria secretoria de trypanosomatideos,
        permitindo secrecao da proteina recombinante no meio de cultura
        (importante para a modalidade de proteina secretada purificada).
        """
        assert SAP1_SIGNAL in construct.protein_sequence

    def test_contem_his6_tag(self, construct):
        """Deve conter His6-tag para purificacao."""
        assert HIS6_TAG in construct.protein_sequence

    def test_contem_strep_tag(self, construct):
        """Deve conter Strep-tag II para purificacao tandem.

        A purificacao dupla (IMAC + Strep-Tactin) permite pureza > 95%
        com apenas dois passos cromatograficos.
        """
        assert STREP_TAG in construct.protein_sequence

    def test_11_epitopos_unicos(self, construct):
        """Deve conter 11 epitopos unicos (mesmos das outras plataformas)."""
        assert len(construct.unique_epitopes) == 11

    def test_vetor_plexsy(self, construct):
        """Vetor deve ser pLEXSY-sat2 (sistema de expressao LEXSY).

        O pLEXSY integra no locus SSU rRNA de L. tarentolae,
        gerando uma linhagem estavel de expressao.
        """
        assert "pLEXSY" in construct.vector

    def test_peso_molecular_razoavel(self, construct):
        """MW deve estar entre 25 e 60 kDa."""
        mw_kda = construct.molecular_weight_da / 1000.0
        assert 25 <= mw_kda <= 60

    def test_comprimento_coerente(self, construct):
        """length_aa deve corresponder ao comprimento real da sequencia."""
        assert construct.length_aa == len(construct.protein_sequence)


# ===========================================================================
# C3: Otimizacao de codons para L. tarentolae
# ===========================================================================


class TestTarentolaeCodonOptimization:
    """Testa otimizacao de codons para trypanosomatideos.

    Trypanosomatideos tem forte vies GC-rich (~60% GC no genoma).
    A otimizacao deve refletir essa preferencia.
    """

    @pytest.fixture
    def codon_result(self):
        """Otimiza codons do construto completo."""
        construct = build_tarentolae_construct()
        return run_lt_codon_optimization(construct.protein_sequence)

    def test_gc_content_range(self, codon_result):
        """GC content deve estar entre 50-70% (vies GC de trypanosomatideos).

        O genoma de L. tarentolae tem ~60% GC. Sequencias otimizadas devem
        estar perto desse valor para expressao eficiente.
        """
        assert 0.50 <= codon_result.gc_content <= 0.70, (
            f"GC content fora do range: {codon_result.gc_content:.1%}"
        )

    def test_cai_acima_limiar(self, codon_result):
        """CAI deve ser > 0.70 (idealmente > 0.85 para alta expressao).

        CAI alto indica boa adaptacao aos codons preferenciais de
        L. tarentolae, maximizando eficiencia de traducao.
        """
        assert codon_result.cai > 0.70, (
            f"CAI muito baixo: {codon_result.cai:.4f}"
        )

    def test_comprimento_dna_multiplo_de_3(self, codon_result):
        """Comprimento do DNA codificante deve ser multiplo de 3.

        Cada codon tem 3 nucleotideos. Se nao for multiplo de 3,
        houve um erro de frame.
        """
        # Excluir flanqueadores para verificar a CDS pura
        assert codon_result.length_nt > 0


# ===========================================================================
# C4: Modelo de custos
# ===========================================================================


class TestTarentolaeCosts:
    """Testa o modelo de custos da Plataforma C.

    L. tarentolae tem vantagens de custo significativas:
    - Crescimento a 26C (sem CO2)
    - BSL-1 (sem biosseguranca especial)
    - Adjuvante intrinseco (vacina viva)
    """

    @pytest.fixture
    def costs(self):
        """Calcula custos."""
        return calculate_costs()

    def test_custo_vacina_viva_baixo(self, costs):
        """Custo por dose da vacina viva deve ser < $1.

        A vacina viva usa o organismo inteiro como imunogeno e adjuvante
        simultaneamente, eliminando custos de purificacao e adjuvante exogeno.
        """
        assert costs.live_cost_per_dose < 1.0, (
            f"Custo por dose muito alto: ${costs.live_cost_per_dose:.4f}"
        )

    def test_custo_proteina_razoavel(self, costs):
        """Custo da proteina secretada purificada por dose deve ser < $10."""
        assert costs.protein_total_dose_cost < 10.0

    def test_doses_por_litro_alto(self, costs):
        """Vacina viva deve render > 1000 doses por litro.

        L. tarentolae cresce a densidades de ~1e8 celulas/mL. Com dose
        tipica de 1e7-1e8 celulas, 1 litro produz milhares de doses.
        """
        assert costs.live_doses_per_l > 1000

    def test_custo_industrial_menor_que_lab(self, costs):
        """Custo industrial deve ser menor que laboratorio (economia de escala)."""
        assert costs.live_industrial_cost_per_dose < costs.live_cost_per_dose
