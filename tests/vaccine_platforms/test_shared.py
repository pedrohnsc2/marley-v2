"""Testes para os modulos compartilhados: epitopes.py e validators.py.

Garante a integridade da SSOT (Single Source of Truth) dos epitopos e
a corretude dos validadores biofisicos usados por todas as plataformas.
"""

import pytest

from vaccine_platforms.shared.epitopes import (
    ADJUVANT_EPITOPE_LINKER,
    ADJUVANT_SEQUENCE,
    EPITOPES,
    LINKERS,
    SIGNAL_PEPTIDE_TPA,
    Epitope,
    get_construct_sequence,
    get_epitope_by_allele,
    get_epitope_by_gene,
    get_epitope_sequences,
    get_unique_source_genes,
)
from vaccine_platforms.shared.validators import (
    calculate_aromaticity,
    calculate_cai,
    calculate_gc_content,
    calculate_gravy,
    calculate_instability_index,
    calculate_isoelectric_point,
    calculate_molecular_weight,
    check_restriction_sites,
)


# ===========================================================================
# Testes de SSOT: epitopos imutaveis
# ===========================================================================


class TestEpitopesSST:
    """Verifica a integridade da tuple imutavel de epitopos canonicos."""

    def test_exatamente_11_epitopos(self):
        """Devem existir exatamente 11 epitopos unicos na SSOT.

        Motivo biologico: o pipeline de vacinologia reversa selecionou 11
        peptideos unicos com IC50 < 200 nM para alelos DLA caninos.
        """
        assert len(EPITOPES) == 11

    def test_epitopos_tupla_imutavel(self):
        """A SSOT deve ser uma tupla (imutavel), nao uma lista.

        Impedir modificacoes acidentais e critico -- se alguem adicionar ou
        remover epitopos, todas as plataformas seriam afetadas silenciosamente.
        """
        assert isinstance(EPITOPES, tuple)

    def test_epitopo_frozen_dataclass(self):
        """Cada epitopo deve ser um dataclass frozen (imutavel).

        Frozen impede que se altere atributos apos criacao, garantindo que
        os dados vindos do pipeline de vacinologia reversa nao sejam corrompidos.
        """
        ep = EPITOPES[0]
        assert isinstance(ep, Epitope)
        with pytest.raises(AttributeError):
            ep.peptide = "XXXXXXXXX"  # type: ignore[misc]

    def test_peptideos_sao_strings_validas(self):
        """Todos os peptideos devem conter apenas aminoacidos validos.

        Aminoacidos padrao do codigo genetico (20 letras). Epitopos CTL
        sao tipicamente 8-11 meros para apresentacao via MHC classe I.
        """
        aminoacidos_validos = set("ACDEFGHIKLMNPQRSTVWY")
        for ep in EPITOPES:
            assert len(ep.peptide) >= 8, (
                f"Epitopo muito curto: {ep.peptide} ({len(ep.peptide)} aa)"
            )
            assert len(ep.peptide) <= 15, (
                f"Epitopo muito longo: {ep.peptide} ({len(ep.peptide)} aa)"
            )
            assert set(ep.peptide).issubset(aminoacidos_validos), (
                f"Caracteres invalidos em {ep.peptide}"
            )

    def test_ic50_positivo_e_menor_500(self):
        """IC50 de todos os epitopos deve ser positivo e < 500 nM.

        IC50 < 500 nM e o limiar padrao para ligantes de MHC classe I
        no NetMHCpan. Nossos epitopos sao todos < 200 nM (forte).
        """
        for ep in EPITOPES:
            assert 0 < ep.ic50 < 500, (
                f"IC50 invalido para {ep.peptide}: {ep.ic50}"
            )

    def test_get_epitope_sequences_retorna_11(self):
        """get_epitope_sequences() deve retornar lista com 11 peptideos."""
        seqs = get_epitope_sequences()
        assert isinstance(seqs, list)
        assert len(seqs) == 11

    def test_get_unique_source_genes(self):
        """Deve agrupar epitopos em pelo menos 4 genes distintos.

        Cobertura antigenica: 4+ proteinas diferentes do parasita garante
        que resposta imune atinge multiplos alvos (redundancia protetora).
        """
        genes = get_unique_source_genes()
        assert len(genes) >= 4

    def test_get_epitope_by_gene(self):
        """Busca por gene LINF_240013900 deve retornar 3 epitopos (emp24/GOLD)."""
        eps = get_epitope_by_gene("LINF_240013900")
        assert len(eps) == 3
        for ep in eps:
            assert ep.gene_id == "LINF_240013900"

    def test_get_epitope_by_allele(self):
        """Busca por alelo DLA-8850101 deve retornar epitopos validos."""
        eps = get_epitope_by_allele("DLA-8850101")
        assert len(eps) >= 1
        for ep in eps:
            assert ep.allele == "DLA-8850101"

    def test_construct_sequence_inclui_adjuvante(self):
        """Construto vacinal deve conter a sequencia do adjuvante L7/L12.

        O adjuvante L7/L12 ribosomal bacteriano ativa TLR4, induzindo
        resposta Th1 -- essencial contra Leishmania (patogeno intracelular).
        """
        construct = get_construct_sequence()
        assert ADJUVANT_SEQUENCE in construct

    def test_construct_sequence_inclui_sinal_tpa(self):
        """Construto com sinal deve iniciar com peptideo tPA.

        O tPA direciona a proteina para a via secretoria, aumentando a
        apresentacao via MHC-I por cross-priming.
        """
        construct = get_construct_sequence(include_signal=True)
        assert construct.startswith(SIGNAL_PEPTIDE_TPA)

    def test_construct_sequence_sem_sinal(self):
        """Construto sem sinal nao deve comecar com tPA."""
        construct = get_construct_sequence(include_signal=False)
        assert not construct.startswith(SIGNAL_PEPTIDE_TPA)

    def test_construct_cpb_repeats_invalido(self):
        """cpb_repeats fora do range 1-5 deve levantar ValueError."""
        with pytest.raises(ValueError):
            get_construct_sequence(cpb_repeats=0)
        with pytest.raises(ValueError):
            get_construct_sequence(cpb_repeats=6)


# ===========================================================================
# Testes dos validadores biofisicos
# ===========================================================================


class TestValidatorsMolecularWeight:
    """Testes para o calculo de peso molecular de proteinas."""

    def test_mw_retorna_float(self):
        """Peso molecular deve ser um float positivo."""
        mw = calculate_molecular_weight("MKLLVV")
        assert isinstance(mw, float)
        assert mw > 0

    def test_mw_sequencia_curta_conhecida(self):
        """Peso de um tripeptideo simples (GGG) deve ser razoavel.

        Glicina: 75.032 Da. Para GGG: 3*75.032 - 2*18.015 = 189.066 Da.
        """
        mw = calculate_molecular_weight("GGG")
        assert abs(mw - 189.07) < 0.1

    def test_mw_sequencia_vazia(self):
        """Sequencia vazia deve retornar 0."""
        mw = calculate_molecular_weight("")
        assert mw == 0.0

    def test_mw_caractere_invalido(self):
        """Aminoacido invalido (X, B, Z, etc.) deve levantar ValueError."""
        with pytest.raises(ValueError):
            calculate_molecular_weight("MKXLVV")


class TestValidatorsIsoelectricPoint:
    """Testes para o calculo do ponto isoeletrico (pI)."""

    def test_pi_retorna_float(self):
        """pI deve ser um float entre 0 e 14."""
        pi = calculate_isoelectric_point("MKLLVV")
        assert isinstance(pi, float)
        assert 0 < pi < 14

    def test_pi_peptideo_basico(self):
        """Peptideo rico em Lys/Arg deve ter pI alto (basico).

        Lisina e arginina sao aminoacidos basicos; um peptideo rico
        neles deve ter pI > 9.
        """
        pi = calculate_isoelectric_point("KKKRRR")
        assert pi > 9.0

    def test_pi_peptideo_acido(self):
        """Peptideo rico em Asp/Glu deve ter pI baixo (acido).

        Aspartato e glutamato sao aminoacidos acidos; um peptideo rico
        neles deve ter pI < 5.
        """
        pi = calculate_isoelectric_point("DDDEEE")
        assert pi < 5.0


class TestValidatorsInstabilityIndex:
    """Testes para o indice de instabilidade de Guruprasad."""

    def test_ii_retorna_float(self):
        """Indice de instabilidade deve ser um float."""
        ii = calculate_instability_index("MKLLVV")
        assert isinstance(ii, float)

    def test_ii_sequencia_curta(self):
        """Sequencia com < 2 residuos deve retornar 0."""
        assert calculate_instability_index("M") == 0.0
        assert calculate_instability_index("") == 0.0


class TestValidatorsGRAVY:
    """Testes para o calculo de GRAVY (hidrofobicidade media)."""

    def test_gravy_retorna_float(self):
        """GRAVY deve ser um float."""
        gravy = calculate_gravy("MKLLVV")
        assert isinstance(gravy, float)

    def test_gravy_peptideo_hidrofobico(self):
        """Peptideo com residuos hidrofobicos (I, L, V) deve ter GRAVY > 0.

        Isoleucina, leucina e valina tem altos valores na escala de
        Kyte-Doolittle, indicando hidrofobicidade.
        """
        gravy = calculate_gravy("IIILLLLVVV")
        assert gravy > 0

    def test_gravy_peptideo_hidrofilico(self):
        """Peptideo com residuos carregados (D, E, K, R) deve ter GRAVY < 0.

        Aminoacidos carregados sao altamente hidrofilicos.
        """
        gravy = calculate_gravy("DDEEKKRR")
        assert gravy < 0


class TestValidatorsAromaticity:
    """Testes para o calculo de aromaticidade."""

    def test_aromaticity_retorna_float(self):
        """Aromaticidade deve ser float entre 0 e 1."""
        arom = calculate_aromaticity("FWYAAA")
        assert isinstance(arom, float)
        assert 0.0 <= arom <= 1.0

    def test_aromaticity_somente_aromaticos(self):
        """Sequencia 100% aromatica (F, W, Y) deve ter aromaticidade = 1.0."""
        arom = calculate_aromaticity("FWY")
        assert abs(arom - 1.0) < 0.001

    def test_aromaticity_sem_aromaticos(self):
        """Sequencia sem aminoacidos aromaticos deve ter aromaticidade = 0."""
        arom = calculate_aromaticity("AAALLLDDD")
        assert arom == 0.0


class TestValidatorsGCContent:
    """Testes para o calculo de conteudo GC de DNA/mRNA."""

    def test_gc_retorna_float(self):
        """GC content deve ser float entre 0 e 1."""
        gc = calculate_gc_content("ATGCATGC")
        assert isinstance(gc, float)
        assert 0 <= gc <= 1.0

    def test_gc_sequencia_50_pct(self):
        """ATGC deve ter 50% de GC."""
        gc = calculate_gc_content("ATGC")
        assert abs(gc - 0.5) < 0.01

    def test_gc_somente_gc(self):
        """Sequencia 100% GC deve retornar 1.0."""
        gc = calculate_gc_content("GGGCCC")
        assert abs(gc - 1.0) < 0.001

    def test_gc_somente_at(self):
        """Sequencia 100% AT deve retornar 0.0."""
        gc = calculate_gc_content("AAATTT")
        assert abs(gc - 0.0) < 0.001


class TestValidatorsCAI:
    """Testes para o Codon Adaptation Index."""

    def test_cai_retorna_float(self):
        """CAI deve ser float entre 0 e 1."""
        # Sequencia codificante minima (ATG = Met)
        cai = calculate_cai("ATGAAAGCC", organism="canis")
        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0

    def test_cai_organismo_invalido(self):
        """Organismo desconhecido deve levantar ValueError."""
        with pytest.raises(ValueError):
            calculate_cai("ATGAAAGCC", organism="mars_bacteria")

    def test_cai_ecoli(self):
        """CAI para sequencia otimizada de E. coli deve ser calculavel."""
        # CTG = Leu preferido em E. coli, AAA = Lys preferido, GCG = Ala preferido
        cai = calculate_cai("CTGAAAGCG", organism="ecoli")
        assert isinstance(cai, float)
        assert cai > 0.0


class TestValidatorsRestrictionSites:
    """Testes para busca de sitios de restricao."""

    def test_encontra_ecori(self):
        """Deve encontrar sitio EcoRI (GAATTC) numa sequencia que o contem."""
        sites = check_restriction_sites("AAAGAATTCAAA")
        assert "EcoRI" in sites

    def test_nao_encontra_em_sequencia_limpa(self):
        """Sequencia sem sitios de restricao deve retornar lista vazia."""
        # Sequencia simples de poli-A: sem nenhum sitio de restricao
        sites = check_restriction_sites("AAAAAAAAAAAA")
        assert len(sites) == 0
