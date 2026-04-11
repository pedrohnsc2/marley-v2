"""Testes para a Plataforma B -- Vacina recombinante em E. coli.

Testa os sub-modulos B1 (construct.py) e B2 (codon_optimizer.py) usando
dados sinteticos ou derivados da SSOT de epitopos.
"""

import pytest

from vaccine_platforms.platform_b_ecoli.construct import (
    EcoliConstruct,
    HIS6_TAG,
    TEV_SITE,
    build_ecoli_construct,
    calculate_instability_index,
    calculate_isoelectric_point,
    calculate_molecular_weight,
    extract_unique_epitopes,
    validate_construct,
)
from vaccine_platforms.platform_b_ecoli.codon_optimizer import (
    add_flanking_elements,
    calculate_cai,
    calculate_gc_content,
    find_homopolymer_runs,
    find_rare_codons,
    find_restriction_sites,
    find_shine_dalgarno,
    get_best_codon,
    optimize_codons,
    run_codon_optimization,
)
from vaccine_platforms.shared.epitopes import EPITOPES


# ---------------------------------------------------------------------------
# Dados de teste: epitopos no formato do construct_card.json
# ---------------------------------------------------------------------------

# Simula o formato do construct_card.json (lista de dicts com 'peptide')
_MOCK_EPITOPES_RAW = [
    {"peptide": "RMMRSLTPF", "gene_id": "LINF_240013900"},
    {"peptide": "YIYETFASM", "gene_id": "LINF_240013900"},
    {"peptide": "SLMCVFYFK", "gene_id": "LINF_230010300"},
    {"peptide": "YLAALVPAL", "gene_id": "LINF_230010300"},
    {"peptide": "LIIEDLSLV", "gene_id": "LINF_160022200"},
    {"peptide": "FAFSVSARR", "gene_id": "LINF_230010300"},
    {"peptide": "MILGTFVRL", "gene_id": "LINF_240013900"},
    {"peptide": "MQNVTFVPK", "gene_id": "LINF_160022200"},
    {"peptide": "RILESISNV", "gene_id": "LINF_100010400"},
    {"peptide": "ILYNKISGL", "gene_id": "LINF_160022200"},
    {"peptide": "LLTANVCYK", "gene_id": "LINF_080014800"},
    # Duplicata intencional -- CPB tem 5 copias no genoma
    {"peptide": "LLTANVCYK", "gene_id": "LINF_080014900"},
    {"peptide": "LLTANVCYK", "gene_id": "LINF_080015000"},
]


# ===========================================================================
# B1: Construcao do constructo proteico
# ===========================================================================


class TestExtractUniqueEpitopes:
    """Testa a extracao de epitopos unicos do formato raw."""

    def test_remove_duplicatas(self):
        """Deve extrair exatamente 11 epitopos unicos de 13 entradas.

        LLTANVCYK aparece 3x nos dados mock (representando as 5 copias
        do gene CPB no genoma de L. infantum). Deve manter apenas 1.
        """
        unique = extract_unique_epitopes(_MOCK_EPITOPES_RAW)
        assert len(unique) == 11

    def test_preserva_ordem(self):
        """Epitopos devem manter a ordem de aparicao original.

        A ordem foi otimizada pelo modulo 04 para minimizar energia livre
        de juncao entre epitopos adjacentes.
        """
        unique = extract_unique_epitopes(_MOCK_EPITOPES_RAW)
        assert unique[0] == "RMMRSLTPF"
        assert unique[-1] == "LLTANVCYK"


class TestBuildEcoliConstruct:
    """Testa a montagem do constructo para expressao em E. coli."""

    @pytest.fixture
    def construct(self) -> EcoliConstruct:
        """Constroi o constructo uma unica vez para os testes da classe."""
        return build_ecoli_construct(_MOCK_EPITOPES_RAW)

    def test_inicia_com_met(self, construct: EcoliConstruct):
        """Proteina deve iniciar com Met (inicio universal da traducao).

        O Met inicial e codificado pelo ATG do sitio NdeI no vetor pET-28a.
        """
        assert construct.protein_sequence.startswith("M")

    def test_contem_his6_tag(self, construct: EcoliConstruct):
        """Constructo deve conter His6-tag para purificacao por IMAC (Ni-NTA).

        A His6-tag coordena ions niquel na resina de afinidade, permitindo
        purificacao em um unico passo com pureza >90%.
        """
        assert HIS6_TAG in construct.protein_sequence

    def test_contem_sitio_tev(self, construct: EcoliConstruct):
        """Constructo deve conter sitio de clivagem TEV (ENLYFQS).

        A protease TEV permite remover a His6-tag apos purificacao,
        restaurando a proteina quasi-nativa.
        """
        assert TEV_SITE in construct.protein_sequence

    def test_11_epitopos_unicos(self, construct: EcoliConstruct):
        """Constructo deve conter exatamente 11 epitopos unicos."""
        assert len(construct.unique_epitopes) == 11

    def test_peso_molecular_positivo(self, construct: EcoliConstruct):
        """Peso molecular deve ser positivo e razoavel (20-60 kDa).

        Proteinas E. coli expressas em BL21 tipicamente ficam abaixo de
        100 kDa para boa expressao soluvel.
        """
        assert construct.molecular_weight_da > 20_000
        assert construct.molecular_weight_da < 60_000

    def test_pi_entre_0_e_14(self, construct: EcoliConstruct):
        """Ponto isoeletrico deve estar entre 0 e 14."""
        assert 0 < construct.isoelectric_point < 14

    def test_comprimento_coerente(self, construct: EcoliConstruct):
        """length_aa deve corresponder ao comprimento real da sequencia."""
        assert construct.length_aa == len(construct.protein_sequence)


class TestValidateConstruct:
    """Testa as validacoes de integridade do constructo E. coli."""

    def test_constructo_valido_sem_erros(self):
        """Constructo padrao nao deve gerar erros (pode ter avisos).

        Os erros fatais sao prefixados com 'ERRO:'; avisos com 'AVISO:'.
        """
        construct = build_ecoli_construct(_MOCK_EPITOPES_RAW)
        warnings = validate_construct(construct)
        erros = [w for w in warnings if w.startswith("ERRO:")]
        assert len(erros) == 0, f"Erros encontrados: {erros}"


# ===========================================================================
# B2: Otimizacao de codons para E. coli
# ===========================================================================


class TestCodonOptimizer:
    """Testa o otimizador de codons para E. coli K12."""

    def test_get_best_codon_leucina(self):
        """Codon preferido para Leu em E. coli deve ser CTG (~50%).

        Diferente de eucariotos onde CTC domina, em E. coli o CTG e o codon
        de leucina mais abundante (reflexo da composicao genomica).
        """
        assert get_best_codon("L") == "CTG"

    def test_get_best_codon_met(self):
        """Metionina tem codon unico: ATG."""
        assert get_best_codon("M") == "ATG"

    def test_get_best_codon_invalido(self):
        """Aminoacido invalido deve levantar ValueError."""
        with pytest.raises(ValueError):
            get_best_codon("X")

    def test_optimize_codons_comprimento_correto(self):
        """DNA otimizado deve ter exatamente 3x o comprimento da proteina.

        Cada aminoacido e codificado por um triplet (3 nucleotideos).
        """
        protein = "MKLLVVRR"
        dna = optimize_codons(protein)
        assert len(dna) == len(protein) * 3

    def test_gc_content_range(self):
        """GC content do DNA otimizado para E. coli deve estar entre 40-60%.

        O genoma de E. coli K12 tem ~50.8% GC. Sequencias otimizadas devem
        ficar perto deste valor para estabilidade do mRNA e eficiencia de traducao.
        """
        protein = "MKLLVVRRAADEEFGHIKLMNPQRSTVWY"
        dna = optimize_codons(protein)
        gc = calculate_gc_content(dna)
        assert 0.35 <= gc <= 0.65, f"GC content fora do range: {gc:.1%}"

    def test_cai_alto(self):
        """CAI do DNA otimizado deve ser > 0.70 (idealmente > 0.80).

        O Codon Adaptation Index mede a adaptacao ao hospedeiro. Como usamos
        estrategia de max-frequency, o CAI deve ser naturalmente alto.
        """
        protein = "MKLLVVRRAADEEFGHIKLMNPQRSTVWY"
        dna = optimize_codons(protein)
        cai = calculate_cai(dna)
        assert cai > 0.70, f"CAI muito baixo: {cai:.4f}"

    def test_add_flanking_elements(self):
        """Inserto completo deve ter NdeI no inicio e XhoI no final.

        Estrutura: NdeI(CATATG) + CDS + TAA(stop) + XhoI(CTCGAG).
        O ATG do NdeI serve como codon de inicio.
        """
        coding = "ATGAAAGCC"  # M-K-A
        insert = add_flanking_elements(coding)
        assert insert.startswith("CATATG")
        assert insert.endswith("CTCGAG")
        assert "TAA" in insert  # Stop codon antes do XhoI

    def test_run_codon_optimization_completo(self):
        """Pipeline completo de otimizacao deve retornar resultado valido.

        Testa a sequencia completa: otimizacao + flanqueamento + metricas.
        """
        protein = "MKLLVVRRAADEEFGHIKLMNPQRSTVWY"
        result = run_codon_optimization(protein)

        # Metricas basicas
        assert result.length_nt > 0
        assert 0 < result.gc_content < 1.0
        assert 0 < result.cai <= 1.0
        assert isinstance(result.rare_codons, list)
        assert isinstance(result.warnings, list)

        # DNA deve comecar com CATATG (NdeI) e terminar com CTCGAG (XhoI)
        assert result.dna_sequence.startswith("CATATG")
        assert result.dna_sequence.endswith("CTCGAG")


class TestCodonQualityChecks:
    """Testa as verificacoes de qualidade pos-otimizacao."""

    def test_find_homopolymer_runs_detecta(self):
        """Deve detectar corridas homopolimericas > 8 nt.

        Homopolimeros longos causam erros de slippage na DNA polimerase
        durante sintese genica.
        """
        # 10 As consecutivos -- deve ser detectado
        runs = find_homopolymer_runs("GGGAAAAAAAAAAAGGG")
        assert len(runs) >= 1

    def test_find_homopolymer_runs_ignora_curtos(self):
        """Nao deve reportar corridas de 6 nt ou menos."""
        runs = find_homopolymer_runs("GGGAAAAAAGGG")
        assert len(runs) == 0

    def test_find_shine_dalgarno(self):
        """Deve detectar sequencias Shine-Dalgarno internas (AAGGAG).

        SD internas podem causar iniciacao espuria de traducao em E. coli,
        produzindo fragmentos truncados.
        """
        positions = find_shine_dalgarno("ATGAAGGAGATG")
        assert len(positions) >= 1

    def test_find_restriction_sites_ndei_interno(self):
        """Deve detectar sitio NdeI interno (CATATG) na regiao codificante.

        Sitios de restricao internos impossibilitam a clonagem direcional.
        """
        # Inserto com NdeI flanqueadores e um NdeI interno
        insert = "CATATG" + "AAACATATGAAA" + "TAA" + "CTCGAG"
        sites = find_restriction_sites(insert)
        ndei_sites = [s for s in sites if s["enzyme"] == "NdeI"]
        assert len(ndei_sites) >= 1
