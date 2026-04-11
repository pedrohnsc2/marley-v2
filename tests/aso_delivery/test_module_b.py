"""Testes unitarios para aso_delivery/module_b_membrane.

Valida os modelos de interacao ASO-membrana e propriedades fisico-quimicas:
1. Peso molecular e propriedades de ASOs com modificacoes PS e LNA
2. Carga liquida dependente do pH (Henderson-Hasselbalch para PS)
3. Interacao com bicamada lipidica (eletrostatica, hidrofobica, Born)
4. Coeficiente de particao membrana/agua

A conclusao biologica: ASOs sao macromoleculas altamente carregadas
que NAO podem cruzar membranas por difusao passiva. Endocitose ativa
e o UNICO mecanismo viavel — e macrofagos sao os melhores captadores.
"""

import math

import pytest

from aso_delivery.module_b_membrane.models import (
    MW_NUCLEOTIDE_PO,
    MW_PS_ADDITION_PER_LINK,
    MW_LNA_ADDITION_PER_RESIDUE,
    LOGP_PER_LINK_PO,
    LOGP_PER_LINK_PS,
    PKA_PS_LINKAGE,
    KB_KCAL,
    T_PHYSIOLOGICAL_K,
    EPSILON_WATER,
    EPSILON_MEMBRANE,
    compute_molecular_weight,
    compute_logp,
    compute_logd,
    compute_net_charge,
    compute_hydrodynamic_radius,
    compute_lipinski_violations,
    compute_physicochemical_properties,
    compute_electrostatic_repulsion,
    compute_hydrophobic_insertion,
    compute_born_solvation_energy,
    compute_membrane_interaction,
    compute_partition_coefficient,
)


# ---------------------------------------------------------------------------
# 1. Peso molecular
# ---------------------------------------------------------------------------


class TestMolecularWeight:
    """Testes para compute_molecular_weight."""

    def test_mrl_aso_001_mw(self):
        """MW do MRL-ASO-001 (25 nt, 24 PS, 10 LNA) deve ser ~8774 Da.
        Calculo: 25*330 + 24*16 + 10*14 = 8250 + 384 + 140 = 8774.
        """
        mw = compute_molecular_weight(
            length_nt=25, n_ps_linkages=24, n_lna_residues=10,
        )
        assert mw == pytest.approx(8774.0, abs=0.1)

    def test_po_only_no_modifications(self):
        """ASO PO sem modificacoes: MW = n * 330 Da apenas.
        Backbone nativo mais simples.
        """
        mw = compute_molecular_weight(length_nt=25, n_ps_linkages=0, n_lna_residues=0)
        assert mw == pytest.approx(25 * MW_NUCLEOTIDE_PO, abs=0.1)

    def test_ps_adds_mass(self):
        """Cada ligacao PS adiciona 16 Da (S-O substituicao, net +16).
        PS deve ser mais pesado que PO para mesmo comprimento.
        """
        mw_po = compute_molecular_weight(25, 0, 0)
        mw_ps = compute_molecular_weight(25, 24, 0)
        assert mw_ps > mw_po
        assert mw_ps - mw_po == pytest.approx(24 * MW_PS_ADDITION_PER_LINK, abs=0.1)

    def test_lna_adds_mass(self):
        """Cada residuo LNA adiciona 14 Da (ponte metileno CH2).
        Gapmer com LNA e mais pesado que PS puro.
        """
        mw_ps = compute_molecular_weight(25, 24, 0)
        mw_gapmer = compute_molecular_weight(25, 24, 10)
        assert mw_gapmer > mw_ps
        assert mw_gapmer - mw_ps == pytest.approx(
            10 * MW_LNA_ADDITION_PER_RESIDUE, abs=0.1,
        )

    def test_zero_length_zero_weight(self):
        """Comprimento zero deve dar MW zero (caso limite)."""
        mw = compute_molecular_weight(0, 0, 0)
        assert mw == 0.0


# ---------------------------------------------------------------------------
# 2. logP e logD
# ---------------------------------------------------------------------------


class TestLogPLogD:
    """Testes para compute_logp, compute_logd e compute_net_charge."""

    def test_logp_ps_less_negative_than_po(self):
        """PS e ligeiramente menos polar que PO (enxofre vs oxigenio).
        logP menos negativo = mais lipofilico. Biologicamente, isso
        permite interacao com proteinas plasmaticas e scavenger receptors.
        """
        logp_po = compute_logp(24, 0, backbone_type="PO")
        logp_ps = compute_logp(24, 0, backbone_type="PS")
        assert logp_ps > logp_po  # menos negativo = mais lipofilico

    def test_logp_extremely_negative(self):
        """Para ASOs, logP deve ser extremamente negativo (< -10).
        ASOs sao macromoleculas polianioncas — muito mais polares
        que qualquer farmaco convencional.
        """
        logp = compute_logp(24, 10, backbone_type="PS")
        assert logp < -10

    def test_logd_more_negative_than_logp(self):
        """logD e ainda mais negativo que logP quando a molecula esta
        carregada. As 24 cargas negativas do PS-ASO fazem logD
        disparar para valores extremos (< -20).
        """
        logp = compute_logp(24, 10, backbone_type="PS")
        logd = compute_logd(logp, -24.0)
        assert logd < logp

    def test_net_charge_physiological_ph(self):
        """A pH 7.4, todas as ligacoes PS estao desprotonadas (pKa ~1.5).
        Carga liquida = -(numero de ligacoes PS) = -24 para 25-mer.
        """
        charge = compute_net_charge(7.4, 24)
        assert charge == pytest.approx(-24.0, abs=0.1)

    def test_net_charge_phagolysosome_ph(self):
        """Mesmo no fagolisossomo (pH 4.5), a carga permanece ~-24.
        pH 4.5 >> pKa 1.5 -> quase todas as ligacoes desprotonadas.
        A carga NÃO muda significativamente ao longo do pathway endocitico.
        """
        charge = compute_net_charge(4.5, 24)
        assert charge == pytest.approx(-24.0, abs=0.1)

    def test_net_charge_at_pka(self):
        """No pH = pKa (1.5), metade das ligacoes esta desprotonada.
        Carga = -24 * 0.5 = -12 (cenario extremo nao fisiologico).
        """
        charge = compute_net_charge(PKA_PS_LINKAGE, 24)
        assert charge == pytest.approx(-12.0, abs=0.5)

    def test_net_charge_zero_linkages(self):
        """Sem ligacoes PS, carga e zero (caso limite)."""
        charge = compute_net_charge(7.4, 0)
        assert charge == 0.0


# ---------------------------------------------------------------------------
# 3. Raio hidrodinamico e violacoes de Lipinski
# ---------------------------------------------------------------------------


class TestHydrodynamicAndLipinski:
    """Testes para compute_hydrodynamic_radius e compute_lipinski_violations."""

    def test_rh_increases_with_length(self):
        """Moleculas maiores tem raio hidrodinamico maior.
        Modelo de cadeia semiflexivel: R_h ~ sqrt(comprimento).
        """
        rh_20 = compute_hydrodynamic_radius(20)
        rh_25 = compute_hydrodynamic_radius(25)
        rh_30 = compute_hydrodynamic_radius(30)
        assert rh_20 < rh_25 < rh_30

    def test_rh_mrl_aso_001(self):
        """Raio do MRL-ASO-001 (25 nt) deve ser ~2.72 nm.
        Calculado pelo modelo worm-like chain.
        """
        rh = compute_hydrodynamic_radius(25)
        assert rh == pytest.approx(2.72, abs=0.05)

    def test_lipinski_all_violations(self):
        """ASO com MW 8774, logP -18, carga -24 viola TODAS as regras.
        ASOs nao sao moleculas pequenas — isso justifica a necessidade
        de mecanismos de uptake ativos (endocitose).
        """
        violations = compute_lipinski_violations(
            mw=8774.0, logp=-18.0, net_charge=-24.0,
        )
        assert violations == 4  # MW>500, |logP|>5, |charge|>3, MW>1000

    def test_lipinski_small_molecule_no_violations(self):
        """Molecula pequena ideal (aspirina-like) nao viola nenhuma regra."""
        violations = compute_lipinski_violations(
            mw=180.0, logp=1.2, net_charge=0.0,
        )
        assert violations == 0


# ---------------------------------------------------------------------------
# 4. Interacao com bicamada lipidica
# ---------------------------------------------------------------------------


class TestMembraneInteraction:
    """Testes para modelos de interacao ASO-membrana."""

    def test_electrostatic_repulsion_positive(self):
        """A repulsao eletrostatica entre ASO(-24) e membrana(-30 mV)
        deve ser POSITIVA (desfavoravel). Cargas de mesmo sinal se repelem.
        """
        e = compute_electrostatic_repulsion(-24.0)
        assert e > 0

    def test_electrostatic_zero_charge(self):
        """Molecula neutra nao sofre repulsao eletrostatica.
        So cargas interagem com o potencial de superficie.
        """
        e = compute_electrostatic_repulsion(0.0)
        assert e == pytest.approx(0.0, abs=1e-6)

    def test_hydrophobic_ps_favorable(self):
        """PS tem contribuicao hidrofobica levemente FAVORAVEL (negativa).
        O enxofre e parcialmente apolar, favorecendo insercao na membrana.
        """
        e = compute_hydrophobic_insertion(24, backbone_type="PS")
        assert e < 0

    def test_hydrophobic_po_unfavorable(self):
        """PO tem contribuicao hidrofobica DESFAVORAVEL (positiva).
        O fosfato diester e muito polar — penalidade para entrar no lipideo.
        """
        e = compute_hydrophobic_insertion(24, backbone_type="PO")
        assert e > 0

    def test_born_solvation_always_positive(self):
        """Energia de Born e sempre POSITIVA (custo de desolvatacao).
        Transferir cargas de agua (eps=74) para lipideo (eps=2) e
        termodinamicamente proibitivo.
        """
        e = compute_born_solvation_energy(-24.0, 2.72)
        assert e > 0

    def test_born_solvation_increases_with_charge(self):
        """Mais carga = maior custo de Born (proporcional a z^2).
        O termo quadratico faz a barreira crescer muito rapidamente.
        """
        e_low = compute_born_solvation_energy(-5.0, 2.72)
        e_high = compute_born_solvation_energy(-24.0, 2.72)
        assert e_high > e_low

    def test_total_barrier_positive(self):
        """A barreira total de membrana deve ser fortemente positiva.
        Born + eletrostatica >> contribuicao hidrofobica do PS.
        """
        interaction = compute_membrane_interaction(
            name="test", net_charge=-24.0, n_linkages=24,
            radius_nm=2.72, backbone_type="PS",
        )
        assert interaction.total_barrier_kcal > 0

    def test_passive_diffusion_infeasible_for_aso(self):
        """Difusao passiva e impossivel para um ASO (barreira >> 20 kT).
        Esta e a conclusao principal do Modulo B — uptake precisa
        ser ativo (endocitose).
        """
        interaction = compute_membrane_interaction(
            name="MRL-ASO-001", net_charge=-24.0, n_linkages=24,
            radius_nm=2.72, backbone_type="PS",
        )
        assert interaction.passive_diffusion_feasible is False
        assert interaction.barrier_in_kt > 20.0


# ---------------------------------------------------------------------------
# 5. Coeficiente de particao membrana/agua
# ---------------------------------------------------------------------------


class TestPartitionCoefficient:
    """Testes para compute_partition_coefficient."""

    def test_log_k_extremely_negative(self):
        """logK deve ser extremamente negativo para ASOs.
        Para barreira de ~50 kcal/mol: logK ~ -35. Nenhuma molecula
        de ASO particionaria na membrana no equilibrio.
        """
        partition = compute_partition_coefficient(
            name="test", backbone_type="PS", barrier_kcal=50.0,
        )
        assert partition.log_k_partition < -30

    def test_k_value_underflow_protection(self):
        """Para barreiras muito altas, K numerico deve ser 0 (sem underflow).
        exp(-700) e o limite pratico de float64.
        """
        partition = compute_partition_coefficient(
            name="test", backbone_type="PS", barrier_kcal=500.0,
        )
        assert partition.k_partition == 0.0

    def test_permeability_zero_when_k_zero(self):
        """Quando K = 0, permeabilidade tambem e zero.
        Nenhum transporte passivo e possivel.
        """
        partition = compute_partition_coefficient(
            name="test", backbone_type="PS", barrier_kcal=500.0,
        )
        assert partition.passive_permeability_cm_s == 0.0

    def test_lower_barrier_higher_k(self):
        """Barreira menor = maior coeficiente de particao.
        Relacao exponencial: K = exp(-dG/RT).
        """
        p_low = compute_partition_coefficient("a", "PS", 10.0)
        p_high = compute_partition_coefficient("b", "PS", 50.0)
        assert p_low.log_k_partition > p_high.log_k_partition


# ---------------------------------------------------------------------------
# 6. Perfil completo (integracao)
# ---------------------------------------------------------------------------


class TestPhysicochemicalProperties:
    """Testes de integracao para compute_physicochemical_properties."""

    def test_mrl_aso_001_full_profile(self):
        """Perfil completo do MRL-ASO-001 deve ter todos os campos preenchidos
        e valores fisicamente razoaveis.
        """
        props = compute_physicochemical_properties(
            name="MRL-ASO-001",
            length_nt=25,
            n_ps_linkages=24,
            n_lna_residues=10,
            backbone_type="PS",
        )
        assert props.molecular_weight_da > 8000
        assert props.logp < -10
        assert props.net_charge_ph74 < -20
        assert props.net_charge_ph45 < -20
        assert props.hydrodynamic_radius_nm > 2.0
        assert props.lipinski_violations == 4

    def test_charge_similar_at_both_phs(self):
        """Carga deve ser similar em pH 7.4 e 4.5 (pKa PS = 1.5).
        O ASO permanece fortemente anionico ao longo de todo o
        pathway endocitico.
        """
        props = compute_physicochemical_properties(
            name="test", length_nt=25, n_ps_linkages=24,
            n_lna_residues=10, backbone_type="PS",
        )
        diff = abs(props.net_charge_ph74 - props.net_charge_ph45)
        assert diff < 0.1
