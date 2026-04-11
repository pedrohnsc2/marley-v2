"""Testes unitarios para aso_delivery/module_f_immune_sde.

Valida os modelos ODE deterministico e SDE estocastico de macrofago
infectado por L. infantum sob tratamento com MRL-ASO-001.

Testes cobrem:
- Modelo ODE: RHS do sistema de 5 variaveis, integracao, cenarios
- Modelo SDE: Euler-Maruyama, ensemble Monte Carlo, distribuicao
- Dose-resposta: curva EC50/EC90, interpolacao
- Propriedades biologicas: clearance parasitaria, funcao dual

Todos os calculos cientificos sao executados sem mocks.
Para testes SDE, usamos semente fixa para reprodutibilidade.
"""

import math

import numpy as np
import pytest

from aso_delivery.module_f_immune_sde.ode_model import (
    IDX_A,
    IDX_I,
    IDX_N,
    IDX_P,
    IDX_T,
    N_VARS,
    MacrophageParams,
    ScenarioResult,
    SimulationConfig,
    compute_clearance_time,
    integrate_ode,
    ode_rhs,
    run_scenario,
)
from aso_delivery.module_f_immune_sde.sde_model import (
    DEFAULT_N_SIMS,
    DEFAULT_SEED,
    DEFAULT_SIGMA,
    DoseResponsePoint,
    SDEEnsembleResult,
    _interpolate_ec,
    euler_maruyama_step,
    estimate_ec_values,
    run_dose_response,
    run_sde_ensemble,
    run_single_sde,
)


# ---------------------------------------------------------------------------
# Fixtures comuns — parametros e config reduzidos para testes rapidos
# ---------------------------------------------------------------------------


@pytest.fixture
def default_params() -> MacrophageParams:
    """Parametros padrao do modelo de macrofago."""
    return MacrophageParams()


@pytest.fixture
def short_config() -> SimulationConfig:
    """Configuracao reduzida (24h, dt=0.5) para testes rapidos."""
    return SimulationConfig(t_end=24.0, dt=0.5)


@pytest.fixture
def full_config() -> SimulationConfig:
    """Configuracao completa (168h = 7 dias) para testes de clearance."""
    return SimulationConfig(t_end=168.0, dt=0.1)


# =========================================================================
# 1. Modelo ODE — lado direito (RHS)
# =========================================================================


class TestODERHS:
    """Testes do lado direito do sistema de ODEs (dy/dt)."""

    def test_five_state_variables(self):
        """Sistema deve ter 5 variaveis de estado: P, A, I, T, N."""
        assert N_VARS == 5

    def test_rhs_output_shape(self, default_params):
        """RHS deve retornar vetor de 5 derivadas."""
        y = np.array([30.0, 5.0, 1.0, 0.5, 0.3])
        dydt = ode_rhs(y, default_params)
        assert dydt.shape == (5,)

    def test_parasite_growth_without_treatment(self):
        """Sem ASO e sem NO, parasitas devem crescer logisticamente.

        dP/dt = r_P * P * (1 - P/K_P) quando k_kill=0, k_aso=0.
        Crescimento positivo quando P < K_P (capacidade de carga).
        """
        params = MacrophageParams(k_kill=0.0, k_aso=0.0, k_tlr9=0.0,
                                  k_tnf=0.0, k_uptake=0.0, A_ext=0.0)
        y = np.array([30.0, 0.0, 0.0, 0.0, 0.0])
        dydt = ode_rhs(y, params)
        # dP/dt = 0.02 * 30 * (1 - 30/200) > 0
        assert dydt[IDX_P] > 0

    def test_aso_uptake_positive(self, default_params):
        """Captacao de ASO deve ser positiva quando A_ext > 0 e A = 0.

        dA/dt = k_uptake * A_ext - k_deg * A
        Com A = 0: dA/dt = k_uptake * A_ext > 0
        """
        y = np.array([30.0, 0.0, 0.0, 0.0, 0.0])
        dydt = ode_rhs(y, default_params)
        assert dydt[IDX_A] > 0

    def test_ifn_gamma_induction_by_tlr9(self, default_params):
        """IFN-gamma deve ser induzido quando ASO ativa TLR9 (A > 0).

        dI/dt = k_tlr9 * A - d_I * I
        Com I = 0 e A > 0: dI/dt > 0
        """
        y = np.array([30.0, 10.0, 0.0, 0.0, 0.0])
        dydt = ode_rhs(y, default_params)
        assert dydt[IDX_I] > 0

    def test_no_production_from_ifn_gamma(self, default_params):
        """NO deve ser produzido quando IFN-gamma > 0.

        dN/dt = k_no * I - d_N * N
        Via IFN-gamma -> iNOS -> NO (mecanismo leishmanicida principal).
        """
        y = np.array([30.0, 10.0, 5.0, 2.0, 0.0])
        dydt = ode_rhs(y, default_params)
        assert dydt[IDX_N] > 0

    def test_parasite_killing_by_no(self):
        """NO deve reduzir a taxa de crescimento do parasita.

        Com NO alto, k_kill * N * P domina -> dP/dt < 0 (parasita morre).
        """
        params = MacrophageParams()
        # Estado com muito NO -> matanca
        y = np.array([30.0, 10.0, 5.0, 2.0, 10.0])
        dydt = ode_rhs(y, params)
        # Com NO = 10: k_kill * 10 * 30 = 30 (matanca)
        # Crescimento: 0.02 * 30 * (1 - 30/200) = 0.51
        # ASO kill: 0.05 * 10 * 30 = 15
        # Total dP/dt = 0.51 - 30 - 15 < 0
        assert dydt[IDX_P] < 0

    def test_rhs_handles_zero_state(self, default_params):
        """RHS deve funcionar com estado todo zero (sem crash)."""
        y = np.zeros(5)
        dydt = ode_rhs(y, default_params)
        assert dydt.shape == (5,)
        # Unico termo nao-zero: dA/dt = k_uptake * A_ext > 0
        assert dydt[IDX_A] > 0


# =========================================================================
# 2. Integracao ODE e cenarios deterministicos
# =========================================================================


class TestODEIntegration:
    """Testes da integracao numerica do sistema ODE."""

    def test_integration_returns_correct_shapes(self, default_params, short_config):
        """Saida deve ter shapes corretos: t=(n+1,), sol=(n+1, 5)."""
        t, sol = integrate_ode(default_params, short_config)
        n_steps = short_config.n_steps
        assert t.shape == (n_steps + 1,)
        assert sol.shape == (n_steps + 1, N_VARS)

    def test_initial_conditions_preserved(self, default_params, short_config):
        """Primeira linha da solucao deve corresponder as condicoes iniciais."""
        t, sol = integrate_ode(default_params, short_config)
        np.testing.assert_array_almost_equal(sol[0], short_config.y0)

    def test_non_negative_concentrations(self, default_params, short_config):
        """Todas as concentracoes devem ser >= 0 ao longo da simulacao.

        Concentracoes biologicas (parasitas, citocinas, NO) nao podem ser
        negativas. O integrador impoe floor de zero.
        """
        t, sol = integrate_ode(default_params, short_config)
        assert np.all(sol >= 0.0)

    def test_aso_accumulates_over_time(self, default_params, short_config):
        """ASO deve acumular no fagolisossomo ao longo do tempo.

        A(t) cresce pela captacao do meio extracelular e e degradado
        lentamente (t1/2 = 1083h >> duracao da simulacao).
        """
        t, sol = integrate_ode(default_params, short_config)
        A_initial = sol[0, IDX_A]
        A_final = sol[-1, IDX_A]
        assert A_final > A_initial

    def test_dual_function_reduces_parasites(self, default_params):
        """Com ambos os mecanismos ativos, parasitas devem ser eliminados.

        A funcao dual de MRL-ASO-001 (antisense + TLR9) deve resultar
        em clearance parasitaria completa em 7 dias de simulacao.
        """
        config = SimulationConfig(t_end=168.0, dt=0.1)
        t, sol = integrate_ode(default_params, config)
        P_final = sol[-1, IDX_P]
        P_initial = sol[0, IDX_P]
        # Reducao de pelo menos 90% em 7 dias
        reduction_pct = (1.0 - P_final / P_initial) * 100.0
        assert reduction_pct > 90.0

    def test_run_scenario_returns_result(self, default_params, short_config):
        """run_scenario deve retornar ScenarioResult com metricas."""
        result = run_scenario(
            name="test",
            description="Test scenario",
            params=default_params,
            config=short_config,
        )
        assert isinstance(result, ScenarioResult)
        assert result.initial_parasite_load == short_config.P0
        assert result.final_parasite_load >= 0.0

    def test_antisense_only_reduces_parasites(self):
        """Cenario antisense-only (sem TLR9): reducao parcial de parasitas.

        Bloqueio de splicing sozinho reduz parasitas mas nao ativa
        imunidade inata do hospedeiro.
        """
        params = MacrophageParams(k_tlr9=0.0, k_tnf=0.0)
        config = SimulationConfig(t_end=168.0, dt=0.1)
        result = run_scenario(
            name="antisense_only",
            description="Only antisense",
            params=params,
            config=config,
        )
        assert result.parasite_reduction_pct > 0.0

    def test_dual_faster_clearance_than_antisense_only(self):
        """Funcao dual deve atingir clearance MAIS RAPIDO que antisense sozinho.

        Ambos os cenarios atingem ~100% de reducao em 168h, mas a funcao
        dual (antisense + TLR9) deve chegar la mais cedo porque ativa
        a matanca via NO alem do bloqueio de splicing.
        """
        config = SimulationConfig(t_end=168.0, dt=0.1)

        # Antisense only
        params_aso = MacrophageParams(k_tlr9=0.0, k_tnf=0.0)
        result_aso = run_scenario("aso", "ASO only", params_aso, config)

        # Dual function
        params_dual = MacrophageParams()
        result_dual = run_scenario("dual", "Dual", params_dual, config)

        # Dual deve atingir 90% de clearance mais rapido
        assert result_dual.time_to_90pct_clearance_hours < result_aso.time_to_90pct_clearance_hours


# =========================================================================
# 3. Funcoes auxiliares
# =========================================================================


class TestClearanceTime:
    """Testes da funcao de calculo de tempo de clearance."""

    def test_clearance_time_for_decaying_trajectory(self):
        """Trajetoria com decaimento deve ter tempo de clearance positivo."""
        t = np.linspace(0, 100, 101)
        P = 30.0 * np.exp(-0.05 * t)  # decaimento exponencial
        ct = compute_clearance_time(t, P, threshold_fraction=0.1)
        # 90% clearance: 30*exp(-0.05*t) = 3 -> t = ln(10)/0.05 = 46.05
        assert ct > 0.0
        assert abs(ct - math.log(10) / 0.05) < 2.0  # tolerancia de 2h

    def test_clearance_time_not_achieved(self):
        """Trajetoria crescente nunca atinge clearance -> retorna -1."""
        t = np.linspace(0, 100, 101)
        P = 30.0 + 0.5 * t  # crescendo
        ct = compute_clearance_time(t, P, threshold_fraction=0.1)
        assert ct == -1.0

    def test_clearance_time_zero_initial(self):
        """Se P0 = 0, clearance ja foi atingida -> retorna 0."""
        t = np.linspace(0, 10, 11)
        P = np.zeros(11)
        ct = compute_clearance_time(t, P)
        assert ct == 0.0


# =========================================================================
# 4. Modelo SDE — Euler-Maruyama
# =========================================================================


class TestSDE:
    """Testes do modelo estocastico (SDE com ruido Wiener multiplicativo)."""

    def test_euler_maruyama_step_shape(self, default_params):
        """Um passo Euler-Maruyama deve retornar vetor de 5 componentes."""
        y = np.array([30.0, 5.0, 1.0, 0.5, 0.3])
        dW = np.array([0.1, -0.05, 0.02, 0.03, -0.01])
        y_new = euler_maruyama_step(y, default_params, dt=0.1, sigma=0.1, dW=dW)
        assert y_new.shape == (5,)

    def test_euler_maruyama_non_negative(self, default_params):
        """Passo Euler-Maruyama deve impor nao-negatividade.

        Ruido multiplicativo pode gerar valores negativos — o integrador
        impoe floor de zero para concentracoes biologicas.
        """
        y = np.array([0.1, 0.01, 0.001, 0.001, 0.001])
        # dW fortemente negativo para tentar forcar valores negativos
        dW = np.array([-5.0, -5.0, -5.0, -5.0, -5.0])
        y_new = euler_maruyama_step(y, default_params, dt=0.1, sigma=0.5, dW=dW)
        assert np.all(y_new >= 0.0)

    def test_zero_sigma_equals_ode(self, default_params, short_config):
        """Com sigma=0 (sem ruido), SDE deve reproduzir a solucao ODE.

        Euler-Maruyama com sigma=0 equivale a Euler explicito.
        """
        rng = np.random.default_rng(42)
        sol_sde = run_single_sde(default_params, short_config, sigma=0.0, rng=rng)
        t, sol_ode = integrate_ode(default_params, short_config)
        # Comparar trajetorias ponto a ponto (devem ser identicas)
        np.testing.assert_array_almost_equal(sol_sde, sol_ode, decimal=6)

    def test_single_sde_trajectory_shape(self, default_params, short_config):
        """Trajetoria SDE deve ter shape (n_steps+1, 5)."""
        rng = np.random.default_rng(42)
        sol = run_single_sde(default_params, short_config, sigma=0.1, rng=rng)
        expected_rows = short_config.n_steps + 1
        assert sol.shape == (expected_rows, N_VARS)


# =========================================================================
# 5. Ensemble Monte Carlo
# =========================================================================


class TestMonteCarloEnsemble:
    """Testes do ensemble Monte Carlo para quantificar variabilidade."""

    def test_ensemble_result_type(self, default_params, short_config):
        """Ensemble deve retornar SDEEnsembleResult."""
        result = run_sde_ensemble(
            name="test",
            params=default_params,
            config=short_config,
            n_sims=10,  # poucos para teste rapido
            sigma=0.1,
            seed=42,
        )
        assert isinstance(result, SDEEnsembleResult)

    def test_ensemble_n_simulations(self, default_params, short_config):
        """Numero de simulacoes registrado deve corresponder ao solicitado."""
        result = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=50, sigma=0.1, seed=42,
        )
        assert result.n_simulations == 50

    def test_ensemble_mean_close_to_ode(self, default_params, short_config):
        """Media do ensemble SDE deve estar proxima da solucao ODE.

        Com muitas simulacoes, a media do ensemble converge para a
        solucao deterministica (lei dos grandes numeros).
        Verificamos usando concentracao de ASO (que nao decai para zero)
        para evitar problemas de divisao com valores muito pequenos.
        """
        result = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=200, sigma=0.1, seed=42,
        )
        t, sol_ode = integrate_ode(default_params, short_config)

        # Concentracao de ASO (variavel estavel, nao decai para zero)
        A_ode_final = sol_ode[-1, IDX_A]
        A_sde_mean_final = result.mean[-1, IDX_A]
        # Tolerancia relativa de 15% para concentracao de ASO
        assert abs(A_sde_mean_final - A_ode_final) / max(A_ode_final, 1e-6) < 0.15

    def test_ensemble_std_positive(self, default_params, short_config):
        """Desvio padrao do ensemble deve ser > 0 (ruido gera dispersao).

        Se std = 0, o ruido nao esta sendo aplicado corretamente.
        """
        result = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=50, sigma=0.1, seed=42,
        )
        # STD da carga parasitaria no ultimo passo
        assert result.std[-1, IDX_P] > 0.0

    def test_ensemble_ci_contains_mean(self, default_params, short_config):
        """Intervalo de confianca 95% deve conter a media.

        percentile_2.5 <= mean <= percentile_97.5 para cada variavel.
        """
        result = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=100, sigma=0.1, seed=42,
        )
        # Verificar no ultimo passo temporal
        for var_idx in range(N_VARS):
            mean_val = result.mean[-1, var_idx]
            low = result.percentile_2_5[-1, var_idx]
            high = result.percentile_97_5[-1, var_idx]
            assert low <= mean_val <= high

    def test_ensemble_clearance_probability_between_0_and_1(
        self, default_params, short_config
    ):
        """Probabilidade de clearance deve estar entre 0 e 1."""
        result = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=50, sigma=0.1, seed=42,
        )
        assert 0.0 <= result.clearance_probability <= 1.0

    def test_ensemble_reproducibility(self, default_params, short_config):
        """Mesma semente deve produzir resultados identicos.

        Reprodutibilidade e essencial para validacao cientifica.
        """
        result1 = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=20, sigma=0.1, seed=42,
        )
        result2 = run_sde_ensemble(
            "test", default_params, short_config,
            n_sims=20, sigma=0.1, seed=42,
        )
        np.testing.assert_array_equal(result1.mean, result2.mean)
        assert result1.clearance_probability == result2.clearance_probability

    def test_higher_sigma_more_variance(self, default_params, short_config):
        """Maior sigma (ruido) deve produzir maior variancia no ensemble.

        sigma=0.2 deve ter dispersao maior que sigma=0.05.
        """
        result_low = run_sde_ensemble(
            "low", default_params, short_config,
            n_sims=100, sigma=0.05, seed=42,
        )
        result_high = run_sde_ensemble(
            "high", default_params, short_config,
            n_sims=100, sigma=0.2, seed=42,
        )
        # Desvio padrao da carga parasitaria final
        assert result_high.std_final_parasite_load > result_low.std_final_parasite_load


# =========================================================================
# 6. Dose-resposta e EC50/EC90
# =========================================================================


class TestDoseResponse:
    """Testes da curva dose-resposta e estimativa de EC50/EC90."""

    def test_higher_dose_more_killing(self, default_params):
        """Dose mais alta deve produzir maior reducao parasitaria.

        Relacao dose-resposta monotonicamente crescente para o ASO.
        """
        config = SimulationConfig(t_end=72.0, dt=0.5)
        doses = [1.0, 50.0]
        results = run_dose_response(
            doses_uM=doses,
            base_params=default_params,
            config=config,
            n_sims=20,  # poucos para teste rapido
            sigma=0.1,
            seed=42,
        )
        assert results[1].mean_parasite_reduction_pct > results[0].mean_parasite_reduction_pct

    def test_interpolate_ec_basic(self):
        """Interpolacao linear deve encontrar valor correto.

        Para reductions [20, 60] em doses [1, 3]:
        EC50 deve ser interpolado em dose = 1 + (50-20)/(60-20) * (3-1) = 2.5
        """
        doses = [1.0, 3.0]
        reductions = [20.0, 60.0]
        ec50 = _interpolate_ec(doses, reductions, 50.0)
        assert abs(ec50 - 2.5) < 0.01

    def test_interpolate_ec_not_found(self):
        """Se nenhuma dose atinge o alvo, deve retornar -1."""
        doses = [1.0, 5.0]
        reductions = [10.0, 30.0]
        ec90 = _interpolate_ec(doses, reductions, 90.0)
        assert ec90 == -1.0

    def test_interpolate_ec_first_exceeds(self):
        """Se a primeira dose ja supera o alvo, retornar a primeira dose."""
        doses = [0.5, 1.0, 5.0]
        reductions = [60.0, 80.0, 95.0]
        ec50 = _interpolate_ec(doses, reductions, 50.0)
        assert ec50 == 0.5

    def test_estimate_ec_values_returns_dict(self):
        """estimate_ec_values deve retornar dict com EC50_uM e EC90_uM."""
        # Criar dados ficticios de dose-resposta
        points = [
            DoseResponsePoint(1.0, 0.1, 100.0, 20.0, 30.0),
            DoseResponsePoint(5.0, 0.3, 60.0, 10.0, 60.0),
            DoseResponsePoint(10.0, 0.6, 40.0, 5.0, 80.0),
            DoseResponsePoint(50.0, 0.9, 20.0, 1.0, 95.0),
        ]
        ec = estimate_ec_values(points)
        assert "EC50_uM" in ec
        assert "EC90_uM" in ec
