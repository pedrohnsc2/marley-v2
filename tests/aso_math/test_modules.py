"""Testes de resposta conhecida para os 5 modulos do aso_math.

Cada teste executa um modulo e verifica que os resultados-chave
correspondem aos valores esperados para L. infantum (MRL-ASO-001).

Estes testes executam os modulos reais — sao mais lentos que unit tests
mas verificam que o pipeline completo funciona.
"""

import importlib

import pytest


# ---------------------------------------------------------------------------
# Helpers — importar modulos com nomes iniciando por digito via importlib
# ---------------------------------------------------------------------------

def _load_module(module_path: str):
    """Importa modulo pelo caminho completo via importlib."""
    return importlib.import_module(module_path)


# ---------------------------------------------------------------------------
# Cache de resultados — cada modulo e executado apenas uma vez por sessao
# de teste para evitar re-execucao desnecessaria.
# ---------------------------------------------------------------------------

_CACHE: dict[str, dict] = {}


def _run_module(module_path: str, **kwargs) -> dict:
    """Executa main() do modulo e faz cache do resultado."""
    if module_path not in _CACHE:
        mod = _load_module(module_path)
        _CACHE[module_path] = mod.main(**kwargs)
    return _CACHE[module_path]


# =========================================================================
# Modulo 01: Paisagem termodinamica
# =========================================================================


class TestModule01Thermodynamic:
    """Modulo 01: Paisagem termodinamica."""

    MODULE = "aso_math.01_thermodynamic_landscape.run"

    def _result(self) -> dict:
        return _run_module(self.MODULE)

    def test_status_success(self):
        """Modulo deve completar com status success."""
        result = self._result()
        assert result["status"] == "success"

    def test_wildtype_dg(self):
        """dG do wildtype MRL-ASO-001 deve ser -27.97."""
        result = self._result()
        wt = result["data"]["wildtype"]
        assert wt["dg_binding_kcal"] == -27.97

    def test_wildtype_tm(self):
        """Tm do wildtype deve ser 68.48."""
        result = self._result()
        wt = result["data"]["wildtype"]
        assert wt["tm_celsius"] == 68.48

    def test_75_point_mutants(self):
        """Deve gerar exatamente 75 mutantes pontuais (25 pos x 3 bases)."""
        result = self._result()
        assert result["data"]["point_mutants"]["total"] == 75

    def test_2700_double_mutants(self):
        """Deve gerar exatamente 2700 mutantes duplos (C(25,2) x 9)."""
        result = self._result()
        assert result["data"]["double_mutants"]["total"] == 2700

    def test_175_length_variants(self):
        """Deve avaliar 175 variantes de comprimento."""
        result = self._result()
        assert result["data"]["length_variants"]["total_evaluated"] == 175

    def test_heatmap_dimensions(self):
        """Heatmap deve ter 25 linhas (posicoes) x 4 colunas (bases)."""
        result = self._result()
        heatmap = result["data"]["heatmap"]
        assert len(heatmap["dg_matrix"]) == 25
        assert len(heatmap["dg_matrix"][0]) == 4
        assert heatmap["bases"] == ["A", "C", "G", "T"]


# =========================================================================
# Modulo 02: Prova de seletividade
# =========================================================================


class TestModule02Selectivity:
    """Modulo 02: Prova de seletividade."""

    MODULE = "aso_math.02_selectivity_proof.run"

    def _result(self) -> dict:
        return _run_module(self.MODULE)

    def test_passes(self):
        """Modulo deve reportar status pass."""
        result = self._result()
        assert result["status"] == "pass"

    def test_max_complementarity_below_threshold(self):
        """Max complementaridade off-target deve ser < 14 bp."""
        result = self._result()
        comp = result["data"]["complementarity_screen"]
        assert comp["max_complementarity_found"] < 14

    def test_complementarity_passes(self):
        """Teste de complementaridade deve passar."""
        result = self._result()
        comp = result["data"]["complementarity_screen"]
        assert comp["passes"] is True

    def test_off_target_risk_negligible(self):
        """Risco off-target deve ser negligible."""
        result = self._result()
        combined = result["data"]["combined_assessment"]
        assert combined["off_target_risk"] == "negligible"

    def test_conservation_analysis_present(self):
        """Dados de conservacao do SL RNA devem estar presentes."""
        result = self._result()
        conservation = result["data"]["conservation_analysis"]
        assert conservation["species_count"] >= 5


# =========================================================================
# Modulo 03: Conservacao evolutiva
# =========================================================================


class TestModule03Conservation:
    """Modulo 03: Conservacao evolutiva."""

    MODULE = "aso_math.03_evolutionary_conservation.run"

    def _result(self) -> dict:
        return _run_module(self.MODULE)

    def test_status_success(self):
        """Modulo deve completar com status success."""
        result = self._result()
        assert result["status"] == "success"

    def test_leishmania_invariant(self):
        """Regiao-alvo deve ser invariante em Leishmania."""
        result = self._result()
        target = result["data"]["target_region_assessment"]
        assert target["leishmania_invariant"] is True

    def test_purifying_selection(self):
        """Omega deve ser < 0.5 (selecao purificadora forte)."""
        result = self._result()
        div = result["data"]["divergence_analysis"]
        assert div["purifying_selection_ratio"] < 0.5

    def test_eleven_species(self):
        """Deve analisar 11 especies de kinetoplastideos."""
        result = self._result()
        species_data = result["data"]["species_data"]
        assert species_data["count"] == 11

    def test_conservation_score_high(self):
        """Score de conservacao da regiao-alvo deve ser >= 0.60.

        Nota: score e calculado sobre TODAS as 11 especies (incluindo
        Trypanosoma que divergiu ha ~350 Mya). Posicoes variaveis entre
        generos reduzem o score para ~0.68. O importante e que dentro
        de Leishmania a regiao e 100% invariante.
        """
        result = self._result()
        target = result["data"]["target_region_assessment"]
        assert target["conservation_score"] >= 0.60


# =========================================================================
# Modulo 04: Otimizacao exaustiva
# =========================================================================


class TestModule04Optimization:
    """Modulo 04: Otimizacao exaustiva."""

    MODULE = "aso_math.04_exhaustive_optimization.run"

    def _result(self) -> dict:
        return _run_module(self.MODULE)

    def test_status_success(self):
        """Modulo deve completar com status success."""
        result = self._result()
        assert result["status"] == "success"

    def test_pareto_includes_mrl_aso_001(self):
        """MRL-ASO-001 deve estar na frente de Pareto."""
        result = self._result()
        pareto = result["data"]["pareto_front"]
        assert pareto["mrl_in_pareto"] is True

    def test_evaluates_2800_designs(self):
        """Deve avaliar 2800 combinacoes (sequencia x LNA)."""
        result = self._result()
        total = result["data"]["full_enumeration"]["total_evaluated"]
        assert total == 2800

    def test_top_10_present(self):
        """Envelope deve conter top 10 designs rankeados."""
        result = self._result()
        top_10 = result["data"]["full_enumeration"]["top_10"]
        assert len(top_10) == 10

    def test_pareto_front_non_empty(self):
        """Frente de Pareto deve ter pelo menos 1 candidato."""
        result = self._result()
        pareto = result["data"]["pareto_front"]
        assert pareto["total_non_dominated"] >= 1


# =========================================================================
# Modulo 05: Modelo de resistencia
# =========================================================================


class TestModule05Resistance:
    """Modulo 05: Modelo de resistencia."""

    MODULE = "aso_math.05_resistance_model.run"

    def _result(self) -> dict:
        return _run_module(self.MODULE)

    def test_status_success(self):
        """Modulo deve completar com status success."""
        result = self._result()
        assert result["status"] == "success"

    def test_zero_viable_escape(self):
        """Nenhuma mutacao de escape viavel no modelo realista."""
        result = self._result()
        escape = result["data"]["escape_mutations"]
        assert escape["functionally_viable"] == 0

    def test_infinite_resistance_time(self):
        """Tempo de resistencia deve ser infinito para populacao 1e8."""
        result = self._result()
        time_data = result["data"]["time_to_resistance"]
        # Chave usa notacao cientifica: "population_1e+08"
        pop_1e8 = time_data.get("population_1e+08", {})
        years = pop_1e8.get("expected_years", None)
        # Deve ser "infinity" (string) porque JSON nao suporta inf
        assert years == "infinity"

    def test_four_population_sizes(self):
        """Deve modelar 4 tamanhos de populacao (1e3, 1e6, 1e8, 1e10)."""
        result = self._result()
        time_data = result["data"]["time_to_resistance"]
        assert len(time_data) == 4

    def test_sensitivity_analysis_present(self):
        """Dados da analise de sensibilidade devem estar presentes."""
        result = self._result()
        sensitivity = result["data"]["sensitivity_matrix"]
        assert sensitivity is not None
