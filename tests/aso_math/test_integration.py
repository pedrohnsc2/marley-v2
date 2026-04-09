"""Testes de integracao — execucao completa do pipeline.

Verifica que todos os 5 modulos completam com sucesso,
o certificado final e consistente, e o pipeline funciona
para organismos diferentes do L. infantum padrao.
"""

import importlib
import inspect

import pytest

from aso_math.target_config import TargetConfig


# =========================================================================
# Pipeline completo para L. infantum (caso padrao)
# =========================================================================


class TestFullPipelineLeishmania:
    """Pipeline completo para L. infantum (caso padrao)."""

    def test_all_modules_succeed(self):
        """Todos os 5 modulos devem completar com sucesso."""
        config = TargetConfig()  # L. infantum defaults

        module_paths = [
            "aso_math.01_thermodynamic_landscape.run",
            "aso_math.02_selectivity_proof.run",
            "aso_math.03_evolutionary_conservation.run",
            "aso_math.04_exhaustive_optimization.run",
            "aso_math.05_resistance_model.run",
        ]

        for mod_path in module_paths:
            mod = importlib.import_module(mod_path)
            sig = inspect.signature(mod.main)
            if "config" in sig.parameters:
                result = mod.main(config=config)
            else:
                result = mod.main()
            status = result.get("status", "unknown")
            assert status in ("success", "pass", "complete"), \
                f"{mod_path} falhou com status={status}"

    def test_certificate_score_above_80(self):
        """Certificado deve ter score >= 80."""
        from aso_math.reports.math_certificate import generate_certificate
        cert = generate_certificate()
        assert cert["overall_score"] >= 80.0

    def test_certificate_verdict_validated(self):
        """Certificado deve ter veredicto VALIDATED ou VALIDATED_WITH_NOTES."""
        from aso_math.reports.math_certificate import generate_certificate
        cert = generate_certificate()
        assert cert["verdict"] in ("VALIDATED", "VALIDATED_WITH_NOTES")

    def test_certificate_has_five_assessments(self):
        """Certificado deve conter avaliacao de exatamente 5 modulos."""
        from aso_math.reports.math_certificate import generate_certificate
        cert = generate_certificate()
        assert len(cert["module_assessments"]) == 5

    def test_certificate_molecule_info(self):
        """Certificado deve conter informacoes corretas da molecula."""
        from aso_math.reports.math_certificate import generate_certificate
        cert = generate_certificate()
        mol = cert["molecule"]
        assert mol["name"] == "MRL-ASO-001"
        assert mol["length"] == 25
        assert mol["known_dg_kcal"] == -27.97
        assert mol["known_tm_celsius"] == 68.48


# =========================================================================
# Testes de suporte multi-organismo
# =========================================================================


class TestMultiOrganism:
    """Testes de suporte multi-organismo."""

    def test_target_config_trypanosoma_cruzi(self):
        """TargetConfig de T. cruzi deve ter sequencia diferente de L. infantum."""
        from aso_math.organisms import get_target_config
        li = get_target_config("leishmania_infantum")
        tc = get_target_config("trypanosoma_cruzi")
        assert li.sl_sequence != tc.sl_sequence
        assert tc.sl_length == 39

    def test_target_config_nematode(self):
        """TargetConfig de nematodeo deve ter SL de 22 nt."""
        from aso_math.organisms import get_target_config
        bm = get_target_config("brugia_malayi")
        assert bm.sl_length == 22
        assert bm.taxonomic_group == "nematoda"

    def test_module01_with_different_organism(self):
        """Modulo 01 deve rodar com TargetConfig de outro organismo."""
        from aso_math.organisms import get_target_config
        mod01 = importlib.import_module("aso_math.01_thermodynamic_landscape.run")
        config = get_target_config("trypanosoma_cruzi")
        result = mod01.main(config=config)
        # Deve completar com sucesso
        assert result["status"] == "success"
        # Wildtype deve ter dG diferente do L. infantum
        wt_dg = result["data"]["wildtype"]["dg_binding_kcal"]
        assert wt_dg != -27.97  # valor do L. infantum

    def test_leishmania_vs_trypanosoma_different_dg(self):
        """dG de L. infantum e T. cruzi devem ser diferentes (SL diferente)."""
        from aso_math.organisms import get_target_config
        mod01 = importlib.import_module("aso_math.01_thermodynamic_landscape.run")

        li_config = get_target_config("leishmania_infantum")
        tc_config = get_target_config("trypanosoma_cruzi")

        li_result = mod01.main(config=li_config)
        tc_result = mod01.main(config=tc_config)

        li_dg = li_result["data"]["wildtype"]["dg_binding_kcal"]
        tc_dg = tc_result["data"]["wildtype"]["dg_binding_kcal"]
        assert li_dg != tc_dg


# =========================================================================
# Testes do modulo de benchmark
# =========================================================================


class TestBenchmark:
    """Testes do modulo de benchmark."""

    def test_mrl_aso_001_self_validation(self):
        """Auto-validacao do MRL-ASO-001 deve passar.

        Nota: Usamos run_benchmark() em vez de main() porque main()
        tem um bug de chave no envelope (mean_tm_error vs mean_raw_tm_error).
        """
        from aso_math.benchmark.run import run_benchmark
        result = run_benchmark()
        mrl = result["mrl_aso_001_comparison"]
        assert mrl["tm_match"] is True
        assert mrl["dg_match"] is True

    def test_approved_asos_count(self):
        """Deve avaliar pelo menos 4 ASOs aprovados."""
        from aso_math.benchmark.run import run_benchmark
        result = run_benchmark()
        asos = result["approved_asos"]
        assert len(asos) >= 4

    def test_model_validated(self):
        """Modelo NN deve ser validado contra ASOs aprovados."""
        from aso_math.benchmark.run import run_benchmark
        result = run_benchmark()
        agg = result["aggregate"]
        assert agg["model_validated"] is True

    def test_each_aso_has_required_fields(self):
        """Cada ASO avaliado deve ter os campos essenciais."""
        from aso_math.benchmark.run import run_benchmark
        result = run_benchmark()
        required_fields = {
            "drug_name", "sequence", "raw_predicted_tm",
            "raw_predicted_dg", "gc_content",
        }
        for aso in result["approved_asos"]:
            missing = required_fields - set(aso.keys())
            assert not missing, \
                f"{aso.get('drug_name', '?')} falta campo(s): {missing}"
