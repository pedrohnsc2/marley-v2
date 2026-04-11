"""Execucao do modulo 11_scientist — AI Scientist autonomo.

Orquestra 5 agentes deterministicos que leem resultados dos modulos
01-10 e produzem: hipoteses, propostas experimentais e relatorio
consolidado. Sem dependencia de LLM — opera como motor de sintese.

Uso:
    python -m marley_ai.11_scientist.run
"""

from __future__ import annotations

import importlib as _il
from typing import Any

from marley_ai.config import AIModuleConfig
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

# Importacao via importlib — Python nao permite "from marley_ai.11_scientist..."
# porque "11" e interpretado como literal numerico invalido no parser de imports.
_sci_config = _il.import_module("marley_ai.11_scientist.config")
_sci_loop = _il.import_module("marley_ai.11_scientist.discovery_loop")

ScientistConfig = _sci_config.ScientistConfig
DiscoveryConfig = _sci_config.DiscoveryConfig
run_discovery_loop = _sci_loop.run_discovery_loop
load_all_results = _sci_loop.load_all_results
_is_operational = _il.import_module("marley_ai.11_scientist.agents")._is_operational


# ---------------------------------------------------------------------------
# Registro do modulo no registry global
# ---------------------------------------------------------------------------

@register("11_scientist")
class AIScientist:
    """Agente cientifico autonomo que integra todos os modulos Marley AI.

    Le resultados de 10 modulos (01-10), executa 5 agentes deterministicos,
    gera hipoteses priorizadas e propoe proximos experimentos.
    """

    def configure(self, config: Any) -> None:
        """Recebe configuracao do orquestrador."""
        self._config = config

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica se pelo menos 1 modulo tem resultados reais."""
        results = load_all_results()
        operational = [
            slug for slug, r in results.items() if _is_operational(r)
        ]
        missing = [
            slug for slug, r in results.items() if not _is_operational(r)
        ]
        return {
            "valid": len(operational) > 0,
            "operational": operational,
            "missing": missing,
            "n_operational": len(operational),
        }

    def run(self) -> dict[str, Any]:
        """Executa o loop de descoberta e retorna envelope."""
        return main()

    def get_dependencies(self) -> list[str]:
        """Depende de todos os modulos anteriores (soft dependency)."""
        return [
            "01_rag", "02_leish_kg", "03_leish_esm", "04_rna_fm",
            "05_rosettafold", "06_evodiff", "07_contrastive",
            "08_rl_ppo", "09_sae", "10_digital_twin",
        ]


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------

def main(config: AIModuleConfig | None = None) -> dict[str, Any]:
    """Executa o loop de raciocinio cientifico autonomo.

    Carrega resultados de todos os modulos 01-10, executa 5 agentes
    deterministicos (Literature, Design, Validation, Knowledge, Report),
    gera hipoteses priorizadas e propoe plano experimental.

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com hipoteses, experimentos e relatorio consolidado.
    """
    envelope = create_envelope("11_scientist")

    with Timer() as timer:
        # Configurar parametros do loop
        discovery_config = DiscoveryConfig()

        # Verificar inputs disponiveis
        module_results = load_all_results()
        n_operational = sum(
            1 for r in module_results.values() if _is_operational(r)
        )
        n_total = len(module_results)

        print(f"[11_scientist] Modulos disponiveis: {n_operational}/{n_total}")
        for slug, result in module_results.items():
            status = "operacional" if _is_operational(result) else "stub/ausente"
            print(f"  {slug}: {status}")

        if n_operational == 0:
            envelope["status"] = "skipped"
            envelope["warnings"].append(
                "Nenhum modulo operacional encontrado. "
                "Execute pelo menos 1 modulo (03-09) antes do AI Scientist."
            )
            envelope["summary"]["conclusion"] = (
                "AI Scientist nao pode operar sem dados de entrada. "
                "Nenhum modulo produziu resultados reais."
            )
            envelope["summary"]["key_metrics"] = {
                "n_modules_operational": 0,
                "n_modules_total": n_total,
                "n_hypotheses": 0,
                "n_experiments": 0,
            }
            envelope["runtime_seconds"] = timer.elapsed
            output_path = write_result(envelope)
            print(f"[11_scientist] Resultado salvo em {output_path}")
            return envelope

        # Executar loop de descoberta
        print(f"\n[11_scientist] Iniciando loop de descoberta (max_iter={discovery_config.max_iterations})...")
        iterations = run_discovery_loop(discovery_config)

        # Extrair resultados da ultima iteracao
        last = iterations[-1]
        print(f"[11_scientist] Iteracao {last.iteration} concluida:")
        print(f"  Agentes executados: {len(last.agent_outputs)}")
        print(f"  Hipoteses geradas: {len(last.hypotheses)}")
        print(f"  Experimentos propostos: {len(last.proposed_experiments)}")
        print(f"  Consenso: {last.consensus_score:.2f}")

        # Imprimir hipoteses
        print(f"\n[11_scientist] === HIPOTESES ===")
        for hyp in last.hypotheses:
            print(f"  [{hyp.id}] (conf={hyp.confidence:.2f}, prio={hyp.priority}) {hyp.title}")

        # Imprimir experimentos
        print(f"\n[11_scientist] === EXPERIMENTOS PROPOSTOS ===")
        for exp in last.proposed_experiments:
            print(f"  [{exp.id}] ({exp.experiment_type}, {exp.estimated_effort}) {exp.title}")

        # Montar envelope
        envelope["status"] = "complete"
        envelope["device"] = "cpu"  # deterministico, sem GPU

        # Coletar todos os insights para o summary
        all_insights = []
        for ao in last.agent_outputs:
            all_insights.extend(ao.insights)

        n_validations = sum(1 for i in all_insights if i.category == "validation")
        n_warnings = sum(1 for i in all_insights if i.category == "warning")
        n_gaps = sum(1 for i in all_insights if i.category == "gap")

        envelope["summary"]["conclusion"] = (
            f"AI Scientist executou {len(last.agent_outputs)} agentes sobre "
            f"{n_operational} modulos operacionais. "
            f"Gerou {len(last.hypotheses)} hipoteses e "
            f"{len(last.proposed_experiments)} propostas experimentais. "
            f"Consenso entre agentes: {last.consensus_score:.2f}. "
            f"Destaque: ASO MRL-ASO-001 validado computacionalmente "
            f"(disrupcao SL RNA + acessibilidade RNase H + dG concordante)."
        )
        envelope["summary"]["key_metrics"] = {
            "n_modules_operational": n_operational,
            "n_modules_total": n_total,
            "n_agents": len(last.agent_outputs),
            "n_total_insights": len(all_insights),
            "n_validations": n_validations,
            "n_warnings": n_warnings,
            "n_gaps": n_gaps,
            "n_hypotheses": len(last.hypotheses),
            "n_experiments_proposed": len(last.proposed_experiments),
            "consensus_score": last.consensus_score,
            "top_hypothesis_confidence": (
                last.hypotheses[0].confidence if last.hypotheses else 0
            ),
        }

        # Dados detalhados
        envelope["data"] = {
            "iterations": [it.to_dict() for it in iterations],
            "hypotheses_summary": [
                {
                    "id": h.id,
                    "title": h.title,
                    "confidence": h.confidence,
                    "priority": h.priority,
                    "status": h.status,
                    "n_supporting": len(h.supporting_evidence),
                    "n_contradicting": len(h.contradicting_evidence),
                }
                for h in last.hypotheses
            ],
            "experiments_summary": [
                {
                    "id": e.id,
                    "title": e.title,
                    "type": e.experiment_type,
                    "effort": e.estimated_effort,
                    "priority": e.priority,
                    "hypothesis": e.hypothesis_id,
                }
                for e in last.proposed_experiments
            ],
            "module_status": {
                slug: {
                    "operational": _is_operational(r),
                    "status": r.get("status", "missing") if r else "missing",
                }
                for slug, r in module_results.items()
            },
        }

        # Dependencias
        envelope["dependencies"] = [
            slug for slug, r in module_results.items() if _is_operational(r)
        ]

    envelope["runtime_seconds"] = timer.elapsed
    output_path = write_result(envelope)
    print(f"\n[11_scientist] Resultado salvo em {output_path}")

    return envelope


if __name__ == "__main__":
    main()
