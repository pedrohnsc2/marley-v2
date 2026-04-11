"""Loop de descoberta cientifica do AI Scientist.

Orquestra os 5 agentes em um ciclo de raciocinio:
    1. Cada agente analisa estado atual
    2. Cruzar achados entre agentes (consenso/divergencia)
    3. Gerar hipoteses priorizadas
    4. Propor plano de validacao experimental
    5. Escrever entrada no log de descoberta
    6. (Se Claude API disponivel) Gerar sintese narrativa com LLM

Agentes sao deterministicos. A sintese LLM e opcional (tier 2).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

import importlib as _il

# Imports via importlib — "11_scientist" começa com digito, invalido no parser.
_agents_mod = _il.import_module("marley_ai.11_scientist.agents")
_config_mod = _il.import_module("marley_ai.11_scientist.config")

AgentInsight = _agents_mod.AgentInsight
AgentOutput = _agents_mod.AgentOutput
DesignAgent = _agents_mod.DesignAgent
KnowledgeAgent = _agents_mod.KnowledgeAgent
LiteratureAgent = _agents_mod.LiteratureAgent
ReportAgent = _agents_mod.ReportAgent
ValidationAgent = _agents_mod.ValidationAgent
_is_operational = _agents_mod._is_operational
_load_module_result = _agents_mod._load_module_result

AGENT_ROLES = _config_mod.AGENT_ROLES
DiscoveryConfig = _config_mod.DiscoveryConfig
MODULE_RESULT_PATHS = _config_mod.MODULE_RESULT_PATHS


# ---------------------------------------------------------------------------
# Estruturas de dados do loop
# ---------------------------------------------------------------------------

@dataclass
class Hypothesis:
    """Uma hipotese cientifica gerada pelo loop de descoberta."""

    id: str
    title: str
    description: str
    confidence: float               # 0-1
    supporting_evidence: list[str]
    contradicting_evidence: list[str]
    testable_with: list[str]        # modulos que podem testar esta hipotese
    priority: int                   # 1 = mais urgente
    status: str = "proposed"        # proposed, supported, refuted

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "title": self.title,
            "description": self.description,
            "confidence": self.confidence,
            "supporting_evidence": self.supporting_evidence,
            "contradicting_evidence": self.contradicting_evidence,
            "testable_with": self.testable_with,
            "priority": self.priority,
            "status": self.status,
        }


@dataclass
class ProposedExperiment:
    """Um experimento proposto para validar uma hipotese."""

    id: str
    hypothesis_id: str
    title: str
    description: str
    experiment_type: str      # computational, wet_lab
    modules_required: list[str]
    estimated_effort: str     # low, medium, high
    expected_outcome: str
    priority: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "hypothesis_id": self.hypothesis_id,
            "title": self.title,
            "description": self.description,
            "experiment_type": self.experiment_type,
            "modules_required": self.modules_required,
            "estimated_effort": self.estimated_effort,
            "expected_outcome": self.expected_outcome,
            "priority": self.priority,
        }


@dataclass
class DiscoveryIteration:
    """Resultado de uma iteracao do loop de descoberta."""

    iteration: int
    timestamp: str
    agent_outputs: list[AgentOutput]
    hypotheses: list[Hypothesis]
    proposed_experiments: list[ProposedExperiment]
    consensus_score: float
    summary: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "iteration": self.iteration,
            "timestamp": self.timestamp,
            "agent_outputs": [ao.to_dict() for ao in self.agent_outputs],
            "n_hypotheses": len(self.hypotheses),
            "hypotheses": [h.to_dict() for h in self.hypotheses],
            "n_proposed_experiments": len(self.proposed_experiments),
            "proposed_experiments": [e.to_dict() for e in self.proposed_experiments],
            "consensus_score": self.consensus_score,
            "summary": self.summary,
        }


# ---------------------------------------------------------------------------
# Funcoes do loop
# ---------------------------------------------------------------------------

def load_all_results() -> dict[str, dict[str, Any] | None]:
    """Carrega resultados de todos os modulos 01-10.

    Retorna dict slug -> resultado (ou None se indisponivel).
    """
    results: dict[str, dict[str, Any] | None] = {}
    for slug, path in MODULE_RESULT_PATHS.items():
        results[slug] = _load_module_result(path)
    return results


def _run_agents(
    module_results: dict[str, dict[str, Any] | None],
) -> list[AgentOutput]:
    """Executa os 5 agentes em sequencia e coleta outputs.

    Ordem importa: ReportAgent roda por ultimo para consumir
    outputs dos outros 4 agentes.
    """
    # Fase 1: agentes independentes
    literature = LiteratureAgent(module_results).analyze()
    design = DesignAgent(module_results).analyze()
    validation = ValidationAgent(module_results).analyze()
    knowledge = KnowledgeAgent(module_results).analyze()

    # Fase 2: agente de sintese (depende dos anteriores)
    prior_outputs = [literature, design, validation, knowledge]
    report = ReportAgent(module_results, agent_outputs=prior_outputs).analyze()

    return [literature, design, validation, knowledge, report]


def _calculate_consensus(agent_outputs: list[AgentOutput]) -> float:
    """Calcula score de consenso entre agentes (0-1).

    Baseado em: proporcao de insights que sao validacoes vs warnings,
    ponderado pela confianca media.
    """
    all_insights: list[AgentInsight] = []
    for ao in agent_outputs:
        all_insights.extend(ao.insights)

    if not all_insights:
        return 0.0

    # Insights positivos: validation, cross_validation, synthesis
    positive_categories = {"validation", "cross_validation", "synthesis", "pattern"}
    n_positive = sum(1 for i in all_insights if i.category in positive_categories)

    # Insights negativos: warning, divergence
    negative_categories = {"warning", "divergence"}
    n_negative = sum(1 for i in all_insights if i.category in negative_categories)

    # Score: proporcao positiva, ajustada pela confianca media
    total_evaluative = n_positive + n_negative
    if total_evaluative == 0:
        return 0.5  # neutro se nao ha insights avaliativos

    positive_ratio = n_positive / total_evaluative
    mean_confidence = sum(i.confidence for i in all_insights) / len(all_insights)

    return round(positive_ratio * mean_confidence, 4)


def _generate_hypotheses(
    agent_outputs: list[AgentOutput],
    module_results: dict[str, dict[str, Any] | None],
    config: DiscoveryConfig,
) -> list[Hypothesis]:
    """Gera hipoteses baseadas nos insights dos agentes.

    Cada hipotese e derivada de padroes observados nos dados e
    gaps identificados pelos agentes.
    """
    hypotheses: list[Hypothesis] = []
    h_id = 0

    # Hipotese 1: sempre — o ASO MRL-ASO-001 e eficaz contra o SL RNA
    rna_fm = module_results.get("04_rna_fm")
    rosetta = module_results.get("05_rosettafold")
    supporting = []
    contradicting = []

    if _is_operational(rna_fm):
        disrupted = rna_fm.get("data", {}).get("structure", {}).get(
            "structural_impact", {}
        ).get("pairs_disrupted", 0)
        conservation = rna_fm.get("data", {}).get("conservation", {}).get(
            "mean_identity", 0
        )
        supporting.append(f"Disrupta {disrupted}/13 pares de base")
        supporting.append(f"SL RNA conservado {conservation:.1%} cross-species")

    if _is_operational(rosetta):
        dg = rosetta.get("metrics", {}).get("dg_predicted", 0)
        rnase = rosetta.get("metrics", {}).get("rnase_h_accessible", False)
        supporting.append(f"dG={dg:.2f} kcal/mol (strong binder)")
        if rnase:
            supporting.append("RNase H acessivel — mecanismo de clivagem confirmado")

    h_id += 1
    hypotheses.append(Hypothesis(
        id=f"H{h_id:03d}",
        title="ASO MRL-ASO-001 inibe trans-splicing em L. infantum via RNase H",
        description=(
            "O ASO gapmer MRL-ASO-001 liga o SL RNA de L. infantum com alta "
            "afinidade, disrupta sua estrutura secundaria nativa, e ativa "
            "clivagem pelo mecanismo RNase H1. Isso impede o trans-splicing "
            "de todos os mRNAs do parasita, resultando em morte celular."
        ),
        confidence=0.85 if supporting else 0.50,
        supporting_evidence=supporting,
        contradicting_evidence=contradicting,
        testable_with=["04_rna_fm", "05_rosettafold"],
        priority=1,
    ))

    # Hipotese 2: variantes RL tem perfil termodinamico superior
    rl = module_results.get("08_rl_ppo")
    if _is_operational(rl):
        rl_best = rl.get("data", {}).get("best_variants", [])
        baseline = rl.get("data", {}).get("baseline", {})
        if rl_best:
            best = rl_best[0]
            h_id += 1
            hypotheses.append(Hypothesis(
                id=f"H{h_id:03d}",
                title="Variantes RL-otimizadas tem binding mais forte mas risco off-target",
                description=(
                    f"A variante {best.get('sequence', 'N/A')} encontrada por RL "
                    f"tem dG={best.get('dg_binding', 0):.2f} kcal/mol vs baseline "
                    f"{baseline.get('dg_binding', 0):.2f}. Porem, GC={best.get('gc_content', 0):.2f} "
                    f"aumenta risco de off-target. Trade-off a ser resolvido."
                ),
                confidence=0.75,
                supporting_evidence=[
                    f"dg_improvement={best.get('dg_binding', 0) - baseline.get('dg_binding', 0):.2f} kcal/mol",
                    f"reward_improvement={rl.get('summary', {}).get('key_metrics', {}).get('improvement_pct', 0):.1f}%",
                ],
                contradicting_evidence=[
                    f"gc_content={best.get('gc_content', 0):.2f} > 0.55 (risco off-target)",
                    "Sequencia nao testada para impacto estrutural no SL RNA",
                ],
                testable_with=["04_rna_fm", "05_rosettafold", "08_rl_ppo"],
                priority=2,
            ))

    # Hipotese 3: perfil hidrofobico e essencial para epitopos vacinais
    sae = module_results.get("09_sae")
    if _is_operational(sae):
        features = sae.get("data", {}).get("interpreted_features", [])
        hydro_features = [
            f for f in features if "hydrophobicity" in f.get("label", "")
        ]
        if len(hydro_features) >= 3:
            h_id += 1
            hypotheses.append(Hypothesis(
                id=f"H{h_id:03d}",
                title="Hidrofobicidade N-terminal e preditor principal de antigenicidade em Leishmania",
                description=(
                    f"SAE identificou {len(hydro_features)} features hidrofobicas "
                    f"entre as top-{len(features)}. Features N-terminal dominam "
                    f"(neuronios 190, 170, 207). Hipotese: residuos hidrofobicos "
                    f"nas posicoes 1-3 sao necessarios para ancoragem no sulco "
                    f"MHC-I canino (alelos DLA-88)."
                ),
                confidence=0.70,
                supporting_evidence=[
                    f"n_hydro_features={len(hydro_features)}/{len(features)}",
                    "Consistente com ancoragem MHC-I pos 2/9",
                    "Top epitopos (YLAALVPAL, MILGTFVRL) iniciam com hidrofobicos",
                ],
                contradicting_evidence=[
                    "SAE treinado com 66 sequencias (dataset pequeno)",
                    f"R2={sae.get('metrics', {}).get('reconstruction_r2', 0):.3f} (reconstrucao moderada)",
                ],
                testable_with=["09_sae", "06_evodiff"],
                priority=3,
            ))

    # Hipotese 4: proteina hipotetica LINF_230010300 e fator de virulencia
    esm = module_results.get("03_leish_esm")
    if _is_operational(esm):
        hyp_data = esm.get("data", {}).get("epitope_analysis", {}).get(
            "source_hypothetical_conserved", {}
        )
        if hyp_data:
            h_id += 1
            hypotheses.append(Hypothesis(
                id=f"H{h_id:03d}",
                title="LINF_230010300 (proteina hipotetica conservada) e fator de virulencia",
                description=(
                    "Proteina hipotetica conservada (gene LINF_230010300) gera 3 "
                    "epitopos preditos com IC50 entre 33-65 nM. A divergencia ESM-2 "
                    f"e baixa ({hyp_data.get('divergence_score', 0):.4f}), sugerindo "
                    "funcao essencial conservada. Candidata a fator de virulencia "
                    "nao anotado."
                ),
                confidence=0.60,
                supporting_evidence=[
                    "3 epitopos preditos: SLMCVFYFK, YLAALVPAL, FAFSVSARR",
                    f"ESM divergence_score={hyp_data.get('divergence_score', 0):.4f}",
                    "Gene conservado entre Leishmania spp.",
                ],
                contradicting_evidence=[
                    "Funcao desconhecida — sem anotacao funcional",
                    "Nenhum dado experimental disponivel",
                ],
                testable_with=["02_leish_kg", "03_leish_esm"],
                priority=4,
            ))

    # Hipotese 5: EvoDiff + RL convergem em regiao quimica otima
    evodiff = module_results.get("06_evodiff")
    if _is_operational(evodiff) and _is_operational(rl):
        h_id += 1
        evo_best = evodiff.get("data", {}).get("aso", {}).get("top_candidates", [{}])
        rl_best_list = rl.get("data", {}).get("best_variants", [{}])
        hypotheses.append(Hypothesis(
            id=f"H{h_id:03d}",
            title="Regiao otima de design ASO esta em GC 0.40-0.50 e dG -30 a -33",
            description=(
                "EvoDiff (conservativo, GC~0.32-0.40) e RL (agressivo, GC~0.60) "
                "indicam uma regiao intermediaria de GC 0.40-0.50 como ideal: "
                "suficiente para binding forte (dG < -30) sem risco excessivo "
                "de off-target. O ASO ACAGAAACGGCTACTTTTGTAGCTT do EvoDiff "
                "(GC=0.40, score=0.89) exemplifica este trade-off."
            ),
            confidence=0.65,
            supporting_evidence=[
                f"EvoDiff top: score={evo_best[0].get('score', 'N/A') if evo_best else 'N/A'}",
                f"RL top: reward={rl_best_list[0].get('reward', 'N/A') if rl_best_list else 'N/A'}",
                "GC intermediario reduz off-target sem sacrificar dG",
            ],
            contradicting_evidence=[
                "Nao testado experimentalmente",
                "Trade-off GC/selectividade e especulativo",
            ],
            testable_with=["06_evodiff", "08_rl_ppo", "04_rna_fm"],
            priority=2,
        ))

    # Ordenar por prioridade e filtrar por confianca
    hypotheses.sort(key=lambda h: (h.priority, -h.confidence))
    return hypotheses[:config.max_hypotheses_per_iteration]


def _propose_experiments(
    hypotheses: list[Hypothesis],
    module_results: dict[str, dict[str, Any] | None],
    config: DiscoveryConfig,
) -> list[ProposedExperiment]:
    """Propoe experimentos para validar hipoteses."""
    experiments: list[ProposedExperiment] = []
    exp_id = 0

    for hyp in hypotheses[:config.max_proposed_experiments]:
        if hyp.id == "H001":
            exp_id += 1
            experiments.append(ProposedExperiment(
                id=f"EXP{exp_id:03d}",
                hypothesis_id=hyp.id,
                title="Validar ASO in vitro: RT-qPCR de trans-splicing em promastigotas",
                description=(
                    "Transfectar promastigotas de L. infantum com MRL-ASO-001 "
                    "(LNA gapmer, 1-10 uM). Medir nivel de mRNAs processados "
                    "por trans-splicing via RT-qPCR (alvo: GP63 mRNA). "
                    "Controles: scrambled ASO, mock."
                ),
                experiment_type="wet_lab",
                modules_required=[],
                estimated_effort="high",
                expected_outcome="Reducao >= 50% nos niveis de mRNA processado",
                priority=1,
            ))
        elif hyp.id == "H002":
            exp_id += 1
            experiments.append(ProposedExperiment(
                id=f"EXP{exp_id:03d}",
                hypothesis_id=hyp.id,
                title="Avaliar estrutura das variantes RL com RNA-FM e RoseTTAFold",
                description=(
                    "Submeter top-5 variantes RL ao pipeline RNA-FM + RoseTTAFold "
                    "para verificar se mantem disrupcao do SL RNA e acessibilidade "
                    "RNase H. Filtrar por BLAST contra transcriptoma canino."
                ),
                experiment_type="computational",
                modules_required=["04_rna_fm", "05_rosettafold"],
                estimated_effort="low",
                expected_outcome="Identificar 1-2 variantes com binding forte e baixo off-target",
                priority=2,
            ))
        elif hyp.id == "H003":
            exp_id += 1
            experiments.append(ProposedExperiment(
                id=f"EXP{exp_id:03d}",
                hypothesis_id=hyp.id,
                title="Validar epitopos preditos com ELISPOT canino",
                description=(
                    "Sintetizar top-5 epitopos (RMMRSLTPF, SLMCVFYFK, YLAALVPAL, "
                    "LIIEDLSLV, MILGTFVRL) e testar reatividade de celulas T caninas "
                    "de animais expostos a L. infantum via IFN-gamma ELISPOT."
                ),
                experiment_type="wet_lab",
                modules_required=[],
                estimated_effort="high",
                expected_outcome="Pelo menos 3 de 5 epitopos induzem resposta IFN-gamma",
                priority=3,
            ))
        elif hyp.id == "H004":
            exp_id += 1
            experiments.append(ProposedExperiment(
                id=f"EXP{exp_id:03d}",
                hypothesis_id=hyp.id,
                title="Anotar LINF_230010300 via homologia estrutural e KG",
                description=(
                    "Usar ESM-2 embeddings para buscar proteinas homologas anotadas "
                    "em outros organismos. Popular KG com triplas (LINF_230010300, "
                    "homolog_of, X) e inferir funcao."
                ),
                experiment_type="computational",
                modules_required=["02_leish_kg", "03_leish_esm"],
                estimated_effort="medium",
                expected_outcome="Anotacao funcional tentativa via transferencia de funcao",
                priority=4,
            ))
        elif hyp.id == "H005":
            exp_id += 1
            experiments.append(ProposedExperiment(
                id=f"EXP{exp_id:03d}",
                hypothesis_id=hyp.id,
                title="Gerar variantes ASO na faixa GC 0.40-0.50 com EvoDiff constrained",
                description=(
                    "Configurar EvoDiff com constraint de GC em [0.40, 0.50] e "
                    "gerar 100 variantes. Avaliar com metricas termodinamicas e "
                    "comparar com variantes RL na mesma faixa."
                ),
                experiment_type="computational",
                modules_required=["06_evodiff", "08_rl_ppo"],
                estimated_effort="low",
                expected_outcome="Variantes com score > 0.90 e GC balanceado",
                priority=2,
            ))

    experiments.sort(key=lambda e: e.priority)
    return experiments


# ---------------------------------------------------------------------------
# Loop principal
# ---------------------------------------------------------------------------

def run_discovery_loop(
    config: DiscoveryConfig | None = None,
) -> list[DiscoveryIteration]:
    """Executa o loop de descoberta cientifica.

    Carrega resultados de todos os modulos, executa os 5 agentes,
    gera hipoteses e propoe experimentos.

    Args:
        config: Parametros do loop. Se None, usa defaults.

    Returns:
        Lista de iteracoes executadas.
    """
    if config is None:
        config = DiscoveryConfig()

    # Carregar todos os resultados
    module_results = load_all_results()

    iterations: list[DiscoveryIteration] = []

    for i in range(config.max_iterations):
        timestamp = datetime.now(tz=timezone.utc).isoformat()

        # Passo 1: executar agentes
        agent_outputs = _run_agents(module_results)

        # Passo 2: calcular consenso
        consensus = _calculate_consensus(agent_outputs)

        # Passo 3: gerar hipoteses
        hypotheses = _generate_hypotheses(agent_outputs, module_results, config)

        # Passo 4: propor experimentos
        experiments = _propose_experiments(hypotheses, module_results, config)

        # Passo 5: criar entrada no log
        iteration = DiscoveryIteration(
            iteration=i + 1,
            timestamp=timestamp,
            agent_outputs=agent_outputs,
            hypotheses=hypotheses,
            proposed_experiments=experiments,
            consensus_score=consensus,
            summary=(
                f"Iteracao {i + 1}: {len(agent_outputs)} agentes executados, "
                f"{len(hypotheses)} hipoteses geradas, "
                f"{len(experiments)} experimentos propostos. "
                f"Consenso: {consensus:.2f}."
            ),
        )

        iterations.append(iteration)

        # Passo 6: sintese LLM (se Claude API disponivel)
        llm_synthesis = _generate_llm_synthesis(iteration)
        if llm_synthesis:
            iteration.summary += f"\n\n=== SINTESE CLAUDE ===\n{llm_synthesis}"

        # Verificar convergencia (para loops com mais de 1 iteracao)
        if consensus >= config.consensus_threshold and i > 0:
            break

    return iterations


# ---------------------------------------------------------------------------
# Tier 2: Sintese com Claude API (opcional)
# ---------------------------------------------------------------------------

def _detect_claude() -> bool:
    """Verifica se Claude API esta disponivel."""
    try:
        import anthropic  # noqa: F401
        import os
        from dotenv import load_dotenv
        load_dotenv()
        return bool(os.environ.get("ANTHROPIC_API_KEY"))
    except ImportError:
        return False


def _generate_llm_synthesis(iteration: "DiscoveryIteration") -> str | None:
    """Gera sintese narrativa usando Claude a partir dos resultados do loop.

    Envia hipoteses + experimentos + insights como contexto e pede
    uma sintese cientifica concisa. Custo: ~$0.003 por chamada (Haiku).

    Returns:
        Texto da sintese ou None se API indisponivel.
    """
    if not _detect_claude():
        return None

    try:
        import anthropic
        import os
        from dotenv import load_dotenv
        load_dotenv()

        # Montar contexto compacto dos achados
        context_parts = []

        # Hipoteses
        context_parts.append("HYPOTHESES:")
        for h in iteration.hypotheses:
            context_parts.append(
                f"  [{h.id}] (confidence={h.confidence:.2f}) {h.title}"
            )
            context_parts.append(f"    {h.description}")
            if h.supporting_evidence:
                context_parts.append(f"    Evidence: {'; '.join(h.supporting_evidence[:3])}")
            if h.contradicting_evidence:
                context_parts.append(f"    Caveats: {'; '.join(h.contradicting_evidence[:2])}")

        # Experimentos
        context_parts.append("\nPROPOSED EXPERIMENTS:")
        for e in iteration.proposed_experiments:
            context_parts.append(
                f"  [{e.id}] ({e.experiment_type}) {e.title}"
            )

        # Insights chave dos agentes
        context_parts.append("\nKEY AGENT INSIGHTS:")
        for ao in iteration.agent_outputs:
            for ins in ao.insights[:3]:  # top 3 por agente
                context_parts.append(
                    f"  [{ao.agent_name}] {ins.title} (conf={ins.confidence:.2f})"
                )

        context = "\n".join(context_parts)

        system_prompt = (
            "You are an AI scientist analyzing computational drug discovery results "
            "for MRL-ASO-001, an antisense oligonucleotide targeting the spliced leader "
            "RNA of Leishmania infantum (canine leishmaniasis). "
            "Write a concise scientific synthesis (max 300 words) that: "
            "1) Summarizes the key findings across all computational modules, "
            "2) Highlights the strongest hypothesis and its evidence, "
            "3) Identifies the most critical next experiment, "
            "4) Notes any contradictions or risks. "
            "Be scientifically rigorous. Cite specific numbers from the data."
        )

        client = anthropic.Anthropic()
        response = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=500,
            system=system_prompt,
            messages=[{"role": "user", "content": context}],
        )
        return response.content[0].text

    except Exception as e:
        print(f"[11_scientist] AVISO: Claude synthesis falhou ({e}). Usando modo deterministico.")
        return None
