import { safeLoadModuleJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import DonutChart from "@/components/charts/donut-chart";

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

interface ModuleEnvelope {
  module: string;
  version: string;
  generated_at: string;
  runtime_seconds: number;
  device: string;
  status: string;
  warnings: string[];
  summary: {
    conclusion: string;
    key_metrics: Record<string, number | string>;
  };
  dependencies: string[];
  data: Record<string, unknown>;
}

interface AgentOutput {
  agent_name: string;
  n_insights: number;
  status: string;
  summary: string;
}

interface Hypothesis {
  id: string;
  title: string;
  description: string;
  confidence: number;
  supporting_evidence: string[];
  contradicting_evidence: string[];
  priority: number;
  status: string;
}

interface ScientistIteration {
  iteration: number;
  agent_outputs: AgentOutput[];
  hypotheses: Hypothesis[];
}

interface ScientistData extends ModuleEnvelope {
  summary: {
    conclusion: string;
    key_metrics: {
      n_modules_operational: number;
      n_modules_total: number;
      n_agents: number;
      n_total_insights: number;
      n_validations: number;
      n_warnings: number;
      n_gaps: number;
      n_hypotheses: number;
      n_experiments_proposed: number;
      consensus_score: number;
      top_hypothesis_confidence: number;
    };
  };
  data: {
    iterations: ScientistIteration[];
  };
}

/* ------------------------------------------------------------------ */
/*  Module metadata                                                    */
/* ------------------------------------------------------------------ */

const AI_MODULES = [
  { id: "01_rag", name: "RAG", desc: "Retrieval-Augmented Generation" },
  { id: "02_leish_kg", name: "Knowledge Graph", desc: "Leishmania Knowledge Graph" },
  { id: "03_leish_esm", name: "ESM Protein", desc: "Protein Language Model" },
  { id: "04_rna_fm", name: "RNA Foundation", desc: "RNA Structure Prediction" },
  { id: "05_rosettafold", name: "RoseTTAFold", desc: "Protein Structure" },
  { id: "06_evodiff", name: "EvoDiff", desc: "Evolutionary Protein Design" },
  { id: "07_contrastive", name: "Contrastive", desc: "Molecular Similarity" },
  { id: "08_rl_ppo", name: "RL-PPO", desc: "Reinforcement Learning" },
  { id: "09_sae", name: "SAE", desc: "Sparse Autoencoder Interpretability" },
  { id: "10_digital_twin", name: "Digital Twin", desc: "In-silico Patient Simulation" },
  { id: "11_scientist", name: "AI Scientist", desc: "Autonomous Hypothesis Generation" },
];

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

function truncate(text: string, max: number): string {
  if (text.length <= max) return text;
  return text.slice(0, max) + "...";
}

function statusBadge(status: string | null) {
  if (!status) {
    return (
      <span className="rounded-full bg-gray-100 px-2.5 py-0.5 text-xs font-semibold text-gray-500">
        No data
      </span>
    );
  }
  const s = status.toLowerCase();
  if (s === "complete" || s === "success") {
    return (
      <span className="app-badge-success rounded-full px-2.5 py-0.5 text-xs font-semibold">
        {status}
      </span>
    );
  }
  if (s === "stub") {
    return (
      <span className="rounded-full bg-amber-50 px-2.5 py-0.5 text-xs font-semibold text-amber-700">
        stub
      </span>
    );
  }
  return (
    <span className="rounded-full bg-gray-100 px-2.5 py-0.5 text-xs font-semibold text-gray-600">
      {status}
    </span>
  );
}

function confidenceBadge(confidence: number) {
  const pct = (confidence * 100).toFixed(0);
  if (confidence >= 0.8) {
    return (
      <span className="rounded-full bg-emerald-50 px-2.5 py-0.5 text-xs font-semibold text-emerald-700">
        {pct}%
      </span>
    );
  }
  if (confidence >= 0.5) {
    return (
      <span className="rounded-full bg-amber-50 px-2.5 py-0.5 text-xs font-semibold text-amber-700">
        {pct}%
      </span>
    );
  }
  return (
    <span className="rounded-full bg-red-50 px-2.5 py-0.5 text-xs font-semibold text-red-700">
      {pct}%
    </span>
  );
}

/* ------------------------------------------------------------------ */
/*  Page                                                               */
/* ------------------------------------------------------------------ */

export default function AiPage() {
  /* Load all module JSONs */
  const moduleData: Record<string, ModuleEnvelope | null> = {};
  for (const mod of AI_MODULES) {
    const filename = `${mod.id}.json`;
    moduleData[mod.id] = safeLoadModuleJson("marley_ai", filename) as ModuleEnvelope | null;
  }

  /* Load AI Scientist data specifically */
  const scientist = moduleData["11_scientist"] as ScientistData | null;
  const metrics = scientist?.summary?.key_metrics;

  /* Compute KPI values */
  const operationalCount = AI_MODULES.filter((m) => {
    const d = moduleData[m.id];
    return d !== null && d.status !== undefined;
  }).length;

  const hypothesesCount = metrics?.n_hypotheses ?? 0;
  const consensusScore = metrics?.consensus_score ?? 0;
  const topConfidence = metrics?.top_hypothesis_confidence ?? 0;

  /* Donut chart: status distribution */
  let completeCount = 0;
  let stubCount = 0;
  let missingCount = 0;

  for (const mod of AI_MODULES) {
    const d = moduleData[mod.id];
    if (!d) {
      missingCount++;
    } else if (d.status === "complete" || d.status === "success") {
      completeCount++;
    } else if (d.status === "stub") {
      stubCount++;
    } else {
      missingCount++;
    }
  }

  const donutSeries = [completeCount, stubCount, missingCount].filter((v) => v > 0);
  const donutLabels: string[] = [];
  const donutColors: string[] = [];
  if (completeCount > 0) {
    donutLabels.push(`Complete (${completeCount})`);
    donutColors.push("#10B981");
  }
  if (stubCount > 0) {
    donutLabels.push(`Stub (${stubCount})`);
    donutColors.push("#F59E0B");
  }
  if (missingCount > 0) {
    donutLabels.push(`No data (${missingCount})`);
    donutColors.push("#9CA3AF");
  }

  /* Hypotheses from AI Scientist */
  const iteration = scientist?.data?.iterations?.[0];
  const hypotheses = (iteration?.hypotheses ?? [])
    .slice()
    .sort((a, b) => b.confidence - a.confidence);

  /* Agent outputs from AI Scientist */
  const agents = iteration?.agent_outputs ?? [];

  return (
    <div>
      {/* ---- Page header ---- */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-sky-100 px-2.5 py-1 text-xs font-bold text-sky-600">
          v0.1
        </span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">AI / ML Engine</h1>
          <p className="text-sm text-gray-500">
            11 computational modules for drug discovery and vaccine design against{" "}
            <span className="italic">L.&nbsp;infantum</span>
          </p>
        </div>
      </div>

      {/* ---- KPI row ---- */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Modules Loaded"
          value={`${operationalCount} / ${AI_MODULES.length}`}
          subtitle={`${completeCount} complete, ${stubCount} stubs`}
          accentColor="bg-sky-500"
        />
        <KpiCard
          title="Hypotheses"
          value={hypothesesCount}
          subtitle="Generated by AI Scientist"
          accentColor="bg-blue-500"
        />
        <KpiCard
          title="Consensus Score"
          value={consensusScore.toFixed(2)}
          subtitle="Inter-agent agreement"
          accentColor="bg-sky-500"
        />
        <KpiCard
          title="Top Confidence"
          value={`${(topConfidence * 100).toFixed(0)}%`}
          subtitle="Highest hypothesis confidence"
          accentColor="bg-blue-500"
        />
      </div>

      {/* ---- Module grid + Donut chart row ---- */}
      <div className="mb-6 grid gap-6 lg:grid-cols-3">
        {/* Module grid: 2 columns in the wider area */}
        <div className="lg:col-span-2">
          <div className="grid gap-4 sm:grid-cols-2">
            {AI_MODULES.map((mod) => {
              const d = moduleData[mod.id];
              const status = d?.status ?? null;
              const conclusion = d?.summary?.conclusion ?? null;
              const runtime = d?.runtime_seconds ?? null;
              const device = d?.device ?? null;

              return (
                <div
                  key={mod.id}
                  className="rounded-xl bg-white shadow-card p-5"
                  data-testid={`module-card-${mod.id}`}
                >
                  <div className="mb-2 flex items-start justify-between gap-2">
                    <div className="min-w-0">
                      <p className="text-sm font-bold text-gray-900">{mod.name}</p>
                      <p className="text-xs text-gray-400">{mod.desc}</p>
                    </div>
                    {statusBadge(status)}
                  </div>

                  {(runtime !== null && runtime > 0) || (device && device.length > 0) ? (
                    <div className="mb-2 flex items-center gap-3 text-xs text-gray-400">
                      {runtime !== null && runtime > 0 && (
                        <span>{runtime.toFixed(1)}s</span>
                      )}
                      {device && device.length > 0 && (
                        <span className="rounded bg-gray-100 px-1.5 py-0.5 font-mono text-xs text-gray-500">
                          {device}
                        </span>
                      )}
                    </div>
                  ) : null}

                  {conclusion && (
                    <p className="text-xs leading-relaxed text-gray-500">
                      {truncate(conclusion, 100)}
                    </p>
                  )}
                </div>
              );
            })}
          </div>
        </div>

        {/* Donut chart */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Module Status</h2>
          <p className="mt-0.5 mb-4 text-xs text-gray-400">
            Distribution across {AI_MODULES.length} AI/ML modules
          </p>
          <DonutChart
            series={donutSeries}
            labels={donutLabels}
            colors={donutColors}
            height={280}
          />
        </div>
      </div>

      {/* ---- Hypotheses card ---- */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">AI Scientist Hypotheses</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {hypotheses.length} hypotheses generated by autonomous multi-agent reasoning, sorted by confidence
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="hypotheses-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">Hypothesis</th>
                <th className="px-4 py-3 text-center">Confidence</th>
                <th className="px-4 py-3">Description</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {hypotheses.map((h, i) => (
                <tr
                  key={h.id}
                  className="transition-colors hover:bg-gray-50"
                >
                  <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                  <td className="max-w-xs px-4 py-2.5 text-sm font-medium text-gray-800">
                    {h.title}
                  </td>
                  <td className="px-4 py-2.5 text-center">
                    {confidenceBadge(h.confidence)}
                  </td>
                  <td className="max-w-md px-4 py-2.5 text-xs text-gray-500">
                    {truncate(h.description, 120)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* ---- Agents card ---- */}
      <div className="rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">AI Scientist Agents</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {agents.length} autonomous agents contributing insights to hypothesis generation
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="agents-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">Agent</th>
                <th className="px-4 py-3 text-right">Insights</th>
                <th className="px-4 py-3">Status</th>
                <th className="px-4 py-3">Summary</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {agents.map((agent) => (
                <tr
                  key={agent.agent_name}
                  className="transition-colors hover:bg-gray-50"
                >
                  <td className="px-4 py-2.5">
                    <span className="rounded-md bg-sky-50 px-2 py-0.5 font-mono text-xs font-semibold text-sky-700">
                      {agent.agent_name}
                    </span>
                  </td>
                  <td className="px-4 py-2.5 text-right font-mono text-sm font-bold text-gray-700">
                    {agent.n_insights}
                  </td>
                  <td className="px-4 py-2.5">
                    <span className="app-badge-success rounded-full px-2.5 py-0.5 text-xs font-semibold">
                      {agent.status}
                    </span>
                  </td>
                  <td className="max-w-md px-4 py-2.5 text-xs text-gray-500">
                    {truncate(agent.summary, 120)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}
