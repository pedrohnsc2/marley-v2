import { safeLoadModuleJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

interface QaoaConfig {
  bitstring: number[];
  fval_qubo: number;
  dg_total: number;
  probability: number;
  configuration: string;
  sp_count: number;
  rp_count: number;
  rp_triples: number;
}

interface QaoaTop10Data {
  metadata: { dg_base_kcal_mol: number; aso_sequence: string; aso_length: number };
  top10_configurations: QaoaConfig[];
  energy_unit: string;
  configuration_encoding: string;
}

interface ComparisonData {
  metadata: { dg_base_kcal_mol: number };
  comparison: {
    qaoa: { best_energy: number; method: string; sp_count: number; rp_count: number; wall_time_seconds: number };
    greedy: { best_energy: number; sp_count: number; rp_count: number };
  };
  detailed_comparison: {
    exhaustive_n10: { dg_optimal: number; configuration: string; sp_count: number };
    greedy_n10: { dg_found: number; configuration: string; sp_count: number; gap_vs_exhaustive: number };
    greedy_n25: { dg_found: number; configuration: string; sp_count: number };
    qaoa: { dg_found: number; configuration: string; sp_count: number; n_positions: number; n_layers: number; optimizer: string };
    qaoa_vs_exhaustive: { gap_kcal: number; gap_percent: number; qaoa_found_optimal: boolean };
  };
}

interface VqeData {
  mulliken_charges: { symbols: string[]; charges: number[] };
  total_charge: number;
  spin_multiplicity: number;
  notes: string;
}

function PendingCard({ label }: { label: string }) {
  return (
    <div className="rounded-xl border border-dashed border-gray-200 bg-gray-50 p-8 text-center">
      <div className="mx-auto mb-3 flex h-10 w-10 items-center justify-center rounded-full bg-indigo-100">
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.5} className="h-5 w-5 text-indigo-500">
          <path strokeLinecap="round" strokeLinejoin="round" d="M12 6v6l4 2m6-2a10 10 0 11-20 0 10 10 0 0120 0z" />
        </svg>
      </div>
      <p className="text-sm font-semibold text-gray-600">{label}</p>
      <p className="mt-1 text-xs text-gray-400">Computação em andamento — resultados disponíveis em breve</p>
    </div>
  );
}

export default function QuantumPage() {
  const qaoa = safeLoadModuleJson("qaoa", "qaoa_top10_configurations.json") as QaoaTop10Data | null;
  const comparison = safeLoadModuleJson("qaoa", "qaoa_vs_classical_comparison.json") as ComparisonData | null;
  const vqe = safeLoadModuleJson("vqe", "vqe_mulliken_charges.json") as VqeData | null;

  const bestConfig = qaoa?.top10_configurations?.[0];
  const dgBase = qaoa?.metadata?.dg_base_kcal_mol ?? 0;
  const bestEnergy = bestConfig?.dg_total ?? 0;
  const improvement = Math.abs(bestEnergy - dgBase);
  const qaoaFoundOptimal = comparison?.detailed_comparison?.qaoa_vs_exhaustive?.qaoa_found_optimal ?? false;

  const top5 = qaoa?.top10_configurations?.slice(0, 5) ?? [];
  const barCategories = top5.map((_, i) => `Config ${i + 1}`);
  const barSeries = [{ name: "dG (kcal/mol)", data: top5.map((c) => Math.abs(c.dg_total)) }];

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-indigo-100 px-2.5 py-1 text-xs font-bold text-indigo-600">QAOA+VQE</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Quantum Computing</h1>
          <p className="text-sm text-gray-500">
            QAOA and VQE algorithms for ASO phosphorothioate stereochemistry optimization
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Best Energy"
          value={qaoa ? `${bestEnergy.toFixed(2)} kcal/mol` : "—"}
          subtitle="Optimal QAOA configuration"
          accentColor="bg-indigo-500"
        />
        <KpiCard
          title="Optimal Config"
          value={bestConfig ? "All-Sp" : "—"}
          subtitle={bestConfig ? `${bestConfig.sp_count} Sp / ${bestConfig.rp_count} Rp` : "Pending"}
          accentColor="bg-indigo-400"
        />
        <KpiCard
          title="QAOA Accuracy"
          value={comparison
            ? qaoaFoundOptimal
              ? "100%"
              : `${(100 - comparison.detailed_comparison.qaoa_vs_exhaustive.gap_percent).toFixed(1)}%`
            : "—"}
          subtitle={qaoaFoundOptimal ? "Found global optimum" : "vs exhaustive search"}
          accentColor="bg-violet-500"
        />
        <KpiCard
          title="Improvement vs Base"
          value={qaoa ? `+${improvement.toFixed(2)} kcal/mol` : "—"}
          subtitle={qaoa ? `${bestEnergy.toFixed(2)} vs ${dgBase.toFixed(2)} base` : "Pending"}
          accentColor="bg-purple-500"
        />
      </div>

      {/* Main content */}
      {qaoa ? (
        <div className="grid gap-6 xl:grid-cols-3">
          {/* QAOA Top 10 table */}
          <div className="xl:col-span-2">
            <div className="rounded-xl bg-white shadow-card overflow-hidden">
              <div className="border-b border-gray-100 px-5 py-4">
                <h2 className="text-sm font-semibold text-gray-900">QAOA Top 10 Configurations</h2>
                <p className="text-xs text-gray-400 mt-0.5">
                  Ranked by free energy ({qaoa.energy_unit}) &mdash; {qaoa.configuration_encoding}
                </p>
              </div>
              <div className="overflow-x-auto">
                <table className="w-full text-left text-sm" data-testid="qaoa-top10-table">
                  <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                    <tr>
                      <th className="px-4 py-3">Rank</th>
                      <th className="px-4 py-3">Configuration</th>
                      <th className="px-4 py-3 text-right">dG (kcal/mol)</th>
                      <th className="px-4 py-3 text-right">Probability</th>
                      <th className="px-4 py-3 text-right">Sp</th>
                      <th className="px-4 py-3 text-right">Rp</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-gray-50">
                    {qaoa.top10_configurations.map((cfg, i) => (
                      <tr key={`${cfg.configuration}-${i}`} className="transition-colors hover:bg-gray-50">
                        <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                        <td className="px-4 py-2.5">
                          <span className="rounded-md bg-indigo-50 px-2 py-0.5 font-mono text-xs font-semibold text-indigo-700">
                            {cfg.configuration}
                          </span>
                        </td>
                        <td className="px-4 py-2.5 text-right font-mono text-xs font-semibold text-indigo-600">
                          {cfg.dg_total.toFixed(2)}
                        </td>
                        <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">
                          {(cfg.probability * 100).toFixed(2)}%
                        </td>
                        <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">{cfg.sp_count}</td>
                        <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">{cfg.rp_count}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          </div>

          {/* Right column */}
          <div className="space-y-6">
            <div className="rounded-xl bg-white shadow-card p-5">
              <h2 className="text-sm font-semibold text-gray-900">Top 5 Configurations</h2>
              <p className="text-xs text-gray-400 mt-0.5 mb-4">Free energy magnitude (|dG|, kcal/mol)</p>
              <BarChart categories={barCategories} series={barSeries} colors={["#6366F1"]} height={240} horizontal />
            </div>

            {comparison && (
              <div className="rounded-xl bg-gradient-to-br from-indigo-500 to-violet-500 p-5 text-white shadow-card">
                <div className="mb-3 flex items-center gap-2">
                  <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">
                    {qaoaFoundOptimal ? "Optimal Found" : "Near-Optimal"}
                  </span>
                </div>
                <p className="text-base font-bold">QAOA Result</p>
                <p className="mt-1 text-xs text-indigo-100">
                  {comparison.detailed_comparison.qaoa.n_positions} qubits, {comparison.detailed_comparison.qaoa.n_layers} layer, {comparison.detailed_comparison.qaoa.optimizer} optimizer
                </p>
                <div className="mt-4 flex items-baseline gap-2">
                  <span className="text-3xl font-bold">{bestEnergy.toFixed(2)}</span>
                  <span className="text-sm text-indigo-200">kcal/mol</span>
                </div>
                <p className="mt-1 text-xs text-indigo-100">
                  Gap vs exhaustive: {comparison.detailed_comparison.qaoa_vs_exhaustive.gap_kcal.toFixed(2)} kcal/mol ({comparison.detailed_comparison.qaoa_vs_exhaustive.gap_percent.toFixed(1)}%)
                </p>
                <div className="mt-4 border-t border-white/20 pt-4 grid grid-cols-2 gap-2 text-xs">
                  <div>
                    <p className="text-indigo-200">Wall Time</p>
                    <p className="font-semibold">{comparison.comparison.qaoa.wall_time_seconds.toFixed(1)}s</p>
                  </div>
                  <div>
                    <p className="text-indigo-200">Probability</p>
                    <p className="font-semibold">{(bestConfig!.probability * 100).toFixed(1)}%</p>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
      ) : (
        <PendingCard label="QAOA Top 10" />
      )}

      {/* Comparison section */}
      {comparison ? (
        <div className="mb-6 mt-6">
          <div className="rounded-xl bg-white shadow-card overflow-hidden">
            <div className="border-b border-gray-100 px-5 py-4">
              <h2 className="text-sm font-semibold text-gray-900">QAOA vs Classical Comparison</h2>
              <p className="text-xs text-gray-400 mt-0.5">Benchmark across exhaustive, greedy, and quantum approaches</p>
            </div>
            <div className="grid gap-4 p-5 sm:grid-cols-3" data-testid="comparison-cards">
              <div className="rounded-lg border border-gray-100 p-4">
                <div className="flex items-center gap-2">
                  <span className="inline-flex h-2 w-2 rounded-full bg-gray-400" />
                  <p className="text-xs font-semibold uppercase tracking-wider text-gray-400">Exhaustive (N=10)</p>
                </div>
                <p className="mt-3 text-2xl font-bold text-gray-900">{comparison.detailed_comparison.exhaustive_n10.dg_optimal.toFixed(2)}</p>
                <p className="text-xs text-gray-400">kcal/mol</p>
                <p className="mt-2 text-xs text-gray-500">{comparison.detailed_comparison.exhaustive_n10.sp_count} Sp &mdash; {comparison.detailed_comparison.exhaustive_n10.configuration}</p>
              </div>
              <div className="rounded-lg border border-gray-100 p-4">
                <div className="flex items-center gap-2">
                  <span className="inline-flex h-2 w-2 rounded-full bg-emerald-400" />
                  <p className="text-xs font-semibold uppercase tracking-wider text-gray-400">Greedy (N=25)</p>
                </div>
                <p className="mt-3 text-2xl font-bold text-gray-900">{comparison.detailed_comparison.greedy_n25.dg_found.toFixed(2)}</p>
                <p className="text-xs text-gray-400">kcal/mol</p>
                <p className="mt-2 text-xs text-gray-500">{comparison.detailed_comparison.greedy_n25.sp_count} Sp &mdash; full-length ASO</p>
              </div>
              <div className="rounded-lg border-2 border-indigo-200 bg-indigo-50/50 p-4">
                <div className="flex items-center gap-2">
                  <span className="inline-flex h-2 w-2 rounded-full bg-indigo-500" />
                  <p className="text-xs font-semibold uppercase tracking-wider text-indigo-500">QAOA (N=10)</p>
                </div>
                <p className="mt-3 text-2xl font-bold text-indigo-700">{comparison.detailed_comparison.qaoa.dg_found.toFixed(2)}</p>
                <p className="text-xs text-indigo-400">kcal/mol</p>
                <p className="mt-2 text-xs text-indigo-600">
                  {qaoaFoundOptimal ? "Matched exhaustive optimum" : `Gap: ${comparison.detailed_comparison.qaoa_vs_exhaustive.gap_kcal.toFixed(2)} kcal/mol`}
                </p>
              </div>
            </div>
          </div>
        </div>
      ) : (
        <div className="mt-6 mb-6"><PendingCard label="QAOA vs Classical Comparison" /></div>
      )}

      {/* VQE Mulliken Charges */}
      {vqe ? (
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">VQE Active Site Analysis &mdash; RNase H1 Catalytic Site</h2>
            <p className="text-xs text-gray-400 mt-0.5">Mulliken population analysis of the two-metal catalytic center</p>
          </div>
          <div className="p-5">
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm" data-testid="vqe-charges-table">
                <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                  <tr>
                    <th className="px-4 py-3">Atom</th>
                    <th className="px-4 py-3">Symbol</th>
                    <th className="px-4 py-3 text-right">Mulliken Charge</th>
                    <th className="px-4 py-3 text-center">Character</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-50">
                  {vqe.mulliken_charges.symbols.map((symbol, i) => {
                    const charge = vqe.mulliken_charges.charges[i];
                    const character = charge > 0 ? "Electrophilic" : charge < 0 ? "Nucleophilic" : "Neutral";
                    const badgeClass = charge > 0
                      ? "inline-flex items-center rounded-full bg-indigo-50 px-2 py-0.5 text-xs font-semibold text-indigo-600"
                      : charge < 0
                        ? "inline-flex items-center rounded-full bg-rose-50 px-2 py-0.5 text-xs font-semibold text-rose-600"
                        : "inline-flex items-center rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-400";
                    return (
                      <tr key={`${symbol}-${i}`} className="transition-colors hover:bg-gray-50">
                        <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                        <td className="px-4 py-2.5 font-semibold text-gray-900">{symbol}</td>
                        <td className="px-4 py-2.5 text-right font-mono text-xs font-semibold text-indigo-600">
                          {charge > 0 ? "+" : ""}{charge.toFixed(4)}
                        </td>
                        <td className="px-4 py-2.5 text-center"><span className={badgeClass}>{character}</span></td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
            <div className="mt-4 rounded-lg bg-indigo-50 p-4">
              <p className="text-xs font-semibold text-indigo-900">About this analysis</p>
              <p className="mt-1 text-xs text-indigo-700">
                VQE computes electronic structure of the RNase H1 active site. Total charge: {vqe.total_charge.toFixed(1)} | Spin multiplicity: {vqe.spin_multiplicity}.
              </p>
            </div>
          </div>
        </div>
      ) : (
        <PendingCard label="VQE Mulliken Charges" />
      )}
    </div>
  );
}
