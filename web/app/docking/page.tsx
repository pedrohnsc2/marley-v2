import { loadCsv } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

export default function DockingPage() {
  const allScores = loadCsv("docking_scores.csv");
  const top10 = allScores.slice(0, 10);

  const bestAffinity =
    allScores.length > 0
      ? Math.min(...allScores.map((s) => parseFloat(s.binding_affinity ?? "0")))
      : 0;

  const uniqueTargets = new Set(allScores.map((s) => s.target_gene_name)).size;
  const uniqueCompounds = new Set(allScores.map((s) => s.compound_id)).size;

  // Chart data: top 10 binding affinities (negative = better, show as positive for visual clarity)
  const chartCategories = top10.map((r) => r.compound_name ?? r.compound_id ?? "");
  const chartSeries = [
    {
      name: "Binding Affinity",
      data: top10.map((r) => Math.abs(parseFloat(r.binding_affinity ?? "0"))),
    },
  ];

  const selectivity = [
    { compound: "MRL-003", tryr: -7.32, humanGr: -8.68, delta: -1.35, selective: false },
    { compound: "Pemetrexed", tryr: -6.92, humanGr: -9.28, delta: -2.36, selective: false },
    { compound: "Menadione", tryr: -4.99, humanGr: -6.5, delta: -1.51, selective: false },
  ];

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-emerald-100 px-2.5 py-1 text-xs font-bold text-emerald-600">v4</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Molecular Docking</h1>
          <p className="text-sm text-gray-500">
            AutoDock Vina virtual screening of repurposed and custom compounds
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Total Dockings"
          value={allScores.length}
          subtitle="Target-compound pairs"
          accentColor="bg-emerald-500"
        />
        <KpiCard
          title="Unique Targets"
          value={uniqueTargets}
          subtitle="Receptor structures"
          accentColor="bg-teal-500"
        />
        <KpiCard
          title="Compounds"
          value={uniqueCompounds}
          subtitle="Tested molecules"
          accentColor="bg-green-500"
        />
        <KpiCard
          title="Best Affinity"
          value={`${bestAffinity.toFixed(2)} kcal/mol`}
          subtitle={top10[0]?.compound_name ?? ""}
          accentColor="bg-lime-500"
        />
      </div>

      {/* Table + Chart */}
      <div className="mb-6 grid gap-6 xl:grid-cols-5">
        {/* Top hits table */}
        <div className="xl:col-span-3 rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">Top 10 Docking Hits</h2>
            <p className="text-xs text-gray-400 mt-0.5">Ranked by binding affinity (kcal/mol)</p>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-left text-sm" data-testid="docking-table">
              <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                <tr>
                  <th className="px-4 py-3">#</th>
                  <th className="px-4 py-3">Target</th>
                  <th className="px-4 py-3">Compound</th>
                  <th className="px-4 py-3 text-right">Affinity</th>
                  <th className="px-4 py-3 text-right">RMSD LB</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-50">
                {top10.map((row, i) => (
                  <tr
                    key={`${row.target_gene_id}-${row.compound_id}-${i}`}
                    className="transition-colors hover:bg-gray-50"
                  >
                    <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                    <td className="px-4 py-2.5 font-semibold text-gray-900">{row.target_gene_name}</td>
                    <td className="px-4 py-2.5">
                      <span className="rounded-md bg-emerald-50 px-2 py-0.5 font-mono text-xs text-emerald-700">
                        {row.compound_name}
                      </span>
                    </td>
                    <td className="px-4 py-2.5 text-right font-mono text-sm font-bold text-emerald-600">
                      {parseFloat(row.binding_affinity ?? "0").toFixed(2)}
                    </td>
                    <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-400">
                      {parseFloat(row.rmsd_lb ?? "0").toFixed(1)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        {/* Affinity chart */}
        <div className="xl:col-span-2 rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Binding Affinities</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">|kcal/mol| — top 10 compounds</p>
          <BarChart
            categories={chartCategories}
            series={chartSeries}
            colors={["#10B981"]}
            height={280}
            horizontal={true}
          />
        </div>
      </div>

      {/* Selectivity table */}
      <div className="rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Selectivity Analysis</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            L. infantum TryR vs human glutathione reductase (GR). Selective threshold: 1.5 kcal/mol difference.
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-5 py-3">Compound</th>
                <th className="px-5 py-3 text-right">TryR (kcal/mol)</th>
                <th className="px-5 py-3 text-right">Human GR (kcal/mol)</th>
                <th className="px-5 py-3 text-right">Delta</th>
                <th className="px-5 py-3 text-center">Selective</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {selectivity.map((s) => (
                <tr key={s.compound} className="hover:bg-gray-50 transition-colors">
                  <td className="px-5 py-3 font-semibold text-gray-900">{s.compound}</td>
                  <td className="px-5 py-3 text-right font-mono text-sm text-emerald-600">{s.tryr}</td>
                  <td className="px-5 py-3 text-right font-mono text-sm text-gray-500">{s.humanGr}</td>
                  <td className="px-5 py-3 text-right font-mono text-sm text-amber-600">{s.delta}</td>
                  <td className="px-5 py-3 text-center">
                    <span
                      className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-semibold ${
                        s.selective
                          ? "bg-emerald-50 text-emerald-700"
                          : "bg-red-50 text-red-600"
                      }`}
                    >
                      {s.selective ? "Yes" : "No"}
                    </span>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <div className="px-5 py-3 bg-amber-50 border-t border-amber-100">
          <p className="text-xs text-amber-700">
            None of the tested compounds showed parasite selectivity. Further structural optimization is needed.
          </p>
        </div>
      </div>

      {/* Docking pose images */}
      <div className="mt-6 grid gap-6 lg:grid-cols-2">
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">GMPS + Methotrexate Pose</h2>
            <p className="text-xs text-gray-400 mt-0.5">Best docking pose in the GMPS binding pocket · Rendered with PyMOL</p>
          </div>
          {/* eslint-disable-next-line @next/next/no-img-element */}
          <img src="/images/docking_3d.png" alt="Docking pose of methotrexate in the GMPS active site" className="w-full" />
        </div>
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">TryR + MRL-003 Docking</h2>
            <p className="text-xs text-gray-400 mt-0.5">Custom molecule MRL-003 (-7.74 kcal/mol) in the TryR binding site</p>
          </div>
          {/* eslint-disable-next-line @next/next/no-img-element */}
          <img src="/images/tryr_mrl003.png" alt="MRL-003 docked into trypanothione reductase active site" className="w-full" />
        </div>
      </div>
    </div>
  );
}
