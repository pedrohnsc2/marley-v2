import { loadCsv } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

export default function DrugPage() {
  const targets = loadCsv("drug_targets_top20.csv");

  // Count targets by pathway
  const pathwayCounts: Record<string, number> = {};
  for (const t of targets) {
    const pw = t.pathway ?? "unknown";
    pathwayCounts[pw] = (pathwayCounts[pw] ?? 0) + 1;
  }
  const pathways = Object.entries(pathwayCounts).sort((a, b) => b[1] - a[1]);

  const priorityCount = targets.filter((t) => t.priority === "True").length;

  // Chart data
  const barCategories = pathways.map(([pw]) => pw.replace(/_/g, " "));
  const barSeries = [{ name: "Targets", data: pathways.map(([, c]) => c) }];

  const statusBadge = (val: string | undefined) =>
    val === "True"
      ? "inline-flex items-center rounded-full bg-orange-50 px-2 py-0.5 text-xs font-semibold text-orange-600"
      : "inline-flex items-center rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-400";

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-orange-100 px-2.5 py-1 text-xs font-bold text-orange-600">v3</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Drug Targets</h1>
          <p className="text-sm text-gray-500">
            Parasite-specific druggable targets identified through comparative genomics
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Total Targets"
          value={targets.length}
          subtitle="Top ranked"
          accentColor="bg-orange-500"
        />
        <KpiCard
          title="Priority Targets"
          value={priorityCount}
          subtitle="Highest druggability"
          accentColor="bg-orange-400"
        />
        <KpiCard
          title="Pathways"
          value={pathways.length}
          subtitle="Distinct metabolic pathways"
          accentColor="bg-amber-500"
        />
        <KpiCard
          title="Best Score"
          value={targets[0]?.druggability_score ?? "N/A"}
          subtitle={targets[0]?.gene_name ?? ""}
          accentColor="bg-red-500"
        />
      </div>

      {/* Main content */}
      <div className="grid gap-6 xl:grid-cols-3">
        {/* Table */}
        <div className="xl:col-span-2">
          <div className="rounded-xl bg-white shadow-card overflow-hidden">
            <div className="border-b border-gray-100 px-5 py-4">
              <h2 className="text-sm font-semibold text-gray-900">Top Drug Targets</h2>
              <p className="text-xs text-gray-400 mt-0.5">Ranked by druggability score</p>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm" data-testid="drug-targets-table">
                <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                  <tr>
                    <th className="px-4 py-3">#</th>
                    <th className="px-4 py-3">Gene</th>
                    <th className="px-4 py-3">Pathway</th>
                    <th className="px-4 py-3 text-right">Identity</th>
                    <th className="px-4 py-3 text-right">Druggability</th>
                    <th className="px-4 py-3 text-center">Priority</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-50">
                  {targets.map((t, i) => (
                    <tr key={t.gene_id} className="transition-colors hover:bg-gray-50">
                      <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                      <td className="px-4 py-2.5 font-semibold text-gray-900">{t.gene_name}</td>
                      <td className="px-4 py-2.5 text-xs text-gray-500">
                        {(t.pathway ?? "").replace(/_/g, " ")}
                      </td>
                      <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">
                        {parseFloat(t.identity_score ?? "0").toFixed(2)}
                      </td>
                      <td className="px-4 py-2.5 text-right font-mono text-xs font-semibold text-orange-600">
                        {parseFloat(t.druggability_score ?? "0").toFixed(2)}
                      </td>
                      <td className="px-4 py-2.5 text-center">
                        <span className={statusBadge(t.priority)}>
                          {t.priority === "True" ? "Priority" : "Standard"}
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>

        {/* Right column */}
        <div className="space-y-6">
          {/* Pathway chart */}
          <div className="rounded-xl bg-white shadow-card p-5">
            <h2 className="text-sm font-semibold text-gray-900">Pathway Distribution</h2>
            <p className="text-xs text-gray-400 mt-0.5 mb-4">Targets per metabolic pathway</p>
            <BarChart
              categories={barCategories}
              series={barSeries}
              colors={["#F97316"]}
              height={240}
              horizontal={true}
            />
          </div>

          {/* MRL-003 highlight card */}
          <div className="rounded-xl bg-gradient-to-br from-orange-500 to-amber-500 p-5 text-white shadow-card">
            <div className="mb-3 flex items-center gap-2">
              <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">Lead Candidate</span>
            </div>
            <p className="text-base font-bold">MRL-003</p>
            <p className="mt-1 text-xs text-orange-100">Custom molecule — amide tail variant</p>
            <p className="mt-1 text-xs text-orange-200">Requires selectivity optimization — binds human GR with higher affinity than TryR</p>
            <div className="mt-4 flex items-baseline gap-2">
              <span className="text-3xl font-bold">-7.74</span>
              <span className="text-sm text-orange-200">kcal/mol</span>
            </div>
            <p className="mt-1 text-xs text-orange-100">Trypanothione Reductase (TryR)</p>
            <div className="mt-4 border-t border-white/20 pt-4">
              <div className="grid grid-cols-2 gap-2 text-xs">
                <div>
                  <p className="text-orange-200">Mol. Weight</p>
                  <p className="font-semibold">425.4 Da</p>
                </div>
                <div>
                  <p className="text-orange-200">logP</p>
                  <p className="font-semibold">-0.53</p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Data Sources */}
      <div className="mt-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">Data Sources</h2>
        <p className="mt-2 text-xs text-gray-500 leading-relaxed">
          Targets identified from L. infantum JPCM5 proteome (TriTrypDB). Human orthologs from UniProt.
          Identity Score = % sequence identity to nearest human ortholog via BLASTp (lower = more parasite-specific = safer target).
          Druggability scoring based on target essentiality, structural druggability, and pathway context.
        </p>
      </div>
    </div>
  );
}
