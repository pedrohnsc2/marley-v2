import { getTranslations } from "next-intl/server";
import { loadCsv, safeLoadRunCsv, getRunCompletedDate } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";
import RunBanner from "@/components/runs/run-banner";

export default async function DrugPage({
  searchParams,
}: {
  searchParams: Promise<{ run?: string }>;
}) {
  const { run: runId } = await searchParams;
  const t = await getTranslations("drug");

  const runTargets = runId ? safeLoadRunCsv(runId, "drug_targets_top20.csv") : [];
  const targets = runTargets.length > 0 ? runTargets : loadCsv("drug_targets_top20.csv");
  const runCompletedAt = runId ? getRunCompletedDate(runId) : null;

  // Count targets by pathway
  const pathwayCounts: Record<string, number> = {};
  for (const tgt of targets) {
    const pw = tgt.pathway ?? "unknown";
    pathwayCounts[pw] = (pathwayCounts[pw] ?? 0) + 1;
  }
  const pathways = Object.entries(pathwayCounts).sort((a, b) => b[1] - a[1]);

  const priorityCount = targets.filter((tgt) => tgt.priority === "True").length;

  // Chart data
  const barCategories = pathways.map(([pw]) => pw.replace(/_/g, " "));
  const barSeries = [{ name: t("pathwayChart.seriesName"), data: pathways.map(([, c]) => c) }];

  const statusBadge = (val: string | undefined) =>
    val === "True"
      ? "inline-flex items-center rounded-full bg-orange-50 px-2 py-0.5 text-xs font-semibold text-orange-600"
      : "inline-flex items-center rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-400";

  return (
    <div>
      <RunBanner runId={runId ?? null} pipeline="drug" completedAt={runCompletedAt} />

      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-orange-100 px-2.5 py-1 text-xs font-bold text-orange-600">{t("badge")}</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">{t("title")}</h1>
          <p className="text-sm text-gray-500">
            {t("subtitle")}
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title={t("kpi.totalTargets.title")}
          value={targets.length}
          subtitle={t("kpi.totalTargets.subtitle")}
          accentColor="bg-orange-500"
        />
        <KpiCard
          title={t("kpi.priorityTargets.title")}
          value={priorityCount}
          subtitle={t("kpi.priorityTargets.subtitle")}
          accentColor="bg-orange-400"
        />
        <KpiCard
          title={t("kpi.pathways.title")}
          value={pathways.length}
          subtitle={t("kpi.pathways.subtitle")}
          accentColor="bg-amber-500"
        />
        <KpiCard
          title={t("kpi.bestScore.title")}
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
              <h2 className="text-sm font-semibold text-gray-900">{t("table.title")}</h2>
              <p className="text-xs text-gray-400 mt-0.5">{t("table.subtitle")}</p>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm" data-testid="drug-targets-table">
                <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                  <tr>
                    <th className="px-4 py-3">#</th>
                    <th className="px-4 py-3">{t("table.gene")}</th>
                    <th className="px-4 py-3">{t("table.pathway")}</th>
                    <th className="px-4 py-3 text-right">{t("table.identity")}</th>
                    <th className="px-4 py-3 text-right">{t("table.druggability")}</th>
                    <th className="px-4 py-3 text-center">{t("table.priority")}</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-50">
                  {targets.map((tgt, i) => (
                    <tr key={tgt.gene_id} className="transition-colors hover:bg-gray-50">
                      <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                      <td className="px-4 py-2.5 font-semibold text-gray-900">{tgt.gene_name}</td>
                      <td className="px-4 py-2.5 text-xs text-gray-500">
                        {(tgt.pathway ?? "").replace(/_/g, " ")}
                      </td>
                      <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">
                        {parseFloat(tgt.identity_score ?? "0").toFixed(2)}
                      </td>
                      <td className="px-4 py-2.5 text-right font-mono text-xs font-semibold text-orange-600">
                        {parseFloat(tgt.druggability_score ?? "0").toFixed(2)}
                      </td>
                      <td className="px-4 py-2.5 text-center">
                        <span className={statusBadge(tgt.priority)}>
                          {tgt.priority === "True" ? t("table.priorityLabel") : t("table.standardLabel")}
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
            <h2 className="text-sm font-semibold text-gray-900">{t("pathwayChart.title")}</h2>
            <p className="text-xs text-gray-400 mt-0.5 mb-4">{t("pathwayChart.subtitle")}</p>
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
              <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">{t("leadCandidate.badge")}</span>
            </div>
            <p className="text-base font-bold">{t("leadCandidate.name")}</p>
            <p className="mt-1 text-xs text-orange-100">{t("leadCandidate.description")}</p>
            <p className="mt-1 text-xs text-orange-200">{t("leadCandidate.warning")}</p>
            <div className="mt-4 flex items-baseline gap-2">
              <span className="text-3xl font-bold">-7.74</span>
              <span className="text-sm text-orange-200">kcal/mol</span>
            </div>
            <p className="mt-1 text-xs text-orange-100">{t("leadCandidate.target")}</p>
            <div className="mt-4 border-t border-white/20 pt-4">
              <div className="grid grid-cols-2 gap-2 text-xs">
                <div>
                  <p className="text-orange-200">{t("leadCandidate.molWeight")}</p>
                  <p className="font-semibold">425.4 Da</p>
                </div>
                <div>
                  <p className="text-orange-200">{t("leadCandidate.logP")}</p>
                  <p className="font-semibold">-0.53</p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Data Sources */}
      <div className="mt-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">{t("dataSources.title")}</h2>
        <p className="mt-2 text-xs text-gray-500 leading-relaxed">
          {t("dataSources.text")}
        </p>
      </div>
    </div>
  );
}
