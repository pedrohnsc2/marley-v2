import dynamic from "next/dynamic";
import { getTranslations } from "next-intl/server";
import { loadCsv, safeLoadPdb, safeLoadRunCsv, safeLoadRunPdb, getRunCompletedDate } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";
import RunBanner from "@/components/runs/run-banner";

const MolViewer = dynamic(() => import("@/components/mol-viewer/MolViewer"), { ssr: false });

export default async function DockingPage({
  searchParams,
}: {
  searchParams: Promise<{ run?: string }>;
}) {
  const { run: runId } = await searchParams;
  const t = await getTranslations("docking");

  const runScores = runId ? safeLoadRunCsv(runId, "docking_scores.csv") : [];
  const allScores = runScores.length > 0 ? runScores : loadCsv("docking_scores.csv");

  // Load 3D structure data for interactive viewers
  const gmpsPdb = runId
    ? safeLoadRunPdb(runId, "data/structures/GMPS_A4IBM8.pdb") ?? safeLoadPdb("data/structures/GMPS_A4IBM8.pdb")
    : safeLoadPdb("data/structures/GMPS_A4IBM8.pdb");
  const gmpsLigand = runId
    ? safeLoadRunPdb(runId, "data/docking/GMPS/CHEMBL36_out.pdbqt") ?? safeLoadPdb("data/docking/GMPS/CHEMBL36_out.pdbqt")
    : safeLoadPdb("data/docking/GMPS/CHEMBL36_out.pdbqt");
  const tryrPdb = runId
    ? safeLoadRunPdb(runId, "data/structures/TryR_Q4Q457.pdb") ?? safeLoadPdb("data/structures/TryR_Q4Q457.pdb")
    : safeLoadPdb("data/structures/TryR_Q4Q457.pdb");
  const mrl003Ligand = runId
    ? safeLoadRunPdb(runId, "data/docking/TryR/MRL-003_amide_tail_out.pdbqt") ?? safeLoadPdb("data/docking/TryR/MRL-003_amide_tail_out.pdbqt")
    : safeLoadPdb("data/docking/TryR/MRL-003_amide_tail_out.pdbqt");
  const top10 = allScores.slice(0, 10);
  const runCompletedAt = runId ? getRunCompletedDate(runId) : null;

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
      name: t("affinityChart.seriesName"),
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
      <RunBanner runId={runId ?? null} pipeline="docking" completedAt={runCompletedAt} />

      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-emerald-100 px-2.5 py-1 text-xs font-bold text-emerald-600">{t("badge")}</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">{t("title")}</h1>
          <p className="text-sm text-gray-500">
            {t("subtitle", { compounds: uniqueCompounds, targets: uniqueTargets })}
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title={t("kpi.totalDockings.title")}
          value={allScores.length}
          subtitle={t("kpi.totalDockings.subtitle")}
          accentColor="bg-emerald-500"
        />
        <KpiCard
          title={t("kpi.uniqueTargets.title")}
          value={uniqueTargets}
          subtitle={t("kpi.uniqueTargets.subtitle")}
          accentColor="bg-teal-500"
        />
        <KpiCard
          title={t("kpi.compounds.title")}
          value={uniqueCompounds}
          subtitle={t("kpi.compounds.subtitle")}
          accentColor="bg-green-500"
        />
        <KpiCard
          title={t("kpi.bestAffinity.title")}
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
            <h2 className="text-sm font-semibold text-gray-900">{t("table.title")}</h2>
            <p className="text-xs text-gray-400 mt-0.5">{t("table.subtitle")}</p>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-left text-sm" data-testid="docking-table">
              <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                <tr>
                  <th className="px-4 py-3">#</th>
                  <th className="px-4 py-3">{t("table.target")}</th>
                  <th className="px-4 py-3">{t("table.compound")}</th>
                  <th className="px-4 py-3 text-right">{t("table.affinity")}</th>
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
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        {/* Affinity chart */}
        <div className="xl:col-span-2 rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">{t("affinityChart.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">{t("affinityChart.subtitle")}</p>
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
          <h2 className="text-sm font-semibold text-gray-900">{t("selectivity.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            {t("selectivity.subtitle")}
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-5 py-3">{t("selectivity.compound")}</th>
                <th className="px-5 py-3 text-right">{t("selectivity.tryr")}</th>
                <th className="px-5 py-3 text-right">{t("selectivity.humanGr")}</th>
                <th className="px-5 py-3 text-right">{t("selectivity.delta")}</th>
                <th className="px-5 py-3 text-center">{t("selectivity.selective")}</th>
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
                      {s.selective ? t("selectivity.yes") : t("selectivity.no")}
                    </span>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <div className="px-5 py-3 bg-amber-50 border-t border-amber-100">
          <p className="text-xs text-amber-700">
            {t("selectivity.warning")}
          </p>
        </div>
      </div>

      {/* Docking pose images */}
      <div className="mt-6 grid gap-6 lg:grid-cols-2">
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">{t("poses.gmps.title")}</h2>
            <p className="text-xs text-gray-400 mt-0.5">{t("poses.gmps.subtitle")}</p>
          </div>
          {gmpsPdb ? (
            <MolViewer
              pdbData={gmpsPdb}
              pdbqtData={gmpsLigand ?? undefined}
              preset="docking"
              height={400}
            />
          ) : (
            /* eslint-disable-next-line @next/next/no-img-element */
            <img src="/images/docking_3d.png" alt={t("poses.gmps.altText")} className="w-full" />
          )}
        </div>
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">{t("poses.tryr.title")}</h2>
            <p className="text-xs text-gray-400 mt-0.5">{t("poses.tryr.subtitle")}</p>
          </div>
          {tryrPdb ? (
            <MolViewer
              pdbData={tryrPdb}
              pdbqtData={mrl003Ligand ?? undefined}
              preset="docking"
              height={400}
            />
          ) : (
            /* eslint-disable-next-line @next/next/no-img-element */
            <img src="/images/tryr_mrl003.png" alt={t("poses.tryr.altText")} className="w-full" />
          )}
        </div>
      </div>
    </div>
  );
}
