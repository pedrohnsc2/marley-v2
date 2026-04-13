import { getTranslations } from "next-intl/server";
import { loadModuleJson, loadCsv, safeLoadCsv, safeLoadRunModuleJson, safeLoadRunCsv, getRunCompletedDate } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";
import RunBanner from "@/components/runs/run-banner";

interface SlRnaSummary {
  total_transcripts: number;
  sl_positive_count: number;
  sl_gc_content: number;
  sl_entropy: number;
  absent_in_human: boolean;
  sl_sequence: string;
  max_mismatches_allowed: number;
  search_window_nt: number;
}

export default async function RnaPage({
  searchParams,
}: {
  searchParams: Promise<{ run?: string }>;
}) {
  const { run: runId } = await searchParams;
  const t = await getTranslations("rna");

  const summary = (
    runId
      ? safeLoadRunModuleJson(runId, "results", "rna/sl_rna_summary.json") ?? loadModuleJson("results", "rna/sl_rna_summary.json")
      : loadModuleJson("results", "rna/sl_rna_summary.json")
  ) as SlRnaSummary;

  const runEntropy = runId ? safeLoadRunCsv(runId, "rna/shannon_entropy_profile.csv") : [];
  const entropyProfile = runEntropy.length > 0 ? runEntropy : loadCsv("rna/shannon_entropy_profile.csv");

  const runTopTargets = runId ? safeLoadRunCsv(runId, "rna/top_rna_targets.csv") : [];
  const topTargets = runTopTargets.length > 0 ? runTopTargets : loadCsv("rna/top_rna_targets.csv");

  const runCodonUsage = runId ? safeLoadRunCsv(runId, "rna/codon_usage_comparison.csv") : [];
  const codonUsage = runCodonUsage.length > 0 ? runCodonUsage : safeLoadCsv("rna/codon_usage_comparison.csv");

  const runCompletedAt = runId ? getRunCompletedDate(runId) : null;

  // Sort entropy profile by shannon_entropy descending, take top 15
  const topEntropy = [...entropyProfile]
    .sort((a, b) => parseFloat(b.shannon_entropy ?? "0") - parseFloat(a.shannon_entropy ?? "0"))
    .slice(0, 15);

  // Shannon entropy bar chart data
  const entropyCategories = topEntropy.map((r) => r.gene_name ?? r.gene_id ?? "");
  const entropySeries = [
    { name: t("entropyChart.seriesName"), data: topEntropy.map((r) => parseFloat(r.shannon_entropy ?? "0")) },
  ];

  // Entropy comparison chart: parasite vs human for top 10
  const topComparison = [...entropyProfile]
    .sort((a, b) => parseFloat(b.shannon_entropy ?? "0") - parseFloat(a.shannon_entropy ?? "0"))
    .slice(0, 10);

  const comparisonCategories = topComparison.map((r) => r.gene_name ?? r.gene_id ?? "");
  const comparisonSeries = [
    { name: t("comparisonChart.seriesParasite"), data: topComparison.map((r) => parseFloat(r.shannon_entropy ?? "0")) },
    { name: t("comparisonChart.seriesHuman"), data: topComparison.map((r) => parseFloat(r.human_entropy ?? "0")) },
  ];

  // Sort top targets by entropy_delta descending
  const sortedTargets = [...topTargets].sort(
    (a, b) => parseFloat(b.entropy_delta ?? "0") - parseFloat(a.entropy_delta ?? "0"),
  );

  // Codon usage chart: top 20 codons by absolute RSCU difference
  const codonChart =
    codonUsage.length > 0
      ? [...codonUsage]
          .map((row) => ({
            codon: row.codon ?? "",
            amino_acid: row.amino_acid ?? "",
            linf_rscu: parseFloat(row.linf_rscu ?? "0"),
            human_rscu: parseFloat(row.human_rscu ?? "0"),
            delta: Math.abs(parseFloat(row.linf_rscu ?? "0") - parseFloat(row.human_rscu ?? "0")),
          }))
          .sort((a, b) => b.delta - a.delta)
          .slice(0, 20)
      : [];

  const codonCategories = codonChart.map((c) => `${c.codon} (${c.amino_acid})`);
  const codonSeries = [
    { name: t("codonChart.seriesLinf"), data: codonChart.map((c) => c.linf_rscu) },
    { name: t("codonChart.seriesHuman"), data: codonChart.map((c) => c.human_rscu) },
  ];

  const deltaBadge = (val: string | undefined) => {
    const delta = parseFloat(val ?? "0");
    if (delta > 0.3)
      return "inline-flex items-center rounded-full bg-teal-50 px-2 py-0.5 text-xs font-semibold text-teal-700";
    if (delta > 0)
      return "inline-flex items-center rounded-full bg-emerald-50 px-2 py-0.5 text-xs font-semibold text-emerald-600";
    return "inline-flex items-center rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-400";
  };

  return (
    <div>
      <RunBanner runId={runId ?? null} pipeline="rna" completedAt={runCompletedAt} />

      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-teal-100 px-2.5 py-1 text-xs font-bold text-teal-600">{t("badge")}</span>
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
          title={t("kpi.totalTranscripts.title")}
          value={summary.total_transcripts}
          subtitle={t("kpi.totalTranscripts.subtitle")}
          accentColor="bg-teal-500"
        />
        <KpiCard
          title={t("kpi.slConservation.title")}
          value={t("kpi.slConservation.value")}
          subtitle={t("kpi.slConservation.subtitle")}
          accentColor="bg-teal-400"
        />
        <KpiCard
          title={t("kpi.offTargetSafety.title")}
          value={t("kpi.offTargetSafety.value")}
          subtitle={t("kpi.offTargetSafety.subtitle")}
          accentColor="bg-cyan-500"
        />
        <KpiCard
          title={t("kpi.humanAbsent.title")}
          value={summary.absent_in_human ? t("kpi.humanAbsent.yes") : t("kpi.humanAbsent.no")}
          subtitle={t("kpi.humanAbsent.subtitleAbsent")}
          accentColor={summary.absent_in_human ? "bg-emerald-500" : "bg-red-500"}
        />
      </div>

      {/* SL RNA highlight card */}
      <div className="mb-6 rounded-xl bg-gradient-to-br from-teal-500 to-cyan-500 p-5 text-white shadow-card">
        <div className="mb-3 flex items-center gap-2">
          <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">{t("slCard.badgeSL")}</span>
          <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">{t("slCard.badgeTarget")}</span>
        </div>
        <p className="text-base font-bold">{t("slCard.title")}</p>
        <p className="mt-1 text-xs text-teal-100">
          {t("slCard.description")}
        </p>
        <div className="mt-4 rounded-lg bg-white/10 px-4 py-3">
          <p
            className="break-all font-mono text-lg font-bold tracking-widest"
            data-testid="sl-rna-sequence"
          >
            {summary.sl_sequence}
          </p>
        </div>
        <div className="mt-4 grid grid-cols-3 gap-4 border-t border-white/20 pt-4 text-xs">
          <div>
            <p className="text-teal-200">{t("slCard.entropy")}</p>
            <p className="font-semibold">{summary.sl_entropy.toFixed(4)} bits</p>
          </div>
          <div>
            <p className="text-teal-200">{t("slCard.gcContent")}</p>
            <p className="font-semibold">{(summary.sl_gc_content * 100).toFixed(1)}%</p>
          </div>
          <div>
            <p className="text-teal-200">{t("slCard.length")}</p>
            <p className="font-semibold">{summary.sl_sequence.length} nt</p>
          </div>
        </div>
      </div>

      {/* Charts row */}
      <div className="mb-6 grid gap-6 lg:grid-cols-2">
        {/* Shannon entropy chart */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">{t("entropyChart.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            {t("entropyChart.subtitle")}
          </p>
          <BarChart
            categories={entropyCategories}
            series={entropySeries}
            colors={["#14B8A6"]}
            height={280}
          />
        </div>

        {/* Entropy comparison chart */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">{t("comparisonChart.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            {t("comparisonChart.subtitle")}
          </p>
          <BarChart
            categories={comparisonCategories}
            series={comparisonSeries}
            colors={["#14B8A6", "#9CA3AF"]}
            height={280}
          />
        </div>
      </div>

      {/* Top targets table */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">{t("targetsTable.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            {t("targetsTable.subtitle", { count: sortedTargets.length })}
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="rna-targets-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">{t("targetsTable.geneId")}</th>
                <th className="px-4 py-3">{t("targetsTable.geneName")}</th>
                <th className="px-4 py-3 text-right">{t("targetsTable.shannonEntropy")}</th>
                <th className="px-4 py-3 text-right">{t("targetsTable.humanEntropy")}</th>
                <th className="px-4 py-3 text-right">{t("targetsTable.delta")}</th>
                <th className="px-4 py-3 text-center">{t("targetsTable.priority")}</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {sortedTargets.map((tgt, i) => (
                <tr key={`${tgt.gene_id}-${i}`} className="transition-colors hover:bg-gray-50">
                  <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                  <td className="px-4 py-2.5 font-mono text-xs text-gray-600">{tgt.gene_id}</td>
                  <td className="px-4 py-2.5 font-semibold text-gray-900">{tgt.gene_name}</td>
                  <td className="px-4 py-2.5 text-right font-mono text-xs text-teal-600">
                    {parseFloat(tgt.shannon_entropy ?? "0").toFixed(4)}
                  </td>
                  <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">
                    {parseFloat(tgt.human_entropy ?? "0").toFixed(4)}
                  </td>
                  <td className="px-4 py-2.5 text-right">
                    <span className={deltaBadge(tgt.entropy_delta)}>
                      {parseFloat(tgt.entropy_delta ?? "0") > 0 ? "+" : ""}
                      {parseFloat(tgt.entropy_delta ?? "0").toFixed(4)}
                    </span>
                  </td>
                  <td className="px-4 py-2.5 text-center">
                    {tgt.is_priority === "True" ? (
                      <span className="inline-flex items-center rounded-full bg-teal-50 px-2 py-0.5 text-xs font-semibold text-teal-700">
                        {t("targetsTable.priorityLabel")}
                      </span>
                    ) : (
                      <span className="inline-flex items-center rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-400">
                        {t("targetsTable.standardLabel")}
                      </span>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Codon usage section — only rendered if data exists */}
      {codonChart.length > 0 && (
        <div className="mb-6 rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">{t("codonChart.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            {t("codonChart.subtitle")}
          </p>
          <BarChart
            categories={codonCategories}
            series={codonSeries}
            colors={["#14B8A6", "#D1D5DB"]}
            height={320}
          />
        </div>
      )}
    </div>
  );
}
