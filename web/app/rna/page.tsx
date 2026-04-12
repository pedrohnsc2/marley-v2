import { safeLoadModuleJson, safeLoadCsv } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

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

export default function RnaPage() {
  const summary = safeLoadModuleJson("results", "rna/sl_rna_summary.json") as SlRnaSummary | null;
  const entropyProfile = safeLoadCsv("rna/shannon_entropy_profile.csv");
  const topTargets = safeLoadCsv("rna/top_rna_targets.csv");
  const codonUsage = safeLoadCsv("rna/codon_usage_comparison.csv");

  // Sort entropy profile by shannon_entropy descending, take top 15
  const topEntropy = [...entropyProfile]
    .sort((a, b) => parseFloat(b.shannon_entropy ?? "0") - parseFloat(a.shannon_entropy ?? "0"))
    .slice(0, 15);

  // Shannon entropy bar chart data
  const entropyCategories = topEntropy.map((r) => r.gene_name ?? r.gene_id ?? "");
  const entropySeries = [
    { name: "Shannon Entropy", data: topEntropy.map((r) => parseFloat(r.shannon_entropy ?? "0")) },
  ];

  // Entropy comparison chart: parasite vs human for top 10
  const topComparison = [...entropyProfile]
    .sort((a, b) => parseFloat(b.shannon_entropy ?? "0") - parseFloat(a.shannon_entropy ?? "0"))
    .slice(0, 10);

  const comparisonCategories = topComparison.map((r) => r.gene_name ?? r.gene_id ?? "");
  const comparisonSeries = [
    { name: "Parasite", data: topComparison.map((r) => parseFloat(r.shannon_entropy ?? "0")) },
    { name: "Human", data: topComparison.map((r) => parseFloat(r.human_entropy ?? "0")) },
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
    { name: "L. infantum RSCU", data: codonChart.map((c) => c.linf_rscu) },
    { name: "Human RSCU", data: codonChart.map((c) => c.human_rscu) },
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
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-teal-100 px-2.5 py-1 text-xs font-bold text-teal-600">v1</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">RNA Entropy</h1>
          <p className="text-sm text-gray-500">
            Information theory analysis of <em>L. infantum</em> transcriptome and SL RNA targets
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Total Transcripts"
          value={summary?.total_transcripts ?? "—"}
          subtitle="L. infantum transcriptome"
          accentColor="bg-teal-500"
        />
        <KpiCard
          title="SL RNA Entropy"
          value={summary ? `${summary.sl_entropy.toFixed(3)} bits` : "—"}
          subtitle="Shannon entropy of SL sequence"
          accentColor="bg-teal-400"
        />
        <KpiCard
          title="SL GC Content"
          value={summary ? `${(summary.sl_gc_content * 100).toFixed(1)}%` : "—"}
          subtitle="Low GC = high AU content"
          accentColor="bg-cyan-500"
        />
        <KpiCard
          title="Human Absent"
          value={summary ? (summary.absent_in_human ? "✓ Yes" : "✗ No") : "—"}
          subtitle="Not found in human transcriptome"
          accentColor={summary?.absent_in_human ? "bg-emerald-500" : "bg-teal-500"}
        />
      </div>

      {/* SL RNA highlight card */}
      {summary ? (
        <div className="mb-6 rounded-xl bg-gradient-to-br from-teal-500 to-cyan-500 p-5 text-white shadow-card">
          <div className="mb-3 flex items-center gap-2">
            <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">SL RNA</span>
            <span className="rounded-md bg-white/20 px-2 py-0.5 text-xs font-semibold">MRL-ASO-001 Target</span>
          </div>
          <p className="text-base font-bold">Spliced Leader RNA Sequence</p>
          <p className="mt-1 text-xs text-teal-100">
            39-nt conserved sequence, absent in human transcriptome — primary target for antisense oligonucleotide therapy
          </p>
          <div className="mt-4 rounded-lg bg-white/10 px-4 py-3">
            <p className="break-all font-mono text-lg font-bold tracking-widest" data-testid="sl-rna-sequence">
              {summary.sl_sequence}
            </p>
          </div>
          <div className="mt-4 grid grid-cols-3 gap-4 border-t border-white/20 pt-4 text-xs">
            <div>
              <p className="text-teal-200">Entropy</p>
              <p className="font-semibold">{summary.sl_entropy.toFixed(4)} bits</p>
            </div>
            <div>
              <p className="text-teal-200">GC Content</p>
              <p className="font-semibold">{(summary.sl_gc_content * 100).toFixed(1)}%</p>
            </div>
            <div>
              <p className="text-teal-200">Length</p>
              <p className="font-semibold">{summary.sl_sequence.length} nt</p>
            </div>
          </div>
        </div>
      ) : (
        <div className="mb-6 rounded-xl border border-dashed border-gray-200 bg-gray-50 p-8 text-center">
          <div className="mx-auto mb-3 flex h-10 w-10 items-center justify-center rounded-full bg-teal-100">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.5} className="h-5 w-5 text-teal-500">
              <path strokeLinecap="round" strokeLinejoin="round" d="M12 6v6l4 2m6-2a10 10 0 11-20 0 10 10 0 0120 0z" />
            </svg>
          </div>
          <p className="text-sm font-semibold text-gray-600">SL RNA Sequence</p>
          <p className="mt-1 text-xs text-gray-400">Computação em andamento — resultados disponíveis em breve</p>
        </div>
      )}

      {/* Charts row */}
      {topEntropy.length > 0 && (
        <div className="mb-6 grid gap-6 lg:grid-cols-2">
          <div className="rounded-xl bg-white shadow-card p-5">
            <h2 className="text-sm font-semibold text-gray-900">Shannon Entropy Profile</h2>
            <p className="text-xs text-gray-400 mt-0.5 mb-4">
              Top 15 genes by Shannon entropy (bits) — higher entropy indicates more sequence complexity
            </p>
            <BarChart categories={entropyCategories} series={entropySeries} colors={["#14B8A6"]} height={280} />
          </div>
          <div className="rounded-xl bg-white shadow-card p-5">
            <h2 className="text-sm font-semibold text-gray-900">Parasite vs Human Entropy</h2>
            <p className="text-xs text-gray-400 mt-0.5 mb-4">
              Comparing Shannon entropy between <em>L. infantum</em> and human orthologs
            </p>
            <BarChart categories={comparisonCategories} series={comparisonSeries} colors={["#14B8A6", "#9CA3AF"]} height={280} />
          </div>
        </div>
      )}

      {/* Top targets table */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Top RNA Targets</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            {sortedTargets.length > 0 ? `${sortedTargets.length} ranked targets sorted by entropy delta` : "Computação em andamento"}
          </p>
        </div>
        {sortedTargets.length > 0 ? (
          <div className="overflow-x-auto">
            <table className="w-full text-left text-sm" data-testid="rna-targets-table">
              <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                <tr>
                  <th className="px-4 py-3">#</th>
                  <th className="px-4 py-3">Gene ID</th>
                  <th className="px-4 py-3">Gene Name</th>
                  <th className="px-4 py-3 text-right">Shannon Entropy</th>
                  <th className="px-4 py-3 text-right">Human Entropy</th>
                  <th className="px-4 py-3 text-right">Delta</th>
                  <th className="px-4 py-3 text-center">Priority</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-50">
                {sortedTargets.map((t, i) => (
                  <tr key={`${t.gene_id}-${i}`} className="transition-colors hover:bg-gray-50">
                    <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                    <td className="px-4 py-2.5 font-mono text-xs text-gray-600">{t.gene_id}</td>
                    <td className="px-4 py-2.5 font-semibold text-gray-900">{t.gene_name}</td>
                    <td className="px-4 py-2.5 text-right font-mono text-xs text-teal-600">
                      {parseFloat(t.shannon_entropy ?? "0").toFixed(4)}
                    </td>
                    <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-600">
                      {parseFloat(t.human_entropy ?? "0").toFixed(4)}
                    </td>
                    <td className="px-4 py-2.5 text-right">
                      <span className={deltaBadge(t.entropy_delta)}>
                        {parseFloat(t.entropy_delta ?? "0") > 0 ? "+" : ""}
                        {parseFloat(t.entropy_delta ?? "0").toFixed(4)}
                      </span>
                    </td>
                    <td className="px-4 py-2.5 text-center">
                      {t.is_priority === "True" ? (
                        <span className="inline-flex items-center rounded-full bg-teal-50 px-2 py-0.5 text-xs font-semibold text-teal-700">Priority</span>
                      ) : (
                        <span className="inline-flex items-center rounded-full bg-gray-100 px-2 py-0.5 text-xs text-gray-400">Standard</span>
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        ) : (
          <div className="px-5 py-10 text-center">
            <p className="text-sm text-gray-400">Resultados disponíveis em breve</p>
          </div>
        )}
      </div>

      {/* Codon usage section */}
      {codonChart.length > 0 && (
        <div className="mb-6 rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Codon Usage Comparison</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            Top 20 codons with largest RSCU difference between <em>L. infantum</em> and human
          </p>
          <BarChart categories={codonCategories} series={codonSeries} colors={["#14B8A6", "#D1D5DB"]} height={320} />
        </div>
      )}
    </div>
  );
}
