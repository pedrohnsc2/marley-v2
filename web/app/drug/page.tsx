import { loadCsv } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";

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

  return (
    <div className="px-8 py-10">
      <header className="mb-8">
        <div className="flex items-center gap-2">
          <span className="rounded bg-orange-900/50 px-2 py-0.5 text-xs font-mono text-orange-400">
            v3
          </span>
          <h1 className="text-2xl font-bold tracking-tight text-white">
            Drug Targets
          </h1>
        </div>
        <p className="mt-2 text-sm text-zinc-400">
          Parasite-specific druggable targets identified through comparative
          genomics and pathway analysis
        </p>
      </header>

      <section className="mb-10">
        <div className="grid grid-cols-2 gap-4 lg:grid-cols-4">
          <KpiCard
            title="Total Targets"
            value={targets.length}
            subtitle="Top 20 ranked"
            accentColor="border-orange-500"
          />
          <KpiCard
            title="Priority Targets"
            value={priorityCount}
            subtitle="Highest druggability"
            accentColor="border-orange-400"
          />
          <KpiCard
            title="Pathways"
            value={pathways.length}
            subtitle="Distinct metabolic pathways"
            accentColor="border-orange-500"
          />
          <KpiCard
            title="Best Score"
            value={targets[0]?.druggability_score ?? "N/A"}
            subtitle={targets[0]?.gene_name ?? ""}
            accentColor="border-orange-400"
          />
        </div>
      </section>

      <div className="grid gap-8 xl:grid-cols-3">
        <section className="xl:col-span-2">
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            Top 20 Drug Targets
          </h2>
          <div className="overflow-x-auto rounded-lg border border-zinc-800">
            <table
              className="w-full text-left text-sm"
              data-testid="drug-targets-table"
            >
              <thead className="border-b border-zinc-800 bg-zinc-900/60 text-xs uppercase tracking-wider text-zinc-400">
                <tr>
                  <th className="px-4 py-3">Rank</th>
                  <th className="px-4 py-3">Gene</th>
                  <th className="px-4 py-3">Pathway</th>
                  <th className="px-4 py-3 text-right">Identity</th>
                  <th className="px-4 py-3 text-right">Druggability</th>
                  <th className="px-4 py-3 text-center">Priority</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-zinc-800">
                {targets.map((t, i) => (
                  <tr
                    key={t.gene_id}
                    className="transition-colors hover:bg-zinc-900/40"
                  >
                    <td className="px-4 py-2 font-mono text-xs text-zinc-600">
                      {i + 1}
                    </td>
                    <td className="px-4 py-2 font-semibold text-white">
                      {t.gene_name}
                    </td>
                    <td className="px-4 py-2 text-xs text-zinc-400">
                      {(t.pathway ?? "").replace(/_/g, " ")}
                    </td>
                    <td className="px-4 py-2 text-right font-mono text-xs text-zinc-300">
                      {parseFloat(t.identity_score ?? "0").toFixed(2)}
                    </td>
                    <td className="px-4 py-2 text-right font-mono text-xs text-orange-300">
                      {parseFloat(t.druggability_score ?? "0").toFixed(2)}
                    </td>
                    <td className="px-4 py-2 text-center">
                      {t.priority === "True" ? (
                        <span className="inline-block rounded-full bg-orange-900/60 px-2 py-0.5 text-xs font-medium text-orange-300">
                          Priority
                        </span>
                      ) : (
                        <span className="text-xs text-zinc-600">--</span>
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>

        <section>
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            Pathway Distribution
          </h2>
          <div className="space-y-3 rounded-lg border border-zinc-800 bg-zinc-900/40 p-5">
            {pathways.map(([pathway, count]) => {
              const pct = (count / targets.length) * 100;
              return (
                <div key={pathway}>
                  <div className="mb-1 flex items-center justify-between text-xs">
                    <span className="text-zinc-300">
                      {pathway.replace(/_/g, " ")}
                    </span>
                    <span className="font-mono text-zinc-500">
                      {count} ({pct.toFixed(0)}%)
                    </span>
                  </div>
                  <div className="h-2 w-full rounded-full bg-zinc-800">
                    <div
                      className="h-2 rounded-full bg-orange-600"
                      style={{ width: `${pct}%` }}
                    />
                  </div>
                </div>
              );
            })}
          </div>

          <h2 className="mb-4 mt-8 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            MRL-003 Custom Molecule
          </h2>
          <div className="rounded-lg border border-orange-800 bg-orange-900/20 p-5">
            <p className="text-sm font-semibold text-white">
              MRL-003 (amide tail variant)
            </p>
            <p className="mt-1 text-xs text-zinc-400">
              Best custom molecule targeting trypanothione reductase (TryR)
            </p>
            <div className="mt-3 flex items-baseline gap-2">
              <span className="text-2xl font-bold text-orange-300">
                -7.74
              </span>
              <span className="text-xs text-zinc-500">kcal/mol</span>
            </div>
            <p className="mt-2 text-xs text-zinc-500">
              Designed via scaffold-hopping from methotrexate. MW 425.4 Da,
              logP -0.53, 1 Lipinski violation.
            </p>
          </div>
        </section>
      </div>
    </div>
  );
}
