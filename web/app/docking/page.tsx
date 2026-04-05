import { loadCsv } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";

export default function DockingPage() {
  const allScores = loadCsv("docking_scores.csv");
  const top10 = allScores.slice(0, 10);

  const bestAffinity = allScores.length > 0
    ? Math.min(...allScores.map((s) => parseFloat(s.binding_affinity ?? "0")))
    : 0;

  const uniqueTargets = new Set(allScores.map((s) => s.target_gene_name)).size;
  const uniqueCompounds = new Set(allScores.map((s) => s.compound_id)).size;

  // Selectivity data (MRL-003 vs human GR)
  const selectivity = [
    {
      compound: "MRL-003",
      tryr: -7.32,
      humanGr: -8.68,
      delta: -1.35,
      selective: false,
    },
    {
      compound: "Pemetrexed",
      tryr: -6.92,
      humanGr: -9.28,
      delta: -2.36,
      selective: false,
    },
    {
      compound: "Menadione",
      tryr: -4.99,
      humanGr: -6.5,
      delta: -1.51,
      selective: false,
    },
  ];

  return (
    <div className="px-8 py-10">
      <header className="mb-8">
        <div className="flex items-center gap-2">
          <span className="rounded bg-green-900/50 px-2 py-0.5 text-xs font-mono text-green-400">
            v4
          </span>
          <h1 className="text-2xl font-bold tracking-tight text-white">
            Molecular Docking
          </h1>
        </div>
        <p className="mt-2 text-sm text-zinc-400">
          AutoDock Vina virtual screening of repurposed and custom compounds
          against parasite drug targets
        </p>
      </header>

      <section className="mb-10">
        <div className="grid grid-cols-2 gap-4 lg:grid-cols-4">
          <KpiCard
            title="Total Dockings"
            value={allScores.length}
            subtitle="Target-compound pairs"
            accentColor="border-green-500"
          />
          <KpiCard
            title="Unique Targets"
            value={uniqueTargets}
            subtitle="Receptor structures"
            accentColor="border-green-400"
          />
          <KpiCard
            title="Compounds"
            value={uniqueCompounds}
            subtitle="Tested molecules"
            accentColor="border-green-500"
          />
          <KpiCard
            title="Best Affinity"
            value={`${bestAffinity.toFixed(2)} kcal/mol`}
            subtitle={top10[0]?.compound_name ?? ""}
            accentColor="border-green-400"
          />
        </div>
      </section>

      <section className="mb-10">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          Top 10 Docking Hits
        </h2>
        <div className="overflow-x-auto rounded-lg border border-zinc-800">
          <table
            className="w-full text-left text-sm"
            data-testid="docking-table"
          >
            <thead className="border-b border-zinc-800 bg-zinc-900/60 text-xs uppercase tracking-wider text-zinc-400">
              <tr>
                <th className="px-4 py-3">Rank</th>
                <th className="px-4 py-3">Target</th>
                <th className="px-4 py-3">Compound</th>
                <th className="px-4 py-3 text-right">Affinity (kcal/mol)</th>
                <th className="px-4 py-3 text-right">RMSD LB</th>
                <th className="px-4 py-3 text-right">RMSD UB</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-zinc-800">
              {top10.map((row, i) => (
                <tr
                  key={`${row.target_gene_id}-${row.compound_id}-${i}`}
                  className="transition-colors hover:bg-zinc-900/40"
                >
                  <td className="px-4 py-2 font-mono text-xs text-zinc-600">
                    {i + 1}
                  </td>
                  <td className="px-4 py-2 font-semibold text-white">
                    {row.target_gene_name}
                  </td>
                  <td className="px-4 py-2 font-mono text-xs text-green-300">
                    {row.compound_name}
                  </td>
                  <td className="px-4 py-2 text-right font-mono text-sm text-zinc-200">
                    {parseFloat(row.binding_affinity ?? "0").toFixed(3)}
                  </td>
                  <td className="px-4 py-2 text-right font-mono text-xs text-zinc-500">
                    {parseFloat(row.rmsd_lb ?? "0").toFixed(1)}
                  </td>
                  <td className="px-4 py-2 text-right font-mono text-xs text-zinc-500">
                    {parseFloat(row.rmsd_ub ?? "0").toFixed(1)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      <div className="grid gap-8 xl:grid-cols-2">
        <section>
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            Selectivity Analysis: Parasite TryR vs Human GR
          </h2>
          <div className="rounded-lg border border-zinc-800 bg-zinc-900/40 p-5">
            <p className="mb-4 text-xs text-zinc-400">
              Comparison of binding affinities against L. infantum TryR versus
              human glutathione reductase (GR). A selective compound should bind
              more tightly to the parasite target (threshold: 1.5 kcal/mol
              difference).
            </p>
            <div className="overflow-x-auto">
              <table className="w-full text-left text-sm">
                <thead className="border-b border-zinc-800 text-xs uppercase text-zinc-500">
                  <tr>
                    <th className="pb-2 pr-4">Compound</th>
                    <th className="pb-2 pr-4 text-right">TryR</th>
                    <th className="pb-2 pr-4 text-right">Human GR</th>
                    <th className="pb-2 pr-4 text-right">Delta</th>
                    <th className="pb-2 text-center">Selective</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-zinc-800/50">
                  {selectivity.map((s) => (
                    <tr key={s.compound}>
                      <td className="py-2 pr-4 font-medium text-white">
                        {s.compound}
                      </td>
                      <td className="py-2 pr-4 text-right font-mono text-xs text-green-300">
                        {s.tryr}
                      </td>
                      <td className="py-2 pr-4 text-right font-mono text-xs text-zinc-400">
                        {s.humanGr}
                      </td>
                      <td className="py-2 pr-4 text-right font-mono text-xs text-yellow-400">
                        {s.delta}
                      </td>
                      <td className="py-2 text-center">
                        {s.selective ? (
                          <span className="text-xs text-green-400">Yes</span>
                        ) : (
                          <span className="text-xs text-red-400">No</span>
                        )}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <p className="mt-4 text-xs text-zinc-600">
              Note: None of the tested compounds showed parasite selectivity.
              This highlights the need for further structural optimization.
            </p>
          </div>
        </section>

        <section>
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            Docking Pose: GMPS + Methotrexate
          </h2>
          <div className="overflow-hidden rounded-lg border border-zinc-800">
            {/* eslint-disable-next-line @next/next/no-img-element */}
            <img
              src="/images/docking_3d.png"
              alt="Docking pose of methotrexate in the GMPS active site"
              className="w-full"
            />
          </div>
          <p className="mt-2 text-xs text-zinc-500">
            Best docking pose of the top hit in the GMPS binding pocket.
            Rendered with PyMOL.
          </p>

          <h2 className="mb-4 mt-8 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            TryR + MRL-003 Docking
          </h2>
          <div className="overflow-hidden rounded-lg border border-zinc-800">
            {/* eslint-disable-next-line @next/next/no-img-element */}
            <img
              src="/images/tryr_mrl003.png"
              alt="MRL-003 docked into trypanothione reductase active site"
              className="w-full"
            />
          </div>
          <p className="mt-2 text-xs text-zinc-500">
            MRL-003 custom molecule (-7.74 kcal/mol) in the TryR binding site.
          </p>
        </section>
      </div>
    </div>
  );
}
