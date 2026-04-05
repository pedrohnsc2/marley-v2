import { loadJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";

interface ConstructEpitope {
  peptide: string;
  gene_id: string;
  gene_name: string;
  allele: string;
  ic50: number;
  rank: number;
  position: number;
}

interface ConstructCard {
  generated: string;
  signal_peptide: string;
  adjuvant: string;
  epitope_count: number;
  protein_length_aa: number;
  mrna_length_nt: number;
  gc_content: number;
  molecular_weight_da: number;
  isoelectric_point: number;
  instability_index: number;
  gravy: number;
  aromaticity: number;
  vaxijen_score: number | null;
  vaxijen_threshold: number;
  restriction_sites_removed: number;
  epitopes: ConstructEpitope[];
}

function extractGeneName(rawName: string): string {
  const match = rawName.match(/gene_product=([^|]+)/);
  if (match) return match[1].trim();
  return rawName.length > 40 ? rawName.slice(0, 40) + "..." : rawName;
}

export default function VaccinePage() {
  const construct = loadJson("construct/construct_card.json") as ConstructCard;

  const marleyScore = 0.3235;
  const leishTecScore = 0.234;
  const maxScore = Math.max(marleyScore, leishTecScore);

  return (
    <div className="px-8 py-10">
      <header className="mb-8">
        <div className="flex items-center gap-2">
          <span className="rounded bg-cyan-900/50 px-2 py-0.5 text-xs font-mono text-cyan-400">
            v1-v2
          </span>
          <h1 className="text-2xl font-bold tracking-tight text-white">
            Vaccine Construct
          </h1>
        </div>
        <p className="mt-2 text-sm text-zinc-400">
          Multi-epitope mRNA vaccine design for canine visceral leishmaniasis
        </p>
      </header>

      <section className="mb-10">
        <div className="grid grid-cols-2 gap-4 lg:grid-cols-4">
          <KpiCard
            title="Epitopes"
            value={construct.epitope_count}
            subtitle="MHC-I binders"
            accentColor="border-cyan-500"
          />
          <KpiCard
            title="Protein Length"
            value={`${construct.protein_length_aa} aa`}
            subtitle={`${construct.mrna_length_nt} nt mRNA`}
            accentColor="border-cyan-400"
          />
          <KpiCard
            title="Mol. Weight"
            value={`${(construct.molecular_weight_da / 1000).toFixed(1)} kDa`}
            subtitle={`pI ${construct.isoelectric_point}`}
            accentColor="border-cyan-500"
          />
          <KpiCard
            title="Instability"
            value={construct.instability_index.toFixed(1)}
            subtitle={construct.instability_index < 40 ? "Stable (< 40)" : "Unstable"}
            accentColor="border-cyan-400"
          />
        </div>
      </section>

      <section className="mb-10">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          Epitope Table
        </h2>
        <div className="overflow-x-auto rounded-lg border border-zinc-800">
          <table
            className="w-full text-left text-sm"
            data-testid="epitope-table"
          >
            <thead className="border-b border-zinc-800 bg-zinc-900/60 text-xs uppercase tracking-wider text-zinc-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">Sequence</th>
                <th className="px-4 py-3">Source Gene</th>
                <th className="px-4 py-3">Allele</th>
                <th className="px-4 py-3 text-right">IC50 (nM)</th>
                <th className="px-4 py-3 text-right">Rank</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-zinc-800">
              {construct.epitopes.map((ep, i) => (
                <tr
                  key={`${ep.peptide}-${ep.gene_id}-${i}`}
                  className="transition-colors hover:bg-zinc-900/40"
                >
                  <td className="px-4 py-2 font-mono text-xs text-zinc-600">
                    {i + 1}
                  </td>
                  <td className="px-4 py-2 font-mono text-sm font-medium text-cyan-300">
                    {ep.peptide}
                  </td>
                  <td className="max-w-xs px-4 py-2 text-xs text-zinc-400">
                    {extractGeneName(ep.gene_name)}
                  </td>
                  <td className="px-4 py-2 text-xs text-zinc-400">
                    {ep.allele}
                  </td>
                  <td className="px-4 py-2 text-right font-mono text-xs text-zinc-300">
                    {ep.ic50.toFixed(2)}
                  </td>
                  <td className="px-4 py-2 text-right font-mono text-xs text-zinc-500">
                    {(ep.rank * 100).toFixed(1)}%
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      <section className="mb-10">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          VaxiJen Antigenicity Comparison
        </h2>
        <div className="rounded-lg border border-zinc-800 bg-zinc-900/40 p-6">
          <div className="space-y-4">
            <div>
              <div className="mb-1 flex items-center justify-between text-sm">
                <span className="font-medium text-white">Marley Vaccine</span>
                <span className="font-mono text-cyan-400">{marleyScore}</span>
              </div>
              <div className="h-6 w-full rounded-full bg-zinc-800">
                <div
                  className="flex h-6 items-center justify-end rounded-full bg-gradient-to-r from-cyan-700 to-cyan-500 pr-2 text-xs font-medium text-white"
                  style={{ width: `${(marleyScore / maxScore) * 100}%` }}
                >
                  {marleyScore}
                </div>
              </div>
            </div>
            <div>
              <div className="mb-1 flex items-center justify-between text-sm">
                <span className="font-medium text-white">Leish-Tec (reference)</span>
                <span className="font-mono text-zinc-400">{leishTecScore}</span>
              </div>
              <div className="h-6 w-full rounded-full bg-zinc-800">
                <div
                  className="flex h-6 items-center justify-end rounded-full bg-gradient-to-r from-zinc-700 to-zinc-500 pr-2 text-xs font-medium text-zinc-300"
                  style={{ width: `${(leishTecScore / maxScore) * 100}%` }}
                >
                  {leishTecScore}
                </div>
              </div>
            </div>
          </div>
          <p className="mt-4 text-xs text-zinc-500">
            VaxiJen score measures predicted antigenicity. Higher values indicate
            stronger predicted immune response. Threshold: {construct.vaxijen_threshold}.
          </p>
        </div>
      </section>

      <section className="mb-10">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          3D Structure Prediction (ESMFold)
        </h2>
        <div className="overflow-hidden rounded-lg border border-zinc-800">
          {/* eslint-disable-next-line @next/next/no-img-element */}
          <img
            src="/images/vaccine_3d.png"
            alt="Predicted 3D structure of the multi-epitope vaccine construct via ESMFold"
            className="w-full max-w-2xl"
          />
        </div>
        <p className="mt-2 text-xs text-zinc-500">
          ESMFold-predicted structure of the {construct.protein_length_aa}-residue
          vaccine construct containing {construct.epitope_count} epitopes linked
          by GPGPG/AAY linkers with {construct.adjuvant} adjuvant and{" "}
          {construct.signal_peptide} signal peptide.
        </p>
      </section>
    </div>
  );
}
