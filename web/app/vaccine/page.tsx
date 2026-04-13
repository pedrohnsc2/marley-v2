import dynamic from "next/dynamic";
import { loadJson, safeLoadPdb } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

const MolViewer = dynamic(() => import("@/components/mol-viewer/MolViewer"), { ssr: false });

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
  const vaccinePdb = safeLoadPdb("data/structures/marley_vaccine_construct.pdb");

  // Deduplicate epitopes by peptide sequence
  const uniqueEpitopes = [...new Map(construct.epitopes.map(e => [e.peptide, e])).values()];

  // IC50 distribution chart (unique epitopes only)
  const ic50Categories = uniqueEpitopes.map((e) => e.peptide);
  const ic50Series = [{ name: "IC50 (nM)", data: uniqueEpitopes.map((e) => e.ic50) }];

  const isStable = construct.instability_index < 40;

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-blue-100 px-2.5 py-1 text-xs font-bold text-blue-600">v1-v2</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Vaccine Construct</h1>
          <p className="text-sm text-gray-500">
            Multi-epitope mRNA vaccine design for canine visceral leishmaniasis
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Epitopes"
          value={uniqueEpitopes.length}
          subtitle={`unique MHC-I binders (from ${construct.epitope_count} allele-specific pairs)`}
          accentColor="bg-blue-500"
        />
        <KpiCard
          title="Protein Length"
          value={`${construct.protein_length_aa} aa`}
          subtitle={`${construct.mrna_length_nt} nt mRNA`}
          accentColor="bg-cyan-500"
        />
        <KpiCard
          title="Mol. Weight"
          value={`${(construct.molecular_weight_da / 1000).toFixed(1)} kDa`}
          subtitle={`pI ${construct.isoelectric_point} · GC ${construct.gc_content.toFixed(1)}%`}
          accentColor="bg-indigo-500"
        />
        <KpiCard
          title="Instability Index"
          value={construct.instability_index.toFixed(1)}
          subtitle={isStable ? "Stable (< 40)" : "Unstable (≥ 40)"}
          accentColor={isStable ? "bg-emerald-500" : "bg-red-500"}
        />
      </div>

      {/* Charts row */}
      <div className="mb-6 grid gap-6 lg:grid-cols-2">
        {/* Antigenicity Assessment */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Antigenicity Assessment</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-3">VaxiJen 2.0 (Doytchinova & Flower, 2007)</p>
          <div className="rounded-lg bg-amber-50 p-4">
            <p className="text-sm text-amber-800 font-medium">VaxiJen score not computed for this construct</p>
            <p className="text-xs text-amber-600 mt-1">The VaxiJen server has size limitations for constructs &gt;100 amino acids. Individual epitopes were selected based on predicted immunogenicity via NetMHCpan 4.1. Full-construct antigenicity requires experimental validation (ELISA/ELISPOT).</p>
          </div>
        </div>

        {/* IC50 distribution */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Epitope IC50 Distribution</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            Binding affinity per epitope (nM) — lower is better
          </p>
          <BarChart
            categories={ic50Categories}
            series={ic50Series}
            colors={["#06B6D4"]}
            height={220}
          />
        </div>
      </div>

      {/* Epitope table */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Epitope Table</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            {uniqueEpitopes.length} unique MHC-I binders linked by GPGPG/AAY linkers
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="epitope-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">Sequence</th>
                <th className="px-4 py-3">Source Gene</th>
                <th className="px-4 py-3">Allele</th>
                <th className="px-4 py-3 text-right">IC50 (nM)</th>
                <th className="px-4 py-3 text-right">Rank %</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {uniqueEpitopes.map((ep, i) => (
                <tr
                  key={`${ep.peptide}-${ep.gene_id}-${i}`}
                  className="transition-colors hover:bg-gray-50"
                >
                  <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                  <td className="px-4 py-2.5">
                    <span className="rounded-md bg-blue-50 px-2 py-0.5 font-mono text-xs font-semibold text-blue-700">
                      {ep.peptide}
                    </span>
                  </td>
                  <td className="max-w-xs px-4 py-2.5 text-xs text-gray-500">
                    {extractGeneName(ep.gene_name)}
                  </td>
                  <td className="px-4 py-2.5 text-xs text-gray-500">{ep.allele}</td>
                  <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-700">
                    {ep.ic50.toFixed(2)}
                  </td>
                  <td className="px-4 py-2.5 text-right font-mono text-xs text-gray-400">
                    {(ep.rank * 100).toFixed(1)}%
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* 3D Structure */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">3D Structure Prediction (ESMFold)</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            ESMFold-predicted structure of the {construct.protein_length_aa}-residue construct with {uniqueEpitopes.length} unique epitopes
          </p>
        </div>
        {vaccinePdb ? (
          <MolViewer pdbData={vaccinePdb} preset="protein-plddt" height={450} />
        ) : (
          /* eslint-disable-next-line @next/next/no-img-element */
          <img
            src="/images/vaccine_3d.png"
            alt="Predicted 3D structure of the multi-epitope vaccine construct via ESMFold"
            className="w-full max-w-2xl"
          />
        )}
      </div>

      {/* Construct details card */}
      <div className="mb-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="mb-4 text-sm font-semibold text-gray-900">Construct Properties</h2>
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-6">
          {[
            { label: "Signal Peptide", value: construct.signal_peptide },
            { label: "Adjuvant", value: construct.adjuvant },
            { label: "GC Content", value: `${construct.gc_content.toFixed(1)}%` },
            { label: "GRAVY", value: construct.gravy.toFixed(3) },
            { label: "Aromaticity", value: construct.aromaticity.toFixed(3) },
            { label: "Restriction Sites Removed", value: construct.restriction_sites_removed },
          ].map((prop) => (
            <div key={prop.label} className="rounded-lg bg-gray-50 p-3">
              <p className="text-xs text-gray-400">{prop.label}</p>
              <p className="mt-1 text-sm font-semibold text-gray-900">{prop.value}</p>
            </div>
          ))}
        </div>
      </div>

      {/* Methods card */}
      <div className="rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">Methods</h2>
        <div className="mt-3 grid grid-cols-2 gap-4 sm:grid-cols-3">
          {[
            { label: "Epitope Prediction", value: "NetMHCpan 4.1" },
            { label: "Alleles", value: "DLA-88*03401, *50101, *50801" },
            { label: "IC50 Threshold", value: "< 500 nM" },
            { label: "Adjuvant", value: "L7/L12 ribosomal" },
            { label: "Structure", value: "ESMFold v1" },
            { label: "Codon Optimization", value: "L. infantum codon table" },
          ].map((m) => (
            <div key={m.label} className="rounded-lg bg-gray-50 p-3">
              <p className="text-xs text-gray-400">{m.label}</p>
              <p className="mt-1 text-sm font-semibold text-gray-900">{m.value}</p>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
