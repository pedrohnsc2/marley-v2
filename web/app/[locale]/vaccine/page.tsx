import dynamic from "next/dynamic";
import { getTranslations } from "next-intl/server";
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

export default async function VaccinePage() {
  const t = await getTranslations("vaccine");
  const construct = loadJson("construct/construct_card.json") as ConstructCard;
  const vaccinePdb = safeLoadPdb("data/structures/marley_vaccine_construct.pdb");

  // Deduplicate epitopes by peptide sequence
  const uniqueEpitopes = [...new Map(construct.epitopes.map(e => [e.peptide, e])).values()];

  // IC50 distribution chart (unique epitopes only)
  const ic50Categories = uniqueEpitopes.map((e) => e.peptide);
  const ic50Series = [{ name: t("ic50Chart.seriesName"), data: uniqueEpitopes.map((e) => e.ic50) }];

  const isStable = construct.instability_index < 40;

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-blue-100 px-2.5 py-1 text-xs font-bold text-blue-600">{t("badge")}</span>
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
          title={t("kpi.epitopes.title")}
          value={uniqueEpitopes.length}
          subtitle={t("kpi.epitopes.subtitle", { total: construct.epitope_count })}
          accentColor="bg-blue-500"
        />
        <KpiCard
          title={t("kpi.proteinLength.title")}
          value={`${construct.protein_length_aa} aa`}
          subtitle={t("kpi.proteinLength.subtitle", { length: construct.mrna_length_nt })}
          accentColor="bg-cyan-500"
        />
        <KpiCard
          title={t("kpi.molWeight.title")}
          value={`${(construct.molecular_weight_da / 1000).toFixed(1)} kDa`}
          subtitle={t("kpi.molWeight.subtitle", { pi: construct.isoelectric_point, gc: construct.gc_content.toFixed(1) })}
          accentColor="bg-indigo-500"
        />
        <KpiCard
          title={t("kpi.instabilityIndex.title")}
          value={construct.instability_index.toFixed(1)}
          subtitle={isStable ? t("kpi.instabilityIndex.subtitleStable") : t("kpi.instabilityIndex.subtitleUnstable")}
          accentColor={isStable ? "bg-emerald-500" : "bg-red-500"}
        />
      </div>

      {/* Charts row */}
      <div className="mb-6 grid gap-6 lg:grid-cols-2">
        {/* Antigenicity Assessment */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">{t("antigenicity.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-3">{t("antigenicity.source")}</p>
          <div className="rounded-lg bg-amber-50 p-4">
            <p className="text-sm text-amber-800 font-medium">{t("antigenicity.notComputed")}</p>
            <p className="text-xs text-amber-600 mt-1">{t("antigenicity.disclaimer")}</p>
          </div>
        </div>

        {/* IC50 distribution */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">{t("ic50Chart.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            {t("ic50Chart.subtitle")}
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
          <h2 className="text-sm font-semibold text-gray-900">{t("epitopeTable.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            {t("epitopeTable.subtitle", { count: uniqueEpitopes.length })}
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="epitope-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">{t("epitopeTable.sequence")}</th>
                <th className="px-4 py-3">{t("epitopeTable.sourceGene")}</th>
                <th className="px-4 py-3">{t("epitopeTable.allele")}</th>
                <th className="px-4 py-3 text-right">{t("epitopeTable.ic50")}</th>
                <th className="px-4 py-3 text-right">{t("epitopeTable.rank")}</th>
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
          <h2 className="text-sm font-semibold text-gray-900">{t("structure3d.title")}</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            {t("structure3d.subtitle", { length: construct.protein_length_aa, count: uniqueEpitopes.length })}
          </p>
        </div>
        {vaccinePdb ? (
          <MolViewer pdbData={vaccinePdb} preset="protein-plddt" height={450} />
        ) : (
          /* eslint-disable-next-line @next/next/no-img-element */
          <img
            src="/images/vaccine_3d.png"
            alt={t("structure3d.altText")}
            className="w-full max-w-2xl"
          />
        )}
      </div>

      {/* Construct details card */}
      <div className="mb-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="mb-4 text-sm font-semibold text-gray-900">{t("constructProperties.title")}</h2>
        <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-6">
          {[
            { label: t("constructProperties.signalPeptide"), value: construct.signal_peptide },
            { label: t("constructProperties.adjuvant"), value: construct.adjuvant },
            { label: t("constructProperties.gcContent"), value: `${construct.gc_content.toFixed(1)}%` },
            { label: t("constructProperties.gravy"), value: construct.gravy.toFixed(3) },
            { label: t("constructProperties.aromaticity"), value: construct.aromaticity.toFixed(3) },
            { label: t("constructProperties.restrictionSites"), value: construct.restriction_sites_removed },
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
        <h2 className="text-sm font-semibold text-gray-900">{t("methods.title")}</h2>
        <div className="mt-3 grid grid-cols-2 gap-4 sm:grid-cols-3">
          {[
            { label: t("methods.epitopePrediction"), value: t("methods.epitopePredictionValue") },
            { label: t("methods.alleles"), value: t("methods.allelesValue") },
            { label: t("methods.ic50Threshold"), value: t("methods.ic50ThresholdValue") },
            { label: t("methods.adjuvant"), value: t("methods.adjuvantValue") },
            { label: t("methods.structure"), value: t("methods.structureValue") },
            { label: t("methods.codonOptimization"), value: t("methods.codonOptimizationValue") },
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
