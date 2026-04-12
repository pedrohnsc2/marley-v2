import { loadModuleJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

/* ------------------------------------------------------------------ */
/* Types for the comparison_results.json payload                      */
/* ------------------------------------------------------------------ */

interface MatrixEntry<T = number> {
  label: string;
  platform_a: T;
  platform_b: T;
  platform_c: T;
  note?: string;
  platform_c_protein?: number;
}

interface RadarPlatform {
  name: string;
  values: number[];
}

interface ComparisonResults {
  metadata: {
    generated: string;
    pipeline_version: string;
    module: string;
    platforms_compared: string[];
  };
  comparison_matrix: {
    construct_properties: {
      protein_length_aa: MatrixEntry;
      molecular_weight_kda: MatrixEntry;
      instability_index: MatrixEntry;
      is_stable: MatrixEntry<boolean>;
      epitope_count: MatrixEntry;
      antigenicity_proxy: MatrixEntry;
    };
    production: {
      expression_system: MatrixEntry<string>;
      purification: MatrixEntry<string>;
      adjuvant: MatrixEntry<string>;
      post_translational_mods: MatrixEntry<string>;
    };
    logistics: {
      cold_chain: MatrixEntry<string>;
      infrastructure: MatrixEntry<string>;
    };
    cost: {
      cost_per_dose_lab_usd: MatrixEntry;
      cost_per_dose_industrial_usd: MatrixEntry;
      cost_per_animal_lab_usd: MatrixEntry;
      cost_per_animal_industrial_usd: MatrixEntry;
    };
    compatibility: {
      aso_compatible: MatrixEntry<boolean>;
      aso_note: MatrixEntry<string>;
    };
    strategic: {
      funding_score: MatrixEntry;
      funding_note: MatrixEntry<string>;
    };
  };
  radar_chart: {
    axes: string[];
    axis_labels: string[];
    platforms: {
      platform_a: RadarPlatform;
      platform_b: RadarPlatform;
      platform_c: RadarPlatform;
    };
  };
}

/* ------------------------------------------------------------------ */
/* Helpers                                                            */
/* ------------------------------------------------------------------ */

function fmtUsd(v: number): string {
  if (v >= 1) return `$${v.toFixed(2)}`;
  if (v >= 0.01) return `$${v.toFixed(4)}`;
  return `$${v.toExponential(2)}`;
}

const PLATFORM_LABELS: Record<string, string> = {
  platform_a: "mRNA-LNP",
  platform_b: "E. coli",
  platform_c: "L. tarentolae",
};

const PLATFORM_COLORS = {
  a: "#3B82F6",
  b: "#10B981",
  c: "#F59E0B",
};

/* ------------------------------------------------------------------ */
/* Page component (Server Component)                                  */
/* ------------------------------------------------------------------ */

export default function PlatformsPage() {
  const data = loadModuleJson("platforms", "comparison_results.json") as ComparisonResults;
  const cm = data.comparison_matrix;
  const radar = data.radar_chart;

  /* ---- Radar-as-grouped-bars chart data ---- */
  const radarCategories = radar.axis_labels;
  const radarSeries = [
    { name: radar.platforms.platform_a.name, data: radar.platforms.platform_a.values },
    { name: radar.platforms.platform_b.name, data: radar.platforms.platform_b.values },
    { name: radar.platforms.platform_c.name, data: radar.platforms.platform_c.values },
  ];

  /* ---- Cost chart data (lab + industrial per platform) ---- */
  const costCategories = [PLATFORM_LABELS.platform_a, PLATFORM_LABELS.platform_b, PLATFORM_LABELS.platform_c];
  const costSeries = [
    {
      name: "Lab Scale (USD)",
      data: [
        cm.cost.cost_per_dose_lab_usd.platform_a,
        cm.cost.cost_per_dose_lab_usd.platform_b,
        cm.cost.cost_per_dose_lab_usd.platform_c,
      ],
    },
    {
      name: "Industrial Scale (USD)",
      data: [
        cm.cost.cost_per_dose_industrial_usd.platform_a,
        cm.cost.cost_per_dose_industrial_usd.platform_b,
        cm.cost.cost_per_dose_industrial_usd.platform_c,
      ],
    },
  ];

  /* ---- Comparison table rows ---- */
  const comparisonRows: { property: string; a: string; b: string; c: string }[] = [
    {
      property: "Protein Length (aa)",
      a: String(cm.construct_properties.protein_length_aa.platform_a),
      b: String(cm.construct_properties.protein_length_aa.platform_b),
      c: String(cm.construct_properties.protein_length_aa.platform_c),
    },
    {
      property: "Molecular Weight (kDa)",
      a: String(cm.construct_properties.molecular_weight_kda.platform_a),
      b: String(cm.construct_properties.molecular_weight_kda.platform_b),
      c: String(cm.construct_properties.molecular_weight_kda.platform_c),
    },
    {
      property: "Instability Index",
      a: cm.construct_properties.instability_index.platform_a.toFixed(2),
      b: cm.construct_properties.instability_index.platform_b.toFixed(2),
      c: cm.construct_properties.instability_index.platform_c.toFixed(2),
    },
    {
      property: "Epitopes",
      a: String(cm.construct_properties.epitope_count.platform_a),
      b: String(cm.construct_properties.epitope_count.platform_b),
      c: String(cm.construct_properties.epitope_count.platform_c),
    },
    {
      property: "Antigenicity Proxy",
      a: cm.construct_properties.antigenicity_proxy.platform_a.toFixed(4),
      b: cm.construct_properties.antigenicity_proxy.platform_b.toFixed(4),
      c: cm.construct_properties.antigenicity_proxy.platform_c.toFixed(4),
    },
    {
      property: "Expression System",
      a: cm.production.expression_system.platform_a,
      b: cm.production.expression_system.platform_b,
      c: cm.production.expression_system.platform_c,
    },
    {
      property: "Purification",
      a: cm.production.purification.platform_a,
      b: cm.production.purification.platform_b,
      c: cm.production.purification.platform_c,
    },
    {
      property: "Adjuvant",
      a: cm.production.adjuvant.platform_a,
      b: cm.production.adjuvant.platform_b,
      c: cm.production.adjuvant.platform_c,
    },
    {
      property: "Cold Chain",
      a: cm.logistics.cold_chain.platform_a,
      b: cm.logistics.cold_chain.platform_b,
      c: cm.logistics.cold_chain.platform_c,
    },
    {
      property: "Infrastructure",
      a: cm.logistics.infrastructure.platform_a,
      b: cm.logistics.infrastructure.platform_b,
      c: cm.logistics.infrastructure.platform_c,
    },
  ];

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-emerald-100 px-2.5 py-1 text-xs font-bold text-emerald-600">
          v{data.metadata.pipeline_version}
        </span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Vaccine Platforms</h1>
          <p className="text-sm text-gray-500">
            Cross-platform comparison of 3 production systems for the Marley multi-epitope vaccine
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Platforms Compared"
          value={data.metadata.platforms_compared.length}
          subtitle="mRNA-LNP, E. coli, L. tarentolae"
          accentColor="bg-emerald-500"
        />
        <KpiCard
          title="Best Cost / Dose"
          value={fmtUsd(cm.cost.cost_per_dose_industrial_usd.platform_c)}
          subtitle="Platform C (live, industrial)"
          accentColor="bg-emerald-500"
        />
        <KpiCard
          title="Most Stable"
          value={`II = ${cm.construct_properties.instability_index.platform_c.toFixed(2)}`}
          subtitle="Platform C (L. tarentolae)"
          accentColor="bg-emerald-500"
        />
        <KpiCard
          title="Best Funding Score"
          value={cm.strategic.funding_score.platform_a.toFixed(1)}
          subtitle="Platform A (mRNA-LNP)"
          accentColor="bg-emerald-500"
        />
      </div>

      {/* Radar-as-grouped-bars chart */}
      <div className="mb-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">Multi-Axis Platform Comparison</h2>
        <p className="text-xs text-gray-400 mt-0.5 mb-4">
          Normalized scores (0-1) across 6 strategic axes for each production platform
        </p>
        <BarChart
          categories={radarCategories}
          series={radarSeries}
          colors={[PLATFORM_COLORS.a, PLATFORM_COLORS.b, PLATFORM_COLORS.c]}
          height={340}
        />
      </div>

      {/* Comparison table */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Detailed Platform Comparison</h2>
          <p className="text-xs text-gray-400 mt-0.5">
            Construct properties, production parameters, and logistics across all 3 platforms
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="platforms-comparison-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-5 py-3">Property</th>
                <th className="px-5 py-3">
                  <span className="inline-block h-2 w-2 rounded-full mr-1.5" style={{ backgroundColor: PLATFORM_COLORS.a }} />
                  Platform A (mRNA-LNP)
                </th>
                <th className="px-5 py-3">
                  <span className="inline-block h-2 w-2 rounded-full mr-1.5" style={{ backgroundColor: PLATFORM_COLORS.b }} />
                  Platform B (E. coli)
                </th>
                <th className="px-5 py-3">
                  <span className="inline-block h-2 w-2 rounded-full mr-1.5" style={{ backgroundColor: PLATFORM_COLORS.c }} />
                  Platform C (L. tarentolae)
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {comparisonRows.map((row) => (
                <tr key={row.property} className="transition-colors hover:bg-gray-50">
                  <td className="px-5 py-3 font-semibold text-gray-900 whitespace-nowrap">{row.property}</td>
                  <td className="px-5 py-3 text-gray-600 max-w-xs">{row.a}</td>
                  <td className="px-5 py-3 text-gray-600 max-w-xs">{row.b}</td>
                  <td className="px-5 py-3 text-gray-600 max-w-xs">{row.c}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* Cost comparison chart */}
      <div className="mb-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">Cost per Dose Comparison</h2>
        <p className="text-xs text-gray-400 mt-0.5 mb-4">
          Lab-scale vs industrial-scale cost per dose (USD) across production platforms
        </p>
        <BarChart
          categories={costCategories}
          series={costSeries}
          colors={["#6366F1", "#10B981"]}
          height={280}
        />
      </div>

      {/* ASO Compatibility + Strategic Funding — two-column grid */}
      <div className="mb-6 grid gap-6 lg:grid-cols-2">
        {/* ASO Compatibility card */}
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">ASO Compatibility</h2>
            <p className="text-xs text-gray-400 mt-0.5">
              Compatibility with MRL-ASO-001 antisense oligonucleotide therapy
            </p>
          </div>
          <div className="p-5 grid gap-4">
            {(["platform_a", "platform_b", "platform_c"] as const).map((key) => {
              const isCompatible = cm.compatibility.aso_compatible[key];
              const note = cm.compatibility.aso_note[key];
              return (
                <div
                  key={key}
                  className={`rounded-lg border p-4 ${
                    isCompatible
                      ? "border-emerald-200 bg-emerald-50"
                      : "border-red-200 bg-red-50"
                  }`}
                  data-testid={`aso-compat-${key}`}
                >
                  <div className="flex items-center gap-3 mb-2">
                    {isCompatible ? (
                      <span className="flex h-7 w-7 items-center justify-center rounded-full bg-emerald-100 text-emerald-600">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2.5} className="h-4 w-4">
                          <polyline points="20 6 9 17 4 12" />
                        </svg>
                      </span>
                    ) : (
                      <span className="flex h-7 w-7 items-center justify-center rounded-full bg-red-100 text-red-600">
                        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2.5} className="h-4 w-4">
                          <line x1="18" y1="6" x2="6" y2="18" />
                          <line x1="6" y1="6" x2="18" y2="18" />
                        </svg>
                      </span>
                    )}
                    <div>
                      <p className="text-sm font-semibold text-gray-900">{PLATFORM_LABELS[key]}</p>
                      <p className={`text-xs font-semibold ${isCompatible ? "text-emerald-700" : "text-red-600"}`}>
                        {isCompatible ? "Compatible" : "Incompatible"}
                      </p>
                    </div>
                  </div>
                  <p className="text-xs text-gray-600 leading-relaxed">{note}</p>
                </div>
              );
            })}
          </div>
        </div>

        {/* Strategic Funding card */}
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">Strategic Funding Scores</h2>
            <p className="text-xs text-gray-400 mt-0.5">
              Funding feasibility scores (0-10) with strategic rationale
            </p>
          </div>
          <div className="p-5 grid gap-4">
            {(["platform_a", "platform_b", "platform_c"] as const).map((key) => {
              const score = cm.strategic.funding_score[key];
              const note = cm.strategic.funding_note[key];
              const colorKey = key === "platform_a" ? "a" : key === "platform_b" ? "b" : "c";
              const barWidth = (score / 10) * 100;
              return (
                <div key={key} className="rounded-lg bg-gray-50 p-4" data-testid={`funding-${key}`}>
                  <div className="flex items-center justify-between mb-2">
                    <p className="text-sm font-semibold text-gray-900">{PLATFORM_LABELS[key]}</p>
                    <span className="text-lg font-bold text-gray-900">{score.toFixed(1)}</span>
                  </div>
                  <div className="h-2 w-full rounded-full bg-gray-200 mb-2">
                    <div
                      className="h-2 rounded-full transition-all"
                      style={{
                        width: `${barWidth}%`,
                        backgroundColor: PLATFORM_COLORS[colorKey],
                      }}
                    />
                  </div>
                  <p className="text-xs text-gray-500 leading-relaxed">{note}</p>
                </div>
              );
            })}
          </div>
        </div>
      </div>
    </div>
  );
}
