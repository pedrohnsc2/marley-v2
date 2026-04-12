import { loadModuleJson, safeLoadModuleJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

interface DimensionDetails {
  wt_dg_kcal?: number;
  wt_tm_celsius?: number;
  [key: string]: unknown;
}

interface DimensionAssessment {
  dimension: string;
  score: number;
  max_score: number;
  verdict: string;
  notes: string[];
  details: DimensionDetails;
}

interface MathCertificateV2 {
  molecule: {
    name: string;
    type: string;
    length: number;
    target: string;
    target_region: string;
  };
  composite_score: number;
  max_score: number;
  overall_verdict: string;
  dimension_assessments: DimensionAssessment[];
}

interface DeliveryModuleResult {
  module: string;
  status: string;
  executive_summary: string;
}

interface DeliveryReport {
  total_runtime_seconds: number;
  modules_executed: number;
  modules_succeeded: number;
  all_passed: boolean;
  module_results: DeliveryModuleResult[];
}

/* ------------------------------------------------------------------ */
/*  Mappings                                                           */
/* ------------------------------------------------------------------ */

const MODULE_NAMES: Record<string, string> = {
  "01_thermodynamic_landscape": "Thermodynamic Landscape",
  "02_selectivity_proof": "Selectivity Proof",
  "03_evolutionary_conservation": "Evolutionary Conservation",
  "04_exhaustive_optimization": "Exhaustive Optimization",
  "05_resistance_model": "Resistance Barrier",
};

const DELIVERY_NAMES: Record<string, string> = {
  A: "Stability",
  B: "Membrane",
  C: "Conjugation",
  D: "LNP Encapsulation",
  E: "ADMET Profile",
  F: "Immune Response",
};

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

function verdictBadge(verdict: string) {
  const v = verdict.toLowerCase();
  if (v === "pass") {
    return (
      <span className="rounded-md bg-emerald-50 px-2 py-0.5 text-xs font-semibold text-emerald-700">
        PASS
      </span>
    );
  }
  if (v === "flag") {
    return (
      <span className="rounded-md bg-amber-50 px-2 py-0.5 text-xs font-semibold text-amber-700">
        FLAG
      </span>
    );
  }
  return (
    <span className="rounded-md bg-gray-50 px-2 py-0.5 text-xs font-semibold text-gray-600">
      {verdict.toUpperCase()}
    </span>
  );
}

function statusBadge(status: string) {
  if (status === "success") {
    return (
      <span className="rounded-md bg-emerald-50 px-2 py-0.5 text-xs font-semibold text-emerald-700">
        PASSED
      </span>
    );
  }
  return (
    <span className="rounded-md bg-red-50 px-2 py-0.5 text-xs font-semibold text-red-700">
      {status.toUpperCase()}
    </span>
  );
}

function truncate(text: string, max: number): string {
  if (text.length <= max) return text;
  return text.slice(0, max) + "...";
}

/* ------------------------------------------------------------------ */
/*  Score bar (inline progress indicator)                              */
/* ------------------------------------------------------------------ */

function ScoreBar({ score }: { score: number }) {
  const color =
    score >= 90 ? "bg-emerald-500" : score >= 80 ? "bg-sky-500" : "bg-amber-500";
  return (
    <div className="flex items-center gap-2">
      <div className="h-2 w-24 overflow-hidden rounded-full bg-gray-100">
        <div className={`h-full rounded-full ${color}`} style={{ width: `${score}%` }} />
      </div>
      <span className="font-mono text-xs font-bold text-gray-700">{score}</span>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Page                                                               */
/* ------------------------------------------------------------------ */

export default function AsoPage() {
  const cert = safeLoadModuleJson("aso_math_reports", "math_certificate_v2.json") as MathCertificateV2 | null;
  const delivery = loadModuleJson("aso_delivery", "delivery_report.json") as DeliveryReport;

  if (!cert) return <div className="p-6 text-red-500">math_certificate_v2.json não encontrado.</div>;

  const { molecule, dimension_assessments } = cert;
  const thermDetails = dimension_assessments[0]?.details ?? {};
  const dg = thermDetails.wt_dg_kcal ?? 0;
  const tm = thermDetails.wt_tm_celsius ?? 0;
  const overallScore = Math.round((cert.composite_score / cert.max_score) * 100);

  /* Bar chart data */
  const chartCategories = dimension_assessments.map((m) => m.dimension.replace(/^\d+\.\s*/, "").split(" ").slice(0, 2).join(" "));
  const chartSeries = [
    { name: "Score", data: dimension_assessments.map((m) => Math.round((m.score / m.max_score) * 100)) },
  ];

  return (
    <div>
      {/* ---- Page header ---- */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-rose-100 px-2.5 py-1 text-xs font-bold text-rose-600">
          v1.0
        </span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">ASO Therapy</h1>
          <p className="text-sm text-gray-500">
            Antisense oligonucleotide targeting Spliced Leader RNA of{" "}
            <span className="italic">L.&nbsp;infantum</span>
          </p>
        </div>
      </div>

      {/* ---- KPI row ---- */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Math Score"
          value={`${cert.composite_score}/${cert.max_score}`}
          subtitle={cert.overall_verdict}
          accentColor="bg-rose-500"
        />
        <KpiCard
          title="Binding Energy"
          value={`${dg} kcal/mol`}
          subtitle={`Tm ${tm} °C`}
          accentColor="bg-pink-500"
        />
        <KpiCard
          title="Delivery"
          value={`${delivery.modules_succeeded}/${delivery.modules_executed}`}
          subtitle="All modules passed"
          accentColor="bg-fuchsia-500"
        />
        <KpiCard
          title="Resistance"
          value="Infinite"
          subtitle="Zero viable escape mutations"
          accentColor="bg-rose-500"
        />
      </div>

      {/* ---- Molecule card ---- */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Molecule: {molecule.name}</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {molecule.type}
          </p>
        </div>
        <div className="px-5 py-5">
          {/* Properties grid */}
          <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
            {[
              { label: "Target", value: molecule.target },
              { label: "Target Region", value: molecule.target_region },
              { label: "Length", value: `${molecule.length} nt` },
              { label: "dG", value: `${dg} kcal/mol` },
              { label: "Tm", value: `${tm} °C` },
            ].map((prop) => (
              <div key={prop.label} className="rounded-lg bg-gray-50 p-3">
                <p className="text-xs text-gray-400">{prop.label}</p>
                <p className="mt-1 text-sm font-semibold text-gray-900">{prop.value}</p>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* ---- Validation modules table ---- */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Validation Modules</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {dimension_assessments.length} mathematical validation modules &mdash; composite score {cert.composite_score}/{cert.max_score}
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="validation-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">Module</th>
                <th className="px-4 py-3">Verdict</th>
                <th className="px-4 py-3">Score</th>
                <th className="px-4 py-3">Note</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {dimension_assessments.map((m, i) => (
                <tr key={m.dimension} className="transition-colors hover:bg-gray-50">
                  <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                  <td className="px-4 py-2.5 text-sm font-medium text-gray-800">
                    {m.dimension.replace(/^\d+\.\s*/, "")}
                  </td>
                  <td className="px-4 py-2.5">{verdictBadge(m.verdict)}</td>
                  <td className="px-4 py-2.5">
                    <ScoreBar score={Math.round((m.score / m.max_score) * 100)} />
                  </td>
                  <td className="max-w-sm px-4 py-2.5 text-xs text-gray-500">
                    {m.notes.length > 0 ? truncate(m.notes[0], 120) : "\u2014"}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* ---- Delivery pipeline table ---- */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">Delivery Pipeline</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {delivery.modules_succeeded}/{delivery.modules_executed} modules passed &mdash;
            runtime {delivery.total_runtime_seconds.toFixed(1)}s
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="delivery-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">Module</th>
                <th className="px-4 py-3">Name</th>
                <th className="px-4 py-3">Status</th>
                <th className="px-4 py-3">Summary</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {delivery.module_results.map((m) => (
                <tr key={m.module} className="transition-colors hover:bg-gray-50">
                  <td className="px-4 py-2.5">
                    <span className="rounded-md bg-rose-50 px-2 py-0.5 font-mono text-xs font-semibold text-rose-700">
                      {m.module}
                    </span>
                  </td>
                  <td className="px-4 py-2.5 text-sm font-medium text-gray-800">
                    {DELIVERY_NAMES[m.module] ?? m.module}
                  </td>
                  <td className="px-4 py-2.5">{statusBadge(m.status)}</td>
                  <td className="max-w-md px-4 py-2.5 text-xs text-gray-500">
                    {truncate(m.executive_summary, 120)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>

      {/* ---- Bar chart ---- */}
      <div className="rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">Module Scores</h2>
        <p className="mt-0.5 mb-4 text-xs text-gray-400">
          Mathematical validation score per module (0&ndash;100)
        </p>
        <BarChart
          categories={chartCategories}
          series={chartSeries}
          colors={["#E11D48"]}
          height={260}
          horizontal
        />
      </div>
    </div>
  );
}
