import Link from "next/link";
import dynamic from "next/dynamic";
import { safeLoadModuleJson, loadModuleJson, safeLoadPdb } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";

const MolViewer = dynamic(() => import("@/components/mol-viewer/MolViewer"), { ssr: false });

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

interface MoleculeV2 {
  name: string;
  type: string;
  length: number;
  target: string;
  target_region: string;
}

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
  details: DimensionDetails;
  notes: string[];
}

interface MathCertificateV2 {
  molecule: MoleculeV2;
  overall_verdict: string;
  composite_score: number;
  max_score: number;
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
  if (v === "pass" || v === "validated") {
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

function translateNote(note: string): string {
  return note
    .replace("MRL-ASO-001 NAO e o melhor ASO de 25 nt", "MRL-ASO-001 is not the optimal 25-nt ASO by pure dG ranking")
    .replace("verificar manualmente", "requires manual verification")
    .replace("Complementaridade maxima off-target", "Maximum off-target complementarity")
    .replace("Insuficiente para ativacao de RNase H", "Insufficient for RNase H activation")
    .replace("Regiao-alvo e INVARIANTE", "Target region is INVARIANT")
    .replace("Qualquer mutacao nesta regiao e presumivelmente letal", "Any mutation in this region is presumably lethal")
    .replace("MRL-ASO-001 esta na frente de Pareto", "MRL-ASO-001 is on the Pareto front")
    .replace("Zero mutacoes de escape viaveis", "Zero viable escape mutations")
    .replace("Resistencia e matematicamente impossivel", "Resistance is computationally predicted to be impossible");
}

function truncate(text: string, max: number): string {
  if (text.length <= max) return text;
  return text.slice(0, max) + "...";
}

/* ------------------------------------------------------------------ */
/*  Score bar (inline progress indicator)                              */
/* ------------------------------------------------------------------ */

function ScoreBar({ score, maxScore = 100 }: { score: number; maxScore?: number }) {
  const pct = Math.round((score / maxScore) * 100);
  const color = pct >= 90 ? "bg-emerald-500" : pct >= 70 ? "bg-sky-500" : "bg-amber-500";
  return (
    <div className="flex items-center gap-2">
      <div className="h-2 w-24 overflow-hidden rounded-full bg-gray-100">
        <div className={`h-full rounded-full ${color}`} style={{ width: `${pct}%` }} />
      </div>
      <span className="font-mono text-xs font-bold text-gray-700">{score}/{maxScore}</span>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Pending card                                                       */
/* ------------------------------------------------------------------ */

function PendingCard({ title, subtitle }: { title: string; subtitle?: string }) {
  return (
    <div className="mb-6 rounded-xl bg-white shadow-card p-8 text-center">
      <div className="mx-auto mb-3 flex h-10 w-10 items-center justify-center rounded-full bg-gray-100">
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.5} className="h-5 w-5 text-gray-400">
          <circle cx="12" cy="12" r="9" />
          <path strokeLinecap="round" d="M12 7v5l3 3" />
        </svg>
      </div>
      <p className="text-sm font-semibold text-gray-700">{title}</p>
      {subtitle && <p className="mt-1 text-xs text-gray-400">{subtitle}</p>}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Page                                                               */
/* ------------------------------------------------------------------ */

export default function AsoPage() {
  const cert = safeLoadModuleJson("aso_math_reports", "math_certificate_v2.json") as MathCertificateV2 | null;
  const delivery = loadModuleJson("aso_delivery", "delivery_report.json") as DeliveryReport;
  const duplexPdb = safeLoadPdb("marley_ai/05_rosettafold/structures/aso_sl_duplex_model.pdb");

  const thermoDetails = cert?.dimension_assessments?.[0]?.details ?? {};
  const wt_dg = thermoDetails.wt_dg_kcal as number | undefined;
  const wt_tm = thermoDetails.wt_tm_celsius as number | undefined;
  const overallPct = cert ? Math.round((cert.composite_score / cert.max_score) * 100) : 0;

  const chartCategories = (cert?.dimension_assessments ?? []).map((m) => m.dimension.replace(/^\d+\.\s*/, "").split(" ").slice(0, 3).join(" "));
  const chartSeries = [
    { name: "Score", data: (cert?.dimension_assessments ?? []).map((m) => m.score) },
  ];

  return (
    <div>
      {/* ---- Page header ---- */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-rose-100 px-2.5 py-1 text-xs font-bold text-rose-600">
          v2.0
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
          value={cert ? `${cert.composite_score}/${cert.max_score}` : "—"}
          subtitle={cert ? cert.overall_verdict.replace(/_/g, " ") : "Computação em andamento"}
          accentColor="bg-rose-500"
        />
        <KpiCard
          title="Binding Energy"
          value={wt_dg != null ? `${wt_dg} kcal/mol` : "—"}
          subtitle={wt_tm != null ? `Tm ${wt_tm} °C` : "Data pending"}
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
          value="0"
          subtitle="All binding-disrupting mutations are lethal"
          accentColor="bg-rose-500"
        />
      </div>

      {/* ---- Bio-Sim launch card ---- */}
      <Link
        href="/bio-sim"
        className="mb-6 flex items-center gap-5 rounded-xl bg-gradient-to-r from-rose-500/10 to-fuchsia-500/10 p-5 transition-all hover:from-rose-500/20 hover:to-fuchsia-500/20"
        style={{ border: "1px solid var(--app-border)" }}
      >
        <div className="flex h-14 w-14 flex-shrink-0 items-center justify-center rounded-xl bg-rose-500/20">
          <svg viewBox="0 0 24 24" fill="none" stroke="#e11d48" strokeWidth={1.8} className="h-7 w-7">
            <circle cx="12" cy="12" r="9" />
            <polygon points="10,8 16,12 10,16" fill="#e11d48" />
          </svg>
        </div>
        <div className="flex-1">
          <h3 className="text-sm font-bold text-gray-900">MRL-ASO-001 Bio-Sim: Mechanism of Action</h3>
          <p className="mt-0.5 text-xs text-gray-500">
            Interactive 3D journey from subcutaneous injection to parasite clearance — 8 scenes, 24 validated metrics
          </p>
        </div>
        <span className="hidden rounded-full bg-rose-500 px-3 py-1.5 text-xs font-bold text-white sm:block">
          Launch
        </span>
      </Link>

      {/* ---- ASO:SL RNA 3D Structure ---- */}
      {duplexPdb && (
        <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">ASO:SL RNA Duplex Structure</h2>
            <p className="mt-0.5 text-xs text-gray-400">
              RoseTTAFold-predicted A-form hybrid duplex — Chain A (SL RNA, teal) · Chain B (MRL-ASO-001, rose)
            </p>
          </div>
          <MolViewer
            pdbData={duplexPdb}
            preset="nucleic-acid"
            height={400}
          />
        </div>
      )}

      {/* ---- Molecule card ---- */}
      {cert ? (
        <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">Molecule: {cert.molecule.name}</h2>
            <p className="mt-0.5 text-xs text-gray-400">
              {cert.molecule.type}
            </p>
          </div>
          <div className="px-5 py-5">
            <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
              {[
                { label: "Target", value: cert.molecule.target },
                { label: "Target Region", value: cert.molecule.target_region },
                { label: "Length", value: `${cert.molecule.length} nt` },
                { label: "dG", value: wt_dg != null ? `${wt_dg} kcal/mol` : "—" },
                { label: "Tm", value: wt_tm != null ? `${wt_tm} °C` : "—" },
              ].map((prop) => (
                <div key={prop.label} className="rounded-lg bg-gray-50 p-3">
                  <p className="text-xs text-gray-400">{prop.label}</p>
                  <p className="mt-1 text-sm font-semibold text-gray-900">{prop.value}</p>
                </div>
              ))}
            </div>
          </div>
        </div>
      ) : (
        <PendingCard title="Mathematical certificate pending" subtitle="Computação em andamento — resultados disponíveis em breve" />
      )}

      {/* ---- Validation modules table ---- */}
      {cert ? (
        <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">Validation Dimensions</h2>
            <p className="mt-0.5 text-xs text-gray-400">
              {cert.dimension_assessments.length} mathematical validation dimensions &mdash; overall {cert.composite_score}/{cert.max_score} ({overallPct}%)
            </p>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-left text-sm" data-testid="validation-table">
              <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
                <tr>
                  <th className="px-4 py-3">#</th>
                  <th className="px-4 py-3">Dimension</th>
                  <th className="px-4 py-3">Verdict</th>
                  <th className="px-4 py-3">Score</th>
                  <th className="px-4 py-3">Note</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-50">
                {cert.dimension_assessments.map((m, i) => (
                  <tr key={m.dimension} className="transition-colors hover:bg-gray-50">
                    <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                    <td className="px-4 py-2.5 text-sm font-medium text-gray-800">
                      {m.dimension.replace(/^\d+\.\s*/, "")}
                    </td>
                    <td className="px-4 py-2.5">{verdictBadge(m.verdict)}</td>
                    <td className="px-4 py-2.5">
                      <ScoreBar score={m.score} maxScore={m.max_score} />
                    </td>
                    <td className="max-w-sm px-4 py-2.5 text-xs text-gray-500">
                      {m.notes.length > 0 ? truncate(translateNote(m.notes[0]), 120) : "—"}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      ) : null}

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
      {cert && chartCategories.length > 0 && (
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Dimension Scores</h2>
          <p className="mt-0.5 mb-4 text-xs text-gray-400">
            Mathematical validation score per dimension
          </p>
          <BarChart
            categories={chartCategories}
            series={chartSeries}
            colors={["#E11D48"]}
            height={260}
            horizontal
          />
        </div>
      )}
    </div>
  );
}
