import Link from "next/link";
import dynamic from "next/dynamic";
import { getTranslations } from "next-intl/server";
import { loadModuleJson, safeLoadModuleJson, safeLoadPdb, safeLoadRunModuleJson, safeLoadRunPdb, getRunCompletedDate } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import BarChart from "@/components/charts/bar-chart";
import RunBanner from "@/components/runs/run-banner";

const MolViewer = dynamic(() => import("@/components/mol-viewer/MolViewer"), { ssr: false });

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

interface Molecule {
  name: string;
  sequence: string;
  length: number;
  target: string;
  target_sequence: string;
  target_region: string;
  known_dg_kcal: number;
  known_tm_celsius: number;
}

interface ModuleAssessment {
  module: string;
  verdict: string;
  score: number;
  notes: string[];
}

interface MathCertificate {
  molecule: Molecule;
  verdict: string;
  confidence: string;
  overall_score: number;
  module_assessments: ModuleAssessment[];
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

export default async function AsoPage({
  searchParams,
}: {
  searchParams: Promise<{ run?: string }>;
}) {
  const { run: runId } = await searchParams;
  const t = await getTranslations("aso");

  const cert = (
    runId
      ? safeLoadRunModuleJson(runId, "aso_math", "math_certificate.json") ?? safeLoadModuleJson("aso_math", "math_certificate.json")
      : safeLoadModuleJson("aso_math", "math_certificate.json")
  ) as MathCertificate | null;

  const delivery = (
    runId
      ? safeLoadRunModuleJson(runId, "aso_delivery", "delivery_report.json") ?? safeLoadModuleJson("aso_delivery", "delivery_report.json")
      : safeLoadModuleJson("aso_delivery", "delivery_report.json")
  ) as DeliveryReport | null;

  const duplexPdb = runId
    ? safeLoadRunPdb(runId, "marley_ai/05_rosettafold/structures/aso_sl_duplex_model.pdb") ?? safeLoadPdb("marley_ai/05_rosettafold/structures/aso_sl_duplex_model.pdb")
    : safeLoadPdb("marley_ai/05_rosettafold/structures/aso_sl_duplex_model.pdb");

  const runCompletedAt = runId ? getRunCompletedDate(runId) : null;

  if (!cert || !delivery) {
    return (
      <div className="flex min-h-[60vh] items-center justify-center">
        <div className="text-center" style={{ color: "var(--app-text-2)" }}>
          <h2 className="text-xl font-bold mb-2">Data not available</h2>
          <p className="text-sm">Pipeline results could not be loaded. Run the ASO pipeline first.</p>
          <Link href="/runs" className="mt-4 inline-block text-sm underline" style={{ color: "var(--app-primary)" }}>
            Go to Runs
          </Link>
        </div>
      </div>
    );
  }

  const { molecule, module_assessments } = cert;

  /* Module name mapping from translations */
  const MODULE_NAMES: Record<string, string> = {
    "01_thermodynamic_landscape": t("moduleNames.01_thermodynamic_landscape"),
    "02_selectivity_proof": t("moduleNames.02_selectivity_proof"),
    "03_evolutionary_conservation": t("moduleNames.03_evolutionary_conservation"),
    "04_exhaustive_optimization": t("moduleNames.04_exhaustive_optimization"),
    "05_resistance_model": t("moduleNames.05_resistance_model"),
  };

  /* Delivery name mapping from translations */
  const DELIVERY_NAMES: Record<string, string> = {
    A: t("deliveryNames.A"),
    B: t("deliveryNames.B"),
    C: t("deliveryNames.C"),
    D: t("deliveryNames.D"),
    E: t("deliveryNames.E"),
    F: t("deliveryNames.F"),
  };

  /* Bar chart data */
  const chartCategories = module_assessments.map(
    (m) => MODULE_NAMES[m.module] ?? m.module,
  );
  const chartSeries = [
    { name: t("moduleScoresChart.seriesName"), data: module_assessments.map((m) => m.score) },
  ];

  return (
    <div>
      <RunBanner runId={runId ?? null} pipeline="aso" completedAt={runCompletedAt} />

      {/* ---- Page header ---- */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-rose-100 px-2.5 py-1 text-xs font-bold text-rose-600">
          {t("badge")}
        </span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">{t("title")}</h1>
          <p className="text-sm text-gray-500">
            {t.raw("subtitle")}
          </p>
        </div>
      </div>

      {/* ---- KPI row ---- */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title={t("kpi.mathScore.title")}
          value={`${cert.overall_score}/100`}
          subtitle={cert.verdict.replace(/_/g, " ")}
          accentColor="bg-rose-500"
        />
        <KpiCard
          title={t("kpi.bindingEnergy.title")}
          value={`${molecule.known_dg_kcal} kcal/mol`}
          subtitle={t("kpi.bindingEnergy.subtitle", { tm: molecule.known_tm_celsius })}
          accentColor="bg-pink-500"
        />
        <KpiCard
          title={t("kpi.delivery.title")}
          value={`${delivery.modules_succeeded}/${delivery.modules_executed}`}
          subtitle={t("kpi.delivery.subtitle")}
          accentColor="bg-fuchsia-500"
        />
        <KpiCard
          title={t("kpi.resistance.title")}
          value={t("kpi.resistance.value")}
          subtitle={t("kpi.resistance.subtitle")}
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
          <h3 className="text-sm font-bold text-gray-900">{t("bioSim.title")}</h3>
          <p className="mt-0.5 text-xs text-gray-500">
            {t("bioSim.subtitle")}
          </p>
        </div>
        <span className="hidden rounded-full bg-rose-500 px-3 py-1.5 text-xs font-bold text-white sm:block">
          {t("bioSim.launch")}
        </span>
      </Link>

      {/* ---- ASO:SL RNA 3D Structure ---- */}
      {duplexPdb && (
        <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">{t("duplexStructure.title")}</h2>
            <p className="mt-0.5 text-xs text-gray-400">
              {t("duplexStructure.subtitle")}
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
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="border-b border-gray-100 px-5 py-4">
          <h2 className="text-sm font-semibold text-gray-900">{t("molecule.title", { name: molecule.name })}</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {t("molecule.subtitle", { length: molecule.length })}
          </p>
        </div>
        <div className="px-5 py-5">
          {/* Sequence display */}
          <div className="mb-5 rounded-lg bg-rose-50 px-4 py-3">
            <p className="mb-1 text-xs font-semibold uppercase tracking-wider text-rose-400">
              {t("molecule.sequenceLabel")}
            </p>
            <p
              className="font-mono text-lg font-bold tracking-widest text-rose-700"
              data-testid="aso-sequence"
            >
              {molecule.sequence}
            </p>
          </div>

          {/* Properties grid */}
          <div className="grid grid-cols-2 gap-4 sm:grid-cols-3 lg:grid-cols-5">
            {[
              { label: t("molecule.target"), value: molecule.target },
              { label: t("molecule.targetRegion"), value: molecule.target_region },
              { label: t("molecule.length"), value: `${molecule.length} nt` },
              { label: t("molecule.dG"), value: `${molecule.known_dg_kcal} kcal/mol` },
              { label: t("molecule.tm"), value: `${molecule.known_tm_celsius} C` },
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
          <h2 className="text-sm font-semibold text-gray-900">{t("validationTable.title")}</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {t("validationTable.subtitle", { score: cert.overall_score })}
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="validation-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">#</th>
                <th className="px-4 py-3">{t("validationTable.module")}</th>
                <th className="px-4 py-3">{t("validationTable.verdict")}</th>
                <th className="px-4 py-3">{t("validationTable.score")}</th>
                <th className="px-4 py-3">{t("validationTable.note")}</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-50">
              {module_assessments.map((m, i) => (
                <tr key={m.module} className="transition-colors hover:bg-gray-50">
                  <td className="px-4 py-2.5 font-mono text-xs text-gray-300">{i + 1}</td>
                  <td className="px-4 py-2.5 text-sm font-medium text-gray-800">
                    {MODULE_NAMES[m.module] ?? m.module}
                  </td>
                  <td className="px-4 py-2.5">{verdictBadge(m.verdict)}</td>
                  <td className="px-4 py-2.5">
                    <ScoreBar score={m.score} />
                  </td>
                  <td className="max-w-sm px-4 py-2.5 text-xs text-gray-500">
                    {m.notes.length > 0 ? truncate(translateNote(m.notes[0]), 120) : "\u2014"}
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
          <h2 className="text-sm font-semibold text-gray-900">{t("deliveryTable.title")}</h2>
          <p className="mt-0.5 text-xs text-gray-400">
            {t("deliveryTable.subtitle", {
              passed: delivery.modules_succeeded,
              total: delivery.modules_executed,
              runtime: delivery.total_runtime_seconds.toFixed(1),
            })}
          </p>
        </div>
        <div className="overflow-x-auto">
          <table className="w-full text-left text-sm" data-testid="delivery-table">
            <thead className="border-b border-gray-100 bg-gray-50 text-xs uppercase tracking-wider text-gray-400">
              <tr>
                <th className="px-4 py-3">{t("deliveryTable.module")}</th>
                <th className="px-4 py-3">{t("deliveryTable.name")}</th>
                <th className="px-4 py-3">{t("deliveryTable.status")}</th>
                <th className="px-4 py-3">{t("deliveryTable.summary")}</th>
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
        <h2 className="text-sm font-semibold text-gray-900">{t("moduleScoresChart.title")}</h2>
        <p className="mt-0.5 mb-4 text-xs text-gray-400">
          {t("moduleScoresChart.subtitle")}
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
