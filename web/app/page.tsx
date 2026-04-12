import Link from "next/link";
import KpiCard from "@/components/kpi-card";
import DonutChart from "@/components/charts/donut-chart";

/* ------------------------------------------------------------------ */
/*  KPI data                                                           */
/* ------------------------------------------------------------------ */

const kpisRow1 = [
  {
    title: "Proteins Analyzed",
    value: "8,527",
    subtitle: "L. infantum proteome",
    accentColor: "bg-cyan-500",
  },
  {
    title: "ASO Validation",
    value: "90/100",
    subtitle: "MRL-ASO-001 certificate",
    accentColor: "bg-rose-500",
  },
  {
    title: "Drug Targets",
    value: 52,
    subtitle: "Prioritized by druggability",
    accentColor: "bg-orange-500",
  },
  {
    title: "AI Modules",
    value: 11,
    subtitle: "6 operational",
    accentColor: "bg-sky-500",
  },
];

const kpisRow2 = [
  {
    title: "Epitopes",
    value: 11,
    subtitle: "MHC-I high-affinity",
    accentColor: "bg-blue-500",
  },
  {
    title: "Best Docking",
    value: "-8.07 kcal/mol",
    subtitle: "CHEMBL92 on GMPS",
    accentColor: "bg-emerald-500",
  },
  {
    title: "Quantum Optimum",
    value: "-30.42 kcal/mol",
    subtitle: "QAOA all-Sp config",
    accentColor: "bg-indigo-500",
  },
  {
    title: "Platforms",
    value: 3,
    subtitle: "mRNA, E. coli, L. tarentolae",
    accentColor: "bg-teal-500",
  },
];

/* ------------------------------------------------------------------ */
/*  Pipeline modules by track                                          */
/* ------------------------------------------------------------------ */

interface PipelineItem {
  title: string;
  href: string;
  description: string;
  track: string;
  badgeColor: string;
}

const discoveryModules: PipelineItem[] = [
  {
    title: "Vaccine Construct",
    href: "/vaccine",
    description:
      "Multi-epitope mRNA vaccine with 11 epitopes, VaxiJen 0.3235",
    track: "Discovery",
    badgeColor: "bg-cyan-100 text-cyan-700",
  },
  {
    title: "Drug Targets",
    href: "/drug",
    description:
      "52 druggable parasite-specific targets across 4 pathways",
    track: "Discovery",
    badgeColor: "bg-cyan-100 text-cyan-700",
  },
  {
    title: "Molecular Docking",
    href: "/docking",
    description:
      "77 compounds screened, best hit CHEMBL92 at -8.07 kcal/mol",
    track: "Discovery",
    badgeColor: "bg-cyan-100 text-cyan-700",
  },
  {
    title: "Immune Simulation",
    href: "/simulation",
    description:
      "Th1 dominance 82.4%, memory protection >693 days",
    track: "Discovery",
    badgeColor: "bg-cyan-100 text-cyan-700",
  },
];

const therapeuticsModules: PipelineItem[] = [
  {
    title: "ASO Therapy",
    href: "/aso",
    description:
      "MRL-ASO-001 validated 90/100, resistance mathematically impossible",
    track: "Therapeutics",
    badgeColor: "bg-rose-100 text-rose-700",
  },
  {
    title: "RNA Entropy",
    href: "/rna",
    description:
      "Information theory analysis of 500 transcripts, SL RNA target confirmed",
    track: "Therapeutics",
    badgeColor: "bg-rose-100 text-rose-700",
  },
  {
    title: "Quantum Computing",
    href: "/quantum",
    description:
      "QAOA+VQE optimization, -30.42 kcal/mol optimal stereochemistry",
    track: "Therapeutics",
    badgeColor: "bg-rose-100 text-rose-700",
  },
];

const intelligenceModules: PipelineItem[] = [
  {
    title: "AI/ML Engine",
    href: "/ai",
    description:
      "11 modules, 5 autonomous agents, 27 insights, consensus 0.60",
    track: "Intelligence",
    badgeColor: "bg-sky-100 text-sky-700",
  },
  {
    title: "Vaccine Platforms",
    href: "/platforms",
    description:
      "3 production systems compared across 6 dimensions",
    track: "Intelligence",
    badgeColor: "bg-sky-100 text-sky-700",
  },
];

/* ------------------------------------------------------------------ */
/*  Chart data                                                         */
/* ------------------------------------------------------------------ */

const donutSeries = [4, 3, 2];
const donutLabels = ["Discovery", "Therapeutics", "Intelligence"];
const donutColors = ["#06B6D4", "#F43F5E", "#0EA5E9"];

/* ------------------------------------------------------------------ */
/*  Pipeline card component                                            */
/* ------------------------------------------------------------------ */

function PipelineCard({ item }: { item: PipelineItem }) {
  return (
    <Link
      href={item.href}
      className="rounded-xl bg-white shadow-card p-5 transition-all hover:shadow-md"
      data-testid={`pipeline-${item.href.slice(1)}`}
    >
      <div className="mb-2 flex items-center gap-2">
        <span
          className={`rounded-lg ${item.badgeColor} px-2 py-0.5 text-xs font-bold`}
        >
          {item.track}
        </span>
        <h3 className="text-sm font-semibold text-gray-900">{item.title}</h3>
      </div>
      <p className="text-xs leading-relaxed text-gray-500">
        {item.description}
      </p>
    </Link>
  );
}

/* ------------------------------------------------------------------ */
/*  Page                                                               */
/* ------------------------------------------------------------------ */

export default function Home() {
  return (
    <div>
      {/* Page header */}
      <div className="mb-6">
        <h1 className="text-2xl font-bold text-gray-900">Marley Dashboard</h1>
        <p className="mt-1 text-sm text-gray-500">
          Reverse vaccinology and drug discovery for canine visceral
          leishmaniasis.
        </p>
      </div>

      {/* KPI Grid - Row 1 */}
      <div
        className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4"
        data-testid="kpi-grid"
      >
        {kpisRow1.map((kpi) => (
          <KpiCard
            key={kpi.title}
            title={kpi.title}
            value={kpi.value}
            subtitle={kpi.subtitle}
            accentColor={kpi.accentColor}
          />
        ))}
      </div>

      {/* KPI Grid - Row 2 */}
      <div
        className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4"
        data-testid="kpi-grid-2"
      >
        {kpisRow2.map((kpi) => (
          <KpiCard
            key={kpi.title}
            title={kpi.title}
            value={kpi.value}
            subtitle={kpi.subtitle}
            accentColor={kpi.accentColor}
          />
        ))}
      </div>

      {/* Track overview chart */}
      <div className="mb-6 rounded-xl bg-white shadow-card overflow-hidden">
        <div className="p-5">
          <h2 className="mb-1 text-sm font-semibold text-gray-900">
            Track Overview
          </h2>
          <p className="mb-4 text-xs text-gray-400">
            Modules per research track
          </p>
          <DonutChart
            series={donutSeries}
            labels={donutLabels}
            colors={donutColors}
            height={260}
          />
        </div>
      </div>

      {/* Pipeline Modules */}
      <div className="mb-6">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-gray-400">
          Pipeline Modules
        </h2>

        {/* Discovery Pipeline */}
        <p className="mb-2 px-1 text-xs font-semibold uppercase tracking-wider text-gray-400">
          Discovery Pipeline
        </p>
        <div className="mb-6 grid gap-4 md:grid-cols-2 xl:grid-cols-3">
          {discoveryModules.map((item) => (
            <PipelineCard key={item.href} item={item} />
          ))}
        </div>

        {/* Therapeutics Pipeline */}
        <p className="mb-2 px-1 text-xs font-semibold uppercase tracking-wider text-gray-400">
          Therapeutics Pipeline
        </p>
        <div className="mb-6 grid gap-4 md:grid-cols-2 xl:grid-cols-3">
          {therapeuticsModules.map((item) => (
            <PipelineCard key={item.href} item={item} />
          ))}
        </div>

        {/* Intelligence Pipeline */}
        <p className="mb-2 px-1 text-xs font-semibold uppercase tracking-wider text-gray-400">
          Intelligence Pipeline
        </p>
        <div className="mb-6 grid gap-4 md:grid-cols-2 xl:grid-cols-3">
          {intelligenceModules.map((item) => (
            <PipelineCard key={item.href} item={item} />
          ))}
        </div>
      </div>
    </div>
  );
}
