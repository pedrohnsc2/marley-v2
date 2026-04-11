import Link from "next/link";
import KpiCard from "@/components/kpi-card";
import DonutChart from "@/components/charts/donut-chart";
import BarChart from "@/components/charts/bar-chart";

const kpis = [
  {
    title: "Proteins Analyzed",
    value: "8,527",
    subtitle: "L. infantum proteome",
    accentColor: "bg-blue-500",
  },
  {
    title: "Epitopes Selected",
    value: 11,
    subtitle: "MHC-I high-affinity binders",
    accentColor: "bg-cyan-500",
  },
  {
    title: "Drug Targets",
    value: 52,
    subtitle: "Prioritized by druggability",
    accentColor: "bg-orange-500",
  },
  {
    title: "Compounds Docked",
    value: 77,
    subtitle: "Repurposed + custom molecules",
    accentColor: "bg-emerald-500",
  },
  {
    title: "Custom Molecules",
    value: 20,
    subtitle: "MRL-series designed in silico",
    accentColor: "bg-amber-500",
  },
  {
    title: "Best Affinity",
    value: "-7.74 kcal/mol",
    subtitle: "MRL-003 on TryR",
    accentColor: "bg-green-500",
  },
];

const pipelineSteps = [
  {
    id: "v1-v2",
    title: "Vaccine Construct",
    description:
      "Proteome screening + multi-epitope mRNA design with L7/L12 adjuvant, tPA signal peptide, and ESMFold 3D prediction.",
    color: "border-l-blue-500",
    badge: "bg-blue-50 text-blue-600",
    href: "/vaccine",
  },
  {
    id: "v3",
    title: "Drug Target Discovery",
    description:
      "52 druggable parasite-specific targets across trypanothione, sterol, folate, and purine salvage pathways.",
    color: "border-l-orange-500",
    badge: "bg-orange-50 text-orange-600",
    href: "/drug",
  },
  {
    id: "v4",
    title: "Molecular Docking",
    description:
      "77 compounds screened via AutoDock Vina. Top hit: CHEMBL92 at -8.07 kcal/mol. 20 custom MRL-series molecules.",
    color: "border-l-emerald-500",
    badge: "bg-emerald-50 text-emerald-600",
    href: "/docking",
  },
  {
    id: "v5",
    title: "Immune Simulation",
    description:
      "ODE-based kinetics predicting Th1 dominance (82.4%), memory formation, and protection beyond 693 days.",
    color: "border-l-purple-500",
    badge: "bg-purple-50 text-purple-600",
    href: "/simulation",
  },
];

// Chart data
const donutSeries = [8527, 11, 52, 77, 20];
const donutLabels = ["Proteins", "Epitopes", "Drug Targets", "Compounds", "Custom Mol."];
const donutColors = ["#3B82F6", "#06B6D4", "#F97316", "#10B981", "#F59E0B"];

const barCategories = ["Proteome\nScreening", "Vaccine\nConstruct", "Drug\nTargets", "Mol.\nDocking", "Immune\nSim"];
const barSeries = [
  { name: "Proteins Analyzed", data: [8527, 0, 0, 0, 0] },
  { name: "Drug Targets / Compounds", data: [0, 11, 52, 77, 0] },
  { name: "Sim Days", data: [0, 0, 0, 0, 365] },
];
const barColors = ["#3B82F6", "#F97316", "#8B5CF6"];

export default function Home() {
  return (
    <div>
      {/* Page header */}
      <div className="mb-6">
        <h1 className="text-2xl font-bold text-gray-900">Dashboard</h1>
        <p className="mt-1 text-sm text-gray-500">
          Reverse vaccinology and drug discovery pipeline for canine visceral leishmaniasis.
        </p>
      </div>

      {/* KPI Grid */}
      <div
        className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-3 xl:grid-cols-6"
        data-testid="kpi-grid"
      >
        {kpis.map((kpi) => (
          <KpiCard
            key={kpi.title}
            title={kpi.title}
            value={kpi.value}
            subtitle={kpi.subtitle}
            accentColor={kpi.accentColor}
          />
        ))}
      </div>

      {/* Charts row */}
      <div className="mb-6 grid gap-4 lg:grid-cols-5">
        {/* Donut chart - pipeline overview */}
        <div className="rounded-xl bg-white p-5 shadow-card lg:col-span-2">
          <h2 className="mb-1 text-sm font-semibold text-gray-900">Pipeline Overview</h2>
          <p className="mb-4 text-xs text-gray-400">Entities per module</p>
          <DonutChart
            series={donutSeries}
            labels={donutLabels}
            colors={donutColors}
            height={280}
          />
        </div>

        {/* Bar chart - pipeline scale */}
        <div className="rounded-xl bg-white p-5 shadow-card lg:col-span-3">
          <h2 className="mb-1 text-sm font-semibold text-gray-900">Pipeline Scale</h2>
          <p className="mb-4 text-xs text-gray-400">Key quantities per module</p>
          <BarChart
            categories={barCategories}
            series={barSeries}
            colors={barColors}
            height={280}
            stacked={false}
          />
        </div>
      </div>

      {/* Pipeline Modules */}
      <div>
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-gray-400">
          Pipeline Modules
        </h2>
        <div className="grid gap-4 md:grid-cols-2 xl:grid-cols-4">
          {pipelineSteps.map((step) => (
            <Link
              key={step.id}
              href={step.href}
              className={`group rounded-xl border border-gray-100 border-l-4 ${step.color} bg-white p-5 shadow-card transition-shadow hover:shadow-card-hover`}
              data-testid={`pipeline-${step.id}`}
            >
              <div className="mb-3 flex items-center gap-2">
                <span className={`rounded-md px-2 py-0.5 text-xs font-semibold ${step.badge}`}>
                  {step.id}
                </span>
                <h3 className="text-sm font-semibold text-gray-900 group-hover:text-blue-600 transition-colors">
                  {step.title}
                </h3>
              </div>
              <p className="text-xs leading-relaxed text-gray-500">{step.description}</p>
              <div className="mt-4 flex items-center gap-1 text-xs font-medium text-blue-500 opacity-0 group-hover:opacity-100 transition-opacity">
                View details
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} className="h-3.5 w-3.5">
                  <path d="M5 12h14M12 5l7 7-7 7" />
                </svg>
              </div>
            </Link>
          ))}
        </div>
      </div>
    </div>
  );
}
