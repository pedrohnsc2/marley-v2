import Link from "next/link";
import KpiCard from "@/components/kpi-card";

const kpis = [
  {
    title: "Proteins Analyzed",
    value: "8,527",
    subtitle: "L. infantum proteome",
    accent: "border-cyan-500",
  },
  {
    title: "Epitopes Selected",
    value: 11,
    subtitle: "MHC-I high-affinity binders",
    accent: "border-cyan-400",
  },
  {
    title: "Drug Targets",
    value: 52,
    subtitle: "Prioritized by druggability",
    accent: "border-orange-500",
  },
  {
    title: "Compounds Docked",
    value: 77,
    subtitle: "Repurposed + custom molecules",
    accent: "border-green-500",
  },
  {
    title: "Custom Molecules",
    value: 20,
    subtitle: "MRL-series designed in silico",
    accent: "border-orange-400",
  },
  {
    title: "Best Affinity",
    value: "-7.74 kcal/mol",
    subtitle: "MRL-003 on TryR",
    accent: "border-green-400",
  },
];

const pipelineSteps = [
  {
    id: "v1",
    title: "Proteome Screening",
    description:
      "Filtered 8,527 proteins through signal peptide detection, conservation analysis, and immunogenicity scoring to identify top vaccine candidates.",
    color: "bg-cyan-900/30 border-cyan-800",
    href: "/vaccine",
  },
  {
    id: "v2",
    title: "Vaccine Construct",
    description:
      "Designed multi-epitope mRNA vaccine with 15 epitopes, L7/L12 adjuvant, tPA signal peptide, and codon-optimized sequence validated via ESMFold 3D prediction.",
    color: "bg-cyan-900/30 border-cyan-800",
    href: "/vaccine",
  },
  {
    id: "v3",
    title: "Drug Target Discovery",
    description:
      "Identified 52 druggable parasite-specific targets across trypanothione, sterol, folate, and purine salvage pathways with low human homology.",
    color: "bg-orange-900/30 border-orange-800",
    href: "/drug",
  },
  {
    id: "v4",
    title: "Molecular Docking",
    description:
      "Screened 77 compounds via AutoDock Vina against GMPS and TryR. Top hit: CHEMBL92 at -8.07 kcal/mol. Designed 20 custom MRL-series molecules.",
    color: "bg-green-900/30 border-green-800",
    href: "/docking",
  },
  {
    id: "v5",
    title: "Immune Simulation",
    description:
      "ODE-based immune kinetics model predicting Th1 dominance (82.4%), robust memory cell formation, and estimated protection beyond 693 days.",
    color: "bg-purple-900/30 border-purple-800",
    href: "/simulation",
  },
];

export default function Home() {
  return (
    <div className="px-8 py-10">
      <header className="mb-10">
        <h1 className="text-3xl font-bold tracking-tight text-white">
          Marley Dashboard
        </h1>
        <p className="mt-2 max-w-2xl text-sm text-zinc-400">
          Reverse vaccinology and drug discovery pipeline for canine visceral
          leishmaniasis. Five computational modules spanning proteome analysis
          through immune simulation.
        </p>
      </header>

      <section className="mb-12">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          Key Metrics
        </h2>
        <div
          className="grid grid-cols-2 gap-4 lg:grid-cols-3 xl:grid-cols-6"
          data-testid="kpi-grid"
        >
          {kpis.map((kpi) => (
            <KpiCard
              key={kpi.title}
              title={kpi.title}
              value={kpi.value}
              subtitle={kpi.subtitle}
              accentColor={kpi.accent}
            />
          ))}
        </div>
      </section>

      <section>
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          Pipeline Modules
        </h2>
        <div className="grid gap-4 md:grid-cols-2 xl:grid-cols-3">
          {pipelineSteps.map((step) => (
            <Link
              key={step.id}
              href={step.href}
              className={`rounded-lg border ${step.color} p-5 transition-colors hover:brightness-125`}
              data-testid={`pipeline-${step.id}`}
            >
              <div className="mb-2 flex items-center gap-2">
                <span className="rounded bg-zinc-800 px-2 py-0.5 text-xs font-mono text-zinc-500">
                  {step.id}
                </span>
                <h3 className="text-sm font-semibold text-white">
                  {step.title}
                </h3>
              </div>
              <p className="text-xs leading-relaxed text-zinc-400">
                {step.description}
              </p>
            </Link>
          ))}
        </div>
      </section>
    </div>
  );
}
