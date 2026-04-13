import Link from "next/link";
import { getTranslations } from "next-intl/server";
import KpiCard from "@/components/kpi-card";

/* ------------------------------------------------------------------ */
/*  Pipeline card component                                            */
/* ------------------------------------------------------------------ */

interface PipelineItem {
  title: string;
  href: string;
  description: string;
  track: string;
  badgeColor: string;
}

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

export default async function Home() {
  const t = await getTranslations("dashboard");

  /* ---- KPI data ---- */

  const kpisRow1 = [
    {
      title: t("kpi.proteinsAnalyzed.title"),
      value: "8,527",
      subtitle: t("kpi.proteinsAnalyzed.subtitle", { count: 52 }),
      accentColor: "bg-cyan-500",
    },
    {
      title: t("kpi.asoValidation.title"),
      value: "90/100",
      subtitle: t("kpi.asoValidation.subtitle"),
      accentColor: "bg-rose-500",
    },
    {
      title: t("kpi.drugTargets.title"),
      value: 52,
      subtitle: t("kpi.drugTargets.subtitle"),
      accentColor: "bg-orange-500",
    },
    {
      title: t("kpi.epitopes.title"),
      value: 11,
      subtitle: t("kpi.epitopes.subtitle", { count: 11 }),
      accentColor: "bg-blue-500",
    },
  ];

  const kpisRow2 = [
    {
      title: t("kpi.dockingHits.title"),
      value: 10,
      subtitle: t("kpi.dockingHits.subtitle"),
      accentColor: "bg-emerald-500",
    },
    {
      title: t("kpi.asoTarget.title"),
      value: "SL RNA",
      subtitle: t("kpi.asoTarget.subtitle"),
      accentColor: "bg-teal-500",
    },
    {
      title: t("kpi.bestDocking.title"),
      value: "-8.07 kcal/mol",
      subtitle: t("kpi.bestDocking.subtitle"),
      accentColor: "bg-emerald-500",
    },
    {
      title: t("kpi.platforms.title"),
      value: 3,
      subtitle: t("kpi.platforms.subtitle"),
      accentColor: "bg-teal-500",
    },
  ];

  /* ---- Pipeline modules by track ---- */

  const discoveryModules: PipelineItem[] = [
    {
      title: t("pipeline.vaccineConstruct.title"),
      href: "/vaccine",
      description: t("pipeline.vaccineConstruct.description"),
      track: t("pipeline.vaccineConstruct.track"),
      badgeColor: "bg-cyan-100 text-cyan-700",
    },
    {
      title: t("pipeline.drugTargets.title"),
      href: "/drug",
      description: t("pipeline.drugTargets.description"),
      track: t("pipeline.drugTargets.track"),
      badgeColor: "bg-cyan-100 text-cyan-700",
    },
    {
      title: t("pipeline.molecularDocking.title"),
      href: "/docking",
      description: t("pipeline.molecularDocking.description"),
      track: t("pipeline.molecularDocking.track"),
      badgeColor: "bg-cyan-100 text-cyan-700",
    },
  ];

  const therapeuticsModules: PipelineItem[] = [
    {
      title: t("pipeline.asoTherapy.title"),
      href: "/aso",
      description: t("pipeline.asoTherapy.description"),
      track: t("pipeline.asoTherapy.track"),
      badgeColor: "bg-rose-100 text-rose-700",
    },
    {
      title: t("pipeline.targetValidation.title"),
      href: "/rna",
      description: t("pipeline.targetValidation.description"),
      track: t("pipeline.targetValidation.track"),
      badgeColor: "bg-rose-100 text-rose-700",
    },
    {
      title: t("pipeline.vaccinePlatforms.title"),
      href: "/platforms",
      description: t("pipeline.vaccinePlatforms.description"),
      track: t("pipeline.vaccinePlatforms.track"),
      badgeColor: "bg-rose-100 text-rose-700",
    },
  ];

  return (
    <div>
      {/* Page header */}
      <div className="mb-6">
        <h1 className="text-2xl font-bold text-gray-900">{t("title")}</h1>
        <p className="mt-1 text-sm text-gray-500">
          {t("subtitle")}
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

      {/* Pipeline Modules */}
      <div className="mb-6">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-gray-400">
          {t("sectionPipelineModules")}
        </h2>

        {/* Discovery Pipeline */}
        <p className="mb-2 px-1 text-xs font-semibold uppercase tracking-wider text-gray-400">
          {t("sectionDiscovery")}
        </p>
        <div className="mb-6 grid gap-4 md:grid-cols-2 xl:grid-cols-3">
          {discoveryModules.map((item) => (
            <PipelineCard key={item.href} item={item} />
          ))}
        </div>

        {/* Therapeutics */}
        <p className="mb-2 px-1 text-xs font-semibold uppercase tracking-wider text-gray-400">
          {t("sectionTherapeutics")}
        </p>
        <div className="mb-6 grid gap-4 md:grid-cols-2 xl:grid-cols-3">
          {therapeuticsModules.map((item) => (
            <PipelineCard key={item.href} item={item} />
          ))}
        </div>
      </div>
    </div>
  );
}
