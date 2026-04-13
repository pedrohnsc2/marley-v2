import { getTranslations } from "next-intl/server";
import { listRuns } from "@/lib/data-loader";
import { Link } from "@/i18n/routing";
import KpiCard from "@/components/kpi-card";
import PipelineFilter from "./pipeline-filter";

interface RunStage {
  stage_id: string;
  name: string;
  status: string;
  duration_s: number;
}

interface RunMeta {
  run_id: string;
  pipeline: string;
  status: string;
  created_at: string;
  git_sha: string;
  parameters: Record<string, unknown>;
  stages: RunStage[];
  total_duration_s: number;
  tags: string[];
  notes: string;
}

const PIPELINE_ROUTES: Record<string, string> = {
  vaccine: "/vaccine",
  drug: "/drug",
  docking: "/docking",
  rna: "/rna",
  aso_math: "/aso",
};

function truncateRunId(runId: string): string {
  if (runId.length <= 20) return runId;
  return runId.slice(0, 8) + "..." + runId.slice(-8);
}

function formatDuration(seconds: number): string {
  if (seconds < 60) return `${seconds.toFixed(1)}s`;
  const min = Math.floor(seconds / 60);
  const sec = seconds % 60;
  return `${min}m ${sec.toFixed(0)}s`;
}

function formatDate(iso: string): string {
  try {
    const d = new Date(iso);
    return d.toLocaleDateString(undefined, {
      year: "numeric",
      month: "short",
      day: "numeric",
      hour: "2-digit",
      minute: "2-digit",
    });
  } catch {
    return iso;
  }
}

export default async function RunsPage({
  searchParams,
}: {
  searchParams: Promise<{ pipeline?: string }>;
}) {
  const t = await getTranslations("runs");
  const params = await searchParams;
  const pipelineFilter = params.pipeline ?? undefined;

  const allRuns = listRuns() as RunMeta[];
  const filteredRuns = pipelineFilter
    ? (listRuns(pipelineFilter) as RunMeta[])
    : allRuns;

  // KPI calculations
  const totalRuns = allRuns.length;
  const completedRuns = allRuns.filter((r) => r.status === "completed").length;
  const failedRuns = allRuns.filter((r) => r.status === "failed").length;
  const pipelinesUsed = new Set(allRuns.map((r) => r.pipeline)).size;

  const statusColor = (status: string): string => {
    switch (status) {
      case "completed":
        return "var(--app-status-success, #16a34a)";
      case "failed":
        return "var(--app-status-error, #dc2626)";
      case "running":
        return "var(--app-status-running, #2563eb)";
      case "cancelled":
        return "var(--app-status-cancelled, #9ca3af)";
      default:
        return "var(--app-text-3)";
    }
  };

  const statusBg = (status: string): string => {
    switch (status) {
      case "completed":
        return "var(--app-status-success-bg, #f0fdf4)";
      case "failed":
        return "var(--app-status-error-bg, #fef2f2)";
      case "running":
        return "var(--app-status-running-bg, #eff6ff)";
      case "cancelled":
        return "var(--app-status-cancelled-bg, #f9fafb)";
      default:
        return "var(--app-surface-2)";
    }
  };

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span
          className="rounded-lg px-2.5 py-1 text-xs font-bold"
          style={{
            backgroundColor: "var(--app-surface-2)",
            color: "var(--app-text-2)",
          }}
        >
          runs
        </span>
        <div>
          <h1
            className="text-2xl font-bold"
            style={{ color: "var(--app-text)" }}
          >
            {t("title")}
          </h1>
          <p className="text-sm" style={{ color: "var(--app-text-3)" }}>
            {t("subtitle")}
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title={t("kpi.totalRuns")}
          value={totalRuns}
          subtitle={t("subtitle")}
          accentColor="bg-blue-500"
        />
        <KpiCard
          title={t("kpi.completed")}
          value={completedRuns}
          subtitle={t("kpi.completedSub")}
          accentColor="bg-emerald-500"
        />
        <KpiCard
          title={t("kpi.failed")}
          value={failedRuns}
          subtitle={t("kpi.failedSub")}
          accentColor="bg-red-500"
        />
        <KpiCard
          title={t("kpi.pipelinesUsed")}
          value={pipelinesUsed}
          subtitle={t("kpi.pipelinesUsedSub")}
          accentColor="bg-purple-500"
        />
      </div>

      {/* Table card */}
      <div
        className="rounded-xl overflow-hidden"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        {/* Table header with filter */}
        <div
          className="flex items-center justify-between px-5 py-4"
          style={{ borderBottom: "1px solid var(--app-border)" }}
        >
          <div>
            <h2
              className="text-sm font-semibold"
              style={{ color: "var(--app-text)" }}
            >
              {t("table.title")}
            </h2>
            <p
              className="text-xs mt-0.5"
              style={{ color: "var(--app-text-3)" }}
            >
              {t("table.subtitle", { count: filteredRuns.length })}
            </p>
          </div>
          <PipelineFilter current={pipelineFilter} />
        </div>

        {filteredRuns.length > 0 ? (
          <div className="overflow-x-auto">
            <table
              className="w-full text-left text-sm"
              data-testid="runs-table"
            >
              <thead
                className="text-xs uppercase tracking-wider"
                style={{
                  backgroundColor: "var(--app-surface-2)",
                  color: "var(--app-text-3)",
                  borderBottom: "1px solid var(--app-border)",
                }}
              >
                <tr>
                  <th className="px-4 py-3">{t("table.runId")}</th>
                  <th className="px-4 py-3">{t("table.pipeline")}</th>
                  <th className="px-4 py-3">{t("table.status")}</th>
                  <th className="px-4 py-3 text-right">
                    {t("table.duration")}
                  </th>
                  <th className="px-4 py-3">{t("table.date")}</th>
                  <th className="px-4 py-3">{t("table.tags")}</th>
                </tr>
              </thead>
              <tbody>
                {filteredRuns.map((run) => {
                  const route = PIPELINE_ROUTES[run.pipeline] ?? "/";
                  return (
                    <tr
                      key={run.run_id}
                      className="transition-colors"
                      style={{
                        borderBottom: "1px solid var(--app-border)",
                      }}
                    >
                      <td className="px-4 py-2.5">
                        <Link
                          href={`${route}?run=${run.run_id}`}
                          className="font-mono text-xs font-semibold hover:underline"
                          style={{ color: "var(--app-text)" }}
                          title={run.run_id}
                          data-testid="run-link"
                        >
                          {truncateRunId(run.run_id)}
                        </Link>
                      </td>
                      <td className="px-4 py-2.5">
                        <span
                          className="inline-flex items-center rounded-md px-2 py-0.5 text-xs font-medium capitalize"
                          style={{
                            backgroundColor: "var(--app-surface-2)",
                            color: "var(--app-text-2)",
                          }}
                        >
                          {t(`pipelines.${run.pipeline}`)}
                        </span>
                      </td>
                      <td className="px-4 py-2.5">
                        <span
                          className="inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-semibold"
                          style={{
                            backgroundColor: statusBg(run.status),
                            color: statusColor(run.status),
                          }}
                        >
                          {t(`status.${run.status}`)}
                        </span>
                      </td>
                      <td
                        className="px-4 py-2.5 text-right font-mono text-xs"
                        style={{ color: "var(--app-text-2)" }}
                      >
                        {formatDuration(run.total_duration_s)}
                      </td>
                      <td
                        className="px-4 py-2.5 text-xs"
                        style={{ color: "var(--app-text-3)" }}
                      >
                        {formatDate(run.created_at)}
                      </td>
                      <td className="px-4 py-2.5">
                        <div className="flex flex-wrap gap-1">
                          {run.tags.map((tag) => (
                            <span
                              key={tag}
                              className="inline-flex rounded-md px-1.5 py-0.5 text-xs"
                              style={{
                                backgroundColor: "var(--app-surface-2)",
                                color: "var(--app-text-3)",
                              }}
                            >
                              {tag}
                            </span>
                          ))}
                        </div>
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        ) : (
          /* Empty state */
          <div className="px-5 py-16 text-center">
            <svg
              viewBox="0 0 24 24"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              className="mx-auto mb-4 h-12 w-12"
              style={{ color: "var(--app-text-3)" }}
            >
              <path d="M3 7v10a2 2 0 002 2h14a2 2 0 002-2V9a2 2 0 00-2-2h-6l-2-2H5a2 2 0 00-2 2z" />
            </svg>
            <p
              className="text-sm font-semibold"
              style={{ color: "var(--app-text-2)" }}
            >
              {t("empty.title")}
            </p>
            <p
              className="mt-1 text-xs"
              style={{ color: "var(--app-text-3)" }}
            >
              {t("empty.subtitle")}
            </p>
            <div
              className="mx-auto mt-4 max-w-md rounded-lg p-3"
              style={{
                backgroundColor: "var(--app-surface-2)",
                border: "1px solid var(--app-border)",
              }}
            >
              <code
                className="text-xs font-mono"
                style={{ color: "var(--app-text-2)" }}
              >
                python -m marley run vaccine --organism leishmania_infantum
              </code>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
