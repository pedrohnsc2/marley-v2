"use client";

import { useState } from "react";
import { useParams } from "next/navigation";
import { useTranslations } from "next-intl";
import { Link } from "@/i18n/routing";
import { useRunRealtime } from "@/lib/hooks/use-run-realtime";
import RunStatusBadge from "@/components/runs/run-status-badge";
import StageProgress from "@/components/runs/stage-progress";

const PIPELINE_ROUTES: Record<string, string> = {
  vaccine: "/vaccine",
  drug: "/drug",
  docking: "/docking",
  rna: "/rna",
  aso_math: "/aso",
};

function formatDate(iso: string | null): string {
  if (!iso) return "--";
  try {
    return new Date(iso).toLocaleString(undefined, {
      month: "short",
      day: "numeric",
      year: "numeric",
      hour: "2-digit",
      minute: "2-digit",
    });
  } catch {
    return iso;
  }
}

export default function RunDetailPage() {
  const t = useTranslations("runs");
  const params = useParams();
  const runId = params.runId as string;

  const { run, stages, isLive, error } = useRunRealtime(runId);
  const [paramsOpen, setParamsOpen] = useState(false);

  if (error && !run) {
    return (
      <div className="space-y-4" data-testid="run-detail-error">
        <Link
          href="/runs"
          className="text-sm font-medium underline"
          style={{ color: "var(--app-text-2)" }}
        >
          {t("detail.backToRuns")}
        </Link>
        <div
          className="rounded-xl p-5"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
          }}
        >
          <p className="text-sm" style={{ color: "rgb(244 63 94)" }}>
            {error}
          </p>
        </div>
      </div>
    );
  }

  if (!run) {
    return (
      <div className="flex items-center justify-center py-20" data-testid="run-detail-loading">
        <div
          className="h-6 w-6 animate-spin rounded-full border-2 border-t-transparent"
          style={{ borderColor: "var(--app-border)", borderTopColor: "transparent" }}
        />
      </div>
    );
  }

  const pipelineRoute = PIPELINE_ROUTES[run.pipeline];
  const paramEntries = Object.entries(run.parameters);

  return (
    <div className="space-y-6" data-testid="run-detail-page">
      {/* Back link */}
      <Link
        href="/runs"
        className="inline-flex items-center gap-1 text-sm font-medium"
        style={{ color: "var(--app-text-2)" }}
        data-testid="back-to-runs"
      >
        <svg
          viewBox="0 0 16 16"
          className="h-4 w-4"
          fill="currentColor"
        >
          <path d="M10.354 3.354a.5.5 0 0 0-.708-.708l-5 5a.5.5 0 0 0 0 .708l5 5a.5.5 0 0 0 .708-.708L5.707 8l4.647-4.646z" />
        </svg>
        {t("detail.backToRuns")}
      </Link>

      {/* Header */}
      <div
        className="rounded-xl p-5"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
        data-testid="run-detail-header"
      >
        <div className="flex flex-wrap items-start justify-between gap-3">
          <div className="min-w-0 flex-1">
            <div className="flex items-center gap-3">
              <p
                className="text-xs font-semibold uppercase tracking-wider"
                style={{ color: "var(--app-text-2)" }}
              >
                {t("detail.title")}
              </p>
              {isLive && (
                <span
                  className="inline-flex items-center gap-1.5 rounded-full px-2 py-0.5 text-xs font-medium"
                  style={{
                    backgroundColor: "rgb(16 185 129 / 0.15)",
                    color: "rgb(16 185 129)",
                  }}
                  data-testid="live-indicator"
                >
                  <span className="relative flex h-2 w-2">
                    <span
                      className="absolute inline-flex h-full w-full animate-ping rounded-full opacity-75"
                      style={{ backgroundColor: "rgb(16 185 129)" }}
                    />
                    <span
                      className="relative inline-flex h-2 w-2 rounded-full"
                      style={{ backgroundColor: "rgb(16 185 129)" }}
                    />
                  </span>
                  {t("detail.live")}
                </span>
              )}
            </div>
            <p
              className="mt-1 truncate font-mono text-sm font-medium"
              style={{ color: "var(--app-text)" }}
              title={run.run_id}
            >
              {run.run_id}
            </p>
          </div>
          <RunStatusBadge status={run.status} />
        </div>

        <div className="mt-4 grid grid-cols-2 gap-x-6 gap-y-2 sm:grid-cols-4">
          <DetailField label="Pipeline" value={run.pipeline} />
          <DetailField label="Created" value={formatDate(run.created_at)} />
          <DetailField label="Started" value={formatDate(run.started_at)} />
          <DetailField label="Completed" value={formatDate(run.completed_at)} />
        </div>

        <div
          className="mt-4 h-1 w-full rounded-full"
          style={{ backgroundColor: "var(--app-accent-bar)", opacity: 0.7 }}
        />
      </div>

      {/* Success banner */}
      {run.status === "completed" && (
        <div
          className="flex flex-wrap items-center justify-between gap-3 rounded-xl p-4"
          style={{
            backgroundColor: "rgb(16 185 129 / 0.1)",
            border: "1px solid rgb(16 185 129 / 0.3)",
          }}
          data-testid="run-completed-banner"
        >
          <p className="text-sm font-medium" style={{ color: "rgb(16 185 129)" }}>
            {t("detail.completed")}
          </p>
          {pipelineRoute && (
            <Link
              href={`${pipelineRoute}?run=${run.run_id}`}
              className="rounded-lg px-3 py-1.5 text-xs font-medium"
              style={{
                backgroundColor: "rgb(16 185 129)",
                color: "white",
              }}
              data-testid="view-results-link"
            >
              {t("detail.viewResults")}
            </Link>
          )}
        </div>
      )}

      {/* Failed banner */}
      {run.status === "failed" && (
        <div
          className="rounded-xl p-4"
          style={{
            backgroundColor: "rgb(244 63 94 / 0.1)",
            border: "1px solid rgb(244 63 94 / 0.3)",
          }}
          data-testid="run-failed-banner"
        >
          <p className="text-sm font-medium" style={{ color: "rgb(244 63 94)" }}>
            {t("detail.failed")}
          </p>
        </div>
      )}

      {/* Stalled warning */}
      {run.status === "running" && !isLive && (
        <div
          className="rounded-xl p-4"
          style={{
            backgroundColor: "rgb(245 158 11 / 0.1)",
            border: "1px solid rgb(245 158 11 / 0.3)",
          }}
          data-testid="run-stalled-warning"
        >
          <p className="text-sm font-medium" style={{ color: "rgb(245 158 11)" }}>
            {t("detail.stalled")}
          </p>
        </div>
      )}

      {/* Stage progress */}
      <div
        className="rounded-xl p-5"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
        data-testid="run-stages-section"
      >
        <StageProgress stages={stages} />
      </div>

      {/* Parameters (collapsible) */}
      {paramEntries.length > 0 && (
        <details
          open={paramsOpen}
          onToggle={(e) =>
            setParamsOpen((e.target as HTMLDetailsElement).open)
          }
          className="rounded-xl"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
          }}
          data-testid="run-parameters-section"
        >
          <summary
            className="cursor-pointer px-5 py-3 text-sm font-medium"
            style={{ color: "var(--app-text-2)" }}
          >
            {t("detail.parameters")}
          </summary>
          <div className="px-5 pb-4">
            <div className="grid grid-cols-1 gap-2 sm:grid-cols-2 lg:grid-cols-3">
              {paramEntries.map(([key, value]) => (
                <div key={key} className="min-w-0">
                  <p
                    className="text-xs"
                    style={{ color: "var(--app-text-3)" }}
                  >
                    {key.replace(/_/g, " ")}
                  </p>
                  <p
                    className="truncate text-sm font-mono font-medium"
                    style={{ color: "var(--app-text)" }}
                    title={String(value)}
                  >
                    {String(value)}
                  </p>
                </div>
              ))}
            </div>
          </div>
        </details>
      )}
    </div>
  );
}

function DetailField({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0">
      <p className="text-xs" style={{ color: "var(--app-text-3)" }}>
        {label}
      </p>
      <p
        className="truncate text-sm font-medium"
        style={{ color: "var(--app-text)" }}
        title={value}
      >
        {value}
      </p>
    </div>
  );
}
