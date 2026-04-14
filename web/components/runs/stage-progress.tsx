"use client";

import { useTranslations } from "next-intl";
import type { StageRecord, StageStatus } from "@/lib/types/run";
import FriendlyError from "@/components/runs/friendly-error";

interface StageProgressProps {
  stages: StageRecord[];
  totalExpected?: number;
  runId?: string;
}

function formatDuration(seconds: number): string {
  if (seconds < 60) return `${seconds.toFixed(1)}s`;
  const min = Math.floor(seconds / 60);
  const sec = (seconds % 60).toFixed(0);
  return `${min}m ${sec}s`;
}

function StatusIndicator({ status }: { status: StageStatus }) {
  switch (status) {
    case "running":
      return (
        <span
          className="relative flex h-4 w-4 items-center justify-center"
          data-testid="stage-indicator-running"
        >
          <span
            className="absolute inline-flex h-full w-full animate-ping rounded-full opacity-75"
            style={{ backgroundColor: "var(--app-accent-bar)" }}
          />
          <span
            className="relative inline-flex h-2.5 w-2.5 rounded-full"
            style={{ backgroundColor: "var(--app-accent-bar)" }}
          />
        </span>
      );

    case "success":
      return (
        <span
          className="flex h-4 w-4 items-center justify-center rounded-full"
          style={{ backgroundColor: "rgb(16 185 129)" }}
          data-testid="stage-indicator-success"
        >
          <svg
            viewBox="0 0 12 12"
            className="h-2.5 w-2.5"
            fill="none"
            stroke="white"
            strokeWidth="2"
            strokeLinecap="round"
            strokeLinejoin="round"
          >
            <path d="M2.5 6.5L5 9L9.5 3.5" />
          </svg>
        </span>
      );

    case "failed":
      return (
        <span
          className="flex h-4 w-4 items-center justify-center rounded-full"
          style={{ backgroundColor: "rgb(244 63 94)" }}
          data-testid="stage-indicator-failed"
        >
          <svg
            viewBox="0 0 12 12"
            className="h-2.5 w-2.5"
            fill="none"
            stroke="white"
            strokeWidth="2"
            strokeLinecap="round"
            strokeLinejoin="round"
          >
            <path d="M3 3L9 9M9 3L3 9" />
          </svg>
        </span>
      );

    case "skipped":
      return (
        <span
          className="flex h-4 w-4 items-center justify-center rounded-full border-2 border-dashed"
          style={{ borderColor: "var(--app-text-3)" }}
          data-testid="stage-indicator-skipped"
        />
      );

    case "pending":
    default:
      return (
        <span
          className="flex h-4 w-4 items-center justify-center rounded-full border-2"
          style={{ borderColor: "var(--app-text-3)" }}
          data-testid="stage-indicator-pending"
        />
      );
  }
}

export default function StageProgress({
  stages,
  totalExpected,
  runId = "",
}: StageProgressProps) {
  const t = useTranslations("runs");

  const completedCount = stages.filter(
    (s) => s.status === "success" || s.status === "failed" || s.status === "skipped",
  ).length;
  const total = totalExpected ?? stages.length;

  return (
    <div data-testid="stage-progress">
      {/* Overall progress bar */}
      <div className="mb-4">
        <div className="mb-1.5 flex items-center justify-between">
          <p
            className="text-xs font-semibold uppercase tracking-wider"
            style={{ color: "var(--app-text-2)" }}
          >
            {t("detail.stages")}
          </p>
          <p className="text-xs" style={{ color: "var(--app-text-3)" }}>
            {t("detail.overallProgress", {
              completed: completedCount,
              total,
            })}
          </p>
        </div>
        <div
          className="h-2 w-full overflow-hidden rounded-full"
          style={{ backgroundColor: "var(--app-surface-2)" }}
        >
          <div
            className="h-full rounded-full transition-all duration-500"
            style={{
              width: total > 0 ? `${(completedCount / total) * 100}%` : "0%",
              backgroundColor: "var(--app-accent-bar)",
            }}
          />
        </div>
      </div>

      {/* Stage list */}
      <div className="relative">
        {stages.map((stage, index) => {
          const isLast = index === stages.length - 1;
          const hasError = stage.status === "failed" && stage.error;

          return (
            <div key={stage.stage_id} className="relative" data-testid={`stage-row-${stage.stage_id}`}>
              {/* Connecting line */}
              {!isLast && (
                <div
                  className="absolute left-2 top-5 w-px"
                  style={{
                    backgroundColor: "var(--app-border)",
                    height: "calc(100% - 4px)",
                  }}
                />
              )}

              {/* Stage row */}
              <div className="flex items-start gap-3 pb-4">
                {/* Status indicator */}
                <div className="flex-shrink-0 pt-0.5">
                  <StatusIndicator status={stage.status} />
                </div>

                {/* Stage info */}
                <div className="min-w-0 flex-1">
                  <div className="flex items-center justify-between gap-2">
                    <p
                      className="truncate text-sm font-medium"
                      style={{ color: "var(--app-text)" }}
                    >
                      {stage.name}
                    </p>
                    <div className="flex items-center gap-2">
                      <span
                        className="whitespace-nowrap text-xs"
                        style={{
                          color:
                            stage.status === "running"
                              ? "var(--app-accent-bar)"
                              : "var(--app-text-3)",
                        }}
                      >
                        {t(`detail.stageStatus.${stage.status}`)}
                      </span>
                      {stage.duration_s > 0 && (
                        <span
                          className="whitespace-nowrap text-xs font-mono"
                          style={{ color: "var(--app-text-3)" }}
                        >
                          {formatDuration(stage.duration_s)}
                        </span>
                      )}
                    </div>
                  </div>

                  {/* Friendly error component for failed stages */}
                  {hasError && (
                    <FriendlyError
                      errorInfo={stage.error_info ?? null}
                      rawError={stage.error}
                      runId={runId}
                      stageId={stage.stage_id}
                    />
                  )}
                </div>
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
}
