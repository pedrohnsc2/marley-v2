import type { RunMetadata } from "@/lib/types/run";
import RunStatusBadge from "@/components/runs/run-status-badge";

interface RunMetadataCardProps {
  run: RunMetadata;
}

function formatDuration(seconds: number): string {
  if (seconds < 60) return `${seconds.toFixed(1)}s`;
  const min = Math.floor(seconds / 60);
  const sec = (seconds % 60).toFixed(0);
  return `${min}m ${sec}s`;
}

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

export default function RunMetadataCard({ run }: RunMetadataCardProps) {
  const paramCount = Object.keys(run.parameters).length;
  const stageCount = run.stages.length;
  const successCount = run.stages.filter((s) => s.status === "success").length;

  return (
    <div
      className="rounded-xl p-5"
      style={{
        backgroundColor: "var(--app-surface)",
        border: "1px solid var(--app-border)",
        boxShadow: "var(--app-card-shadow)",
      }}
      data-testid="run-metadata-card"
    >
      {/* Header row */}
      <div className="flex items-start justify-between gap-3">
        <div className="min-w-0 flex-1">
          <p
            className="text-xs font-semibold uppercase tracking-wider"
            style={{ color: "var(--app-text-2)" }}
          >
            Pipeline Run
          </p>
          <p
            className="mt-1 truncate text-sm font-mono font-medium"
            style={{ color: "var(--app-text)" }}
            title={run.run_id}
          >
            {run.run_id}
          </p>
        </div>
        <RunStatusBadge status={run.status} />
      </div>

      {/* Details grid */}
      <div className="mt-4 grid grid-cols-2 gap-x-6 gap-y-2">
        <Detail label="Pipeline" value={run.pipeline} />
        <Detail label="Created" value={formatDate(run.created_at)} />
        <Detail label="Duration" value={formatDuration(run.total_duration_s)} />
        <Detail label="Git SHA" value={run.git_sha || "--"} mono />
        <Detail label="Parameters" value={String(paramCount)} />
        <Detail label="Stages" value={`${successCount}/${stageCount}`} />
      </div>

      {/* Tags */}
      {run.tags.length > 0 && (
        <div className="mt-3 flex flex-wrap gap-1.5">
          {run.tags.map((tag) => (
            <span
              key={tag}
              className="rounded-md px-2 py-0.5 text-xs font-medium"
              style={{
                backgroundColor: "var(--app-surface-2)",
                color: "var(--app-text-2)",
              }}
            >
              {tag}
            </span>
          ))}
        </div>
      )}

      {/* Accent bar matching KpiCard */}
      <div
        className="mt-4 h-1 w-full rounded-full"
        style={{ backgroundColor: "var(--app-accent-bar)", opacity: 0.7 }}
      />
    </div>
  );
}

function Detail({
  label,
  value,
  mono = false,
}: {
  label: string;
  value: string;
  mono?: boolean;
}) {
  return (
    <div className="min-w-0">
      <p
        className="text-xs"
        style={{ color: "var(--app-text-3)" }}
      >
        {label}
      </p>
      <p
        className={`truncate text-sm font-medium ${mono ? "font-mono" : ""}`}
        style={{ color: "var(--app-text)" }}
        title={value}
      >
        {value}
      </p>
    </div>
  );
}
