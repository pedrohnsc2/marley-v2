"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

interface RunBannerProps {
  runId: string | null;
  pipeline: string;
  completedAt?: string | null;
}

/**
 * Thin banner displayed at the top of pipeline pages when viewing
 * a specific run instead of the latest data. Shows run ID, completion
 * date, and links to switch back to latest or compare runs.
 *
 * Renders nothing when runId is null (viewing latest data).
 */
export default function RunBanner({ runId, pipeline, completedAt }: RunBannerProps) {
  const pathname = usePathname();

  if (!runId) return null;

  // Strip the ?run= param to get back to the latest view
  const latestHref = pathname;
  const compareHref = `${pathname}?compare=${runId}`;

  const formattedDate = completedAt
    ? new Date(completedAt).toLocaleDateString(undefined, {
        year: "numeric",
        month: "short",
        day: "numeric",
        hour: "2-digit",
        minute: "2-digit",
      })
    : null;

  return (
    <div
      className="mb-4 flex flex-wrap items-center gap-3 rounded-lg px-4 py-2.5 text-sm"
      style={{
        backgroundColor: "var(--app-surface)",
        border: "1px solid var(--app-border)",
        color: "var(--app-text-2)",
      }}
      data-testid="run-banner"
    >
      <span
        className="rounded-md px-2 py-0.5 text-xs font-bold uppercase tracking-wider"
        style={{
          backgroundColor: "color-mix(in srgb, var(--app-border) 60%, transparent)",
          color: "var(--app-text)",
        }}
      >
        {pipeline}
      </span>

      <span style={{ color: "var(--app-text)" }}>
        Viewing run{" "}
        <span className="font-mono font-semibold">{runId}</span>
      </span>

      {formattedDate && (
        <span style={{ color: "var(--app-text-3)" }}>
          (completed {formattedDate})
        </span>
      )}

      <span className="flex-1" />

      <Link
        href={latestHref}
        className="rounded-md px-3 py-1 text-xs font-semibold transition-opacity hover:opacity-80"
        style={{
          backgroundColor: "color-mix(in srgb, var(--app-border) 80%, transparent)",
          color: "var(--app-text)",
        }}
        data-testid="run-banner-view-latest"
      >
        View Latest
      </Link>

      <Link
        href={compareHref}
        className="rounded-md px-3 py-1 text-xs font-semibold transition-opacity hover:opacity-80"
        style={{
          backgroundColor: "color-mix(in srgb, var(--app-border) 80%, transparent)",
          color: "var(--app-text)",
        }}
        data-testid="run-banner-compare"
      >
        Compare
      </Link>
    </div>
  );
}
