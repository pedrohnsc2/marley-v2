"use client";

import { useRouter, useSearchParams } from "next/navigation";
import { useCallback } from "react";
import type { RunMetadata } from "@/lib/types/run";
import RunStatusBadge from "@/components/runs/run-status-badge";

interface RunSelectorProps {
  pipeline: string;
  runs: RunMetadata[];
  currentRunId?: string;
}

function formatDate(iso: string): string {
  try {
    const d = new Date(iso);
    return d.toLocaleDateString(undefined, {
      month: "short",
      day: "numeric",
      hour: "2-digit",
      minute: "2-digit",
    });
  } catch {
    return iso;
  }
}

function truncateId(runId: string): string {
  if (runId.length <= 20) return runId;
  return runId.slice(0, 17) + "...";
}

export default function RunSelector({ pipeline, runs, currentRunId }: RunSelectorProps) {
  const router = useRouter();
  const searchParams = useSearchParams();

  const handleChange = useCallback(
    (e: React.ChangeEvent<HTMLSelectElement>) => {
      const value = e.target.value;
      const params = new URLSearchParams(searchParams.toString());

      if (value === "") {
        params.delete("run");
      } else {
        params.set("run", value);
      }

      const qs = params.toString();
      router.push(qs ? `?${qs}` : window.location.pathname);
    },
    [router, searchParams],
  );

  if (runs.length === 0) {
    return null;
  }

  return (
    <div className="flex items-center gap-3" data-testid="run-selector">
      <label
        htmlFor={`run-select-${pipeline}`}
        className="text-xs font-semibold uppercase tracking-wider"
        style={{ color: "var(--app-text-2)" }}
      >
        Run
      </label>

      <div className="relative">
        <select
          id={`run-select-${pipeline}`}
          value={currentRunId ?? ""}
          onChange={handleChange}
          className="appearance-none rounded-lg px-3 py-1.5 pr-8 text-sm font-medium outline-none"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
            color: "var(--app-text)",
          }}
        >
          <option value="">Latest</option>
          {runs.map((run) => (
            <option key={run.run_id} value={run.run_id}>
              {formatDate(run.created_at)} - {truncateId(run.run_id)} [{run.status}]
            </option>
          ))}
        </select>

        {/* Dropdown chevron */}
        <svg
          className="pointer-events-none absolute right-2 top-1/2 h-4 w-4 -translate-y-1/2"
          style={{ color: "var(--app-text-3)" }}
          viewBox="0 0 24 24"
          fill="none"
          stroke="currentColor"
          strokeWidth={2}
        >
          <path d="M6 9l6 6 6-6" />
        </svg>
      </div>

      {/* Show badge for current selection */}
      {currentRunId && runs.length > 0 && (() => {
        const selected = runs.find((r) => r.run_id === currentRunId);
        return selected ? <RunStatusBadge status={selected.status} /> : null;
      })()}
    </div>
  );
}
