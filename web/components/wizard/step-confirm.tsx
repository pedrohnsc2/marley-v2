"use client";

import { useState, useMemo } from "react";
import { useRouter } from "@/i18n/routing";
import {
  PIPELINE_METADATA,
  PARAMETER_METADATA,
} from "@/lib/pipeline-metadata";
import type { WizardState } from "./pipeline-wizard";

/* ------------------------------------------------------------------ */
/*  Props                                                              */
/* ------------------------------------------------------------------ */

interface StepConfirmProps {
  state: WizardState;
  onBack: () => void;
}

/* ------------------------------------------------------------------ */
/*  Step component                                                     */
/* ------------------------------------------------------------------ */

export function StepConfirm({ state, onBack }: StepConfirmProps) {
  const router = useRouter();
  const [isLaunching, setIsLaunching] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const pipelineMeta = PIPELINE_METADATA[state.pipeline!];
  const presetLabel = state.presetDisplayName
    ?? (state.preset
      ? state.preset
          .split("_")
          .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
          .join(" ")
      : "Custom");

  // Detect modified parameters
  const modifiedParams = useMemo(() => {
    const mods: { key: string; label: string; oldValue: unknown; newValue: unknown }[] = [];
    for (const key of Object.keys(state.parameters)) {
      const original = state.originalParameters[key];
      const current = state.parameters[key];
      if (JSON.stringify(original) !== JSON.stringify(current)) {
        const meta = PARAMETER_METADATA[key];
        mods.push({
          key,
          label: meta?.label ?? key.replace(/_/g, " "),
          oldValue: original,
          newValue: current,
        });
      }
    }
    return mods;
  }, [state.parameters, state.originalParameters]);

  const handleLaunch = async () => {
    setIsLaunching(true);
    setError(null);

    try {
      const res = await fetch("/api/runs/start", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          pipeline: state.pipeline,
          preset: state.preset || undefined,
          parameters: state.parameters,
        }),
      });

      if (!res.ok) {
        const body = await res.json().catch(() => ({}));
        throw new Error(
          (body as Record<string, string>).error ??
            `Request failed: ${res.status}`,
        );
      }

      const data = (await res.json()) as { run_id: string };
      router.push(`/runs/${data.run_id}`);
    } catch (err) {
      const message =
        err instanceof Error ? err.message : "Failed to start pipeline";
      setError(message);
    } finally {
      setIsLaunching(false);
    }
  };

  if (!pipelineMeta) {
    return (
      <div
        className="rounded-xl p-8 text-center"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
        }}
      >
        <p style={{ color: "var(--app-text-3)" }}>
          No pipeline selected. Please go back and select a pipeline.
        </p>
      </div>
    );
  }

  return (
    <div data-testid="step-confirm">
      {/* Header */}
      <div className="mb-6">
        <h2
          className="text-lg font-bold"
          style={{ color: "var(--app-text)" }}
        >
          Review and Launch
        </h2>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
          Review your configuration before starting the analysis.
        </p>
      </div>

      {/* Summary card */}
      <div
        className="mb-6 rounded-xl p-6"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
        data-testid="confirm-summary"
      >
        {/* Pipeline info */}
        <div className="mb-5 flex items-start gap-4">
          <span className="text-3xl" role="img" aria-hidden="true">
            {pipelineMeta.icon}
          </span>
          <div className="flex-1">
            <h3
              className="text-base font-bold"
              style={{ color: "var(--app-text)" }}
            >
              {pipelineMeta.displayName}
            </h3>
            <p
              className="mt-0.5 text-sm"
              style={{ color: "var(--app-text-2)" }}
            >
              {pipelineMeta.description}
            </p>
          </div>
        </div>

        {/* Details grid */}
        <div
          className="grid grid-cols-1 gap-4 border-t pt-5 sm:grid-cols-3"
          style={{ borderColor: "var(--app-border)" }}
        >
          {/* Preset */}
          <div>
            <p
              className="text-xs font-semibold uppercase tracking-wider"
              style={{ color: "var(--app-text-3)" }}
            >
              Preset
            </p>
            <p
              className="mt-1 text-sm font-medium"
              style={{ color: "var(--app-text)" }}
            >
              {presetLabel}
            </p>
          </div>

          {/* Estimated duration */}
          <div>
            <p
              className="text-xs font-semibold uppercase tracking-wider"
              style={{ color: "var(--app-text-3)" }}
            >
              Estimated Duration
            </p>
            <p
              className="mt-1 text-sm font-medium"
              style={{ color: "var(--app-text)" }}
            >
              {pipelineMeta.estimatedDuration}
            </p>
          </div>

          {/* Stages */}
          <div>
            <p
              className="text-xs font-semibold uppercase tracking-wider"
              style={{ color: "var(--app-text-3)" }}
            >
              Stages
            </p>
            <p
              className="mt-1 text-sm font-medium"
              style={{ color: "var(--app-text)" }}
            >
              {pipelineMeta.stages.length} stages
            </p>
          </div>
        </div>

        {/* Expected outputs */}
        <div
          className="mt-5 border-t pt-5"
          style={{ borderColor: "var(--app-border)" }}
        >
          <p
            className="mb-2 text-xs font-semibold uppercase tracking-wider"
            style={{ color: "var(--app-text-3)" }}
          >
            Expected Outputs
          </p>
          <ul className="space-y-1.5">
            {pipelineMeta.expectedOutputs.map((output) => (
              <li
                key={output}
                className="flex items-start gap-2 text-sm"
                style={{ color: "var(--app-text-2)" }}
              >
                <svg
                  viewBox="0 0 12 12"
                  className="mt-0.5 h-3 w-3 flex-shrink-0"
                  fill="none"
                  stroke="currentColor"
                  strokeWidth="2"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  style={{ color: "var(--app-primary)" }}
                >
                  <path d="M2.5 6.5L5 9L9.5 3.5" />
                </svg>
                {output}
              </li>
            ))}
          </ul>
        </div>

        {/* Modified parameters */}
        {modifiedParams.length > 0 && (
          <div
            className="mt-5 border-t pt-5"
            style={{ borderColor: "var(--app-border)" }}
          >
            <p
              className="mb-2 text-xs font-semibold uppercase tracking-wider"
              style={{ color: "rgb(245 158 11)" }}
            >
              Modified Parameters ({modifiedParams.length})
            </p>
            <div className="space-y-2">
              {modifiedParams.map((mod) => (
                <div
                  key={mod.key}
                  className="flex items-center justify-between rounded-lg px-3 py-2 text-sm"
                  style={{
                    backgroundColor: "var(--app-surface-2)",
                    border: "1px solid var(--app-border)",
                  }}
                  data-testid={`modified-param-${mod.key}`}
                >
                  <span
                    className="font-medium"
                    style={{ color: "var(--app-text)" }}
                  >
                    {mod.label}
                  </span>
                  <div className="flex items-center gap-2 text-xs">
                    <span
                      className="line-through"
                      style={{ color: "var(--app-text-3)" }}
                    >
                      {formatValue(mod.oldValue)}
                    </span>
                    <svg
                      viewBox="0 0 24 24"
                      fill="none"
                      stroke="currentColor"
                      strokeWidth="2"
                      className="h-3 w-3"
                      style={{ color: "var(--app-text-3)" }}
                    >
                      <path d="M5 12h14M12 5l7 7-7 7" />
                    </svg>
                    <span
                      className="font-semibold"
                      style={{ color: "rgb(245 158 11)" }}
                    >
                      {formatValue(mod.newValue)}
                    </span>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Pipeline stages preview */}
        <div
          className="mt-5 border-t pt-5"
          style={{ borderColor: "var(--app-border)" }}
        >
          <p
            className="mb-3 text-xs font-semibold uppercase tracking-wider"
            style={{ color: "var(--app-text-3)" }}
          >
            Pipeline Stages
          </p>
          <div className="space-y-2">
            {pipelineMeta.stages.map((stage, idx) => (
              <div
                key={stage.stageId}
                className="flex items-center gap-3 text-sm"
              >
                <span
                  className="flex h-6 w-6 flex-shrink-0 items-center justify-center rounded-full text-xs font-bold"
                  style={{
                    backgroundColor: "var(--app-surface-2)",
                    color: "var(--app-text-3)",
                    border: "1px solid var(--app-border)",
                  }}
                >
                  {idx + 1}
                </span>
                <div>
                  <span
                    className="font-medium"
                    style={{ color: "var(--app-text)" }}
                  >
                    {stage.friendlyName}
                  </span>
                  <span
                    className="ml-2 text-xs"
                    style={{ color: "var(--app-text-3)" }}
                  >
                    {stage.description}
                  </span>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>

      {/* Error display */}
      {error && (
        <div
          className="mb-4 rounded-lg p-4 text-sm"
          style={{
            backgroundColor: "rgb(244 63 94 / 0.1)",
            color: "rgb(244 63 94)",
            border: "1px solid rgb(244 63 94 / 0.2)",
          }}
          data-testid="launch-error"
        >
          {error}
        </div>
      )}

      {/* Navigation */}
      <div className="flex items-center justify-between">
        <button
          type="button"
          onClick={onBack}
          disabled={isLaunching}
          data-testid="wizard-back-step-4"
          className="rounded-lg px-6 py-2.5 text-sm font-semibold transition-opacity disabled:opacity-40"
          style={{
            backgroundColor: "var(--app-surface-2)",
            color: "var(--app-text-2)",
            border: "1px solid var(--app-border)",
          }}
        >
          Back
        </button>

        <button
          type="button"
          onClick={handleLaunch}
          disabled={isLaunching}
          data-testid="launch-analysis-button"
          className="rounded-lg px-8 py-3 text-sm font-bold transition-opacity disabled:opacity-60"
          style={{
            backgroundColor: "var(--app-primary)",
            color: "var(--app-primary-tx)",
          }}
        >
          {isLaunching ? (
            <span className="inline-flex items-center gap-2">
              <svg
                className="h-4 w-4 animate-spin"
                viewBox="0 0 24 24"
                fill="none"
              >
                <circle
                  className="opacity-25"
                  cx="12"
                  cy="12"
                  r="10"
                  stroke="currentColor"
                  strokeWidth="4"
                />
                <path
                  className="opacity-75"
                  fill="currentColor"
                  d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4z"
                />
              </svg>
              Launching...
            </span>
          ) : (
            "Launch Analysis"
          )}
        </button>
      </div>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

function formatValue(value: unknown): string {
  if (Array.isArray(value)) {
    if (value.length <= 3) return value.join(", ");
    return `${value.slice(0, 2).join(", ")} +${value.length - 2} more`;
  }
  if (typeof value === "boolean") return value ? "true" : "false";
  if (value === null || value === undefined) return "--";
  return String(value);
}
