"use client";

import { useState, useEffect, useCallback } from "react";
import { PIPELINE_METADATA } from "@/lib/pipeline-metadata";
import type { PresetMeta } from "@/lib/presets";

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

interface PresetListEntry {
  name: string;
  meta: PresetMeta | null;
}

/* ------------------------------------------------------------------ */
/*  Props                                                              */
/* ------------------------------------------------------------------ */

interface StepPresetSelectProps {
  pipeline: string;
  selectedPreset: string | null;
  onSelect: (presetName: string, presetParams: Record<string, unknown>, displayName: string | null) => void;
  onBack: () => void;
  onNext: () => void;
}

/* ------------------------------------------------------------------ */
/*  Step component                                                     */
/* ------------------------------------------------------------------ */

export function StepPresetSelect({
  pipeline,
  selectedPreset,
  onSelect,
  onBack,
  onNext,
}: StepPresetSelectProps) {
  const [presets, setPresets] = useState<PresetListEntry[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [loadingPreset, setLoadingPreset] = useState<string | null>(null);

  const pipelineMeta = PIPELINE_METADATA[pipeline];

  // Fetch preset list for the selected pipeline
  useEffect(() => {
    let cancelled = false;
    setIsLoading(true);
    setError(null);

    fetch(`/api/presets/${pipeline}`)
      .then((res) => {
        if (!res.ok) throw new Error(`Failed to load presets (${res.status})`);
        return res.json();
      })
      .then((data: PresetListEntry[]) => {
        if (!cancelled) {
          setPresets(data);
          setIsLoading(false);
        }
      })
      .catch((err: Error) => {
        if (!cancelled) {
          setError(err.message);
          setIsLoading(false);
        }
      });

    return () => {
      cancelled = true;
    };
  }, [pipeline]);

  // Load full preset parameters when a preset is clicked
  const handlePresetClick = useCallback(
    async (presetName: string) => {
      setLoadingPreset(presetName);
      setError(null);

      try {
        const res = await fetch(
          `/api/presets/${pipeline}?name=${presetName}`,
        );
        if (!res.ok) {
          throw new Error(`Failed to load preset details (${res.status})`);
        }
        const data = (await res.json()) as {
          meta: PresetMeta | null;
          params: Record<string, unknown>;
        };
        onSelect(presetName, data.params, data.meta?.display_name ?? null);
      } catch (err) {
        const message =
          err instanceof Error ? err.message : "Failed to load preset";
        setError(message);
      } finally {
        setLoadingPreset(null);
      }
    },
    [pipeline, onSelect],
  );

  // Humanize preset name: use _meta from API, fallback to formatted filename
  function getPresetLabel(entry: PresetListEntry): {
    title: string;
    description: string;
    recommended: boolean;
  } {
    if (entry.meta) {
      return {
        title: entry.meta.display_name,
        description: entry.meta.description,
        recommended: entry.meta.recommended === true,
      };
    }
    // Fallback: capitalize and replace underscores
    return {
      title: entry.name
        .split("_")
        .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
        .join(" "),
      description: `Configuration preset for ${entry.name.replace(/_/g, " ")}`,
      recommended: false,
    };
  }

  return (
    <div data-testid="step-preset-select">
      {/* Header */}
      <div className="mb-6">
        <h2
          className="text-lg font-bold"
          style={{ color: "var(--app-text)" }}
        >
          Choose a Configuration
        </h2>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
          Select a preset for{" "}
          <span className="font-semibold" style={{ color: "var(--app-text)" }}>
            {pipelineMeta?.icon} {pipelineMeta?.displayName}
          </span>
          . You can customize parameters in the next step.
        </p>
      </div>

      {/* Loading state */}
      {isLoading && (
        <div
          className="flex items-center justify-center rounded-xl p-12"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
          }}
          data-testid="preset-loading"
        >
          <div className="flex items-center gap-3">
            <svg
              className="h-5 w-5 animate-spin"
              viewBox="0 0 24 24"
              fill="none"
              style={{ color: "var(--app-text-3)" }}
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
            <span
              className="text-sm"
              style={{ color: "var(--app-text-3)" }}
            >
              Loading presets...
            </span>
          </div>
        </div>
      )}

      {/* Error state */}
      {error && (
        <div
          className="mb-4 rounded-lg p-4 text-sm"
          style={{
            backgroundColor: "rgb(244 63 94 / 0.1)",
            color: "rgb(244 63 94)",
            border: "1px solid rgb(244 63 94 / 0.2)",
          }}
          data-testid="preset-error"
        >
          {error}
        </div>
      )}

      {/* Empty state */}
      {!isLoading && !error && presets.length === 0 && (
        <div
          className="rounded-xl p-8 text-center"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
          }}
          data-testid="preset-empty"
        >
          <p
            className="text-sm"
            style={{ color: "var(--app-text-3)" }}
          >
            No presets found for this pipeline. You can still proceed with
            default parameters.
          </p>
        </div>
      )}

      {/* Preset cards */}
      {!isLoading && presets.length > 0 && (
        <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
          {presets.map((entry) => {
            const label = getPresetLabel(entry);
            const name = entry.name;
            const isSelected = selectedPreset === name;
            const isLoadingThis = loadingPreset === name;

            return (
              <button
                key={name}
                type="button"
                onClick={() => handlePresetClick(name)}
                disabled={isLoadingThis}
                data-testid={`preset-card-${name}`}
                className="w-full rounded-xl p-5 text-left transition-all"
                style={{
                  backgroundColor: "var(--app-surface)",
                  border: isSelected
                    ? "2px solid var(--app-primary)"
                    : "1px solid var(--app-border)",
                  boxShadow: isSelected
                    ? "0 0 0 3px color-mix(in srgb, var(--app-primary) 20%, transparent)"
                    : "var(--app-card-shadow)",
                  opacity: isLoadingThis ? 0.7 : 1,
                }}
              >
                <div className="flex items-start justify-between gap-2">
                  <div>
                    <h3
                      className="text-sm font-bold"
                      style={{ color: "var(--app-text)" }}
                    >
                      {label.title}
                    </h3>
                    <p
                      className="mt-1 text-xs leading-relaxed"
                      style={{ color: "var(--app-text-2)" }}
                    >
                      {label.description}
                    </p>
                  </div>

                  {isSelected && (
                    <div
                      className="flex h-5 w-5 flex-shrink-0 items-center justify-center rounded-full"
                      style={{ backgroundColor: "var(--app-primary)" }}
                    >
                      <svg
                        viewBox="0 0 12 12"
                        className="h-3 w-3"
                        fill="none"
                        stroke="currentColor"
                        strokeWidth="2"
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        style={{ color: "var(--app-primary-tx)" }}
                      >
                        <path d="M2.5 6.5L5 9L9.5 3.5" />
                      </svg>
                    </div>
                  )}

                  {isLoadingThis && (
                    <svg
                      className="h-5 w-5 flex-shrink-0 animate-spin"
                      viewBox="0 0 24 24"
                      fill="none"
                      style={{ color: "var(--app-text-3)" }}
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
                  )}
                </div>

                {/* Badges */}
                <div className="mt-3 flex items-center gap-2">
                  {label.recommended && (
                    <span
                      className="rounded-md px-2 py-0.5 text-xs font-medium"
                      style={{
                        backgroundColor: "color-mix(in srgb, var(--app-primary) 15%, transparent)",
                        color: "var(--app-primary)",
                      }}
                    >
                      Recommended
                    </span>
                  )}
                  <span
                    className="rounded-md px-2 py-0.5 font-mono text-xs"
                    style={{
                      backgroundColor: "var(--app-surface-2)",
                      color: "var(--app-text-3)",
                    }}
                  >
                    {name}.json
                  </span>
                </div>
              </button>
            );
          })}
        </div>
      )}

      {/* Navigation */}
      <div className="mt-8 flex items-center justify-between">
        <button
          type="button"
          onClick={onBack}
          data-testid="wizard-back-step-2"
          className="rounded-lg px-6 py-2.5 text-sm font-semibold transition-opacity"
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
          onClick={onNext}
          disabled={!selectedPreset}
          data-testid="wizard-next-step-2"
          className="rounded-lg px-6 py-2.5 text-sm font-semibold transition-opacity disabled:cursor-not-allowed disabled:opacity-40"
          style={{
            backgroundColor: "var(--app-primary)",
            color: "var(--app-primary-tx)",
          }}
        >
          Continue
        </button>
      </div>
    </div>
  );
}
