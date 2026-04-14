"use client";

import { useState, useMemo } from "react";
import {
  PARAMETER_METADATA,
  type ParameterMeta,
} from "@/lib/pipeline-metadata";

/* ------------------------------------------------------------------ */
/*  Props                                                              */
/* ------------------------------------------------------------------ */

interface StepParametersProps {
  pipeline: string;
  parameters: Record<string, unknown>;
  originalParameters: Record<string, unknown>;
  onUpdate: (key: string, value: unknown) => void;
  onBack: () => void;
  onNext: () => void;
}

/* ------------------------------------------------------------------ */
/*  Tooltip                                                            */
/* ------------------------------------------------------------------ */

function TooltipIcon({ text }: { text: string }) {
  const [visible, setVisible] = useState(false);

  return (
    <span className="relative inline-flex">
      <button
        type="button"
        className="flex h-4 w-4 items-center justify-center rounded-full text-xs"
        style={{
          backgroundColor: "var(--app-surface-2)",
          color: "var(--app-text-3)",
          border: "1px solid var(--app-border)",
        }}
        onMouseEnter={() => setVisible(true)}
        onMouseLeave={() => setVisible(false)}
        onFocus={() => setVisible(true)}
        onBlur={() => setVisible(false)}
        aria-label="More information"
      >
        ?
      </button>
      {visible && (
        <div
          className="absolute bottom-full left-1/2 z-50 mb-2 w-56 -translate-x-1/2 rounded-lg p-2.5 text-xs leading-relaxed shadow-lg"
          style={{
            backgroundColor: "var(--app-surface)",
            color: "var(--app-text-2)",
            border: "1px solid var(--app-border)",
          }}
          role="tooltip"
        >
          {text}
        </div>
      )}
    </span>
  );
}

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

function isModified(
  key: string,
  current: Record<string, unknown>,
  original: Record<string, unknown>,
): boolean {
  return JSON.stringify(current[key]) !== JSON.stringify(original[key]);
}

function inferInputType(value: unknown): "text" | "number" | "boolean" | "array" {
  if (typeof value === "boolean") return "boolean";
  if (typeof value === "number") return "number";
  if (Array.isArray(value)) return "array";
  return "text";
}

function getMeta(key: string): ParameterMeta {
  return (
    PARAMETER_METADATA[key] ?? {
      label: key
        .split("_")
        .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
        .join(" "),
      tooltip: `Parameter: ${key}`,
      group: "Other",
      advanced: false,
    }
  );
}

/* ------------------------------------------------------------------ */
/*  Parameter input component                                          */
/* ------------------------------------------------------------------ */

function ParameterInput({
  paramKey,
  value,
  modified,
  meta,
  onUpdate,
}: {
  paramKey: string;
  value: unknown;
  modified: boolean;
  meta: ParameterMeta;
  onUpdate: (key: string, value: unknown) => void;
}) {
  const inputType = inferInputType(value);

  return (
    <div data-testid={`param-field-${paramKey}`}>
      {/* Label row */}
      <div className="mb-1.5 flex items-center gap-2">
        <label
          htmlFor={`param-${paramKey}`}
          className="text-xs font-medium"
          style={{ color: "var(--app-text-2)" }}
        >
          {meta.label}
        </label>
        <TooltipIcon text={meta.tooltip} />
        {modified ? (
          <span
            className="rounded px-1.5 py-0.5 text-xs font-semibold"
            style={{
              backgroundColor: "rgb(245 158 11 / 0.15)",
              color: "rgb(245 158 11)",
            }}
          >
            modified
          </span>
        ) : (
          <span
            className="rounded px-1.5 py-0.5 text-xs"
            style={{
              backgroundColor: "var(--app-surface-2)",
              color: "var(--app-text-3)",
            }}
          >
            default
          </span>
        )}
      </div>

      {/* Input */}
      {inputType === "boolean" && (
        <button
          type="button"
          id={`param-${paramKey}`}
          onClick={() => onUpdate(paramKey, !value)}
          className="rounded-lg border px-3 py-2 text-sm font-medium"
          style={{
            backgroundColor: value
              ? "rgb(16 185 129 / 0.1)"
              : "var(--app-surface-2)",
            borderColor: value
              ? "rgb(16 185 129 / 0.3)"
              : "var(--app-border)",
            color: "var(--app-text)",
          }}
          data-testid={`param-toggle-${paramKey}`}
        >
          {value ? "true" : "false"}
        </button>
      )}

      {inputType === "number" && (
        <input
          type="number"
          id={`param-${paramKey}`}
          value={value as number}
          step="any"
          onChange={(e) =>
            onUpdate(paramKey, parseFloat(e.target.value) || 0)
          }
          className="w-full rounded-lg border px-3 py-2 text-sm"
          style={{
            backgroundColor: "var(--app-surface-2)",
            borderColor: "var(--app-border)",
            color: "var(--app-text)",
          }}
          data-testid={`param-input-${paramKey}`}
        />
      )}

      {inputType === "text" && (
        <input
          type="text"
          id={`param-${paramKey}`}
          value={String(value ?? "")}
          onChange={(e) => onUpdate(paramKey, e.target.value)}
          className="w-full rounded-lg border px-3 py-2 text-sm"
          style={{
            backgroundColor: "var(--app-surface-2)",
            borderColor: "var(--app-border)",
            color: "var(--app-text)",
          }}
          data-testid={`param-input-${paramKey}`}
        />
      )}

      {inputType === "array" && (
        <textarea
          id={`param-${paramKey}`}
          value={Array.isArray(value) ? value.join(", ") : String(value ?? "")}
          onChange={(e) => {
            const items = e.target.value
              .split(",")
              .map((s) => s.trim())
              .filter(Boolean);
            onUpdate(paramKey, items);
          }}
          rows={2}
          className="w-full rounded-lg border px-3 py-2 text-sm"
          style={{
            backgroundColor: "var(--app-surface-2)",
            borderColor: "var(--app-border)",
            color: "var(--app-text)",
            resize: "vertical",
          }}
          data-testid={`param-input-${paramKey}`}
        />
      )}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Step component                                                     */
/* ------------------------------------------------------------------ */

export function StepParameters({
  pipeline: _pipeline,
  parameters,
  originalParameters,
  onUpdate,
  onBack,
  onNext,
}: StepParametersProps) {
  // _pipeline reserved for future pipeline-specific validation
  void _pipeline;
  const [showAdvanced, setShowAdvanced] = useState(false);

  // Group parameters by category
  const grouped = useMemo(() => {
    const groups: Record<string, { key: string; meta: ParameterMeta }[]> = {};

    for (const key of Object.keys(parameters)) {
      const meta = getMeta(key);
      const groupName = meta.group;
      if (!groups[groupName]) groups[groupName] = [];
      groups[groupName].push({ key, meta });
    }

    return groups;
  }, [parameters]);

  // Separate standard and advanced parameters
  const standardGroups = useMemo(() => {
    const result: Record<string, { key: string; meta: ParameterMeta }[]> = {};
    for (const [group, items] of Object.entries(grouped)) {
      const standard = items.filter((i) => !i.meta.advanced);
      if (standard.length > 0) result[group] = standard;
    }
    return result;
  }, [grouped]);

  const advancedGroups = useMemo(() => {
    const result: Record<string, { key: string; meta: ParameterMeta }[]> = {};
    for (const [group, items] of Object.entries(grouped)) {
      const advanced = items.filter((i) => i.meta.advanced);
      if (advanced.length > 0) result[group] = advanced;
    }
    return result;
  }, [grouped]);

  const hasAdvanced = Object.keys(advancedGroups).length > 0;

  const modifiedCount = Object.keys(parameters).filter((key) =>
    isModified(key, parameters, originalParameters),
  ).length;

  return (
    <div data-testid="step-parameters">
      {/* Header */}
      <div className="mb-6">
        <h2
          className="text-lg font-bold"
          style={{ color: "var(--app-text)" }}
        >
          Adjust Parameters
        </h2>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
          Fine-tune the analysis parameters. Default values are from the
          selected preset.
          {modifiedCount > 0 && (
            <span className="ml-2 font-medium" style={{ color: "rgb(245 158 11)" }}>
              {modifiedCount} parameter{modifiedCount !== 1 ? "s" : ""} modified
            </span>
          )}
        </p>
      </div>

      {/* Standard parameter groups */}
      {Object.entries(standardGroups).map(([group, items]) => (
        <div
          key={group}
          className="mb-6 rounded-xl p-5"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
            boxShadow: "var(--app-card-shadow)",
          }}
        >
          <h3
            className="mb-4 text-xs font-semibold uppercase tracking-wider"
            style={{ color: "var(--app-text-2)" }}
          >
            {group}
          </h3>
          <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3">
            {items.map(({ key, meta }) => (
              <ParameterInput
                key={key}
                paramKey={key}
                value={parameters[key]}
                modified={isModified(key, parameters, originalParameters)}
                meta={meta}
                onUpdate={onUpdate}
              />
            ))}
          </div>
        </div>
      ))}

      {/* Advanced parameters (collapsible) */}
      {hasAdvanced && (
        <div className="mb-6">
          <button
            type="button"
            onClick={() => setShowAdvanced((v) => !v)}
            data-testid="toggle-advanced-params"
            className="mb-4 flex items-center gap-2 text-sm font-medium"
            style={{ color: "var(--app-text-2)" }}
          >
            <svg
              viewBox="0 0 24 24"
              fill="none"
              stroke="currentColor"
              strokeWidth="2"
              className="h-4 w-4 transition-transform"
              style={{
                transform: showAdvanced ? "rotate(90deg)" : "rotate(0deg)",
              }}
            >
              <path d="M9 18l6-6-6-6" />
            </svg>
            Advanced Parameters
          </button>

          {showAdvanced &&
            Object.entries(advancedGroups).map(([group, items]) => (
              <div
                key={group}
                className="mb-4 rounded-xl p-5"
                style={{
                  backgroundColor: "var(--app-surface)",
                  border: "1px solid var(--app-border)",
                  boxShadow: "var(--app-card-shadow)",
                }}
              >
                <h3
                  className="mb-4 text-xs font-semibold uppercase tracking-wider"
                  style={{ color: "var(--app-text-2)" }}
                >
                  {group}
                </h3>
                <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3">
                  {items.map(({ key, meta }) => (
                    <ParameterInput
                      key={key}
                      paramKey={key}
                      value={parameters[key]}
                      modified={isModified(
                        key,
                        parameters,
                        originalParameters,
                      )}
                      meta={meta}
                      onUpdate={onUpdate}
                    />
                  ))}
                </div>
              </div>
            ))}
        </div>
      )}

      {/* Navigation */}
      <div className="mt-8 flex items-center justify-between">
        <button
          type="button"
          onClick={onBack}
          data-testid="wizard-back-step-3"
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
          data-testid="wizard-next-step-3"
          className="rounded-lg px-6 py-2.5 text-sm font-semibold transition-opacity"
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
