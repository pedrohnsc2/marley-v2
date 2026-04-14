"use client";

import { useState, useEffect } from "react";
import { useTranslations } from "next-intl";
import { useRouter } from "@/i18n/routing";

interface ParamField {
  key: string;
  value: string | number | boolean;
  type: "text" | "number" | "boolean";
}

const PIPELINE_OPTIONS = [
  { id: "vaccine", label: "Reverse Vaccinology" },
  { id: "drug", label: "Drug Target Discovery" },
  { id: "docking", label: "Molecular Docking" },
  { id: "rna", label: "RNA Entropy Analysis" },
  { id: "aso_math", label: "ASO Mathematical Validation" },
];

function inferType(value: unknown): "text" | "number" | "boolean" {
  if (typeof value === "boolean") return "boolean";
  if (typeof value === "number") return "number";
  return "text";
}

function toParamFields(obj: Record<string, unknown>): ParamField[] {
  return Object.entries(obj)
    .filter(([, v]) => typeof v !== "object" || v === null)
    .map(([key, value]) => ({
      key,
      value: value as string | number | boolean,
      type: inferType(value),
    }));
}

function buildCliCommand(pipeline: string, params: ParamField[]): string {
  const parts = [`python3 -m core.run_cli --pipeline ${pipeline}`];
  for (const p of params) {
    parts.push(`--param ${p.key}=${p.value}`);
  }
  return parts.join(" \\\n  ");
}

export default function NewRunPage() {
  const t = useTranslations("runs");
  const router = useRouter();

  const [pipeline, setPipeline] = useState("vaccine");
  const [presets, setPresets] = useState<string[]>([]);
  const [selectedPreset, setSelectedPreset] = useState("");
  const [params, setParams] = useState<ParamField[]>([]);
  const [copied, setCopied] = useState(false);
  const [loadError, setLoadError] = useState("");
  const [isRunning, setIsRunning] = useState(false);
  const [runError, setRunError] = useState("");

  // Load presets when pipeline changes
  useEffect(() => {
    setLoadError("");
    fetch(`/api/presets/${pipeline}`)
      .then((r) => r.json())
      .then((data: string[]) => {
        setPresets(data);
        if (data.length > 0) {
          setSelectedPreset(data[0]);
        }
      })
      .catch(() => {
        setPresets([]);
        setLoadError("Failed to load presets");
      });
  }, [pipeline]);

  // Load preset params when preset changes
  useEffect(() => {
    if (!selectedPreset) return;
    setLoadError("");
    fetch(`/api/presets/${pipeline}?name=${selectedPreset}`)
      .then((r) => r.json())
      .then((data: Record<string, unknown>) => {
        setParams(toParamFields(data));
      })
      .catch(() => {
        setParams([]);
        setLoadError("Failed to load preset parameters");
      });
  }, [pipeline, selectedPreset]);

  const updateParam = (index: number, value: string | number | boolean) => {
    setParams((prev) => {
      const next = [...prev];
      next[index] = { ...next[index], value };
      return next;
    });
  };

  const handleRun = async () => {
    setIsRunning(true);
    setRunError("");
    try {
      const paramObj = Object.fromEntries(params.map((p) => [p.key, p.value]));
      const res = await fetch("/api/runs/start", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          pipeline,
          preset: selectedPreset || undefined,
          parameters: paramObj,
        }),
      });
      if (!res.ok) {
        const body = await res.json().catch(() => ({}));
        throw new Error(body.error || `Request failed: ${res.status}`);
      }
      const data = await res.json();
      router.push(`/runs/${data.run_id}`);
    } catch (err) {
      const message =
        err instanceof Error ? err.message : t("newRun.runError");
      setRunError(message);
    } finally {
      setIsRunning(false);
    }
  };

  const cliCommand = buildCliCommand(pipeline, params);

  const handleCopy = () => {
    navigator.clipboard.writeText(cliCommand.replace(/\\\n\s*/g, " "));
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <div className="space-y-8">
      {/* Header */}
      <div>
        <h1 className="text-2xl font-bold" style={{ color: "var(--app-text)" }}>
          {t("newRun.title")}
        </h1>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
          {t("newRun.subtitle")}
        </p>
      </div>

      {/* Pipeline + Preset selectors */}
      <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
        <div>
          <label className="mb-1.5 block text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            Pipeline
          </label>
          <select
            value={pipeline}
            onChange={(e) => setPipeline(e.target.value)}
            data-testid="pipeline-select"
            className="w-full rounded-lg border px-3 py-2 text-sm"
            style={{
              backgroundColor: "var(--app-surface)",
              borderColor: "var(--app-border)",
              color: "var(--app-text)",
            }}
          >
            {PIPELINE_OPTIONS.map((p) => (
              <option key={p.id} value={p.id}>{p.label}</option>
            ))}
          </select>
        </div>

        <div>
          <label className="mb-1.5 block text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            Preset
          </label>
          <select
            value={selectedPreset}
            onChange={(e) => setSelectedPreset(e.target.value)}
            data-testid="preset-select"
            className="w-full rounded-lg border px-3 py-2 text-sm"
            style={{
              backgroundColor: "var(--app-surface)",
              borderColor: "var(--app-border)",
              color: "var(--app-text)",
            }}
          >
            {presets.map((p) => (
              <option key={p} value={p}>{p.replace(/_/g, " ")}</option>
            ))}
          </select>
        </div>
      </div>

      {/* Error display */}
      {loadError && (
        <p className="text-sm" style={{ color: "rgb(244,63,94)" }}>{loadError}</p>
      )}

      {/* Parameters form */}
      {params.length > 0 && (
        <div
          className="rounded-xl p-5"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
            boxShadow: "var(--app-card-shadow)",
          }}
        >
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            {t("newRun.parameters")}
          </h2>
          <div className="grid grid-cols-1 gap-3 sm:grid-cols-2 lg:grid-cols-3">
            {params.map((p, i) => (
              <div key={p.key}>
                <label className="mb-1 block text-xs" style={{ color: "var(--app-text-3)" }}>
                  {p.key.replace(/_/g, " ")}
                </label>
                {p.type === "boolean" ? (
                  <button
                    onClick={() => updateParam(i, !p.value)}
                    className="rounded-lg border px-3 py-1.5 text-sm"
                    style={{
                      backgroundColor: p.value ? "rgba(16,185,129,0.1)" : "var(--app-surface-2)",
                      borderColor: p.value ? "rgba(16,185,129,0.3)" : "var(--app-border)",
                      color: "var(--app-text)",
                    }}
                  >
                    {p.value ? "true" : "false"}
                  </button>
                ) : (
                  <input
                    type={p.type}
                    value={p.value as string | number}
                    step={p.type === "number" ? "any" : undefined}
                    onChange={(e) =>
                      updateParam(
                        i,
                        p.type === "number" ? parseFloat(e.target.value) || 0 : e.target.value,
                      )
                    }
                    className="w-full rounded-lg border px-3 py-1.5 text-sm"
                    style={{
                      backgroundColor: "var(--app-surface-2)",
                      borderColor: "var(--app-border)",
                      color: "var(--app-text)",
                    }}
                  />
                )}
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Run Pipeline button */}
      <div
        className="rounded-xl p-5"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <button
          onClick={handleRun}
          disabled={isRunning}
          data-testid="run-pipeline-button"
          className="w-full rounded-lg px-4 py-3 text-sm font-semibold text-white transition-opacity disabled:opacity-60 sm:w-auto sm:min-w-[200px]"
          style={{ backgroundColor: "var(--app-accent-bar)" }}
        >
          {isRunning ? (
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
              {t("newRun.running")}
            </span>
          ) : (
            t("newRun.runButton")
          )}
        </button>
        {runError && (
          <p
            className="mt-3 text-sm"
            style={{ color: "rgb(244 63 94)" }}
            data-testid="run-error-message"
          >
            {runError}
          </p>
        )}
      </div>

      {/* CLI command */}
      <div
        className="rounded-xl p-5"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <div className="mb-3 flex items-center justify-between">
          <h2 className="text-sm font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            {t("newRun.orCli")}
          </h2>
          <button
            onClick={handleCopy}
            data-testid="copy-command-button"
            className="rounded-lg border px-3 py-1 text-xs font-medium transition-colors"
            style={{
              borderColor: "var(--app-border)",
              color: copied ? "rgb(16,185,129)" : "var(--app-text-3)",
              backgroundColor: "var(--app-surface-2)",
            }}
          >
            {copied ? t("newRun.copied") : t("newRun.copyCommand")}
          </button>
        </div>

        <pre
          className="overflow-x-auto rounded-lg p-4 text-sm leading-relaxed"
          style={{
            backgroundColor: "var(--app-surface-2)",
            color: "var(--app-text)",
            border: "1px solid var(--app-border)",
          }}
        >
          {cliCommand}
        </pre>

        <p className="mt-3 text-xs" style={{ color: "var(--app-text-3)" }}>
          {t("newRun.cliNote")}
        </p>
      </div>

      {/* JSON preview */}
      <details
        className="rounded-xl"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
        }}
      >
        <summary
          className="cursor-pointer px-5 py-3 text-sm font-medium"
          style={{ color: "var(--app-text-2)" }}
        >
          {t("newRun.jsonPreview")}
        </summary>
        <pre
          className="overflow-x-auto px-5 pb-4 text-xs"
          style={{ color: "var(--app-text-3)" }}
        >
          {JSON.stringify(
            Object.fromEntries(params.map((p) => [p.key, p.value])),
            null,
            2,
          )}
        </pre>
      </details>
    </div>
  );
}
