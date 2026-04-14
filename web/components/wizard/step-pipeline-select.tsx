"use client";

import { PIPELINE_METADATA, type PipelineMetadata } from "@/lib/pipeline-metadata";

/* ------------------------------------------------------------------ */
/*  Props                                                              */
/* ------------------------------------------------------------------ */

interface StepPipelineSelectProps {
  selectedPipeline: string | null;
  onSelect: (pipelineId: string) => void;
  onNext: () => void;
}

/* ------------------------------------------------------------------ */
/*  Pipeline card                                                      */
/* ------------------------------------------------------------------ */

function PipelineCard({
  meta,
  isSelected,
  onSelect,
}: {
  meta: PipelineMetadata;
  isSelected: boolean;
  onSelect: () => void;
}) {
  return (
    <button
      type="button"
      onClick={onSelect}
      data-testid={`pipeline-card-${meta.id}`}
      className="w-full rounded-xl p-5 text-left transition-all"
      style={{
        backgroundColor: "var(--app-surface)",
        border: isSelected
          ? "2px solid var(--app-primary)"
          : "1px solid var(--app-border)",
        boxShadow: isSelected
          ? "0 0 0 3px color-mix(in srgb, var(--app-primary) 20%, transparent)"
          : "var(--app-card-shadow)",
      }}
    >
      {/* Header row */}
      <div className="mb-3 flex items-start justify-between gap-3">
        <div className="flex items-center gap-3">
          <span className="text-2xl" role="img" aria-hidden="true">
            {meta.icon}
          </span>
          <h3
            className="text-base font-bold"
            style={{ color: "var(--app-text)" }}
          >
            {meta.displayName}
          </h3>
        </div>

        {/* Duration badge */}
        <span
          className="flex-shrink-0 rounded-full px-2.5 py-0.5 text-xs font-semibold"
          style={{
            backgroundColor: "var(--app-surface-2)",
            color: "var(--app-text-2)",
            border: "1px solid var(--app-border)",
          }}
        >
          {meta.estimatedDuration}
        </span>
      </div>

      {/* Description */}
      <p
        className="mb-3 text-sm leading-relaxed"
        style={{ color: "var(--app-text-2)" }}
      >
        {meta.description}
      </p>

      {/* Expected outputs */}
      <div className="mb-2">
        <p
          className="mb-1.5 text-xs font-semibold uppercase tracking-wider"
          style={{ color: "var(--app-text-3)" }}
        >
          Expected Outputs
        </p>
        <ul className="space-y-1">
          {meta.expectedOutputs.map((output) => (
            <li
              key={output}
              className="flex items-start gap-2 text-xs"
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

      {/* Stage count */}
      <p className="text-xs" style={{ color: "var(--app-text-3)" }}>
        {meta.stages.length} stages
      </p>
    </button>
  );
}

/* ------------------------------------------------------------------ */
/*  Step component                                                     */
/* ------------------------------------------------------------------ */

export function StepPipelineSelect({
  selectedPipeline,
  onSelect,
  onNext,
}: StepPipelineSelectProps) {
  const pipelines = Object.values(PIPELINE_METADATA);

  return (
    <div data-testid="step-pipeline-select">
      {/* Header */}
      <div className="mb-6">
        <h2
          className="text-lg font-bold"
          style={{ color: "var(--app-text)" }}
        >
          Select a Pipeline
        </h2>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
          Choose the analysis pipeline you want to run. Each pipeline is
          optimized for a specific type of bioinformatics investigation.
        </p>
      </div>

      {/* Pipeline grid */}
      <div className="grid grid-cols-1 gap-4 md:grid-cols-2 xl:grid-cols-3">
        {pipelines.map((meta) => (
          <PipelineCard
            key={meta.id}
            meta={meta}
            isSelected={selectedPipeline === meta.id}
            onSelect={() => onSelect(meta.id)}
          />
        ))}
      </div>

      {/* Navigation */}
      <div className="mt-8 flex justify-end">
        <button
          type="button"
          onClick={onNext}
          disabled={!selectedPipeline}
          data-testid="wizard-next-step-1"
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
