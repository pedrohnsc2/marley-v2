"use client";

import { useState, useCallback } from "react";
import { StepPipelineSelect } from "./step-pipeline-select";
import { StepPresetSelect } from "./step-preset-select";
import { StepParameters } from "./step-parameters";
import { StepConfirm } from "./step-confirm";

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

export interface WizardState {
  pipeline: string | null;
  preset: string | null;
  parameters: Record<string, unknown>;
  /** Snapshot of the original preset parameters, used to detect modifications */
  originalParameters: Record<string, unknown>;
}

const INITIAL_STATE: WizardState = {
  pipeline: null,
  preset: null,
  parameters: {},
  originalParameters: {},
};

const STEPS = [
  { number: 1, label: "Pipeline" },
  { number: 2, label: "Configuration" },
  { number: 3, label: "Parameters" },
  { number: 4, label: "Launch" },
] as const;

/* ------------------------------------------------------------------ */
/*  Stepper indicator                                                  */
/* ------------------------------------------------------------------ */

function WizardStepper({ currentStep }: { currentStep: number }) {
  return (
    <nav
      aria-label="Wizard progress"
      className="mb-8"
      data-testid="wizard-stepper"
    >
      <ol className="flex items-center gap-2">
        {STEPS.map((step, idx) => {
          const isActive = step.number === currentStep;
          const isCompleted = step.number < currentStep;
          const isFuture = step.number > currentStep;

          return (
            <li
              key={step.number}
              className="flex items-center gap-2"
              style={{ flex: idx < STEPS.length - 1 ? 1 : "none" }}
            >
              {/* Circle */}
              <div
                className="flex h-8 w-8 flex-shrink-0 items-center justify-center rounded-full text-xs font-bold transition-colors"
                style={{
                  backgroundColor: isCompleted
                    ? "var(--app-primary)"
                    : isActive
                      ? "var(--app-primary)"
                      : "var(--app-surface-2)",
                  color: isCompleted || isActive
                    ? "var(--app-primary-tx)"
                    : "var(--app-text-3)",
                  border: isFuture
                    ? "2px solid var(--app-border)"
                    : "none",
                }}
                data-testid={`wizard-step-${step.number}`}
              >
                {isCompleted ? (
                  <svg
                    viewBox="0 0 12 12"
                    className="h-3.5 w-3.5"
                    fill="none"
                    stroke="currentColor"
                    strokeWidth="2"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  >
                    <path d="M2.5 6.5L5 9L9.5 3.5" />
                  </svg>
                ) : (
                  step.number
                )}
              </div>

              {/* Label */}
              <span
                className="hidden text-sm font-medium sm:inline"
                style={{
                  color: isActive
                    ? "var(--app-text)"
                    : isCompleted
                      ? "var(--app-text-2)"
                      : "var(--app-text-3)",
                }}
              >
                {step.label}
              </span>

              {/* Connector line */}
              {idx < STEPS.length - 1 && (
                <div
                  className="mx-2 hidden h-0.5 flex-1 rounded-full sm:block"
                  style={{
                    backgroundColor: isCompleted
                      ? "var(--app-primary)"
                      : "var(--app-border)",
                  }}
                />
              )}
            </li>
          );
        })}
      </ol>
    </nav>
  );
}

/* ------------------------------------------------------------------ */
/*  Main wizard                                                        */
/* ------------------------------------------------------------------ */

export function PipelineWizard() {
  const [step, setStep] = useState(1);
  const [state, setState] = useState<WizardState>(INITIAL_STATE);

  const goNext = useCallback(() => setStep((s) => Math.min(s + 1, 4)), []);
  const goBack = useCallback(() => setStep((s) => Math.max(s - 1, 1)), []);

  const selectPipeline = useCallback(
    (pipelineId: string) => {
      setState((prev) => ({
        ...prev,
        pipeline: pipelineId,
        preset: null,
        parameters: {},
        originalParameters: {},
      }));
    },
    [],
  );

  const selectPreset = useCallback(
    (presetName: string, presetParams: Record<string, unknown>) => {
      setState((prev) => ({
        ...prev,
        preset: presetName,
        parameters: { ...presetParams },
        originalParameters: { ...presetParams },
      }));
    },
    [],
  );

  const updateParameter = useCallback(
    (key: string, value: unknown) => {
      setState((prev) => ({
        ...prev,
        parameters: { ...prev.parameters, [key]: value },
      }));
    },
    [],
  );

  return (
    <div data-testid="pipeline-wizard">
      <WizardStepper currentStep={step} />

      {step === 1 && (
        <StepPipelineSelect
          selectedPipeline={state.pipeline}
          onSelect={selectPipeline}
          onNext={goNext}
        />
      )}

      {step === 2 && (
        <StepPresetSelect
          pipeline={state.pipeline!}
          selectedPreset={state.preset}
          onSelect={selectPreset}
          onBack={goBack}
          onNext={goNext}
        />
      )}

      {step === 3 && (
        <StepParameters
          pipeline={state.pipeline!}
          parameters={state.parameters}
          originalParameters={state.originalParameters}
          onUpdate={updateParameter}
          onBack={goBack}
          onNext={goNext}
        />
      )}

      {step === 4 && (
        <StepConfirm
          state={state}
          onBack={goBack}
        />
      )}
    </div>
  );
}
