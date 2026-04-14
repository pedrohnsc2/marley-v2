"use client";

import { PipelineWizard } from "@/components/wizard/pipeline-wizard";

export default function NewRunPage() {
  return (
    <div className="mx-auto max-w-5xl">
      {/* Page header */}
      <div className="mb-8">
        <h1
          className="text-2xl font-bold"
          style={{ color: "var(--app-text)" }}
        >
          New Pipeline Run
        </h1>
        <p
          className="mt-1 text-sm"
          style={{ color: "var(--app-text-3)" }}
        >
          Configure and launch a new bioinformatics analysis pipeline.
        </p>
      </div>

      {/* Wizard */}
      <PipelineWizard />
    </div>
  );
}
