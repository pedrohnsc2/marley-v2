"use client";

import { useState, useCallback } from "react";
import type { ErrorInfo } from "@/lib/types/run";

interface FriendlyErrorProps {
  errorInfo: ErrorInfo | null;
  rawError: string | null;
  runId: string;
  stageId: string;
}

const ERROR_DISPLAY: Record<string, { title: string; body: string; suggestion: string }> = {
  "data.missing_input": {
    title: "Input Data Not Found",
    body: "A required input file was not found. This usually means a previous step did not complete.",
    suggestion: "Re-run the pipeline from the beginning.",
  },
  "data.invalid_format": {
    title: "Invalid Data Format",
    body: "An input file could not be read because it has an unexpected format.",
    suggestion: "Verify that input files have not been modified manually.",
  },
  "data.empty_result": {
    title: "No Results Found",
    body: "This step completed but found no matching results. The filtering criteria may be too restrictive.",
    suggestion: "Try adjusting the pipeline parameters.",
  },
  "network.timeout": {
    title: "Server Timeout",
    body: "A request to an external server timed out.",
    suggestion: "Check your internet connection and try again.",
  },
  "network.connection_refused": {
    title: "Connection Failed",
    body: "Could not connect to an external service required by this step.",
    suggestion: "Verify your internet connection.",
  },
  "network.dns_failure": {
    title: "DNS Resolution Failed",
    body: "The server address could not be resolved.",
    suggestion: "Check your network configuration.",
  },
  "network.api_error": {
    title: "External Service Error",
    body: "An external bioinformatics service returned an error.",
    suggestion: "The service may be temporarily unavailable. Try again later.",
  },
  "network.ssl_error": {
    title: "Secure Connection Failed",
    body: "Could not establish a secure connection to an external service.",
    suggestion: "This may be caused by a corporate proxy or firewall.",
  },
  "compute.out_of_memory": {
    title: "Insufficient Memory",
    body: "The system ran out of memory while processing this step.",
    suggestion: "Try reducing the dataset size.",
  },
  "compute.numerical_error": {
    title: "Computation Error",
    body: "A numerical calculation produced an invalid result.",
    suggestion: "Check your input data for extreme or missing values.",
  },
  "compute.timeout": {
    title: "Computation Timeout",
    body: "This step exceeded its time limit.",
    suggestion: "Try reducing the dataset size.",
  },
  "compute.subprocess_failed": {
    title: "External Tool Error",
    body: "An external bioinformatics tool exited with an error.",
    suggestion: "Check that the tool is properly installed.",
  },
  "config.invalid_parameter": {
    title: "Invalid Parameter",
    body: "One of the pipeline parameters has an incorrect value.",
    suggestion: "Review the parameters and correct the invalid value.",
  },
  "config.missing_parameter": {
    title: "Missing Parameter",
    body: "A required pipeline parameter was not provided.",
    suggestion: "Check the preset configuration.",
  },
  "config.preset_not_found": {
    title: "Preset Not Found",
    body: "The selected preset does not exist.",
    suggestion: "Use the pipeline wizard to see available presets.",
  },
  "dependency.module_missing": {
    title: "Missing Software Package",
    body: "A required Python package is not installed.",
    suggestion: "Run: pip install -r requirements.txt",
  },
  "dependency.tool_missing": {
    title: "Missing External Tool",
    body: "A required bioinformatics tool is not installed.",
    suggestion: "See the Methods page for required software.",
  },
  "permission.file_access": {
    title: "File Permission Denied",
    body: "The pipeline could not read or write a required file.",
    suggestion: "Check directory permissions.",
  },
  "permission.api_auth": {
    title: "Authentication Failed",
    body: "An API key was rejected by an external service.",
    suggestion: "Verify API keys in the .env file.",
  },
  "internal.assertion": {
    title: "Internal Error",
    body: "An internal consistency check failed. This is likely a software bug.",
    suggestion: "Please report this issue with the run ID.",
  },
  "internal.unknown": {
    title: "Unexpected Error",
    body: "An unexpected error occurred during this step.",
    suggestion: "Try again. If it persists, report the issue.",
  },
};

const FALLBACK = {
  title: "Error",
  body: "This step encountered an error.",
  suggestion: "Check the technical details for more information.",
};

function getDisplay(errorInfo: ErrorInfo | null): { title: string; body: string; suggestion: string } {
  if (!errorInfo) return FALLBACK;
  return ERROR_DISPLAY[errorInfo.code] ?? FALLBACK;
}

export default function FriendlyError({
  errorInfo,
  rawError,
  runId,
  stageId,
}: FriendlyErrorProps) {
  const [techOpen, setTechOpen] = useState(false);
  const [copied, setCopied] = useState(false);

  const display = getDisplay(errorInfo);
  const suggestion = errorInfo?.suggestion || display.suggestion;

  const handleCopyDiagnostics = useCallback(async () => {
    const diagnostics = [
      `Run ID: ${runId}`,
      `Stage ID: ${stageId}`,
      errorInfo ? `Error code: ${errorInfo.code}` : null,
      errorInfo ? `Category: ${errorInfo.category}` : null,
      errorInfo ? `Message: ${errorInfo.message}` : null,
      rawError ? `\nRaw error:\n${rawError}` : null,
    ]
      .filter(Boolean)
      .join("\n");

    try {
      await navigator.clipboard.writeText(diagnostics);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch {
      // Fallback: select a hidden textarea (not needed in most modern browsers)
    }
  }, [runId, stageId, errorInfo, rawError]);

  return (
    <div
      className="mt-2 rounded-lg p-4"
      style={{
        backgroundColor: "rgb(244 63 94 / 0.06)",
        border: "1px solid rgb(244 63 94 / 0.2)",
      }}
      data-testid={`friendly-error-${stageId}`}
    >
      {/* Title */}
      <p
        className="text-sm font-semibold"
        style={{ color: "rgb(244 63 94)" }}
        data-testid={`friendly-error-title-${stageId}`}
      >
        {display.title}
      </p>

      {/* Description */}
      <p
        className="mt-1 text-sm"
        style={{ color: "var(--app-text-2)" }}
        data-testid={`friendly-error-body-${stageId}`}
      >
        {display.body}
      </p>

      {/* Suggestion */}
      <div
        className="mt-3 flex items-start gap-2 rounded-md p-2.5"
        style={{
          backgroundColor: "var(--app-surface-2)",
          border: "1px solid var(--app-border)",
        }}
      >
        <svg
          viewBox="0 0 16 16"
          className="mt-0.5 h-3.5 w-3.5 flex-shrink-0"
          fill="currentColor"
          style={{ color: "var(--app-text-3)" }}
        >
          <path d="M8 1.5a6.5 6.5 0 1 0 0 13 6.5 6.5 0 0 0 0-13zM0 8a8 8 0 1 1 16 0A8 8 0 0 1 0 8zm6.5-.25A.75.75 0 0 1 7.25 7h1a.75.75 0 0 1 .75.75v2.75h.25a.75.75 0 0 1 0 1.5h-2a.75.75 0 0 1 0-1.5h.25v-2h-.25a.75.75 0 0 1-.75-.75zM8 6a1 1 0 1 1 0-2 1 1 0 0 1 0 2z" />
        </svg>
        <p
          className="text-xs"
          style={{ color: "var(--app-text-2)" }}
          data-testid={`friendly-error-suggestion-${stageId}`}
        >
          {suggestion}
        </p>
      </div>

      {/* Actions */}
      <div className="mt-3 flex flex-wrap items-center gap-2">
        <button
          onClick={handleCopyDiagnostics}
          className="rounded-md px-3 py-1.5 text-xs font-medium transition-opacity hover:opacity-80"
          style={{
            backgroundColor: "var(--app-surface-2)",
            color: "var(--app-text-2)",
            border: "1px solid var(--app-border)",
          }}
          data-testid={`friendly-error-report-${stageId}`}
        >
          {copied ? "Copied!" : "Report Issue"}
        </button>

        {rawError && (
          <button
            onClick={() => setTechOpen((prev) => !prev)}
            className="text-xs font-medium underline"
            style={{ color: "var(--app-text-3)" }}
            data-testid={`friendly-error-tech-toggle-${stageId}`}
          >
            {techOpen ? "Hide Technical Details" : "Technical Details"}
          </button>
        )}
      </div>

      {/* Technical details (collapsible) */}
      {techOpen && rawError && (
        <pre
          className="mt-2 overflow-x-auto rounded-lg p-3 text-xs"
          style={{
            backgroundColor: "var(--app-surface-2)",
            color: "var(--app-text-3)",
            border: "1px solid var(--app-border)",
          }}
          data-testid={`friendly-error-raw-${stageId}`}
        >
          {rawError}
        </pre>
      )}
    </div>
  );
}
