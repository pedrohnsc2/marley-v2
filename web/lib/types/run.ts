/**
 * TypeScript interfaces for pipeline run data.
 * Mirrors the Python PipelineRun / StageRecord models in core/run.py.
 */

export type RunStatus =
  | "created"
  | "running"
  | "completed"
  | "failed"
  | "cancelled";

export type StageStatus =
  | "pending"
  | "running"
  | "success"
  | "failed"
  | "skipped";

export type ErrorCategory =
  | "data"
  | "network"
  | "compute"
  | "config"
  | "dependency"
  | "permission"
  | "internal";

export interface ErrorInfo {
  category: ErrorCategory;
  code: string;
  message: string;
  stage_id: string;
  suggestion: string;
}

export interface StageRecord {
  stage_id: string;
  name: string;
  status: StageStatus;
  started_at: string | null;
  completed_at: string | null;
  duration_s: number;
  error: string | null;
  error_info: ErrorInfo | null;
  output_file: string | null;
  key_metrics: Record<string, unknown>;
}

export interface RunMetadata {
  run_id: string;
  pipeline: string;
  status: RunStatus;
  created_at: string;
  started_at: string | null;
  completed_at: string | null;
  git_sha: string;
  parameters: RunParameters;
  stages: StageRecord[];
  total_duration_s: number;
  output_dir: string;
  notes: string;
  tags: string[];
}

export type RunParameters = Record<string, unknown>;
