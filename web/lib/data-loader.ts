import fs from "fs";
import path from "path";
import { safePath, validateId } from "@/lib/safe-path";

const PROJECT_ROOT = path.join(process.cwd(), "..");
const RESULTS_DIR = path.join(PROJECT_ROOT, "results");
const LATEST_DIR = path.join(RESULTS_DIR, "latest");
const RUNS_DIR = path.join(RESULTS_DIR, "runs");

const MODULE_DIRS: Record<string, string> = {
  results: RESULTS_DIR,
  aso_math: path.join(PROJECT_ROOT, "aso_math", "results"),
  aso_math_reports: path.join(PROJECT_ROOT, "aso_math", "reports", "results"),
  aso_delivery: path.join(PROJECT_ROOT, "aso_delivery", "results"),
  marley_ai: path.join(PROJECT_ROOT, "marley_ai", "results"),
  qaoa: path.join(PROJECT_ROOT, "mrl_quantum", "results", "qaoa"),
  vqe: path.join(PROJECT_ROOT, "mrl_quantum", "results", "vqe"),
  platforms: path.join(PROJECT_ROOT, "vaccine_platforms", "reports", "results"),
};

// ---------------------------------------------------------------------------
// Run-aware resolution: if results/latest/<pipeline> exists as a symlink,
// use it. Otherwise fall back to the hardcoded MODULE_DIRS paths.
// This enables dynamic pipeline runs without breaking existing static data.
// ---------------------------------------------------------------------------

function resolveModuleDir(module: string): string {
  validateId(module);
  // Check if latest/ symlink exists for this module
  const latestLink = path.join(LATEST_DIR, module);
  try {
    const stat = fs.lstatSync(latestLink);
    if (stat.isSymbolicLink() || stat.isDirectory()) {
      return fs.realpathSync(latestLink);
    }
  } catch {
    // No latest symlink — fall through to static dirs
  }

  const dir = MODULE_DIRS[module];
  if (!dir) throw new Error(`Unknown module: ${module}`);
  return dir;
}

// ---------------------------------------------------------------------------
// Run-specific data loading
// ---------------------------------------------------------------------------

export function getRunDir(runId: string): string {
  validateId(runId);
  return path.join(RUNS_DIR, runId, "outputs");
}

export function loadRunMetadata(runId: string): unknown {
  validateId(runId);
  const metaPath = path.join(RUNS_DIR, runId, "metadata.json");
  return JSON.parse(fs.readFileSync(metaPath, "utf-8"));
}

export function listRuns(pipeline?: string): unknown[] {
  const runsPath = RUNS_DIR;
  if (!fs.existsSync(runsPath)) return [];

  const entries = fs.readdirSync(runsPath, { withFileTypes: true });
  const runs: unknown[] = [];

  for (const entry of entries) {
    if (!entry.isDirectory()) continue;
    const metaPath = path.join(runsPath, entry.name, "metadata.json");
    try {
      const meta = JSON.parse(fs.readFileSync(metaPath, "utf-8"));
      if (!pipeline || meta.pipeline === pipeline) {
        runs.push(meta);
      }
    } catch {
      continue;
    }
  }

  return runs.sort((a, b) => {
    const dateA = (a as Record<string, string>).created_at ?? "";
    const dateB = (b as Record<string, string>).created_at ?? "";
    return dateB.localeCompare(dateA);
  });
}

// ---------------------------------------------------------------------------
// Module-based loading (existing API, now run-aware)
// ---------------------------------------------------------------------------

export function loadModuleJson(module: keyof typeof MODULE_DIRS, filename: string): unknown {
  const dir = resolveModuleDir(module);
  const filePath = safePath(dir, filename);
  return JSON.parse(fs.readFileSync(filePath, "utf-8"));
}

export function safeLoadModuleJson(module: keyof typeof MODULE_DIRS, filename: string): unknown | null {
  try {
    return loadModuleJson(module, filename);
  } catch {
    return null;
  }
}

export function loadModuleCsv(module: keyof typeof MODULE_DIRS, filename: string): Record<string, string>[] {
  const dir = resolveModuleDir(module);
  const filePath = safePath(dir, filename);
  return parseCsv(fs.readFileSync(filePath, "utf-8"));
}

function parseCsv(content: string): Record<string, string>[] {
  const lines = content.trim().split("\n");
  if (lines.length === 0) return [];
  const headers = lines[0].split(",");
  return lines.slice(1).map((line) => {
    const values = line.split(",");
    const row: Record<string, string> = {};
    headers.forEach((header, i) => {
      row[header.trim()] = (values[i] ?? "").trim();
    });
    return row;
  });
}

export function loadJson(filename: string): unknown {
  const filePath = safePath(RESULTS_DIR, filename);
  const content = fs.readFileSync(filePath, "utf-8");
  return JSON.parse(content);
}

export function loadCsv(filename: string): Record<string, string>[] {
  const filePath = safePath(RESULTS_DIR, filename);
  const content = fs.readFileSync(filePath, "utf-8");
  const lines = content.trim().split("\n");
  if (lines.length === 0) return [];
  const headers = lines[0].split(",");
  return lines.slice(1).map((line) => {
    const values = line.split(",");
    const row: Record<string, string> = {};
    headers.forEach((header, i) => {
      row[header.trim()] = (values[i] ?? "").trim();
    });
    return row;
  });
}

export function loadMarkdown(filename: string): string {
  const filePath = safePath(RESULTS_DIR, filename);
  return fs.readFileSync(filePath, "utf-8");
}

export function loadPdb(relativePath: string): string {
  const filePath = safePath(PROJECT_ROOT, relativePath);
  return fs.readFileSync(filePath, "utf-8");
}

export function safeLoadPdb(relativePath: string): string | null {
  try {
    return loadPdb(relativePath);
  } catch {
    return null;
  }
}

export function safeLoadJson(filename: string): unknown | null {
  try {
    return loadJson(filename);
  } catch {
    return null;
  }
}

export function safeLoadCsv(filename: string): Record<string, string>[] {
  try {
    return loadCsv(filename);
  } catch {
    return [];
  }
}

// ---------------------------------------------------------------------------
// Run-specific data loading by filename (mirrors loadJson/loadCsv but reads
// from the run's output directory instead of the global results dir)
// ---------------------------------------------------------------------------

export function loadRunJson(runId: string, filename: string): unknown {
  const runDir = getRunDir(runId);
  const filePath = safePath(runDir, filename);
  return JSON.parse(fs.readFileSync(filePath, "utf-8"));
}

export function safeLoadRunJson(runId: string, filename: string): unknown | null {
  try {
    return loadRunJson(runId, filename);
  } catch {
    return null;
  }
}

export function loadRunCsv(runId: string, filename: string): Record<string, string>[] {
  const runDir = getRunDir(runId);
  const filePath = safePath(runDir, filename);
  return parseCsv(fs.readFileSync(filePath, "utf-8"));
}

export function safeLoadRunCsv(runId: string, filename: string): Record<string, string>[] {
  try {
    return loadRunCsv(runId, filename);
  } catch {
    return [];
  }
}

export function loadRunModuleJson(runId: string, module: string, filename: string): unknown {
  validateId(module);
  // For run-specific data, try to find the module's output within the run dir
  const runDir = getRunDir(runId);
  const filePath = safePath(runDir, module, filename);
  return JSON.parse(fs.readFileSync(filePath, "utf-8"));
}

export function safeLoadRunModuleJson(runId: string, module: string, filename: string): unknown | null {
  try {
    return loadRunModuleJson(runId, module, filename);
  } catch {
    return null;
  }
}

export function safeLoadRunPdb(runId: string, relativePath: string): string | null {
  try {
    const runDir = getRunDir(runId);
    const filePath = safePath(runDir, relativePath);
    return fs.readFileSync(filePath, "utf-8");
  } catch {
    return null;
  }
}

export function getRunCompletedDate(runId: string): string | null {
  try {
    const meta = loadRunMetadata(runId) as Record<string, string>;
    return meta.completed_at ?? meta.created_at ?? null;
  } catch {
    return null;
  }
}
