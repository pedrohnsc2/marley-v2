import fs from "fs";
import path from "path";

const RESULTS_DIR = path.join(process.cwd(), "..", "results");

const MODULE_DIRS: Record<string, string> = {
  results: RESULTS_DIR,
  aso_math: path.join(process.cwd(), "..", "aso_math", "results"),
  aso_delivery: path.join(process.cwd(), "..", "aso_delivery", "results"),
  marley_ai: path.join(process.cwd(), "..", "marley_ai", "results"),
  qaoa: path.join(process.cwd(), "..", "mrl_quantum", "results", "qaoa"),
  vqe: path.join(process.cwd(), "..", "mrl_quantum", "results", "vqe"),
  platforms: path.join(process.cwd(), "..", "vaccine_platforms", "reports", "results"),
};

export function loadModuleJson(module: keyof typeof MODULE_DIRS, filename: string): unknown {
  const dir = MODULE_DIRS[module];
  if (!dir) throw new Error(`Unknown module: ${module}`);
  const filePath = path.join(dir, filename);
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
  const dir = MODULE_DIRS[module];
  if (!dir) throw new Error(`Unknown module: ${module}`);
  const filePath = path.join(dir, filename);
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
  const filePath = path.join(RESULTS_DIR, filename);
  const content = fs.readFileSync(filePath, "utf-8");
  return JSON.parse(content);
}

export function loadCsv(filename: string): Record<string, string>[] {
  const filePath = path.join(RESULTS_DIR, filename);
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
  const filePath = path.join(RESULTS_DIR, filename);
  return fs.readFileSync(filePath, "utf-8");
}

const PROJECT_ROOT = path.join(process.cwd(), "..");

export function loadPdb(relativePath: string): string {
  const filePath = path.join(PROJECT_ROOT, relativePath);
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
