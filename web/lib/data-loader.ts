import fs from "fs";
import path from "path";

const RESULTS_DIR = path.join(process.cwd(), "..", "results");

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
