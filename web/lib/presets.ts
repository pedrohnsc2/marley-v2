import fs from "fs";
import path from "path";
import { validateId } from "@/lib/safe-path";

const PRESETS_DIR = path.join(process.cwd(), "..", "configs", "presets");

/* ------------------------------------------------------------------ */
/*  Types                                                              */
/* ------------------------------------------------------------------ */

export interface PresetMeta {
  display_name: string;
  description: string;
  recommended?: boolean;
}

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

/**
 * Separate the `_meta` block from the actual pipeline parameters.
 * Returns null for meta if the preset has no `_meta` key (backwards compat).
 */
export function splitPreset(raw: Record<string, unknown>): {
  meta: PresetMeta | null;
  params: Record<string, unknown>;
} {
  const { _meta, ...params } = raw;
  return {
    meta: (_meta as PresetMeta) ?? null,
    params,
  };
}

/**
 * List available preset names for a given pipeline.
 * Returns an array of preset names (filename without .json extension).
 */
export function listPresets(pipeline: string): string[] {
  validateId(pipeline);
  const dir = path.join(PRESETS_DIR, pipeline);
  try {
    const entries = fs.readdirSync(dir, { withFileTypes: true });
    return entries
      .filter((e) => e.isFile() && e.name.endsWith(".json"))
      .map((e) => e.name.replace(/\.json$/, ""))
      .sort();
  } catch {
    return [];
  }
}

/**
 * Load a specific preset's configuration as a JSON object.
 * Returns null if the preset does not exist.
 */
export function loadPreset(
  pipeline: string,
  name: string,
): Record<string, unknown> | null {
  validateId(pipeline);
  validateId(name);
  const filePath = path.join(PRESETS_DIR, pipeline, `${name}.json`);
  try {
    const content = fs.readFileSync(filePath, "utf-8");
    return JSON.parse(content) as Record<string, unknown>;
  } catch {
    return null;
  }
}

/**
 * List all pipelines that have at least one preset defined.
 */
export function listPipelinesWithPresets(): string[] {
  try {
    const entries = fs.readdirSync(PRESETS_DIR, { withFileTypes: true });
    return entries
      .filter((e) => e.isDirectory())
      .map((e) => e.name)
      .sort();
  } catch {
    return [];
  }
}
