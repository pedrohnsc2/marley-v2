import path from "path";

/**
 * Validates that a resolved path stays within the base directory.
 * Throws if path traversal is detected.
 */
export function safePath(base: string, ...segments: string[]): string {
  const resolved = path.resolve(base, ...segments);
  const realBase = path.resolve(base);
  if (!resolved.startsWith(realBase + path.sep) && resolved !== realBase) {
    throw new Error("Path traversal detected");
  }
  return resolved;
}

/**
 * Validates an ID (runId, pipeline name, preset name) contains only safe characters.
 * Allows: alphanumeric, underscore, hyphen, dot (no ..)
 */
export function validateId(id: string): void {
  if (!id || /[^a-zA-Z0-9_\-.]/.test(id) || id.includes("..")) {
    throw new Error("Invalid identifier");
  }
}
