/**
 * Storage abstraction for loading pipeline run results.
 * Uses local filesystem in dev, Supabase Storage in production.
 */

import { createClient } from "@/lib/supabase/server";

const STORAGE_BUCKET = process.env.STORAGE_BUCKET || "pipeline-results";

export interface StorageLoader {
  loadFile(runId: string, filename: string): Promise<Buffer | null>;
  loadJson<T = unknown>(runId: string, filename: string): Promise<T | null>;
  loadText(runId: string, filename: string): Promise<string | null>;
  getPublicUrl(runId: string, filename: string): string | null;
}

// --- Local filesystem loader (dev) ---
import fs from "fs";
import path from "path";
import { safePath } from "@/lib/safe-path";

const RUNS_DIR = path.join(process.cwd(), "..", "results", "runs");

class LocalStorageLoader implements StorageLoader {
  async loadFile(runId: string, filename: string): Promise<Buffer | null> {
    try {
      const filePath = safePath(path.join(RUNS_DIR, runId), filename);
      return fs.readFileSync(filePath);
    } catch {
      return null;
    }
  }

  async loadJson<T = unknown>(runId: string, filename: string): Promise<T | null> {
    const buf = await this.loadFile(runId, filename);
    if (!buf) return null;
    return JSON.parse(buf.toString("utf-8")) as T;
  }

  async loadText(runId: string, filename: string): Promise<string | null> {
    const buf = await this.loadFile(runId, filename);
    return buf ? buf.toString("utf-8") : null;
  }

  getPublicUrl(): string | null {
    return null; // No public URLs in local mode
  }
}

// --- Supabase Storage loader (production) ---

class SupabaseStorageLoader implements StorageLoader {
  async loadFile(runId: string, filename: string): Promise<Buffer | null> {
    try {
      const supabase = await createClient();
      const storagePath = `${runId}/${filename}`;
      const { data, error } = await supabase.storage
        .from(STORAGE_BUCKET)
        .download(storagePath);
      if (error || !data) return null;
      const arrayBuffer = await data.arrayBuffer();
      return Buffer.from(arrayBuffer);
    } catch {
      return null;
    }
  }

  async loadJson<T = unknown>(runId: string, filename: string): Promise<T | null> {
    const buf = await this.loadFile(runId, filename);
    if (!buf) return null;
    return JSON.parse(buf.toString("utf-8")) as T;
  }

  async loadText(runId: string, filename: string): Promise<string | null> {
    const buf = await this.loadFile(runId, filename);
    return buf ? buf.toString("utf-8") : null;
  }

  getPublicUrl(runId: string, filename: string): string | null {
    const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL;
    if (!supabaseUrl) return null;
    return `${supabaseUrl}/storage/v1/object/public/${STORAGE_BUCKET}/${runId}/${filename}`;
  }
}

// --- Factory ---

let _instance: StorageLoader | null = null;

export function getStorageLoader(): StorageLoader {
  if (_instance) return _instance;

  const useSupabase =
    process.env.NEXT_PUBLIC_SUPABASE_URL &&
    process.env.NODE_ENV === "production";

  _instance = useSupabase ? new SupabaseStorageLoader() : new LocalStorageLoader();
  return _instance;
}
