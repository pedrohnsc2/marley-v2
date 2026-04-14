import { NextRequest, NextResponse } from "next/server";
import { execFileSync } from "child_process";
import path from "path";
import { requireWriteAuth } from "@/lib/api-auth";
import { rateLimit, acquireRunSlot, releaseRunSlot } from "@/lib/rate-limit";

const PROJECT_ROOT = path.join(process.cwd(), "..");
const PYTHON_BIN = path.join(PROJECT_ROOT, ".venv", "bin", "python");
const ALLOWED_PIPELINES = new Set([
  "vaccine",
  "drug",
  "docking",
  "rna",
  "aso_math",
]);
const PRESET_RE = /^[a-z0-9_]{1,50}$/;

export async function POST(request: NextRequest) {
  // Auth check
  const authErr = requireWriteAuth(request);
  if (authErr) return authErr;

  // Rate limit: 3 requests per minute
  const limited = rateLimit(request, { maxRequests: 3, windowMs: 60_000 });
  if (limited) return limited;

  // Concurrency check
  if (!acquireRunSlot()) {
    return NextResponse.json(
      { error: "Too many concurrent runs. Maximum 3 allowed." },
      { status: 429 },
    );
  }

  try {
    const body = await request.json();

    // Validate pipeline
    const pipeline = body.pipeline;
    if (!pipeline || !ALLOWED_PIPELINES.has(pipeline)) {
      releaseRunSlot();
      return NextResponse.json(
        {
          error:
            "Invalid pipeline. Allowed: vaccine, drug, docking, rna, aso_math",
        },
        { status: 400 },
      );
    }

    // Validate preset (optional)
    if (body.preset && !PRESET_RE.test(body.preset)) {
      releaseRunSlot();
      return NextResponse.json(
        { error: "Invalid preset name" },
        { status: 400 },
      );
    }

    // Validate parameters -- reject dangerous keys/values
    if (body.parameters && typeof body.parameters === "object") {
      for (const [key, value] of Object.entries(body.parameters)) {
        if (typeof key !== "string" || /[^a-zA-Z0-9_]/.test(key)) {
          releaseRunSlot();
          return NextResponse.json(
            { error: `Invalid parameter key: ${key}` },
            { status: 400 },
          );
        }
        if (
          typeof value === "string" &&
          (value.includes("..") || value.includes("/") || value.includes("\\"))
        ) {
          releaseRunSlot();
          return NextResponse.json(
            { error: "Invalid parameter value" },
            { status: 400 },
          );
        }
      }
    }

    // NEVER accept config_path from HTTP input
    const sanitizedBody = {
      pipeline: body.pipeline,
      preset: body.preset || undefined,
      parameters: body.parameters || {},
      tags: Array.isArray(body.tags)
        ? body.tags.map(String).slice(0, 10)
        : [],
      notes: typeof body.notes === "string" ? body.notes.slice(0, 500) : "",
    };

    // Call Python launcher via execFileSync (no shell involved)
    const result = execFileSync(PYTHON_BIN, ["-m", "core.launcher"], {
      input: JSON.stringify(sanitizedBody),
      cwd: PROJECT_ROOT,
      timeout: 10_000,
      encoding: "utf-8",
    });

    const parsed = JSON.parse(result.trim());

    // Release run slot after spawn since the Python process runs independently
    releaseRunSlot();

    return NextResponse.json(parsed, { status: 201 });
  } catch {
    releaseRunSlot();
    return NextResponse.json(
      { error: "Failed to start pipeline" },
      { status: 500 },
    );
  }
}
