import { NextRequest, NextResponse } from "next/server";
import { requireWriteAuth } from "@/lib/api-auth";
import { rateLimit, acquireRunSlot, releaseRunSlot } from "@/lib/rate-limit";

const WORKER_API_URL = process.env.NEXT_PUBLIC_WORKER_API_URL || "";
const API_SECRET_KEY = process.env.API_SECRET_KEY || "";

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
  const { response: authErr, userId } = await requireWriteAuth(request);
  if (authErr) return authErr;

  // Rate limit: 3 requests per minute
  const limited = rateLimit(request, { maxRequests: 3, windowMs: 60_000 });
  if (limited) return limited;

  try {
    const body = await request.json();

    // Validate pipeline
    const pipeline = body.pipeline;
    if (!pipeline || !ALLOWED_PIPELINES.has(pipeline)) {
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
      return NextResponse.json(
        { error: "Invalid preset name" },
        { status: 400 },
      );
    }

    // Validate parameters -- reject dangerous keys/values
    if (body.parameters && typeof body.parameters === "object") {
      for (const [key, value] of Object.entries(body.parameters)) {
        if (typeof key !== "string" || /[^a-zA-Z0-9_]/.test(key)) {
          return NextResponse.json(
            { error: `Invalid parameter key: ${key}` },
            { status: 400 },
          );
        }
        if (
          typeof value === "string" &&
          (value.includes("..") || value.includes("/") || value.includes("\\"))
        ) {
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
      user_id: userId ?? undefined,
    };

    // -----------------------------------------------------------------------
    // Production: forward to Worker API (FastAPI + rq queue)
    // Development: fall back to local subprocess (execFileSync)
    // -----------------------------------------------------------------------
    if (WORKER_API_URL) {
      const workerResp = await fetch(`${WORKER_API_URL}/api/runs/start`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          ...(API_SECRET_KEY ? { "x-api-key": API_SECRET_KEY } : {}),
        },
        body: JSON.stringify(sanitizedBody),
      });

      if (!workerResp.ok) {
        const err = await workerResp.json().catch(() => ({}));
        return NextResponse.json(
          { error: err.detail || "Worker API error" },
          { status: workerResp.status },
        );
      }

      const parsed = await workerResp.json();
      return NextResponse.json(parsed, { status: 201 });
    }

    // --- Local development fallback: direct subprocess ---
    if (!acquireRunSlot()) {
      return NextResponse.json(
        { error: "Too many concurrent runs. Maximum 3 allowed." },
        { status: 429 },
      );
    }

    try {
      const { execFileSync } = await import("child_process");
      const path = await import("path");

      const PROJECT_ROOT = path.join(process.cwd(), "..");
      const PYTHON_BIN = path.join(PROJECT_ROOT, ".venv", "bin", "python");

      const result = execFileSync(PYTHON_BIN, ["-m", "core.launcher"], {
        input: JSON.stringify(sanitizedBody),
        cwd: PROJECT_ROOT,
        timeout: 10_000,
        encoding: "utf-8",
        stdio: ["pipe", "pipe", "pipe"],
      });

      const lines = result.trim().split("\n");
      const jsonLine = lines.reverse().find((l: string) => l.startsWith("{"));
      if (!jsonLine) {
        throw new Error("No JSON output from launcher");
      }
      const parsed = JSON.parse(jsonLine);

      releaseRunSlot();
      return NextResponse.json(parsed, { status: 201 });
    } catch {
      releaseRunSlot();
      return NextResponse.json(
        { error: "Failed to start pipeline" },
        { status: 500 },
      );
    }
  } catch {
    return NextResponse.json(
      { error: "Failed to start pipeline" },
      { status: 500 },
    );
  }
}
