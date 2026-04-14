import { NextRequest, NextResponse } from "next/server";
import { loadRunMetadata } from "@/lib/data-loader";
import { requireAuth } from "@/lib/api-auth";
import { rateLimit } from "@/lib/rate-limit";
import type { RunMetadata } from "@/lib/types/run";

export async function GET(
  request: NextRequest,
  { params }: { params: { runId: string } },
) {
  const { response: authErr } = await requireAuth(request);
  if (authErr) return authErr;

  const limited = rateLimit(request);
  if (limited) return limited;

  const { runId } = params;

  if (!runId) {
    return NextResponse.json(
      { error: "Missing runId parameter" },
      { status: 400 },
    );
  }

  // Security: reject path traversal attempts
  if (runId.includes("..") || runId.includes("/") || runId.includes("\\")) {
    return NextResponse.json(
      { error: "Invalid runId" },
      { status: 400 },
    );
  }

  try {
    const run = loadRunMetadata(runId) as RunMetadata;
    return NextResponse.json(run);
  } catch (err) {
    const message = err instanceof Error ? err.message : "";

    if (message.includes("ENOENT") || message.includes("no such file")) {
      return NextResponse.json(
        { error: "Run not found" },
        { status: 404 },
      );
    }

    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500 },
    );
  }
}
