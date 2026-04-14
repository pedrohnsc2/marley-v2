import { NextRequest, NextResponse } from "next/server";
import { listRuns } from "@/lib/data-loader";
import { rateLimit } from "@/lib/rate-limit";
import type { RunMetadata } from "@/lib/types/run";

export async function GET(request: NextRequest) {
  const limited = rateLimit(request);
  if (limited) return limited;

  const { searchParams } = new URL(request.url);
  const pipeline = searchParams.get("pipeline") ?? undefined;

  try {
    const runs = listRuns(pipeline) as RunMetadata[];
    return NextResponse.json(runs);
  } catch {
    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500 },
    );
  }
}
