import { NextRequest, NextResponse } from "next/server";
import { listPresets, loadPreset } from "@/lib/presets";
import { rateLimit } from "@/lib/rate-limit";

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ pipeline: string }> },
) {
  const limited = rateLimit(request);
  if (limited) return limited;

  try {
    const { pipeline } = await params;

    const url = new URL(request.url);
    const name = url.searchParams.get("name");

    if (name) {
      const preset = loadPreset(pipeline, name);
      if (!preset) {
        return NextResponse.json(
          { error: "Preset not found" },
          { status: 404 },
        );
      }
      return NextResponse.json(preset);
    }

    const presets = listPresets(pipeline);
    return NextResponse.json(presets);
  } catch {
    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500 },
    );
  }
}
