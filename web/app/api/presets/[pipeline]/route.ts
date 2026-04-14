import { NextRequest, NextResponse } from "next/server";
import { listPresets, loadPreset, splitPreset } from "@/lib/presets";
import { requireAuth } from "@/lib/api-auth";
import { rateLimit } from "@/lib/rate-limit";

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ pipeline: string }> },
) {
  const { response: authErr } = await requireAuth(request);
  if (authErr) return authErr;

  const limited = rateLimit(request);
  if (limited) return limited;

  try {
    const { pipeline } = await params;

    const url = new URL(request.url);
    const name = url.searchParams.get("name");

    // Single preset: return { meta, params } with _meta separated
    if (name) {
      const raw = loadPreset(pipeline, name);
      if (!raw) {
        return NextResponse.json(
          { error: "Preset not found" },
          { status: 404 },
        );
      }
      const { meta, params: presetParams } = splitPreset(raw);
      return NextResponse.json({ meta, params: presetParams });
    }

    // List presets: return enriched array with meta for each preset
    const names = listPresets(pipeline);
    const enriched = names.map((presetName) => {
      const raw = loadPreset(pipeline, presetName);
      const meta = raw ? splitPreset(raw).meta : null;
      return { name: presetName, meta };
    });
    return NextResponse.json(enriched);
  } catch {
    return NextResponse.json(
      { error: "Internal server error" },
      { status: 500 },
    );
  }
}
