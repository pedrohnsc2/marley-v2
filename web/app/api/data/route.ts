import { NextRequest, NextResponse } from "next/server";
import fs from "fs";
import path from "path";

const RESULTS_DIR = path.join(process.cwd(), "..", "results");

function parseCsv(content: string): Record<string, string>[] {
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

export async function GET(request: NextRequest) {
  const { searchParams } = new URL(request.url);
  const file = searchParams.get("file");

  if (!file) {
    return NextResponse.json(
      { error: "Missing 'file' query parameter" },
      { status: 400 },
    );
  }

  // Security: prevent path traversal
  const normalized = path.normalize(file);
  if (normalized.includes("..") || path.isAbsolute(normalized)) {
    return NextResponse.json(
      { error: "Invalid file path" },
      { status: 400 },
    );
  }

  const filePath = path.join(RESULTS_DIR, normalized);

  // Ensure resolved path is within RESULTS_DIR
  if (!filePath.startsWith(RESULTS_DIR)) {
    return NextResponse.json(
      { error: "Access denied" },
      { status: 403 },
    );
  }

  if (!fs.existsSync(filePath)) {
    return NextResponse.json(
      { error: "File not found" },
      { status: 404 },
    );
  }

  try {
    const content = fs.readFileSync(filePath, "utf-8");
    const ext = path.extname(filePath).toLowerCase();

    if (ext === ".csv") {
      return NextResponse.json(parseCsv(content));
    }

    if (ext === ".json") {
      return NextResponse.json(JSON.parse(content));
    }

    // Return raw content for markdown and other text files
    return new NextResponse(content, {
      headers: { "Content-Type": "text/plain; charset=utf-8" },
    });
  } catch (err) {
    const message = err instanceof Error ? err.message : "Unknown error";
    return NextResponse.json(
      { error: `Failed to read file: ${message}` },
      { status: 500 },
    );
  }
}
