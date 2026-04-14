import { NextRequest, NextResponse } from "next/server";

const API_KEY = process.env.API_SECRET_KEY;

/**
 * Validates API key from X-API-Key header.
 * Returns null if valid, or a 401 NextResponse if invalid.
 * If API_SECRET_KEY env var is not set, auth is skipped (dev mode).
 */
export function requireAuth(request: NextRequest): NextResponse | null {
  // Skip auth if no key configured (local dev)
  if (!API_KEY) return null;

  const provided = request.headers.get("x-api-key");
  if (provided !== API_KEY) {
    return NextResponse.json({ error: "Unauthorized" }, { status: 401 });
  }
  return null;
}

/**
 * Same as requireAuth but for write operations (POST/PUT/DELETE).
 * Always enforces auth when API_SECRET_KEY is set.
 */
export function requireWriteAuth(request: NextRequest): NextResponse | null {
  return requireAuth(request);
}
