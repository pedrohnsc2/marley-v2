import { createClient } from "@/lib/supabase/server";
import { NextRequest, NextResponse } from "next/server";
import { timingSafeEqual } from "crypto";

interface AuthResult {
  response: NextResponse | null;
  userId: string | null;
}

/**
 * Validates the request against API key or Supabase session.
 *
 * Priority:
 *  1. Service-to-service: X-API-Key header matches API_SECRET_KEY (constant-time)
 *  2. User session: valid Supabase auth cookie
 *  3. Dev fallback: if API_SECRET_KEY is unset AND Supabase fails, allow through
 *
 * Returns { response: null } when auth succeeds and the route should proceed.
 * Returns { response: NextResponse(401) } when auth fails.
 */
export async function requireAuth(
  request: NextRequest
): Promise<AuthResult> {
  const apiKey = process.env.API_SECRET_KEY;

  // 1. Service-to-service: constant-time API key check
  const providedKey = request.headers.get("x-api-key");
  if (apiKey && providedKey) {
    try {
      const a = Buffer.from(apiKey, "utf-8");
      const b = Buffer.from(providedKey, "utf-8");
      if (a.length === b.length && timingSafeEqual(a, b)) {
        return { response: null, userId: null };
      }
    } catch {
      // length mismatch or encoding issue — fall through
    }
  }

  // 2. User auth: check Supabase session via cookies
  try {
    const supabase = await createClient();
    const {
      data: { user },
      error,
    } = await supabase.auth.getUser();

    if (!error && user) {
      return { response: null, userId: user.id };
    }
  } catch {
    // Supabase client creation failed (e.g., no URL configured)
  }

  // 3. Dev fallback: skip auth only when API_SECRET_KEY is not configured
  if (!apiKey) {
    return { response: null, userId: null };
  }

  // Production: both API key and session failed
  return {
    response: NextResponse.json(
      { error: "Unauthorized" },
      { status: 401 }
    ),
    userId: null,
  };
}

/**
 * Same as requireAuth but for write operations (POST/PUT/DELETE).
 */
export async function requireWriteAuth(
  request: NextRequest
): Promise<AuthResult> {
  return requireAuth(request);
}
