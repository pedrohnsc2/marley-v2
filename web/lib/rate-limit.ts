import { NextRequest, NextResponse } from "next/server";

interface RateLimitEntry {
  count: number;
  resetAt: number;
}

const store = new Map<string, RateLimitEntry>();

// Clean up expired entries every 60s
if (typeof setInterval !== "undefined") {
  setInterval(() => {
    const now = Date.now();
    for (const [key, entry] of store) {
      if (entry.resetAt <= now) store.delete(key);
    }
  }, 60_000);
}

function getClientIp(request: NextRequest): string {
  return (
    request.headers.get("x-forwarded-for")?.split(",")[0]?.trim() ??
    request.headers.get("x-real-ip") ??
    "unknown"
  );
}

/**
 * Rate limit check. Returns null if allowed, or 429 response if exceeded.
 */
export function rateLimit(
  request: NextRequest,
  {
    maxRequests = 60,
    windowMs = 60_000,
  }: { maxRequests?: number; windowMs?: number } = {},
): NextResponse | null {
  const ip = getClientIp(request);
  const key = `${ip}:${request.nextUrl.pathname}`;
  const now = Date.now();

  const entry = store.get(key);
  if (!entry || entry.resetAt <= now) {
    store.set(key, { count: 1, resetAt: now + windowMs });
    return null;
  }

  entry.count++;
  if (entry.count > maxRequests) {
    return NextResponse.json(
      { error: "Too many requests" },
      {
        status: 429,
        headers: {
          "Retry-After": String(Math.ceil((entry.resetAt - now) / 1000)),
        },
      },
    );
  }
  return null;
}

// Track active pipeline runs for concurrency control
let activeRuns = 0;
const MAX_CONCURRENT_RUNS = 3;

export function acquireRunSlot(): boolean {
  if (activeRuns >= MAX_CONCURRENT_RUNS) return false;
  activeRuns++;
  return true;
}

export function releaseRunSlot(): void {
  if (activeRuns > 0) activeRuns--;
}
