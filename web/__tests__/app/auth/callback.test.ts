import { describe, it, expect, vi, beforeEach } from "vitest";

const mockExchangeCodeForSession = vi.fn();
const mockSupabaseClient = {
  auth: {
    exchangeCodeForSession: mockExchangeCodeForSession,
  },
};

vi.mock("@/lib/supabase/server", () => ({
  createClient: vi.fn(() => Promise.resolve(mockSupabaseClient)),
}));

import { GET } from "@/app/auth/callback/route";

function createCallbackRequest(
  params: Record<string, string> = {}
): Request {
  const url = new URL("http://localhost:3000/auth/callback");
  for (const [key, value] of Object.entries(params)) {
    url.searchParams.set(key, value);
  }
  return new Request(url.toString());
}

describe("auth/callback route - GET", () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it("exchanges code for session and redirects to / on success", async () => {
    mockExchangeCodeForSession.mockResolvedValue({ error: null });
    const request = createCallbackRequest({ code: "valid-auth-code" });

    const response = await GET(request);

    expect(mockExchangeCodeForSession).toHaveBeenCalledWith("valid-auth-code");
    expect(response.status).toBe(307); // NextResponse.redirect uses 307
    expect(response.headers.get("location")).toBe("http://localhost:3000/");
  });

  it("redirects to custom next path after successful exchange", async () => {
    mockExchangeCodeForSession.mockResolvedValue({ error: null });
    const request = createCallbackRequest({
      code: "valid-code",
      next: "/dashboard",
    });

    const response = await GET(request);

    expect(response.headers.get("location")).toBe(
      "http://localhost:3000/dashboard"
    );
  });

  it("redirects to login with error when code is missing", async () => {
    const request = createCallbackRequest();

    const response = await GET(request);

    expect(mockExchangeCodeForSession).not.toHaveBeenCalled();
    expect(response.status).toBe(307);
    expect(response.headers.get("location")).toBe(
      "http://localhost:3000/login?error=auth_callback_failed"
    );
  });

  it("redirects to login with error when session exchange fails", async () => {
    mockExchangeCodeForSession.mockResolvedValue({
      error: { message: "Invalid code" },
    });
    const request = createCallbackRequest({ code: "expired-code" });

    const response = await GET(request);

    expect(mockExchangeCodeForSession).toHaveBeenCalledWith("expired-code");
    expect(response.status).toBe(307);
    expect(response.headers.get("location")).toBe(
      "http://localhost:3000/login?error=auth_callback_failed"
    );
  });

  it("does not call createClient when code is missing", async () => {
    const { createClient } = await import("@/lib/supabase/server");
    vi.mocked(createClient).mockClear();

    const request = createCallbackRequest();
    await GET(request);

    expect(createClient).not.toHaveBeenCalled();
  });

  // Open redirect prevention (C3 security fix)
  it("rejects protocol-relative redirect (//evil.com)", async () => {
    mockExchangeCodeForSession.mockResolvedValue({ error: null });
    const request = createCallbackRequest({
      code: "valid-code",
      next: "//evil.com",
    });

    const response = await GET(request);

    // Should redirect to / instead of //evil.com
    expect(response.headers.get("location")).toBe("http://localhost:3000/");
  });

  it("rejects absolute URL redirect", async () => {
    mockExchangeCodeForSession.mockResolvedValue({ error: null });
    const request = createCallbackRequest({
      code: "valid-code",
      next: "https://evil.com/steal",
    });

    const response = await GET(request);

    // Should redirect to / since it doesn't start with /
    expect(response.headers.get("location")).toBe("http://localhost:3000/");
  });
});
