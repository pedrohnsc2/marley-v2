import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { NextRequest } from "next/server";

// Mock the Supabase server client before importing the module under test
const mockGetUser = vi.fn();
const mockSupabaseClient = {
  auth: {
    getUser: mockGetUser,
  },
};

vi.mock("@/lib/supabase/server", () => ({
  createClient: vi.fn(() => Promise.resolve(mockSupabaseClient)),
}));

// Import after mocks are set up
import { requireAuth, requireWriteAuth } from "@/lib/api-auth";

function createMockRequest(headers: Record<string, string> = {}): NextRequest {
  const reqHeaders = new Headers(headers);
  return new NextRequest("http://localhost:3000/api/test", {
    headers: reqHeaders,
  });
}

describe("requireAuth", () => {
  const originalEnv = process.env;

  beforeEach(() => {
    vi.clearAllMocks();
    // Reset mock implementations (clearAllMocks only clears call history)
    mockGetUser.mockReset();
    process.env = { ...originalEnv };
    delete process.env.API_SECRET_KEY;
  });

  afterEach(() => {
    process.env = originalEnv;
  });

  it("falls back to dev mode when no API key and no Supabase session", async () => {
    mockGetUser.mockResolvedValue({
      data: { user: null },
      error: null,
    });
    const request = createMockRequest();

    const result = await requireAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBeNull();
  });

  it("returns userId from Supabase even in dev mode (no API key)", async () => {
    mockGetUser.mockResolvedValue({
      data: { user: { id: "dev-user-1" } },
      error: null,
    });
    const request = createMockRequest();

    const result = await requireAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBe("dev-user-1");
  });

  it("succeeds with valid API key (constant-time comparison)", async () => {
    process.env.API_SECRET_KEY = "test-secret-key";
    const request = createMockRequest({ "x-api-key": "test-secret-key" });

    const result = await requireAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBeNull();
    // Supabase should NOT be called when API key matches
    expect(mockGetUser).not.toHaveBeenCalled();
  });

  it("rejects mismatched API key and falls through to Supabase", async () => {
    process.env.API_SECRET_KEY = "test-secret-key";
    mockGetUser.mockResolvedValue({
      data: { user: null },
      error: { message: "No session" },
    });
    const request = createMockRequest({ "x-api-key": "wrong-key" });

    const result = await requireAuth(request);

    expect(result.response).not.toBeNull();
    expect(result.response!.status).toBe(401);
    expect(result.userId).toBeNull();
  });

  it("returns 401 when API key is missing and no Supabase session", async () => {
    process.env.API_SECRET_KEY = "test-secret-key";
    mockGetUser.mockResolvedValue({
      data: { user: null },
      error: null,
    });
    const request = createMockRequest();

    const result = await requireAuth(request);

    expect(result.response).not.toBeNull();
    expect(result.response!.status).toBe(401);
    expect(result.userId).toBeNull();
  });

  it("returns 401 when Supabase returns error even with user data", async () => {
    process.env.API_SECRET_KEY = "test-secret-key";
    mockGetUser.mockResolvedValue({
      data: { user: { id: "user-123" } },
      error: { message: "Token expired" },
    });
    const request = createMockRequest();

    const result = await requireAuth(request);

    expect(result.response).not.toBeNull();
    expect(result.response!.status).toBe(401);
    expect(result.userId).toBeNull();
  });

  it("succeeds with valid Supabase session and returns userId", async () => {
    process.env.API_SECRET_KEY = "test-secret-key";
    mockGetUser.mockResolvedValue({
      data: { user: { id: "user-abc-123" } },
      error: null,
    });
    const request = createMockRequest();

    const result = await requireAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBe("user-abc-123");
  });

  it("handles Supabase client creation failure in dev mode", async () => {
    mockGetUser.mockRejectedValue(new Error("Connection refused"));
    const request = createMockRequest();

    const result = await requireAuth(request);

    // Dev mode fallback: allow through
    expect(result.response).toBeNull();
    expect(result.userId).toBeNull();
  });

  it("returns 401 when Supabase fails in production", async () => {
    process.env.API_SECRET_KEY = "test-secret-key";
    mockGetUser.mockRejectedValue(new Error("Connection refused"));
    const request = createMockRequest();

    const result = await requireAuth(request);

    expect(result.response).not.toBeNull();
    expect(result.response!.status).toBe(401);
  });
});

describe("requireWriteAuth", () => {
  const originalEnv = process.env;

  beforeEach(() => {
    vi.clearAllMocks();
    mockGetUser.mockReset();
    process.env = { ...originalEnv };
    delete process.env.API_SECRET_KEY;
  });

  afterEach(() => {
    process.env = originalEnv;
  });

  it("delegates to requireAuth (dev mode, no session)", async () => {
    mockGetUser.mockResolvedValue({
      data: { user: null },
      error: null,
    });
    const request = createMockRequest();

    const result = await requireWriteAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBeNull();
  });

  it("delegates to requireAuth (valid API key)", async () => {
    process.env.API_SECRET_KEY = "write-secret";
    const request = createMockRequest({ "x-api-key": "write-secret" });

    const result = await requireWriteAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBeNull();
  });

  it("delegates to requireAuth (no key, valid session)", async () => {
    process.env.API_SECRET_KEY = "write-secret";
    mockGetUser.mockResolvedValue({
      data: { user: { id: "writer-456" } },
      error: null,
    });
    const request = createMockRequest();

    const result = await requireWriteAuth(request);

    expect(result.response).toBeNull();
    expect(result.userId).toBe("writer-456");
  });
});
