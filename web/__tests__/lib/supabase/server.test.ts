import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";

const { mockCookieStore, mockCreateServerClient } = vi.hoisted(() => {
  const mockCookieStore = {
    getAll: vi.fn(() => [{ name: "sb-token", value: "abc123" }]),
    set: vi.fn(),
  };

  const mockCreateServerClient = vi.fn(() => ({
    auth: { getUser: vi.fn() },
  }));

  return { mockCookieStore, mockCreateServerClient };
});

vi.mock("next/headers", () => ({
  cookies: vi.fn(() => Promise.resolve(mockCookieStore)),
}));

vi.mock("@supabase/ssr", () => ({
  createServerClient: mockCreateServerClient,
}));

import { createClient } from "@/lib/supabase/server";

describe("supabase/server - createClient", () => {
  const originalEnv = process.env;

  beforeEach(() => {
    vi.clearAllMocks();
    process.env = {
      ...originalEnv,
      NEXT_PUBLIC_SUPABASE_URL: "https://test.supabase.co",
      NEXT_PUBLIC_SUPABASE_ANON_KEY: "test-anon-key-456",
    };
  });

  afterEach(() => {
    process.env = originalEnv;
  });

  it("calls createServerClient with correct env vars", async () => {
    await createClient();

    expect(mockCreateServerClient).toHaveBeenCalledOnce();
    expect(mockCreateServerClient).toHaveBeenCalledWith(
      "https://test.supabase.co",
      "test-anon-key-456",
      expect.objectContaining({
        cookies: expect.objectContaining({
          getAll: expect.any(Function),
          setAll: expect.any(Function),
        }),
      })
    );
  });

  it("cookie getAll delegates to the cookie store", async () => {
    await createClient();

    // Extract the cookies config passed to createServerClient
    const cookiesConfig = mockCreateServerClient.mock.calls[0][2].cookies;
    const result = cookiesConfig.getAll();

    expect(mockCookieStore.getAll).toHaveBeenCalled();
    expect(result).toEqual([{ name: "sb-token", value: "abc123" }]);
  });

  it("cookie setAll delegates to the cookie store", async () => {
    await createClient();

    const cookiesConfig = mockCreateServerClient.mock.calls[0][2].cookies;
    const cookiesToSet = [
      { name: "sb-token", value: "new-value", options: { path: "/" } },
      { name: "sb-refresh", value: "refresh-val", options: { httpOnly: true } },
    ];

    cookiesConfig.setAll(cookiesToSet);

    expect(mockCookieStore.set).toHaveBeenCalledTimes(2);
    expect(mockCookieStore.set).toHaveBeenCalledWith(
      "sb-token",
      "new-value",
      { path: "/" }
    );
    expect(mockCookieStore.set).toHaveBeenCalledWith(
      "sb-refresh",
      "refresh-val",
      { httpOnly: true }
    );
  });

  it("returns the client instance", async () => {
    const mockClient = { auth: { getUser: vi.fn() }, from: vi.fn() };
    mockCreateServerClient.mockReturnValue(mockClient);

    const client = await createClient();

    expect(client).toBe(mockClient);
  });
});
