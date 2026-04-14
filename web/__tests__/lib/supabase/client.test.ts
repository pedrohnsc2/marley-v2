import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";

const { mockCreateBrowserClient } = vi.hoisted(() => ({
  mockCreateBrowserClient: vi.fn(() => ({ from: vi.fn() })),
}));

vi.mock("@supabase/ssr", () => ({
  createBrowserClient: mockCreateBrowserClient,
}));

import { createClient } from "@/lib/supabase/client";

describe("supabase/client - createClient", () => {
  const originalEnv = process.env;

  beforeEach(() => {
    vi.clearAllMocks();
    process.env = {
      ...originalEnv,
      NEXT_PUBLIC_SUPABASE_URL: "https://test.supabase.co",
      NEXT_PUBLIC_SUPABASE_ANON_KEY: "test-anon-key-123",
    };
  });

  afterEach(() => {
    process.env = originalEnv;
  });

  it("calls createBrowserClient with correct env vars", () => {
    createClient();

    expect(mockCreateBrowserClient).toHaveBeenCalledOnce();
    expect(mockCreateBrowserClient).toHaveBeenCalledWith(
      "https://test.supabase.co",
      "test-anon-key-123"
    );
  });

  it("returns the client instance from createBrowserClient", () => {
    const mockClient = { from: vi.fn(), auth: {} };
    mockCreateBrowserClient.mockReturnValue(mockClient);

    const client = createClient();

    expect(client).toBe(mockClient);
  });
});
