import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";

const { mockGetUser, mockCreateServerClient, mockNextResponse } = vi.hoisted(
  () => {
    const mockGetUser = vi.fn(() =>
      Promise.resolve({ data: { user: null }, error: null })
    );

    const mockCreateServerClient = vi.fn(() => ({
      auth: { getUser: mockGetUser },
    }));

    // Build a fake NextResponse that has a .cookies helper
    const createFakeResponse = () => {
      const cookieMap = new Map<string, { value: string; options?: unknown }>();
      return {
        cookies: {
          set: (name: string, value: string, options?: unknown) => {
            cookieMap.set(name, { value, options });
          },
          get: (name: string) => cookieMap.get(name),
          getAll: () =>
            Array.from(cookieMap.entries()).map(([name, { value }]) => ({
              name,
              value,
            })),
        },
        headers: new Headers(),
      };
    };

    const mockNextResponse = {
      next: vi.fn(() => createFakeResponse()),
      _createFakeResponse: createFakeResponse,
    };

    return { mockGetUser, mockCreateServerClient, mockNextResponse };
  }
);

vi.mock("@supabase/ssr", () => ({
  createServerClient: mockCreateServerClient,
}));

// Mock next/server to provide both NextRequest and our controlled NextResponse
vi.mock("next/server", async () => {
  const actual = await vi.importActual<typeof import("next/server")>(
    "next/server"
  );
  return {
    ...actual,
    NextResponse: {
      ...actual.NextResponse,
      next: mockNextResponse.next,
    },
  };
});

import { updateSession } from "@/lib/supabase/middleware";
import { NextRequest } from "next/server";

function createMockRequest(url = "http://localhost:3000/"): NextRequest {
  return new NextRequest(url);
}

describe("supabase/middleware - updateSession", () => {
  const originalEnv = process.env;

  beforeEach(() => {
    vi.clearAllMocks();
    process.env = {
      ...originalEnv,
      NEXT_PUBLIC_SUPABASE_URL: "https://test.supabase.co",
      NEXT_PUBLIC_SUPABASE_ANON_KEY: "test-anon-key-789",
    };
  });

  afterEach(() => {
    process.env = originalEnv;
  });

  it("creates a Supabase client with correct env vars", async () => {
    const request = createMockRequest();

    await updateSession(request);

    expect(mockCreateServerClient).toHaveBeenCalledOnce();
    expect(mockCreateServerClient).toHaveBeenCalledWith(
      "https://test.supabase.co",
      "test-anon-key-789",
      expect.objectContaining({
        cookies: expect.objectContaining({
          getAll: expect.any(Function),
          setAll: expect.any(Function),
        }),
      })
    );
  });

  it("calls getUser to refresh the session", async () => {
    const request = createMockRequest();

    await updateSession(request);

    expect(mockGetUser).toHaveBeenCalledOnce();
  });

  it("returns supabaseResponse and user", async () => {
    const request = createMockRequest();

    const result = await updateSession(request);

    expect(result).toHaveProperty("supabaseResponse");
    expect(result).toHaveProperty("user");
    expect(result.supabaseResponse.cookies).toBeDefined();
  });

  it("returns null user when not authenticated", async () => {
    const request = createMockRequest();

    const { user } = await updateSession(request);

    expect(user).toBeNull();
  });

  it("returns user when authenticated", async () => {
    const fakeUser = { id: "user-123", email: "test@example.com" };
    mockGetUser.mockResolvedValueOnce({
      data: { user: fakeUser },
      error: null,
    });

    const request = createMockRequest();

    const { user } = await updateSession(request);

    expect(user).toEqual(fakeUser);
  });

  it("cookie getAll reads from the request cookies", async () => {
    const request = createMockRequest();
    request.cookies.set("sb-auth", "token-value");

    await updateSession(request);

    const cookiesConfig = mockCreateServerClient.mock.calls[0][2].cookies;
    const allCookies = cookiesConfig.getAll();

    // Should include the cookie we set on the request
    expect(allCookies).toEqual(
      expect.arrayContaining([
        expect.objectContaining({ name: "sb-auth", value: "token-value" }),
      ])
    );
  });

  it("cookie setAll updates the request cookies", async () => {
    const request = createMockRequest();

    await updateSession(request);

    const cookiesConfig = mockCreateServerClient.mock.calls[0][2].cookies;

    // Simulate Supabase setting cookies during getUser
    cookiesConfig.setAll([
      {
        name: "sb-access-token",
        value: "new-token",
        options: { path: "/", httpOnly: true },
      },
    ]);

    // The request cookie should have been updated
    expect(request.cookies.get("sb-access-token")?.value).toBe("new-token");
  });

  it("cookie setAll creates a new NextResponse.next and sets cookies on it", async () => {
    const request = createMockRequest();

    await updateSession(request);

    const cookiesConfig = mockCreateServerClient.mock.calls[0][2].cookies;

    // setAll should trigger a new NextResponse.next() call
    const callCountBefore = mockNextResponse.next.mock.calls.length;
    cookiesConfig.setAll([
      {
        name: "sb-token",
        value: "refreshed",
        options: { path: "/" },
      },
    ]);
    const callCountAfter = mockNextResponse.next.mock.calls.length;

    // setAll creates a new NextResponse.next({ request })
    expect(callCountAfter).toBe(callCountBefore + 1);
  });
});
