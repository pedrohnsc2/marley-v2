import { NextResponse, type NextRequest } from "next/server";
import createIntlMiddleware from "next-intl/middleware";
import { routing } from "./i18n/routing";
import { updateSession } from "./lib/supabase/middleware";

const intlMiddleware = createIntlMiddleware(routing);

const AUTH_ROUTES = ["/login", "/signup"];

function stripLocale(pathname: string): string {
  // Match locale prefix like /pt-BR/ or /en/ (must be followed by / or end)
  return pathname.replace(/^\/[a-z]{2}(-[A-Z]{2})?(\/|$)/, "/") || "/";
}

export async function middleware(request: NextRequest) {
  const { pathname } = request.nextUrl;

  // Skip auth for API routes, static files, and auth callback
  if (
    pathname.startsWith("/api") ||
    pathname.startsWith("/_next") ||
    pathname.startsWith("/auth")
  ) {
    return NextResponse.next();
  }

  // Refresh Supabase session and get user
  const { supabaseResponse, user } = await updateSession(request);

  const pathWithoutLocale = stripLocale(pathname);

  // Authenticated user on auth page → redirect to /runs
  if (user && AUTH_ROUTES.includes(pathWithoutLocale)) {
    const url = request.nextUrl.clone();
    url.pathname = "/runs";
    return NextResponse.redirect(url);
  }

  // Unauthenticated user on protected page → redirect to /login
  // Public pages: /, /login, /signup
  const isPublic =
    pathWithoutLocale === "/" ||
    AUTH_ROUTES.includes(pathWithoutLocale);

  if (!user && !isPublic) {
    const url = request.nextUrl.clone();
    url.pathname = "/login";
    return NextResponse.redirect(url);
  }

  // Run intl middleware (handles locale detection, redirects, rewrites)
  const intlResponse = intlMiddleware(request);

  // Copy Supabase auth cookies into the intl response so the browser
  // receives both the locale-related headers AND the refreshed session.
  supabaseResponse.cookies.getAll().forEach((cookie) => {
    intlResponse.cookies.set(cookie.name, cookie.value);
  });

  return intlResponse;
}

export const config = {
  matcher: ["/((?!_next/static|_next/image|favicon.ico|.*\\..*).*)"],
};
