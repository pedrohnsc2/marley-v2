"use client";

import { useState } from "react";
import { createClient } from "@/lib/supabase/client";
import { Link, useRouter } from "@/i18n/routing";
import LoginParticles from "@/components/login-particles";

type AuthMode = "password" | "magic-link";

export default function LoginPage() {
  const router = useRouter();
  const supabase = createClient();

  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [mode, setMode] = useState<AuthMode>("password");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [magicLinkSent, setMagicLinkSent] = useState(false);

  async function handlePasswordLogin(e: React.FormEvent) {
    e.preventDefault();
    setLoading(true);
    setError(null);

    const { error: signInError } = await supabase.auth.signInWithPassword({
      email,
      password,
    });

    if (signInError) {
      setError("Invalid email or password.");
      setLoading(false);
      return;
    }

    router.push("/runs");
  }

  async function handleMagicLink(e: React.FormEvent) {
    e.preventDefault();
    setLoading(true);
    setError(null);

    const { error: otpError } = await supabase.auth.signInWithOtp({
      email,
      options: {
        emailRedirectTo: `${window.location.origin}/auth/callback?next=/runs`,
      },
    });

    if (otpError) {
      setError("Could not send magic link. Please try again.");
      setLoading(false);
      return;
    }

    setMagicLinkSent(true);
    setLoading(false);
  }

  return (
    <div className="login-bg relative flex min-h-screen flex-col items-center justify-center px-4">
      <LoginParticles />

      {/* ── Card ── */}
      <div className="login-card relative z-10 w-full max-w-[420px] px-8 py-10 sm:px-10">
        <div className="login-stagger">
          {/* Wordmark */}
          <div className="mb-8 text-center">
            <h1
              className="text-3xl font-bold"
              style={{ color: "var(--app-text)", letterSpacing: "-0.02em" }}
            >
              Marley
            </h1>
            <p
              className="mt-1.5 text-[11px] font-medium uppercase"
              style={{ color: "var(--app-text-3)", letterSpacing: "0.1em" }}
            >
              Reverse Vaccinology Platform
            </p>
          </div>

          {/* Error */}
          {error && (
            <div
              className="mb-4 rounded-lg px-4 py-3 text-sm"
              style={{
                backgroundColor: "var(--app-badge-warn-bg)",
                color: "var(--app-badge-warn-tx)",
              }}
              role="alert"
              data-testid="login-error"
            >
              {error}
            </div>
          )}

          {/* Magic link sent */}
          {magicLinkSent ? (
            <div className="text-center">
              <div
                className="mb-6 rounded-xl px-4 py-8 text-sm"
                style={{
                  backgroundColor: "var(--app-badge-success-bg)",
                  color: "var(--app-badge-success-tx)",
                }}
                data-testid="magic-link-sent"
              >
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.5} className="mx-auto mb-3 h-10 w-10">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M21.75 6.75v10.5a2.25 2.25 0 01-2.25 2.25h-15a2.25 2.25 0 01-2.25-2.25V6.75m19.5 0A2.25 2.25 0 0019.5 4.5h-15a2.25 2.25 0 00-2.25 2.25m19.5 0v.243a2.25 2.25 0 01-1.07 1.916l-7.5 4.615a2.25 2.25 0 01-2.36 0L3.32 8.91a2.25 2.25 0 01-1.07-1.916V6.75" />
                </svg>
                <p className="font-medium">Check your email</p>
                <p className="mt-1 opacity-80">
                  We sent a login link to <strong>{email}</strong>
                </p>
              </div>
              <button
                onClick={() => { setMagicLinkSent(false); setMode("password"); }}
                className="text-sm font-medium underline"
                style={{ color: "var(--app-primary)" }}
              >
                Back to login
              </button>
            </div>
          ) : (
            <>
              {/* Tabs */}
              <div
                className="mb-6 flex rounded-xl p-1"
                style={{ backgroundColor: "var(--app-surface-2)" }}
              >
                {(["password", "magic-link"] as const).map((tab) => (
                  <button
                    key={tab}
                    onClick={() => setMode(tab)}
                    className="flex-1 rounded-lg py-2 text-[13px] font-medium transition-all"
                    style={{
                      backgroundColor: mode === tab ? "var(--app-surface)" : "transparent",
                      color: mode === tab ? "var(--app-text)" : "var(--app-text-3)",
                      boxShadow: mode === tab ? "var(--app-card-shadow)" : "none",
                    }}
                    data-testid={`tab-${tab}`}
                  >
                    {tab === "password" ? "Password" : "Magic Link"}
                  </button>
                ))}
              </div>

              {/* Password form */}
              {mode === "password" && (
                <form onSubmit={handlePasswordLogin} className="space-y-4">
                  <div>
                    <label htmlFor="email" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                      Email
                    </label>
                    <input
                      id="email" type="email" required value={email}
                      onChange={(e) => setEmail(e.target.value)}
                      placeholder="you@university.edu"
                      className="login-input app-input w-full rounded-lg px-4 py-2.5 text-sm outline-none"
                      data-testid="input-email"
                    />
                  </div>
                  <div>
                    <label htmlFor="password" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                      Password
                    </label>
                    <input
                      id="password" type="password" required value={password}
                      onChange={(e) => setPassword(e.target.value)}
                      placeholder="Your password"
                      className="login-input app-input w-full rounded-lg px-4 py-2.5 text-sm outline-none"
                      data-testid="input-password"
                    />
                  </div>
                  <button
                    type="submit" disabled={loading}
                    className="login-btn app-btn-primary w-full rounded-lg py-2.5 text-sm font-semibold disabled:opacity-50"
                    style={{ letterSpacing: "0.02em" }}
                    data-testid="btn-login"
                  >
                    {loading ? "Signing in..." : "Sign In"}
                  </button>
                </form>
              )}

              {/* Magic link form */}
              {mode === "magic-link" && (
                <form onSubmit={handleMagicLink} className="space-y-4">
                  <div>
                    <label htmlFor="magic-email" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                      Email
                    </label>
                    <input
                      id="magic-email" type="email" required value={email}
                      onChange={(e) => setEmail(e.target.value)}
                      placeholder="you@university.edu"
                      className="login-input app-input w-full rounded-lg px-4 py-2.5 text-sm outline-none"
                      data-testid="input-magic-email"
                    />
                  </div>
                  <p className="text-xs" style={{ color: "var(--app-text-3)" }}>
                    We&apos;ll email you a link to sign in — no password needed.
                  </p>
                  <button
                    type="submit" disabled={loading}
                    className="login-btn app-btn-primary w-full rounded-lg py-2.5 text-sm font-semibold disabled:opacity-50"
                    style={{ letterSpacing: "0.02em" }}
                    data-testid="btn-magic-link"
                  >
                    {loading ? "Sending link..." : "Send Magic Link"}
                  </button>
                </form>
              )}

              {/* Divider + Signup */}
              <div className="my-6 flex items-center gap-3">
                <div className="h-px flex-1" style={{ backgroundColor: "var(--app-border)" }} />
                <span className="text-[11px]" style={{ color: "var(--app-text-3)" }}>New to Marley?</span>
                <div className="h-px flex-1" style={{ backgroundColor: "var(--app-border)" }} />
              </div>

              <Link
                href="/signup"
                className="login-btn block w-full rounded-lg border py-2.5 text-center text-sm font-medium"
                style={{
                  borderColor: "var(--app-border)",
                  color: "var(--app-text)",
                  backgroundColor: "transparent",
                }}
                data-testid="link-signup"
              >
                Create an account
              </Link>
            </>
          )}
        </div>
      </div>

      {/* ── Institution bar ── */}
      <div
        className="relative z-10 mt-8 pb-6 text-center"
        style={{ opacity: 0.4 }}
      >
        <p
          className="text-[11px] font-medium uppercase"
          style={{ color: "var(--app-text-3)", letterSpacing: "0.1em" }}
        >
          UFMG &nbsp;&middot;&nbsp; Fiocruz &nbsp;&middot;&nbsp; CTVacinas
        </p>
      </div>
    </div>
  );
}
