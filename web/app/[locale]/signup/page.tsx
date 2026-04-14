"use client";

import { useState } from "react";
import { createClient } from "@/lib/supabase/client";
import { Link } from "@/i18n/routing";
import LoginBackground from "@/components/login-background";

export default function SignupPage() {
  const supabase = createClient();

  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [fullName, setFullName] = useState("");
  const [institution, setInstitution] = useState("");
  const [orcid, setOrcid] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [success, setSuccess] = useState(false);

  async function handleSignup(e: React.FormEvent) {
    e.preventDefault();
    setLoading(true);
    setError(null);

    if (password.length < 8) {
      setError("Password must be at least 8 characters.");
      setLoading(false);
      return;
    }

    const { error: signUpError } = await supabase.auth.signUp({
      email,
      password,
      options: {
        data: {
          full_name: fullName,
          institution: institution || undefined,
          orcid: orcid || undefined,
        },
        emailRedirectTo: `${window.location.origin}/auth/callback?next=/runs`,
      },
    });

    if (signUpError) {
      setError("Could not create account. Please try again.");
      setLoading(false);
      return;
    }

    setSuccess(true);
    setLoading(false);
  }

  return (
    <div className="login-bg relative flex min-h-screen flex-col items-center justify-center px-4">
      <LoginBackground />

      {/* ── Card ── */}
      <div className="login-card relative z-10 w-full max-w-[420px] overflow-hidden px-8 py-10 sm:px-10">
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
              Computational Biology Platform
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
              data-testid="signup-error"
            >
              {error}
            </div>
          )}

          {/* Success */}
          {success ? (
            <div className="text-center">
              <div
                className="mb-6 rounded-xl px-4 py-8 text-sm"
                style={{
                  backgroundColor: "var(--app-badge-success-bg)",
                  color: "var(--app-badge-success-tx)",
                }}
                data-testid="signup-success"
              >
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.5} className="mx-auto mb-3 h-10 w-10">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75L11.25 15 15 9.75M21 12a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
                <p className="font-medium">Account created</p>
                <p className="mt-1 opacity-80">Check your email for a confirmation link.</p>
              </div>
              <Link
                href="/login"
                className="text-sm font-medium underline"
                style={{ color: "var(--app-primary)" }}
              >
                Back to login
              </Link>
            </div>
          ) : (
            <>
              <form onSubmit={handleSignup} className="space-y-4">
                <div>
                  <label htmlFor="fullName" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                    Full Name
                  </label>
                  <input
                    id="fullName" type="text" required value={fullName}
                    onChange={(e) => setFullName(e.target.value)}
                    placeholder="Jane Doe"
                    className="login-input app-input w-full rounded-lg px-4 py-2.5 text-sm outline-none"
                    data-testid="input-full-name"
                  />
                </div>

                <div>
                  <label htmlFor="signup-email" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                    Email
                  </label>
                  <input
                    id="signup-email" type="email" required value={email}
                    onChange={(e) => setEmail(e.target.value)}
                    placeholder="you@university.edu"
                    className="login-input app-input w-full rounded-lg px-4 py-2.5 text-sm outline-none"
                    data-testid="input-email"
                  />
                </div>

                <div>
                  <label htmlFor="signup-password" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                    Password
                  </label>
                  <input
                    id="signup-password" type="password" required minLength={8}
                    value={password}
                    onChange={(e) => setPassword(e.target.value)}
                    placeholder="At least 8 characters"
                    className="login-input app-input w-full rounded-lg px-4 py-2.5 text-sm outline-none"
                    data-testid="input-password"
                  />
                </div>

                {/* Institution + ORCID side by side */}
                <div className="grid grid-cols-2 gap-3">
                  <div>
                    <label htmlFor="institution" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                      Institution <span className="font-normal" style={{ color: "var(--app-text-3)" }}>(opt.)</span>
                    </label>
                    <input
                      id="institution" type="text" value={institution}
                      onChange={(e) => setInstitution(e.target.value)}
                      placeholder="e.g. UFMG"
                      className="login-input app-input w-full rounded-lg px-3 py-2.5 text-sm outline-none"
                      data-testid="input-institution"
                    />
                  </div>
                  <div>
                    <label htmlFor="orcid" className="mb-1.5 block text-[13px] font-medium" style={{ color: "var(--app-text-2)" }}>
                      ORCID <span className="font-normal" style={{ color: "var(--app-text-3)" }}>(opt.)</span>
                    </label>
                    <input
                      id="orcid" type="text" value={orcid}
                      onChange={(e) => setOrcid(e.target.value)}
                      placeholder="0000-0000-..."
                      pattern="\d{4}-\d{4}-\d{4}-\d{3}[\dX]"
                      title="ORCID format: 0000-0000-0000-0000"
                      className="login-input app-input w-full rounded-lg px-3 py-2.5 text-sm outline-none"
                      data-testid="input-orcid"
                    />
                  </div>
                </div>

                <button
                  type="submit" disabled={loading}
                  className="login-btn app-btn-primary w-full rounded-lg py-2.5 text-sm font-semibold disabled:opacity-50"
                  style={{ letterSpacing: "0.02em" }}
                  data-testid="btn-signup"
                >
                  {loading ? "Creating account..." : "Create Account"}
                </button>
              </form>

              {/* Divider + Login */}
              <div className="my-6 flex items-center gap-3">
                <div className="h-px flex-1" style={{ backgroundColor: "var(--app-border)" }} />
                <span className="text-[11px]" style={{ color: "var(--app-text-3)" }}>Already have an account?</span>
                <div className="h-px flex-1" style={{ backgroundColor: "var(--app-border)" }} />
              </div>

              <Link
                href="/login"
                className="login-btn block w-full rounded-lg border py-2.5 text-center text-sm font-medium"
                style={{
                  borderColor: "var(--app-border)",
                  color: "var(--app-text)",
                  backgroundColor: "transparent",
                }}
                data-testid="link-login"
              >
                Sign in
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
