"use client";

import { useEffect, useRef } from "react";
import { createPortal } from "react-dom";
import { useTranslations } from "next-intl";
import Lottie from "lottie-react";
import dogNoseAnimation from "@/public/dog-nose.json";
import { useOnboarding } from "./OnboardingProvider";

const FEATURE_KEYS = ["pipeline", "drugTargets", "vaccine", "docking", "aso", "bioSim"] as const;

const FEATURE_ICONS: Record<string, React.ReactNode> = {
  pipeline: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
      <path d="M4 6h16M4 12h16M4 18h16" />
      <path d="M8 6l4 3-4 3" fill="currentColor" stroke="none" />
    </svg>
  ),
  drugTargets: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
      <circle cx="12" cy="12" r="9" /><path d="M12 3v18M3 12h18" />
    </svg>
  ),
  vaccine: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
      <path d="M9 3h6v4H9z" /><path d="M12 7v14" />
      <path d="M8 11h8" /><path d="M10 15h4" />
    </svg>
  ),
  docking: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
      <circle cx="12" cy="12" r="3" />
      <circle cx="5" cy="6" r="2" /><circle cx="19" cy="6" r="2" />
      <circle cx="5" cy="18" r="2" /><circle cx="19" cy="18" r="2" />
      <path d="M7 7l3 3M14 14l3 3M17 7l-3 3M10 14l-3 3" />
    </svg>
  ),
  aso: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
      <path d="M4 4l16 16M4 20L20 4" />
      <circle cx="4" cy="4" r="2" /><circle cx="20" cy="20" r="2" />
      <circle cx="4" cy="20" r="2" /><circle cx="20" cy="4" r="2" />
    </svg>
  ),
  bioSim: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
      <circle cx="12" cy="12" r="9" />
      <polygon points="10,8 16,12 10,16" fill="currentColor" />
    </svg>
  ),
};

export default function WelcomeModal() {
  const { showWelcome, completeWelcome } = useOnboarding();
  const t = useTranslations("onboarding");
  const startBtnRef = useRef<HTMLButtonElement>(null);

  useEffect(() => {
    if (showWelcome) {
      startBtnRef.current?.focus();
    }
  }, [showWelcome]);

  useEffect(() => {
    if (!showWelcome) return;

    function handleEsc(e: KeyboardEvent) {
      if (e.key === "Escape") completeWelcome(false);
    }
    document.addEventListener("keydown", handleEsc);
    return () => document.removeEventListener("keydown", handleEsc);
  }, [showWelcome, completeWelcome]);

  if (!showWelcome) return null;

  return createPortal(
    <div
      className="onboarding-backdrop"
      data-testid="onboarding-welcome"
      onClick={(e) => {
        if (e.target === e.currentTarget) completeWelcome(false);
      }}
    >
      <div
        className="onboarding-modal"
        role="dialog"
        aria-modal="true"
        aria-labelledby="onboarding-welcome-title"
      >
        <div className="flex items-center gap-3 mb-4">
          <div style={{ width: 48, height: 48, filter: "var(--app-logo-filter)" }}>
            <Lottie animationData={dogNoseAnimation} loop autoplay style={{ width: 48, height: 48 }} />
          </div>
          <div>
            <h2
              id="onboarding-welcome-title"
              className="text-xl font-bold"
              style={{ color: "var(--app-text)" }}
            >
              {t("welcome.title")}
            </h2>
          </div>
        </div>

        <p className="text-sm mb-6" style={{ color: "var(--app-text-2)" }}>
          {t("welcome.subtitle")}
        </p>

        <div className="onboarding-features">
          {FEATURE_KEYS.map((key) => (
            <div key={key} className="onboarding-feature-card">
              <div className="mb-2" style={{ color: "var(--app-primary)" }}>
                {FEATURE_ICONS[key]}
              </div>
              <p className="text-sm font-semibold mb-1" style={{ color: "var(--app-text)" }}>
                {t(`welcome.features.${key}.title`)}
              </p>
              <p className="text-xs leading-relaxed" style={{ color: "var(--app-text-2)" }}>
                {t(`welcome.features.${key}.desc`)}
              </p>
            </div>
          ))}
        </div>

        <div className="flex items-center gap-3 mt-6">
          <button
            ref={startBtnRef}
            onClick={() => completeWelcome(true)}
            className="app-btn-primary"
            data-testid="onboarding-start-tour"
          >
            {t("welcome.startTour")}
          </button>
          <button
            onClick={() => completeWelcome(false)}
            className="text-sm font-medium"
            style={{ color: "var(--app-text-3)", background: "none", border: "none", cursor: "pointer" }}
            data-testid="onboarding-skip-welcome"
          >
            {t("welcome.skip")}
          </button>
        </div>
      </div>
    </div>,
    document.body,
  );
}
