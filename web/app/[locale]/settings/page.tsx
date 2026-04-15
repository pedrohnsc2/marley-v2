"use client";

import ThemeSelector from "@/components/theme-selector";
import { useOnboarding } from "@/components/onboarding/OnboardingProvider";
import { useTranslations } from "next-intl";

const FEATURE_KEYS = ["dashboard", "runs", "drug", "vaccine", "docking", "aso", "bioSim", "settings"] as const;

export default function SettingsPage() {
  const { startTour, resetOnboarding, hintsEnabled, setHintsEnabled, state } = useOnboarding();
  const t = useTranslations("onboarding");

  const exploredCount = FEATURE_KEYS.filter((k) => state.featuresVisited[k]).length;
  const totalFeatures = FEATURE_KEYS.length;

  return (
    <div>
      <div className="mb-6">
        <h1 className="text-2xl font-bold" style={{ color: "var(--app-text)" }}>
          Settings
        </h1>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-2)" }}>
          Manage your dashboard preferences
        </p>
      </div>

      {/* Theme section */}
      <div
        className="rounded-xl p-6 mb-6"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <div className="mb-6">
          <h2 className="text-base font-semibold" style={{ color: "var(--app-text)" }}>
            Appearance
          </h2>
          <p className="mt-1 text-sm" style={{ color: "var(--app-text-2)" }}>
            Choose a visual theme for the dashboard
          </p>
        </div>
        <ThemeSelector />
      </div>

      {/* Onboarding section */}
      <div
        className="rounded-xl p-6 mb-6"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
        data-testid="onboarding-settings"
      >
        <div className="mb-4">
          <h2 className="text-base font-semibold" style={{ color: "var(--app-text)" }}>
            {t("settings.title")}
          </h2>
          <p className="mt-1 text-sm" style={{ color: "var(--app-text-2)" }}>
            {t("settings.subtitle")}
          </p>
        </div>

        {/* Progress bar */}
        <div className="mb-4">
          <div className="flex justify-between text-sm mb-2">
            <span style={{ color: "var(--app-text-2)" }}>
              {exploredCount >= totalFeatures
                ? t("settings.progressComplete")
                : t("settings.progress", { count: exploredCount, total: totalFeatures })}
            </span>
          </div>
          <div className="h-2 rounded-full" style={{ backgroundColor: "var(--app-surface-2)" }}>
            <div
              className="h-2 rounded-full transition-all duration-500"
              style={{
                backgroundColor: "var(--app-primary)",
                width: `${(exploredCount / totalFeatures) * 100}%`,
              }}
              data-testid="onboarding-settings-progress"
            />
          </div>
        </div>

        {/* Feature checklist */}
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2 mb-4">
          {FEATURE_KEYS.map((key) => (
            <div
              key={key}
              className="flex items-center gap-2 text-sm"
              style={{ color: state.featuresVisited[key] ? "var(--app-text)" : "var(--app-text-3)" }}
            >
              {state.featuresVisited[key] ? (
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} className="h-4 w-4 flex-shrink-0" style={{ color: "var(--app-primary)" }}>
                  <path d="M20 6L9 17l-5-5" />
                </svg>
              ) : (
                <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.5} className="h-4 w-4 flex-shrink-0">
                  <circle cx="12" cy="12" r="9" />
                </svg>
              )}
              <span>{t(`settings.features.${key}`)}</span>
            </div>
          ))}
        </div>

        {/* Actions */}
        <div className="flex flex-wrap items-center gap-3">
          <button
            onClick={() => {
              resetOnboarding();
              setTimeout(() => startTour(), 100);
            }}
            className="app-btn-primary"
            data-testid="onboarding-replay-tour"
          >
            {t("settings.replayTour")}
          </button>
          <button
            onClick={() => setHintsEnabled(!hintsEnabled)}
            className="rounded-lg px-4 py-2 text-sm font-medium transition-colors"
            style={{
              backgroundColor: "var(--app-surface-2)",
              color: "var(--app-text-2)",
              border: "1px solid var(--app-border)",
            }}
            data-testid="onboarding-hints-toggle"
          >
            {t("settings.hintsToggle")}: {hintsEnabled ? "ON" : "OFF"}
          </button>
        </div>
      </div>

      {/* About section */}
      <div
        className="rounded-xl p-6"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <h2 className="mb-4 text-base font-semibold" style={{ color: "var(--app-text)" }}>
          About
        </h2>
        <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
          {[
            { label: "Project", value: "Marley v2.0" },
            { label: "Pipeline", value: "Reverse Vaccinology" },
            { label: "Target", value: "L. infantum" },
            { label: "Species", value: "Canine (CVL)" },
          ].map((item) => (
            <div
              key={item.label}
              className="rounded-lg p-3"
              style={{ backgroundColor: "var(--app-surface-2)", border: "1px solid var(--app-border)" }}
            >
              <p className="text-xs" style={{ color: "var(--app-text-3)" }}>{item.label}</p>
              <p className="mt-1 text-sm font-semibold" style={{ color: "var(--app-text)" }}>{item.value}</p>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
