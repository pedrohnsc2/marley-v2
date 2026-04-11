"use client";

import { useTheme, type AppTheme } from "@/contexts/theme-context";

const themes: {
  id: AppTheme;
  label: string;
  description: string;
  preview: { bg: string; surface: string; sidebar: string; accent: string; text: string; border: string };
  group?: string;
}[] = [
  {
    id: "theme1",
    label: "Tailhub — Light",
    description: "Clean light interface with blue accents",
    group: "Theme 1",
    preview: {
      bg: "#F3F4F6",
      surface: "#FFFFFF",
      sidebar: "#FFFFFF",
      accent: "#3B82F6",
      text: "#111827",
      border: "#E5E7EB",
    },
  },
  {
    id: "theme2-dark",
    label: "LimeDock — Dark",
    description: "Dark navy with lime green accents",
    group: "Theme 2",
    preview: {
      bg: "#080D16",
      surface: "#0D1526",
      sidebar: "#080D16",
      accent: "#B4F740",
      text: "#F3F4F6",
      border: "#1C2840",
    },
  },
  {
    id: "theme2-light",
    label: "LimeDock — Light",
    description: "Light mode with lime green accents",
    group: "Theme 2",
    preview: {
      bg: "#F3F4F6",
      surface: "#FFFFFF",
      sidebar: "#FFFFFF",
      accent: "#B4F740",
      text: "#080D16",
      border: "#E5E7EB",
    },
  },
];

function ThemePreview({ preview, active }: {
  preview: typeof themes[0]["preview"];
  active: boolean;
}) {
  return (
    <div
      className="relative overflow-hidden rounded-lg"
      style={{
        backgroundColor: preview.bg,
        border: `2px solid ${active ? preview.accent : preview.border}`,
        height: 120,
      }}
    >
      {/* Mini sidebar */}
      <div
        className="absolute left-0 top-0 bottom-0 flex flex-col gap-1 p-2"
        style={{ width: 36, backgroundColor: preview.sidebar, borderRight: `1px solid ${preview.border}` }}
      >
        {/* Logo dot */}
        <div className="h-5 w-5 rounded-md mb-1" style={{ backgroundColor: preview.accent }} />
        {/* Nav items */}
        {[0, 1, 2].map((i) => (
          <div
            key={i}
            className="h-3 rounded-sm"
            style={{
              backgroundColor: i === 0 ? preview.accent : preview.border,
              opacity: i === 0 ? 1 : 0.6,
            }}
          />
        ))}
      </div>

      {/* Main area */}
      <div className="absolute left-9 top-0 right-0 bottom-0 p-2">
        {/* Header bar */}
        <div
          className="mb-2 h-5 rounded"
          style={{ backgroundColor: preview.surface, border: `1px solid ${preview.border}` }}
        />
        {/* Cards grid */}
        <div className="grid grid-cols-2 gap-1 mb-2">
          {[0, 1].map((i) => (
            <div
              key={i}
              className="h-8 rounded"
              style={{ backgroundColor: preview.surface, border: `1px solid ${preview.border}` }}
            >
              <div className="mx-1 mt-1 h-1.5 w-6 rounded" style={{ backgroundColor: preview.border }} />
              <div className="mx-1 mt-1 h-2 w-8 rounded" style={{ backgroundColor: preview.text, opacity: 0.6 }} />
            </div>
          ))}
        </div>
        {/* Chart placeholder */}
        <div
          className="h-12 rounded"
          style={{ backgroundColor: preview.surface, border: `1px solid ${preview.border}` }}
        >
          <div className="flex items-end gap-1 h-full p-2 pb-1">
            {[40, 70, 55, 90, 65, 80].map((h, i) => (
              <div
                key={i}
                className="flex-1 rounded-sm"
                style={{ height: `${h}%`, backgroundColor: preview.accent, opacity: i === 3 ? 1 : 0.5 }}
              />
            ))}
          </div>
        </div>
      </div>

      {/* Active checkmark */}
      {active && (
        <div
          className="absolute right-2 top-2 flex h-5 w-5 items-center justify-center rounded-full"
          style={{ backgroundColor: preview.accent }}
        >
          <svg viewBox="0 0 24 24" fill="none" stroke={preview.bg || "#080D16"} strokeWidth={3} className="h-3 w-3">
            <path d="M20 6L9 17l-5-5" />
          </svg>
        </div>
      )}
    </div>
  );
}

export default function ThemeSelector() {
  const { theme, setTheme } = useTheme();

  // Group themes
  const groups = [
    { label: "Theme 1 — Tailhub", items: themes.filter((t) => t.group === "Theme 1") },
    { label: "Theme 2 — LimeDock", items: themes.filter((t) => t.group === "Theme 2") },
  ];

  return (
    <div className="space-y-8">
      {groups.map((group) => (
        <div key={group.label}>
          <p className="mb-3 text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-3)" }}>
            {group.label}
          </p>
          <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
            {group.items.map((t) => {
              const active = theme === t.id;
              return (
                <button
                  key={t.id}
                  onClick={() => setTheme(t.id)}
                  className="text-left rounded-xl p-1 transition-all"
                  style={{
                    outline: active ? `2px solid ${t.preview.accent}` : "2px solid transparent",
                    outlineOffset: 2,
                  }}
                >
                  <ThemePreview preview={t.preview} active={active} />
                  <div className="mt-2 px-1">
                    <div className="flex items-center gap-2">
                      <p className="text-sm font-semibold" style={{ color: "var(--app-text)" }}>
                        {t.label}
                      </p>
                      {active && (
                        <span
                          className="rounded-full px-2 py-0.5 text-xs font-semibold"
                          style={{
                            backgroundColor: "var(--app-primary-lt)",
                            color: "var(--app-primary)",
                          }}
                        >
                          Active
                        </span>
                      )}
                    </div>
                    <p className="mt-0.5 text-xs" style={{ color: "var(--app-text-3)" }}>
                      {t.description}
                    </p>
                  </div>
                </button>
              );
            })}
          </div>
        </div>
      ))}

      {/* Theme 2 dark/light quick toggle */}
      {(theme === "theme2-dark" || theme === "theme2-light") && (
        <div
          className="flex items-center justify-between rounded-xl p-4"
          style={{ backgroundColor: "var(--app-surface-2)", border: "1px solid var(--app-border)" }}
        >
          <div>
            <p className="text-sm font-semibold" style={{ color: "var(--app-text)" }}>
              Dark / Light Mode
            </p>
            <p className="text-xs mt-0.5" style={{ color: "var(--app-text-3)" }}>
              Toggle between dark and light for LimeDock
            </p>
          </div>
          <button
            onClick={() => setTheme(theme === "theme2-dark" ? "theme2-light" : "theme2-dark")}
            className="relative inline-flex h-6 w-11 flex-shrink-0 cursor-pointer rounded-full transition-colors duration-200"
            style={{ backgroundColor: theme === "theme2-dark" ? "#1C2840" : "var(--app-primary)" }}
            role="switch"
            aria-checked={theme === "theme2-dark"}
          >
            <span
              className="inline-block h-5 w-5 transform rounded-full shadow transition duration-200"
              style={{
                backgroundColor: "white",
                transform: theme === "theme2-dark" ? "translateX(1px) translateY(2px)" : "translateX(22px) translateY(2px)",
              }}
            />
          </button>
        </div>
      )}
    </div>
  );
}
