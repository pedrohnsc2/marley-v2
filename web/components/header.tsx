"use client";

import { usePathname } from "next/navigation";
import { useTheme } from "@/contexts/theme-context";

const breadcrumbMap: Record<string, string> = {
  "/": "Dashboard",
  "/vaccine": "Vaccine Construct",
  "/drug": "Drug Targets",
  "/docking": "Molecular Docking",
  "/simulation": "Immune Simulation",
  "/settings": "Settings",
};

export default function Header() {
  const pathname = usePathname();
  const { isTheme2, isDark, setTheme, theme } = useTheme();
  const pageTitle = breadcrumbMap[pathname] ?? "Page";

  const toggleDarkLight = () => {
    setTheme(isDark ? "theme2-light" : "theme2-dark");
  };

  return (
    <header
      className="sticky top-0 z-20 flex h-16 items-center gap-4 px-6"
      style={{
        backgroundColor: "var(--app-header-bg)",
        borderBottom: "1px solid var(--app-header-border)",
        boxShadow: "0 1px 3px rgb(0 0 0 / 0.04)",
      }}
    >
      {/* Breadcrumb */}
      <div className="flex items-center gap-2 text-sm">
        <span style={{ color: "var(--app-text-3)" }}>Marley</span>
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} className="h-3.5 w-3.5" style={{ color: "var(--app-text-3)" }}>
          <path d="M9 18l6-6-6-6" />
        </svg>
        <span className="font-semibold" style={{ color: "var(--app-text)" }}>{pageTitle}</span>
      </div>

      {/* Search */}
      <div className="ml-4 hidden flex-1 max-w-xs sm:block">
        <div className="relative">
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2" style={{ color: "var(--app-text-3)" }}>
            <circle cx="11" cy="11" r="7" />
            <path d="m21 21-4.35-4.35" />
          </svg>
          <input
            type="text"
            placeholder="Search..."
            className="app-input w-full rounded-lg py-2 pl-9 pr-4 text-sm outline-none focus:ring-2"
            style={{ "--tw-ring-color": "var(--app-primary)" } as React.CSSProperties}
          />
        </div>
      </div>

      <div className="ml-auto flex items-center gap-2">
        {/* Dark/Light toggle — only for Theme 2 */}
        {isTheme2 && (
          <button
            onClick={toggleDarkLight}
            className="flex h-9 w-9 items-center justify-center rounded-lg transition-colors"
            style={{
              backgroundColor: "var(--app-surface-2)",
              border: "1px solid var(--app-border)",
              color: "var(--app-text-2)",
            }}
            title={isDark ? "Switch to light mode" : "Switch to dark mode"}
          >
            {isDark ? (
              <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-4 w-4">
                <circle cx="12" cy="12" r="5" />
                <path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42" />
              </svg>
            ) : (
              <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-4 w-4">
                <path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z" />
              </svg>
            )}
          </button>
        )}

        {/* Notifications */}
        <button
          className="relative flex h-9 w-9 items-center justify-center rounded-lg transition-colors"
          style={{ color: "var(--app-text-2)" }}
        >
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={1.8} className="h-5 w-5">
            <path d="M18 8A6 6 0 0 0 6 8c0 7-3 9-3 9h18s-3-2-3-9" />
            <path d="M13.73 21a2 2 0 0 1-3.46 0" />
          </svg>
          <span
            className="absolute right-1.5 top-1.5 h-2 w-2 rounded-full"
            style={{ backgroundColor: "#EF4444" }}
          />
        </button>

        {/* Avatar */}
        <div className="flex items-center gap-2.5">
          <div
            className="flex h-9 w-9 items-center justify-center rounded-full text-sm font-semibold shadow-sm"
            style={{
              backgroundColor: "var(--app-primary)",
              color: "var(--app-primary-tx)",
            }}
          >
            R
          </div>
          <div className="hidden sm:block">
            <p className="text-sm font-semibold leading-tight" style={{ color: "var(--app-text)" }}>
              Researcher
            </p>
            <p className="text-xs" style={{ color: "var(--app-text-3)" }}>Admin</p>
          </div>
        </div>
      </div>
    </header>
  );
}
