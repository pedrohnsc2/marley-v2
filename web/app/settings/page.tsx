import ThemeSelector from "@/components/theme-selector";

export default function SettingsPage() {
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
