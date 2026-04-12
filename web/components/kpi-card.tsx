interface KpiCardProps {
  title: string;
  value: string | number;
  subtitle?: string;
  accentColor?: string;
  icon?: React.ReactNode;
}

export default function KpiCard({ title, value, subtitle, accentColor = "bg-blue-500", icon }: KpiCardProps) {
  return (
    <div
      className="rounded-xl p-5"
      style={{
        backgroundColor: "var(--app-surface)",
        border: "1px solid var(--app-border)",
        boxShadow: "var(--app-card-shadow)",
      }}
      data-testid="kpi-card"
    >
      <div className="flex items-start justify-between">
        <div className="flex-1 min-w-0">
          <p className="text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            {title}
          </p>
          <p className="mt-2 text-2xl font-bold leading-tight truncate" style={{ color: "var(--app-text)" }}>
            {value}
          </p>
          {subtitle && (
            <p className="mt-1 text-xs truncate" style={{ color: "var(--app-text-3)" }}>{subtitle}</p>
          )}
        </div>
        {icon && (
          <div
            className="ml-4 flex-shrink-0 flex h-11 w-11 items-center justify-center rounded-xl"
            style={{ backgroundColor: "var(--app-surface-2)", color: "var(--app-text-3)" }}
          >
            {icon}
          </div>
        )}
      </div>
      <div
        className={`mt-4 h-1 w-full rounded-full ${accentColor} opacity-70`}
        style={{ backgroundColor: "var(--app-accent-bar)", opacity: 0.7 }}
      />
    </div>
  );
}
