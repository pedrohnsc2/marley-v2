interface KpiCardProps {
  title: string;
  value: string | number;
  subtitle?: string;
  accentColor?: string;
}

export default function KpiCard({
  title,
  value,
  subtitle,
  accentColor = "border-cyan-500",
}: KpiCardProps) {
  return (
    <div
      className={`rounded-lg border-l-4 ${accentColor} border border-zinc-800 bg-zinc-900/50 px-5 py-4`}
      data-testid="kpi-card"
    >
      <p className="text-xs font-medium uppercase tracking-wider text-zinc-500">
        {title}
      </p>
      <p className="mt-1 text-2xl font-bold text-white">{value}</p>
      {subtitle && (
        <p className="mt-1 text-xs text-zinc-500">{subtitle}</p>
      )}
    </div>
  );
}
