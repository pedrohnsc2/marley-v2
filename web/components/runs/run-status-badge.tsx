import type { RunStatus } from "@/lib/types/run";

interface RunStatusBadgeProps {
  status: string;
}

const STATUS_STYLES: Record<RunStatus, { bg: string; text: string }> = {
  completed: { bg: "rgb(16 185 129 / 0.15)", text: "rgb(16 185 129)" },
  failed:    { bg: "rgb(244 63 94 / 0.15)",   text: "rgb(244 63 94)" },
  running:   { bg: "rgb(245 158 11 / 0.15)",  text: "rgb(245 158 11)" },
  created:   { bg: "rgb(161 161 170 / 0.15)", text: "rgb(161 161 170)" },
  cancelled: { bg: "rgb(161 161 170 / 0.15)", text: "rgb(161 161 170)" },
};

export default function RunStatusBadge({ status }: RunStatusBadgeProps) {
  const style = STATUS_STYLES[status as RunStatus] ?? STATUS_STYLES.created;

  return (
    <span
      className="inline-flex items-center rounded-full px-2 py-0.5 text-xs font-medium"
      style={{ backgroundColor: style.bg, color: style.text }}
      data-testid="run-status-badge"
    >
      {status}
    </span>
  );
}
