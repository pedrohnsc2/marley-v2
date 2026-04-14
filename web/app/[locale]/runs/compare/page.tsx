import { getTranslations } from "next-intl/server";
import { listRuns, loadRunMetadata } from "@/lib/data-loader";
import type { RunMetadata } from "@/lib/types/run";

function formatDuration(seconds: number): string {
  if (seconds < 1) return "--";
  if (seconds < 60) return `${seconds.toFixed(1)}s`;
  return `${(seconds / 60).toFixed(1)}m`;
}

function formatDate(iso: string): string {
  if (!iso) return "--";
  return iso.slice(0, 16).replace("T", " ");
}

export default async function ComparePage({
  searchParams,
}: {
  searchParams: Promise<{ a?: string; b?: string }>;
}) {
  const t = await getTranslations("runs");
  const { a: runIdA, b: runIdB } = await searchParams;

  // If no runs selected, show selector
  if (!runIdA || !runIdB) {
    const allRuns = listRuns() as RunMetadata[];
    return (
      <div className="space-y-6">
        <div>
          <h1 className="text-2xl font-bold" style={{ color: "var(--app-text)" }}>
            {t("compare.title")}
          </h1>
          <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
            {t("compare.selectPrompt")}
          </p>
        </div>

        <form className="grid grid-cols-1 gap-4 sm:grid-cols-2">
          <div>
            <label className="mb-1.5 block text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
              Run A
            </label>
            <select
              name="a"
              defaultValue={runIdA ?? ""}
              data-testid="compare-run-a"
              className="w-full rounded-lg border px-3 py-2 text-sm"
              style={{
                backgroundColor: "var(--app-surface)",
                borderColor: "var(--app-border)",
                color: "var(--app-text)",
              }}
            >
              <option value="">-- Select --</option>
              {allRuns.map((r) => (
                <option key={r.run_id} value={r.run_id}>
                  {r.run_id} ({r.pipeline})
                </option>
              ))}
            </select>
          </div>

          <div>
            <label className="mb-1.5 block text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
              Run B
            </label>
            <select
              name="b"
              defaultValue={runIdB ?? ""}
              data-testid="compare-run-b"
              className="w-full rounded-lg border px-3 py-2 text-sm"
              style={{
                backgroundColor: "var(--app-surface)",
                borderColor: "var(--app-border)",
                color: "var(--app-text)",
              }}
            >
              <option value="">-- Select --</option>
              {allRuns.map((r) => (
                <option key={r.run_id} value={r.run_id}>
                  {r.run_id} ({r.pipeline})
                </option>
              ))}
            </select>
          </div>

          <div className="sm:col-span-2">
            <button
              type="submit"
              data-testid="compare-submit"
              className="rounded-lg px-4 py-2 text-sm font-medium"
              style={{
                backgroundColor: "var(--app-accent-bar)",
                color: "white",
              }}
            >
              {t("compare.compareButton")}
            </button>
          </div>
        </form>
      </div>
    );
  }

  // Load both runs
  let runA: RunMetadata;
  let runB: RunMetadata;
  try {
    runA = loadRunMetadata(runIdA) as RunMetadata;
    runB = loadRunMetadata(runIdB) as RunMetadata;
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Unknown error";
    return (
      <div className="p-8 text-center" style={{ color: "var(--app-text-3)" }}>
        {t("compare.notFound")}: {message}
      </div>
    );
  }

  // Parameter diff
  const allParamKeys = new Set([
    ...Object.keys(runA.parameters ?? {}),
    ...Object.keys(runB.parameters ?? {}),
  ]);
  const paramDiffs: { key: string; a: unknown; b: unknown; changed: boolean }[] = [];
  for (const key of [...allParamKeys].sort()) {
    const a = (runA.parameters ?? {})[key];
    const b = (runB.parameters ?? {})[key];
    paramDiffs.push({ key, a, b, changed: JSON.stringify(a) !== JSON.stringify(b) });
  }

  // Stage comparison
  const stagesA = new Map((runA.stages ?? []).map((s) => [s.stage_id, s]));
  const stagesB = new Map((runB.stages ?? []).map((s) => [s.stage_id, s]));
  const allStageIds = new Set([...stagesA.keys(), ...stagesB.keys()]);

  return (
    <div className="space-y-8">
      {/* Header */}
      <div>
        <h1 className="text-2xl font-bold" style={{ color: "var(--app-text)" }}>
          {t("compare.title")}
        </h1>
        <p className="mt-1 text-sm" style={{ color: "var(--app-text-3)" }}>
          {runIdA} vs {runIdB}
        </p>
      </div>

      {/* Summary cards */}
      <div className="grid grid-cols-2 gap-4">
        {[
          { label: "Run A", run: runA },
          { label: "Run B", run: runB },
        ].map(({ label, run }) => (
          <div
            key={label}
            className="rounded-xl p-4"
            style={{
              backgroundColor: "var(--app-surface)",
              border: "1px solid var(--app-border)",
              boxShadow: "var(--app-card-shadow)",
            }}
          >
            <p className="text-xs font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
              {label}
            </p>
            <p className="mt-1 text-sm font-mono truncate" style={{ color: "var(--app-text)" }}>
              {run.run_id}
            </p>
            <div className="mt-2 flex gap-4 text-xs" style={{ color: "var(--app-text-3)" }}>
              <span>{run.pipeline}</span>
              <span>{run.status}</span>
              <span>{formatDuration(run.total_duration_s ?? 0)}</span>
              <span>{formatDate(run.created_at ?? "")}</span>
            </div>
          </div>
        ))}
      </div>

      {/* Parameter diff */}
      <div
        className="rounded-xl overflow-hidden"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <div className="px-5 py-3" style={{ borderBottom: "1px solid var(--app-border)" }}>
          <h2 className="text-sm font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            {t("compare.parameters")}
          </h2>
        </div>
        <table className="w-full text-sm">
          <thead>
            <tr style={{ borderBottom: "1px solid var(--app-border)" }}>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>Parameter</th>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>Run A</th>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>Run B</th>
            </tr>
          </thead>
          <tbody>
            {paramDiffs.map((d) => (
              <tr
                key={d.key}
                style={{
                  borderBottom: "1px solid var(--app-border)",
                  backgroundColor: d.changed ? "rgba(251,191,36,0.06)" : "transparent",
                }}
              >
                <td className="px-5 py-2 font-mono text-xs" style={{ color: "var(--app-text)" }}>
                  {d.key}
                </td>
                <td className="px-5 py-2 font-mono text-xs" style={{ color: d.changed ? "rgb(251,191,36)" : "var(--app-text-3)" }}>
                  {JSON.stringify(d.a) ?? "--"}
                </td>
                <td className="px-5 py-2 font-mono text-xs" style={{ color: d.changed ? "rgb(251,191,36)" : "var(--app-text-3)" }}>
                  {JSON.stringify(d.b) ?? "--"}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Stage comparison */}
      <div
        className="rounded-xl overflow-hidden"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <div className="px-5 py-3" style={{ borderBottom: "1px solid var(--app-border)" }}>
          <h2 className="text-sm font-semibold uppercase tracking-wider" style={{ color: "var(--app-text-2)" }}>
            {t("compare.stages")}
          </h2>
        </div>
        <table className="w-full text-sm">
          <thead>
            <tr style={{ borderBottom: "1px solid var(--app-border)" }}>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>Stage</th>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>A Status</th>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>A Duration</th>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>B Status</th>
              <th className="px-5 py-2 text-left text-xs font-semibold uppercase" style={{ color: "var(--app-text-3)" }}>B Duration</th>
            </tr>
          </thead>
          <tbody>
            {[...allStageIds].map((id) => {
              const sA = stagesA.get(id);
              const sB = stagesB.get(id);
              return (
                <tr key={id} style={{ borderBottom: "1px solid var(--app-border)" }}>
                  <td className="px-5 py-2 text-xs" style={{ color: "var(--app-text)" }}>
                    {sA?.name ?? sB?.name ?? id}
                  </td>
                  <td className="px-5 py-2 text-xs" style={{ color: "var(--app-text-3)" }}>
                    {sA?.status ?? "--"}
                  </td>
                  <td className="px-5 py-2 text-xs font-mono" style={{ color: "var(--app-text-3)" }}>
                    {sA ? formatDuration(sA.duration_s) : "--"}
                  </td>
                  <td className="px-5 py-2 text-xs" style={{ color: "var(--app-text-3)" }}>
                    {sB?.status ?? "--"}
                  </td>
                  <td className="px-5 py-2 text-xs font-mono" style={{ color: "var(--app-text-3)" }}>
                    {sB ? formatDuration(sB.duration_s) : "--"}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}
