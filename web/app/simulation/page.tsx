import { loadJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";

interface SimPeak {
  peak_level: number;
  day_of_peak: number;
  final_level: number;
}

interface SimData {
  generated: string;
  vaccine_params: {
    epitope_count: number;
    avg_ic50: number;
    ctl_count: number;
    htl_count: number;
    adjuvant: string;
  };
  dose_schedule: { day: number; amount: number }[];
  simulation_days: number;
  peaks: Record<string, SimPeak>;
  th1_dominance_pct: number;
  th1_peak_day: number;
  memory_after_doses: number[];
  memory_stability_pct: number;
  estimated_protection: string;
}

export default function SimulationPage() {
  const sim = loadJson("immune_sim/immune_sim_data.json") as SimData;

  const cellTypes = [
    { key: "Th", label: "T-helper (Th)", color: "text-blue-400" },
    { key: "Tc", label: "Cytotoxic T (Tc)", color: "text-red-400" },
    { key: "B", label: "B cells", color: "text-yellow-400" },
    { key: "Ab", label: "Antibodies", color: "text-green-400" },
    { key: "M", label: "Memory cells", color: "text-purple-400" },
  ];

  return (
    <div className="px-8 py-10">
      <header className="mb-8">
        <div className="flex items-center gap-2">
          <span className="rounded bg-purple-900/50 px-2 py-0.5 text-xs font-mono text-purple-400">
            v5
          </span>
          <h1 className="text-2xl font-bold tracking-tight text-white">
            Immune Simulation
          </h1>
        </div>
        <p className="mt-2 text-sm text-zinc-400">
          ODE-based immune kinetics model simulating the vaccine response over{" "}
          {sim.simulation_days} days with {sim.dose_schedule.length}-dose
          schedule
        </p>
      </header>

      <section className="mb-10">
        <div className="grid grid-cols-2 gap-4 lg:grid-cols-4">
          <KpiCard
            title="Th1 Dominance"
            value={`${sim.th1_dominance_pct}%`}
            subtitle="Pro-inflammatory bias"
            accentColor="border-purple-500"
          />
          <KpiCard
            title="Memory Stability"
            value={`${sim.memory_stability_pct}%`}
            subtitle="Post-vaccination"
            accentColor="border-purple-400"
          />
          <KpiCard
            title="Est. Protection"
            value={sim.estimated_protection.replace("~", "").replace(">", ">")}
            subtitle="Based on memory decay model"
            accentColor="border-purple-500"
          />
          <KpiCard
            title="Dose Schedule"
            value={`${sim.dose_schedule.length} doses`}
            subtitle={sim.dose_schedule.map((d) => `Day ${d.day}`).join(", ")}
            accentColor="border-purple-400"
          />
        </div>
      </section>

      <section className="mb-10">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          Peak Immune Response
        </h2>
        <div className="rounded-lg border border-zinc-800 bg-zinc-900/40 p-5">
          <div className="grid gap-4 md:grid-cols-5">
            {cellTypes.map(({ key, label, color }) => {
              const peak = sim.peaks[key];
              if (!peak) return null;
              return (
                <div key={key} className="text-center">
                  <p className="text-xs text-zinc-500">{label}</p>
                  <p className={`mt-1 text-xl font-bold ${color}`}>
                    {peak.peak_level.toFixed(2)}
                  </p>
                  <p className="text-xs text-zinc-600">
                    Day {peak.day_of_peak.toFixed(0)}
                  </p>
                </div>
              );
            })}
          </div>
        </div>
      </section>

      <section className="mb-10">
        <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
          Memory Cell Formation
        </h2>
        <div className="rounded-lg border border-zinc-800 bg-zinc-900/40 p-5">
          <div className="flex items-end gap-6">
            {sim.memory_after_doses.map((val, i) => {
              const maxVal = Math.max(...sim.memory_after_doses);
              const heightPct = maxVal > 0 ? (val / maxVal) * 100 : 0;
              return (
                <div key={i} className="flex flex-col items-center gap-2">
                  <span className="text-xs font-mono text-purple-300">
                    {val.toFixed(4)}
                  </span>
                  <div
                    className="w-16 rounded-t bg-gradient-to-t from-purple-800 to-purple-500"
                    style={{ height: `${Math.max(heightPct, 5)}px`, minHeight: "8px" }}
                  />
                  <span className="text-xs text-zinc-500">Dose {i + 1}</span>
                </div>
              );
            })}
          </div>
          <p className="mt-4 text-xs text-zinc-500">
            Memory cell levels after each vaccine dose, showing progressive
            accumulation with each boost.
          </p>
        </div>
      </section>

      <div className="grid gap-8 xl:grid-cols-2">
        <section>
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            Immune Kinetics
          </h2>
          <div className="overflow-hidden rounded-lg border border-zinc-800">
            {/* eslint-disable-next-line @next/next/no-img-element */}
            <img
              src="/images/immune_kinetics.png"
              alt="Immune kinetics plot showing cell populations over 365 days"
              className="w-full"
            />
          </div>
          <p className="mt-2 text-xs text-zinc-500">
            Time course of immune cell populations (Th, Tc, B, Ab, Memory) over
            the {sim.simulation_days}-day simulation with boosts at days{" "}
            {sim.dose_schedule.map((d) => d.day).join(", ")}.
          </p>
        </section>

        <section>
          <h2 className="mb-4 text-sm font-semibold uppercase tracking-wider text-zinc-500">
            Th1/Th2 Balance
          </h2>
          <div className="overflow-hidden rounded-lg border border-zinc-800">
            {/* eslint-disable-next-line @next/next/no-img-element */}
            <img
              src="/images/th1_th2.png"
              alt="Th1/Th2 balance plot showing dominant Th1 response"
              className="w-full"
            />
          </div>
          <p className="mt-2 text-xs text-zinc-500">
            Th1/Th2 polarization showing {sim.th1_dominance_pct}% Th1 dominance,
            consistent with the pro-inflammatory response needed for
            intracellular pathogen clearance in leishmaniasis.
          </p>
        </section>
      </div>
    </div>
  );
}
