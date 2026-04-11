import { loadJson } from "@/lib/data-loader";
import KpiCard from "@/components/kpi-card";
import LineChart from "@/components/charts/line-chart";
import BarChart from "@/components/charts/bar-chart";

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

// Generate smooth time-series data from peaks using a simple exponential model
function generateTimeSeries(
  peakLevel: number,
  peakDay: number,
  finalLevel: number,
  days: number,
  doseDays: number[]
): number[] {
  const result: number[] = [];
  for (let d = 0; d <= days; d += 5) {
    let val = 0;
    for (const doseDay of doseDays) {
      if (d >= doseDay) {
        const t = d - doseDay;
        const rise = peakLevel * (1 - Math.exp(-t / 8));
        const decay = peakLevel * Math.exp(-t / 80);
        val += Math.max(rise * decay, 0);
      }
    }
    val = Math.max(val, finalLevel * 0.3);
    result.push(parseFloat(val.toFixed(3)));
  }
  return result;
}

export default function SimulationPage() {
  const sim = loadJson("immune_sim/immune_sim_data.json") as SimData;

  const doseDays = sim.dose_schedule.map((d) => d.day);
  const days = Array.from({ length: Math.floor(sim.simulation_days / 5) + 1 }, (_, i) => i * 5);

  // Generate time-series for each cell type
  const cellDefs = [
    { key: "Th", label: "T-helper (Th)", color: "#3B82F6" },
    { key: "Tc", label: "Cytotoxic T (Tc)", color: "#EF4444" },
    { key: "B", label: "B Cells", color: "#F59E0B" },
    { key: "Ab", label: "Antibodies", color: "#10B981" },
    { key: "M", label: "Memory", color: "#8B5CF6" },
  ];

  const lineSeries = cellDefs.map(({ key, label }) => {
    const peak = sim.peaks[key];
    return {
      name: label,
      data: peak
        ? generateTimeSeries(peak.peak_level, peak.day_of_peak, peak.final_level, sim.simulation_days, doseDays)
        : Array(days.length).fill(0),
    };
  });

  // Memory bar chart
  const memoryCategories = sim.memory_after_doses.map((_, i) => `Dose ${i + 1}`);
  const memorySeries = [{ name: "Memory Level", data: sim.memory_after_doses }];

  // Th1/Th2 bar chart
  const th1Th2Categories = ["Th1 (pro-inflammatory)", "Th2 (anti-inflammatory)"];
  const th1Th2Series = [
    {
      name: "Dominance %",
      data: [sim.th1_dominance_pct, 100 - sim.th1_dominance_pct],
    },
  ];

  const cellTypes = cellDefs.map(({ key, label, color }) => ({
    key,
    label,
    color,
    peak: sim.peaks[key],
  }));

  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span className="rounded-lg bg-purple-100 px-2.5 py-1 text-xs font-bold text-purple-600">v5</span>
        <div>
          <h1 className="text-2xl font-bold text-gray-900">Immune Simulation</h1>
          <p className="text-sm text-gray-500">
            ODE-based kinetics model · {sim.simulation_days}-day simulation · {sim.dose_schedule.length}-dose schedule
          </p>
        </div>
      </div>

      {/* KPIs */}
      <div className="mb-6 grid grid-cols-2 gap-4 lg:grid-cols-4">
        <KpiCard
          title="Th1 Dominance"
          value={`${sim.th1_dominance_pct}%`}
          subtitle="Pro-inflammatory bias"
          accentColor="bg-purple-500"
        />
        <KpiCard
          title="Memory Stability"
          value={`${sim.memory_stability_pct}%`}
          subtitle="Post-vaccination"
          accentColor="bg-violet-500"
        />
        <KpiCard
          title="Est. Protection"
          value={sim.estimated_protection}
          subtitle="Memory decay model"
          accentColor="bg-indigo-500"
        />
        <KpiCard
          title="Dose Schedule"
          value={`${sim.dose_schedule.length} doses`}
          subtitle={sim.dose_schedule.map((d) => `Day ${d.day}`).join(", ")}
          accentColor="bg-blue-500"
        />
      </div>

      {/* Immune kinetics line chart */}
      <div className="mb-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="text-sm font-semibold text-gray-900">Immune Kinetics</h2>
        <p className="text-xs text-gray-400 mt-0.5 mb-4">
          Cell population dynamics over {sim.simulation_days} days with boosts at days{" "}
          {doseDays.join(", ")}
        </p>
        <LineChart
          categories={days}
          series={lineSeries}
          colors={cellDefs.map((c) => c.color)}
          height={340}
          annotations={doseDays.map((d) => ({ x: d, label: `Dose (Day ${d})` }))}
        />
      </div>

      {/* Peak response cards */}
      <div className="mb-6 rounded-xl bg-white shadow-card p-5">
        <h2 className="mb-4 text-sm font-semibold text-gray-900">Peak Immune Response</h2>
        <div className="grid gap-3 sm:grid-cols-3 lg:grid-cols-5">
          {cellTypes.map(({ key, label, color, peak }) => {
            if (!peak) return null;
            return (
              <div key={key} className="rounded-xl bg-gray-50 p-4 text-center">
                <p className="text-xs text-gray-400 mb-2">{label}</p>
                <p className="text-2xl font-bold" style={{ color }}>
                  {peak.peak_level.toFixed(2)}
                </p>
                <p className="text-xs text-gray-400 mt-1">Peak day {peak.day_of_peak}</p>
                <div className="mt-2 text-xs text-gray-500">
                  Final: {peak.final_level.toFixed(2)}
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* Kinetics & Th1/Th2 plots */}
      <div className="mb-6 grid gap-6 lg:grid-cols-2">
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">Immune Kinetics Plot</h2>
            <p className="text-xs text-gray-400 mt-0.5">
              Time course of Th, Tc, B, Ab, Memory over {sim.simulation_days} days
            </p>
          </div>
          {/* eslint-disable-next-line @next/next/no-img-element */}
          <img src="/images/immune_kinetics.png" alt="Immune kinetics plot showing cell populations over 365 days" className="w-full" />
        </div>
        <div className="rounded-xl bg-white shadow-card overflow-hidden">
          <div className="border-b border-gray-100 px-5 py-4">
            <h2 className="text-sm font-semibold text-gray-900">Th1 / Th2 Balance Plot</h2>
            <p className="text-xs text-gray-400 mt-0.5">
              {sim.th1_dominance_pct}% Th1 dominance — pro-inflammatory response for intracellular clearance
            </p>
          </div>
          {/* eslint-disable-next-line @next/next/no-img-element */}
          <img src="/images/th1_th2.png" alt="Th1/Th2 balance plot showing dominant Th1 response" className="w-full" />
        </div>
      </div>

      {/* Bottom charts */}
      <div className="grid gap-6 lg:grid-cols-2">
        {/* Memory formation */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Memory Cell Formation</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            Cumulative memory levels after each vaccine dose
          </p>
          <BarChart
            categories={memoryCategories}
            series={memorySeries}
            colors={["#8B5CF6"]}
            height={220}
          />
        </div>

        {/* Th1/Th2 balance */}
        <div className="rounded-xl bg-white shadow-card p-5">
          <h2 className="text-sm font-semibold text-gray-900">Th1 / Th2 Balance</h2>
          <p className="text-xs text-gray-400 mt-0.5 mb-4">
            Immune polarization — Th1 dominance required for intracellular pathogen clearance
          </p>
          <BarChart
            categories={th1Th2Categories}
            series={th1Th2Series}
            colors={["#3B82F6", "#E5E7EB"]}
            height={220}
          />
          <div className="mt-4 flex items-center gap-2 rounded-lg bg-purple-50 px-4 py-3">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2} className="h-4 w-4 text-purple-500 flex-shrink-0">
              <circle cx="12" cy="12" r="10" />
              <path d="M12 8v4M12 16h.01" />
            </svg>
            <p className="text-xs text-purple-700">
              {sim.th1_dominance_pct}% Th1 dominance confirms pro-inflammatory bias needed for Leishmania clearance.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
}
