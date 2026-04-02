import { createServerClient } from "@/lib/supabase";
import { Candidate } from "@/lib/types";

function PriorityBadge({ priority }: { priority: string }) {
  const colors: Record<string, string> = {
    high: "bg-green-900 text-green-300 border-green-700",
    medium: "bg-yellow-900 text-yellow-300 border-yellow-700",
    low: "bg-zinc-800 text-zinc-400 border-zinc-600",
  };

  return (
    <span
      className={`inline-block rounded-full border px-2.5 py-0.5 text-xs font-medium ${colors[priority] ?? colors.low}`}
    >
      {priority}
    </span>
  );
}

async function fetchCandidates(): Promise<Candidate[]> {
  const supabase = createServerClient();

  const { data, error } = await supabase
    .from("candidates")
    .select("*")
    .order("priority", { ascending: false })
    .order("final_score", { ascending: false });

  if (error) {
    console.error("Failed to fetch candidates:", error.message);
    return [];
  }

  return (data as Candidate[]) ?? [];
}

export default async function Home() {
  const candidates = await fetchCandidates();

  return (
    <main className="mx-auto max-w-6xl px-6 py-12">
      <header className="mb-10">
        <h1 className="text-3xl font-bold tracking-tight text-white">
          Marley
        </h1>
        <p className="mt-1 text-sm text-zinc-400">
          Vaccine candidate antigen rankings for canine visceral leishmaniasis
        </p>
      </header>

      <section>
        <h2 className="mb-4 text-lg font-semibold text-zinc-200">
          Candidate Rankings
        </h2>

        {candidates.length === 0 ? (
          <div
            className="rounded-lg border border-zinc-800 bg-zinc-900/40 px-6 py-16 text-center"
            data-testid="empty-state"
          >
            <p className="text-sm text-zinc-400">
              No candidates found. Run the pipeline to populate results.
            </p>
          </div>
        ) : (
          <>
            <div className="overflow-x-auto rounded-lg border border-zinc-800">
              <table
                className="w-full text-left text-sm"
                data-testid="candidates-table"
              >
                <thead className="border-b border-zinc-800 bg-zinc-900/60 text-xs uppercase tracking-wider text-zinc-400">
                  <tr>
                    <th className="px-4 py-3">Rank</th>
                    <th className="px-4 py-3">Gene ID</th>
                    <th className="px-4 py-3">Gene Name</th>
                    <th className="px-4 py-3 text-right">Final Score</th>
                    <th className="px-4 py-3">Priority</th>
                    <th className="px-4 py-3">Evidence</th>
                    <th className="px-4 py-3">Source</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-zinc-800">
                  {candidates.map((candidate, index) => (
                    <tr
                      key={candidate.gene_id}
                      className="transition-colors hover:bg-zinc-900/40"
                      data-testid={`candidate-row-${index}`}
                    >
                      <td className="px-4 py-3 font-mono text-zinc-500">
                        {index + 1}
                      </td>
                      <td className="px-4 py-3 font-mono text-zinc-300">
                        {candidate.gene_id}
                      </td>
                      <td className="px-4 py-3 font-semibold text-white">
                        {candidate.gene_name}
                      </td>
                      <td className="px-4 py-3 text-right font-mono text-zinc-200">
                        {candidate.final_score.toFixed(2)}
                      </td>
                      <td className="px-4 py-3">
                        <PriorityBadge priority={candidate.priority} />
                      </td>
                      <td className="max-w-xs px-4 py-3 text-zinc-400">
                        {Array.isArray(candidate.evidence) &&
                        candidate.evidence.length > 0 ? (
                          <ul className="list-inside list-disc space-y-0.5">
                            {candidate.evidence.map((item) => (
                              <li key={item} className="truncate text-xs">
                                {item}
                              </li>
                            ))}
                          </ul>
                        ) : (
                          <span className="text-xs text-zinc-600">--</span>
                        )}
                      </td>
                      <td className="px-4 py-3 text-zinc-400">
                        {candidate.source}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            <p className="mt-3 text-xs text-zinc-600">
              Showing {candidates.length} candidate
              {candidates.length !== 1 ? "s" : ""} from Supabase
            </p>
          </>
        )}
      </section>
    </main>
  );
}
