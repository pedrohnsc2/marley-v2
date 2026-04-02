import { Candidate } from "@/lib/types";

const MOCK_CANDIDATES: Candidate[] = [
  {
    gene_id: "LinJ.28.0310",
    gene_name: "LiHyp1",
    sequence: "MKLF...",
    has_signal_peptide: true,
    conservation_score: 0.95,
    immunogenicity_score: 0.91,
    final_score: 0.93,
    priority: "high",
    source: "TriTrypDB",
    evidence: ["Validated in clinical trial", "Strong Th1 response"],
    status: "validated",
    filters_passed: 6,
  },
  {
    gene_id: "LinJ.22.0680",
    gene_name: "A2",
    sequence: "MAST...",
    has_signal_peptide: false,
    conservation_score: 0.92,
    immunogenicity_score: 0.89,
    final_score: 0.9,
    priority: "high",
    source: "TriTrypDB",
    evidence: ["Used in Leish-Tec vaccine", "Amastigote-specific"],
    status: "validated",
    filters_passed: 5,
  },
  {
    gene_id: "LinJ.25.1680",
    gene_name: "KMP-11",
    sequence: "MATK...",
    has_signal_peptide: false,
    conservation_score: 0.88,
    immunogenicity_score: 0.85,
    final_score: 0.86,
    priority: "medium",
    source: "TriTrypDB",
    evidence: ["Kinetoplast membrane protein", "Cross-reactive"],
    status: "validated",
    filters_passed: 5,
  },
  {
    gene_id: "LinJ.12.0945",
    gene_name: "LACK",
    sequence: "MRVL...",
    has_signal_peptide: false,
    conservation_score: 0.84,
    immunogenicity_score: 0.82,
    final_score: 0.83,
    priority: "medium",
    source: "TriTrypDB",
    evidence: ["Leishmania homolog of RACK1", "Th1 adjuvant candidate"],
    status: "validated",
    filters_passed: 4,
  },
];

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

export default function Home() {
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
                <th className="px-4 py-3">Source</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-zinc-800">
              {MOCK_CANDIDATES.map((candidate, index) => (
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
                  <td className="px-4 py-3 text-zinc-400">
                    {candidate.source}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>

        <p className="mt-3 text-xs text-zinc-600">
          Showing {MOCK_CANDIDATES.length} candidates (mock data)
        </p>
      </section>
    </main>
  );
}
