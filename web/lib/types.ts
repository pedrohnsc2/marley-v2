export interface Candidate {
  gene_id: string;
  gene_name: string;
  sequence: string;
  has_signal_peptide: boolean;
  conservation_score: number;
  immunogenicity_score: number;
  final_score: number;
  priority: "high" | "medium" | "low";
  source: string;
  evidence: string[];
  status: string;
  filters_passed: number;
}
