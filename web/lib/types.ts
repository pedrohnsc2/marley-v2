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

export interface DrugTarget {
  gene_id: string;
  gene_name: string;
  pathway: string;
  identity_score: number;
  druggability_score: number;
  priority: boolean;
}

export interface DockingResult {
  target_gene_name: string;
  compound_id: string;
  compound_name: string;
  binding_affinity: number;
  composite_score: number;
  is_approved_drug: boolean;
}

export interface Epitope {
  sequence: string;
  source_gene_name: string;
  epitope_type: string;
  allele: string;
  ic50: number;
}

export interface KpiData {
  proteinsAnalyzed: number;
  epitopesSelected: number;
  compoundsDocked: number;
  bestAffinity: number;
  drugTargets: number;
  customMolecules: number;
}
