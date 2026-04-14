/**
 * Static metadata about each pipeline, used by the wizard and other UI components.
 * Provides humanized names, descriptions, expected outputs, and stage information.
 */

export interface PipelineMetadata {
  id: string;
  displayName: string;
  description: string;
  expectedOutputs: string[];
  estimatedDuration: string;
  icon: string;
  stages: { stageId: string; friendlyName: string; description: string }[];
}

export const PIPELINE_METADATA: Record<string, PipelineMetadata> = {
  vaccine: {
    id: "vaccine",
    displayName: "Vaccine Candidate Design",
    description:
      "Designs a multi-epitope mRNA vaccine construct by predicting MHC-I binding epitopes from the parasite proteome.",
    expectedOutputs: [
      "Ranked epitope list with binding affinities",
      "Full vaccine construct sequence",
      "3D structure prediction",
      "Codon-optimized mRNA sequence",
    ],
    estimatedDuration: "3-8 min",
    icon: "\u{1F9EC}",
    stages: [
      {
        stageId: "01_fetch_genome",
        friendlyName: "Downloading proteome data",
        description: "Fetching protein sequences from TriTrypDB",
      },
      {
        stageId: "02_filter_surface",
        friendlyName: "Filtering surface proteins",
        description: "Identifying membrane and secreted proteins",
      },
      {
        stageId: "03_conservation",
        friendlyName: "Analyzing conservation",
        description: "Comparing across Leishmania species",
      },
      {
        stageId: "04_immunogenicity",
        friendlyName: "Predicting immunogenicity",
        description: "Running epitope prediction and scoring",
      },
      {
        stageId: "05_report",
        friendlyName: "Generating report",
        description: "Compiling results and visualizations",
      },
    ],
  },
  drug: {
    id: "drug",
    displayName: "Drug Target Discovery",
    description:
      "Identifies parasite-specific proteins that could serve as drug targets by comparing parasite and host proteomes.",
    expectedOutputs: [
      "Prioritized target list with druggability scores",
      "Pathway analysis",
      "Human ortholog comparison",
    ],
    estimatedDuration: "2-5 min",
    icon: "\u{1F48A}",
    stages: [
      {
        stageId: "01_fetch_enzymes",
        friendlyName: "Fetching enzyme data",
        description: "Downloading parasite enzyme database",
      },
      {
        stageId: "02_human_comparison",
        friendlyName: "Comparing with human proteome",
        description: "Finding parasite-specific proteins",
      },
      {
        stageId: "03_essentiality",
        friendlyName: "Assessing essentiality",
        description: "Evaluating target importance for parasite survival",
      },
      {
        stageId: "04_druggability",
        friendlyName: "Scoring druggability",
        description: "Evaluating each target's drug potential",
      },
      {
        stageId: "05_report",
        friendlyName: "Generating report",
        description: "Compiling results",
      },
    ],
  },
  docking: {
    id: "docking",
    displayName: "Molecular Docking",
    description:
      "Screens small-molecule compounds against drug targets using AutoDock Vina to find potential inhibitors.",
    expectedOutputs: [
      "Binding affinity rankings",
      "Top compound-target pairs",
      "3D binding pose visualization",
    ],
    estimatedDuration: "5-15 min",
    icon: "\u{1F52C}",
    stages: [
      {
        stageId: "01_fetch_structures",
        friendlyName: "Fetching structures",
        description: "Downloading 3D protein structures",
      },
      {
        stageId: "02_compound_library",
        friendlyName: "Preparing compounds",
        description: "Loading and processing compound library",
      },
      {
        stageId: "03_docking",
        friendlyName: "Running docking simulation",
        description: "AutoDock Vina virtual screening",
      },
      {
        stageId: "04_admet",
        friendlyName: "ADMET filtering",
        description: "Evaluating drug-likeness properties",
      },
      {
        stageId: "05_report",
        friendlyName: "Generating report",
        description: "Compiling docking results",
      },
    ],
  },
  rna: {
    id: "rna",
    displayName: "RNA Target Validation",
    description:
      "Analyzes RNA sequence conservation and structural entropy to validate therapeutic targets across kinetoplastid species.",
    expectedOutputs: [
      "Conservation heatmap",
      "Entropy scores per transcript",
      "Cross-species comparison",
    ],
    estimatedDuration: "1-3 min",
    icon: "\u{1F9EA}",
    stages: [
      {
        stageId: "01_transcriptome",
        friendlyName: "Loading transcriptome",
        description: "Fetching RNA sequence data",
      },
      {
        stageId: "02_codon_usage",
        friendlyName: "Analyzing codon usage",
        description: "Computing codon frequency tables",
      },
      {
        stageId: "03_entropy",
        friendlyName: "Calculating entropy",
        description: "Measuring sequence variability",
      },
      {
        stageId: "04_sl_rna",
        friendlyName: "SL RNA analysis",
        description: "Analyzing spliced leader RNA",
      },
      {
        stageId: "05_human_comparison",
        friendlyName: "Human comparison",
        description: "Cross-referencing with human transcripts",
      },
      {
        stageId: "06_structure",
        friendlyName: "Structure prediction",
        description: "Predicting RNA secondary structure",
      },
      {
        stageId: "07_report",
        friendlyName: "Generating report",
        description: "Compiling analysis results",
      },
    ],
  },
  aso_math: {
    id: "aso_math",
    displayName: "ASO Mathematical Validation",
    description:
      "Runs mathematical validation of antisense oligonucleotide designs, testing thermodynamic stability, off-target risk, and escape mutation resistance.",
    expectedOutputs: [
      "Validation certificate with score",
      "Thermodynamic analysis",
      "Off-target assessment",
      "Escape mutation prediction",
    ],
    estimatedDuration: "1-2 min",
    icon: "\u{1F4D0}",
    stages: [
      {
        stageId: "01_thermodynamic",
        friendlyName: "Thermodynamic validation",
        description: "Checking ASO binding stability",
      },
      {
        stageId: "02_selectivity",
        friendlyName: "Selectivity analysis",
        description: "Screening for off-target binding",
      },
      {
        stageId: "03_conservation",
        friendlyName: "Conservation check",
        description: "Verifying target conservation across species",
      },
      {
        stageId: "04_optimization",
        friendlyName: "Optimization",
        description: "Finding optimal ASO modifications",
      },
      {
        stageId: "05_resistance",
        friendlyName: "Resistance analysis",
        description: "Predicting escape mutations",
      },
    ],
  },
};

/** Parameter humanization */
export interface ParameterMeta {
  label: string;
  tooltip: string;
  group: string;
  advanced?: boolean;
}

export const PARAMETER_METADATA: Record<string, ParameterMeta> = {
  organism: {
    label: "Organism",
    tooltip: "Target organism for analysis",
    group: "Input Data",
  },
  species_name: {
    label: "Species Name",
    tooltip: "Full species name",
    group: "Input Data",
  },
  proteome_source: {
    label: "Proteome Source",
    tooltip: "Database to fetch proteome from",
    group: "Input Data",
  },
  proteome_id: {
    label: "Proteome ID",
    tooltip: "Identifier in the proteome database",
    group: "Input Data",
  },
  conservation_weight: {
    label: "Conservation Weight",
    tooltip:
      "Weight given to evolutionary conservation in scoring (0-1)",
    group: "Scoring",
  },
  immunogenicity_weight: {
    label: "Immunogenicity Weight",
    tooltip:
      "Weight given to predicted immune response in scoring (0-1)",
    group: "Scoring",
  },
  mhc_alleles: {
    label: "MHC Alleles",
    tooltip: "Target MHC-I alleles for epitope prediction",
    group: "Analysis",
    advanced: true,
  },
  identity_threshold: {
    label: "Max Human Identity (%)",
    tooltip:
      "Proteins above this identity to human orthologs are excluded. Lower = stricter.",
    group: "Filtering",
  },
  top_n_targets: {
    label: "Top N Targets",
    tooltip: "Number of top-scoring targets to include in results",
    group: "Output",
  },
  host_organism: {
    label: "Host Organism",
    tooltip: "Host species for ortholog comparison",
    group: "Input Data",
    advanced: true,
  },
  force_rerun: {
    label: "Force Re-run",
    tooltip: "Re-run all stages even if cached results exist",
    group: "Execution",
    advanced: true,
  },
  organism_slug: {
    label: "Organism Slug",
    tooltip: "Machine-readable organism identifier",
    group: "Input Data",
    advanced: true,
  },
  taxonomic_group: {
    label: "Taxonomic Group",
    tooltip: "Taxonomic classification of the target organism",
    group: "Input Data",
    advanced: true,
  },
  disease_name: {
    label: "Disease Name",
    tooltip: "Associated disease for this analysis",
    group: "Input Data",
  },
  sl_sequence: {
    label: "SL Sequence",
    tooltip: "Spliced leader RNA sequence",
    group: "Analysis",
    advanced: true,
  },
  mutation_rate: {
    label: "Mutation Rate",
    tooltip: "Estimated mutation rate per base per generation",
    group: "Analysis",
    advanced: true,
  },
  generation_time_hours: {
    label: "Generation Time (h)",
    tooltip: "Organism generation time in hours",
    group: "Analysis",
    advanced: true,
  },
  sl_copy_number: {
    label: "SL Copy Number",
    tooltip: "Number of SL RNA gene copies",
    group: "Analysis",
    advanced: true,
  },
  dg_functional_threshold: {
    label: "dG Functional Threshold",
    tooltip: "Minimum free energy threshold for functional binding (kcal/mol)",
    group: "Analysis",
  },
  length_scan_min: {
    label: "Min ASO Length",
    tooltip: "Minimum ASO length for scanning",
    group: "Analysis",
  },
  length_scan_max: {
    label: "Max ASO Length",
    tooltip: "Maximum ASO length for scanning",
    group: "Analysis",
  },
  conservation_threshold: {
    label: "Conservation Threshold",
    tooltip: "Minimum conservation score for candidate selection (0-1)",
    group: "Filtering",
  },
};

/** Helper to get stage friendly name, falling back to raw name */
export function getStageFriendlyName(
  pipelineId: string,
  stageId: string,
): string | null {
  const pipeline = PIPELINE_METADATA[pipelineId];
  if (!pipeline) return null;
  const stage = pipeline.stages.find((s) => s.stageId === stageId);
  return stage?.friendlyName ?? null;
}
