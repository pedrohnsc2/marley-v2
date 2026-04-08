"""
Domain models for vaccine candidate data.

The central entity is ``Candidate``, which tracks a single gene through
every stage of the reverse vaccinology pipeline -- from raw sequence to
final scored and filtered result.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Scoring weights
CONSERVATION_WEIGHT: float = 0.4
IMMUNOGENICITY_WEIGHT: float = 0.6

# Valid status values
STATUS_PENDING: str = "pending"
STATUS_APPROVED: str = "approved"
STATUS_REJECTED: str = "rejected"
VALID_STATUSES: set[str] = {STATUS_PENDING, STATUS_APPROVED, STATUS_REJECTED}


@dataclass
class Candidate:
    """A vaccine candidate gene progressing through the pipeline.

    Attributes:
        gene_id: Unique identifier for the gene (e.g. NCBI accession).
        gene_name: Human-readable gene name.
        sequence: Amino-acid or nucleotide sequence.
        has_signal_peptide: Whether a signal peptide was detected.
        conservation_score: Conservation score in the range [0.0, 1.0].
        immunogenicity_score: Predicted immunogenicity in [0.0, 1.0].
        final_score: Weighted average of conservation and immunogenicity.
        filters_passed: Names of pipeline filters this candidate cleared.
        status: Current pipeline status (pending, approved, rejected).
        priority: Whether this is a pre-validated antigen with priority status.
        source: Origin of the antigen data (e.g. researcher/institution).
        evidence: Summary of experimental evidence supporting this candidate.
    """

    gene_id: str
    gene_name: str
    sequence: str
    has_signal_peptide: bool = False
    conservation_score: float = 0.0
    immunogenicity_score: float = 0.0
    final_score: float = 0.0
    filters_passed: list[str] = field(default_factory=list)
    status: str = STATUS_PENDING
    priority: bool = False
    source: str = ""
    evidence: str = ""

    def compute_final_score(self) -> None:
        """Set ``final_score`` as the weighted average of sub-scores.

        Formula: 40% conservation + 60% immunogenicity.
        """
        self.final_score = (
            CONSERVATION_WEIGHT * self.conservation_score
            + IMMUNOGENICITY_WEIGHT * self.immunogenicity_score
        )

    def to_dict(self) -> dict:
        """Serialize to a plain dict suitable for Supabase upsert.

        Returns:
            Dictionary with all candidate fields.
        """
        return {
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "sequence": self.sequence,
            "has_signal_peptide": self.has_signal_peptide,
            "conservation_score": self.conservation_score,
            "immunogenicity_score": self.immunogenicity_score,
            "final_score": self.final_score,
            "filters_passed": self.filters_passed,
            "status": self.status,
            "priority": self.priority,
            "source": self.source,
            "evidence": self.evidence,
        }

    @classmethod
    def from_dict(cls, data: dict) -> Candidate:
        """Construct a ``Candidate`` from a dictionary (e.g. Supabase row).

        Args:
            data: Dictionary containing candidate fields.  Unknown keys
                  are silently ignored.

        Returns:
            A new ``Candidate`` instance.
        """
        return cls(
            gene_id=data["gene_id"],
            gene_name=data["gene_name"],
            sequence=data["sequence"],
            has_signal_peptide=data.get("has_signal_peptide", False),
            conservation_score=data.get("conservation_score", 0.0),
            immunogenicity_score=data.get("immunogenicity_score", 0.0),
            final_score=data.get("final_score", 0.0),
            filters_passed=data.get("filters_passed", []),
            status=data.get("status", STATUS_PENDING),
            priority=data.get("priority", False),
            source=data.get("source", ""),
            evidence=data.get("evidence", ""),
        )


# ---------------------------------------------------------------------------
# Module 06 -- mRNA vaccine construct models
# ---------------------------------------------------------------------------


@dataclass
class Epitope:
    """A single predicted epitope selected for the vaccine construct."""

    sequence: str
    source_gene_id: str
    source_gene_name: str
    epitope_type: str  # "CTL" or "HTL"
    allele: str
    ic50: float
    rank: float = 0.0
    start_position: int = 0

    def to_dict(self) -> dict:
        """Serialize to a plain dict suitable for Supabase upsert.

        Returns:
            Dictionary with all epitope fields.
        """
        return {
            "sequence": self.sequence,
            "source_gene_id": self.source_gene_id,
            "source_gene_name": self.source_gene_name,
            "epitope_type": self.epitope_type,
            "allele": self.allele,
            "ic50": self.ic50,
            "rank": self.rank,
            "start_position": self.start_position,
        }

    @classmethod
    def from_dict(cls, data: dict) -> Epitope:
        """Construct an ``Epitope`` from a dictionary (e.g. Supabase row).

        Args:
            data: Dictionary containing epitope fields.

        Returns:
            A new ``Epitope`` instance.
        """
        return cls(
            sequence=data["sequence"],
            source_gene_id=data["source_gene_id"],
            source_gene_name=data.get("source_gene_name", ""),
            epitope_type=data["epitope_type"],
            allele=data["allele"],
            ic50=float(data["ic50"]),
            rank=float(data.get("rank", 0.0)),
            start_position=int(data.get("start_position", 0)),
        )


@dataclass
class VaccineConstruct:
    """A complete multi-epitope mRNA vaccine design."""

    construct_id: str
    protein_sequence: str
    mrna_sequence: str
    signal_peptide_name: str
    adjuvant_name: str
    epitopes: list[Epitope] = field(default_factory=list)
    molecular_weight: float = 0.0
    isoelectric_point: float = 0.0
    instability_index: float = 0.0
    gravy: float = 0.0
    gc_content: float = 0.0
    vaxijen_score: float | None = None
    allergenicity: str | None = None
    created_at: str = ""

    def to_dict(self) -> dict:
        """Serialize to a plain dict suitable for Supabase upsert.

        Returns:
            Dictionary with all construct fields.
        """
        return {
            "construct_id": self.construct_id,
            "protein_sequence": self.protein_sequence,
            "mrna_sequence": self.mrna_sequence,
            "signal_peptide": self.signal_peptide_name,
            "adjuvant_name": self.adjuvant_name,
            "epitope_count": len(self.epitopes),
            "molecular_weight": self.molecular_weight,
            "isoelectric_point": self.isoelectric_point,
            "instability_index": self.instability_index,
            "gravy": self.gravy,
            "gc_content": self.gc_content,
            "vaxijen_score": self.vaxijen_score,
            "allergenicity": self.allergenicity,
            "created_at": self.created_at,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "VaccineConstruct":
        """Construct a ``VaccineConstruct`` from a dictionary (e.g. Supabase row).

        Args:
            data: Dictionary containing construct fields.

        Returns:
            A new ``VaccineConstruct`` instance.
        """
        return cls(
            construct_id=data["construct_id"],
            protein_sequence=data["protein_sequence"],
            mrna_sequence=data["mrna_sequence"],
            signal_peptide_name=data.get("signal_peptide", "tPA"),
            adjuvant_name=data.get("adjuvant_name", "L7L12"),
            molecular_weight=float(data.get("molecular_weight", 0.0)),
            isoelectric_point=float(data.get("isoelectric_point", 0.0)),
            instability_index=float(data.get("instability_index", 0.0)),
            gravy=float(data.get("gravy", 0.0)),
            gc_content=float(data.get("gc_content", 0.0)),
            vaxijen_score=data.get("vaxijen_score"),
            allergenicity=data.get("allergenicity"),
            created_at=data.get("created_at", ""),
        )


# ---------------------------------------------------------------------------
# Module v2 -- Drug target discovery models
# ---------------------------------------------------------------------------

# Druggability score weights
DIVERGENCE_WEIGHT: float = 0.40
ACTIVE_SITE_WEIGHT: float = 0.35
ESSENTIALITY_WEIGHT: float = 0.25


@dataclass
class DrugTarget:
    """An enzymatic drug target from *L. infantum* evaluated for selective inhibitor design.

    Attributes:
        gene_id: Unique identifier for the gene (e.g. TriTrypDB accession).
        gene_name: Human-readable gene/enzyme name.
        sequence: Amino-acid sequence.
        enzyme_class: Enzyme classification (e.g. "phosphoribosyltransferase").
        pathway: Metabolic pathway (e.g. "purine_salvage", "sterol_biosynthesis").
        human_homolog_id: UniProt ID of the closest human homolog.
        identity_score: Percent identity with human homolog in [0.0, 1.0].
        active_site_diff: Percent difference in active-site residues in [0.0, 1.0].
        is_essential: Whether the gene is essential for parasite survival.
        druggability_score: Composite score in [0.0, 1.0] (higher = better target).
        alphafold_url: AlphaFold structure URL for the target.
        evidence: Bibliographic reference or experimental evidence summary.
        priority: True for targets with prior experimental validation.
        status: Pipeline status (pending, approved, rejected, priority_validated).
    """

    gene_id: str
    gene_name: str
    sequence: str
    enzyme_class: str = ""
    pathway: str = ""
    human_homolog_id: str = ""
    identity_score: float = 0.0
    active_site_diff: float = 0.0
    is_essential: bool = False
    druggability_score: float = 0.0
    alphafold_url: str = ""
    evidence: str = ""
    priority: bool = False
    status: str = STATUS_PENDING

    def compute_druggability_score(self) -> None:
        """Calculate the composite druggability score.

        Formula:
            (1 - identity_score) * 0.40
            + active_site_diff * 0.35
            + essentiality_factor * 0.25

        Where essentiality_factor is 1.0 if essential, 0.3 otherwise.
        """
        essentiality_factor = 1.0 if self.is_essential else 0.3
        self.druggability_score = (
            (1 - self.identity_score) * DIVERGENCE_WEIGHT
            + self.active_site_diff * ACTIVE_SITE_WEIGHT
            + essentiality_factor * ESSENTIALITY_WEIGHT
        )

    def to_dict(self) -> dict:
        """Serialize to a plain dict suitable for Supabase upsert.

        Returns:
            Dictionary with all drug target fields.
        """
        return {
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "sequence": self.sequence,
            "enzyme_class": self.enzyme_class,
            "pathway": self.pathway,
            "human_homolog_id": self.human_homolog_id,
            "identity_score": self.identity_score,
            "active_site_diff": self.active_site_diff,
            "is_essential": self.is_essential,
            "druggability_score": self.druggability_score,
            "alphafold_url": self.alphafold_url,
            "evidence": self.evidence,
            "priority": self.priority,
            "status": self.status,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "DrugTarget":
        """Construct a ``DrugTarget`` from a dictionary (e.g. Supabase row).

        Args:
            data: Dictionary containing drug target fields. Unknown keys
                  are silently ignored.

        Returns:
            A new ``DrugTarget`` instance.
        """
        return cls(
            gene_id=data["gene_id"],
            gene_name=data["gene_name"],
            sequence=data.get("sequence", ""),
            enzyme_class=data.get("enzyme_class", ""),
            pathway=data.get("pathway", ""),
            human_homolog_id=data.get("human_homolog_id", ""),
            identity_score=float(data.get("identity_score", 0.0)),
            active_site_diff=float(data.get("active_site_diff", 0.0)),
            is_essential=data.get("is_essential", False),
            druggability_score=float(data.get("druggability_score", 0.0)),
            alphafold_url=data.get("alphafold_url", ""),
            evidence=data.get("evidence", ""),
            priority=data.get("priority", False),
            status=data.get("status", STATUS_PENDING),
        )


# ---------------------------------------------------------------------------
# Module v3 -- Molecular docking models
# ---------------------------------------------------------------------------


@dataclass
class DockingCompound:
    """A small molecule compound prepared for molecular docking.

    Attributes:
        compound_id: Unique identifier (ChEMBL ID, ZINC ID, or DrugBank ID).
        name: Common drug name (if approved drug).
        smiles: Canonical SMILES string.
        source: Origin database ("chembl", "drugbank", "zinc", "repurposing_hub").
        is_approved_drug: Whether this compound is an approved drug.
        molecular_weight: Molecular weight in Daltons.
        logp: Calculated LogP (lipophilicity).
        inchi_key: InChIKey for deduplication.
    """

    compound_id: str
    name: str
    smiles: str
    source: str = ""
    is_approved_drug: bool = False
    molecular_weight: float = 0.0
    logp: float = 0.0
    inchi_key: str = ""

    def to_dict(self) -> dict:
        return {
            "compound_id": self.compound_id,
            "name": self.name,
            "smiles": self.smiles,
            "source": self.source,
            "is_approved": self.is_approved_drug,
            "mol_weight": self.molecular_weight,
            "logp": self.logp,
            "inchi_key": self.inchi_key,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "DockingCompound":
        return cls(
            compound_id=data["compound_id"],
            name=data.get("name", ""),
            smiles=data["smiles"],
            source=data.get("source", ""),
            is_approved_drug=data.get("is_approved", False),
            molecular_weight=float(data.get("mol_weight", 0.0)),
            logp=float(data.get("logp", 0.0)),
            inchi_key=data.get("inchi_key", ""),
        )


# Docking score weights
AFFINITY_WEIGHT: float = 0.50
ADMET_WEIGHT: float = 0.25
REPURPOSING_WEIGHT: float = 0.15
SELECTIVITY_WEIGHT: float = 0.10


@dataclass
class DockingResult:
    """Result of a molecular docking simulation.

    Attributes:
        target_gene_id: Gene ID of the drug target.
        target_gene_name: Human-readable target name.
        compound_id: ID of the docked compound.
        compound_name: Common name of the compound.
        smiles: SMILES of the compound.
        binding_affinity: Vina binding affinity in kcal/mol (more negative = better).
        rmsd_lb: RMSD lower bound from Vina.
        rmsd_ub: RMSD upper bound from Vina.
        lipinski_violations: Number of Lipinski Rule of 5 violations (0-4).
        is_approved_drug: Whether compound is an approved drug.
        admet_score: ADMET composite score in [0.0, 1.0].
        composite_score: Final ranking score combining docking + ADMET.
        pdbqt_path: Path to output PDBQT with docked pose.
        source: Origin of the compound.
        status: Pipeline status.
    """

    target_gene_id: str
    target_gene_name: str
    compound_id: str
    compound_name: str = ""
    smiles: str = ""
    binding_affinity: float = 0.0
    rmsd_lb: float = 0.0
    rmsd_ub: float = 0.0
    lipinski_violations: int = 0
    is_approved_drug: bool = False
    admet_score: float = 0.0
    composite_score: float = 0.0
    pdbqt_path: str = ""
    source: str = ""
    status: str = STATUS_PENDING

    def compute_composite_score(self) -> None:
        """Calculate composite score from docking affinity and ADMET.

        Normalizes binding affinity from [-12, 0] kcal/mol to [0, 1],
        then combines with ADMET score and repurposing bonus.
        """
        norm_affinity = min(1.0, max(0.0, -self.binding_affinity / 12.0))
        repurposing_bonus = 1.0 if self.is_approved_drug else 0.0
        self.composite_score = (
            norm_affinity * AFFINITY_WEIGHT
            + self.admet_score * ADMET_WEIGHT
            + repurposing_bonus * REPURPOSING_WEIGHT
            + SELECTIVITY_WEIGHT
        )

    def to_dict(self) -> dict:
        return {
            "target_gene_id": self.target_gene_id,
            "target_gene_name": self.target_gene_name,
            "compound_id": self.compound_id,
            "compound_name": self.compound_name,
            "smiles": self.smiles,
            "binding_affinity": self.binding_affinity,
            "rmsd_lb": self.rmsd_lb,
            "rmsd_ub": self.rmsd_ub,
            "lipinski_violations": self.lipinski_violations,
            "is_approved_drug": self.is_approved_drug,
            "admet_score": self.admet_score,
            "composite_score": self.composite_score,
            "pdbqt_path": self.pdbqt_path,
            "source": self.source,
            "status": self.status,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "DockingResult":
        return cls(
            target_gene_id=data["target_gene_id"],
            target_gene_name=data.get("target_gene_name", ""),
            compound_id=data["compound_id"],
            compound_name=data.get("compound_name", ""),
            smiles=data.get("smiles", ""),
            binding_affinity=float(data.get("binding_affinity", 0.0)),
            rmsd_lb=float(data.get("rmsd_lb", 0.0)),
            rmsd_ub=float(data.get("rmsd_ub", 0.0)),
            lipinski_violations=int(data.get("lipinski_violations", 0)),
            is_approved_drug=data.get("is_approved_drug", False),
            admet_score=float(data.get("admet_score", 0.0)),
            composite_score=float(data.get("composite_score", 0.0)),
            pdbqt_path=data.get("pdbqt_path", ""),
            source=data.get("source", ""),
            status=data.get("status", "pending"),
        )


# ---------------------------------------------------------------------------
# Module v4 -- Vaccine optimization models
# ---------------------------------------------------------------------------


@dataclass
class OptimizedEpitope:
    """An epitope variant generated by APL (Altered Peptide Ligand) design."""

    original_sequence: str
    optimized_sequence: str
    source_gene_id: str
    epitope_type: str
    allele: str
    original_ic50: float
    optimized_ic50: float
    improvement_ratio: float
    mutation_positions: list[int] = field(default_factory=list)
    is_safe: bool = True
    blast_identity: float = 0.0
    blast_hit_gene: str = ""

    def to_dict(self) -> dict:
        return {
            "original_sequence": self.original_sequence,
            "optimized_sequence": self.optimized_sequence,
            "source_gene_id": self.source_gene_id,
            "epitope_type": self.epitope_type,
            "allele": self.allele,
            "original_ic50": self.original_ic50,
            "optimized_ic50": self.optimized_ic50,
            "improvement_ratio": self.improvement_ratio,
            "mutation_positions": self.mutation_positions,
            "is_safe": self.is_safe,
            "blast_identity": self.blast_identity,
            "blast_hit_gene": self.blast_hit_gene,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "OptimizedEpitope":
        return cls(
            original_sequence=data["original_sequence"],
            optimized_sequence=data["optimized_sequence"],
            source_gene_id=data.get("source_gene_id", ""),
            epitope_type=data.get("epitope_type", ""),
            allele=data.get("allele", ""),
            original_ic50=float(data.get("original_ic50", 0.0)),
            optimized_ic50=float(data.get("optimized_ic50", 0.0)),
            improvement_ratio=float(data.get("improvement_ratio", 0.0)),
            mutation_positions=data.get("mutation_positions", []),
            is_safe=data.get("is_safe", True),
            blast_identity=float(data.get("blast_identity", 0.0)),
            blast_hit_gene=data.get("blast_hit_gene", ""),
        )


@dataclass
class AdjuvantResult:
    """Result of adjuvant screening for a vaccine construct."""

    adjuvant_name: str
    adjuvant_sequence: str
    vaxijen_score: float | None = None
    is_antigenic: bool = False
    th_bias: str = ""
    construct_length: int = 0
    notes: str = ""

    def to_dict(self) -> dict:
        return {
            "adjuvant_name": self.adjuvant_name,
            "adjuvant_sequence": self.adjuvant_sequence,
            "vaxijen_score": self.vaxijen_score,
            "is_antigenic": self.is_antigenic,
            "th_bias": self.th_bias,
            "construct_length": self.construct_length,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "AdjuvantResult":
        return cls(
            adjuvant_name=data["adjuvant_name"],
            adjuvant_sequence=data.get("adjuvant_sequence", ""),
            vaxijen_score=data.get("vaxijen_score"),
            is_antigenic=data.get("is_antigenic", False),
            th_bias=data.get("th_bias", ""),
            construct_length=int(data.get("construct_length", 0)),
            notes=data.get("notes", ""),
        )


@dataclass
class OptimizedConstruct:
    """A v4 optimized vaccine construct for comparison."""

    construct_id: str
    base_construct_id: str
    adjuvant_name: str
    epitope_order: str
    protein_sequence: str
    mrna_sequence: str
    epitope_count: int = 0
    optimized_epitope_count: int = 0
    safe_epitope_count: int = 0
    molecular_weight: float = 0.0
    instability_index: float = 0.0
    plddt_mean: float = 0.0
    plddt_min: float = 0.0
    vaxijen_score: float | None = None
    gc_content: float = 0.0
    improvement_summary: str = ""

    def to_dict(self) -> dict:
        return {
            "construct_id": self.construct_id,
            "base_construct_id": self.base_construct_id,
            "adjuvant_name": self.adjuvant_name,
            "epitope_order": self.epitope_order,
            "protein_sequence": self.protein_sequence,
            "mrna_sequence": self.mrna_sequence,
            "epitope_count": self.epitope_count,
            "optimized_epitope_count": self.optimized_epitope_count,
            "safe_epitope_count": self.safe_epitope_count,
            "molecular_weight": self.molecular_weight,
            "instability_index": self.instability_index,
            "plddt_mean": self.plddt_mean,
            "plddt_min": self.plddt_min,
            "vaxijen_score": self.vaxijen_score,
            "gc_content": self.gc_content,
            "improvement_summary": self.improvement_summary,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "OptimizedConstruct":
        return cls(
            construct_id=data["construct_id"],
            base_construct_id=data.get("base_construct_id", ""),
            adjuvant_name=data.get("adjuvant_name", ""),
            epitope_order=data.get("epitope_order", ""),
            protein_sequence=data.get("protein_sequence", ""),
            mrna_sequence=data.get("mrna_sequence", ""),
            epitope_count=int(data.get("epitope_count", 0)),
            optimized_epitope_count=int(data.get("optimized_epitope_count", 0)),
            safe_epitope_count=int(data.get("safe_epitope_count", 0)),
            molecular_weight=float(data.get("molecular_weight", 0.0)),
            instability_index=float(data.get("instability_index", 0.0)),
            plddt_mean=float(data.get("plddt_mean", 0.0)),
            plddt_min=float(data.get("plddt_min", 0.0)),
            vaxijen_score=data.get("vaxijen_score"),
            gc_content=float(data.get("gc_content", 0.0)),
            improvement_summary=data.get("improvement_summary", ""),
        )


# ---------------------------------------------------------------------------
# Module v4-RNA -- RNA information theory models
# ---------------------------------------------------------------------------

# Information score weights
ENTROPY_DELTA_WEIGHT: float = 0.35
CONSERVATION_WEIGHT_RNA: float = 0.25
CODON_BIAS_WEIGHT: float = 0.20
SL_RNA_WEIGHT: float = 0.10
STRUCTURE_WEIGHT: float = 0.10


@dataclass
class RNATarget:
    """An RNA target identified by Shannon entropy analysis.

    Regions of low entropy in *L. infantum* that have high entropy in
    humans are mathematically ideal targets — conserved in the parasite,
    variable in the host.

    Attributes:
        gene_id: Unique identifier for the gene.
        gene_name: Human-readable gene name.
        sequence_rna: RNA sequence.
        gc_content: Percentage of G+C nucleotides (expected ~60-63% in L. infantum).
        shannon_entropy: Shannon entropy H(X) in bits (lower = more conserved).
        human_entropy: Equivalent Shannon entropy in human homolog.
        entropy_delta: Difference: human_entropy - shannon_entropy (higher = better target).
        codon_bias_score: RSCU distance from human codon usage, normalized to [0.0, 1.0].
        has_sl_rna: Whether the transcript has the 39nt spliced leader sequence.
        min_free_energy: Minimum free energy of RNA secondary structure (kcal/mol).
        conservation_score: Conservation across Leishmania strains in [0.0, 1.0].
        information_score: Composite score in [0.0, 1.0] (higher = better target).
        is_priority: True for targets with prior experimental validation.
        evidence: Bibliographic reference or evidence summary.
        status: Pipeline status (pending, approved, rejected, priority_validated).
    """

    gene_id: str
    gene_name: str
    sequence_rna: str = ""
    gc_content: float = 0.0
    shannon_entropy: float = 0.0
    human_entropy: float = 0.0
    entropy_delta: float = 0.0
    codon_bias_score: float = 0.0
    has_sl_rna: bool = False
    min_free_energy: float = 0.0
    conservation_score: float = 0.0
    information_score: float = 0.0
    is_priority: bool = False
    evidence: str = ""
    status: str = STATUS_PENDING

    def compute_information_score(
        self, max_delta: float = 2.0, max_mfe: float = 50.0
    ) -> None:
        """Calculate composite information score.

        Formula:
            (entropy_delta / max_delta) * 0.35
            + (1 - shannon_entropy) * 0.25
            + codon_bias_score * 0.20
            + (1.0 if has_sl_rna else 0.0) * 0.10
            + (abs(min_free_energy) / max_mfe) * 0.10
        """
        delta_norm = min(1.0, max(0.0, self.entropy_delta / max_delta))
        conservation = min(1.0, max(0.0, 1.0 - self.shannon_entropy))
        mfe_norm = min(1.0, abs(self.min_free_energy) / max_mfe) if max_mfe > 0 else 0.0
        sl_bonus = 1.0 if self.has_sl_rna else 0.0

        self.information_score = (
            delta_norm * ENTROPY_DELTA_WEIGHT
            + conservation * CONSERVATION_WEIGHT_RNA
            + self.codon_bias_score * CODON_BIAS_WEIGHT
            + sl_bonus * SL_RNA_WEIGHT
            + mfe_norm * STRUCTURE_WEIGHT
        )

    def to_dict(self) -> dict:
        """Serialize to a plain dict suitable for Supabase upsert."""
        return {
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "sequence_rna": self.sequence_rna,
            "gc_content": self.gc_content,
            "shannon_entropy": self.shannon_entropy,
            "human_entropy": self.human_entropy,
            "entropy_delta": self.entropy_delta,
            "codon_bias_score": self.codon_bias_score,
            "has_sl_rna": self.has_sl_rna,
            "min_free_energy": self.min_free_energy,
            "conservation_score": self.conservation_score,
            "information_score": self.information_score,
            "is_priority": self.is_priority,
            "evidence": self.evidence,
            "status": self.status,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "RNATarget":
        """Construct an ``RNATarget`` from a dictionary."""
        return cls(
            gene_id=data["gene_id"],
            gene_name=data["gene_name"],
            sequence_rna=data.get("sequence_rna", ""),
            gc_content=float(data.get("gc_content", 0.0)),
            shannon_entropy=float(data.get("shannon_entropy", 0.0)),
            human_entropy=float(data.get("human_entropy", 0.0)),
            entropy_delta=float(data.get("entropy_delta", 0.0)),
            codon_bias_score=float(data.get("codon_bias_score", 0.0)),
            has_sl_rna=data.get("has_sl_rna", False),
            min_free_energy=float(data.get("min_free_energy", 0.0)),
            conservation_score=float(data.get("conservation_score", 0.0)),
            information_score=float(data.get("information_score", 0.0)),
            is_priority=data.get("is_priority", False),
            evidence=data.get("evidence", ""),
            status=data.get("status", STATUS_PENDING),
        )
