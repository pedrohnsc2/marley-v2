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
