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
    def from_dict(cls, data: dict) -> VaccineConstruct:
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
