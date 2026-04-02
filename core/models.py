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
