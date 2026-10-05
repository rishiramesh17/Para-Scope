"""
Component ontology and evidence record schemas.

Core principle from CLAUDE.md and SPEC §6:
- Only compare like with like. Heavy chain against heavy chain, CDR-H3 against CDR-H3.
- Never CDR-H3 against CDR-L3; never VH against VL.
- Related-class comparisons (VH↔VHH, kappa↔lambda) must be explicitly labeled and down-ranked.
"""

from dataclasses import dataclass, field
from typing import Optional
from enum import Enum


class ComponentType(str, Enum):
    """Supported component types from SPEC §6.2-6.5."""
    # Variable domains
    VH = "VH"
    VL_KAPPA = "VL-kappa"
    VL_LAMBDA = "VL-lambda"
    VHH = "VHH"

    # CDRs
    CDR_H1 = "CDR-H1"
    CDR_H2 = "CDR-H2"
    CDR_H3 = "CDR-H3"
    CDR_L1 = "CDR-L1"
    CDR_L2 = "CDR-L2"
    CDR_L3 = "CDR-L3"

    # Framework regions
    FR_H1 = "FR-H1"
    FR_H2 = "FR-H2"
    FR_H3 = "FR-H3"
    FR_H4 = "FR-H4"
    FR_L1 = "FR-L1"
    FR_L2 = "FR-L2"
    FR_L3 = "FR-L3"
    FR_L4 = "FR-L4"

    # Constant regions (SPEC S6.5), as actually recoverable from the sources.
    #
    # SPEC S6.5 names CH1 / Hinge / CH2 / CH3 / CL-kappa / CL-lambda separately.
    # These two types are deliberately coarser: what the corpus yields is the
    # residues C-terminal to a numbered V-domain, and deciding which constant
    # domains that span contains needs a constant-region reference set that is
    # not bundled. Calling a ~100-residue remainder "CH1" because its length is
    # about right would be a claim, not a measurement -- so it is stored
    # undecomposed and named that way. Decomposition is future work; these types
    # are what honest ingestion can produce today.
    #
    # They are NOT full chains: a Fab heavy chain is V + CH1 and carries no
    # CH2/CH3, so FULL_HEAVY would overstate it.
    CONSTANT_HEAVY = "Constant region (heavy, undecomposed)"
    CONSTANT_LIGHT = "Constant region (light, undecomposed)"

    # Full chains
    FULL_HEAVY = "Full heavy chain"
    FULL_LIGHT = "Full light chain"


class ComponentRelationship(str, Enum):
    """Relationship between input and matched component (SPEC §20)."""
    EXACT_SAME_TYPE = "exact_same_type"
    RELATED_CLASS = "related_class"  # VH vs VHH, kappa vs lambda
    EXCLUDED = "excluded"  # Never compare these


class EvidenceStatus(str, Enum):
    """
    Outcome of attempting to gather one evidence dimension (SPEC §22).

    Exists so that "we looked and found nothing", "this does not apply here", and
    "we could not look" are three different statements rather than three zeroes.
    SPEC §22: *do not force inapplicable dimensions into a zero*; CLAUDE.md:
    *missing evidence is not zero evidence*. Any dimension carrying a status other
    than AVAILABLE must leave its value unset rather than defaulting to 0, to a
    midpoint, or to a neutral -- a failed database call must never lower a
    biological score.
    """

    AVAILABLE = "available"
    NOT_FOUND = "not_found"
    NOT_APPLICABLE = "not_applicable"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    ANALYSIS_FAILED = "analysis_failed"
    DATABASE_UNAVAILABLE = "database_unavailable"
    ACCESS_RESTRICTED = "access_restricted"
    UNSUPPORTED_COMPONENT = "unsupported_component"


class EvidenceClass(str, Enum):
    """Evidence classes from SPEC §17. Must be reported separately, never combined."""
    GERMLINE = "germline"
    OBSERVED_REPERTOIRE = "observed_repertoire"
    STRUCTURE_LINKED = "structure_linked_experimental"
    THERAPEUTIC_PRECEDENT = "therapeutic_precedent"
    DIRECT_PAIRING = "direct_paired_chain"
    USER_SUPPLIED = "user_supplied"
    COMPUTATIONAL_PREDICTION = "computational_prediction"


def compare_eligibility(input_type: ComponentType, matched_type: ComponentType) -> ComponentRelationship:
    """
    Determine whether two components are eligible to compare under SPEC §20.

    Returns:
        - EXACT_SAME_TYPE if they match exactly
        - RELATED_CLASS if they're related but not identical (e.g., VH vs VHH)
        - EXCLUDED if they should never be compared (e.g., CDR-H3 vs CDR-L3)
    """
    if input_type == matched_type:
        return ComponentRelationship.EXACT_SAME_TYPE

    # VH vs VHH: related class
    if (input_type == ComponentType.VH and matched_type == ComponentType.VHH) or \
       (input_type == ComponentType.VHH and matched_type == ComponentType.VH):
        return ComponentRelationship.RELATED_CLASS

    # Kappa vs lambda: related class
    if (input_type == ComponentType.VL_KAPPA and matched_type == ComponentType.VL_LAMBDA) or \
       (input_type == ComponentType.VL_LAMBDA and matched_type == ComponentType.VL_KAPPA):
        return ComponentRelationship.RELATED_CLASS

    # Heavy vs light: excluded
    heavy_cdr = {ComponentType.CDR_H1, ComponentType.CDR_H2, ComponentType.CDR_H3}
    light_cdr = {ComponentType.CDR_L1, ComponentType.CDR_L2, ComponentType.CDR_L3}
    heavy_fr = {ComponentType.FR_H1, ComponentType.FR_H2, ComponentType.FR_H3, ComponentType.FR_H4}
    light_fr = {ComponentType.FR_L1, ComponentType.FR_L2, ComponentType.FR_L3, ComponentType.FR_L4}

    if (input_type in heavy_cdr and matched_type in light_cdr) or \
       (input_type in light_cdr and matched_type in heavy_cdr) or \
       (input_type in heavy_fr and matched_type in light_fr) or \
       (input_type in light_fr and matched_type in heavy_fr):
        return ComponentRelationship.EXCLUDED

    # Heavy vs light constant region: excluded, same rule as heavy vs light
    # variable domains. A CH1 and a CL are different components, and their
    # mutual similarity would be a comparison across chain classes.
    if (input_type == ComponentType.CONSTANT_HEAVY and matched_type == ComponentType.CONSTANT_LIGHT) or \
       (input_type == ComponentType.CONSTANT_LIGHT and matched_type == ComponentType.CONSTANT_HEAVY):
        return ComponentRelationship.EXCLUDED

    # VH/VL comparisons: excluded
    vl_types = {ComponentType.VL_KAPPA, ComponentType.VL_LAMBDA}
    if (input_type == ComponentType.VH and matched_type in vl_types) or \
       (input_type in vl_types and matched_type == ComponentType.VH):
        return ComponentRelationship.EXCLUDED

    # Default: excluded (safe default per SPEC §3)
    return ComponentRelationship.EXCLUDED


@dataclass
class SourceSnapshot:
    """
    Immutable reference to a versioned source snapshot (SPEC §13, §14).

    Every evidence record points to a snapshot so that:
    - Old analyses remain reproducible after a source update.
    - Score changes across versions are explainable.
    """
    source_name: str
    snapshot_id: str
    retrieved_at: str  # ISO 8601 UTC
    license_id: str
    raw_checksum: str  # sha256 hex of raw source file
    record_count: int
    parser_version: str
    schema_version: str


@dataclass
class EvidenceRecord:
    """
    Immutable evidence record for one match (SPEC §38).

    Core principle: Raw measurements first (SPEC §18), layers 2-4 deferred.
    Every record carries provenance, version, and confidence so that:
    - Every score is reproducible from raw evidence + formula version.
    - Score changes across software versions are explainable.
    - Missing evidence is explicitly distinguished from zero evidence (SPEC §3).
    """
    # Identity
    evidence_id: str  # UUID
    analysis_id: str

    # Input and matched components
    input_component_id: str
    input_component_type: ComponentType
    matched_component_type: ComponentType
    component_relationship: ComponentRelationship

    # Sequence identifiers
    input_sequence_checksum: str  # sha256 hex
    matched_sequence_checksum: str

    # Source provenance (SPEC §13)
    source_database: str
    source_snapshot: SourceSnapshot
    source_accession: str

    # Raw alignment measurements (SPEC §18)
    alignment_method: str  # e.g., "blast", "exact_match"
    sequence_identity: float  # 0.0-1.0
    query_coverage: float  # 0.0-1.0
    subject_coverage: float  # 0.0-1.0
    alignment_length: int
    insertions: int
    deletions: int

    # Alignment metadata
    query_length: int
    subject_length: int

    # Evidence classification (SPEC §17 — never combined)
    evidence_class: EvidenceClass
    evidence_directness: str  # e.g., "direct_sequence_observation"

    # Quality and completeness
    metadata_completeness: float  # 0.0-1.0

    # Versioning (SPEC §19)
    formula_version: str
    calibration_version: str

    # Optional fields (must come after non-default fields)
    e_value: Optional[float] = None
    bit_score: Optional[float] = None
    warnings: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate constraint: component types determine eligibility."""
        eligibility = compare_eligibility(self.input_component_type, self.matched_component_type)
        if eligibility != self.component_relationship:
            raise ValueError(
                f"Component relationship mismatch: "
                f"input={self.input_component_type}, matched={self.matched_component_type}, "
                f"declared={self.component_relationship}, expected={eligibility}"
            )
