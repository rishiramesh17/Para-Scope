"""
Evidence strength, confidence, and completeness (Version 2, Gate 10).

Implements SPEC §28 (completeness must be split) and SPEC §19 Layers 3-4
(evidence support grade, sufficiency and uncertainty).

The one rule this module exists to enforce
------------------------------------------
SPEC §28: *"Do not use one ambiguous completeness score."* Four things get
conflated by a single number, and they have different causes, different owners,
and different remedies:

  28.1 **Input completeness**     did the submitter give us the components?
  28.2 **Analysis completion**    did our modules run?
  28.3 **Metadata completeness**  are the matched source records well annotated?
  28.4 **Evidence sufficiency**   is there enough independent evidence to interpret?

A candidate submitted without a light chain (28.1), a corpus snapshot that failed
to load (28.2), a matched PDB entry with no resolution recorded (28.3), and a
match resting on one record (28.4) are four different situations. Averaged into
"completeness: 0.62" they are indistinguishable, and the user cannot tell whether
to supply more input, retry, distrust the annotation, or distrust the conclusion.

So there is no aggregate completeness field here, and a test asserts none exists.

The distinction that matters most
--------------------------------
SPEC §28: *"A database outage affects analysis completion. It must not be treated
as biological evidence against the candidate."*

A failed module and a module that ran and found nothing produce the same *absence*
of evidence and must never produce the same grade. Infrastructure failure yields
`INSUFFICIENT` -- we cannot grade -- never `D` -- the evidence is weak. `D` is a
statement about a candidate; `INSUFFICIENT` is a statement about us. This is
`_INFRASTRUCTURAL_STATUSES` below, and it short-circuits grading entirely.

What a grade is not
-------------------
A grade describes **the evidence**, not the candidate. Grade A means "direct,
curated, same-component experimental evidence with strong provenance" exists --
not that the candidate is good, likely to work, or safe. A candidate with grade A
evidence of similarity to a discontinued programme has excellent evidence of
similarity to a discontinued programme (CLAUDE.md, prohibited claims).
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from para_scope.models import ComponentType, EvidenceClass, EvidenceStatus

# Bump on any change to how the four fields or the counts are computed.
#   1.0.0 -- initial: SPEC §28 four-field split, SPEC §19 Layer 4 counts
EVIDENCE_STRENGTH_VERSION = "1.0.0"

# Versioned separately from the engine, because SPEC §19 Layer 3 requires it:
# "The exact grade thresholds must be benchmarked and versioned." A threshold
# change alters every stored grade without altering any measurement, so the two
# versions must be independently traceable.
#
# PROVISIONAL -- see reports/GATE10_EVIDENCE_STRENGTH.md. These thresholds are
# chosen to match SPEC §19's prose definitions and then checked against the real
# corpus for a sensible spread. That is a weaker claim than "benchmarked": no
# external ground truth of evidence quality exists to benchmark against, so what
# was verified is that the bands separate real cases rather than collapsing onto
# one grade.
GRADE_THRESHOLD_VERSION = "1.1.0-provisional"


class EvidenceGrade(str, Enum):
    """SPEC §19 Layer 3. A statement about the evidence, never about the candidate."""

    A = "A"  # direct, curated, same-component experimental evidence, strong provenance
    B = "B"  # strong curated or independently repeated component precedent
    C = "C"  # valid sequence precedent, limited experimental/structure/pairing context
    D = "D"  # weak, indirect, sparse, or poorly annotated
    INSUFFICIENT = "Insufficient"  # not enough evidence to grade


# Statuses that mean *we* failed, not that the candidate lacks evidence. These can
# never produce a biological grade (SPEC §28's closing rule).
_INFRASTRUCTURAL_STATUSES = frozenset({
    EvidenceStatus.DATABASE_UNAVAILABLE,
    EvidenceStatus.ANALYSIS_FAILED,
    EvidenceStatus.ACCESS_RESTRICTED,
})

# Thresholds for each band. Named rather than inlined so a threshold change is a
# visible, reviewable diff and so GRADE_THRESHOLD_VERSION means something.
GRADE_A_MIN_PARENT_ANTIBODIES = 2      # independently repeated, not one record twice
GRADE_A_MIN_STRUCTURE_RECORDS = 1      # "experimental" requires an experiment
GRADE_A_MIN_IDENTITY = 0.95
GRADE_B_MIN_PARENT_ANTIBODIES = 2
GRADE_B_MIN_IDENTITY = 0.85
GRADE_C_MIN_IDENTITY = 0.70            # matches the admissibility floor
GRADE_MIN_RECORDS_TO_GRADE_AT_ALL = 1

# SPEC S19 defines D as "weak, indirect, predicted, sparse, or poorly annotated".
# An identity threshold alone cannot express that, and relying on one made D
# unreachable: 0 of 120 real candidates graded D, because GRADE_C_MIN_IDENTITY
# coincides with the admissibility coverage floor, so anything admitted at all
# already cleared C. Sparseness and annotation are the other two routes the spec
# names, and they are the ones that actually occur.
GRADE_D_MAX_PARENTS_FOR_SPARSE = 1     # resting on a single parent antibody
GRADE_D_MAX_METADATA_FRACTION = 0.25   # matched records barely annotated

GRADE_DISCLAIMER = (
    "An evidence grade describes the strength and provenance of the database "
    "evidence found, not the candidate. Grade A means direct, curated, "
    "same-component experimental evidence exists -- it is not a statement that "
    "the candidate will fold, express, pair, bind, be stable, be safe, be "
    "effective, or succeed."
)


@dataclass(frozen=True)
class DimensionReport:
    """
    One evidence dimension's contribution, in the shape Gate 10 consumes.

    Deliberately a flat, generic record rather than a union of each gate's own
    result type. Gate 10 would otherwise have to know the field names of every
    gate, and would break whenever one of them changed shape -- which is exactly
    what happened to the pair index when component checksums moved. Adapters
    below translate; this stays stable.
    """

    name: str
    status: EvidenceStatus
    evidence_class: Optional[EvidenceClass] = None
    reason: Optional[str] = None

    # SPEC §19 Layer 4 counts. Deduplicated counts are the caller's
    # responsibility -- this module adds them up, it cannot know that two records
    # describe one antibody (SPEC §15).
    independent_records: int = 0
    parent_antibodies: int = 0
    paired_records: int = 0
    structure_linked_records: int = 0
    therapeutic_records: int = 0

    # Metadata completeness (28.3), as present/expected rather than a ratio, so a
    # dimension with nothing to annotate is distinguishable from one whose records
    # are entirely unannotated.
    metadata_fields_present: int = 0
    metadata_fields_expected: int = 0

    # Identities of the parent antibodies this dimension matched, when known.
    #
    # SPEC S15 requires one antibody to be counted once. Summing
    # `parent_antibodies` across dimensions double-counts an antibody matched via
    # both therapeutic precedent and pairing -- the normal case, since a
    # therapeutic that has also been crystallised appears in both. When every
    # contributing dimension supplies identifiers the union is used; when any does
    # not, the sum is used and the inability to deduplicate is recorded.
    parent_identifiers: frozenset[str] = frozenset()

    best_identity: Optional[float] = None
    calibrated: bool = False
    # Terms this dimension could not compute, and why. Surfaced as missing
    # evidence (Layer 4) rather than silently absent.
    blocked_terms: dict[str, str] = field(default_factory=dict)

    @property
    def is_infrastructural_failure(self) -> bool:
        return self.status in _INFRASTRUCTURAL_STATUSES

    @property
    def produced_evidence(self) -> bool:
        return self.status is EvidenceStatus.AVAILABLE and self.independent_records > 0


@dataclass
class InputCompleteness:
    """SPEC §28.1 -- which expected components the submitter supplied."""

    expected: list[str] = field(default_factory=list)
    supplied: list[str] = field(default_factory=list)
    # Each missing component with why it is missing, because "not supplied" and
    # "not applicable to this format" are different (a VHH has no light chain by
    # design; reporting that as missing input would describe an architecture as a
    # deficiency).
    missing: dict[str, str] = field(default_factory=dict)
    not_applicable: dict[str, str] = field(default_factory=dict)

    @property
    def supplied_count(self) -> int:
        return len(self.supplied)

    @property
    def required_count(self) -> int:
        """Expected components excluding those that do not apply to this format."""
        return len([c for c in self.expected if c not in self.not_applicable])

    @property
    def is_complete(self) -> bool:
        return not self.missing


@dataclass
class AnalysisCompletion:
    """SPEC §28.2 -- which modules ran. Never read as evidence about the candidate."""

    attempted: list[str] = field(default_factory=list)
    completed: list[str] = field(default_factory=list)
    # Module -> (status, reason). Split from `completed` because a module that ran
    # and found nothing completed successfully.
    failed: dict[str, tuple[str, str]] = field(default_factory=dict)
    # Failures attributable to infrastructure rather than to the data. Separated
    # so a consumer can retry rather than conclude.
    infrastructural_failures: dict[str, str] = field(default_factory=dict)

    @property
    def is_complete(self) -> bool:
        return not self.failed

    @property
    def has_infrastructural_failure(self) -> bool:
        return bool(self.infrastructural_failures)


@dataclass
class MetadataCompleteness:
    """SPEC §28.3 -- how well annotated the matched source records are."""

    fields_present: int = 0
    fields_expected: int = 0
    # Per dimension, so a poorly annotated source is attributable rather than
    # averaged into an opaque figure.
    by_dimension: dict[str, tuple[int, int]] = field(default_factory=dict)

    @property
    def fraction(self) -> Optional[float]:
        """None when nothing was expected -- not 0.0, and not 1.0.

        A dimension with no annotatable records has no metadata completeness.
        Reporting 0.0 would say the records are unannotated; 1.0 would say they
        are perfectly annotated. Both are claims about records that do not exist.
        """
        if self.fields_expected == 0:
            return None
        return self.fields_present / self.fields_expected


@dataclass
class EvidenceSufficiency:
    """SPEC §28.4 and §19 Layer 4 -- is there enough independent evidence?"""

    independent_records: int = 0
    parent_antibodies: int = 0
    paired_records: int = 0
    structure_linked_records: int = 0
    therapeutic_records: int = 0
    dimensions_with_evidence: int = 0
    dimensions_attempted: int = 0
    # False means `parent_antibodies` is a sum across dimensions and may count one
    # antibody more than once. Reported rather than silently assumed either way.
    parent_antibodies_deduplicated: bool = False

    # Evidence classes are reported separately and never summed (SPEC §17).
    by_evidence_class: dict[str, int] = field(default_factory=dict)

    # Layer 4's "missing evidence": every term no dimension could compute.
    missing_evidence: dict[str, str] = field(default_factory=dict)

    # Layer 4 also lists study count, subject count, and a bootstrap interval.
    # None of the three is derivable from the current sources, and each is stated
    # as unavailable rather than silently omitted or filled with a count that
    # means something else.
    study_count: Optional[int] = None
    subject_count: Optional[int] = None
    bootstrap_interval: Optional[tuple[float, float]] = None
    unavailable_layer4_fields: dict[str, str] = field(
        default_factory=lambda: {
            "study_count": (
                "Neither SAbDab2 nor Thera-SAbDab records a study identifier "
                "distinct from the structure or therapeutic entry, so independent "
                "studies cannot be counted without conflating them with records."
            ),
            "subject_count": (
                "Subject-level provenance exists only in repertoire sources (OAS). "
                "promb's bundled reference is pre-aggregated at a 10%-of-subjects "
                "threshold and exposes no per-subject counts."
            ),
            "bootstrap_interval": (
                "Requires resampling a reference distribution per query. The "
                "calibration artifact stores full observations so this is "
                "computable, but it is not computed here and an interval has not "
                "been validated."
            ),
        }
    )


@dataclass
class EvidenceStrengthReport:
    """
    The four SPEC §28 fields, plus the graded layer, deliberately never combined.

    There is no `overall_completeness`. A consumer that wants one must decide for
    itself what to trade against what, in the open, rather than inherit a weighting
    from this module.
    """

    input_completeness: InputCompleteness
    analysis_completion: AnalysisCompletion
    metadata_completeness: MetadataCompleteness
    evidence_sufficiency: EvidenceSufficiency

    grade: EvidenceGrade
    grade_rationale: str
    # Set when grading was refused because a module failed infrastructurally.
    grade_blocked_reason: Optional[str] = None

    engine_version: str = EVIDENCE_STRENGTH_VERSION
    grade_threshold_version: str = GRADE_THRESHOLD_VERSION
    disclaimer: str = GRADE_DISCLAIMER
    limitations: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------
# Adapters: each gate's own result type -> DimensionReport
# --------------------------------------------------------------------------

# Display fields SPEC §26 expects on a therapeutic match. Counted for 28.3.
_THERAPEUTIC_METADATA_FIELDS = (
    "target", "format_raw", "ch1_isotype", "highest_clinical_trial",
    "development_status",
)


def from_therapeutic_precedent(result, name: str = "therapeutic_precedent") -> DimensionReport:
    """Adapt a TherapeuticPrecedentResult."""
    present = 0
    expected = 0
    structure_linked = 0
    for match in result.matches:
        for attribute in _THERAPEUTIC_METADATA_FIELDS:
            expected += 1
            if getattr(match, attribute, None):
                present += 1
        if getattr(match, "has_structure_record", None) or getattr(match, "pdb_ids", None):
            structure_linked += 1

    return DimensionReport(
        name=name,
        status=result.status,
        evidence_class=result.evidence_class,
        reason=result.reason,
        independent_records=result.raw_match_count,
        # One therapeutic may contribute several matched records; the therapeutic
        # is the parent entity (SPEC §15).
        parent_antibodies=result.therapeutics_matched,
        parent_identifiers=frozenset(m.therapeutic_name for m in result.matches),
        therapeutic_records=result.therapeutics_matched,
        structure_linked_records=structure_linked,
        metadata_fields_present=present,
        metadata_fields_expected=expected,
        best_identity=result.best_global_identity,
        calibrated=any(
            getattr(m, "calibrated_similarity", None) is not None
            and m.calibrated_similarity.status is EvidenceStatus.AVAILABLE
            for m in result.matches
        ),
        blocked_terms=dict(result.blocked_terms),
    )


def from_pairing_precedent(result, name: str = "pairing_precedent") -> DimensionReport:
    """Adapt a PairingPrecedentResult."""
    structure_linked = sum(1 for m in result.matches if m.structure_supported is True)
    therapeutic_records = len({
        therapeutic
        for match in result.matches
        for therapeutic in match.therapeutic_names
    })
    present = 0
    expected = 0
    for match in result.matches:
        for value in (match.targets, match.development_statuses, match.structure_ids):
            expected += 1
            if value:
                present += 1

    return DimensionReport(
        name=name,
        status=result.status,
        evidence_class=result.evidence_class,
        reason=result.reason,
        independent_records=result.raw_observation_count,
        parent_antibodies=result.pairs_matched,
        # A pair's parents are the therapeutics documenting it where named, and
        # otherwise the pair's own record ids -- a structure-only pair is still a
        # distinct parent antibody even with no drug name.
        parent_identifiers=frozenset(
            name
            for match in result.matches
            for name in (match.therapeutic_names or [f"pair:{match.heavy_record_id}"])
        ),
        paired_records=result.raw_observation_count,
        structure_linked_records=structure_linked,
        therapeutic_records=therapeutic_records,
        metadata_fields_present=present,
        metadata_fields_expected=expected,
        best_identity=(
            min(result.best_heavy_identity, result.best_light_identity)
            if result.best_heavy_identity is not None
            and result.best_light_identity is not None
            else None
        ),
        blocked_terms=dict(result.blocked_dimensions),
    )


def from_germline(result, name: str = "germline") -> DimensionReport:
    """Adapt a GermlineAssignmentResult.

    Germline is a reference-set lookup, not database precedent: a hit is one
    alignment against a curated allele, so it contributes no parent antibodies
    and no independent records in the Layer 4 sense. Counting the allele as an
    "independent record" would let a single germline call inflate sufficiency.
    """
    return DimensionReport(
        name=name,
        status=result.status,
        evidence_class=result.evidence_class,
        reason=result.reason,
        independent_records=0,
        best_identity=result.nearest_germline_identity,
        calibrated=result.calibration_version not in (None, "", "uncalibrated"),
    )


def from_window_coverage(result, name: str = "window_coverage") -> DimensionReport:
    """Adapt a WindowCoverageResult. Same reasoning as germline: a reference lookup."""
    return DimensionReport(
        name=name,
        status=result.status,
        evidence_class=EvidenceClass.OBSERVED_REPERTOIRE,
        reason=result.reason,
        independent_records=0,
        best_identity=result.coverage,
    )


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------


def _grade(
    dimensions: list[DimensionReport],
    sufficiency: EvidenceSufficiency,
    metadata: "MetadataCompleteness",
):
    """Assign SPEC §19 Layer 3's grade, or refuse to.

    Returns (grade, rationale, blocked_reason).
    """
    infrastructural = [d for d in dimensions if d.is_infrastructural_failure]
    if infrastructural:
        names = ", ".join(f"{d.name} ({d.status.value})" for d in infrastructural)
        return (
            EvidenceGrade.INSUFFICIENT,
            "Grading was not attempted.",
            (
                f"A module failed for infrastructural reasons: {names}. A failed "
                f"lookup and a lookup that found nothing produce the same absence "
                f"of evidence and must not produce the same grade -- so this is "
                f"Insufficient (we cannot grade), never D (the evidence is weak). "
                f"Retry rather than interpret."
            ),
        )

    if sufficiency.independent_records < GRADE_MIN_RECORDS_TO_GRADE_AT_ALL:
        return (
            EvidenceGrade.INSUFFICIENT,
            (
                f"No dimension returned an evidence record "
                f"({sufficiency.dimensions_with_evidence} of "
                f"{sufficiency.dimensions_attempted} dimensions produced evidence)."
            ),
            None,
        )

    best = max(
        (d.best_identity for d in dimensions if d.best_identity is not None),
        default=None,
    )
    if best is None:
        return (
            EvidenceGrade.INSUFFICIENT,
            "Evidence records exist but none carries a similarity measurement to grade.",
            None,
        )

    repeated = sufficiency.parent_antibodies >= GRADE_A_MIN_PARENT_ANTIBODIES
    experimental = sufficiency.structure_linked_records >= GRADE_A_MIN_STRUCTURE_RECORDS

    if best >= GRADE_A_MIN_IDENTITY and repeated and experimental:
        return (
            EvidenceGrade.A,
            (
                f"Identity {best:.3f} to curated records, repeated across "
                f"{sufficiency.parent_antibodies} parent antibodies, with "
                f"{sufficiency.structure_linked_records} structure-linked record(s)."
            ),
            None,
        )
    if best >= GRADE_B_MIN_IDENTITY and repeated:
        return (
            EvidenceGrade.B,
            (
                f"Identity {best:.3f} repeated across "
                f"{sufficiency.parent_antibodies} parent antibodies, but "
                f"{'no' if not experimental else 'limited'} structure-linked "
                f"experimental context."
            ),
            None,
        )
    # D's other two routes, from SPEC S19's own wording: sparse, or poorly
    # annotated. Checked before C, because evidence resting on one barely
    # annotated record is weak regardless of how well it aligns.
    sparse = sufficiency.parent_antibodies <= GRADE_D_MAX_PARENTS_FOR_SPARSE
    fraction = metadata.fraction
    poorly_annotated = fraction is not None and fraction <= GRADE_D_MAX_METADATA_FRACTION

    if best >= GRADE_C_MIN_IDENTITY and (sparse and poorly_annotated):
        return (
            EvidenceGrade.D,
            (
                f"Sparse and poorly annotated: identity {best:.3f} rests on "
                f"{sufficiency.parent_antibodies} parent antibody with "
                f"{fraction:.0%} of expected metadata fields populated. The "
                f"alignment is valid; the evidence behind it is thin."
            ),
            None,
        )
    if best >= GRADE_C_MIN_IDENTITY:
        return (
            EvidenceGrade.C,
            (
                f"Valid sequence precedent at identity {best:.3f}, resting on "
                f"{sufficiency.parent_antibodies} parent antibody(s) with limited "
                f"experimental, structural, or pairing context."
            ),
            None,
        )
    return (
        EvidenceGrade.D,
        (
            f"Weak precedent: best identity {best:.3f} is below the "
            f"{GRADE_C_MIN_IDENTITY} threshold for valid sequence precedent."
        ),
        None,
    )


def assess_evidence(
    dimensions: list[DimensionReport],
    expected_components: Optional[dict[str, Optional[str]]] = None,
    not_applicable_components: Optional[dict[str, str]] = None,
) -> EvidenceStrengthReport:
    """
    Assemble the four SPEC §28 completeness fields and the SPEC §19 Layer 3 grade.

    `expected_components` maps a component name to the sequence supplied for it, or
    None when it was not supplied. `not_applicable_components` maps a component to
    why it does not apply to this construct format -- a VHH has no light chain by
    design, and reporting that as missing input would describe an architecture as a
    deficiency.
    """
    expected_components = expected_components or {}
    not_applicable_components = not_applicable_components or {}

    # --- 28.1 input completeness ---
    supplied = [name for name, sequence in expected_components.items() if sequence]
    missing = {
        name: "not supplied"
        for name, sequence in expected_components.items()
        if not sequence and name not in not_applicable_components
    }
    input_completeness = InputCompleteness(
        expected=list(expected_components),
        supplied=supplied,
        missing=missing,
        not_applicable=dict(not_applicable_components),
    )

    # --- 28.2 analysis completion ---
    completion = AnalysisCompletion(attempted=[d.name for d in dimensions])
    for dimension in dimensions:
        if dimension.is_infrastructural_failure:
            completion.infrastructural_failures[dimension.name] = (
                dimension.reason or dimension.status.value
            )
            completion.failed[dimension.name] = (
                dimension.status.value, dimension.reason or ""
            )
        elif dimension.status is EvidenceStatus.ANALYSIS_FAILED:
            completion.failed[dimension.name] = (
                dimension.status.value, dimension.reason or ""
            )
        else:
            # A module that ran and found nothing completed successfully. Counting
            # NOT_FOUND as a failure would turn an answer into an outage.
            completion.completed.append(dimension.name)

    # --- 28.3 metadata completeness ---
    metadata = MetadataCompleteness()
    for dimension in dimensions:
        if dimension.metadata_fields_expected:
            metadata.by_dimension[dimension.name] = (
                dimension.metadata_fields_present, dimension.metadata_fields_expected
            )
            metadata.fields_present += dimension.metadata_fields_present
            metadata.fields_expected += dimension.metadata_fields_expected

    # --- 28.4 evidence sufficiency ---
    sufficiency = EvidenceSufficiency(
        dimensions_attempted=len(dimensions),
        dimensions_with_evidence=sum(1 for d in dimensions if d.produced_evidence),
    )
    # Deduplicate parent antibodies across dimensions when possible (SPEC S15).
    contributing = [d for d in dimensions if d.parent_antibodies]
    deduplicated = bool(contributing) and all(d.parent_identifiers for d in contributing)
    if deduplicated:
        sufficiency.parent_antibodies = len(
            frozenset().union(*(d.parent_identifiers for d in contributing))
        )
        sufficiency.parent_antibodies_deduplicated = True

    for dimension in dimensions:
        sufficiency.independent_records += dimension.independent_records
        if not deduplicated:
            sufficiency.parent_antibodies += dimension.parent_antibodies
        sufficiency.paired_records += dimension.paired_records
        sufficiency.structure_linked_records += dimension.structure_linked_records
        sufficiency.therapeutic_records += dimension.therapeutic_records
        if dimension.evidence_class is not None and dimension.produced_evidence:
            key = dimension.evidence_class.value
            sufficiency.by_evidence_class[key] = (
                sufficiency.by_evidence_class.get(key, 0) + dimension.independent_records
            )
        for term, why in dimension.blocked_terms.items():
            sufficiency.missing_evidence[f"{dimension.name}.{term}"] = why

    grade, rationale, blocked = _grade(dimensions, sufficiency, metadata)

    report = EvidenceStrengthReport(
        input_completeness=input_completeness,
        analysis_completion=completion,
        metadata_completeness=metadata,
        evidence_sufficiency=sufficiency,
        grade=grade,
        grade_rationale=rationale,
        grade_blocked_reason=blocked,
    )

    if contributing and not deduplicated:
        report.limitations.append(
            "Parent-antibody count is a sum across dimensions and may count one "
            "antibody more than once: at least one contributing dimension supplied "
            "no parent identifiers, so cross-dimension deduplication was not "
            "possible (SPEC S15 wants one antibody counted once)."
        )

    if metadata.fraction is None:
        report.limitations.append(
            "No matched record carried annotatable metadata, so metadata "
            "completeness is undefined rather than zero."
        )
    if not input_completeness.is_complete:
        report.limitations.append(
            f"Input incomplete: {sorted(input_completeness.missing)} not supplied. "
            f"This limits what could be analysed and is not evidence about the "
            f"candidate."
        )
    report.limitations.append(
        f"Grade thresholds are provisional (version {GRADE_THRESHOLD_VERSION}). "
        f"SPEC §19 requires them to be benchmarked; no external ground truth of "
        f"evidence quality exists to benchmark against, so what has been verified "
        f"is that the bands separate real corpus cases rather than collapsing onto "
        f"one grade."
    )
    return report
