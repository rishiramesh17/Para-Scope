"""
Whole-construct profile assembly and the comparative evidence index (Gate 11, SPEC §37).

The component profile is authoritative; the index is a summary of it
-------------------------------------------------------------------
SPEC §37 opens with "The component profile is authoritative." That ordering is the
whole design. The index is one number derived from per-component evidence, and it
exists to be skimmed, not to replace what it summarises. Every field it carries is
there so a reader can get back to the components.

The arithmetic that makes silence look like completeness
--------------------------------------------------------
§37's index is a reliability-weighted geometric mean:

    100 × exp( Σ[wᵢ × rᵢ × ln(max(sᵢ/100, ε))] / Σ[wᵢ × rᵢ] )

Note the denominator. It sums only over components that *contributed*. Drop a
component -- because its data does not exist, because a module failed, because the
corpus has no constant regions -- and the divisor shrinks with the numerator. The
index self-normalises and comes out looking exactly as confident as a complete
one. Nothing in the formula notices.

So coverage reporting is not a nicety here; it is the only thing standing between
this index and a number that overstates itself every time. `evidence_coverage` is
computed against the **full weight set**, not against what contributed:

    evidence_coverage = Σ[wᵢ for contributing i] / Σ[wᵢ for all i in the weight set]

An index resting on 0.80 of its own weight says so. CLAUDE.md states the rule
directly: never renormalise silently.

Scope is declared, and the name is derived from it
--------------------------------------------------
`coverage_scope` is `VARIABLE_DOMAIN_ONLY` or `FULL_CONSTRUCT`, and the displayed
name follows from it rather than being chosen. Today the corpus has constant-region
components but no full chains, so an index renders as "Comparative Evidence Index
(variable domain only)". When constant-region coverage improves the qualifier
disappears on its own -- no rename, no second name for one concept across reports
and git history.

This is deliberately *not* a rename of the index. The term "Whole-Construct
Comparative Evidence Index" names the spec's concept and the platform's target,
not whichever subset of components happens to be populated this quarter
(CLAUDE.md, "The index must state its own scope").

"Critical" is per format, not absolute
--------------------------------------
Constant/Fc is not required to claim variable-domain scope, and is required to
claim full-construct scope. A VHH has no light chain and no VL weight at all;
asking it for one would describe an architecture as a deficiency. So
`REQUIRED_COMPONENTS` is keyed by (format class, scope), and the index may compute
today while only ever claiming the scope it actually covered -- which satisfies
§37's "do not calculate the index when critical evidence is missing" without
declaring Fc non-critical by fiat.
"""

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from para_scope.construct_format import (
    FORMAT_VOCABULARY_VERSION,
    Architecture,
    ConstructFormat,
    Specificity,
)
from para_scope.evidence_strength import (
    EvidenceGrade,
    EvidenceStrengthReport,
)
from para_scope.models import EvidenceStatus

# Bump on any change to assembly, gating, or the derived fields.
#   1.0.0 -- initial: SPEC §37 index, gating, coverage, bottleneck
PROFILE_VERSION = "1.0.0"

# Versioned separately: §37's weights are explicitly "provisional and
# format-specific", and a weight change alters every index without altering any
# measurement. The suffix names the format the set was written for, because using
# the IgG-like set for a VHH is the mistake §37 warns about by name.
WEIGHT_SET_VERSION = "1.0.0-provisional-igg-like"

# Floor inside the logarithm, from §37's own `max(score/100, ε)`. Without it a
# single zero-scoring component sends the whole index to zero via ln(0); with it a
# zero component is very bad rather than annihilating.
SCORE_EPSILON = 0.01

# Verbatim, and asserted verbatim by test. Required on every emitted index by both
# SPEC §37 and CLAUDE.md.
INDEX_DISCLAIMER = (
    "This index summarizes comparative database evidence. It is not a probability "
    "of folding, expression, pairing, binding, stability, manufacturability, "
    "safety, efficacy, or experimental success."
)

# Exact string §37 requires when the index may not be calculated.
INSUFFICIENT_DISPLAY = "Insufficient evidence for a whole-construct index."


class ComponentSlot(str, Enum):
    """The weighted slots of §37's table. Coarser than ComponentType by design:
    §37 weights CDR-H1 and H2 together and all light CDRs together, so a slot is
    not always one component."""

    VH = "VH"
    VL = "VL"
    CDR_H3 = "CDR-H3"
    CDR_H1_H2 = "CDR-H1 and H2 combined"
    CDR_L1_L2_L3 = "CDR-L1, L2, and L3 combined"
    HEAVY_FRAMEWORKS = "Heavy frameworks"
    LIGHT_FRAMEWORKS = "Light frameworks"
    OBSERVED_PAIRING = "Observed heavy-light pairing"
    STRUCTURE_PRECEDENT = "Structure-linked experimental precedent"
    CONSTANT_FC = "Constant/Fc region"


class CoverageScope(str, Enum):
    VARIABLE_DOMAIN_ONLY = "VARIABLE_DOMAIN_ONLY"
    FULL_CONSTRUCT = "FULL_CONSTRUCT"


class FormatClass(str, Enum):
    """Which weight set applies. Not the same as construct_format's axes: several
    formats share a weight set, and some have none written yet."""

    IGG_LIKE = "igg_like"           # paired Fv or IgG-like, §37's table
    SINGLE_DOMAIN = "single_domain"  # VHH: no light chain at all
    UNSUPPORTED = "unsupported"      # scFv, bispecific, fusion: no set written


# SPEC §37's provisional table, verbatim. Sums to 1.00.
_IGG_LIKE_WEIGHTS: dict[ComponentSlot, float] = {
    ComponentSlot.VH: 0.14,
    ComponentSlot.VL: 0.14,
    ComponentSlot.CDR_H3: 0.14,
    ComponentSlot.CDR_H1_H2: 0.10,
    ComponentSlot.CDR_L1_L2_L3: 0.12,
    ComponentSlot.HEAVY_FRAMEWORKS: 0.08,
    ComponentSlot.LIGHT_FRAMEWORKS: 0.08,
    ComponentSlot.OBSERVED_PAIRING: 0.10,
    ComponentSlot.STRUCTURE_PRECEDENT: 0.05,
    ComponentSlot.CONSTANT_FC: 0.05,
}

# A VHH has no light chain and no pairing, by architecture rather than by missing
# data. §37 says explicitly not to reuse its weights for VHH without a separate
# versioned set, so this is one -- the heavy-side weights renormalised over what a
# single-domain construct actually has. It is provisional in exactly the way §37's
# own table is, and is NOT a claim that these proportions are calibrated.
_SINGLE_DOMAIN_WEIGHTS: dict[ComponentSlot, float] = {
    ComponentSlot.VH: 0.26,
    ComponentSlot.CDR_H3: 0.26,
    ComponentSlot.CDR_H1_H2: 0.18,
    ComponentSlot.HEAVY_FRAMEWORKS: 0.15,
    ComponentSlot.STRUCTURE_PRECEDENT: 0.08,
    ComponentSlot.CONSTANT_FC: 0.07,
}

WEIGHT_SETS: dict[FormatClass, dict[ComponentSlot, float]] = {
    FormatClass.IGG_LIKE: _IGG_LIKE_WEIGHTS,
    FormatClass.SINGLE_DOMAIN: _SINGLE_DOMAIN_WEIGHTS,
}

# Slots that must carry a valid score before the index may be calculated, per
# (format class, claimed scope). Constant/Fc appears only under FULL_CONSTRUCT --
# that is the entire difference between the two scopes.
REQUIRED_COMPONENTS: dict[tuple[FormatClass, CoverageScope], frozenset[ComponentSlot]] = {
    (FormatClass.IGG_LIKE, CoverageScope.VARIABLE_DOMAIN_ONLY): frozenset({
        ComponentSlot.VH, ComponentSlot.VL, ComponentSlot.CDR_H3,
    }),
    (FormatClass.IGG_LIKE, CoverageScope.FULL_CONSTRUCT): frozenset({
        ComponentSlot.VH, ComponentSlot.VL, ComponentSlot.CDR_H3,
        ComponentSlot.CONSTANT_FC,
    }),
    (FormatClass.SINGLE_DOMAIN, CoverageScope.VARIABLE_DOMAIN_ONLY): frozenset({
        ComponentSlot.VH, ComponentSlot.CDR_H3,
    }),
    (FormatClass.SINGLE_DOMAIN, CoverageScope.FULL_CONSTRUCT): frozenset({
        ComponentSlot.VH, ComponentSlot.CDR_H3, ComponentSlot.CONSTANT_FC,
    }),
}

# Reliability per evidence grade. §37 weights each component by a reliability but
# does not define one, so this does: a component graded on sparse or
# poorly-annotated evidence contributes less than one graded on repeated
# experimental records. Versioned with the weight set because it multiplies into
# the same product.
RELIABILITY_BY_GRADE: dict[EvidenceGrade, float] = {
    EvidenceGrade.A: 1.00,
    EvidenceGrade.B: 0.80,
    EvidenceGrade.C: 0.60,
    EvidenceGrade.D: 0.30,
    # Zero, so an ungradeable component cannot move the index at all -- it is
    # excluded from both sums and shows up as missing coverage instead.
    EvidenceGrade.INSUFFICIENT: 0.00,
}


def classify_format(construct_format: Optional[ConstructFormat]) -> FormatClass:
    """Which weight set a construct format maps to.

    Returns UNSUPPORTED rather than falling back to the IgG-like set for scFv,
    bispecific and fusion formats. §37 names that reuse as the mistake not to
    make, and an index computed under the wrong weight set is wrong in a way no
    coverage field reveals.
    """
    if construct_format is None:
        return FormatClass.UNSUPPORTED
    if construct_format.architecture is Architecture.SINGLE_DOMAIN:
        return FormatClass.SINGLE_DOMAIN
    if construct_format.specificity not in (
        Specificity.MONOSPECIFIC, Specificity.UNKNOWN
    ):
        return FormatClass.UNSUPPORTED
    if construct_format.architecture in (
        Architecture.WHOLE_MAB, Architecture.FAB, Architecture.MAB_UNSPECIFIED
    ):
        return FormatClass.IGG_LIKE
    return FormatClass.UNSUPPORTED


@dataclass
class ComponentEvidence:
    """One weighted slot's contribution to the index."""

    slot: ComponentSlot
    status: EvidenceStatus
    # 0-100, on the same scale as §26's provisional component score. None whenever
    # status is not AVAILABLE -- never 0.0 as a stand-in, since zero is a score and
    # absence is not.
    score: Optional[float] = None
    grade: EvidenceGrade = EvidenceGrade.INSUFFICIENT
    reason: Optional[str] = None

    @property
    def reliability(self) -> float:
        return RELIABILITY_BY_GRADE[self.grade]

    @property
    def contributes(self) -> bool:
        """Whether this slot enters either sum of the geometric mean."""
        return (
            self.status is EvidenceStatus.AVAILABLE
            and self.score is not None
            and self.reliability > 0.0
        )


@dataclass
class WholeConstructIndex:
    """SPEC §37's index plus every field required to read it honestly."""

    status: EvidenceStatus
    # None whenever the index may not be calculated. The display string is what
    # §37 requires be shown in that case.
    index: Optional[float] = None
    display: str = INSUFFICIENT_DISPLAY

    # --- scope, and the name derived from it ---
    coverage_scope: CoverageScope = CoverageScope.VARIABLE_DOMAIN_ONLY
    displayed_name: str = ""
    components_present: list[str] = field(default_factory=list)
    # Slot -> why it is absent. Status and reason, never a bare omission.
    components_missing: dict[str, str] = field(default_factory=dict)
    weight_set_version: str = WEIGHT_SET_VERSION
    format_class: FormatClass = FormatClass.UNSUPPORTED

    # --- §37's six additional required outputs ---
    bottleneck_slot: Optional[str] = None
    bottleneck_score: Optional[float] = None
    evidence_coverage: float = 0.0
    evidence_confidence: Optional[str] = None
    observed_pairing_status: str = EvidenceStatus.NOT_APPLICABLE.value
    assembly_validation_status: str = "not_validated"
    missing_critical_evidence: dict[str, str] = field(default_factory=dict)

    # --- provenance ---
    profile_version: str = PROFILE_VERSION
    vocabulary_version: str = FORMAT_VOCABULARY_VERSION
    disclaimer: str = INDEX_DISCLAIMER
    gating_failures: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)


def _derive_scope(
    contributions: dict[ComponentSlot, ComponentEvidence],
) -> CoverageScope:
    """Full-construct scope requires constant-region evidence; nothing else does."""
    constant = contributions.get(ComponentSlot.CONSTANT_FC)
    if constant is not None and constant.contributes:
        return CoverageScope.FULL_CONSTRUCT
    return CoverageScope.VARIABLE_DOMAIN_ONLY


def _displayed_name(scope: CoverageScope) -> str:
    """The name follows from the scope rather than being chosen.

    The concept keeps its name; the qualifier appears and disappears with the data.
    """
    if scope is CoverageScope.FULL_CONSTRUCT:
        return "Whole-Construct Comparative Evidence Index"
    return "Comparative Evidence Index (variable domain only)"


def assemble_index(
    components: list[ComponentEvidence],
    construct_format: Optional[ConstructFormat] = None,
    strength: Optional[EvidenceStrengthReport] = None,
    observed_pairing_status: EvidenceStatus = EvidenceStatus.NOT_APPLICABLE,
    assembly_validated: Optional[bool] = None,
) -> WholeConstructIndex:
    """
    Assemble SPEC §37's whole-construct comparative evidence index.

    Returns a result in every case. When §37's gating conditions are not met the
    index is None, `display` is §37's required string, and `gating_failures` says
    which conditions failed -- an index is withheld with reasons, never replaced by
    a lower number.
    """
    by_slot = {c.slot: c for c in components}
    format_class = classify_format(construct_format)

    result = WholeConstructIndex(
        status=EvidenceStatus.INSUFFICIENT_EVIDENCE,
        format_class=format_class,
        observed_pairing_status=observed_pairing_status.value,
        assembly_validation_status=(
            "not_validated" if assembly_validated is None
            else ("passed" if assembly_validated else "failed")
        ),
    )

    # --- no weight set, no index. §37 forbids reusing the IgG-like weights. ---
    if format_class is FormatClass.UNSUPPORTED:
        described = construct_format.raw if construct_format else "undeclared"
        result.gating_failures.append(
            f"No versioned weight set exists for this construct format "
            f"({described!r}). SPEC §37's provisional weights are format-specific "
            f"and must not be reused for VHH, scFv, bispecific or fusion formats; "
            f"an index computed under the wrong weight set is wrong in a way no "
            f"coverage field would reveal."
        )
        result.coverage_scope = _derive_scope(by_slot)
        result.displayed_name = _displayed_name(result.coverage_scope)
        return result

    weights = WEIGHT_SETS[format_class]
    scope = _derive_scope(by_slot)
    result.coverage_scope = scope
    result.displayed_name = _displayed_name(scope)

    # --- present / missing, against the weight set rather than against input ---
    for slot in weights:
        evidence = by_slot.get(slot)
        if evidence is not None and evidence.contributes:
            result.components_present.append(slot.value)
        elif evidence is None:
            result.components_missing[slot.value] = (
                "no evidence was supplied for this slot"
            )
        else:
            result.components_missing[slot.value] = (
                f"{evidence.status.value}"
                + (f": {evidence.reason}" if evidence.reason else "")
                + (
                    f" (grade {evidence.grade.value} gives reliability 0.0, so it "
                    f"cannot contribute)"
                    if evidence.status is EvidenceStatus.AVAILABLE
                    and evidence.reliability == 0.0
                    else ""
                )
            )

    # --- evidence coverage, against the FULL weight set ---
    #
    # This is the field that makes the geometric mean's self-normalisation
    # visible. Computing it over contributing weight only would always yield 1.0.
    total_weight = sum(weights.values())
    present_weight = sum(
        weight for slot, weight in weights.items()
        if (ev := by_slot.get(slot)) is not None and ev.contributes
    )
    result.evidence_coverage = present_weight / total_weight if total_weight else 0.0

    # --- §37 gating ---
    required = REQUIRED_COMPONENTS[(format_class, scope)]
    for slot in sorted(required, key=lambda s: s.value):
        evidence = by_slot.get(slot)
        if evidence is None or not evidence.contributes:
            reason = (
                "not supplied" if evidence is None
                else f"{evidence.status.value} (grade {evidence.grade.value})"
            )
            result.missing_critical_evidence[slot.value] = reason
    if result.missing_critical_evidence:
        result.gating_failures.append(
            f"Required components for {format_class.value} at "
            f"{scope.value} scope are missing or unscored: "
            f"{sorted(result.missing_critical_evidence)}."
        )

    if strength is not None:
        result.evidence_confidence = strength.grade.value
        if strength.analysis_completion.has_infrastructural_failure:
            result.gating_failures.append(
                f"A module failed for infrastructural reasons "
                f"({sorted(strength.analysis_completion.infrastructural_failures)}); "
                f"an index computed over a partial analysis would read as a verdict "
                f"on the candidate rather than an outage. Retry."
            )
        if strength.grade is EvidenceGrade.INSUFFICIENT:
            result.gating_failures.append(
                "Minimum evidence sufficiency is not met: the evidence strength "
                "assessment graded Insufficient."
            )

    if assembly_validated is False:
        result.gating_failures.append(
            "Assembly validation failed; §37 requires the arrangement to pass "
            "basic architecture validation before an index is calculated."
        )

    contributing = [ev for ev in by_slot.values() if ev.contributes and ev.slot in weights]
    if not contributing:
        result.gating_failures.append("No component carries a valid, gradeable score.")

    if result.gating_failures:
        result.display = INSUFFICIENT_DISPLAY
        return result

    # --- the index itself ---
    numerator = 0.0
    denominator = 0.0
    for evidence in contributing:
        weight = weights[evidence.slot]
        product = weight * evidence.reliability
        numerator += product * math.log(max(evidence.score / 100.0, SCORE_EPSILON))
        denominator += product

    result.index = 100.0 * math.exp(numerator / denominator)
    result.status = EvidenceStatus.AVAILABLE

    bottleneck = min(contributing, key=lambda ev: ev.score)
    result.bottleneck_slot = bottleneck.slot.value
    result.bottleneck_score = bottleneck.score

    result.display = (
        f"{result.displayed_name}: {result.index:.1f} "
        f"(evidence coverage {result.evidence_coverage:.0%}, "
        f"weight set {result.weight_set_version})"
    )

    if result.evidence_coverage < 1.0:
        result.limitations.append(
            f"The index rests on {result.evidence_coverage:.0%} of its weight set. "
            f"The geometric mean divides by the weight that contributed, so an "
            f"incomplete index is numerically indistinguishable from a complete "
            f"one -- this coverage figure is the only thing that distinguishes "
            f"them. Absent slots: {sorted(result.components_missing)}."
        )
    if scope is CoverageScope.VARIABLE_DOMAIN_ONLY:
        result.limitations.append(
            "Variable-domain scope: no constant-region evidence contributed, so "
            "this index says nothing about the constant region or the construct as "
            "a whole. The name carries the qualifier for that reason."
        )
    return result


# --------------------------------------------------------------------------
# Orchestration: a submitted candidate -> a whole-construct profile
# --------------------------------------------------------------------------

# Which IMGT component each weighted slot is scored from. A slot that combines
# several components (§37 weights CDR-H1 and H2 together) is scored from the mean
# of its parts, so one absent part lowers the slot rather than voiding it.
_SLOT_COMPONENTS: dict[ComponentSlot, tuple[str, ...]] = {
    ComponentSlot.VH: ("VH",),
    ComponentSlot.VL: ("VL-kappa", "VL-lambda"),
    ComponentSlot.CDR_H3: ("CDR-H3",),
    ComponentSlot.CDR_H1_H2: ("CDR-H1", "CDR-H2"),
    ComponentSlot.CDR_L1_L2_L3: ("CDR-L1", "CDR-L2", "CDR-L3"),
    ComponentSlot.HEAVY_FRAMEWORKS: ("FR-H1", "FR-H2", "FR-H3", "FR-H4"),
    ComponentSlot.LIGHT_FRAMEWORKS: ("FR-L1", "FR-L2", "FR-L3", "FR-L4"),
    ComponentSlot.CONSTANT_FC: (
        "Constant region (heavy, undecomposed)",
        "Constant region (light, undecomposed)",
    ),
}


@dataclass
class ComponentProfileEntry:
    """One component's own evidence, which §37 calls authoritative."""

    component_type: str
    slot: Optional[str]
    status: str
    score: Optional[float]
    score_reason: str
    grade: str
    best_identity: Optional[float]
    therapeutics_matched: int
    blocked_terms: dict[str, str] = field(default_factory=dict)
    limitations: list[str] = field(default_factory=list)


@dataclass
class WholeConstructProfile:
    """
    The full profile: per-component evidence plus the index that summarises it.

    `components` comes first because SPEC §37 opens by saying the component profile
    is authoritative. The index is a derived convenience, and a reader who
    disbelieves it must be able to go straight to what produced it.
    """

    components: list[ComponentProfileEntry] = field(default_factory=list)
    index: Optional[WholeConstructIndex] = None
    strength: Optional[EvidenceStrengthReport] = None
    declared_format: Optional[str] = None
    profile_version: str = PROFILE_VERSION
    disclaimer: str = INDEX_DISCLAIMER
    limitations: list[str] = field(default_factory=list)


def slot_for_component(component_type_value: str) -> Optional[ComponentSlot]:
    """The weighted slot a component contributes to, or None if §37 weights none."""
    for slot, members in _SLOT_COMPONENTS.items():
        if component_type_value in members:
            return slot
    return None


def build_component_evidence(
    entries: list[ComponentProfileEntry],
    pairing_score: Optional[float] = None,
    pairing_grade: EvidenceGrade = EvidenceGrade.INSUFFICIENT,
    structure_score: Optional[float] = None,
    structure_grade: EvidenceGrade = EvidenceGrade.INSUFFICIENT,
) -> list[ComponentEvidence]:
    """
    Collapse per-component entries onto §37's weighted slots.

    A combined slot averages its contributing components. Averaging is safe here
    and would not be safe for the index itself: these are the same measurement on
    the same scale for sibling components of one region, whereas the index's slots
    measure different things and carry different weights for that reason.

    Pairing and structure-linked precedent arrive separately because they are
    properties of the construct, not of any one component.
    """
    scored: dict[ComponentSlot, list[tuple[float, EvidenceGrade]]] = {}
    for entry in entries:
        if entry.score is None or entry.slot is None:
            continue
        scored.setdefault(ComponentSlot(entry.slot), []).append(
            (entry.score, EvidenceGrade(entry.grade))
        )

    evidence: list[ComponentEvidence] = []
    for slot, parts in scored.items():
        mean_score = sum(score for score, _ in parts) / len(parts)
        # The weakest grade among the parts governs: a slot is only as reliable as
        # its least-evidenced component, and taking the best would let one
        # well-evidenced CDR vouch for its siblings.
        weakest = min(
            (grade for _, grade in parts),
            key=lambda g: RELIABILITY_BY_GRADE[g],
        )
        evidence.append(
            ComponentEvidence(
                slot=slot, status=EvidenceStatus.AVAILABLE,
                score=mean_score, grade=weakest,
            )
        )

    for slot, score, grade in (
        (ComponentSlot.OBSERVED_PAIRING, pairing_score, pairing_grade),
        (ComponentSlot.STRUCTURE_PRECEDENT, structure_score, structure_grade),
    ):
        if score is None:
            evidence.append(
                ComponentEvidence(
                    slot=slot, status=EvidenceStatus.NOT_FOUND,
                    reason="no score was produced for this construct-level slot",
                )
            )
        else:
            evidence.append(
                ComponentEvidence(
                    slot=slot, status=EvidenceStatus.AVAILABLE,
                    score=score, grade=grade,
                )
            )
    return evidence
