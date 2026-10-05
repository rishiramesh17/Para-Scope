"""
Short-component chance correction (SPEC §21).

Why an alignment score is not enough here
-----------------------------------------
A 13-residue CDR-H3 and a 120-residue VH are not the same measurement problem.
Two unrelated VH domains share perhaps 70% identity and that number means
something; two unrelated CDR-H3s of the same length can share 40% by accident,
because there are few positions and antibody CDRs draw from a biased residue
pool. Reporting 0.40 identity for a CDR as though it were 0.40 for a domain
overstates it, which is why CDRs are excluded from the alignment index entirely
(`indexes/approximate.ELIGIBLE_TYPES`) and why that exclusion left them with no
score at all until this module.

SPEC §21's answer is an empirical null: compare the observed identity against the
distribution of identities that *unrelated* components of the same type and length
actually achieve, and report where it falls.

    empirical p = (1 + nulls at least as strong as observed) / (1 + total nulls)

The +1 on both sides is the spec's, and it matters: it bounds p away from zero, so
an observation stronger than every null gets a p of 1/(N+1) rather than 0. Claiming
p = 0 from a finite null would assert impossibility from a sample.

The null is real components, not random strings
-----------------------------------------------
The background is drawn from pairs of *distinct, unrelated* real components of the
same type and length band. Random strings would make every real match look
spectacular, because real CDRs share compositional biases that random sequence does
not -- the null has to contain that bias or it is not a null for this question.

Position-aware, and length compatibility reported separately
------------------------------------------------------------
§21 asks for position-aware similarity and length compatibility as separate
outputs. They are separate facts: two CDR-H3s of different lengths may be related
and are not positionally comparable, so identity is computed only over the
compatible span and the length difference is reported rather than folded into the
score. Equal length is not required, but it is stated.

Multiple testing
----------------
A query is compared against every same-type component in the corpus, so the
smallest p-value among thousands of comparisons is not surprising on its own. §21
requires FDR correction across the searched set; Benjamini-Hochberg q-values are
reported alongside the raw p-values, and the q-value is the one to interpret.
"""

from dataclasses import dataclass, field
from typing import Optional

from para_scope.models import ComponentType, EvidenceStatus

# Bump on any change to the null construction, the p-value, or the FDR method.
#   1.0.0 -- initial: empirical null over unrelated same-type same-band components
SHORT_COMPONENT_VERSION = "1.0.0"

# Measurement id of the null distributions in the calibration artifact. Named here
# so producer (scripts/build_short_component_nulls.py) and consumer cannot drift.
NULL_MEASUREMENT_ID = "null.short_component_identity"

# Short component types. These are exactly the types excluded from the alignment
# index, which is the reason this module exists.
SHORT_TYPES = frozenset({
    ComponentType.CDR_H1, ComponentType.CDR_H2, ComponentType.CDR_H3,
    ComponentType.CDR_L1, ComponentType.CDR_L2, ComponentType.CDR_L3,
})

# Fewest null comparisons for a p-value to be reported. Below this the finest
# p-value obtainable is too coarse to mean anything: with 50 nulls the smallest
# non-trivial p is 1/51, so a "significant" result cannot be distinguished from
# the resolution limit.
MINIMUM_NULL_COMPARISONS = 200

# Largest length difference still treated as positionally comparable, as a
# fraction of the longer sequence. Beyond it identity over the shared span
# describes a different pair of objects than the ones submitted.
#
# PROVISIONAL, like the coverage floor and the grade thresholds.
MAX_LENGTH_DIFFERENCE_FRACTION = 0.25

SHORT_COMPONENT_DISCLAIMER = (
    "A short-component p-value states how often unrelated components of the same "
    "type and length reach the observed similarity by chance in this reference "
    "set. It is not a probability that the components are functionally related, "
    "and it says nothing about binding, specificity, or any property of either "
    "construct."
)


@dataclass
class ShortComponentMatch:
    """One short-component comparison with its chance correction."""

    matched_record_id: str
    component_type: ComponentType

    # §21's two required measurements, kept separate.
    position_aware_similarity: float
    length_compatible: bool
    query_length: int
    matched_length: int
    length_difference: int

    p_value: Optional[float] = None
    # Benjamini-Hochberg q-value across the searched set. This is the number to
    # interpret; a raw p-value over thousands of comparisons is not.
    q_value: Optional[float] = None
    null_comparisons: Optional[int] = None


@dataclass
class ShortComponentResult:
    """Chance-corrected short-component evidence for one query."""

    status: EvidenceStatus
    component_type: ComponentType
    query_length: int = 0
    matches: list[ShortComponentMatch] = field(default_factory=list)

    # Percentile of the best observed similarity within the null, in [0, 1]. This
    # is what a calibrated similarity term may legitimately consume for a short
    # component: it is the observed value placed in a chance distribution, which
    # is what "calibrated" means. None whenever no null was available.
    calibrated_similarity: Optional[float] = None

    null_stratum: Optional[str] = None
    null_comparisons: Optional[int] = None
    engine_version: str = SHORT_COMPONENT_VERSION
    reason: Optional[str] = None
    limitations: list[str] = field(default_factory=list)
    disclaimer: str = SHORT_COMPONENT_DISCLAIMER


def position_aware_similarity(query: str, subject: str) -> tuple[float, bool, int]:
    """
    Positional identity over the comparable span, plus length compatibility.

    Returns (similarity, length_compatible, length_difference). Similarity is
    computed over the shorter length -- no gaps are introduced, because inserting
    gaps into a 13-residue CDR to improve its score is exactly the kind of
    flattery an empirical null exists to prevent.
    """
    if not query or not subject:
        return 0.0, False, abs(len(query) - len(subject))

    longer = max(len(query), len(subject))
    shorter = min(len(query), len(subject))
    difference = longer - shorter
    compatible = (difference / longer) <= MAX_LENGTH_DIFFERENCE_FRACTION

    identities = sum(1 for a, b in zip(query, subject) if a == b)
    return identities / shorter, compatible, difference


def empirical_p_value(observed: float, null_values: tuple[float, ...]) -> float:
    """
    SPEC §21's empirical p-value, verbatim.

        (1 + nulls at least as strong as observed) / (1 + total nulls)

    The +1 on both sides bounds p away from zero: an observation stronger than
    every null gets 1/(N+1), never 0. A finite null cannot establish impossibility.
    """
    at_least_as_strong = sum(1 for value in null_values if value >= observed)
    return (1 + at_least_as_strong) / (1 + len(null_values))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """
    BH step-up q-values, order preserved.

    A query is compared against every same-type component, so the smallest p among
    thousands is unremarkable by itself. Without this a short-component result
    would report a "significant" match for essentially every query.
    """
    count = len(p_values)
    if count == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda pair: pair[1])
    q_values = [0.0] * count
    previous = 1.0
    # Walk from the largest p downward, enforcing monotonicity.
    for rank in range(count, 0, -1):
        original_index, p = indexed[rank - 1]
        previous = min(previous, p * count / rank)
        q_values[original_index] = min(1.0, previous)
    return q_values


def compute_short_component_evidence(
    component_type: ComponentType,
    query: str,
    subjects: list[tuple[str, str]],
    null_values: Optional[tuple[float, ...]] = None,
    null_stratum: Optional[str] = None,
) -> ShortComponentResult:
    """
    Compare a short component against candidate subjects with §21 chance correction.

    `subjects` is (record_id, sequence). `null_values` is the empirical background
    for this component type and length band, from the calibration artifact. Without
    a null the result is INSUFFICIENT_EVIDENCE: a similarity for a 13-residue
    sequence reported without chance correction is the defect §21 exists to prevent,
    so it is withheld rather than reported bare.
    """
    cleaned = query.strip().upper()
    result = ShortComponentResult(
        status=EvidenceStatus.NOT_FOUND,
        component_type=component_type,
        query_length=len(cleaned),
        null_stratum=null_stratum,
        null_comparisons=len(null_values) if null_values else None,
    )

    if component_type not in SHORT_TYPES:
        result.status = EvidenceStatus.UNSUPPORTED_COMPONENT
        result.reason = (
            f"{component_type.value} is not a short component; it is scored by "
            f"alignment, which needs no length-matched null."
        )
        return result

    if not cleaned:
        result.status = EvidenceStatus.UNSUPPORTED_COMPONENT
        result.reason = "Empty sequence."
        return result

    if not null_values or len(null_values) < MINIMUM_NULL_COMPARISONS:
        result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
        available = len(null_values) if null_values else 0
        result.reason = (
            f"No usable empirical null for stratum {null_stratum!r} "
            f"({available} comparisons, {MINIMUM_NULL_COMPARISONS} required). A "
            f"similarity over {len(cleaned)} residues reported without chance "
            f"correction is the defect SPEC §21 exists to prevent, so it is "
            f"withheld rather than reported uncorrected."
        )
        return result

    measured: list[ShortComponentMatch] = []
    for record_id, subject in subjects:
        similarity, compatible, difference = position_aware_similarity(cleaned, subject)
        if not compatible:
            # Length-incompatible subjects are dropped rather than scored over a
            # shared span: that span describes a different pair of objects.
            continue
        measured.append(
            ShortComponentMatch(
                matched_record_id=record_id,
                component_type=component_type,
                position_aware_similarity=similarity,
                length_compatible=compatible,
                query_length=len(cleaned),
                matched_length=len(subject),
                length_difference=difference,
                p_value=empirical_p_value(similarity, null_values),
                null_comparisons=len(null_values),
            )
        )

    if not measured:
        result.reason = (
            f"No same-type component was within "
            f"{MAX_LENGTH_DIFFERENCE_FRACTION:.0%} of the query's length, so none "
            f"is positionally comparable to it."
        )
        return result

    for match, q_value in zip(measured, benjamini_hochberg([m.p_value for m in measured])):
        match.q_value = q_value

    measured.sort(key=lambda m: (-m.position_aware_similarity, m.matched_record_id))
    result.matches = measured
    result.status = EvidenceStatus.AVAILABLE

    # Percentile of the best observation within the null. 1 - p is exactly that:
    # the fraction of the background the observation exceeds.
    best = measured[0]
    result.calibrated_similarity = 1.0 - best.p_value

    if best.length_difference:
        result.limitations.append(
            f"The best match differs in length by {best.length_difference} "
            f"residue(s); similarity is positional over the shorter span and no "
            f"gaps were introduced."
        )
    return result
