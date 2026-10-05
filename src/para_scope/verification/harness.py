"""
Adversarial audit harness: run the case suite against the real engine and measure it.

The generator proposes inputs; this module decides whether the engine's answers
were defensible. The split is deliberate -- a generator that judged its own
output could not be wrong in a way the harness would notice.

Each expectation is a property, not a golden value. Golden values pin behaviour to
a corpus and break on every legitimate rebuild; properties like "more mutations
must not score closer" hold across corpora and across engine versions, which is
what makes a failure here worth investigating rather than worth re-recording.

Every case also crosses the grounding boundary: its result is rendered as a
GroundedAnswer, so the audit measures prohibited-claim rate and citation coverage
on the engine's actual outputs rather than on hand-written examples.
"""

from dataclasses import dataclass, field
from typing import Optional

from para_scope.calibration import CalibrationSet
from para_scope.corpus import Corpus
from para_scope.germline import SUPPORTED_TYPES, compute_germline_assignment
from para_scope.models import ComponentType, EvidenceStatus
from para_scope.numbering import normalize_query_to_variable_domain
from para_scope.therapeutic_precedent import compute_therapeutic_precedent
from para_scope.verification.adversarial import (
    ADVERSARIAL_VERSION,
    EXPECT_EXPLICIT_STATUS_NOT_A_NUMBER,
    EXPECT_MONOTONE_IN_DISTANCE,
    EXPECT_NO_CONFIDENT_MATCH,
    EXPECT_NORMALIZED_TO_DOMAIN,
    EXPECT_WORSE_THAN_PARENT,
    AdversarialCase,
)
from para_scope.verification.grounding import (
    GROUNDING_VERSION,
    Citation,
    Claim,
    GroundedAnswer,
    ProhibitedClaimError,
    Refusal,
    scan_for_prohibited_claims,
)

HARNESS_VERSION = "1.0.0"

# An identity at or above this, from an input with no relationship to any
# antibody, is an artifact of the scoring method rather than a finding. 0.90 is
# the same threshold Gate 1's regression test used after the fragmentary-alignment
# defect; it is a judgement, recorded as one.
CONFIDENT_MATCH_IDENTITY = 0.90

# Monotonicity is checked with a tolerance: alignment is not a metric, and two
# mutation counts can legitimately produce near-identical global identity when
# substitutions land in regions that align the same way. A *reversal* beyond this
# is what the check is for.
MONOTONE_TOLERANCE = 0.02


@dataclass
class CaseOutcome:
    case_id: str
    family: str
    expectation: str
    passed: bool
    detail: str
    measurements: dict = field(default_factory=dict)


@dataclass
class AuditReport:
    outcomes: list[CaseOutcome] = field(default_factory=list)
    prohibited_claim_count: int = 0
    prohibited_claim_examples: list[str] = field(default_factory=list)
    answers_rendered: int = 0
    answers_refused: int = 0
    citation_coverage: float = 1.0
    uncited_claim_count: int = 0
    crashes: list[str] = field(default_factory=list)

    adversarial_version: str = ADVERSARIAL_VERSION
    grounding_version: str = GROUNDING_VERSION
    harness_version: str = HARNESS_VERSION

    @property
    def total(self) -> int:
        return len(self.outcomes)

    @property
    def passed(self) -> int:
        return sum(1 for o in self.outcomes if o.passed)

    @property
    def failures(self) -> list[CaseOutcome]:
        return [o for o in self.outcomes if not o.passed]

    def by_family(self) -> dict[str, tuple[int, int]]:
        families: dict[str, list[CaseOutcome]] = {}
        for outcome in self.outcomes:
            families.setdefault(outcome.family, []).append(outcome)
        return {
            family: (sum(1 for o in group if o.passed), len(group))
            for family, group in sorted(families.items())
        }


def _best_identity(corpus, case_type, sequence, calibration, top_k=5):
    """Best global_identity the engine will report for a sequence, with its status."""
    result = compute_therapeutic_precedent(
        corpus, case_type, sequence, top_k=top_k,
        shortlist_limit=None, calibration=calibration,
    )
    best = max((m.global_identity for m in result.matches), default=None)
    return result, best


def _answer_for(case: AdversarialCase, result) -> GroundedAnswer:
    """Render the engine's result as a grounded answer.

    Phrased to report evidence and nothing beyond it. If the engine found nothing,
    the answer is an explicit refusal rather than an empty string -- "no therapeutic
    precedent was found" and "we did not look" must not read the same.
    """
    answer = GroundedAnswer(question=f"What therapeutic precedent exists for {case.case_id}?")

    if result.status is not EvidenceStatus.AVAILABLE or not result.matches:
        answer.refusals.append(
            Refusal(
                question=answer.question,
                reason=(
                    f"No therapeutic-precedent evidence is available for this input "
                    f"(status {result.status.value}: {result.reason or 'no matches'})"
                ),
                category=(
                    "no_evidence" if result.status is EvidenceStatus.NOT_FOUND
                    else "insufficient_evidence"
                ),
            )
        )
        return answer

    for match in result.matches[:3]:
        citation = Citation(
            record_id=match.matched_record_id,
            source_database="Thera-SAbDab",
            source_checksum=(
                result.source_snapshots[0].get("raw_checksum", "unknown")
                if getattr(result, "source_snapshots", None) else "unknown"
            ),
            engine_version=match.match_method,
            calibration_version=(
                match.calibrated_similarity.calibration_version
                if match.calibrated_similarity is not None
                and match.calibrated_similarity.status is EvidenceStatus.AVAILABLE
                else None
            ),
        )
        text = (
            f"The submitted {match.matched_component_type.value} aligns to the "
            f"{match.matched_component_type.value} of {match.therapeutic_name} at "
            f"global identity {match.global_identity:.3f} "
            f"(query coverage {match.query_coverage:.2f})."
        )
        if (
            match.calibrated_similarity is not None
            and match.calibrated_similarity.status is EvidenceStatus.AVAILABLE
        ):
            percentile = match.calibrated_similarity
            text += (
                f" That identity sits at the {percentile.percentile:.0f}th percentile "
                f"of {percentile.component_class} similarity to the nearest therapeutic "
                f"in a reference population of {percentile.sample_size} components."
            )
        answer.claims.append(Claim(text=text, citations=(citation,)))

    return answer


def run_suite(
    corpus: Corpus,
    cases: list[AdversarialCase],
    calibration: Optional[CalibrationSet] = None,
) -> AuditReport:
    """Run every case, check its expectation, and measure the grounding boundary."""
    report = AuditReport()
    # Monotone families are judged as a ladder, so results are collected first.
    ladders: dict[tuple[str, str, Optional[str]], list[tuple[int, Optional[float]]]] = {}
    claim_total = 0
    claim_cited = 0

    for case in cases:
        try:
            if case.expectation == EXPECT_EXPLICIT_STATUS_NOT_A_NUMBER:
                report.outcomes.append(_check_degenerate(corpus, case, calibration))
                continue

            result, best = _best_identity(corpus, case.component_type, case.sequence, calibration)

            # --- grounding boundary, measured on real output ---
            answer = _answer_for(case, result)
            claim_total += len(answer.claims)
            claim_cited += sum(1 for c in answer.claims if c.citations)
            try:
                rendered = answer.render()
                report.answers_rendered += 1
                if answer.refusals and not answer.claims:
                    report.answers_refused += 1
                leaked = scan_for_prohibited_claims(rendered)
                if leaked:
                    report.prohibited_claim_count += len(leaked)
                    report.prohibited_claim_examples.append(f"{case.case_id}: {leaked}")
            except ProhibitedClaimError as exc:
                # The boundary did its job -- but an engine output that *would*
                # have asserted something prohibited is itself a finding.
                report.prohibited_claim_count += 1
                report.prohibited_claim_examples.append(f"{case.case_id}: blocked -- {exc}")

            # --- expectation ---
            if case.expectation == EXPECT_NO_CONFIDENT_MATCH:
                report.outcomes.append(_check_no_confident_match(case, best))
            elif case.expectation == EXPECT_WORSE_THAN_PARENT:
                _, parent_best = _best_identity(
                    corpus, case.component_type, case.parent_sequence, calibration
                )
                report.outcomes.append(_check_worse_than_parent(case, best, parent_best))
            elif case.expectation == EXPECT_NORMALIZED_TO_DOMAIN:
                report.outcomes.append(_check_normalized(case, best, corpus, calibration))
            elif case.expectation == EXPECT_MONOTONE_IN_DISTANCE:
                key = (case.family, case.component_type.value, case.parent_record_id)
                ladders.setdefault(key, []).append((case.distance or 0, best))
                report.outcomes.append(
                    CaseOutcome(
                        case_id=case.case_id, family=case.family,
                        expectation=case.expectation, passed=True,
                        detail="recorded for ladder check",
                        measurements={"distance": case.distance, "best_identity": best},
                    )
                )
        except Exception as exc:  # noqa: BLE001 -- a crash is a finding, not an abort
            report.crashes.append(f"{case.case_id}: {type(exc).__name__}: {exc}")
            report.outcomes.append(
                CaseOutcome(
                    case_id=case.case_id, family=case.family,
                    expectation=case.expectation, passed=False,
                    detail=f"raised {type(exc).__name__}: {exc}",
                )
            )

    for (family, component_class, parent), points in sorted(ladders.items(), key=lambda kv: str(kv[0])):
        report.outcomes.append(_check_monotone(family, component_class, parent, points))

    report.citation_coverage = (claim_cited / claim_total) if claim_total else 1.0
    report.uncited_claim_count = claim_total - claim_cited
    return report


# Characters outside the standard 20 make an input not-a-sequence, which is a
# different fact from a sequence the corpus has no match for.
_NON_SEQUENCE_CASES = frozenset({
    "degenerate_digits", "degenerate_lowercase_junk", "degenerate_ambiguity_codes",
    "degenerate_stop_codon_character", "degenerate_gap_characters",
    "degenerate_empty", "degenerate_whitespace_only",
})


def _check_degenerate(corpus, case, calibration) -> CaseOutcome:
    """
    Degenerate input must produce a status, no measurement, and the *right* status.

    The original version of this check only asserted "not AVAILABLE", and passed
    while five of eight degenerate inputs returned NOT_FOUND -- which asserts the
    corpus holds no precedent for a string that was never a sequence. A check that
    accepts any non-success status cannot see that, so it now requires
    UNSUPPORTED_COMPONENT for inputs that are not sequences at all.
    """
    result = compute_therapeutic_precedent(
        corpus, case.component_type, case.sequence,
        top_k=5, shortlist_limit=None, calibration=calibration,
    )
    gave_number = result.status is EvidenceStatus.AVAILABLE and bool(result.matches)
    has_reason = result.reason is not None

    expects_unsupported = case.case_id in _NON_SEQUENCE_CASES
    status_is_right = (
        result.status is EvidenceStatus.UNSUPPORTED_COMPONENT
        if expects_unsupported
        else result.status is not EvidenceStatus.AVAILABLE
    )
    passed = (not gave_number) and has_reason and status_is_right

    detail = (
        f"status={result.status.value} reason={'set' if has_reason else 'MISSING'} "
        f"matches={len(result.matches)}"
    )
    if gave_number:
        detail += " -- scored a degenerate input"
    elif expects_unsupported and not status_is_right:
        detail += (
            " -- not a sequence, so NOT_FOUND wrongly claims the corpus lacks "
            "precedent for it; expected unsupported_component"
        )
    return CaseOutcome(
        case_id=case.case_id, family=case.family, expectation=case.expectation,
        passed=passed, detail=detail,
        measurements={"status": result.status.value, "match_count": len(result.matches),
                      "expected_unsupported": expects_unsupported},
    )


def _check_no_confident_match(case, best) -> CaseOutcome:
    if best is None:
        return CaseOutcome(
            case_id=case.case_id, family=case.family, expectation=case.expectation,
            passed=True, detail="no match reported at all",
            measurements={"best_identity": None},
        )
    passed = best < CONFIDENT_MATCH_IDENTITY
    return CaseOutcome(
        case_id=case.case_id, family=case.family, expectation=case.expectation,
        passed=passed,
        detail=(
            f"best global_identity {best:.3f} "
            f"{'<' if passed else '>='} {CONFIDENT_MATCH_IDENTITY}"
        ),
        measurements={"best_identity": best},
    )


def _check_worse_than_parent(case, best, parent_best) -> CaseOutcome:
    if parent_best is None:
        return CaseOutcome(
            case_id=case.case_id, family=case.family, expectation=case.expectation,
            passed=True, detail="parent itself had no match; nothing to compare",
        )
    if best is None:
        return CaseOutcome(
            case_id=case.case_id, family=case.family, expectation=case.expectation,
            passed=True, detail=f"no match reported (parent reached {parent_best:.3f})",
            measurements={"best_identity": None, "parent_best": parent_best},
        )
    passed = best < parent_best - MONOTONE_TOLERANCE
    return CaseOutcome(
        case_id=case.case_id, family=case.family, expectation=case.expectation,
        passed=passed,
        detail=f"best {best:.3f} vs parent {parent_best:.3f}",
        measurements={"best_identity": best, "parent_best": parent_best},
    )


def _check_normalized(case, best, corpus, calibration) -> CaseOutcome:
    """A doubled domain must be reduced to one before comparison."""
    _, note = normalize_query_to_variable_domain(
        case.component_type.value, case.sequence
    )
    _, parent_best = _best_identity(
        corpus, case.component_type, case.parent_sequence, calibration
    )
    normalized = note is not None
    matches_parent = (
        best is not None and parent_best is not None
        and abs(best - parent_best) <= MONOTONE_TOLERANCE
    )
    passed = normalized and matches_parent
    return CaseOutcome(
        case_id=case.case_id, family=case.family, expectation=case.expectation,
        passed=passed,
        detail=(
            f"normalization {'reported' if normalized else 'NOT reported'}; "
            f"doubled {best if best is None else round(best, 3)} vs "
            f"parent {parent_best if parent_best is None else round(parent_best, 3)}"
        ),
        measurements={"best_identity": best, "parent_best": parent_best,
                      "normalization_note": normalized},
    )


def _check_monotone(family, component_class, parent, points) -> CaseOutcome:
    """Identity must not increase as distance from the parent increases."""
    ordered = sorted(points, key=lambda p: p[0])
    scored = [(d, v) for d, v in ordered if v is not None]
    case_id = f"ladder_{family}_{component_class}_{parent or 'none'}"
    if len(scored) < 2:
        return CaseOutcome(
            case_id=case_id, family=f"{family}_ladder",
            expectation=EXPECT_MONOTONE_IN_DISTANCE, passed=True,
            detail=f"only {len(scored)} scored point(s); nothing to compare",
        )
    reversals = [
        (scored[i][0], scored[i][1], scored[i + 1][0], scored[i + 1][1])
        for i in range(len(scored) - 1)
        if scored[i + 1][1] > scored[i][1] + MONOTONE_TOLERANCE
    ]
    return CaseOutcome(
        case_id=case_id, family=f"{family}_ladder",
        expectation=EXPECT_MONOTONE_IN_DISTANCE,
        passed=not reversals,
        detail=(
            "monotone: " + " -> ".join(f"d{d}:{v:.3f}" for d, v in scored)
            if not reversals
            else f"{len(reversals)} reversal(s): " + "; ".join(
                f"d{a}:{av:.3f} then d{b}:{bv:.3f}" for a, av, b, bv in reversals
            )
        ),
        measurements={"ladder": scored},
    )
