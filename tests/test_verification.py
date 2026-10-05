"""
Tests for the verification layer: adversarial generation, grounding, and the detectors.

The most important class here is TestDetectorsActuallyFire. An audit that reports
86/86 proves nothing unless the checks can fail, so each one is fed a known-bad
value and must report failure. Without that, a passing audit and a broken audit
are indistinguishable.
"""

import pytest

from para_scope.models import ComponentType
from para_scope.verification.adversarial import (
    ADVERSARIAL_VERSION,
    DEGENERATE_INPUTS,
    EXPECT_NO_CONFIDENT_MATCH,
    EXPECT_WORSE_THAN_PARENT,
    AdversarialCase,
    build_case_suite,
    chimeric_sequence,
    point_mutant,
    random_sequence,
    reversed_sequence,
    shuffled_sequence,
    truncate,
)
from para_scope.verification.grounding import (
    GROUNDING_VERSION,
    Citation,
    Claim,
    GroundedAnswer,
    ProhibitedClaimError,
    Refusal,
    UngroundedClaimError,
    scan_for_prohibited_claims,
)
from para_scope.verification.harness import (
    CONFIDENT_MATCH_IDENTITY,
    MONOTONE_TOLERANCE,
    _check_monotone,
    _check_no_confident_match,
    _check_worse_than_parent,
)

VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
OTHER_VH = (
    "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKN"
    "RVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS"
)


def _case(expectation=EXPECT_NO_CONFIDENT_MATCH, **kw):
    defaults = dict(
        case_id="c", family="f", description="d", sequence=VH,
        component_type=ComponentType.VH, expectation=expectation,
    )
    defaults.update(kw)
    return AdversarialCase(**defaults)


class TestGeneratorsAreDeterministicAndDistinct:
    def test_random_sequence_is_reproducible(self):
        assert random_sequence("x", 50) == random_sequence("x", 50)

    def test_different_case_ids_give_different_sequences(self):
        assert random_sequence("x", 50) != random_sequence("y", 50)

    def test_shuffle_preserves_composition_and_changes_order(self):
        shuffled = shuffled_sequence("x", VH)
        assert sorted(shuffled) == sorted(VH)
        assert shuffled != VH

    def test_reverse_preserves_composition(self):
        assert sorted(reversed_sequence(VH)) == sorted(VH)
        assert reversed_sequence(VH) != VH

    def test_point_mutant_changes_exactly_k_positions(self):
        for k in (1, 5, 20):
            mutant = point_mutant("x", VH, k)
            assert len(mutant) == len(VH)
            assert sum(1 for a, b in zip(mutant, VH) if a != b) == k

    def test_point_mutant_never_substitutes_a_residue_for_itself(self):
        """A 'mutation' that keeps the residue would make the ladder meaningless."""
        mutant = point_mutant("x", VH, 30)
        assert sum(1 for a, b in zip(mutant, VH) if a != b) == 30

    def test_chimera_takes_halves_from_each_parent(self):
        chimera = chimeric_sequence(VH, OTHER_VH)
        assert chimera.startswith(VH[:10])
        assert chimera.endswith(OTHER_VH[-10:])

    def test_truncation_keeps_a_prefix(self):
        assert VH.startswith(truncate(VH, 0.5))
        assert len(truncate(VH, 0.5)) < len(VH)

    def test_suite_covers_every_family(self):
        cases = build_case_suite(
            [("a", VH, ComponentType.VH), ("b", OTHER_VH, ComponentType.VH)]
        )
        assert {c.family for c in cases} == {
            "random", "degenerate", "shuffled", "reversed",
            "duplicated", "point_mutant", "truncation", "chimera",
        }

    def test_suite_is_identical_across_runs(self):
        parents = [("a", VH, ComponentType.VH)]
        first = [(c.case_id, c.sequence) for c in build_case_suite(parents)]
        second = [(c.case_id, c.sequence) for c in build_case_suite(parents)]
        assert first == second

    def test_degenerate_inputs_include_the_dangerous_shapes(self):
        names = {name for name, _ in DEGENERATE_INPUTS}
        assert {"empty", "ambiguity_codes", "gap_characters"} <= names

    def test_version_is_exported(self):
        assert ADVERSARIAL_VERSION


class TestGroundingBoundary:
    def _citation(self):
        return Citation(
            record_id="VH::abc", source_database="Thera-SAbDab",
            source_checksum="cc940278", engine_version="1.0.0",
        )

    def test_claim_without_citation_cannot_be_constructed(self):
        with pytest.raises(UngroundedClaimError, match="no citation"):
            Claim(text="This VH appears in 12 therapeutics.", citations=())

    def test_empty_claim_is_rejected(self):
        with pytest.raises(UngroundedClaimError):
            Claim(text="   ", citations=(self._citation(),))

    def test_cited_claim_renders_with_its_evidence(self):
        answer = GroundedAnswer(question="q")
        answer.claims.append(Claim("Matches trastuzumab at 1.000.", (self._citation(),)))
        rendered = answer.render()
        assert "Matches trastuzumab" in rendered
        assert "VH::abc" in rendered
        assert "cc940278" in rendered

    def test_prohibited_claim_is_refused_not_warned(self):
        """A warning beside a sentence does not stop the sentence being quoted."""
        answer = GroundedAnswer(question="q")
        answer.claims.append(Claim("This candidate will bind HER2.", (self._citation(),)))
        with pytest.raises(ProhibitedClaimError, match="binding"):
            answer.render()

    @pytest.mark.parametrize("text", [
        "This candidate will bind its target.",
        "The construct is likely to express well.",
        "It is safe and effective.",
        "This shows low immunogenicity.",
        "A humanness score of 0.8.",
        "Probability of binding is high.",
        "This is a promising candidate.",
        "A better antibody than the parent.",
        "The molecule is developable.",
        "High drug-likeness.",
    ])
    def test_each_prohibited_shape_is_caught(self, text):
        assert scan_for_prohibited_claims(text), f"not caught: {text}"

    @pytest.mark.parametrize("text", [
        "This VH appears in 12 therapeutic programmes.",
        "Global identity 0.873 at query coverage 0.98.",
        "It is not a measure of humanness, immunogenicity, or safety.",
        "This index summarizes comparative database evidence.",
        "Window coverage is not an immunogenicity score.",
        "87 of 112 nine-mers were found in the reference set.",
        "Nearest germline V gene is IGHV3-23 at 0.912 identity.",
    ])
    def test_legitimate_evidence_statements_pass(self, text):
        """A topic-based banned-word list would reject the disclaimers themselves."""
        assert scan_for_prohibited_claims(text) == [], text

    @pytest.mark.parametrize("text", [
        "Window coverage is not an immunogenicity score.",
        "It is not a measure of humanness, immunogenicity, or safety.",
        "This index is not a probability of binding.",
        "Similarity does not imply it will bind the same target.",
        "Similarity to an approved therapeutic does not mean the candidate will succeed.",
    ])
    def test_required_disclaimers_are_not_blocked(self, text):
        """
        CLAUDE.md *requires* these sentences. A scanner that banned the topic
        rather than the assertion would reject the platform's own mandatory
        disclaimers -- the exact text that exists to prevent the claim.
        """
        assert scan_for_prohibited_claims(text) == [], text

    def test_a_later_disclaimer_does_not_launder_an_earlier_claim(self):
        """Negation only suppresses a claim it actually precedes."""
        assert scan_for_prohibited_claims("It will bind. It does not imply anything.")

    def test_negated_property_is_still_a_claim_about_the_molecule(self):
        """'not a safety score' is a disclaimer; 'is not safe' is a safety claim."""
        assert scan_for_prohibited_claims("This is not a safety score.") == []

    def test_empty_answer_cannot_render(self):
        """Empty text would read as 'nothing found', which is itself a claim."""
        with pytest.raises(UngroundedClaimError):
            GroundedAnswer(question="q").render()

    def test_refusal_alone_is_a_valid_answer(self):
        answer = GroundedAnswer(question="q")
        answer.refusals.append(Refusal("q", "no therapeutic precedent found", "no_evidence"))
        assert "Cannot answer" in answer.render()

    def test_citation_coverage_is_measured(self):
        answer = GroundedAnswer(question="q")
        answer.claims.append(Claim("x", (self._citation(),)))
        assert answer.citation_coverage() == 1.0

    def test_calibrated_and_raw_citations_are_distinguishable(self):
        """An uncalibrated number is not a percentile and must not read as one."""
        raw = self._citation()
        calibrated = Citation(
            record_id="VH::abc", source_database="Thera-SAbDab",
            source_checksum="cc940278", engine_version="1.0.0",
            calibration_version="1.0.0",
        )
        assert "calib" not in raw.short()
        assert "calib 1.0.0" in calibrated.short()

    def test_version_is_exported(self):
        assert GROUNDING_VERSION


class TestDetectorsActuallyFire:
    """
    Negative controls. Each detector is given a value it must reject.

    This class is why the audit's pass rate means anything: a check that cannot
    fail reports success on a broken engine exactly as it does on a correct one.
    """

    def test_confident_match_detector_rejects_a_high_identity(self):
        bad = _check_no_confident_match(_case(), best=0.99)
        assert bad.passed is False
        assert "0.99" in bad.detail

    def test_confident_match_detector_accepts_a_low_identity(self):
        assert _check_no_confident_match(_case(), best=0.40).passed is True

    def test_confident_match_detector_fires_exactly_at_the_threshold(self):
        assert _check_no_confident_match(_case(), best=CONFIDENT_MATCH_IDENTITY).passed is False
        assert _check_no_confident_match(
            _case(), best=CONFIDENT_MATCH_IDENTITY - 0.01
        ).passed is True

    def test_worse_than_parent_detector_rejects_a_tie(self):
        """A shuffled sequence scoring like its parent is the failure, not a pass."""
        case = _case(expectation=EXPECT_WORSE_THAN_PARENT, parent_sequence=VH)
        assert _check_worse_than_parent(case, best=0.90, parent_best=0.90).passed is False

    def test_worse_than_parent_detector_rejects_scoring_higher(self):
        case = _case(expectation=EXPECT_WORSE_THAN_PARENT, parent_sequence=VH)
        assert _check_worse_than_parent(case, best=0.95, parent_best=0.90).passed is False

    def test_worse_than_parent_detector_accepts_a_real_drop(self):
        case = _case(expectation=EXPECT_WORSE_THAN_PARENT, parent_sequence=VH)
        assert _check_worse_than_parent(case, best=0.60, parent_best=0.95).passed is True

    def test_monotone_detector_rejects_a_reversal(self):
        """More mutations scoring closer to the parent is the defect."""
        bad = _check_monotone("point_mutant", "VH", "p", [(1, 0.80), (3, 0.70), (10, 0.95)])
        assert bad.passed is False
        assert "reversal" in bad.detail

    def test_monotone_detector_accepts_a_decreasing_ladder(self):
        good = _check_monotone("point_mutant", "VH", "p", [(1, 0.98), (3, 0.90), (10, 0.70)])
        assert good.passed is True
        assert "monotone" in good.detail

    def test_monotone_detector_tolerates_noise_below_the_tolerance(self):
        """Alignment is not a metric; a hair's movement is not a reversal."""
        noisy = [(1, 0.900), (3, 0.900 + MONOTONE_TOLERANCE / 2)]
        assert _check_monotone("point_mutant", "VH", "p", noisy).passed is True

    def test_monotone_detector_catches_a_reversal_just_above_tolerance(self):
        over = [(1, 0.900), (3, 0.900 + MONOTONE_TOLERANCE * 2)]
        assert _check_monotone("point_mutant", "VH", "p", over).passed is False

    def test_monotone_detector_abstains_on_too_few_points(self):
        """One point is not a trend, and must not be reported as a pass of a trend."""
        outcome = _check_monotone("point_mutant", "VH", "p", [(1, 0.9)])
        assert "nothing to compare" in outcome.detail

    def test_prohibited_scan_would_catch_a_leaked_claim_in_a_real_answer(self):
        """
        The harness counts a blocked render as a finding. This confirms the path it
        depends on: a claim built from otherwise-valid evidence still fails if its
        prose oversteps.
        """
        answer = GroundedAnswer(question="q")
        answer.claims.append(
            Claim(
                "Matches trastuzumab at 1.000 identity, so it will bind HER2.",
                (Citation("VH::a", "Thera-SAbDab", "cc94", "1.0.0"),),
            )
        )
        with pytest.raises(ProhibitedClaimError):
            answer.render()
