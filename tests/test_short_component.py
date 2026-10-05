"""
Tests for short-component chance correction (SPEC §21).

The statistical content is small but unforgiving: an off-by-one in the p-value or a
BH implementation that is not monotone produces plausible numbers that are wrong,
and nothing downstream would notice. Each formula is checked against values worked
out by hand.
"""

import pytest

from para_scope.models import ComponentType, EvidenceStatus
from para_scope.short_component import (
    MAX_LENGTH_DIFFERENCE_FRACTION,
    MINIMUM_NULL_COMPARISONS,
    SHORT_TYPES,
    benjamini_hochberg,
    compute_short_component_evidence,
    empirical_p_value,
    position_aware_similarity,
)


def _null(n=1000, value=0.2):
    return tuple([value] * n)


class TestPositionAwareSimilarity:
    def test_identical_sequences_are_one(self):
        similarity, compatible, difference = position_aware_similarity("WGGDGFYAMDY", "WGGDGFYAMDY")
        assert similarity == 1.0
        assert compatible is True
        assert difference == 0

    def test_completely_different_equal_length_is_zero(self):
        assert position_aware_similarity("AAAA", "CCCC")[0] == 0.0

    def test_similarity_is_positional_not_compositional(self):
        """An anagram shares every residue and no position."""
        assert position_aware_similarity("ACDE", "EDCA")[0] == 0.0

    def test_half_matching_is_a_half(self):
        assert position_aware_similarity("ACDE", "ACXX")[0] == pytest.approx(0.5)

    def test_similarity_is_over_the_shorter_length(self):
        similarity, _, difference = position_aware_similarity("ACDEF", "ACDE")
        assert similarity == 1.0
        assert difference == 1

    def test_no_gaps_are_introduced(self):
        """Gapping a 13-residue CDR to improve its score is the flattery §21 prevents."""
        # A one-residue insertion at the front destroys positional identity; a
        # gapped aligner would score this near 1.0.
        assert position_aware_similarity("ACDEFG", "XACDEF")[0] < 0.2

    def test_length_incompatibility_is_reported(self):
        _, compatible, difference = position_aware_similarity("A" * 10, "A" * 20)
        assert compatible is False
        assert difference == 10

    def test_compatibility_threshold_is_a_fraction_of_the_longer(self):
        longer = 20
        allowed = int(longer * MAX_LENGTH_DIFFERENCE_FRACTION)
        assert position_aware_similarity("A" * (longer - allowed), "A" * longer)[1] is True
        assert position_aware_similarity("A" * (longer - allowed - 1), "A" * longer)[1] is False

    def test_empty_input_is_not_compatible(self):
        assert position_aware_similarity("", "ACDE")[1] is False


class TestEmpiricalPValue:
    def test_matches_the_spec_formula_by_hand(self):
        """(1 + nulls at least as strong) / (1 + total)."""
        nulls = (0.1, 0.2, 0.3, 0.4)
        # 0.25 is exceeded-or-equalled by 0.3 and 0.4 -> (1+2)/(1+4)
        assert empirical_p_value(0.25, nulls) == pytest.approx(3 / 5)

    def test_an_observation_stronger_than_every_null_is_bounded_away_from_zero(self):
        """A finite null cannot establish impossibility."""
        nulls = _null(99, 0.2)
        assert empirical_p_value(1.0, nulls) == pytest.approx(1 / 100)
        assert empirical_p_value(1.0, nulls) > 0.0

    def test_an_observation_weaker_than_every_null_approaches_one(self):
        nulls = _null(99, 0.9)
        assert empirical_p_value(0.0, nulls) == pytest.approx(100 / 100)

    def test_ties_count_as_at_least_as_strong(self):
        """A null equal to the observation is evidence against it, not for it."""
        assert empirical_p_value(0.5, (0.5,)) == pytest.approx(2 / 2)

    def test_p_value_decreases_as_the_observation_strengthens(self):
        nulls = tuple(i / 100 for i in range(100))
        values = [empirical_p_value(x, nulls) for x in (0.1, 0.5, 0.9)]
        assert values == sorted(values, reverse=True)

    def test_p_value_is_bounded(self):
        nulls = tuple(i / 100 for i in range(100))
        for observed in (0.0, 0.5, 1.0, 2.0):
            assert 0.0 < empirical_p_value(observed, nulls) <= 1.0


class TestBenjaminiHochberg:
    def test_empty_input(self):
        assert benjamini_hochberg([]) == []

    def test_single_p_value_is_unchanged(self):
        assert benjamini_hochberg([0.03]) == pytest.approx([0.03])

    def test_known_example_by_hand(self):
        """BH on [0.01, 0.02, 0.03] with n=3: q = min over ranks of p*n/rank."""
        q = benjamini_hochberg([0.01, 0.02, 0.03])
        assert q[2] == pytest.approx(0.03)        # 0.03*3/3
        assert q[1] == pytest.approx(0.03)        # min(0.02*3/2, 0.03) = 0.03
        assert q[0] == pytest.approx(0.03)        # min(0.01*3/1, 0.03) = 0.03

    def test_order_is_preserved(self):
        p_values = [0.5, 0.001, 0.2]
        q = benjamini_hochberg(p_values)
        # The smallest p must still map to the smallest q, in its original slot.
        assert q[1] == min(q)

    def test_q_values_are_monotone_in_p(self):
        p_values = [0.001, 0.01, 0.02, 0.2, 0.5]
        q = benjamini_hochberg(p_values)
        assert q == sorted(q)

    def test_q_values_are_bounded_by_one(self):
        assert all(value <= 1.0 for value in benjamini_hochberg([0.9, 0.95, 0.99]))

    def test_correction_makes_a_lone_small_p_less_impressive(self):
        """The reason FDR is required: the smallest of many p-values is not news."""
        alone = benjamini_hochberg([0.01])[0]
        among_many = benjamini_hochberg([0.01] + [0.9] * 99)[0]
        assert among_many > alone


class TestChanceCorrectionRefusesWithoutANull:
    def test_no_null_is_insufficient_evidence_not_an_uncorrected_score(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", [("r1", "WGGDGFYAMDY")],
            null_values=None, null_stratum="CDR-H3|len=9-11",
        )
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.calibrated_similarity is None
        assert "without chance correction" in result.reason

    def test_a_null_below_the_minimum_is_refused(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", [("r1", "WGGDGFYAMDY")],
            null_values=_null(MINIMUM_NULL_COMPARISONS - 1), null_stratum="s",
        )
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert str(MINIMUM_NULL_COMPARISONS) in result.reason

    def test_a_long_component_is_unsupported_here(self):
        """VH is scored by alignment and needs no length-matched null."""
        result = compute_short_component_evidence(
            ComponentType.VH, "EVQLVESGG", [], null_values=_null(),
        )
        assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT

    def test_every_cdr_type_is_in_scope(self):
        for component_type in (
            ComponentType.CDR_H1, ComponentType.CDR_H2, ComponentType.CDR_H3,
            ComponentType.CDR_L1, ComponentType.CDR_L2, ComponentType.CDR_L3,
        ):
            assert component_type in SHORT_TYPES


class TestChanceCorrectedEvidence:
    def test_an_exact_match_exceeds_a_weak_null(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", [("r1", "WGGDGFYAMDY")],
            null_values=_null(999, 0.2), null_stratum="CDR-H3|len=9-11",
        )
        assert result.status is EvidenceStatus.AVAILABLE
        assert result.matches[0].position_aware_similarity == 1.0
        assert result.matches[0].p_value == pytest.approx(1 / 1000)
        assert result.calibrated_similarity == pytest.approx(1 - 1 / 1000)

    def test_an_ordinary_match_against_a_strong_null_is_not_significant(self):
        """The whole point: 0.4 identity on a short CDR can be unremarkable."""
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "AAAAACCCCC", [("r1", "AAAAAXXXXX")],
            null_values=_null(999, 0.9), null_stratum="s",
        )
        assert result.calibrated_similarity < 0.5

    def test_length_incompatible_subjects_are_dropped_not_scored(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "A" * 10, [("r1", "A" * 30)],
            null_values=_null(), null_stratum="s",
        )
        assert result.matches == []
        assert result.status is EvidenceStatus.NOT_FOUND
        assert "positionally comparable" in result.reason

    def test_q_values_are_assigned_across_the_searched_set(self):
        subjects = [(f"r{i}", "WGGDGFYAMDY") for i in range(5)]
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", subjects,
            null_values=_null(), null_stratum="s",
        )
        assert len(result.matches) == 5
        assert all(m.q_value is not None for m in result.matches)
        assert all(m.q_value >= m.p_value for m in result.matches)

    def test_matches_are_ordered_by_similarity(self):
        subjects = [("weak", "WGGDGFYAMDX"), ("exact", "WGGDGFYAMDY")]
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", subjects,
            null_values=_null(), null_stratum="s",
        )
        assert result.matches[0].matched_record_id == "exact"

    def test_length_difference_is_recorded_as_a_limitation(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMD", [("r1", "WGGDGFYAMDY")],
            null_values=_null(), null_stratum="s",
        )
        assert any("differs in length" in lim for lim in result.limitations)

    def test_null_provenance_travels_with_the_result(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", [("r1", "WGGDGFYAMDY")],
            null_values=_null(777), null_stratum="CDR-H3|len=9-11",
        )
        assert result.null_comparisons == 777
        assert result.null_stratum == "CDR-H3|len=9-11"

    def test_disclaimer_refuses_a_functional_reading(self):
        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", [("r1", "WGGDGFYAMDY")],
            null_values=_null(), null_stratum="s",
        )
        assert "not a probability that the components are functionally related" in result.disclaimer
        assert "binding" in result.disclaimer

    def test_result_is_json_serializable(self):
        import json
        from dataclasses import asdict

        result = compute_short_component_evidence(
            ComponentType.CDR_H3, "WGGDGFYAMDY", [("r1", "WGGDGFYAMDY")],
            null_values=_null(), null_stratum="s",
        )
        json.dumps(asdict(result), default=str)
