"""
Tests for whole-construct profile assembly (Gate 11, SPEC §37).

The load-bearing test here is
`test_index_is_identical_at_full_and_partial_coverage`. §37's geometric mean
divides by the weight that contributed, so dropping components shrinks numerator
and denominator together and the index comes out unchanged. An incomplete index is
numerically indistinguishable from a complete one, and `evidence_coverage` is the
only thing that distinguishes them. If that test ever fails because coverage
stopped being reported, every index this platform emits starts overstating itself.
"""

import math

import pytest

from para_scope.construct_format import parse_format
from para_scope.evidence_strength import (
    DimensionReport,
    EvidenceGrade,
    assess_evidence,
)
from para_scope.models import EvidenceClass, EvidenceStatus
from para_scope.profile import (
    INDEX_DISCLAIMER,
    INSUFFICIENT_DISPLAY,
    PROFILE_VERSION,
    RELIABILITY_BY_GRADE,
    REQUIRED_COMPONENTS,
    SCORE_EPSILON,
    WEIGHT_SET_VERSION,
    WEIGHT_SETS,
    ComponentEvidence,
    ComponentSlot,
    CoverageScope,
    FormatClass,
    assemble_index,
    classify_format,
)


def ev(slot, score, grade=EvidenceGrade.A, status=EvidenceStatus.AVAILABLE):
    return ComponentEvidence(slot=slot, status=status, score=score, grade=grade)


def _igg_components(score=90, exclude=()):
    return [
        ev(slot, score) for slot in WEIGHT_SETS[FormatClass.IGG_LIKE]
        if slot not in exclude
    ]


WHOLE_MAB = parse_format("Whole mAb")
NANOBODY = parse_format("Nanobody")


class TestWeightSets:
    def test_igg_weights_match_the_spec_table_and_sum_to_one(self):
        weights = WEIGHT_SETS[FormatClass.IGG_LIKE]
        assert weights[ComponentSlot.VH] == 0.14
        assert weights[ComponentSlot.VL] == 0.14
        assert weights[ComponentSlot.CDR_H3] == 0.14
        assert weights[ComponentSlot.CDR_H1_H2] == 0.10
        assert weights[ComponentSlot.CDR_L1_L2_L3] == 0.12
        assert weights[ComponentSlot.HEAVY_FRAMEWORKS] == 0.08
        assert weights[ComponentSlot.LIGHT_FRAMEWORKS] == 0.08
        assert weights[ComponentSlot.OBSERVED_PAIRING] == 0.10
        assert weights[ComponentSlot.STRUCTURE_PRECEDENT] == 0.05
        assert weights[ComponentSlot.CONSTANT_FC] == 0.05
        assert sum(weights.values()) == pytest.approx(1.00)

    def test_every_weight_set_sums_to_one(self):
        for format_class, weights in WEIGHT_SETS.items():
            assert sum(weights.values()) == pytest.approx(1.00), format_class

    def test_single_domain_set_has_no_light_chain_or_pairing_slots(self):
        """A VHH has no light chain by architecture, not by missing data."""
        weights = WEIGHT_SETS[FormatClass.SINGLE_DOMAIN]
        for absent in (
            ComponentSlot.VL, ComponentSlot.CDR_L1_L2_L3,
            ComponentSlot.LIGHT_FRAMEWORKS, ComponentSlot.OBSERVED_PAIRING,
        ):
            assert absent not in weights

    def test_weight_set_version_names_the_format_it_was_written_for(self):
        """Using the IgG-like set for a VHH is the mistake §37 warns about."""
        assert "igg" in WEIGHT_SET_VERSION.lower()
        assert "provisional" in WEIGHT_SET_VERSION


class TestFormatClassification:
    @pytest.mark.parametrize("raw,expected", [
        ("Whole mAb", FormatClass.IGG_LIKE),
        ("Fab", FormatClass.IGG_LIKE),
        ("Nanobody", FormatClass.SINGLE_DOMAIN),
        ("scFv", FormatClass.UNSUPPORTED),
        ("Bispecific Whole mAb", FormatClass.UNSUPPORTED),
        ("Bispecific Mixed mAb and scFv", FormatClass.UNSUPPORTED),
        ("Fusion Protein", FormatClass.UNSUPPORTED),
    ])
    def test_classification(self, raw, expected):
        assert classify_format(parse_format(raw)) is expected

    def test_undeclared_format_is_unsupported_not_assumed_igg(self):
        assert classify_format(None) is FormatClass.UNSUPPORTED

    def test_unsupported_format_withholds_the_index(self):
        """§37: the provisional weights must not be reused for other formats."""
        result = assemble_index(_igg_components(), parse_format("scFv"))
        assert result.index is None
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert any("weight set" in f for f in result.gating_failures)


class TestCoverageMakesOmissionVisible:
    def test_index_is_identical_at_full_and_partial_coverage(self):
        """
        The reason coverage reporting exists.

        §37's geometric mean divides by the weight that contributed, so removing
        components shrinks numerator and denominator together. With equal
        per-component scores the index is *exactly* unchanged -- an index resting
        on 42% of its weight set reads identically to one resting on 100%.
        """
        full = assemble_index(_igg_components(90), WHOLE_MAB)
        partial = assemble_index(
            [ev(ComponentSlot.VH, 90), ev(ComponentSlot.VL, 90), ev(ComponentSlot.CDR_H3, 90)],
            WHOLE_MAB,
        )
        assert full.index == pytest.approx(partial.index)
        assert full.evidence_coverage == pytest.approx(1.0)
        assert partial.evidence_coverage < 0.5
        # The numbers are the same; only coverage separates them.
        assert full.evidence_coverage != partial.evidence_coverage

    def test_coverage_is_measured_against_the_full_weight_set(self):
        """Measured against contributing weight it would always be 1.0."""
        result = assemble_index(
            _igg_components(exclude={ComponentSlot.CONSTANT_FC}), WHOLE_MAB
        )
        assert result.evidence_coverage == pytest.approx(0.95)

    def test_partial_coverage_is_recorded_as_a_limitation(self):
        result = assemble_index(
            _igg_components(exclude={ComponentSlot.CONSTANT_FC}), WHOLE_MAB
        )
        assert any("indistinguishable" in lim for lim in result.limitations)

    def test_every_missing_slot_carries_a_reason(self):
        result = assemble_index(
            _igg_components(exclude={ComponentSlot.CONSTANT_FC, ComponentSlot.OBSERVED_PAIRING}),
            WHOLE_MAB,
        )
        assert set(result.components_missing) == {
            ComponentSlot.CONSTANT_FC.value, ComponentSlot.OBSERVED_PAIRING.value
        }
        for reason in result.components_missing.values():
            assert reason.strip()

    def test_present_and_missing_partition_the_weight_set(self):
        result = assemble_index(
            _igg_components(exclude={ComponentSlot.CONSTANT_FC}), WHOLE_MAB
        )
        covered = set(result.components_present) | set(result.components_missing)
        assert covered == {s.value for s in WEIGHT_SETS[FormatClass.IGG_LIKE]}


class TestScopeIsDeclaredAndTheNameFollows:
    def test_without_constant_evidence_scope_is_variable_domain_only(self):
        result = assemble_index(
            _igg_components(exclude={ComponentSlot.CONSTANT_FC}), WHOLE_MAB
        )
        assert result.coverage_scope is CoverageScope.VARIABLE_DOMAIN_ONLY
        assert result.displayed_name == "Comparative Evidence Index (variable domain only)"

    def test_with_constant_evidence_the_qualifier_disappears(self):
        result = assemble_index(_igg_components(), WHOLE_MAB)
        assert result.coverage_scope is CoverageScope.FULL_CONSTRUCT
        assert result.displayed_name == "Whole-Construct Comparative Evidence Index"

    def test_the_concept_is_never_renamed_only_qualified(self):
        """CLAUDE.md: the term names the spec's concept, not today's data gap."""
        for result in (
            assemble_index(_igg_components(), WHOLE_MAB),
            assemble_index(_igg_components(exclude={ComponentSlot.CONSTANT_FC}), WHOLE_MAB),
        ):
            assert "Comparative Evidence Index" in result.displayed_name

    def test_constant_is_required_for_full_scope_and_not_for_variable(self):
        """"Critical" is per format and per scope, not absolute."""
        variable = REQUIRED_COMPONENTS[
            (FormatClass.IGG_LIKE, CoverageScope.VARIABLE_DOMAIN_ONLY)
        ]
        full = REQUIRED_COMPONENTS[(FormatClass.IGG_LIKE, CoverageScope.FULL_CONSTRUCT)]
        assert ComponentSlot.CONSTANT_FC not in variable
        assert ComponentSlot.CONSTANT_FC in full


class TestGating:
    def test_missing_required_component_withholds_the_index(self):
        result = assemble_index(
            [ev(ComponentSlot.VH, 90), ev(ComponentSlot.CDR_H3, 90)], WHOLE_MAB
        )
        assert result.index is None
        assert ComponentSlot.VL.value in result.missing_critical_evidence

    def test_the_required_display_string_is_exact(self):
        """§37 specifies this string."""
        result = assemble_index([ev(ComponentSlot.VH, 90)], WHOLE_MAB)
        assert result.display == INSUFFICIENT_DISPLAY
        assert result.display == "Insufficient evidence for a whole-construct index."

    def test_an_infrastructural_failure_withholds_the_index(self):
        """An index over a partial analysis would read as a verdict, not an outage."""
        strength = assess_evidence([
            DimensionReport("therapeutic_precedent",
                            EvidenceStatus.DATABASE_UNAVAILABLE, reason="snapshot gone"),
        ])
        result = assemble_index(_igg_components(), WHOLE_MAB, strength=strength)
        assert result.index is None
        assert any("infrastructural" in f for f in result.gating_failures)
        assert any("Retry" in f for f in result.gating_failures)

    def test_failed_assembly_validation_withholds_the_index(self):
        result = assemble_index(_igg_components(), WHOLE_MAB, assembly_validated=False)
        assert result.index is None
        assert result.assembly_validation_status == "failed"

    def test_insufficient_strength_grade_withholds_the_index(self):
        strength = assess_evidence([
            DimensionReport("a", EvidenceStatus.NOT_FOUND, reason="nothing found"),
        ])
        result = assemble_index(_igg_components(), WHOLE_MAB, strength=strength)
        assert result.index is None
        assert any("sufficiency" in f for f in result.gating_failures)

    def test_a_withheld_index_is_never_replaced_by_a_low_number(self):
        result = assemble_index([ev(ComponentSlot.VH, 90)], WHOLE_MAB)
        assert result.index is None
        assert result.gating_failures


class TestIndexArithmetic:
    def test_uniform_scores_give_that_score(self):
        assert assemble_index(_igg_components(90), WHOLE_MAB).index == pytest.approx(90.0)

    def test_a_weak_critical_component_is_not_hidden_by_strong_routine_ones(self):
        """§37 chooses a geometric mean for exactly this reason."""
        strong = assemble_index(_igg_components(95), WHOLE_MAB)
        with_weak_cdr = assemble_index(
            [ev(s, 20 if s is ComponentSlot.CDR_H3 else 95)
             for s in WEIGHT_SETS[FormatClass.IGG_LIKE]],
            WHOLE_MAB,
        )
        arithmetic_mean = (95 * 0.86 + 20 * 0.14)
        assert with_weak_cdr.index < strong.index
        # The geometric mean punishes the weak component harder than an average.
        assert with_weak_cdr.index < arithmetic_mean

    def test_lower_reliability_pulls_the_index_toward_better_evidenced_slots(self):
        confident = assemble_index(
            [ev(s, 50 if s is ComponentSlot.CDR_H3 else 95)
             for s in WEIGHT_SETS[FormatClass.IGG_LIKE]],
            WHOLE_MAB,
        )
        unreliable = assemble_index(
            [ev(s, 50, EvidenceGrade.D) if s is ComponentSlot.CDR_H3 else ev(s, 95)
             for s in WEIGHT_SETS[FormatClass.IGG_LIKE]],
            WHOLE_MAB,
        )
        assert unreliable.index > confident.index

    def test_zero_reliability_excludes_a_slot_rather_than_zeroing_the_index(self):
        result = assemble_index(
            [ev(s, 95) for s in WEIGHT_SETS[FormatClass.IGG_LIKE]
             if s is not ComponentSlot.CONSTANT_FC]
            + [ev(ComponentSlot.CONSTANT_FC, 95, EvidenceGrade.INSUFFICIENT)],
            WHOLE_MAB,
        )
        assert result.index == pytest.approx(95.0)
        assert ComponentSlot.CONSTANT_FC.value in result.components_missing
        assert result.evidence_coverage == pytest.approx(0.95)

    def test_epsilon_keeps_a_zero_scoring_component_finite(self):
        result = assemble_index(
            [ev(s, 0 if s is ComponentSlot.CDR_H3 else 90)
             for s in WEIGHT_SETS[FormatClass.IGG_LIKE]],
            WHOLE_MAB,
        )
        assert result.index is not None
        assert math.isfinite(result.index)
        assert result.index > 0.0

    def test_index_is_bounded_by_zero_and_one_hundred(self):
        for score in (0, 1, 50, 99, 100):
            result = assemble_index(_igg_components(score), WHOLE_MAB)
            assert 0.0 <= result.index <= 100.0

    def test_reliability_is_monotone_in_grade(self):
        order = [EvidenceGrade.A, EvidenceGrade.B, EvidenceGrade.C,
                 EvidenceGrade.D, EvidenceGrade.INSUFFICIENT]
        values = [RELIABILITY_BY_GRADE[g] for g in order]
        assert values == sorted(values, reverse=True)
        assert RELIABILITY_BY_GRADE[EvidenceGrade.INSUFFICIENT] == 0.0

    def test_epsilon_is_positive(self):
        assert 0.0 < SCORE_EPSILON < 1.0


class TestRequiredOutputsAndFraming:
    def test_all_six_additional_spec_outputs_are_present(self):
        result = assemble_index(
            _igg_components(), WHOLE_MAB,
            observed_pairing_status=EvidenceStatus.AVAILABLE,
            assembly_validated=True,
        )
        assert result.bottleneck_slot is not None
        assert result.bottleneck_score is not None
        assert result.evidence_coverage is not None
        assert result.observed_pairing_status == EvidenceStatus.AVAILABLE.value
        assert result.assembly_validation_status == "passed"
        assert result.missing_critical_evidence == {}

    def test_bottleneck_identifies_the_weakest_contributing_slot(self):
        result = assemble_index(
            [ev(s, 30 if s is ComponentSlot.OBSERVED_PAIRING else 95)
             for s in WEIGHT_SETS[FormatClass.IGG_LIKE]],
            WHOLE_MAB,
        )
        assert result.bottleneck_slot == ComponentSlot.OBSERVED_PAIRING.value
        assert result.bottleneck_score == 30

    def test_pairing_is_reported_separately_from_the_index(self):
        """§37 requires pairing evidence reported separately."""
        result = assemble_index(
            _igg_components(), WHOLE_MAB,
            observed_pairing_status=EvidenceStatus.NOT_FOUND,
        )
        assert result.observed_pairing_status == EvidenceStatus.NOT_FOUND.value

    def test_the_disclaimer_is_verbatim(self):
        """Required word-for-word by SPEC §37 and CLAUDE.md."""
        expected = (
            "This index summarizes comparative database evidence. It is not a "
            "probability of folding, expression, pairing, binding, stability, "
            "manufacturability, safety, efficacy, or experimental success."
        )
        assert INDEX_DISCLAIMER == expected
        assert assemble_index(_igg_components(), WHOLE_MAB).disclaimer == expected

    def test_the_disclaimer_is_present_even_when_the_index_is_withheld(self):
        result = assemble_index([ev(ComponentSlot.VH, 90)], WHOLE_MAB)
        assert result.disclaimer == INDEX_DISCLAIMER

    def test_versions_travel_with_every_index(self):
        result = assemble_index(_igg_components(), WHOLE_MAB)
        assert result.profile_version == PROFILE_VERSION
        assert result.weight_set_version == WEIGHT_SET_VERSION
        assert result.vocabulary_version

    def test_evidence_confidence_comes_from_the_strength_grade(self):
        strength = assess_evidence([
            DimensionReport("therapeutic_precedent", EvidenceStatus.AVAILABLE,
                            EvidenceClass.THERAPEUTIC_PRECEDENT,
                            independent_records=11, parent_antibodies=11,
                            parent_identifiers=frozenset({"a", "b"}),
                            structure_linked_records=3,
                            metadata_fields_present=9, metadata_fields_expected=10,
                            best_identity=0.99),
        ])
        result = assemble_index(_igg_components(), WHOLE_MAB, strength=strength)
        assert result.evidence_confidence == strength.grade.value

    def test_result_is_json_serializable(self):
        import json
        from dataclasses import asdict

        json.dumps(asdict(assemble_index(_igg_components(), WHOLE_MAB)), default=str)


class TestSingleDomainUsesItsOwnWeightSet:
    def test_vhh_does_not_borrow_the_igg_weights(self):
        result = assemble_index(
            [ev(s, 90) for s in WEIGHT_SETS[FormatClass.SINGLE_DOMAIN]], NANOBODY
        )
        assert result.format_class is FormatClass.SINGLE_DOMAIN
        assert result.index == pytest.approx(90.0)

    def test_vhh_is_not_penalised_for_having_no_light_chain(self):
        """Absent-by-architecture must not read as missing evidence."""
        result = assemble_index(
            [ev(s, 90) for s in WEIGHT_SETS[FormatClass.SINGLE_DOMAIN]], NANOBODY
        )
        assert result.evidence_coverage == pytest.approx(1.0)
        assert ComponentSlot.VL.value not in result.components_missing
