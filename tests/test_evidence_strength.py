"""
Tests for evidence strength, confidence, and completeness (Gate 10).

Two properties carry most of the weight here:

  The four SPEC §28 fields stay separate. A test asserts no aggregate
  completeness field exists, because the moment one does, every caller uses it
  and the split becomes decorative.

  An infrastructural failure never becomes a biological grade. A failed lookup
  and a lookup that found nothing produce the same absence of evidence; grading
  them the same would turn an outage into a verdict about a candidate.
"""

from dataclasses import fields

import pytest

from para_scope.evidence_strength import (
    EVIDENCE_STRENGTH_VERSION,
    GRADE_A_MIN_IDENTITY,
    GRADE_B_MIN_IDENTITY,
    GRADE_C_MIN_IDENTITY,
    GRADE_THRESHOLD_VERSION,
    DimensionReport,
    EvidenceGrade,
    EvidenceStrengthReport,
    MetadataCompleteness,
    assess_evidence,
)
from para_scope.models import EvidenceClass, EvidenceStatus


def _strong(**kw):
    defaults = dict(
        name="therapeutic_precedent",
        status=EvidenceStatus.AVAILABLE,
        evidence_class=EvidenceClass.THERAPEUTIC_PRECEDENT,
        independent_records=11,
        parent_antibodies=11,
        structure_linked_records=16,
        therapeutic_records=11,
        metadata_fields_present=40,
        metadata_fields_expected=55,
        best_identity=1.0,
    )
    defaults.update(kw)
    return DimensionReport(**defaults)


class TestTheFourFieldsStaySeparate:
    def test_no_aggregate_completeness_field_exists(self):
        """
        SPEC §28: "Do not use one ambiguous completeness score."

        The moment an aggregate exists every caller uses it, and input
        completeness, analysis completion, metadata completeness, and evidence
        sufficiency -- four different situations with four different remedies --
        become indistinguishable.
        """
        names = {f.name for f in fields(EvidenceStrengthReport)}
        for forbidden in (
            "overall_completeness", "completeness", "completeness_score",
            "total_completeness", "combined_completeness", "confidence_score",
        ):
            assert forbidden not in names

    def test_all_four_fields_are_present(self):
        report = assess_evidence([_strong()])
        assert report.input_completeness is not None
        assert report.analysis_completion is not None
        assert report.metadata_completeness is not None
        assert report.evidence_sufficiency is not None

    def test_the_four_fields_are_independently_populated(self):
        """A single number could not express this combination."""
        report = assess_evidence(
            [_strong(metadata_fields_present=0, metadata_fields_expected=55)],
            expected_components={"VH": "EVQ", "VL-kappa": None},
        )
        assert report.input_completeness.is_complete is False      # input missing
        assert report.analysis_completion.is_complete is True      # modules fine
        assert report.metadata_completeness.fraction == 0.0        # records bare
        assert report.evidence_sufficiency.independent_records == 11  # evidence rich


class TestInputCompleteness:
    def test_supplied_and_missing_are_both_listed(self):
        report = assess_evidence([_strong()], expected_components={"VH": "EVQ", "VL-kappa": None})
        assert report.input_completeness.supplied == ["VH"]
        assert "VL-kappa" in report.input_completeness.missing

    def test_not_applicable_is_not_counted_as_missing(self):
        """A VHH has no light chain by design; calling that missing input would
        describe an architecture as a deficiency."""
        report = assess_evidence(
            [_strong()],
            expected_components={"VHH": "QVQ", "VL-kappa": None},
            not_applicable_components={"VL-kappa": "single-domain construct has no light chain"},
        )
        assert report.input_completeness.missing == {}
        assert report.input_completeness.is_complete is True
        assert "VL-kappa" in report.input_completeness.not_applicable
        assert report.input_completeness.required_count == 1

    def test_incomplete_input_is_recorded_as_a_limitation_not_as_evidence(self):
        report = assess_evidence([_strong()], expected_components={"VH": "E", "VL-kappa": None})
        assert any("not evidence about the candidate" in lim for lim in report.limitations)


class TestAnalysisCompletionIsNotEvidence:
    @pytest.mark.parametrize("status", [
        EvidenceStatus.DATABASE_UNAVAILABLE,
        EvidenceStatus.ANALYSIS_FAILED,
        EvidenceStatus.ACCESS_RESTRICTED,
    ])
    def test_infrastructural_failure_blocks_grading_entirely(self, status):
        """SPEC §28: an outage must not be treated as evidence against the candidate."""
        report = assess_evidence([DimensionReport("therapeutic_precedent", status, reason="x")])
        assert report.grade is EvidenceGrade.INSUFFICIENT
        assert report.grade is not EvidenceGrade.D
        assert report.grade_blocked_reason is not None
        assert "Retry" in report.grade_blocked_reason

    def test_an_outage_does_not_downgrade_otherwise_strong_evidence_to_d(self):
        """The failure mode: strong evidence plus one outage scoring worse than weak."""
        report = assess_evidence([
            _strong(),
            DimensionReport("pairing_precedent", EvidenceStatus.DATABASE_UNAVAILABLE, reason="x"),
        ])
        assert report.grade is EvidenceGrade.INSUFFICIENT
        assert report.grade is not EvidenceGrade.D

    def test_outage_is_listed_separately_from_ordinary_failure(self):
        report = assess_evidence([
            DimensionReport("a", EvidenceStatus.DATABASE_UNAVAILABLE, reason="snapshot gone"),
        ])
        assert "a" in report.analysis_completion.infrastructural_failures
        assert report.analysis_completion.has_infrastructural_failure is True

    @pytest.mark.parametrize("status", [
        EvidenceStatus.NOT_FOUND,
        EvidenceStatus.NOT_APPLICABLE,
        EvidenceStatus.INSUFFICIENT_EVIDENCE,
        EvidenceStatus.UNSUPPORTED_COMPONENT,
    ])
    def test_a_module_that_ran_and_found_nothing_completed_successfully(self, status):
        """Counting NOT_FOUND as a failure would turn an answer into an outage."""
        report = assess_evidence([DimensionReport("a", status, reason="looked, found nothing")])
        assert report.analysis_completion.completed == ["a"]
        assert report.analysis_completion.failed == {}
        assert report.analysis_completion.is_complete is True


class TestMetadataCompleteness:
    def test_fraction_is_none_when_nothing_was_annotatable(self):
        """0.0 would say the records are unannotated; 1.0 that they are perfect.
        Both are claims about records that do not exist."""
        assert MetadataCompleteness().fraction is None

    def test_zero_present_of_many_expected_is_zero_not_none(self):
        report = assess_evidence([_strong(metadata_fields_present=0, metadata_fields_expected=10)])
        assert report.metadata_completeness.fraction == 0.0

    def test_completeness_is_attributable_per_dimension(self):
        report = assess_evidence([
            _strong(name="therapeutic_precedent", metadata_fields_present=9, metadata_fields_expected=10),
            _strong(name="pairing_precedent", metadata_fields_present=1, metadata_fields_expected=10),
        ])
        by_dimension = report.metadata_completeness.by_dimension
        assert by_dimension["therapeutic_precedent"] == (9, 10)
        assert by_dimension["pairing_precedent"] == (1, 10)

    def test_undefined_metadata_is_recorded_as_a_limitation(self):
        report = assess_evidence([
            DimensionReport("germline", EvidenceStatus.AVAILABLE, best_identity=0.9),
        ])
        assert any("undefined rather than zero" in lim for lim in report.limitations)


class TestEvidenceSufficiency:
    def test_counts_are_summed_across_dimensions(self):
        report = assess_evidence([
            _strong(name="a", independent_records=5, parent_antibodies=3),
            _strong(name="b", independent_records=7, parent_antibodies=4),
        ])
        assert report.evidence_sufficiency.independent_records == 12
        assert report.evidence_sufficiency.parent_antibodies == 7

    def test_evidence_classes_are_reported_separately_never_summed(self):
        """SPEC §17: evidence classes must be reported separately."""
        report = assess_evidence([
            _strong(name="a", evidence_class=EvidenceClass.THERAPEUTIC_PRECEDENT, independent_records=5),
            _strong(name="b", evidence_class=EvidenceClass.DIRECT_PAIRING, independent_records=7),
        ])
        by_class = report.evidence_sufficiency.by_evidence_class
        assert by_class["therapeutic_precedent"] == 5
        assert by_class["direct_paired_chain"] == 7

    def test_blocked_terms_surface_as_missing_evidence(self):
        report = assess_evidence([
            _strong(blocked_terms={"format_relevance": "needs a declared format"}),
        ])
        missing = report.evidence_sufficiency.missing_evidence
        assert "therapeutic_precedent.format_relevance" in missing

    def test_unavailable_layer4_fields_say_why_rather_than_being_omitted(self):
        """study_count, subject_count, and a bootstrap interval are not derivable
        from the current sources. Silently omitting them hides a known gap."""
        sufficiency = assess_evidence([_strong()]).evidence_sufficiency
        assert sufficiency.study_count is None
        assert sufficiency.subject_count is None
        assert sufficiency.bootstrap_interval is None
        for field_name in ("study_count", "subject_count", "bootstrap_interval"):
            assert field_name in sufficiency.unavailable_layer4_fields
            assert sufficiency.unavailable_layer4_fields[field_name]

    def test_a_reference_lookup_contributes_no_independent_records(self):
        """A germline allele is one curated reference sequence, not database
        precedent; counting it would let a single lookup inflate sufficiency."""
        report = assess_evidence([
            DimensionReport("germline", EvidenceStatus.AVAILABLE,
                            EvidenceClass.GERMLINE, best_identity=0.91),
        ])
        assert report.evidence_sufficiency.independent_records == 0


class TestGradeBands:
    def test_strong_repeated_experimental_evidence_is_a(self):
        assert assess_evidence([_strong()]).grade is EvidenceGrade.A

    def test_repeated_but_not_experimental_is_b(self):
        report = assess_evidence([_strong(structure_linked_records=0)])
        assert report.grade is EvidenceGrade.B
        assert "structure-linked" in report.grade_rationale

    def test_single_parent_antibody_cannot_reach_a_or_b(self):
        """'Independently repeated' means more than one parent, not one record twice."""
        report = assess_evidence([_strong(parent_antibodies=1)])
        assert report.grade is EvidenceGrade.C

    def test_moderate_identity_is_c(self):
        report = assess_evidence([_strong(best_identity=0.75, parent_antibodies=1,
                                          structure_linked_records=0)])
        assert report.grade is EvidenceGrade.C

    def test_sparse_and_poorly_annotated_is_d_even_at_good_identity(self):
        """
        SPEC §19 defines D as "weak, indirect, predicted, sparse, or poorly
        annotated" -- not merely low identity.

        This combination does not occur in the current corpus (measured:
        metadata fraction never drops below 0.80 and identity never below 0.703,
        so D is unreachable against SAbDab2 + Thera-SAbDab -- see
        reports/GATE10_EVIDENCE_STRENGTH.md). The band is kept and tested anyway:
        it is spec-mandated, it would fire on a sparser or less curated source,
        and an untested band is dead code that silently stops working.
        """
        report = assess_evidence([
            _strong(
                best_identity=0.92,
                parent_antibodies=1,
                metadata_fields_present=1,
                metadata_fields_expected=20,
            )
        ])
        assert report.grade is EvidenceGrade.D
        assert "Sparse and poorly annotated" in report.grade_rationale

    def test_sparse_but_well_annotated_is_c_not_d(self):
        """A single, fully annotated exact match is not weak evidence."""
        report = assess_evidence([
            _strong(best_identity=0.92, parent_antibodies=1,
                    metadata_fields_present=20, metadata_fields_expected=20),
        ])
        assert report.grade is EvidenceGrade.C

    def test_below_the_precedent_floor_is_d(self):
        report = assess_evidence([_strong(best_identity=0.50)])
        assert report.grade is EvidenceGrade.D

    def test_no_records_is_insufficient_not_d(self):
        report = assess_evidence([
            DimensionReport("a", EvidenceStatus.NOT_FOUND, reason="no match"),
        ])
        assert report.grade is EvidenceGrade.INSUFFICIENT

    def test_records_without_a_similarity_measurement_cannot_be_graded(self):
        report = assess_evidence([
            DimensionReport("a", EvidenceStatus.AVAILABLE, independent_records=3),
        ])
        assert report.grade is EvidenceGrade.INSUFFICIENT
        assert "none carries a similarity measurement" in report.grade_rationale

    def test_bands_are_ordered_and_do_not_overlap(self):
        assert GRADE_C_MIN_IDENTITY < GRADE_B_MIN_IDENTITY < GRADE_A_MIN_IDENTITY

    def test_every_grade_carries_a_rationale(self):
        for dimension in (
            _strong(),
            _strong(structure_linked_records=0),
            _strong(best_identity=0.75, parent_antibodies=1),
            _strong(best_identity=0.5),
        ):
            report = assess_evidence([dimension])
            assert report.grade_rationale.strip()


class TestFramingAndVersioning:
    def test_grade_disclaimer_refuses_a_candidate_reading(self):
        report = assess_evidence([_strong()])
        assert "not a statement that" in report.disclaimer
        for word in ("bind", "stable", "safe", "succeed"):
            assert word in report.disclaimer

    def test_grade_thresholds_are_versioned_separately_from_the_engine(self):
        """A threshold change alters every stored grade without altering any
        measurement, so the two versions must be independently traceable."""
        report = assess_evidence([_strong()])
        assert report.engine_version == EVIDENCE_STRENGTH_VERSION
        assert report.grade_threshold_version == GRADE_THRESHOLD_VERSION
        assert report.grade_threshold_version != report.engine_version

    def test_threshold_version_admits_it_is_provisional(self):
        """SPEC §19 requires benchmarked thresholds. These are not, and say so."""
        assert "provisional" in GRADE_THRESHOLD_VERSION
        report = assess_evidence([_strong()])
        assert any("provisional" in lim for lim in report.limitations)

    def test_report_is_json_serializable(self):
        import json
        from dataclasses import asdict

        json.dumps(asdict(assess_evidence([_strong()])), default=str)
