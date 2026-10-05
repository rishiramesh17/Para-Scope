"""
Tests for the construct format vocabulary and format relevance (Gate 10.5).

The properties that matter: an axis the source text does not determine never
becomes agreement or disagreement; an uncomparable format is never a relevance of
zero; and the relevance number orders sensibly in format distance, which is the
part that was wrong twice during development.
"""

import pytest

from para_scope.construct_format import (
    AXIS_WEIGHTS,
    FORMAT_VOCABULARY_VERSION,
    MINIMUM_COMPARABLE_AXES,
    REQUIRED_AXIS,
    Architecture,
    AxisComparison,
    HostSpecies,
    Modification,
    Specificity,
    compute_format_relevance,
    parse_format,
)
from para_scope.models import EvidenceStatus


class TestAxisParsing:
    @pytest.mark.parametrize("raw,expected", [
        ("Whole mAb", Specificity.MONOSPECIFIC),
        ("Bispecific mAb", Specificity.BISPECIFIC),
        ("Bispecific Whole mAb ADC", Specificity.BISPECIFIC),
        ("Trispecific and Pentavalent (sc-5-VH); Bispecific entry", Specificity.TRISPECIFIC),
        ("Biparatopic Nanobody", Specificity.BISPECIFIC),
    ])
    def test_specificity(self, raw, expected):
        assert parse_format(raw).specificity is expected

    @pytest.mark.parametrize("raw,expected", [
        ("Whole mAb", Architecture.WHOLE_MAB),
        ("Whole  mAb", Architecture.WHOLE_MAB),          # double space in the real file
        ("Fab", Architecture.FAB),
        ("scFv", Architecture.SCFV),
        ("Nanobody", Architecture.SINGLE_DOMAIN),
        ("Bispecific Single Domains (VH-VH')", Architecture.SINGLE_DOMAIN),
        ("Fusion Protein", Architecture.FUSION_PROTEIN),
        ("Bispecific Mixed mAb and scFv", Architecture.MIXED),
        ("Bispecific Dual Variable Domain IG", Architecture.DUAL_VARIABLE_DOMAIN),
        ("Bispecific mAb", Architecture.MAB_UNSPECIFIED),
        ("Bispecific (VH-VK-VH'-VL Homodimer)", Architecture.UNKNOWN),
    ])
    def test_architecture(self, raw, expected):
        assert parse_format(raw).architecture is expected

    @pytest.mark.parametrize("raw,expected", [
        ("Whole mAb", Modification.NONE),
        ("Whole mAb ADC", Modification.ADC),
        ("Whole mAb Radiolabelled", Modification.RADIOLABELLED),
        ("Whole mAb Fusion", Modification.FUSION),
    ])
    def test_modification(self, raw, expected):
        assert parse_format(raw).modification is expected

    @pytest.mark.parametrize("raw,expected", [
        ("Canine Whole mAb", HostSpecies.CANINE),
        ("Feline Whole mAb", HostSpecies.FELINE),
        ("Whole mAb (Mouse)", HostSpecies.MOUSE),
        ("Whole mAb", HostSpecies.UNSTATED),
    ])
    def test_species(self, raw, expected):
        assert parse_format(raw).species is expected

    def test_absence_of_a_species_word_is_not_read_as_human(self):
        """Inferring human from silence would be a claim the column does not make."""
        assert parse_format("Whole mAb").species is HostSpecies.UNSTATED

    def test_inferred_axes_are_recorded_as_inferred(self):
        """Monospecific-from-absence is a convention of this column, not a fact."""
        parsed = parse_format("Whole mAb")
        assert "specificity" in parsed.inferred_axes
        assert "modification" in parsed.inferred_axes

    def test_stated_axes_are_not_marked_inferred(self):
        parsed = parse_format("Bispecific Whole mAb ADC")
        assert "specificity" not in parsed.inferred_axes
        assert "modification" not in parsed.inferred_axes

    def test_raw_text_is_preserved_verbatim(self):
        """The vocabulary is lossy by design; the source text must stay readable."""
        raw = "Bispecific Mixed Domains (L-kappa-H-gamma1_[VH-G1(CH1-h)]_L-kappa)"
        assert parse_format(raw).raw == raw

    def test_empty_and_none_parse_to_unknown_without_raising(self):
        for raw in (None, "", "   "):
            parsed = parse_format(raw)
            assert parsed.architecture is Architecture.UNKNOWN
            assert parsed.specificity is Specificity.UNKNOWN

    def test_vocabulary_version_is_stamped(self):
        assert parse_format("Whole mAb").vocabulary_version == FORMAT_VOCABULARY_VERSION


class TestUnknownNeverBecomesAgreementOrDisagreement:
    def test_unparseable_architecture_blocks_a_relevance_number(self):
        """Agreement on modifiers cannot establish format relevance."""
        result = compute_format_relevance(
            parse_format("Bispecific Mixed Domains (L-kappa-H-gamma1)"),
            parse_format("Bispecific Mixed Domains (VH-G1CH1-h)"),
        )
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.relevance is None
        assert "architecture" in result.reason

    def test_mab_unspecified_is_uncomparable_not_a_mismatch(self):
        """
        "Bispecific mAb" says the construct is an antibody and stops. It does not
        say the construct is not a whole mAb, so scoring it as a definite
        architecture mismatch answers a question the source never answered -- it
        produced 0.130 for a pair that may be architecturally identical.
        """
        result = compute_format_relevance(
            parse_format("Whole mAb"), parse_format("Bispecific mAb")
        )
        assert result.by_axis["architecture"] is AxisComparison.UNCOMPARABLE
        assert result.relevance is None
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE

    def test_unstated_species_is_uncomparable_on_both_sides(self):
        result = compute_format_relevance(
            parse_format("Whole mAb"), parse_format("Canine Whole mAb")
        )
        assert result.by_axis["species"] is AxisComparison.UNCOMPARABLE

    def test_uncomparable_axes_do_not_count_toward_the_total(self):
        result = compute_format_relevance(
            parse_format("Whole mAb"), parse_format("Whole mAb")
        )
        assert result.comparable_axes == 3        # species uncomparable both sides
        assert result.comparable_weight == pytest.approx(
            AXIS_WEIGHTS["architecture"] + AXIS_WEIGHTS["specificity"]
            + AXIS_WEIGHTS["modification"]
        )


class TestMissingFormatIsNotZeroRelevance:
    def test_no_declared_format_is_insufficient_not_zero(self):
        """Zero would assert the formats are maximally unrelated."""
        result = compute_format_relevance(None, parse_format("Whole mAb"))
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.relevance is None
        assert "0.10 format-relevance weight is unmet" in result.reason

    def test_no_matched_format_is_insufficient(self):
        result = compute_format_relevance(parse_format("Whole mAb"), None)
        assert result.relevance is None

    def test_the_reason_forbids_substituting_a_value(self):
        result = compute_format_relevance(None, parse_format("Whole mAb"))
        assert "not satisfied by a substituted value" in result.reason


class TestRelevanceOrdersByFormatDistance:
    def _rel(self, declared, matched):
        return compute_format_relevance(parse_format(declared), parse_format(matched)).relevance

    def test_identical_format_is_one(self):
        assert self._rel("Whole mAb", "Whole mAb") == 1.0

    def test_architecture_mismatch_scores_worse_than_two_modifier_mismatches(self):
        """
        The ordering bug this guards: at a first cut architecture carried 0.40
        against the modifiers' 0.60, so "Whole mAb" vs "scFv" (0.556) beat
        "Whole mAb" vs "Bispecific Whole mAb ADC" (0.444) -- disagreeing on the
        axis that decides what kind of molecule a construct is scored better than
        disagreeing on two modifiers.
        """
        architecture_differs = self._rel("Whole mAb", "scFv")
        modifiers_differ = self._rel("Whole mAb", "Bispecific Whole mAb ADC")
        assert architecture_differs < modifiers_differ

    def test_relevance_decreases_monotonically_with_added_differences(self):
        ordered = [
            self._rel("Whole mAb", "Whole mAb"),
            self._rel("Whole mAb", "Whole mAb ADC"),
            self._rel("Whole mAb", "Bispecific Whole mAb"),
            self._rel("Whole mAb", "Bispecific Whole mAb ADC"),
            self._rel("Whole mAb", "scFv"),
        ]
        assert ordered == sorted(ordered, reverse=True), ordered

    def test_architecture_weight_exceeds_the_others_combined(self):
        others = sum(v for k, v in AXIS_WEIGHTS.items() if k != REQUIRED_AXIS)
        assert AXIS_WEIGHTS[REQUIRED_AXIS] > others

    def test_relevance_is_bounded(self):
        for matched in ("Whole mAb", "scFv", "Whole mAb ADC", "Bispecific Whole mAb"):
            value = self._rel("Whole mAb", matched)
            assert 0.0 <= value <= 1.0


class TestFramingAndProvenance:
    def test_disclaimer_refuses_a_quality_reading(self):
        result = compute_format_relevance(parse_format("Whole mAb"), parse_format("scFv"))
        assert "not a judgement of either" in result.disclaimer
        for word in ("worse", "developable", "succeed"):
            assert word in result.disclaimer

    def test_both_sides_are_retained_for_audit(self):
        result = compute_format_relevance(parse_format("Whole mAb"), parse_format("scFv"))
        assert result.declared.raw == "Whole mAb"
        assert result.matched.raw == "scFv"

    def test_per_axis_breakdown_is_always_reported(self):
        result = compute_format_relevance(parse_format("Whole mAb"), parse_format("scFv"))
        assert set(result.by_axis) == set(AXIS_WEIGHTS)

    def test_minimum_comparable_axes_is_enforced(self):
        assert MINIMUM_COMPARABLE_AXES >= 2

    def test_result_is_json_serializable(self):
        import json
        from dataclasses import asdict

        result = compute_format_relevance(parse_format("Whole mAb"), parse_format("scFv"))
        json.dumps(asdict(result), default=str)
