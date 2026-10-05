"""
Tests for empirical reference distributions and calibrated percentiles (Gate 9).

The contract under test is that a percentile always names the population it is
relative to, that insufficient data produces a status rather than a plausible
number, and that nothing silently substitutes a different stratum for the one
asked about.
"""

import json

import pytest

from para_scope.calibration import (
    CALIBRATION_VERSION,
    LENGTH_STRATIFIED_CLASSES,
    MINIMUM_SAMPLE_SIZE,
    CalibrationSet,
    build_distribution,
    length_band,
    load_calibration,
    save_calibration,
    stratum_key,
)
from para_scope.models import EvidenceStatus


def _dist(measurement="m", component_class="VH", stratum=None, values=None, n=100):
    values = values if values is not None else [i / n for i in range(n)]
    return build_distribution(
        measurement=measurement,
        component_class=component_class,
        stratum=stratum or component_class,
        stratification={},
        values=values,
        population_description="test population",
        source_snapshots=[{"source_name": "test", "raw_checksum": "abc"}],
    )


def _set(*distributions):
    return CalibrationSet(
        distributions={(d.measurement, d.stratum): d for d in distributions},
        built_at="2026-10-04T00:00:00+00:00",
        source_snapshots=({"source_name": "test", "raw_checksum": "abc"},),
    )


class TestLengthBanding:
    def test_bands_are_contiguous_and_total(self):
        for length in range(0, 60):
            assert length_band(length)

    def test_similar_lengths_share_a_band(self):
        assert length_band(12) == length_band(14)

    def test_distant_lengths_do_not(self):
        assert length_band(9) != length_band(22)

    def test_long_tail_has_an_open_band(self):
        assert length_band(200).endswith("+")

    def test_cdrs_are_length_stratified_and_domains_are_not(self):
        assert "CDR-H3" in LENGTH_STRATIFIED_CLASSES
        assert "VH" not in LENGTH_STRATIFIED_CLASSES
        assert stratum_key("CDR-H3", 13) == "CDR-H3|len=12-14"
        assert stratum_key("VH", None) == "VH"

    def test_length_stratified_class_requires_a_length(self):
        """Silently percentiling a 22-residue CDR-H3 against all lengths is the bug."""
        with pytest.raises(ValueError, match="length"):
            stratum_key("CDR-H3", None)


class TestPercentileComputation:
    def test_minimum_and_maximum_sit_at_the_edges(self):
        d = _dist(values=[float(i) for i in range(100)])
        assert d.percentile_of(0.0) < 1.0
        assert d.percentile_of(99.0) > 99.0

    def test_median_sits_near_fifty(self):
        d = _dist(values=[float(i) for i in range(101)])
        assert 49.0 <= d.percentile_of(50.0) <= 51.0

    def test_value_below_the_population_is_zero(self):
        d = _dist(values=[float(i) for i in range(10, 100)])
        assert d.percentile_of(0.0) == 0.0

    def test_value_above_the_population_is_one_hundred(self):
        d = _dist(values=[float(i) for i in range(100)])
        assert d.percentile_of(1000.0) == 100.0

    def test_ties_sit_at_the_centre_of_the_tied_block(self):
        """
        A saturated measurement is mostly ties. Taking either edge of the tied
        block would report a confidently extreme percentile for the single most
        ordinary value in the population.
        """
        d = _dist(values=[1.0] * 50 + [0.5] * 50)
        assert d.percentile_of(1.0) == 75.0
        assert d.percentile_of(0.5) == 25.0

    def test_summary_quantiles_are_ordered(self):
        s = _dist(values=[float(i) for i in range(1000)]).summary()
        ordered = [s["min"], s["p05"], s["p25"], s["median"], s["p75"], s["p95"], s["max"]]
        assert ordered == sorted(ordered)


class TestMissingAndInsufficientData:
    def test_absent_measurement_is_not_applicable_not_zero(self):
        """None means nothing was measured; it is not a percentile of zero."""
        result = _set(_dist()).percentile("m", "VH", None)
        assert result.status is EvidenceStatus.NOT_APPLICABLE
        assert result.percentile is None
        assert "not a percentile of zero" in result.reason

    def test_unknown_stratum_is_insufficient_evidence(self):
        result = _set(_dist()).percentile("m", "VL-kappa", 0.5)
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.percentile is None

    def test_unknown_stratum_does_not_fall_back_to_another(self):
        """Widening or substituting a stratum answers a different question."""
        result = _set(_dist(component_class="VH")).percentile("m", "VL-kappa", 0.5)
        assert result.percentile is None
        assert "not widened" in result.reason
        assert "fallback" in result.reason

    def test_small_stratum_is_insufficient_evidence_with_its_size(self):
        small = _dist(values=[0.1, 0.2, 0.3])
        result = _set(small).percentile("m", "VH", 0.2)
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.sample_size == 3
        assert str(MINIMUM_SAMPLE_SIZE) in result.reason

    def test_stratum_at_the_minimum_is_reportable(self):
        d = _dist(values=[i / MINIMUM_SAMPLE_SIZE for i in range(MINIMUM_SAMPLE_SIZE)])
        result = _set(d).percentile("m", "VH", 0.5)
        assert result.status is EvidenceStatus.AVAILABLE


class TestEveryPercentileNamesItsPopulation:
    def test_spec_required_fields_all_travel_with_the_number(self):
        """SPEC S19 Layer 2 requires all six, or the percentile is uninterpretable."""
        result = _set(_dist()).percentile("m", "VH", 0.5)
        assert result.status is EvidenceStatus.AVAILABLE
        assert result.percentile is not None
        assert result.component_class == "VH"        # component class
        assert result.stratum == "VH"                # reference distribution
        assert result.sample_size == 100             # sample size
        assert result.source_snapshots               # database snapshot
        assert result.stratification == {}           # stratification variables
        assert result.calibration_version == CALIBRATION_VERSION

    def test_population_is_described_and_not_claimed_to_be_repertoire(self):
        result = _set(_dist()).percentile("m", "VH", 0.5)
        assert result.population_description
        assert "NOT a sample" in result.disclaimer or "not a sample" in result.disclaimer

    def test_disclaimer_refuses_a_quality_reading(self):
        result = _set(_dist()).percentile("m", "VH", 0.5)
        for word in ("worse", "riskier", "less likely to succeed"):
            assert word in result.disclaimer

    def test_population_summary_accompanies_the_percentile(self):
        """A percentile alone cannot be sanity-checked; the population's shape can."""
        result = _set(_dist()).percentile("m", "VH", 0.5)
        assert set(result.population_summary) >= {"min", "median", "max"}


class TestArtifactRoundTrip:
    def test_save_and_load_preserves_distributions(self, tmp_path):
        path = tmp_path / "calib.json"
        d = _dist(values=[0.1, 0.2, 0.3, 0.4])
        save_calibration([d], "2026-10-04T00:00:00+00:00", [{"source_name": "test"}], path)
        loaded = load_calibration(path)
        assert len(loaded) == 1
        back = loaded.get("m", "VH")
        assert back is not None
        assert back.values == (0.1, 0.2, 0.3, 0.4)
        assert back.sample_size == 4

    def test_artifact_is_json_not_pickle(self, tmp_path):
        """Inspectable and diffable: a calibration is a published claim."""
        path = tmp_path / "calib.json"
        save_calibration([_dist()], "2026-10-04T00:00:00+00:00", [], path)
        payload = json.loads(path.read_text())
        assert payload["calibration_version"] == CALIBRATION_VERSION
        assert payload["minimum_sample_size"] == MINIMUM_SAMPLE_SIZE
        assert "population_disclaimer" in payload

    def test_loaded_percentile_matches_in_memory(self, tmp_path):
        path = tmp_path / "calib.json"
        d = _dist(values=[float(i) for i in range(100)])
        save_calibration([d], "2026-10-04T00:00:00+00:00", [], path)
        assert load_calibration(path).percentile("m", "VH", 42.0).percentile == d.percentile_of(42.0)
