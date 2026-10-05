"""
Gate 9: calibration measured on the real corpus, pinned to the artifact it was built from.

A percentile is only meaningful relative to a population, so these tests assert
the *shape* of the real populations as well as the mechanics. A unit fixture
cannot catch a calibration built from the wrong population: 120 synthetic values
percentile correctly while describing nothing.
"""

from pathlib import Path

import pytest

from para_scope.calibration import MINIMUM_SAMPLE_SIZE, load_calibration
from para_scope.models import ComponentType, EvidenceStatus
from para_scope.persistence import load_corpus_snapshot
from para_scope.therapeutic_precedent import (
    SIMILARITY_MEASUREMENT,
    compute_therapeutic_precedent,
)

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
CALIBRATION_PATH = Path("data/cache/calibration_v1.json")

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)

pytestmark = pytest.mark.skipif(
    not (SNAPSHOT_PATH.exists() and CALIBRATION_PATH.exists()),
    reason=(
        "Needs the combined corpus snapshot and calibration artifact. Run "
        "scripts/build_combined_corpus_cache.py, scripts/build_calibration.py, "
        "then scripts/build_similarity_calibration.py."
    ),
)


@pytest.fixture(scope="module")
def calibration():
    return load_calibration(CALIBRATION_PATH)


@pytest.fixture(scope="module")
def corpus():
    loaded, _stats, _manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    return loaded


class TestArtifactShape:
    def test_every_distribution_names_its_population_and_snapshots(self, calibration):
        for (measurement, stratum) in calibration._distributions:
            d = calibration.get(measurement, stratum)
            assert d.population_description, f"{measurement}/{stratum}"
            assert d.source_snapshots, f"{measurement}/{stratum}"
            assert d.sample_size == len(d.values)

    def test_no_distribution_claims_to_be_natural_repertoire(self, calibration):
        """The population is crystallised/developed antibodies. Saying otherwise is the error."""
        for (measurement, stratum) in calibration._distributions:
            text = calibration.get(measurement, stratum).population_description
            assert "NOT a sample of natural antibody repertoire" in text

    def test_expected_measurement_families_are_present(self, calibration):
        measurements = {m for (m, _s) in calibration._distributions}
        families = {m.split(".")[0] for m in measurements}
        assert {"feature", "motif", "window_coverage", "germline", "similarity"} <= families

    def test_values_are_sorted_so_percentiles_are_valid(self, calibration):
        for (measurement, stratum) in calibration._distributions:
            values = calibration.get(measurement, stratum).values
            assert list(values) == sorted(values), f"{measurement}/{stratum} unsorted"


class TestGermlineDistributionIsBiologicallyPlausible:
    def test_germline_identity_sits_in_the_somatic_mutation_range(self, calibration):
        """
        Antibodies carry somatic mutations, so identity to nearest germline V gene
        should be high but well below 1.0. A population pinned at 1.0 would mean
        germline assignment was returning self-matches; near 0.3 would mean it was
        aligning against the wrong locus.
        """
        for component_class in ("VH", "VL-kappa", "VL-lambda"):
            d = calibration.get("germline.nearest_identity", component_class)
            assert d is not None, component_class
            assert d.sample_size >= MINIMUM_SAMPLE_SIZE
            summary = d.summary()
            assert 0.55 <= summary["median"] <= 0.95, (component_class, summary)
            assert summary["max"] < 1.0, "a germline match should not be identical"
            assert summary["min"] > 0.30, (component_class, summary)


class TestSimilarityDistributionIsTheRightPopulation:
    def test_self_matches_are_excluded(self, calibration):
        """
        Every sampled component is in the corpus, and some *are* therapeutics. If
        self-matches leaked in, the population would pile up at 1.0 and no genuine
        near-match could ever look unusual.
        """
        d = calibration.get(SIMILARITY_MEASUREMENT, "VH")
        assert d is not None
        at_one = sum(1 for v in d.values if v >= 0.999)
        assert at_one / d.sample_size < 0.10, (
            f"{at_one}/{d.sample_size} values at 1.0 suggests self-matches leaked in"
        )

    def test_population_median_is_below_an_exact_match(self, calibration):
        d = calibration.get(SIMILARITY_MEASUREMENT, "VH")
        assert d.summary()["median"] < 0.95


class TestCalibrationChangesWhatGate7CanSay:
    def test_term_is_blocked_without_calibration_and_cleared_with_it(self, corpus, calibration):
        without = compute_therapeutic_precedent(
            corpus, ComponentType.VH, TRASTUZUMAB_VH, top_k=3, shortlist_limit=500,
        )
        assert "calibrated_component_similarity" in without.blocked_terms
        assert without.matches[0].calibrated_similarity is None

        with_calib = compute_therapeutic_precedent(
            corpus, ComponentType.VH, TRASTUZUMAB_VH, top_k=3, shortlist_limit=500,
            calibration=calibration,
        )
        assert "calibrated_component_similarity" not in with_calib.blocked_terms
        assert with_calib.matches[0].calibrated_similarity is not None

    def test_an_exact_therapeutic_match_is_an_extreme_percentile(self, corpus, calibration):
        """The point of Layer 2: 1.000 identity means nothing until placed."""
        result = compute_therapeutic_precedent(
            corpus, ComponentType.VH, TRASTUZUMAB_VH, top_k=3, shortlist_limit=500,
            calibration=calibration,
        )
        top = result.matches[0]
        assert top.global_identity == 1.0
        percentile = top.calibrated_similarity
        assert percentile.status is EvidenceStatus.AVAILABLE
        assert percentile.percentile >= 99.0
        # Required by SPEC S19 -- a bare percentile is uninterpretable.
        assert percentile.sample_size >= MINIMUM_SAMPLE_SIZE
        assert percentile.population_description
        assert percentile.population_summary["median"] < percentile.raw_value

    def test_raw_measurement_is_unchanged_by_calibration(self, corpus, calibration):
        """Calibration adds a layer; it must not alter Layer 1."""
        a = compute_therapeutic_precedent(
            corpus, ComponentType.VH, TRASTUZUMAB_VH, top_k=5, shortlist_limit=500,
        )
        b = compute_therapeutic_precedent(
            corpus, ComponentType.VH, TRASTUZUMAB_VH, top_k=5, shortlist_limit=500,
            calibration=calibration,
        )
        assert [m.global_identity for m in a.matches] == [m.global_identity for m in b.matches]
        assert [m.therapeutic_name for m in a.matches] == [m.therapeutic_name for m in b.matches]
