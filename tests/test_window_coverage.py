"""
Tests for human repertoire window coverage (Version 2, Gate 6).

Three things are under test:

1. The arithmetic -- windows are chopped correctly and coverage equals the
   reference implementation's own score.
2. The failure modes -- every path that cannot produce a number returns an
   explicit status with a reason and leaves `coverage` unset. This is the
   `missing evidence is not zero evidence` constraint, and it is the reason most
   of this file exists.
3. The documented limitation -- a camelid VHH scores near a humanized antibody.
   That fact is what forbids naming this "humanness", so it is asserted rather
   than left in a comment where it could quietly stop being true.
"""

import pytest

from para_scope.models import EvidenceStatus
from para_scope.window_coverage import (
    COVERAGE_DISCLAIMER,
    MINIMUM_WINDOWS_FOR_FRACTION,
    REFERENCE_SET_ID,
    WINDOW_COVERAGE_VERSION,
    WINDOW_LENGTH,
    ReferenceSetUnavailable,
    chop_windows,
    compute_window_coverage,
    load_reference_set,
)

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
MUROMONAB_VH = (
    "QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKD"
    "KATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS"
)
CAPLACIZUMAB_VHH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEG"
    "RFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS"
)


@pytest.fixture(scope="module")
def db():
    """The reference set loads once (~3s, ~5.7M peptides) and is reused."""
    return load_reference_set()


class _RaisingDatabase:
    """Stands in for a reference set whose lookups fail mid-analysis."""

    def contains(self, peptide):
        raise RuntimeError("lookup exploded")


class TestWindowChopping:
    def test_windows_are_overlapping_and_ordered(self):
        assert chop_windows("ABCDEFGHIJK") == ["ABCDEFGHI", "BCDEFGHIJ", "CDEFGHIJK"]

    def test_window_count_is_length_minus_window_plus_one(self):
        assert len(chop_windows(TRASTUZUMAB_VH)) == len(TRASTUZUMAB_VH) - WINDOW_LENGTH + 1

    def test_sequence_shorter_than_window_yields_no_windows(self):
        assert chop_windows("EVQLVE") == []

    def test_sequence_exactly_window_length_yields_one_window(self):
        assert chop_windows("EVQLVESGG") == ["EVQLVESGG"]


class TestCoverageArithmetic:
    def test_matches_reference_implementation(self, db):
        """Our per-window count must equal promb's own aggregate score."""
        result = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert result.coverage == pytest.approx(db.compute_peptide_content(TRASTUZUMAB_VH))

    def test_coverage_is_observed_over_total(self, db):
        result = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert result.coverage == pytest.approx(result.windows_observed / result.windows_total)

    def test_counts_are_consistent(self, db):
        result = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert result.windows_total == len(TRASTUZUMAB_VH) - WINDOW_LENGTH + 1
        assert result.windows_observed + result.windows_not_observed == result.windows_total

    def test_observations_align_with_counts(self, db):
        result = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert len(result.observations) == result.windows_total
        assert sum(1 for o in result.observations if o.observed) == result.windows_observed
        assert [o.sequence_offset for o in result.observations] == list(range(result.windows_total))
        assert result.observations[0].window_index == 1

    def test_observations_can_be_suppressed_without_changing_the_value(self, db):
        full = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        lean = compute_window_coverage(TRASTUZUMAB_VH, database=db, record_observations=False)
        assert lean.observations == []
        assert lean.coverage == full.coverage
        assert lean.windows_observed == full.windows_observed

    def test_is_deterministic(self, db):
        first = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        second = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert first.coverage == second.coverage

    def test_whitespace_and_case_are_normalised(self, db):
        messy = compute_window_coverage(f"  {TRASTUZUMAB_VH.lower()}  ", database=db)
        clean = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert messy.coverage == clean.coverage


class TestDocumentedLimitation:
    """
    The naming constraint depends on these orderings. If they change, the module
    docstring's justification is stale and must be revisited before release.
    """

    def test_humanized_scores_above_mouse(self, db):
        humanized = compute_window_coverage(TRASTUZUMAB_VH, database=db).coverage
        mouse = compute_window_coverage(MUROMONAB_VH, database=db).coverage
        assert humanized > mouse

    def test_camelid_vhh_scores_above_mouse(self, db):
        """The load-bearing result: a llama VHH is not separated as non-human."""
        vhh = compute_window_coverage(CAPLACIZUMAB_VHH, database=db).coverage
        mouse = compute_window_coverage(MUROMONAB_VH, database=db).coverage
        assert vhh > mouse

    def test_camelid_vhh_is_close_to_humanized(self, db):
        """Within 0.15 of humanized -- far too close to read as a humanness call."""
        vhh = compute_window_coverage(CAPLACIZUMAB_VHH, database=db).coverage
        humanized = compute_window_coverage(TRASTUZUMAB_VH, database=db).coverage
        assert abs(humanized - vhh) < 0.15

    def test_recorded_validation_figures_still_hold(self, db):
        """Pins the figures quoted in the module docstring and the Gate 6 report."""
        assert compute_window_coverage(TRASTUZUMAB_VH, database=db).coverage == pytest.approx(
            0.7054, abs=0.001
        )
        assert compute_window_coverage(TRASTUZUMAB_VL, database=db).coverage == pytest.approx(
            0.8283, abs=0.001
        )
        assert compute_window_coverage(MUROMONAB_VH, database=db).coverage == pytest.approx(
            0.5225, abs=0.001
        )
        assert compute_window_coverage(CAPLACIZUMAB_VHH, database=db).coverage == pytest.approx(
            0.6333, abs=0.001
        )


class TestMissingEvidenceIsNotZero:
    """Every non-AVAILABLE path must leave coverage unset and say why."""

    def test_short_sequence_is_not_applicable_not_zero(self, db):
        result = compute_window_coverage("EVQLVE", database=db)
        assert result.status is EvidenceStatus.NOT_APPLICABLE
        assert result.coverage is None
        assert result.reason

    def test_cdr_length_sequence_does_not_raise(self, db):
        """promb raises ZeroDivisionError here; a CDR must not crash a profile."""
        for cdr in ("AAA", "GGSA", "QQHYTTPPT"[:5]):
            result = compute_window_coverage(cdr, database=db)
            assert result.coverage is None
            assert result.status in {
                EvidenceStatus.NOT_APPLICABLE,
                EvidenceStatus.INSUFFICIENT_EVIDENCE,
            }

    def test_few_windows_is_insufficient_evidence_with_counts_retained(self, db):
        short = TRASTUZUMAB_VH[:20]  # 12 windows, below the minimum
        result = compute_window_coverage(short, database=db)
        assert result.status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.coverage is None
        assert result.windows_total == len(short) - WINDOW_LENGTH + 1
        assert result.windows_total < MINIMUM_WINDOWS_FOR_FRACTION
        # Raw counts survive so Gate 9 can build a length-matched null from them.
        assert result.windows_observed >= 0
        assert result.observations

    def test_non_standard_residue_is_unsupported_not_scored(self, db):
        """An X must not deflate coverage into looking like a low score."""
        with_x = TRASTUZUMAB_VH[:60] + "X" + TRASTUZUMAB_VH[61:]
        result = compute_window_coverage(with_x, database=db)
        assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT
        assert result.coverage is None
        assert "X" in result.reason

    def test_ambiguity_codes_are_all_rejected(self, db):
        for residue in ("X", "B", "Z", "J", "U", "*"):
            seq = TRASTUZUMAB_VH[:60] + residue + TRASTUZUMAB_VH[61:]
            result = compute_window_coverage(seq, database=db)
            assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT
            assert result.coverage is None

    def test_input_is_never_modified_to_make_it_scoreable(self, db):
        """The X path must not fall back to scoring a cleaned sequence."""
        with_x = TRASTUZUMAB_VH[:60] + "X" + TRASTUZUMAB_VH[61:]
        stripped = with_x.replace("X", "")
        assert compute_window_coverage(with_x, database=db).coverage is None
        assert compute_window_coverage(stripped, database=db).coverage is not None

    def test_empty_sequence_is_unsupported(self, db):
        for empty in ("", "   "):
            result = compute_window_coverage(empty, database=db)
            assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT
            assert result.coverage is None

    def test_database_unavailable_does_not_lower_the_score(self, monkeypatch):
        def explode(*args, **kwargs):
            raise ReferenceSetUnavailable("promb not installed")

        monkeypatch.setattr("para_scope.window_coverage.load_reference_set", explode)
        result = compute_window_coverage(TRASTUZUMAB_VH)
        assert result.status is EvidenceStatus.DATABASE_UNAVAILABLE
        assert result.coverage is None
        assert "not installed" in result.reason
        # Windows were still countable; only the lookup was impossible.
        assert result.windows_total > 0

    def test_lookup_failure_is_analysis_failed(self):
        result = compute_window_coverage(TRASTUZUMAB_VH, database=_RaisingDatabase())
        assert result.status is EvidenceStatus.ANALYSIS_FAILED
        assert result.coverage is None
        assert "lookup exploded" in result.reason

    def test_no_failure_path_returns_a_neutral_value(self, db):
        """Nothing may substitute 0.0, 0.5, or any placeholder for an absent value."""
        failing = [
            compute_window_coverage("EVQLVE", database=db),
            compute_window_coverage("", database=db),
            compute_window_coverage(TRASTUZUMAB_VH[:20], database=db),
            compute_window_coverage("EVQLVESGGGXVQPGG", database=db),
            compute_window_coverage(TRASTUZUMAB_VH, database=_RaisingDatabase()),
        ]
        for result in failing:
            assert result.status is not EvidenceStatus.AVAILABLE
            assert result.coverage is None
            assert result.reason, f"{result.status} carries no reason"


class TestProvenanceAndFraming:
    def test_every_result_carries_reproduction_metadata(self, db):
        result = compute_window_coverage(TRASTUZUMAB_VH, database=db)
        assert result.engine_version == WINDOW_COVERAGE_VERSION
        assert result.window_length == WINDOW_LENGTH
        assert result.reference_set_id == REFERENCE_SET_ID
        assert result.reference_set_file
        assert result.reference_set_description
        assert result.reference_library_version

    def test_failure_results_also_carry_provenance(self, db):
        result = compute_window_coverage("EVQLVE", database=db)
        assert result.engine_version == WINDOW_COVERAGE_VERSION
        assert result.reference_set_id == REFERENCE_SET_ID

    def test_disclaimer_is_present_on_every_result(self, db):
        for seq in (TRASTUZUMAB_VH, "EVQLVE", ""):
            assert compute_window_coverage(seq, database=db).disclaimer == COVERAGE_DISCLAIMER

    def test_disclaimer_refuses_the_prohibited_readings(self):
        lowered = COVERAGE_DISCLAIMER.lower()
        for claim in ("humanness", "immunogenicity", "safety"):
            assert claim in lowered
        assert "not a measure of" in lowered

    def test_module_exposes_no_humanness_vocabulary(self):
        """The naming decision is enforced, not just documented."""
        import para_scope.window_coverage as module

        public = [name for name in dir(module) if not name.startswith("_")]
        for name in public:
            lowered = name.lower()
            assert "humanness" not in lowered
            assert "naturalness" not in lowered
            assert "immunogenic" not in lowered
