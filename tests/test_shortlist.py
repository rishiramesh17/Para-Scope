"""
Tests for the k-mer shortlist index and two-stage search (Version 2, Gate 3).

The shortlist is an optimization, so the properties that matter are: it must not
change scores, it must not leak comparisons across component types, it must be
deterministic, and it must be possible to turn off. Its *recall* is a measured
quantity, not a unit-testable one -- that lives in
scripts/benchmark_two_stage_search.py against the real corpus.
"""

import pytest

from antigens.indexes import ApproximateIndex
from antigens.indexes.shortlist import (
    DEFAULT_K,
    KmerShortlistIndex,
    extract_kmers,
)
from antigens.models import ComponentType

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
OTHER_VH = (
    "QVQLQQSGAELMKPGASVKISCKATGYTFSSYWIEWVKQRPGHGLEWIGEILPGSGSTYYNEKFKG"
    "KATFTADTSSNTAYMQLSSLTSEDSAVYYCARWDFDYWGQGTTLTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)


class TestKmerExtraction:
    def test_extracts_expected_count(self):
        seq = "ABCDEFGH"
        assert len(extract_kmers(seq, k=5)) == len(seq) - 5 + 1

    def test_sequence_shorter_than_k_yields_nothing(self):
        assert extract_kmers("ABC", k=5) == set()

    def test_kmers_are_distinct_not_a_multiset(self):
        """A repeated motif must not let one subject dominate by repetition alone."""
        assert extract_kmers("AAAAAAAAAA", k=5) == {"AAAAA"}


class TestKmerShortlistIndex:
    def test_identical_sequence_has_containment_one(self):
        index = KmerShortlistIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        candidates = index.shortlist(ComponentType.VH, TRASTUZUMAB_VH)
        assert candidates[0].record_id == "rec1"
        assert candidates[0].containment == 1.0

    def test_ranks_closer_sequence_first(self):
        index = KmerShortlistIndex()
        index.add(ComponentType.VH, "close", TRASTUZUMAB_VH)
        index.add(ComponentType.VH, "far", OTHER_VH)
        candidates = index.shortlist(ComponentType.VH, TRASTUZUMAB_VH)
        assert candidates[0].record_id == "close"
        assert candidates[0].containment > candidates[1].containment

    def test_component_types_never_share_a_posting_list(self):
        """A shortlist that mixed types would push an ineligible comparison into stage 2."""
        index = KmerShortlistIndex()
        index.add(ComponentType.VH, "heavy", TRASTUZUMAB_VH)
        index.add(ComponentType.VL_KAPPA, "light", TRASTUZUMAB_VL)

        heavy_hits = index.shortlist(ComponentType.VH, TRASTUZUMAB_VH)
        assert [c.record_id for c in heavy_hits] == ["heavy"]

        light_hits = index.shortlist(ComponentType.VL_KAPPA, TRASTUZUMAB_VL)
        assert [c.record_id for c in light_hits] == ["light"]

    def test_unrelated_query_yields_no_candidates(self):
        index = KmerShortlistIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        # A sequence over a disjoint set of 5-mers shares no k-mers at all.
        assert index.shortlist(ComponentType.VH, "WWWWWWWWWWWWWWWWWWWW") == []

    def test_respects_limit(self):
        index = KmerShortlistIndex()
        for i in range(20):
            index.add(ComponentType.VH, f"rec{i}", TRASTUZUMAB_VH)
        assert len(index.shortlist(ComponentType.VH, TRASTUZUMAB_VH, limit=5)) == 5

    def test_shortlist_is_deterministic(self):
        """Two identical submissions must return identical evidence."""
        index = KmerShortlistIndex()
        for i in range(30):
            index.add(ComponentType.VH, f"rec{i:02d}", TRASTUZUMAB_VH)
        first = [c.record_id for c in index.shortlist(ComponentType.VH, TRASTUZUMAB_VH, limit=10)]
        second = [c.record_id for c in index.shortlist(ComponentType.VH, TRASTUZUMAB_VH, limit=10)]
        assert first == second

    def test_sequence_shorter_than_k_is_skipped_not_an_error(self):
        index = KmerShortlistIndex(k=5)
        index.add(ComponentType.VH, "tiny", "ABC")
        assert len(index) == 0
        assert index.shortlist(ComponentType.VH, "ABC") == []

    def test_is_picklable(self):
        """Snapshots pickle this index; a lambda-factory defaultdict would break that."""
        import pickle

        index = KmerShortlistIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        restored = pickle.loads(pickle.dumps(index))
        assert restored.shortlist(ComponentType.VH, TRASTUZUMAB_VH)[0].record_id == "rec1"


class TestTwoStageSearch:
    def test_shortlist_does_not_change_scores(self):
        """The shortlist decides which alignments run, never what they score."""
        index = ApproximateIndex()
        index.add(ComponentType.VH, "close", TRASTUZUMAB_VH)
        index.add(ComponentType.VH, "far", OTHER_VH)

        staged = index.search(ComponentType.VH, TRASTUZUMAB_VH, top_k=5)
        brute = index.search(ComponentType.VH, TRASTUZUMAB_VH, top_k=5, shortlist_limit=None)

        assert [h.record_id for h in staged] == [h.record_id for h in brute]
        for a, b in zip(staged, brute):
            assert a.sequence_identity == b.sequence_identity
            assert a.global_identity == b.global_identity
            assert a.query_coverage == b.query_coverage

    def test_brute_force_mode_is_available(self):
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH, shortlist_limit=None)
        assert hits[0].global_identity == 1.0

    def test_search_defaults_to_exhaustive(self):
        """Completeness is the default; speed must be opted into explicitly."""
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        assert index.search(ComponentType.VH, TRASTUZUMAB_VH)[0].shortlisted is False

    def test_shortlisted_hits_are_labelled_as_such(self):
        """An approximate hit set must never be mistakable for a complete one."""
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH, shortlist_limit=200)
        assert hits and all(h.shortlisted for h in hits)

    def test_framework_request_falls_back_to_exhaustive(self):
        """
        Shortlisting frameworks is measurably worse and buys nothing (FR-H3
        recall@10 0.78 mean / 0.40 min, against 0.32s exhaustive), so a request
        is downgraded -- and `shortlisted` must report the truth: False.
        """
        index = ApproximateIndex()
        index.add(ComponentType.FR_H3, "rec1", "YYPDSVKGRFTISDDNAESTLYLQMSSLRSEDTAMYYC")
        hits = index.search(
            ComponentType.FR_H3, "YYPDSVKGRFTISDDNAESTLYLQMSSLRSEDTAMYYC", shortlist_limit=200
        )
        assert hits, "expected the exhaustive fallback to still return the match"
        assert hits[0].shortlisted is False

    def test_whole_domain_request_is_honoured(self):
        """The fallback must not swallow shortlisting where it is beneficial."""
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH, shortlist_limit=200)
        assert hits[0].shortlisted is True

    def test_unalignable_sequence_is_rejected_at_insertion(self):
        """IMGT gap placeholders crashed alignment at query time; reject them at add()."""
        index = ApproximateIndex()
        with pytest.raises(ValueError, match="substitution matrix"):
            index.add(ComponentType.FR_H3, "rec1", "YYPDSVK-GRFTISDDNAESTLYLQMSSLRSEDTAMYYC")

    def test_shortlist_index_stays_in_sync_with_sequences(self):
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        index.add(ComponentType.VL_KAPPA, "rec2", TRASTUZUMAB_VL)
        assert len(index) == 2
        assert len(index.shortlist_index) == 2

    def test_coverage_floor_still_applies_under_shortlisting(self):
        """Stage 1 must not become a way to smuggle a sub-coverage hit through."""
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        fragment = TRASTUZUMAB_VH[40:55]
        assert index.search(ComponentType.VH, fragment, top_k=5) == []
        assert index.search(ComponentType.VH, fragment, top_k=5, shortlist_limit=None) == []

    def test_short_component_still_rejected_from_index(self):
        index = ApproximateIndex()
        with pytest.raises(ValueError):
            index.add(ComponentType.CDR_H3, "rec1", "WGGDGFYAMDY")
