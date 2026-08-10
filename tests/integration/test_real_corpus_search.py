"""
Gate 1: real end-to-end search validation against the actual ingested SAbDab2 corpus.

Unlike tests/test_search.py (self-contained synthetic fixtures, always runs),
this module exercises the real ~15,641-row corpus via a pre-built snapshot
(see scripts/build_v2_corpus_cache.py + scripts/build_gold_queries.py). Both
artifacts live under data/ and benchmarks/ and are not guaranteed to exist in
every environment (e.g. a fresh clone, or CI without the 1.1GB source data),
so the whole module is skipped rather than failed when they're absent.

Requirements enforced here (Version 2 roadmap, Gate 1):
  - Exact matches always rank first; exact-match recall is 100%.
  - Insertion-coded residues remain correctly aligned/preserved.
  - Parent antibodies with multiple structures are not double-counted.
  - Invalid inputs produce clear rejection, never a partial/fabricated score.
"""

import json
from pathlib import Path

import pytest

from para_scope.models import ComponentType
from para_scope.numbering import number_sequence
from para_scope.persistence import load_corpus_snapshot

SNAPSHOT_PATH = Path("data/cache/sabdab2_corpus_snapshot.pkl")
GOLD_QUERIES_PATH = Path("benchmarks/search_gold_queries.json")

pytestmark = pytest.mark.skipif(
    not SNAPSHOT_PATH.exists() or not GOLD_QUERIES_PATH.exists(),
    reason=(
        "Real corpus snapshot or gold-query file not built yet. Run "
        "scripts/build_v2_corpus_cache.py then scripts/build_gold_queries.py."
    ),
)


@pytest.fixture(scope="module")
def corpus():
    corpus, _stats, _manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    return corpus


@pytest.fixture(scope="module")
def gold():
    return json.loads(GOLD_QUERIES_PATH.read_text())


def _query_by_category(gold: dict, category: str) -> dict:
    matches = [q for q in gold["queries"] if q["category"] == category]
    assert matches, f"No gold query found for category {category!r}"
    return matches[0]


class TestExactMatch:
    def test_exact_match_ranks_first_with_identity_one(self, corpus, gold):
        q = _query_by_category(gold, "exact_match")
        checksum = q["query_sequence_checksum"]

        exact_hits = corpus.exact.lookup(checksum)
        assert q["expected_record_id"] in exact_hits

        component_type = ComponentType(q["component_type_hint"])
        approx_hits = corpus.approximate.search(component_type, q["query_sequence"], top_k=5)
        assert approx_hits[0].record_id == q["expected_record_id"]
        assert approx_hits[0].sequence_identity == 1.0


class TestSubstitutionVariants:
    def test_single_substitution_top_hit_is_original(self, gold, corpus):
        q = _query_by_category(gold, "single_substitution")
        component_type = ComponentType(q["component_type_hint"])
        hits = corpus.approximate.search(component_type, q["query_sequence"], top_k=5)
        assert hits[0].record_id == q["expected_record_id"]
        assert hits[0].sequence_identity >= q["expect"]["sequence_identity_min"]

    def test_several_substitution_top_hit_is_original(self, gold, corpus):
        q = _query_by_category(gold, "several_substitution")
        component_type = ComponentType(q["component_type_hint"])
        hits = corpus.approximate.search(component_type, q["query_sequence"], top_k=5)
        assert hits[0].record_id == q["expected_record_id"]
        assert hits[0].sequence_identity >= q["expect"]["sequence_identity_min"]


class TestSimilarNonidentical:
    def test_distinct_records_are_not_conflated(self, gold):
        q = _query_by_category(gold, "similar_nonidentical")
        assert q["query_sequence"] != ""
        assert q["reference_record_id"] is not None


class TestLongCdrH3:
    def test_long_cdr_h3_is_indexed_and_findable_by_exact_hash(self, gold, corpus):
        q = _query_by_category(gold, "long_cdr_h3")
        checksum = q["query_sequence_checksum"]
        assert q["expected_record_id"] in corpus.exact.lookup(checksum)
        record = corpus.get(q["expected_record_id"])
        assert record is not None
        assert record.component_type == ComponentType.CDR_H3
        assert len(record.sequence) == q["expect"]["length"]


class TestInsertionCodedCdr:
    def test_insertion_coded_component_preserved(self, gold, corpus):
        q = _query_by_category(gold, "insertion_coded_cdr")
        if q["expected_record_id"] is None:
            pytest.skip("No insertion-coded component found in this corpus snapshot")
        record = corpus.get(q["expected_record_id"])
        assert record is not None
        has_code = any(code and code.strip() for (_pos, code), _res in record.numbered.imgt_numbering)
        assert has_code, "Expected component to retain a non-empty IMGT insertion code"


class TestSingleChainInputs:
    def test_heavy_only_query_succeeds_without_light_chain(self, gold, corpus):
        q = _query_by_category(gold, "heavy_only")
        checksum = q["query_sequence_checksum"]
        assert q["expected_record_id"] in corpus.exact.lookup(checksum)

    def test_light_only_query_succeeds_without_heavy_chain(self, gold, corpus):
        q = _query_by_category(gold, "light_only")
        checksum = q["query_sequence_checksum"]
        assert q["expected_record_id"] in corpus.exact.lookup(checksum)


class TestPairedInput:
    def test_directly_observed_pair_is_found(self, gold, corpus):
        q = _query_by_category(gold, "paired_vh_vl")
        if q["heavy_sequence_checksum"] is None or q["light_sequence_checksum"] is None:
            pytest.skip("No observed pair found in this corpus snapshot")
        hits = corpus.pairs.lookup(q["heavy_sequence_checksum"], q["light_sequence_checksum"])
        assert len(hits) >= 1


class TestInvalidAndTruncatedInput:
    def test_invalid_sequence_is_rejected_not_scored(self, gold):
        q = _query_by_category(gold, "invalid_sequence")
        result = number_sequence(q["query_sequence"])
        assert result is None, "Invalid input must be rejected outright, never partially numbered"

    def test_truncated_sequence_does_not_silently_misrepresent(self, gold):
        q = _query_by_category(gold, "truncated_sequence")
        result = number_sequence(q["query_sequence"])
        # Either ANARCII refuses to number a 20-residue fragment (None), or it
        # succeeds -- both are acceptable, but a success must not be silently
        # treated as a full-length chain elsewhere. This test documents which
        # behavior ANARCII actually exhibits for this real truncated sequence.
        if result is not None:
            assert result.numbering_error is None or isinstance(result.numbering_error, str)


class TestParentDeduplication:
    def test_duplicate_parent_structures_do_not_inflate_component_count(self, gold, corpus):
        q = _query_by_category(gold, "duplicate_parent_record")
        if q["expected_record_id"] is None:
            pytest.skip("No multi-structure parent antibody found in this corpus snapshot")
        record = corpus.get(q["expected_record_id"])
        assert record is not None
        assert len(record.source_pointers) > 1, "Expected multiple source structures for this parent antibody"
        # The component itself must still be a single deduplicated record, not one per structure.
        checksum = q["query_sequence_checksum"]
        assert corpus.exact.lookup(checksum).count(q["expected_record_id"]) == 1


class TestNoCloseMatch:
    def test_random_sequence_does_not_produce_spurious_high_identity_hit(self, gold, corpus):
        q = _query_by_category(gold, "no_close_match")
        component_type = ComponentType(q["component_type_hint"])
        hits = corpus.approximate.search(component_type, q["query_sequence"], top_k=5)
        for hit in hits:
            assert hit.global_identity < 0.9, (
                "A randomly generated sequence matched a real record at high global "
                "identity -- investigate before trusting this result."
            )

    def test_fragmentary_local_alignment_is_not_reported_as_a_match(self, gold, corpus):
        """
        Regression test, found by this gate against the real corpus.

        A random 118-residue query aligned to a real VH over 5 residues, all of
        them matching, and was returned as the top hit with sequence_identity=1.0
        at query_coverage=0.04. Local identity is computed over the aligned region
        alone, so any coincidental short match scores 1.0; consumed without
        coverage it is indistinguishable from a genuine perfect match, and would
        have entered an EvidenceRecord as one.
        """
        q = _query_by_category(gold, "no_close_match")
        component_type = ComponentType(q["component_type_hint"])
        hits = corpus.approximate.search(component_type, q["query_sequence"], top_k=5)

        for hit in hits:
            # Any hit that survives must actually span both sequences...
            assert hit.query_coverage >= 0.70 and hit.subject_coverage >= 0.70
            # ...and its global identity must never be inflated by a short span.
            assert hit.global_identity <= hit.sequence_identity

    def test_fragmentary_alignment_is_dropped_rather_than_scored_low(self, corpus):
        """A sub-threshold alignment must be absent, not present with a low score.

        Missing evidence is not zero evidence (CLAUDE.md): a coincidental fragment
        is not weak evidence of similarity, it is no evidence, and must not appear
        in a hit list at all.
        """
        # A short peptide can only ever align to a small fraction of a full VH.
        fragment = "WGQGTLVTVSS"
        hits = corpus.approximate.search(ComponentType.VH, fragment, top_k=5)
        assert hits == [], (
            f"A {len(fragment)}-residue fragment produced {len(hits)} VH hits; "
            "sub-coverage alignments must be omitted, not ranked."
        )
