"""
Gate 8: pairing counts measured on the real combined corpus, pinned to both sources.

Pairing spans SAbDab2 and Thera-SAbDab, so both checksums guard these numbers --
an update to either source changes them legitimately, and the module skips with
re-measurement instructions rather than failing as a code regression.

Unit fixtures cannot see the failure mode these cover. A pairing bug that counts
source records as pairs looks perfectly correct on a two-row fixture, and only
diverges at scale: 12,032 records across 5,686 pairs.
"""

from collections import Counter
from pathlib import Path

import pytest

from antigens.models import ComponentType, EvidenceStatus
from antigens.pairing_precedent import compute_pairing_precedent
from antigens.persistence import load_corpus_snapshot

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")

PINNED_CHECKSUMS = {
    "SAbDab2": "f0b0aa6dae7d4e064cd09b843cc426c57368bb0bd7223dc3d910f4efe02596bc",
    "Thera-SAbDab": "cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d",
}

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

pytestmark = pytest.mark.skipif(
    not SNAPSHOT_PATH.exists(),
    reason="Combined corpus snapshot not built. Run scripts/build_combined_corpus_cache.py.",
)


@pytest.fixture(scope="module")
def corpus():
    loaded, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    stored = {s["source_name"]: s["raw_checksum"] for s in manifest.source_snapshots}
    for source, expected in PINNED_CHECKSUMS.items():
        actual = stored.get(source)
        if actual != expected:
            pytest.skip(
                f"{source} snapshot is {str(actual)[:16]}, these counts were measured "
                f"against {expected[:16]}. Re-measure with "
                f"scripts/measure_pairing_precedent.py and update every count here. "
                f"Do not adjust the numbers to make it pass."
            )
    return loaded


class TestTheTwoCountsStayDistinct:
    def test_unique_pairs_and_backing_records(self, corpus):
        """5,686 pairs, 12,032 records. Conflating them is the Gate 8 failure mode."""
        unique_keys = len(corpus.pairs)
        total_records = sum(len(v) for v in corpus.pairs._index.values())

        assert unique_keys == 5686
        assert total_records == 12032
        assert total_records / unique_keys == pytest.approx(2.12, abs=0.01)

    def test_most_pairs_rest_on_a_single_record(self, corpus):
        """3,060 of 5,686 have exactly one backing record.

        Worth pinning because it bounds how much corroboration actually exists:
        the majority of observed pairs are seen once, so any future weighting by
        observation count is weighting a long tail, not a typical case.
        """
        depth = Counter(len(v) for v in corpus.pairs._index.values())
        assert depth[1] == 3060
        assert sum(depth.values()) == 5686


class TestSourcesBarelyOverlap:
    def test_only_two_pairs_appear_in_both_sources(self, corpus):
        """4,562 SAbDab2-only, 1,122 Thera-SAbDab-only, 2 in both.

        The structural and therapeutic databases describe almost disjoint sets of
        chain pairs, so neither is redundant with the other and dropping either
        would lose most of the evidence, not a duplicated subset.
        """
        by_sources = Counter(
            tuple(sorted({r.source_database for r in records}))
            for records in corpus.pairs._index.values()
        )
        assert by_sources[("SAbDab2",)] == 4562
        assert by_sources[("Thera-SAbDab",)] == 1122
        assert by_sources[("SAbDab2", "Thera-SAbDab")] == 2


class TestStructureSupportAcrossPairs:
    def test_structure_backed_and_stated_absent_split(self, corpus):
        states = Counter()
        for records in corpus.pairs._index.values():
            available = [r.structure_available for r in records]
            if all(s is None for s in available):
                states["unknown"] += 1
            elif any(available):
                states["backed"] += 1
            else:
                states["absent"] += 1

        assert states["backed"] == 4852
        assert states["absent"] == 834
        assert states["unknown"] == 0
        assert sum(states.values()) == 5686


class TestRealCandidateBehaviour:
    def test_trastuzumab_pair_resolves_with_deduplicated_counts(self, corpus):
        """The exact pair is backed by 11 records and counts as one pair."""
        result = compute_pairing_precedent(
            corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL,
            top_k=5, shortlist_limit=500,
        )

        assert result.status is EvidenceStatus.AVAILABLE
        assert result.exact_pair_match_count == 1
        assert result.raw_observation_count > result.pairs_matched
        top = result.matches[0]
        assert top.heavy_global_identity == 1.0
        assert top.light_global_identity == 1.0
        assert top.observation_count == 11

    def test_every_emitted_pair_is_in_the_index(self, corpus):
        """At corpus scale, not just on a fixture: no pair is manufactured."""
        result = compute_pairing_precedent(
            corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL,
            top_k=10, shortlist_limit=500,
        )
        assert result.matches
        for match in result.matches:
            heavy = match.heavy_record_id.partition("::")[2]
            light = match.light_record_id.partition("::")[2]
            assert corpus.pairs.lookup(heavy, light)

    def test_single_domain_candidate_is_not_applicable(self, corpus):
        result = compute_pairing_precedent(
            corpus, TRASTUZUMAB_VH, None, heavy_component_type=ComponentType.VHH
        )
        assert result.status is EvidenceStatus.NOT_APPLICABLE


class TestCuratedContentAtCorpusScale:
    def test_most_pairs_are_structural_records_with_no_drug_name(self, corpus):
        """4,562 of 5,686. Any display design assuming a name and target exist
        renders empty fields for four results in five."""
        named = sum(
            1 for records in corpus.pairs._index.values()
            if any(r.metadata.get("therapeutic_name") for r in records)
        )
        assert named == 1124
        assert len(corpus.pairs) - named == 4562

    def test_one_chain_pair_can_carry_many_names_and_many_targets(self, corpus):
        """Trastuzumab's exact pair: 11 programmes, 6 distinct targets.

        They are bispecifics sharing a trastuzumab arm, so each targets HER2 plus
        something different. A single "target" field would show one of six. This is
        why the promoted fields are lists.
        """
        result = compute_pairing_precedent(
            corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL, top_k=5, shortlist_limit=500,
        )
        top = result.matches[0]

        assert top.observation_count == 11
        assert len(top.therapeutic_names) == 11
        assert len(top.targets) == 6
        assert all("HER2" in target for target in top.targets)
        # Mixed statuses across the programmes -- one value would misreport the rest.
        assert len(top.development_statuses) > 1

    def test_structural_matches_carry_method_and_resolution(self, corpus):
        result = compute_pairing_precedent(
            corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL, top_k=10, shortlist_limit=500,
        )
        structural = [m for m in result.matches if m.pdb_ids]
        assert structural, "expected at least one structure-backed match"
        for match in structural:
            assert match.structure_methods
            assert all(r > 0 for r in match.structure_resolutions)
            assert match.display_label.startswith("pdb_")


class TestPairRecordsAreWellFormed:
    def test_no_pair_key_is_a_sentinel_or_malformed_checksum(self, corpus):
        """The Gate 7 shape lesson applied to pairing.

        A truthy sentinel reaching a pair key would fabricate a pair rather than a
        structure flag -- worse, because pair counts drive deduplication. Checksums
        are computed rather than parsed, so this should hold trivially; it is
        asserted because "should hold trivially" is what was believed about the
        structure columns.
        """
        offenders = [
            key for key in corpus.pairs._index
            if not all(len(part) == 64 and part.isalnum() for part in key)
        ]
        assert offenders == []

    def test_every_pair_record_names_a_source_and_accession(self, corpus):
        offenders = [
            (record.source_database, record.source_accession)
            for records in corpus.pairs._index.values()
            for record in records
            if not record.source_database or not record.source_accession
            or str(record.source_accession).strip().lower() in {"none", "na", "nan", "null", ""}
        ]
        assert offenders == []
