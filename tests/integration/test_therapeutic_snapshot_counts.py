"""
Gate 7: corpus-wide therapeutic counts, pinned to one exact source snapshot.

These numbers are properties of a *dataset version*, not of the code. Thera-SAbDab
is updated by its maintainers; when it is, these counts will legitimately change
and this module's expectations must be re-measured rather than patched to pass.

That is why they live here and not in the unit suite. Every assertion below is
guarded by the source file's SHA-256: against a different snapshot the module
skips with instructions, so a source update never presents as a code regression.
Unit tests assert *behaviour* (a blank cell means unknown); these assert *the
measured shape of one known dataset*, which is a different kind of claim with a
different lifetime.

The counts exist at all because two of them were fabricated before Gate 7. Every
therapeutic claimed an experimental structure; the true split is 337/794. A
regression there would be invisible to the unit tests, which use fixtures of a
handful of rows.
"""

from pathlib import Path

import pytest

from antigens.parsers.thera_sabdab import PARSER_VERSION
from antigens.persistence import load_corpus_snapshot

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")

# The exact Thera-SAbDab download every number below was measured against.
PINNED_SOURCE = "Thera-SAbDab"
PINNED_CHECKSUM = "cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d"

pytestmark = pytest.mark.skipif(
    not SNAPSHOT_PATH.exists(),
    reason=(
        "Combined corpus snapshot not built. Run "
        "scripts/build_combined_corpus_cache.py."
    ),
)


@pytest.fixture(scope="module")
def snapshot():
    corpus, stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    source = next(
        (s for s in manifest.source_snapshots if s["source_name"] == PINNED_SOURCE),
        None,
    )
    if source is None:
        pytest.skip(f"snapshot holds no {PINNED_SOURCE} source record")
    if source["raw_checksum"] != PINNED_CHECKSUM:
        pytest.skip(
            f"{PINNED_SOURCE} snapshot is {source['raw_checksum'][:16]}, these counts "
            f"were measured against {PINNED_CHECKSUM[:16]}. The source was updated: "
            f"re-measure and update PINNED_CHECKSUM together with every expected "
            f"count in this module. Do not adjust the numbers to make it pass."
        )
    return corpus, stats, manifest, source


@pytest.fixture(scope="module")
def therapeutics(snapshot):
    corpus = snapshot[0]
    return corpus._therapeutics


class TestSnapshotIsTheOneTheseNumbersDescribe:
    def test_parser_version_matches_the_running_code(self, snapshot):
        """A cache built by an older parser would report that parser's numbers.

        The checksum pins the *input*; this pins the code that read it. Both must
        match or the counts describe something other than what is being tested.
        """
        source = snapshot[3]
        assert source["parser_version"] == PARSER_VERSION, (
            f"snapshot was parsed by {source['parser_version']}, code is "
            f"{PARSER_VERSION}. Rebuild: scripts/build_combined_corpus_cache.py"
        )

    def test_source_licence_is_still_unspecified(self, snapshot):
        """Guards the release blocker, not the science.

        Both sources carry an unspecified licence, which blocks public deployment
        and redistribution. If this ever starts failing the licence changed, and
        that is a fact worth surfacing loudly rather than discovering later.
        """
        assert snapshot[3]["license_id"] == "unspecified"


class TestStructureEvidenceSplit:
    """The 794-record defect, measured on the real file."""

    def test_structure_backed_and_stated_absent_counts(self, therapeutics):
        backed = sum(
            1 for t in therapeutics.values()
            if any(arm.has_structure_record is True for arm in t.arms)
        )
        stated_absent = sum(
            1 for t in therapeutics.values()
            if all(arm.has_structure_record is False for arm in t.arms)
        )

        assert backed == 337
        assert stated_absent == 794
        assert backed + stated_absent == len(therapeutics) == 1131

    def test_not_every_therapeutic_claims_a_structure(self, therapeutics):
        """The shape of the original defect: all 1,131 answered True.

        Kept as a separate assertion from the exact counts because it survives a
        source update -- if this one ever fails again, the sentinel handling broke,
        whatever the current dataset happens to hold.
        """
        backed = sum(
            1 for t in therapeutics.values()
            if any(arm.has_structure_record is True for arm in t.arms)
        )
        assert 0 < backed < len(therapeutics)

    def test_no_arm_reports_unknown_on_this_snapshot(self, therapeutics):
        """3,393/3,393 structure values align to arms on this file.

        Unknown is the correct answer under schema drift, so a nonzero count here
        is not a bug -- it means the source shape changed and the per-arm split can
        no longer be trusted for the affected rows.
        """
        unknown = sum(
            1 for t in therapeutics.values()
            for arm in t.arms
            if arm.has_structure_record is None
        )
        assert unknown == 0


class TestEntityStoreShape:
    def test_entity_arm_and_link_counts(self, therapeutics, snapshot):
        corpus = snapshot[0]
        assert len(therapeutics) == 1131
        assert sum(len(t.arms) for t in therapeutics.values()) == 1288
        assert len(corpus._component_to_therapeutics) == 10186

    def test_heavy_only_arms_remain_reachable(self, therapeutics):
        """These 45 arms had no recoverable metadata before the entity store.

        Every single-domain therapeutic is among them. Metadata that depends on an
        arm happening to be paired is unreachable for exactly the constructs that
        never pair -- so this count dropping to 0 would mean they were silently
        excluded, not that the data improved.
        """
        heavy_only = [
            arm for t in therapeutics.values()
            for arm in t.arms
            if arm.light_record_id is None
        ]
        assert len(heavy_only) == 45
        assert all(arm.heavy_record_id is not None for arm in heavy_only)

    def test_ambiguous_aliases_are_preserved_not_resolved(self, snapshot):
        """70 aliases span co-administered combination products.

        `Evusheld (combination of cilgavimab and tixagevimab)` denotes two
        entities, so resolving it to either would be wrong. This count going to
        zero means something started guessing.
        """
        corpus = snapshot[0]
        assert len(corpus._therapeutic_aliases) == 6252
        assert len(corpus._ambiguous_therapeutic_aliases) == 70
