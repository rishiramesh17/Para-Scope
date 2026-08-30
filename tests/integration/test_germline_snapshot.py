"""
Gate 9 milestone 1: germline assignment against real corpus sequences,
pinned to one exact bundled reference set and one exact corpus snapshot.

Two independent checksums guard these numbers, because either input can move
independently: the bundled OGRDB reference file (re-running
scripts/build_germline_reference.py against a newer OGRDB release) or the
Thera-SAbDab source snapshot the queried sequence comes from. Against a
different pairing of either, the module skips with instructions rather than
failing as a code regression -- these assertions are about the measured shape
of one known (reference set, corpus) pair, not about behaviour.
"""

import hashlib
from pathlib import Path

import pytest

from para_scope.germline import REFERENCE_SET_VERSION, _REFERENCE_FILE, compute_germline_assignment
from para_scope.models import ComponentType
from para_scope.persistence import load_corpus_snapshot

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")

PINNED_REFERENCE_CHECKSUM = (
    "6f9214189bed065646fd4bc30bfd8107b47ada751bc37c115db0feb3ce4c568e"
)
PINNED_THERA_SABDAB_SOURCE = "Thera-SAbDab"
PINNED_THERA_SABDAB_CHECKSUM = (
    "cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d"
)

pytestmark = pytest.mark.skipif(
    not SNAPSHOT_PATH.exists(),
    reason=(
        "Combined corpus snapshot not built. Run "
        "scripts/build_combined_corpus_cache.py."
    ),
)


@pytest.fixture(scope="module")
def corpus():
    reference_checksum = hashlib.sha256(_REFERENCE_FILE.read_bytes()).hexdigest()
    if reference_checksum != PINNED_REFERENCE_CHECKSUM:
        pytest.skip(
            f"Bundled germline reference is {reference_checksum[:16]}, these "
            f"numbers were measured against {PINNED_REFERENCE_CHECKSUM[:16]} "
            f"(REFERENCE_SET_VERSION={REFERENCE_SET_VERSION}). Re-measure and "
            f"update the pinned checksum and expected calls together. Do not "
            f"adjust the numbers to make it pass."
        )

    corpus, stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    source = next(
        (
            s
            for s in manifest.source_snapshots
            if s["source_name"] == PINNED_THERA_SABDAB_SOURCE
        ),
        None,
    )
    if source is None:
        pytest.skip(f"snapshot holds no {PINNED_THERA_SABDAB_SOURCE} source record")
    if source["raw_checksum"] != PINNED_THERA_SABDAB_CHECKSUM:
        pytest.skip(
            f"{PINNED_THERA_SABDAB_SOURCE} snapshot is "
            f"{source['raw_checksum'][:16]}, these numbers were measured "
            f"against {PINNED_THERA_SABDAB_CHECKSUM[:16]}. Re-measure and "
            f"update the pinned checksum and expected calls together."
        )
    return corpus


class TestTrastuzumabGermlineAssignment:
    """Measured once against the pinned (reference set, corpus) pair above.

    Re-derive rather than hand-adjust if either checksum ever changes.
    """

    def test_heavy_chain_nearest_germline(self, corpus):
        arm = corpus.get_therapeutic("Trastuzumab").arms[0]
        vh_sequence = corpus._components[arm.heavy_record_id].sequence

        result = compute_germline_assignment(ComponentType.VH, vh_sequence)

        assert result.nearest_gene_call == "IGHV3-66"
        assert result.nearest_germline_identity == pytest.approx(0.6667, abs=1e-3)

    def test_light_chain_nearest_germline(self, corpus):
        arm = corpus.get_therapeutic("Trastuzumab").arms[0]
        vl_sequence = corpus._components[arm.light_record_id].sequence

        result = compute_germline_assignment(ComponentType.VL_KAPPA, vl_sequence)

        assert result.nearest_gene_call == "IGKV1-39"
        assert result.nearest_germline_identity == pytest.approx(0.7664, abs=1e-3)
