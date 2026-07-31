"""
Tests for the corpus persistence layer (Version 2, Gate 0.5).

The whole point of this layer is to avoid re-paying ANARCII's numbering cost on
every ingest -- measured at 511s of a 512s run over the real 15,641-row SAbDab2
split, versus 6.85s on a cache hit. These tests use small synthetic CSVs (same
pattern as test_search.py) so they run in seconds, but they specifically
verify that a cache hit means ANARCII is never invoked a second time --
not just that files round-trip correctly.
"""

import csv
from pathlib import Path

import pytest

import antigens.corpus as corpus_module
from antigens.corpus import Corpus, IngestStats, build_sabdab2_snapshot
from antigens.models import ComponentType, SourceSnapshot
from antigens.numbering import number_sequence, number_sequences
from antigens.parsers.sabdab2 import REQUIRED_COLUMNS
from antigens.persistence import (
    StaleSnapshotError,
    ingest_sabdab2_cached,
    load_corpus_snapshot,
    load_numbering_cache,
    save_corpus_snapshot,
    save_numbering_cache,
    snapshot_matches_source,
)

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTFTFTISSLQPEDFATYYC"
)


def _write_ab_split_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = sorted(REQUIRED_COLUMNS)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _base_row(**overrides) -> dict:
    row = {
        "INSTANCE": "pdb_00000000_H_L",
        "PDB_ID": "pdb_00000000",
        "SABDAB_ID": "H0001L0001",
        "HEAVY_ID": "H0001",
        "LIGHT_ID": "L0001",
        "method": "XRAY",
        "resolution": "2.0",
        "type": "FAB",
        "construct": "",
        "holo": "False",
        "Hchain": "H",
        "Lchain": "L",
        "Hseq_expected": TRASTUZUMAB_VH,
        "Lseq_expected": TRASTUZUMAB_VL,
        "ab_split": "train",
    }
    row.update(overrides)
    return row


class TestNumberingCachePersistence:
    def test_roundtrip_preserves_numbered_component(self, tmp_path):
        numbered = number_sequence(TRASTUZUMAB_VH)
        cache = {numbered.sequence_checksum: numbered}

        save_numbering_cache(tmp_path, "source-checksum-a", cache)
        loaded = load_numbering_cache(tmp_path, "source-checksum-a")

        assert loaded is not None
        restored = loaded[numbered.sequence_checksum]
        assert restored.original_sequence == numbered.original_sequence
        assert restored.chain_type == numbered.chain_type
        assert restored.imgt_numbering == numbered.imgt_numbering

    def test_missing_cache_returns_none(self, tmp_path):
        assert load_numbering_cache(tmp_path, "never-saved") is None


class TestIngestSabdab2Cached:
    def test_second_ingest_does_not_call_anarci_again(self, tmp_path, monkeypatch):
        """The core Gate 0.5 requirement: a cache hit must skip ANARCII entirely."""
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        cache_dir = tmp_path / "cache"

        call_count = {"n": 0}
        original = corpus_module.number_sequences

        def counting_wrapper(sequences, *args, **kwargs):
            call_count["n"] += 1
            return original(sequences, *args, **kwargs)

        monkeypatch.setattr(corpus_module, "number_sequences", counting_wrapper)

        corpus1, stats1, _ = ingest_sabdab2_cached(csv_path, cache_dir)
        assert call_count["n"] == 1
        assert stats1.heavy_chains_parsed == 1

        corpus2, stats2, _ = ingest_sabdab2_cached(csv_path, cache_dir)
        assert call_count["n"] == 1  # unchanged -- second ingest must not renumber
        assert stats2.heavy_chains_parsed == 1
        assert len(corpus2) == len(corpus1)

    def test_cached_ingest_produces_equivalent_search_results(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        cache_dir = tmp_path / "cache"

        corpus, _, _ = ingest_sabdab2_cached(csv_path, cache_dir)
        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        assert corpus.exact.lookup(heavy_numbered.sequence_checksum) != []


class TestCorpusIngestionWithPartialCache:
    def test_merges_new_sequences_not_in_supplied_cache(self, tmp_path):
        """A partial cache (covering only some sequences) must be extended, not replace-or-fail."""
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])

        # Pre-populate a cache with only the light chain numbered; heavy is absent.
        light_numbered = number_sequence(TRASTUZUMAB_VL)
        partial_cache = {light_numbered.sequence_checksum: light_numbered}

        corpus = Corpus()
        stats = IngestStats()
        corpus.ingest_sabdab2(csv_path, stats=stats, numbering_cache=partial_cache)

        assert stats.heavy_chains_parsed == 1
        assert stats.light_chains_parsed == 1
        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        assert heavy_numbered.sequence_checksum in corpus.last_numbering_cache
        assert light_numbered.sequence_checksum in corpus.last_numbering_cache


class TestCorpusSnapshotPersistence:
    def test_roundtrip_preserves_all_indexes(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])

        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        snapshot = build_sabdab2_snapshot(csv_path)

        snapshot_path = tmp_path / "corpus.pkl"
        manifest = save_corpus_snapshot(corpus, stats, [snapshot], snapshot_path)
        assert manifest.stats["heavy_chains_parsed"] == 1

        loaded_corpus, loaded_stats, loaded_manifest = load_corpus_snapshot(snapshot_path)

        assert len(loaded_corpus) == len(corpus)
        assert loaded_stats.heavy_chains_parsed == stats.heavy_chains_parsed

        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        assert loaded_corpus.exact.lookup(heavy_numbered.sequence_checksum) != []
        record_id = f"{ComponentType.VH.value}::{heavy_numbered.sequence_checksum}"
        assert loaded_corpus.get(record_id) is not None

        light_numbered = number_sequence(TRASTUZUMAB_VL)
        pair_hits = loaded_corpus.pairs.lookup(heavy_numbered.sequence_checksum, light_numbered.sequence_checksum)
        assert len(pair_hits) == 1

    def test_stale_corpus_schema_is_rejected_not_silently_loaded(self, tmp_path):
        """
        A snapshot pickled under an older Corpus shape must fail loudly on load.

        Regression: the real SAbDab2 snapshot was built before the Thera-SAbDab
        bookkeeping attribute was added to Corpus.__init__ (now the therapeutic
        entity store, schema 5.0.0). Unpickling restores the stored __dict__
        verbatim, so the loaded object silently lacked that attribute and only
        failed much later, as an AttributeError deep inside ingest_thera_sabdab.
        """
        import pickle

        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        snapshot = build_sabdab2_snapshot(csv_path)

        snapshot_path = tmp_path / "corpus.pkl"
        save_corpus_snapshot(corpus, stats, [snapshot], snapshot_path)

        # Rewrite the manifest to claim an older Corpus schema.
        with open(snapshot_path, "rb") as f:
            blob = pickle.load(f)
        blob["manifest"].corpus_schema_version = "1.0.0"
        with open(snapshot_path, "wb") as f:
            pickle.dump(blob, f)

        with pytest.raises(StaleSnapshotError, match="1.0.0"):
            load_corpus_snapshot(snapshot_path)

    def test_stale_cache_format_is_rejected(self, tmp_path):
        import pickle

        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        snapshot = build_sabdab2_snapshot(csv_path)

        snapshot_path = tmp_path / "corpus.pkl"
        save_corpus_snapshot(corpus, stats, [snapshot], snapshot_path)

        with open(snapshot_path, "rb") as f:
            blob = pickle.load(f)
        blob["manifest"].cache_format_version = "0.9.0"
        with open(snapshot_path, "wb") as f:
            pickle.dump(blob, f)

        with pytest.raises(StaleSnapshotError):
            load_corpus_snapshot(snapshot_path)

    def test_snapshot_matches_source_true_for_same_file(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        snapshot = build_sabdab2_snapshot(csv_path)

        manifest = save_corpus_snapshot(corpus, stats, [snapshot], tmp_path / "corpus.pkl")
        assert snapshot_matches_source(manifest, snapshot)

    def test_snapshot_matches_source_false_for_different_file(self, tmp_path):
        csv_path_a = tmp_path / "a.csv"
        _write_ab_split_csv(csv_path_a, [_base_row()])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path_a)
        snapshot_a = build_sabdab2_snapshot(csv_path_a)
        manifest = save_corpus_snapshot(corpus, stats, [snapshot_a], tmp_path / "corpus.pkl")

        csv_path_b = tmp_path / "b.csv"
        _write_ab_split_csv(csv_path_b, [_base_row(INSTANCE="pdb_other_H_L", PDB_ID="pdb_other")])
        snapshot_b = build_sabdab2_snapshot(csv_path_b)

        assert not snapshot_matches_source(manifest, snapshot_b)


class TestSourceProvenanceRoundTrip:
    """
    A stored snapshot must be answerable about its own provenance.

    `_snapshot_to_dict` previously persisted four of SourceSnapshot's eight fields
    and silently dropped the rest. license_id was the consequential loss: SPEC
    §13/§14 require a source's licence to travel with evidence derived from it,
    both ingested sources currently carry `unspecified`, and a redistribution
    question cannot be answered by re-reading source files that may have moved or
    changed. These tests fail if any field is dropped again.
    """

    def _saved_manifest(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        snapshot = build_sabdab2_snapshot(csv_path)
        save_corpus_snapshot(corpus, stats, [snapshot], tmp_path / "corpus.pkl")
        _corpus, _stats, manifest = load_corpus_snapshot(tmp_path / "corpus.pkl")
        return manifest, snapshot

    def test_every_source_snapshot_field_survives_a_round_trip(self, tmp_path):
        manifest, snapshot = self._saved_manifest(tmp_path)
        stored = manifest.source_snapshots[0]
        for field_name in (
            "source_name", "snapshot_id", "retrieved_at", "license_id",
            "raw_checksum", "record_count", "parser_version", "schema_version",
        ):
            assert field_name in stored, f"{field_name} dropped at save time"
            assert stored[field_name] == getattr(snapshot, field_name)

    def test_license_is_recoverable_from_a_loaded_snapshot(self, tmp_path):
        """The specific field whose loss motivated this test class."""
        manifest, snapshot = self._saved_manifest(tmp_path)
        assert manifest.source_snapshots[0]["license_id"] == snapshot.license_id

    def test_stored_dict_covers_all_dataclass_fields(self, tmp_path):
        """Adding a field to SourceSnapshot must not silently bypass persistence."""
        from dataclasses import fields

        manifest, _snapshot = self._saved_manifest(tmp_path)
        declared = {f.name for f in fields(SourceSnapshot)}
        assert declared == set(manifest.source_snapshots[0])
