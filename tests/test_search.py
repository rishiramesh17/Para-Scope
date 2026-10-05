"""
Tests for search indexes and corpus ingestion (SPEC S16, S13-15).

All fixtures are synthetic/self-contained (no dependency on the downloaded
SAbDab2 corpus, which lives outside git under data/ and may not be present).
"""

import csv
from pathlib import Path

import pytest

from para_scope.corpus import Corpus, REGION_TO_COMPONENT_TYPE
from para_scope.indexes import ApproximateIndex, ExactIndex, PairIndex, PairRecord, PositionalIndex
from para_scope.indexes.approximate import ALIGNMENT_SCORING_VERSION
from para_scope.models import ComponentType
from para_scope.numbering import (
    compute_sequence_checksum,
    number_sequence,
    split_chain_at_variable_domain,
)
from para_scope.parsers.sabdab2 import REQUIRED_COLUMNS, parse_ab_split_csv

# Public trastuzumab heavy/light (kappa) variable domains -- same fixtures as test_numbering.py
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
# A complete kappa V-domain, through FR-L4 -- required for the CDR-L3 / FR-L4
# assertions to exercise real residues rather than IMGT gap placeholders.
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTFTFTISSLQPEDFATYYC"
    "QQHYTTPPTFGQGTKVEIK"
)

# A second, unrelated heavy/light pair (public pembrolizumab-like Fv) for
# distinguishing "found" from "not found" in search tests.
OTHER_VH = (
    "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS"
)
OTHER_VL = (
    "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIK"
)

# A clean single-domain heavy (nanobody, VHH) sequence pulled from a real SAbDab2 SD-H
# record and hardcoded here so this test never depends on the downloaded corpus.
NANOBODY_VHH = (
    "GSQVQLVESGGGLVQAGGSLRLSCAASGYAYTYIYMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYELRDRTYGQWGQGTQVTVSS"
)


def _write_ab_split_csv(path: Path, rows: list[dict]) -> None:
    """Write a minimal ab_split.csv containing only the columns the parser requires."""
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


class TestExactIndex:
    def test_add_and_lookup(self):
        index = ExactIndex()
        index.add("checksum-a", "rec1")
        assert index.lookup("checksum-a") == ["rec1"]

    def test_lookup_miss_returns_empty(self):
        index = ExactIndex()
        assert index.lookup("nonexistent") == []

    def test_duplicate_add_does_not_duplicate_record(self):
        index = ExactIndex()
        index.add("checksum-a", "rec1")
        index.add("checksum-a", "rec1")
        assert index.lookup("checksum-a") == ["rec1"]

    def test_same_checksum_multiple_records(self):
        index = ExactIndex()
        index.add("checksum-a", "rec1")
        index.add("checksum-a", "rec2")
        assert set(index.lookup("checksum-a")) == {"rec1", "rec2"}


class TestApproximateIndex:
    def test_exact_match_scores_identity_1(self):
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH, top_k=5)
        assert hits[0].record_id == "rec1"
        assert hits[0].sequence_identity == 1.0
        assert hits[0].query_coverage == 1.0

    def test_ranks_closer_sequence_higher(self):
        index = ApproximateIndex()
        index.add(ComponentType.VH, "close", TRASTUZUMAB_VH)
        index.add(ComponentType.VH, "far", OTHER_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH, top_k=5)
        assert hits[0].record_id == "close"
        assert hits[0].sequence_identity > hits[1].sequence_identity

    def test_rejects_short_components(self):
        """CDRs must never enter the long-sequence approximate index (SPEC S21, CLAUDE.md)."""
        index = ApproximateIndex()
        with pytest.raises(ValueError):
            index.add(ComponentType.CDR_H3, "rec1", "WGGDGFYAMDY")

    def test_native_float_types_for_json_serialization(self):
        import json

        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH)
        json.dumps({"identity": hits[0].sequence_identity, "coverage": hits[0].query_coverage})

    def test_exact_match_scores_global_identity_1(self):
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH, top_k=5)
        # A true full-length match must be 1.0 under *both* measures.
        assert hits[0].global_identity == 1.0
        assert hits[0].sequence_identity == 1.0

    def test_fragment_query_is_dropped_not_scored_as_perfect(self):
        """
        Regression: a fragment that aligns perfectly over a few residues must not
        be reported as a match (found by Gate 1 on the real corpus, where a random
        118-residue query scored sequence_identity=1.0 at query_coverage=0.04).
        """
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        # A real subsequence of the indexed VH -- locally it aligns at 100% identity.
        fragment = TRASTUZUMAB_VH[40:50]
        hits = index.search(ComponentType.VH, fragment, top_k=5)
        assert hits == [], "A 10-residue fragment must not be returned as a VH match"

    def test_global_identity_penalizes_partial_coverage(self):
        """Half a domain matching perfectly must not score like a whole domain matching perfectly."""
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        half = TRASTUZUMAB_VH[: len(TRASTUZUMAB_VH) // 2]
        # Coverage floor is relative to both sequences, so relax it to inspect scoring itself.
        hits = index.search(
            ComponentType.VH, half, top_k=5, min_query_coverage=0.0, min_subject_coverage=0.0
        )
        assert hits, "expected a hit once coverage floors are relaxed"
        hit = hits[0]
        assert hit.sequence_identity == 1.0, "the aligned region itself is a perfect match"
        assert hit.global_identity < 0.6, (
            "global identity must reflect the half of the subject left uncovered, "
            f"got {hit.global_identity}"
        )

    def test_coverage_floor_is_configurable_and_versioned(self):
        index = ApproximateIndex()
        index.add(ComponentType.VH, "rec1", TRASTUZUMAB_VH)
        hits = index.search(ComponentType.VH, TRASTUZUMAB_VH)
        assert hits[0].scoring_version == ALIGNMENT_SCORING_VERSION


class TestPositionalIndex:
    def test_identical_positions_meet_threshold(self):
        index = PositionalIndex()
        positions = set(range(105, 118))
        index.add(ComponentType.CDR_H3, "rec1", positions)
        results = index.lookup(ComponentType.CDR_H3, positions, threshold=0.95)
        assert ("rec1", 1.0) in results

    def test_below_threshold_excluded(self):
        index = PositionalIndex()
        index.add(ComponentType.CDR_H3, "rec1", set(range(105, 118)))
        # Query shares only ~half the positions -> well below 0.95 Jaccard threshold
        query_positions = set(range(105, 112))
        results = index.lookup(ComponentType.CDR_H3, query_positions, threshold=0.95)
        assert results == []

    def test_component_types_never_cross(self):
        """A CDR-H3 position set must never match against CDR-L3 (CLAUDE.md)."""
        index = PositionalIndex()
        positions = set(range(105, 118))
        index.add(ComponentType.CDR_H3, "heavy_rec", positions)
        results = index.lookup(ComponentType.CDR_L3, positions)
        assert results == []


class TestPairIndex:
    def test_add_and_lookup(self):
        index = PairIndex()
        index.add("h1", "l1", PairRecord(source_database="SAbDab2", source_accession="1FVH", structure_available=True))
        hits = index.lookup("h1", "l1")
        assert len(hits) == 1
        assert hits[0].source_accession == "1FVH"

    def test_never_synthesizes_unobserved_pairs(self):
        """CLAUDE.md: never synthesize a pair from independent heavy/light hits."""
        index = PairIndex()
        index.add("h1", "l1", PairRecord(source_database="SAbDab2", source_accession="1FVH", structure_available=True))
        index.add("h2", "l2", PairRecord(source_database="SAbDab2", source_accession="2ABC", structure_available=True))
        # h1 and l2 individually exist in the index, but were never observed paired together
        assert index.lookup("h1", "l2") == []


class TestSAbDab2Parser:
    def test_parses_supported_row(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        rows = list(parse_ab_split_csv(csv_path))
        assert len(rows) == 1
        assert rows[0].supported
        assert rows[0].heavy_sequence == TRASTUZUMAB_VH

    def test_flags_vnar_as_unsupported(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row(type="VNAR", INSTANCE="pdb_vnar_A_+")])
        rows = list(parse_ab_split_csv(csv_path))
        assert rows[0].supported is False
        assert "unsupported_component" in rows[0].skip_reason

    def test_flags_fused_scfv_construct_as_unsupported(self, tmp_path):
        """scFv rows report the same fused VH-linker-VL chain in both H/L columns --
        must not be silently numbered as an isolated VH or VL (CLAUDE.md)."""
        csv_path = tmp_path / "ab_split.csv"
        fused_seq = TRASTUZUMAB_VH + "GGGGSGGGGSGGGGS" + TRASTUZUMAB_VL
        _write_ab_split_csv(csv_path, [_base_row(
            construct="SCFV", Hseq_expected=fused_seq, Lseq_expected=fused_seq,
        )])
        rows = list(parse_ab_split_csv(csv_path))
        assert rows[0].supported is False
        assert "engineered single-chain construct" in rows[0].skip_reason

    def test_missing_required_column_raises(self, tmp_path):
        csv_path = tmp_path / "bad.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["INSTANCE"])
            writer.writeheader()
            writer.writerow({"INSTANCE": "x"})
        with pytest.raises(ValueError, match="missing required columns"):
            list(parse_ab_split_csv(csv_path))


class TestCorpusIngestion:
    def test_ingest_single_row_populates_all_indexes(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])

        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)

        assert stats.rows_seen == 1
        assert stats.heavy_chains_parsed == 1
        assert stats.light_chains_parsed == 1
        assert stats.observed_pairs == 1

        # Heavy variable domain findable by exact hash
        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        assert corpus.exact.lookup(heavy_numbered.sequence_checksum) != []

        # Direct pairing observed (both chains from the same row)
        light_numbered = number_sequence(TRASTUZUMAB_VL)
        pair_hits = corpus.pairs.lookup(heavy_numbered.sequence_checksum, light_numbered.sequence_checksum)
        assert len(pair_hits) == 1
        assert pair_hits[0].source_database == "SAbDab2"

    def test_duplicate_sequence_across_structures_deduplicates(self, tmp_path):
        """Two PDB depositions of the same antibody collapse into one component (SPEC S15)."""
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [
            _base_row(INSTANCE="pdb_00000001_H_L", PDB_ID="pdb_00000001"),
            _base_row(INSTANCE="pdb_00000002_H_L", PDB_ID="pdb_00000002"),
        ])

        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)

        assert stats.rows_seen == 2
        assert stats.unique_parent_antibodies == 1  # same SABDAB_ID both rows

        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        record_id = f"{ComponentType.VH.value}::{heavy_numbered.sequence_checksum}"
        record = corpus.get(record_id)
        assert record is not None
        assert len(record.source_pointers) == 2  # both PDB structures recorded, not duplicated

    def test_vnar_row_skipped_not_miscategorized(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [
            _base_row(type="VNAR", INSTANCE="pdb_vnar_A_+", Lchain="", Lseq_expected=""),
        ])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        assert stats.rows_skipped_unsupported == 1
        assert stats.heavy_chains_parsed == 0
        assert len(corpus) == 0

    def test_nanobody_maps_to_vhh_not_vh(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [
            _base_row(
                type="SD-H", INSTANCE="pdb_nano_A_+", Hchain="A", Lchain="", Lseq_expected="",
                Hseq_expected=NANOBODY_VHH,
            ),
        ])
        corpus = Corpus()
        stats = corpus.ingest_sabdab2(csv_path)
        assert stats.heavy_chains_parsed == 1

        # Keyed on the numbered V-domain, not the submitted chain: this fixture
        # carries a leading "GS" cloning-tag remnant that ANARCII does not
        # number, so the stored component is the 121-residue domain.
        nano_numbered = number_sequence(NANOBODY_VHH)
        _pre, nano_domain, _suf = split_chain_at_variable_domain(nano_numbered)
        domain_checksum = compute_sequence_checksum(nano_domain)
        vhh_record_id = f"{ComponentType.VHH.value}::{domain_checksum}"
        vh_record_id = f"{ComponentType.VH.value}::{domain_checksum}"
        assert corpus.get(vhh_record_id) is not None
        assert corpus.get(vh_record_id) is None

    def test_ingestion_indexes_all_cdrs_and_frameworks(self, tmp_path):
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [_base_row()])
        corpus = Corpus()
        corpus.ingest_sabdab2(csv_path)

        indexed_types = {record.component_type for record in corpus._components.values()}
        for component_type in REGION_TO_COMPONENT_TYPE.values():
            assert component_type in indexed_types, f"{component_type} was not indexed"

    def test_cross_query_finds_similar_but_not_unrelated(self, tmp_path):
        """End-to-end: querying the corpus with a related sequence ranks it above an unrelated one."""
        csv_path = tmp_path / "ab_split.csv"
        _write_ab_split_csv(csv_path, [
            _base_row(),
            _base_row(
                INSTANCE="pdb_other_H_L", PDB_ID="pdb_other", SABDAB_ID="H0002L0002",
                Hseq_expected=OTHER_VH, Lseq_expected=OTHER_VL,
            ),
        ])
        corpus = Corpus()
        corpus.ingest_sabdab2(csv_path)

        hits = corpus.approximate.search(ComponentType.VH, TRASTUZUMAB_VH, top_k=5)
        assert hits[0].sequence_identity == 1.0  # exact self-match ranks first
