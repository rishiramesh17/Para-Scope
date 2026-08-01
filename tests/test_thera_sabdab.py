"""
Tests for the Thera-SAbDab parser and corpus ingestion (Version 2, Gate 2).

Fixtures are synthetic, self-contained CSVs matching the real column schema
confirmed by direct inspection of the downloaded file (see
parsers/thera_sabdab.py module docstring) -- no dependency on the actual
downloaded corpus under data/, which is gitignored and may not be present.
"""

import csv
from pathlib import Path

import pytest

from antigens.corpus import Corpus
from antigens.models import ComponentType
from antigens.numbering import number_sequence
from antigens.parsers.thera_sabdab import (
    REQUIRED_COLUMNS,
    is_single_domain_format,
    parse_thera_sabdab_csv,
)

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTFTFTISSLQPEDFATYYC"
)
OTHER_VH = (
    "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS"
)
OTHER_VL = (
    "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIK"
)
NANOBODY_VHH = (
    "GSQVQLVESGGGLVQAGGSLRLSCAASGYAYTYIYMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYELRDRTYGQWGQGTQVTVSS"
)


def _write_thera_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = sorted(REQUIRED_COLUMNS)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _base_row(**overrides) -> dict:
    row = {
        "Therapeutic": "Testuzumab",
        "Format": "Whole mAb",
        "CH1 Isotype": "G1",
        "VD LC": "Kappa",
        "Highest_Clin_Trial (Feb '25)": "Approved",
        "Est. Status": "Active",
        "HeavySequence": TRASTUZUMAB_VH,
        "LightSequence": TRASTUZUMAB_VL,
        "HeavySequence(ifbispec)": "na",
        "LightSequence(ifbispec)": "na",
        # The structure columns write missing data as the literal 'None', NOT as
        # 'na'. This fixture said 'na' and so never exercised the real sentinel,
        # which is how 794 therapeutics came to be flagged structure-backed.
        "100% SI Structure": "None",
        "99% SI Structure": "None",
        "95-98% SI Structure": "None",
        "Target": "HER2",
        "Alternative Therapeutic Names": "na",
    }
    row.update(overrides)
    return row


class TestSingleDomainHeuristic:
    def test_detects_nanobody_keyword(self):
        assert is_single_domain_format("Nanobody")
        assert is_single_domain_format("Mixed Nanobody (VHH-CH2-CH3 dimer)")

    def test_detects_vhh_keyword(self):
        assert is_single_domain_format("Anti-target VHH")

    def test_vh_substring_alone_does_not_false_positive_as_vhh(self):
        assert not is_single_domain_format("VH-IgG")

    def test_conventional_format_not_flagged(self):
        assert not is_single_domain_format("Whole mAb")
        assert not is_single_domain_format("Fab")

    def test_none_format_defaults_false(self):
        assert not is_single_domain_format(None)
        assert not is_single_domain_format("")


class TestThereSAbDabParser:
    def test_parses_conventional_row_as_single_primary_arm(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row()])
        rows = list(parse_thera_sabdab_csv(csv_path))

        assert len(rows) == 1
        assert rows[0].supported
        assert len(rows[0].arms) == 1
        arm = rows[0].arms[0]
        assert arm.arm_index == 1
        assert arm.heavy_sequence == TRASTUZUMAB_VH
        assert arm.light_sequence == TRASTUZUMAB_VL
        assert not arm.is_single_domain

    def test_bispecific_row_yields_two_independent_arms(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(
            Format="Bispecific mAb",
            **{"HeavySequence(ifbispec)": OTHER_VH, "LightSequence(ifbispec)": OTHER_VL},
        )])
        rows = list(parse_thera_sabdab_csv(csv_path))

        assert len(rows[0].arms) == 2
        primary, bispec = rows[0].arms
        assert primary.arm_index == 1
        assert primary.heavy_sequence == TRASTUZUMAB_VH
        assert bispec.arm_index == 2
        assert bispec.heavy_sequence == OTHER_VH
        assert bispec.light_sequence == OTHER_VL

    def test_missing_light_chain_still_supported_as_heavy_only(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(
            Format="Nanobody", LightSequence="na", HeavySequence=NANOBODY_VHH,
        )])
        rows = list(parse_thera_sabdab_csv(csv_path))

        assert rows[0].supported
        assert len(rows[0].arms) == 1
        assert rows[0].arms[0].light_sequence is None
        assert rows[0].arms[0].is_single_domain

    def test_no_heavy_sequence_at_all_marks_row_unsupported(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(HeavySequence="na", LightSequence="na")])
        rows = list(parse_thera_sabdab_csv(csv_path))

        assert not rows[0].supported
        assert len(rows[0].arms) == 0
        assert "unsupported_component" in rows[0].skip_reason

    def test_na_sentinel_parsed_as_none_for_metadata_fields(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(Target="na")])
        rows = list(parse_thera_sabdab_csv(csv_path))
        assert rows[0].target is None

    @pytest.mark.parametrize("empty_value", ["", "   ", "NA", "nA"])
    def test_every_empty_shape_reaches_none_not_a_truthy_string(self, tmp_path, empty_value):
        """Blank, whitespace and case-variant sentinels all mean missing.

        Untested behaviour regresses. The structure columns proved this: they had
        a sentinel nobody had asserted on, and 794 therapeutics acquired a
        fabricated experimental-structure flag as a result.
        """
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(Target=empty_value)])
        rows = list(parse_thera_sabdab_csv(csv_path))
        assert rows[0].target is None

    def test_missing_required_column_raises(self, tmp_path):
        csv_path = tmp_path / "bad.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["Therapeutic"])
            writer.writeheader()
            writer.writerow({"Therapeutic": "x"})
        with pytest.raises(ValueError, match="missing required columns"):
            list(parse_thera_sabdab_csv(csv_path))


class TestCorpusThereSAbDabIngestion:
    def test_ingest_conventional_row_populates_pair_and_therapeutic_evidence(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row()])

        corpus = Corpus()
        stats = corpus.ingest_thera_sabdab(csv_path)

        assert stats.heavy_chains_parsed == 1
        assert stats.light_chains_parsed == 1
        assert stats.observed_pairs == 1
        assert stats.unique_parent_antibodies == 1

        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        record_id = f"{ComponentType.VH.value}::{heavy_numbered.sequence_checksum}"
        record = corpus.get(record_id)
        assert record is not None
        assert record.source_pointers[0][0] == "Thera-SAbDab"

    def test_bispecific_arms_are_not_cross_paired(self, tmp_path):
        """The primary arm's light chain must never pair with the bispecific arm's heavy chain."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(
            Therapeutic="Bispecifab", Format="Bispecific mAb",
            **{"HeavySequence(ifbispec)": OTHER_VH, "LightSequence(ifbispec)": OTHER_VL},
        )])

        corpus = Corpus()
        stats = corpus.ingest_thera_sabdab(csv_path)

        assert stats.observed_pairs == 2  # one pair per arm, never a cross-arm pair
        assert stats.unique_parent_antibodies == 1  # one therapeutic, two arms

        primary_heavy = number_sequence(TRASTUZUMAB_VH)
        primary_light = number_sequence(TRASTUZUMAB_VL)
        bispec_heavy = number_sequence(OTHER_VH)
        bispec_light = number_sequence(OTHER_VL)

        assert len(corpus.pairs.lookup(primary_heavy.sequence_checksum, primary_light.sequence_checksum)) == 1
        assert len(corpus.pairs.lookup(bispec_heavy.sequence_checksum, bispec_light.sequence_checksum)) == 1
        # Cross-arm combinations must never appear as observed pairs.
        assert corpus.pairs.lookup(primary_heavy.sequence_checksum, bispec_light.sequence_checksum) == []
        assert corpus.pairs.lookup(bispec_heavy.sequence_checksum, primary_light.sequence_checksum) == []

    def test_single_domain_format_maps_to_vhh_not_vh(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(
            Format="Nanobody", HeavySequence=NANOBODY_VHH, LightSequence="na",
        )])

        corpus = Corpus()
        stats = corpus.ingest_thera_sabdab(csv_path)
        assert stats.heavy_chains_parsed == 1
        assert stats.light_chains_parsed == 0

        nano_numbered = number_sequence(NANOBODY_VHH)
        vhh_record_id = f"{ComponentType.VHH.value}::{nano_numbered.sequence_checksum}"
        vh_record_id = f"{ComponentType.VH.value}::{nano_numbered.sequence_checksum}"
        assert corpus.get(vhh_record_id) is not None
        assert corpus.get(vh_record_id) is None

    def test_evidence_class_recorded_as_therapeutic_precedent_not_observed_repertoire(self, tmp_path):
        """Sanity check that Thera-SAbDab source pointers are tagged distinctly from SAbDab2's."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row()])
        corpus = Corpus()
        corpus.ingest_thera_sabdab(csv_path)

        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        record_id = f"{ComponentType.VH.value}::{heavy_numbered.sequence_checksum}"
        record = corpus.get(record_id)
        assert all(src == "Thera-SAbDab" for src, _accession in record.source_pointers)

    def test_same_sequence_from_both_sources_merges_source_pointers(self, tmp_path):
        """A therapeutic whose Fv also appears in SAbDab2 must dedupe to one component
        with both source databases recorded -- never double-counted (SPEC S15)."""
        import csv as csv_module

        from antigens.parsers.sabdab2 import REQUIRED_COLUMNS as SABDAB2_COLUMNS

        sabdab2_path = tmp_path / "ab_split.csv"
        with open(sabdab2_path, "w", newline="", encoding="utf-8") as f:
            fieldnames = sorted(SABDAB2_COLUMNS)
            writer = csv_module.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({
                "INSTANCE": "pdb_00000000_H_L", "PDB_ID": "pdb_00000000", "SABDAB_ID": "H0001L0001",
                "HEAVY_ID": "H0001", "LIGHT_ID": "L0001", "method": "XRAY", "resolution": "2.0",
                "type": "FAB", "construct": "", "holo": "False", "Hchain": "H", "Lchain": "L",
                "Hseq_expected": TRASTUZUMAB_VH, "Lseq_expected": TRASTUZUMAB_VL, "ab_split": "train",
            })

        thera_path = tmp_path / "thera.csv"
        _write_thera_csv(thera_path, [_base_row()])

        corpus = Corpus()
        corpus.ingest_sabdab2(sabdab2_path)
        corpus.ingest_thera_sabdab(thera_path)

        heavy_numbered = number_sequence(TRASTUZUMAB_VH)
        record_id = f"{ComponentType.VH.value}::{heavy_numbered.sequence_checksum}"
        record = corpus.get(record_id)
        source_dbs = {src for src, _accession in record.source_pointers}
        assert source_dbs == {"SAbDab2", "Thera-SAbDab"}


class TestStructureSentinelAndPerArmSplit:
    """
    The structure columns use their own missing-data sentinel and are split per arm
    (Gate 7).

    Measured on the real file before this was fixed: 794 of 1,131 therapeutics were
    recorded as having an experimental structure record when the source says they
    have none. SPEC S26 scores structure-record availability directly, so the flag
    was fabricated evidence, not a cosmetic error.
    """

    def test_none_sentinel_means_no_structure_not_a_structure_named_none(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row()])
        row = next(iter(parse_thera_sabdab_csv(csv_path)))
        arm = row.arms[0]
        assert arm.structures_100_si == []
        assert arm.has_structure_record is False

    def test_real_structure_ids_are_split_on_slash(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{"100% SI Structure": "1n8z:BA/6bi0:HL:IM"})])
        arm = next(iter(parse_thera_sabdab_csv(csv_path))).arms[0]
        assert arm.structures_100_si == ["1n8z:BA", "6bi0:HL:IM"]
        assert arm.has_structure_record is True

    def test_bispecific_arms_do_not_share_each_others_structures(self, tmp_path):
        """A structure solved for arm 2 is not evidence about arm 1."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{
            "HeavySequence(ifbispec)": OTHER_VH,
            "LightSequence(ifbispec)": OTHER_VL,
            "100% SI Structure": "None;8z8v:B",
            "99% SI Structure": "None;None",
            "95-98% SI Structure": "None;None",
        })])
        arms = next(iter(parse_thera_sabdab_csv(csv_path))).arms
        first, second = arms[0], arms[1]
        assert first.structures_100_si == []
        assert first.has_structure_record is False
        assert second.structures_100_si == ["8z8v:B"]
        assert second.has_structure_record is True

    def test_unmappable_column_is_unknown_rather_than_misattributed(self, tmp_path):
        """Under schema drift the arm mapping is unknown -- never guessed.

        A wrong arm is worse than an absent one: it would credit one specificity
        with another's experimental record.
        """
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{
            "HeavySequence(ifbispec)": OTHER_VH,
            "LightSequence(ifbispec)": OTHER_VL,
            # Two arms, three tokens -- the source shape this parser expects is broken.
            "100% SI Structure": "1n8z:BA;8z8v:B;7ah1:AA",
            "99% SI Structure": "None;None",
            "95-98% SI Structure": "None;None",
        })])
        arms = next(iter(parse_thera_sabdab_csv(csv_path))).arms
        for arm in arms:
            assert arm.structures_100_si is None
        # Still False, not None: the other two tiers mapped cleanly and say "none".
        assert arms[0].has_structure_record is False

    def test_unknown_is_never_coerced_to_false(self, tmp_path):
        """When every tier is unmappable the answer is None, not 'no structure'."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{
            "HeavySequence(ifbispec)": OTHER_VH,
            "LightSequence(ifbispec)": OTHER_VL,
            "100% SI Structure": "a;b;c",
            "99% SI Structure": "a;b;c",
            "95-98% SI Structure": "a;b;c",
        })])
        arms = next(iter(parse_thera_sabdab_csv(csv_path))).arms
        assert all(arm.has_structure_record is None for arm in arms)

    def test_blank_cell_is_unknown_not_a_statement_that_no_structure_exists(self, tmp_path):
        """A blank column asserts nothing, so it must not answer False.

        This is the 794-record defect inverted. There, a missing sentinel was read
        as evidence present; here, an empty cell was read as evidence of absence.
        Both fabricate a fact the source never stated, and this one is scored:
        has_structure_record False feeds SPEC S26's structure term as 0.0, while
        None leaves the term out of the profile entirely.
        """
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{
            "100% SI Structure": "",
            "99% SI Structure": "   ",
            "95-98% SI Structure": "",
        })])
        arm = next(iter(parse_thera_sabdab_csv(csv_path))).arms[0]

        assert arm.structures_100_si is None
        assert arm.structures_99_si is None
        assert arm.has_structure_record is None

    def test_blank_is_distinguishable_from_the_none_sentinel(self, tmp_path):
        """The two must not collapse to the same value -- they mean opposite things."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{"100% SI Structure": "None"})])
        stated_absent = next(iter(parse_thera_sabdab_csv(csv_path))).arms[0]

        blank_path = tmp_path / "thera_blank.csv"
        _write_thera_csv(blank_path, [_base_row(**{"100% SI Structure": ""})])
        unknown = next(iter(parse_thera_sabdab_csv(blank_path))).arms[0]

        assert stated_absent.structures_100_si == []    # source says: none exists
        assert unknown.structures_100_si is None        # source says: nothing
        assert stated_absent.structures_100_si != unknown.structures_100_si

    def test_blank_arm_token_does_not_silently_absolve_the_other_arm(self, tmp_path):
        """One unreadable arm makes the column's arm mapping unknown, not half-known."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row(**{
            "HeavySequence(ifbispec)": OTHER_VH,
            "LightSequence(ifbispec)": OTHER_VL,
            "100% SI Structure": "None;",
            "99% SI Structure": "None;None",
            "95-98% SI Structure": "None;None",
        })])
        arms = next(iter(parse_thera_sabdab_csv(csv_path))).arms
        assert all(arm.structures_100_si is None for arm in arms)

    def test_pair_record_carries_this_arms_structure_status(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [_base_row()])
        corpus = Corpus()
        corpus.ingest_thera_sabdab(csv_path)

        heavy = number_sequence(TRASTUZUMAB_VH)
        light = number_sequence(TRASTUZUMAB_VL)
        records = corpus.pairs.lookup(heavy.sequence_checksum, light.sequence_checksum)
        assert len(records) == 1
        assert records[0].structure_available is False
        assert records[0].metadata["structure_ids"] == []
