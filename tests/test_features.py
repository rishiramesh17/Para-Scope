"""
Tests for the deterministic feature engine (Version 2, Gate 4).

The contract under test is determinism and honesty: identical input gives
identical output, non-standard residues are refused rather than coerced, and
every approximate value carries its limitation.
"""

import pytest

from para_scope.features import (
    FEATURE_ENGINE_VERSION,
    STANDARD_AMINO_ACIDS,
    NonStandardResidueError,
    compute_features,
    count_n_glycosylation_sequons,
)

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)


class TestDeterminism:
    def test_same_sequence_gives_identical_features(self):
        a = compute_features(TRASTUZUMAB_VH)
        b = compute_features(TRASTUZUMAB_VH)
        assert a == b

    def test_version_is_stamped_on_every_result(self):
        assert compute_features(TRASTUZUMAB_VH).feature_engine_version == FEATURE_ENGINE_VERSION

    def test_whitespace_and_case_are_normalized_without_changing_values(self):
        plain = compute_features(TRASTUZUMAB_VH)
        messy = compute_features(f"  {TRASTUZUMAB_VH.lower()}  ")
        assert plain == messy

    def test_features_are_json_serializable(self):
        import json
        from dataclasses import asdict

        json.dumps(asdict(compute_features(TRASTUZUMAB_VH)))


class TestBasicMeasurements:
    def test_length_matches_input(self):
        assert compute_features(TRASTUZUMAB_VH).sequence_length == len(TRASTUZUMAB_VH)

    def test_composition_sums_to_one_and_covers_all_standard_residues(self):
        f = compute_features(TRASTUZUMAB_VH)
        assert set(f.composition) == set(STANDARD_AMINO_ACIDS)
        assert f.composition["C"] == TRASTUZUMAB_VH.count("C") / len(TRASTUZUMAB_VH)
        assert abs(sum(f.composition.values()) - 1.0) < 1e-9

    def test_net_charge_is_computed_from_charged_residue_counts(self):
        seq = "KKKDD"
        assert compute_features(seq).net_charge_at_ph_7 == 1

    def test_gravy_matches_hand_computation_for_uniform_sequence(self):
        # All-isoleucine: GRAVY must equal the Kyte-Doolittle value for I exactly.
        assert compute_features("IIIIIIIIII").gravy == pytest.approx(4.5)

    def test_cysteine_count_matches_input(self):
        f = compute_features(TRASTUZUMAB_VH)
        assert f.cysteine_count == TRASTUZUMAB_VH.count("C")

    def test_vh_and_vl_produce_distinct_feature_sets(self):
        assert compute_features(TRASTUZUMAB_VH) != compute_features(TRASTUZUMAB_VL)


class TestNonStandardResidues:
    def test_ambiguous_residue_is_refused_not_substituted(self):
        """X must never be silently dropped or coerced -- that would alter the input."""
        with pytest.raises(NonStandardResidueError, match="X"):
            compute_features("EVQLXESGGG")

    def test_empty_sequence_is_refused(self):
        with pytest.raises(NonStandardResidueError):
            compute_features("")

    def test_error_names_every_offending_residue(self):
        with pytest.raises(NonStandardResidueError) as exc:
            compute_features("EVQLBZGGG")
        assert "B" in str(exc.value) and "Z" in str(exc.value)


class TestGlycosylationSequons:
    def test_counts_simple_sequon(self):
        assert count_n_glycosylation_sequons("NAS") == 1

    def test_threonine_also_completes_the_motif(self):
        assert count_n_glycosylation_sequons("NAT") == 1

    def test_proline_at_x_blocks_the_motif(self):
        assert count_n_glycosylation_sequons("NPS") == 0

    def test_counts_multiple_occurrences(self):
        assert count_n_glycosylation_sequons("NASNAT") == 2

    def test_no_sequon_in_unrelated_sequence(self):
        assert count_n_glycosylation_sequons("AAAAAA") == 0

    def test_surfaced_on_feature_set(self):
        f = compute_features("EVQLVESGGGNASLVQPGGS")
        assert f.n_glycosylation_sequon_count == 1


class TestFreeCysteineEstimate:
    def test_odd_cysteine_count_flags_one_unpairable(self):
        f = compute_features("CAAAA")
        assert f.free_cysteine_count_estimate == 1

    def test_even_cysteine_count_estimates_zero(self):
        f = compute_features("CAAAAC")
        assert f.free_cysteine_count_estimate == 0

    def test_no_cysteine_yields_none_not_zero(self):
        """Absent is not the same as zero: there is nothing to estimate here."""
        f = compute_features("AAAAA")
        assert f.free_cysteine_count_estimate is None

    def test_estimate_carries_its_limitation(self):
        """An approximate value must never be reported bare."""
        f = compute_features("CAAAA")
        assert any("parity" in note for note in f.notes)
        assert any("structural" in note for note in f.notes)

    def test_no_limitation_note_when_nothing_was_estimated(self):
        assert compute_features("AAAAA").notes == []
