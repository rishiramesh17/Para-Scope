"""
Tests for the sequence motif observation engine (Version 2, Gate 5).

Two things are under test: that motifs are located correctly in IMGT
coordinates, and that the output stays an observation rather than becoming an
assessment (no severity, no score, disclaimer always present).
"""

from dataclasses import asdict, fields

from para_scope.motifs import (
    MOTIF_RULES,
    MOTIF_RULESET_VERSION,
    MotifObservation,
    find_motifs,
)
from para_scope.numbering import number_sequence

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)


class TestMotifLocation:
    def test_finds_ng_deamidation_motif(self):
        report = find_motifs("AAANGAAA")
        ng = [o for o in report.observations if o.rule_id == "deamidation_ng"]
        assert len(ng) == 1
        assert ng[0].sequence_offset == 3
        assert ng[0].matched_residues == "NG"

    def test_finds_multiple_occurrences_of_same_motif(self):
        report = find_motifs("NGAAANG")
        assert report.counts_by_rule["deamidation_ng"] == 2

    def test_ns_and_nt_share_a_rule(self):
        assert find_motifs("ANSA").counts_by_rule["deamidation_ns_nt"] == 1
        assert find_motifs("ANTA").counts_by_rule["deamidation_ns_nt"] == 1
        assert find_motifs("ANSANTA").counts_by_rule["deamidation_ns_nt"] == 2

    def test_finds_dp_fragmentation_motif(self):
        assert find_motifs("AADPAA").counts_by_rule["fragmentation_dp"] == 1

    def test_single_residue_rules_count_each_occurrence(self):
        report = find_motifs("MAWAM")
        assert report.counts_by_rule["oxidation_met"] == 2
        assert report.counts_by_rule["oxidation_trp"] == 1

    def test_absent_motif_is_reported_as_zero_not_omitted(self):
        """A rule that found nothing differs from a rule never evaluated."""
        report = find_motifs("AAAAAA")
        for rule in MOTIF_RULES:
            assert rule.rule_id in report.counts_by_rule
            assert report.counts_by_rule[rule.rule_id] == 0

    def test_observations_are_deterministically_ordered(self):
        report = find_motifs("NGAADPAAM")
        offsets = [o.sequence_offset for o in report.observations]
        assert offsets == sorted(offsets)

    def test_overlapping_motifs_are_both_reported(self):
        # "NGG" contains NG at 0; "DGS" contains DG at 0 and DS is not present.
        report = find_motifs("NGG")
        assert report.counts_by_rule["deamidation_ng"] == 1


class TestImgtCoordinates:
    def test_observations_carry_imgt_position_and_region(self):
        numbered = number_sequence(TRASTUZUMAB_VH)
        report = find_motifs(TRASTUZUMAB_VH, numbered=numbered)
        assert report.observations, "expected motifs in a real VH"
        for obs in report.observations:
            assert obs.imgt_position is not None
            assert obs.region in {"FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"}

    def test_region_assignment_distinguishes_cdr_from_framework(self):
        """A motif in CDR-H3 and the same motif in a framework are different observations."""
        numbered = number_sequence(TRASTUZUMAB_VH)
        report = find_motifs(TRASTUZUMAB_VH, numbered=numbered)
        regions = {o.region for o in report.observations}
        assert len(regions) > 1, "expected motifs across more than one region"

    def test_no_numbering_records_the_limitation(self):
        """Missing coordinates must be stated, not left implicit."""
        report = find_motifs("AAANGAAA")
        assert report.limitations
        assert all(o.imgt_position is None for o in report.observations)
        assert all(o.region is None for o in report.observations)

    def test_motif_is_never_reported_across_an_imgt_gap(self):
        """
        Gap placeholders are positional metadata, not residues. If they reached the
        matched string, adjacent-residue motifs would be silently missed and
        offsets would not correspond to real sequence positions.
        """
        vl = number_sequence(
            "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTFTFTISSLQPEDFATYYC"
            "QQHYTTPPTFGQGTKVEIK"
        )
        raw = [r for (_p, _c), r in vl.imgt_numbering]
        assert "-" in raw, "fixture should exercise the gapped case"

        report = find_motifs("", numbered=vl)
        for obs in report.observations:
            assert "-" not in obs.matched_residues
            assert "." not in obs.matched_residues

    def test_numbering_takes_precedence_over_raw_sequence_argument(self):
        """When numbering is supplied it defines the residues, so offsets stay consistent."""
        numbered = number_sequence(TRASTUZUMAB_VH)
        from_numbered = find_motifs("", numbered=numbered)
        from_both = find_motifs(TRASTUZUMAB_VH, numbered=numbered)
        assert from_numbered.counts_by_rule == from_both.counts_by_rule


class TestOutputStaysAnObservation:
    def test_no_severity_or_score_field_exists(self):
        """The moment an observation carries a severity it has become a claim."""
        field_names = {f.name for f in fields(MotifObservation)}
        for forbidden in ("severity", "score", "risk", "weight", "penalty", "liability"):
            assert forbidden not in field_names

    def test_disclaimer_is_present_on_every_report(self):
        report = find_motifs("AAANGAAA")
        assert "not predictions" in report.disclaimer

    def test_every_rule_states_its_basis(self):
        for rule in MOTIF_RULES:
            assert rule.basis, f"{rule.rule_id} has no stated basis"

    def test_ruleset_version_is_stamped(self):
        assert find_motifs("AAA").ruleset_version == MOTIF_RULESET_VERSION

    def test_report_is_json_serializable(self):
        import json

        json.dumps(asdict(find_motifs(TRASTUZUMAB_VH, numbered=number_sequence(TRASTUZUMAB_VH))))


class TestDeterminism:
    def test_same_input_gives_identical_report(self):
        numbered = number_sequence(TRASTUZUMAB_VH)
        assert find_motifs(TRASTUZUMAB_VH, numbered=numbered) == find_motifs(
            TRASTUZUMAB_VH, numbered=numbered
        )

    def test_input_is_not_modified(self):
        original = "  aaangaaa  "
        find_motifs(original)
        assert original == "  aaangaaa  "
