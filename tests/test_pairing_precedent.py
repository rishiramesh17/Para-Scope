"""
Tests for observed pairing precedent (Version 2, Gate 8; SPEC S35.1).

Fixtures are synthetic, self-contained CSVs -- no dependency on data/, which is
gitignored and may be absent in a fresh clone.

The central test in this file is
`TestPairsAreNeverSynthesized::test_independent_heavy_and_light_hits_do_not_make_a_pair`.
Everything else is supporting. The failure mode this dimension exists to avoid is
taking a good heavy match and a good light match from two unrelated antibodies
and reporting them as pairing evidence -- which would be confident, plausible,
and entirely fabricated.
"""

import csv
from pathlib import Path

import pytest

from antigens.corpus import Corpus
from antigens.models import (
    ComponentRelationship,
    ComponentType,
    EvidenceClass,
    EvidenceStatus,
)
from antigens.numbering import compute_sequence_checksum
from antigens.pairing_precedent import (
    BLOCKED_DIMENSIONS,
    CALIBRATION_VERSION,
    PAIRING_DISCLAIMER,
    PAIRING_PRECEDENT_VERSION,
    _checksum_of,
    compute_pairing_precedent,
)
from antigens.parsers.thera_sabdab import REQUIRED_COLUMNS

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
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=sorted(REQUIRED_COLUMNS))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _row(**overrides) -> dict:
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
        "100% SI Structure": "None",
        "99% SI Structure": "None",
        "95-98% SI Structure": "None",
        "Target": "HER2",
        "Alternative Therapeutic Names": "na",
    }
    row.update(overrides)
    return row


def _corpus(tmp_path: Path, rows: list[dict]) -> Corpus:
    tmp_path.mkdir(parents=True, exist_ok=True)
    csv_path = tmp_path / "thera.csv"
    _write_thera_csv(csv_path, rows)
    corpus = Corpus()
    corpus.ingest_thera_sabdab(csv_path)
    return corpus


class TestPairsAreNeverSynthesized:
    """CLAUDE.md: never synthesize a pair by joining an independent heavy hit
    with an independent light hit. This is the whole point of the dimension."""

    def test_independent_heavy_and_light_hits_do_not_make_a_pair(self, tmp_path):
        """The fabricated combination must never be emitted.

        The corpus holds Alpha(heavy_A, light_A) and Beta(heavy_B, light_B).
        Querying heavy_A with light_B produces a perfect heavy hit and a perfect
        light hit -- both real, both exact, belonging to different antibodies.
        A naive implementation joins them and reports a confident 1.0/1.0 pair.

        Note what the correct answer is *not*: it is not an empty result. Alpha's
        own pair and Beta's own pair each half-match the query, and both are
        genuine observed pairs, so both are legitimately reported. What must never
        appear is the (heavy_A, light_B) combination that no record contains.
        """
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Alpha"),
            _row(Therapeutic="Beta", HeavySequence=OTHER_VH, LightSequence=OTHER_VL),
        ])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, OTHER_VL)

        fabricated = (
            compute_sequence_checksum(TRASTUZUMAB_VH),
            compute_sequence_checksum(OTHER_VL),
        )
        emitted = {
            (_checksum_of(m.heavy_record_id), _checksum_of(m.light_record_id))
            for m in result.matches
        }
        assert fabricated not in emitted

    def test_every_emitted_pair_exists_in_the_index(self, tmp_path):
        """The general form: nothing is reported that the corpus does not hold.

        Stronger than checking one forbidden combination, because it fails for any
        route by which a pair could be manufactured, not just the obvious one.
        """
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Alpha"),
            _row(Therapeutic="Beta", HeavySequence=OTHER_VH, LightSequence=OTHER_VL),
        ])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, OTHER_VL)

        assert result.matches  # the test is vacuous if nothing was returned
        for match in result.matches:
            linking = corpus.pairs.lookup(
                _checksum_of(match.heavy_record_id), _checksum_of(match.light_record_id)
            )
            assert linking, f"{match.heavy_record_id} + {match.light_record_id} is not an observed pair"

    def test_the_exact_pair_is_found_and_ranks_first(self, tmp_path):
        """Control: the genuine combination resolves and outranks partial matches."""
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Alpha"),
            _row(Therapeutic="Beta", HeavySequence=OTHER_VH, LightSequence=OTHER_VL),
        ])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.status is EvidenceStatus.AVAILABLE
        assert result.exact_pair_match_count == 1
        top = result.matches[0]
        assert top.heavy_global_identity == 1.0
        assert top.light_global_identity == 1.0
        assert top.source_accessions == ["Alpha#arm1"]

    def test_sharing_a_therapeutic_name_does_not_license_a_pair(self, tmp_path):
        """Entity-level identity is not a source-record link.

        Both arms belong to one bispecific therapeutic, so they share a name, a
        target and a programme. Arm 1's heavy and arm 2's light are still not a
        pair: no single record links them, and the drug's identity says nothing
        about which chains were observed together.
        """
        corpus = _corpus(tmp_path, [_row(**{
            "HeavySequence(ifbispec)": OTHER_VH,
            "LightSequence(ifbispec)": OTHER_VL,
        })])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, OTHER_VL)

        cross_arm = (
            compute_sequence_checksum(TRASTUZUMAB_VH),
            compute_sequence_checksum(OTHER_VL),
        )
        emitted = {
            (_checksum_of(m.heavy_record_id), _checksum_of(m.light_record_id))
            for m in result.matches
        }
        assert cross_arm not in emitted


class TestTwoCountsNeverOne:
    """SPEC S15: deduplicate before counting."""

    def test_one_pair_seen_in_many_records_is_one_pair(self, tmp_path):
        """Three records of the same chain pair are one pair, three observations.

        On the real corpus this ratio is ~2.1 (5,686 unique keys across 12,032
        records). Counting records as pairs would treat one antibody's crystal
        structures as independent pairing evidence.
        """
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Alpha"),
            _row(Therapeutic="Alpha-biosimilar"),
            _row(Therapeutic="Alpha-conjugate"),
        ])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.pairs_matched == 1
        assert result.raw_observation_count == 3
        assert result.matches[0].observation_count == 3

    def test_breadth_counts_pairs_not_observations(self, tmp_path):
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Alpha"),
            _row(Therapeutic="Alpha-again"),
        ])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.computable_dimensions["independent_paired_record_breadth"] == 1.0
        assert result.raw_observation_count == 2

    def test_corroborating_records_are_still_reported(self, tmp_path):
        """Deduplicating must not discard the supporting records."""
        corpus = _corpus(tmp_path, [_row(Therapeutic="Alpha"), _row(Therapeutic="Beta")])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert sorted(result.matches[0].source_accessions) == ["Alpha#arm1", "Beta#arm1"]


class TestMissingEvidenceIsNotZeroEvidence:
    def test_a_single_domain_candidate_is_not_applicable_not_not_found(self, tmp_path):
        """A VHH has no light chain by design.

        Reporting "no pairing evidence found" would describe a deliberate
        architecture as a missing result, and downstream that is indistinguishable
        from an antibody whose pair simply was not observed.
        """
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(
            corpus, NANOBODY_VHH, None, heavy_component_type=ComponentType.VHH
        )

        assert result.status is EvidenceStatus.NOT_APPLICABLE
        assert "by design" in result.reason
        assert result.matches == []

    def test_corpus_without_pairs_is_database_unavailable(self, tmp_path):
        """Not not_found: we could not look, rather than looked and saw nothing."""
        corpus = Corpus()
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.status is EvidenceStatus.DATABASE_UNAVAILABLE
        assert "not evidence" in result.reason

    def test_empty_heavy_sequence_is_unsupported(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, "", TRASTUZUMAB_VL)
        assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT

    def test_heavy_and_light_types_may_not_be_swapped(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(
            corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL,
            heavy_component_type=ComponentType.VL_KAPPA,
            light_component_type=ComponentType.VH,
        )
        assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT

    @pytest.mark.parametrize("heavy,light,expected", [
        ("", TRASTUZUMAB_VL, EvidenceStatus.UNSUPPORTED_COMPONENT),
        (TRASTUZUMAB_VH, None, EvidenceStatus.NOT_APPLICABLE),
        # A light chain too short to clear the subject-coverage floor: the heavy
        # matches perfectly and there is still no admissible pair.
        (TRASTUZUMAB_VH, TRASTUZUMAB_VL[:30], EvidenceStatus.NOT_FOUND),
    ])
    def test_no_failure_path_returns_a_neutral_number(self, tmp_path, heavy, light, expected):
        """Every non-available outcome leaves the measurements unset, never 0.0."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, heavy, light)

        assert result.status is expected
        assert result.best_heavy_identity is None
        assert result.best_light_identity is None
        assert result.best_joint_coverage is None
        assert result.computable_dimensions == {}


class TestNoCompositeIsInvented:
    """SPEC S35 lists dimensions but specifies no weights."""

    def test_result_exposes_no_composite_score_field(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert not hasattr(result, "score")
        assert result.calibration_version == CALIBRATION_VERSION == "uncalibrated"

    def test_blocked_dimensions_name_themselves_and_say_why(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert set(result.blocked_dimensions) == {
            "vh_vl_germline_pair_precedent",
            "combined_cdr_pattern_similarity",
            "format_relevance",
        }
        assert all(reason for reason in result.blocked_dimensions.values())

    def test_blocked_dimensions_are_reported_even_when_nothing_matched(self, tmp_path):
        """A future session must see what is missing regardless of the outcome."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL[:30])

        assert result.status is EvidenceStatus.NOT_FOUND
        assert len(result.blocked_dimensions) == 3

    def test_the_five_computable_dimensions_are_raw(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert set(result.computable_dimensions) == {
            "heavy_similarity",
            "light_similarity",
            "joint_coverage",
            "direct_parent_record_evidence",
            "independent_paired_record_breadth",
        }
        assert result.computable_dimensions["heavy_similarity"] == 1.0
        assert result.computable_dimensions["light_similarity"] == 1.0


class TestJointCoverage:
    def test_joint_coverage_is_the_weaker_chain_not_the_average(self, tmp_path):
        """A fully-covered heavy must not hide a poorly-covered light."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)
        match = result.matches[0]

        assert match.joint_coverage == min(
            match.heavy_query_coverage, match.light_query_coverage
        )


class TestStructureSupport:
    def test_unknown_structure_support_is_not_false(self, tmp_path):
        """Optional[bool] throughout: None means the source could not tell us."""
        corpus = _corpus(tmp_path, [_row(**{
            "HeavySequence(ifbispec)": OTHER_VH,
            "LightSequence(ifbispec)": OTHER_VL,
            "100% SI Structure": "a;b;c",   # three tokens, two arms -> unmappable
            "99% SI Structure": "a;b;c",
            "95-98% SI Structure": "a;b;c",
        })])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.matches[0].structure_supported is None
        assert result.computable_dimensions["direct_parent_record_evidence"] is None

    def test_stated_absence_is_false_not_unknown(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.matches[0].structure_supported is False
        assert result.computable_dimensions["direct_parent_record_evidence"] == 0.0

    def test_a_real_structure_is_reported_with_its_ids(self, tmp_path):
        corpus = _corpus(tmp_path, [_row(**{"100% SI Structure": "1n8z:BA"})])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.matches[0].structure_supported is True
        assert "1n8z:BA" in result.matches[0].structure_ids


class TestLikeWithLike:
    def test_a_pair_is_only_as_same_class_as_its_weaker_chain(self, tmp_path):
        """Kappa/lambda substitution downgrades the pair, not just the light chain.

        Otherwise an exact heavy match would carry a related-class light match
        into same-class evidence, which CLAUDE.md forbids merging.
        """
        from antigens.pairing_precedent import _worst_relationship

        assert _worst_relationship(
            ComponentRelationship.EXACT_SAME_TYPE,
            ComponentRelationship.RELATED_CLASS,
        ) is ComponentRelationship.RELATED_CLASS

    def test_related_class_pairs_rank_below_same_class_at_equal_identity(self, tmp_path):
        from antigens.pairing_precedent import PairMatch

        def _match(relationship):
            return PairMatch(
                heavy_record_id="h", light_record_id="l",
                heavy_component_type=ComponentType.VH,
                light_component_type=ComponentType.VL_KAPPA,
                pair_relationship=relationship,
                heavy_global_identity=1.0, light_global_identity=1.0,
                heavy_query_coverage=1.0, light_query_coverage=1.0,
                joint_coverage=1.0,
                heavy_match_method="exact_checksum",
                light_match_method="exact_checksum",
            )

        same = _match(ComponentRelationship.EXACT_SAME_TYPE)
        related = _match(ComponentRelationship.RELATED_CLASS)
        assert same.rank_key > related.rank_key


class TestProvenanceAndFraming:
    def test_disclaimer_is_carried_verbatim_on_every_result(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        for heavy, light in [
            (TRASTUZUMAB_VH, TRASTUZUMAB_VL),
            ("", TRASTUZUMAB_VL),
            (TRASTUZUMAB_VH, None),
            (OTHER_VH, OTHER_VL),
        ]:
            result = compute_pairing_precedent(corpus, heavy, light)
            assert result.disclaimer == PAIRING_DISCLAIMER

    def test_disclaimer_denies_the_prediction_reading(self, tmp_path):
        assert "not a prediction" in PAIRING_DISCLAIMER
        assert "absence of an observed pair is not evidence" in PAIRING_DISCLAIMER

    def test_result_states_its_engine_and_calibration_versions(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)

        assert result.engine_version == PAIRING_PRECEDENT_VERSION
        assert result.calibration_version == "uncalibrated"
        assert result.evidence_class is EvidenceClass.DIRECT_PAIRING

    def test_match_names_the_source_records_that_link_the_chains(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_pairing_precedent(corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL)
        match = result.matches[0]

        assert match.source_databases == ["Thera-SAbDab"]
        assert match.source_accessions == ["Testuzumab#arm1"]
        assert match.observation_count == len(match.source_records)

    def test_module_exposes_no_prohibited_claim_vocabulary(self):
        import antigens.pairing_precedent as module

        banned = ("compatib", "will_pair", "predicted_pair", "stability", "affinity")
        offenders = [
            name for name in dir(module)
            if not name.startswith("_") and any(word in name.lower() for word in banned)
        ]
        assert offenders == []

    def test_35_2_interface_prediction_is_absent(self):
        """SPEC S35.2 defers predicted VH-VL interface compatibility entirely."""
        import antigens.pairing_precedent as module

        source = Path(module.__file__).read_text()
        for banned in ("orientation", "buried_surface", "steric"):
            assert f"def {banned}" not in source
