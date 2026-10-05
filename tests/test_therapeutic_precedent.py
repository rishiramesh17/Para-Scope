"""
Tests for the therapeutic-precedent dimension (Version 2, Gate 7; SPEC S26).

Fixtures are synthetic, self-contained CSVs -- no dependency on data/, which is
gitignored and may be absent in a fresh clone.

The bulk of these tests assert things the code must *refuse* to do: fabricate a
composite score from uncalibrated inputs, merge related-class into same-class
evidence, count one therapeutic several times, or report a failed lookup as a
low score. Those are the failure modes that produce a confident wrong answer,
which is worse here than producing none.
"""

import pytest

import csv
from pathlib import Path

from para_scope.corpus import Corpus, split_alternative_names
from para_scope.models import (
    ComponentRelationship,
    ComponentType,
    EvidenceClass,
    EvidenceStatus,
)
from para_scope.numbering import compute_sequence_checksum, number_sequence
from para_scope.parsers.thera_sabdab import REQUIRED_COLUMNS
from para_scope.therapeutic_precedent import (
    BLOCKED_TERMS,
    CALIBRATION_VERSION,
    THERAPEUTIC_DISCLAIMER,
    THERAPEUTIC_PRECEDENT_VERSION,
    compute_therapeutic_precedent,
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
# IMGT CDR3 opens at position 105, so trastuzumab's CDR-H3 begins with SR.
TRASTUZUMAB_CDR_H3 = "SRWGGDGFYAMDY"


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


class TestTherapeuticEntityStore:
    def test_metadata_is_reachable_for_an_arm_with_no_light_chain(self, tmp_path):
        """The gap that motivated the store.

        Therapeutic metadata used to live only on PairRecord, so a single-domain
        therapeutic -- which by construction has no pair -- carried no recoverable
        name, target, or development status at all. SPEC S26 requires those for any
        supported component.
        """
        corpus = _corpus(tmp_path, [_row(
            Therapeutic="Nanobodimab", Format="Nanobody",
            HeavySequence=NANOBODY_VHH, LightSequence="na", Target="VWF",
        )])
        record = corpus.get_therapeutic("Nanobodimab")
        assert record is not None
        assert record.target == "VWF"
        assert record.arms[0].light_record_id is None
        assert record.arms[0].heavy_record_id is not None

    def test_component_resolves_back_to_its_therapeutics(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        record_id = f"{ComponentType.VH.value}::{number_sequence(TRASTUZUMAB_VH).sequence_checksum}"
        names = [t.therapeutic_name for t in corpus.therapeutics_for_component(record_id)]
        assert names == ["Testuzumab"]

    def test_sub_regions_link_to_the_therapeutic_too(self, tmp_path):
        """A CDR must resolve to its therapeutic, not just the whole variable domain."""
        corpus = _corpus(tmp_path, [_row()])
        # Numbering a bare CDR fails (it is not a chain), so the record id comes from
        # the sequence checksum directly -- the same key the corpus stored it under.
        cdr_id = f"{ComponentType.CDR_H3.value}::{compute_sequence_checksum(TRASTUZUMAB_CDR_H3)}"
        assert [t.therapeutic_name for t in corpus.therapeutics_for_component(cdr_id)] == ["Testuzumab"]

    def test_alias_and_development_code_resolve_to_one_entity(self, tmp_path):
        """SPEC S15.5: names, brand names, and development codes are one entity."""
        corpus = _corpus(tmp_path, [_row(**{
            "Alternative Therapeutic Names": "Brandcept;XYZ-123;XYZ 123",
        })])
        for probe in ("Testuzumab", "testuzumab", "Brandcept", "XYZ-123", "xyz 123"):
            assert corpus.get_therapeutic(probe).therapeutic_name == "Testuzumab"

    def test_unknown_name_resolves_to_nothing_rather_than_a_guess(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        assert corpus.get_therapeutic("Nonexistumab") is None

    def test_alias_claimed_by_two_therapeutics_is_left_unresolved(self, tmp_path):
        """Assigning it to either would silently merge two distinct entities."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [
            _row(Therapeutic="Alphamab", **{"Alternative Therapeutic Names": "SHARED-1"}),
            _row(Therapeutic="Betamab", HeavySequence=OTHER_VH, LightSequence=OTHER_VL,
                 **{"Alternative Therapeutic Names": "SHARED-1"}),
        ])
        corpus = Corpus()
        stats = corpus.ingest_thera_sabdab(csv_path)

        assert stats.therapeutic_alias_collisions == 1
        assert corpus.get_therapeutic("SHARED-1") is None
        # The canonical names themselves are untouched by the collision.
        assert corpus.get_therapeutic("Alphamab").therapeutic_name == "Alphamab"
        assert corpus.get_therapeutic("Betamab").therapeutic_name == "Betamab"

    def test_a_third_claimant_cannot_reclaim_an_ambiguous_alias(self, tmp_path):
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [
            _row(Therapeutic="Alphamab", **{"Alternative Therapeutic Names": "SHARED-1"}),
            _row(Therapeutic="Betamab", HeavySequence=OTHER_VH, LightSequence=OTHER_VL,
                 **{"Alternative Therapeutic Names": "SHARED-1"}),
            _row(Therapeutic="Gammamab", HeavySequence=NANOBODY_VHH, LightSequence="na",
                 Format="Nanobody", **{"Alternative Therapeutic Names": "SHARED-1"}),
        ])
        corpus = Corpus()
        corpus.ingest_thera_sabdab(csv_path)
        assert corpus.get_therapeutic("SHARED-1") is None


class TestAlternativeNameSplitting:
    def test_plain_semicolon_list_splits(self):
        assert split_alternative_names("ACA 125;ACA125;ACA-125") == [
            "ACA 125", "ACA125", "ACA-125",
        ]

    def test_semicolons_inside_parentheses_do_not_delimit(self):
        """Regression: the naive split fragmented combination-product names.

        'REGN5713-5714-5715 (combination of Atisnolerbart; Bremzalerbart; and
        Umesolerbart)' yielded a bare 'Bremzalerbart', which then collided with the
        real therapeutic of that name and left a legitimate alias unresolvable --
        plus 'and Umesolerbart)', which is not a name.
        """
        raw = "REGN5713-5714-5715 (combination of Atisnolerbart; Bremzalerbart; and Umesolerbart)"
        assert split_alternative_names(raw) == [raw]

    def test_placeholders_are_not_registered_as_names(self):
        """'TBC' was claimed by two therapeutics and poisoned itself as an alias."""
        assert split_alternative_names("TBC") == []
        assert split_alternative_names("na;Unknown;ND") == []
        assert split_alternative_names("Realname;TBC") == ["Realname"]

    def test_empty_input_yields_no_names(self):
        assert split_alternative_names(None) == []
        assert split_alternative_names("") == []

    def test_combination_identifier_stays_ambiguous_by_design(self, tmp_path):
        """A name denoting two therapeutics must resolve to neither (SPEC S15.5)."""
        csv_path = tmp_path / "thera.csv"
        _write_thera_csv(csv_path, [
            _row(Therapeutic="Cilgavimab",
                 **{"Alternative Therapeutic Names": "Evusheld (combination of cilgavimab and tixagevimab)"}),
            _row(Therapeutic="Tixagevimab", HeavySequence=OTHER_VH, LightSequence=OTHER_VL,
                 **{"Alternative Therapeutic Names": "Evusheld (combination of cilgavimab and tixagevimab)"}),
        ])
        corpus = Corpus()
        corpus.ingest_thera_sabdab(csv_path)

        assert corpus.get_therapeutic("Evusheld (combination of cilgavimab and tixagevimab)") is None
        assert corpus.get_therapeutic("Cilgavimab").therapeutic_name == "Cilgavimab"


class TestMissingEvidenceIsNotZeroEvidence:
    def test_no_match_is_not_found_and_leaves_the_score_unset(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.CDR_H3, "WWWWWWWWWWWWW")

        assert result.status is EvidenceStatus.NOT_FOUND
        assert result.score is None
        assert result.best_global_identity is None
        assert result.therapeutics_matched == 0

    def test_corpus_without_therapeutics_is_database_unavailable_not_not_found(self, tmp_path):
        """"We have no therapeutic records" and "we looked and found none" differ.

        Reporting the first as the second would let an absent database read as
        evidence that no therapeutic precedent exists.
        """
        result = compute_therapeutic_precedent(Corpus(), ComponentType.VH, TRASTUZUMAB_VH)

        assert result.status is EvidenceStatus.DATABASE_UNAVAILABLE
        assert result.score is None
        assert "not evidence of absence" in result.reason

    def test_empty_sequence_is_unsupported_not_a_zero_score(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, "")
        assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT
        assert result.score is None

    def test_no_failure_path_returns_a_neutral_number(self, tmp_path):
        """Every non-AVAILABLE outcome leaves score None -- never 0.0, 0.5, or a mean."""
        corpus = _corpus(tmp_path, [_row()])
        outcomes = [
            compute_therapeutic_precedent(corpus, ComponentType.VH, ""),
            compute_therapeutic_precedent(corpus, ComponentType.CDR_H3, "WWWWWWWWWWWWW"),
            compute_therapeutic_precedent(Corpus(), ComponentType.VH, TRASTUZUMAB_VH),
        ]
        for result in outcomes:
            assert result.status is not EvidenceStatus.AVAILABLE
            assert result.score is None


class TestCompositeScoreIsWithheld:
    def test_score_is_never_emitted_while_calibration_is_missing(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert result.status is EvidenceStatus.AVAILABLE  # evidence exists...
        assert result.score is None                       # ...but no composite
        assert result.score_status is EvidenceStatus.INSUFFICIENT_EVIDENCE
        assert result.calibration_version == CALIBRATION_VERSION == "uncalibrated"

    def test_blocked_terms_name_themselves_and_say_why(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert set(result.blocked_terms) == {
            "calibrated_component_similarity", "format_relevance",
        }
        assert all(reason for reason in result.blocked_terms.values())

    def test_the_four_computable_terms_are_reported_raw(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert set(result.computable_terms) == {
            "coverage",
            "exact_component_class_confirmation",
            "direct_pair_or_structure_support",
            "metadata_completeness",
        }
        assert result.computable_terms["coverage"] == 1.0
        assert result.computable_terms["exact_component_class_confirmation"] == 1.0

    def test_gate_9_alone_does_not_unblock_the_composite(self, tmp_path):
        """Similarity is 0.45 of SPEC S26; format relevance is a further 0.10.

        This test exists because an earlier version of it was named
        "..._so_gate_9_need_only_add_similarity", which stated the opposite and
        would have been read by a future session as a settled decision. Filling
        only the similarity slot leaves 0.90 of available weight, and emitting a
        composite over that is the silent renormalization CLAUDE.md forbids: the
        omitted term self-normalizes and the index reads as complete.

        Two terms are blocked, so two must be unblocked -- or the composite must
        declare the scope it actually covered.
        """
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert len(result.blocked_terms) == 2
        assert "format_relevance" in result.blocked_terms
        # The withheld reason must name both, not just the calibrated one.
        assert "format relevance" in result.score_reason.lower()
        assert "calibrated component similarity" in result.score_reason.lower()

    def test_weights_are_not_applied_to_the_partial_terms(self, tmp_path):
        """A weighted partial sum would look like a composite that merely scored low."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)
        # SPEC S26 weights coverage at 0.15; the raw proportion must survive intact.
        assert result.computable_terms["coverage"] == 1.0
        assert BLOCKED_TERMS["calibrated_component_similarity"]


class TestLikeWithLike:
    def test_exact_match_of_a_shared_sequence_never_crosses_chain_types(self, tmp_path):
        """The exact index is keyed by checksum alone.

        A short region can appear verbatim in both a heavy and a light chain, so a
        raw lookup can return a CDR-L1 for a CDR-H1 query. Reporting that at
        identity 1.0 would be a textbook like-with-like violation.
        """
        corpus = _corpus(tmp_path, [_row()])
        heavy_fr4 = None
        for record_id, record in corpus._components.items():
            if record.component_type is ComponentType.FR_H4:
                heavy_fr4 = record
                break
        assert heavy_fr4 is not None

        result = compute_therapeutic_precedent(
            corpus, ComponentType.FR_L4, heavy_fr4.sequence
        )
        for match in result.matches:
            assert match.matched_component_type is not ComponentType.FR_H4

    def test_related_class_matches_are_labelled_and_ranked_below_same_class(self, tmp_path):
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Sameclass"),
            _row(Therapeutic="Relatedclass", Format="Nanobody",
                 HeavySequence=NANOBODY_VHH, LightSequence="na"),
        ])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert result.status is EvidenceStatus.AVAILABLE
        assert result.matches[0].component_relationship is ComponentRelationship.EXACT_SAME_TYPE
        related = [
            m for m in result.matches
            if m.component_relationship is ComponentRelationship.RELATED_CLASS
        ]
        if related:
            assert any("related-class" in note for note in result.limitations)
            assert related[0].rank_key < result.matches[0].rank_key


class TestDeduplication:
    def test_entities_and_observations_are_counted_separately(self, tmp_path):
        """SPEC S15: one therapeutic matching twice is one precedent, two observations."""
        corpus = _corpus(tmp_path, [_row(
            Therapeutic="Bispecifab", Format="Bispecific mAb",
            **{"HeavySequence(ifbispec)": OTHER_VH, "LightSequence(ifbispec)": OTHER_VL},
        )])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert result.therapeutics_matched == 1
        assert result.raw_match_count >= result.therapeutics_matched

    def test_two_therapeutics_sharing_one_sequence_count_as_two(self, tmp_path):
        """Deduplicating the component must not deduplicate the therapeutics behind it."""
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Firstumab"),
            _row(Therapeutic="Secondumab"),  # identical sequences, distinct programs
        ])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        assert result.therapeutics_matched == 2
        assert {m.therapeutic_name for m in result.matches} == {"Firstumab", "Secondumab"}

    def test_match_order_is_deterministic_under_ties(self, tmp_path):
        """Identical sequences across therapeutics are the normal case, not an edge one."""
        corpus = _corpus(tmp_path, [
            _row(Therapeutic="Zulumab"),
            _row(Therapeutic="Alphamab"),
            _row(Therapeutic="Mikemab"),
        ])
        orders = [
            [m.therapeutic_name for m in
             compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH).matches]
            for _ in range(3)
        ]
        assert orders[0] == orders[1] == orders[2]
        assert orders[0][:3] == ["Alphamab", "Mikemab", "Zulumab"]


class TestShortComponents:
    def test_short_component_reports_exact_hits_and_declines_alignment(self, tmp_path):
        """A CDR-H3 identical to a therapeutic's is evidence; there is nothing to calibrate."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(
            corpus, ComponentType.CDR_H3, TRASTUZUMAB_CDR_H3
        )

        assert result.status is EvidenceStatus.AVAILABLE
        assert result.exact_match_count == 1
        assert result.matches[0].match_method == "exact_checksum"
        assert result.approximate_search_status is EvidenceStatus.UNSUPPORTED_COMPONENT
        assert "S21" in result.approximate_search_reason

    def test_declining_alignment_is_recorded_as_a_limitation(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(
            corpus, ComponentType.CDR_H3, TRASTUZUMAB_CDR_H3
        )
        assert any("Alignment search not performed" in note for note in result.limitations)


class TestMatchCarriesItsSource:
    """Where a match came from must be readable from the match itself.

    The corpus has always held this on the component record, but reaching it meant
    taking matched_record_id back to corpus.get(...).source_pointers. That hop is
    easy to get wrong in a way that fails silently -- the stored source name is
    "Thera-SAbDab", and filtering on "thera_sabdab" returns zero matches rather
    than raising. SPEC S13/S14 also require a source's licence to travel with the
    evidence derived from it, and licence is recorded per source.
    """

    def test_match_names_the_database_it_was_observed_in(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        match = result.matches[0]
        assert match.source_databases == ["Thera-SAbDab"]
        assert match.source_pointers
        assert all(
            database and accession for database, accession in match.source_pointers
        )

    def test_source_is_not_re_derivable_only_by_a_second_lookup(self, tmp_path):
        """The whole point: no corpus round-trip needed to answer 'from where?'."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)

        match = result.matches[0]
        via_match = set(match.source_pointers)
        via_corpus = set(corpus.get(match.matched_record_id).source_pointers)
        assert via_match == via_corpus

    def test_a_component_in_two_sources_reports_both(self, tmp_path):
        """Provenance is every source the component was seen in, not the first one.

        A VH present in SAbDab2 and reached here through a Thera-SAbDab therapeutic
        link genuinely has two pointers. Reporting one would misstate the record.
        """
        corpus = _corpus(tmp_path, [_row()])
        record_id = f"{ComponentType.VH.value}::{compute_sequence_checksum(TRASTUZUMAB_VH)}"
        corpus.get(record_id).source_pointers.append(("SAbDab2", "1n8z_H"))

        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)
        assert result.matches[0].source_databases == ["SAbDab2", "Thera-SAbDab"]

    def test_engine_version_records_that_provenance_is_present(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)
        assert result.engine_version == THERAPEUTIC_PRECEDENT_VERSION == "1.1.0"


class TestProvenanceAndFraming:
    def test_disclaimer_is_carried_verbatim_on_every_result(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        for component_type, query in [
            (ComponentType.VH, TRASTUZUMAB_VH),
            (ComponentType.VH, ""),
            (ComponentType.CDR_H3, "WWWWWWWWWWWWW"),
        ]:
            result = compute_therapeutic_precedent(corpus, component_type, query)
            assert result.disclaimer == THERAPEUTIC_DISCLAIMER
            assert "does not transfer approval" in result.disclaimer

    def test_result_states_its_engine_and_calibration_versions(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)
        assert result.engine_version == THERAPEUTIC_PRECEDENT_VERSION
        assert result.calibration_version == CALIBRATION_VERSION
        assert result.evidence_class is EvidenceClass.THERAPEUTIC_PRECEDENT

    def test_clinical_stage_is_carried_as_metadata_not_as_a_verdict(self, tmp_path):
        corpus = _corpus(tmp_path, [_row(**{"Highest_Clin_Trial (Feb '25)": "Phase-III"})])
        result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH)
        match = result.matches[0]

        assert match.highest_clinical_trial == "Phase-III"
        # It must not have been folded into any numeric term.
        assert "clinical" not in " ".join(result.computable_terms)
        assert "stage" not in " ".join(result.computable_terms)

    def test_discontinued_status_does_not_lower_any_reported_measurement(self, tmp_path):
        """A failed clinical programme says nothing about sequence similarity."""
        active = _corpus(tmp_path / "a", [_row(**{"Est. Status": "Active"})])
        (tmp_path / "b").mkdir(parents=True, exist_ok=True)
        stopped = _corpus(tmp_path / "b", [_row(**{"Est. Status": "Discontinued"})])

        result_active = compute_therapeutic_precedent(active, ComponentType.VH, TRASTUZUMAB_VH)
        result_stopped = compute_therapeutic_precedent(stopped, ComponentType.VH, TRASTUZUMAB_VH)

        assert result_active.best_global_identity == result_stopped.best_global_identity
        assert result_active.computable_terms == result_stopped.computable_terms
        assert result_stopped.matches[0].development_status == "Discontinued"

    def test_module_exposes_no_prohibited_claim_vocabulary(self):
        """The dimension names comparisons, never conclusions drawn from them."""
        import para_scope.therapeutic_precedent as module

        banned = ("efficacy", "safety_score", "will_bind", "success", "approved_score")
        for name in dir(module):
            if name.startswith("_"):
                continue
            assert not any(word in name.lower() for word in banned), name

    def test_structure_support_reflects_the_source_not_an_assumption(self, tmp_path):
        """The term this gate's parser fix exists to protect."""
        without = _corpus(tmp_path / "n", [_row()])
        (tmp_path / "y").mkdir(parents=True, exist_ok=True)
        with_structure = _corpus(
            tmp_path / "y", [_row(**{"100% SI Structure": "1n8z:BA"})]
        )

        no_struct = compute_therapeutic_precedent(without, ComponentType.VH, TRASTUZUMAB_VH)
        has_struct = compute_therapeutic_precedent(with_structure, ComponentType.VH, TRASTUZUMAB_VH)

        assert no_struct.matches[0].has_structure_record is False
        assert no_struct.matches[0].structure_ids == []
        assert has_struct.matches[0].has_structure_record is True
        assert has_struct.matches[0].structure_ids == ["1n8z:BA"]

class TestNonSequenceInputIsUnsupportedNotNotFound:
    """
    Regression, found by the adversarial audit (reports/ADVERSARIAL_AUDIT.md).

    NOT_FOUND asserts "the corpus holds no precedent for this" -- a claim about the
    database. For a string that is not a sequence the true statement is that it was
    never a valid query. Five of eight degenerate inputs returned NOT_FOUND before
    this check existed, and germline.py already validated the other way, so the
    two modules disagreed about the same input.
    """

    @pytest.mark.parametrize("sequence,label", [
        ("1234567890123456789012345678901234567890", "digits"),
        ("not a protein sequence at all", "prose"),
        ("EVQLVESGGG-LVQPGGSLRLSCAAS", "gap character"),
        ("EVQLVESGGG*LVQPGGSLRLSCAAS", "stop codon character"),
        ("EVQLXESGGGLVQPGGSLRLSCAASBFNIKDTYIHZVRQ", "ambiguity codes"),
    ])
    def test_non_sequence_input_is_unsupported_component(self, tmp_path, sequence, label):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(
            corpus, ComponentType.VH, sequence, shortlist_limit=None,
        )
        assert result.status is EvidenceStatus.UNSUPPORTED_COMPONENT, label
        assert result.reason is not None
        assert result.matches == []

    def test_the_reason_names_the_offending_characters(self, tmp_path):
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(
            corpus, ComponentType.VH, "EVQLVESGGG-LVQ", shortlist_limit=None,
        )
        assert "-" in result.reason

    def test_input_is_not_silently_repaired_to_make_it_searchable(self, tmp_path):
        """Stripping the offending character would search a sequence nobody submitted."""
        corpus = _corpus(tmp_path, [_row()])
        result = compute_therapeutic_precedent(
            corpus, ComponentType.VH, "EVQLVESGGG-LVQ", shortlist_limit=None,
        )
        assert result.matches == []
        assert "not modified" in result.reason
