"""
Tests for human germline V-gene assignment (Gate 9 milestone 1).

Uses a small synthetic reference index built from a handful of real IGHV /
IGKV / IGLV sequences, injected via `reference_index=`, rather than the
bundled OGRDB file -- unit tests stay fast and self-contained, and the
integration test (tests/integration/test_germline_snapshot.py) is what
exercises the real, checksum-pinned bundled reference set.
"""

from para_scope.germline import (
    GERMLINE_VERSION,
    CALIBRATION_VERSION,
    GERMLINE_DISCLAIMER,
    compute_germline_assignment,
)
from para_scope.indexes.approximate import ApproximateIndex
from para_scope.models import ComponentType, EvidenceClass, EvidenceStatus

# Real germline V-gene protein sequences (see src/para_scope/reference_data/germline).
IGHV3_23_01 = (
    "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTIS"
    "RDNSKNTLYLQMNSLRAEDTAVYYCAK"
)
IGKV1_39_01 = (
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTD"
    "FTLTISSLQPEDFATYYCQQSYSTP"
)
IGLV1_40_01 = (
    "QSVLTQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSKSG"
    "TSASLAITGLQAEDEADYYCQSYDSSLSG"
)


def _reference_index() -> ApproximateIndex:
    index = ApproximateIndex()
    index.add(ComponentType.VH, "IGHV3-23*01|IGH", IGHV3_23_01)
    index.add(ComponentType.VL_KAPPA, "IGKV1-39*01|IGK", IGKV1_39_01)
    index.add(ComponentType.VL_LAMBDA, "IGLV1-40*01|IGL", IGLV1_40_01)
    return index


class TestExactMatch:
    def test_identical_vh_sequence_matches_its_germline(self):
        result = compute_germline_assignment(
            ComponentType.VH, IGHV3_23_01, reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.AVAILABLE
        assert result.nearest_gene_call == "IGHV3-23"
        assert result.nearest_germline_identity == 1.0
        assert result.matches[0].allele == "IGHV3-23*01"
        assert result.matches[0].locus == "IGH"
        assert result.evidence_class == EvidenceClass.GERMLINE
        assert result.engine_version == GERMLINE_VERSION
        assert result.calibration_version == CALIBRATION_VERSION
        assert result.disclaimer == GERMLINE_DISCLAIMER

    def test_identical_vl_kappa_sequence_matches_its_germline(self):
        result = compute_germline_assignment(
            ComponentType.VL_KAPPA, IGKV1_39_01, reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.AVAILABLE
        assert result.nearest_gene_call == "IGKV1-39"

    def test_identical_vl_lambda_sequence_matches_its_germline(self):
        result = compute_germline_assignment(
            ComponentType.VL_LAMBDA, IGLV1_40_01, reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.AVAILABLE
        assert result.nearest_gene_call == "IGLV1-40"

    def test_kappa_query_is_never_scored_against_a_lambda_only_index(self):
        # Like-with-like: a kappa query must search only kappa-typed entries.
        # A reference index holding no VL_KAPPA entries at all must return no
        # hits, not silently borrow from the VL_LAMBDA bucket.
        lambda_only = ApproximateIndex()
        lambda_only.add(ComponentType.VL_LAMBDA, "IGLV1-40*01|IGL", IGLV1_40_01)
        result = compute_germline_assignment(
            ComponentType.VL_KAPPA, IGKV1_39_01, reference_index=lambda_only
        )
        assert result.status == EvidenceStatus.NOT_FOUND


class TestMissingEvidenceIsNotZero:
    def test_divergent_query_below_coverage_floor_is_not_found_not_scored_low(self):
        # A short, unrelated peptide should never receive a fabricated low
        # identity -- it must be dropped from the result entirely.
        result = compute_germline_assignment(
            ComponentType.VH, "GGGGGGGGGGGGGGGG", reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.NOT_FOUND
        assert result.nearest_germline_identity is None
        assert result.matches == []
        assert result.reason is not None

    def test_empty_sequence_is_unsupported_not_a_score(self):
        result = compute_germline_assignment(
            ComponentType.VH, "", reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.UNSUPPORTED_COMPONENT
        assert result.nearest_germline_identity is None

    def test_non_standard_residues_are_unsupported_and_input_is_not_modified(self):
        query = "EVQLLXSGGGLVQPGGSLRLSCAAS"
        result = compute_germline_assignment(
            ComponentType.VH, query, reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.UNSUPPORTED_COMPONENT
        assert "X" in result.reason
        # The module must never substitute, strip, or otherwise touch the
        # input to make it scoreable (CLAUDE.md).
        assert query == "EVQLLXSGGGLVQPGGSLRLSCAAS"


class TestVHHIsExplicitlyExcluded:
    def test_vhh_is_unsupported_not_scored_against_human_reference(self):
        # Human-only scope, deliberately: there is no camelid germline set
        # bundled, so a VHH must never be silently aligned against human VH
        # germlines and reported as if it were a real assignment.
        result = compute_germline_assignment(
            ComponentType.VHH, IGHV3_23_01, reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.UNSUPPORTED_COMPONENT
        assert "VHH" in result.reason

    def test_short_components_are_unsupported_not_approximated(self):
        # CDRs/frameworks are not compared against a whole V gene here.
        result = compute_germline_assignment(
            ComponentType.CDR_H3, "ARWGGDGFYAMDY", reference_index=_reference_index()
        )
        assert result.status == EvidenceStatus.UNSUPPORTED_COMPONENT


class TestTopKAndRanking:
    def test_top_k_limits_returned_matches(self):
        result = compute_germline_assignment(
            ComponentType.VH, IGHV3_23_01, top_k=1, reference_index=_reference_index()
        )
        assert len(result.matches) <= 1

    def test_matches_are_never_from_the_wrong_component_type(self):
        result = compute_germline_assignment(
            ComponentType.VH, IGHV3_23_01, top_k=10, reference_index=_reference_index()
        )
        for match in result.matches:
            assert match.locus == "IGH"
