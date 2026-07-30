"""
Tests for IMGT numbering and component comparison eligibility.

Key test: verify that CDR-H3 never matches against CDR-L3 (SPEC §20, CLAUDE.md).
"""

import pytest
from antigens.numbering import (
    IMGT_GAP_CHARACTERS,
    extract_cdr_h3,
    extract_cdr_l3,
    extract_region_sequence,
    number_sequence,
)
from antigens.models import ComponentType, compare_eligibility, ComponentRelationship


# Public trastuzumab heavy-chain variable domain
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)

# Public trastuzumab light-chain variable domain (kappa)
# A complete kappa V-domain, through FR-L4. It previously stopped at C104, so it
# covered no CDR-L3 or FR-L4 positions at all -- and the CDR-L3 assertions below
# passed only because region extraction returned a string of IMGT gap
# placeholders, which is truthy but contains no residues.
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTFTFTISSLQPEDFATYYC"
    "QQHYTTPPTFGQGTKVEIK"
)


class TestNumbering:
    """Test ANARCII numbering wrapper."""

    def test_number_vh(self):
        """Trastuzumab VH should number successfully as heavy chain."""
        result = number_sequence(TRASTUZUMAB_VH)
        assert result is not None
        assert result.chain_type == 'H'
        assert result.scheme == 'imgt'
        assert result.numbering_error is None
        assert result.numbering_score > 0
        assert len(result.imgt_numbering) > 100  # Full VH is ~128 positions

    def test_number_vl(self):
        """Trastuzumab VL (kappa) should number successfully as light chain."""
        result = number_sequence(TRASTUZUMAB_VL)
        assert result is not None
        assert result.chain_type == 'K'  # Kappa light chain
        assert result.scheme == 'imgt'
        assert result.numbering_error is None

    def test_extract_cdr_h3(self):
        """CDR-H3 extraction from a numbered heavy chain."""
        vh = number_sequence(TRASTUZUMAB_VH)
        assert vh is not None
        cdr_h3 = extract_cdr_h3(vh)
        assert cdr_h3 is not None
        assert len(cdr_h3) > 0
        # Trastuzumab CDR-H3 is well-known: WGGDGFYAMDY
        assert 'WGGDG' in cdr_h3 or 'WGGDGFYA' in cdr_h3

    def test_extract_cdr_l3(self):
        """CDR-L3 extraction from a numbered light chain."""
        vl = number_sequence(TRASTUZUMAB_VL)
        assert vl is not None
        cdr_l3 = extract_cdr_l3(vl)
        assert cdr_l3 is not None
        assert len(cdr_l3) > 0
        assert cdr_l3 == "QQHYTTPPT"


class TestImgtGapHandling:
    """
    IMGT is fixed-width: unoccupied position slots are emitted as '-'/'.'
    placeholders. Those are positional metadata, not sequence content.

    Regression: they were being joined into region sequences, which (a) crashed
    every framework alignment, since BLOSUM62 has no such letter, (b) were hashed
    into component checksums, splitting identical components into several, and
    (c) let a region covering zero real residues be indexed as a component made
    entirely of gaps.
    """

    def test_region_sequence_excludes_gap_placeholders(self):
        vl = number_sequence(TRASTUZUMAB_VL)
        for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"):
            seq = extract_region_sequence(vl, region)
            assert not (set(seq) & IMGT_GAP_CHARACTERS), f"{region} retained a gap placeholder"

    def test_gapped_cdr_l3_is_contiguous_residues_only(self):
        """The real trastuzumab CDR-L3 spans a gapped IMGT range (110-113 empty)."""
        vl = number_sequence(TRASTUZUMAB_VL)
        raw_slots = [r for (pos, _c), r in vl.imgt_numbering if 105 <= pos <= 117]
        assert "-" in raw_slots, "fixture should exercise the gapped case"
        assert extract_region_sequence(vl, "CDR3") == "QQHYTTPPT"

    def test_region_with_no_real_residues_yields_empty_not_gaps(self):
        """A region covering zero residues must be absent, not a gap-only string."""
        truncated = number_sequence(
            "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTFTFTISSLQPEDFATYYC"
        )
        assert truncated is not None
        assert extract_region_sequence(truncated, "CDR3") == ""
        assert extract_cdr_l3(truncated) is None

    def test_every_region_sequence_is_alignable(self):
        """Whatever this returns must be scoreable by the substitution matrix."""
        from Bio.Align import substitution_matrices

        alphabet = set(substitution_matrices.load("BLOSUM62").alphabet)
        for seq in (TRASTUZUMAB_VH, TRASTUZUMAB_VL):
            numbered = number_sequence(seq)
            for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"):
                extracted = extract_region_sequence(numbered, region)
                assert not (set(extracted) - alphabet), f"{region}: {extracted!r} not alignable"

    def test_cdr_extraction_wrong_chain_type(self):
        """CDR-H3 extraction should fail on light chains."""
        vl = number_sequence(TRASTUZUMAB_VL)
        assert vl is not None
        cdr_h3 = extract_cdr_h3(vl)  # Wrong: VL is light chain
        assert cdr_h3 is None

    def test_sequence_checksum(self):
        """Sequence checksum should be deterministic."""
        result1 = number_sequence(TRASTUZUMAB_VH)
        result2 = number_sequence(TRASTUZUMAB_VH)
        assert result1.sequence_checksum == result2.sequence_checksum


class TestComponentEligibility:
    """Test comparison eligibility rules (SPEC §20)."""

    def test_exact_same_type(self):
        """Same component types always eligible."""
        assert compare_eligibility(ComponentType.CDR_H3, ComponentType.CDR_H3) \
            == ComponentRelationship.EXACT_SAME_TYPE
        assert compare_eligibility(ComponentType.VH, ComponentType.VH) \
            == ComponentRelationship.EXACT_SAME_TYPE

    def test_vh_vs_vhh_related_class(self):
        """VH vs VHH is related class, not exact match."""
        assert compare_eligibility(ComponentType.VH, ComponentType.VHH) \
            == ComponentRelationship.RELATED_CLASS
        assert compare_eligibility(ComponentType.VHH, ComponentType.VH) \
            == ComponentRelationship.RELATED_CLASS

    def test_kappa_vs_lambda_related_class(self):
        """Kappa vs lambda is related class."""
        assert compare_eligibility(ComponentType.VL_KAPPA, ComponentType.VL_LAMBDA) \
            == ComponentRelationship.RELATED_CLASS
        assert compare_eligibility(ComponentType.VL_LAMBDA, ComponentType.VL_KAPPA) \
            == ComponentRelationship.RELATED_CLASS

    def test_cdr_h3_vs_cdr_l3_excluded(self):
        """CDR-H3 vs CDR-L3 is forbidden (SPEC §20, CLAUDE.md)."""
        assert compare_eligibility(ComponentType.CDR_H3, ComponentType.CDR_L3) \
            == ComponentRelationship.EXCLUDED
        assert compare_eligibility(ComponentType.CDR_L3, ComponentType.CDR_H3) \
            == ComponentRelationship.EXCLUDED

    def test_vh_vs_vl_excluded(self):
        """VH vs VL is forbidden."""
        assert compare_eligibility(ComponentType.VH, ComponentType.VL_KAPPA) \
            == ComponentRelationship.EXCLUDED
        assert compare_eligibility(ComponentType.VL_LAMBDA, ComponentType.VH) \
            == ComponentRelationship.EXCLUDED

    def test_heavy_fr_vs_light_fr_excluded(self):
        """Heavy frameworks vs light frameworks are excluded."""
        assert compare_eligibility(ComponentType.FR_H1, ComponentType.FR_L1) \
            == ComponentRelationship.EXCLUDED
        assert compare_eligibility(ComponentType.FR_L2, ComponentType.FR_H3) \
            == ComponentRelationship.EXCLUDED

    def test_heavy_cdr_vs_light_cdr_excluded(self):
        """Heavy CDRs vs light CDRs are excluded."""
        assert compare_eligibility(ComponentType.CDR_H1, ComponentType.CDR_L1) \
            == ComponentRelationship.EXCLUDED
        assert compare_eligibility(ComponentType.CDR_H2, ComponentType.CDR_L2) \
            == ComponentRelationship.EXCLUDED
