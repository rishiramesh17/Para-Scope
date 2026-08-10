"""
IMGT numbering wrapper around ANARCII (SPEC §9).

Captures tool version, numbering confidence, and original coordinates.
All scoring code uses canonical IMGT coordinates internally.
"""

import hashlib
from dataclasses import dataclass
from typing import Optional
from anarcii import Anarcii
import anarcii as anarcii_module


def get_anarcii_version() -> str:
    """Report ANARCII package version."""
    return getattr(anarcii_module, '__version__', 'unknown')


@dataclass
class NumberedComponent:
    """One component after IMGT numbering."""
    original_sequence: str  # Preserved verbatim
    sequence_checksum: str  # sha256 hex
    imgt_numbering: list[tuple[tuple[int, str], str]]  # [((position, insertion_code), residue), ...]
    chain_type: str  # 'H' or 'L'
    scheme: str  # 'imgt' or alternate
    numbering_score: float  # Confidence metric from ANARCII
    numbering_error: Optional[str]
    tool_name: str
    tool_version: str


def _validate_sequence(sequence: str) -> Optional[str]:
    """Normalize and validate; returns the cleaned sequence, or None if invalid."""
    sequence = sequence.strip().upper()
    if not sequence or not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in sequence):
        return None
    return sequence


def compute_sequence_checksum(sequence: str) -> str:
    """SHA256 of a normalized (stripped, uppercased) sequence -- the corpus dedup key.

    Callers that already hold a normalized sequence (e.g. corpus.py resolving a
    chain against a numbering cache) must use this same normalization, or cache
    lookups will silently miss.
    """
    return hashlib.sha256(sequence.strip().upper().encode()).hexdigest()


def _record_to_component(sequence: str, checksum: str, record: dict) -> Optional[NumberedComponent]:
    if record.get('error'):
        return None
    return NumberedComponent(
        original_sequence=sequence,
        sequence_checksum=checksum,
        imgt_numbering=record['numbering'],
        chain_type=record['chain_type'],
        scheme=record['scheme'],
        numbering_score=float(record['score']),
        numbering_error=record.get('error'),
        tool_name='ANARCII',
        tool_version=get_anarcii_version(),
    )


def number_sequence(sequence: str, cpu: bool = True) -> Optional[NumberedComponent]:
    """
    Number an antibody sequence under IMGT using ANARCII (SPEC §9).

    Args:
        sequence: Amino-acid sequence, uppercase.
        cpu: Use CPU-only mode (set False if GPU available).

    Returns:
        NumberedComponent if numbering succeeded, None if it failed.
    """
    validated = _validate_sequence(sequence)
    if validated is None:
        return None

    checksum = compute_sequence_checksum(validated)

    numberer = Anarcii(
        seq_type='antibody',
        mode='accuracy',
        cpu=cpu,
        ncpu=1,
        verbose=False
    )

    result = numberer.number([validated])
    if not result:
        return None

    record = list(result.values())[0]
    return _record_to_component(validated, checksum, record)


def number_sequences(
    sequences: list[str], cpu: bool = True, batch_size: int = 500
) -> dict[str, Optional[NumberedComponent]]:
    """
    Batch-number many sequences under IMGT, deduplicating by checksum first (SPEC S9, S15).

    A single batched ANARCII call over N unique sequences is substantially faster
    than N repeated single-sequence calls, so corpus ingestion should always
    deduplicate and batch through this function rather than looping number_sequence().

    Returns a dict keyed by sequence checksum -> NumberedComponent, or None for
    checksums whose sequence failed to number. Invalid input sequences (non-residue
    characters, empty) are silently excluded from the result, matching
    number_sequence()'s treatment of the same input.
    """
    unique_by_checksum: dict[str, str] = {}
    for sequence in sequences:
        validated = _validate_sequence(sequence)
        if validated is None:
            continue
        checksum = compute_sequence_checksum(validated)
        unique_by_checksum.setdefault(checksum, validated)

    results: dict[str, Optional[NumberedComponent]] = {}
    checksums = list(unique_by_checksum.keys())
    seqs = list(unique_by_checksum.values())

    for start in range(0, len(seqs), batch_size):
        batch_checksums = checksums[start:start + batch_size]
        batch_seqs = seqs[start:start + batch_size]

        numberer = Anarcii(
            seq_type='antibody',
            mode='accuracy',
            cpu=cpu,
            ncpu=1,
            verbose=False
        )
        raw_results = numberer.number(batch_seqs)

        for checksum, seq, record in zip(batch_checksums, batch_seqs, raw_results.values()):
            results[checksum] = _record_to_component(seq, checksum, record)

    return results


# IMGT region boundaries (numeric position, inclusive) for a V-domain (SPEC S9).
# Insertion-coded positions (e.g. 111A, 111B) belong to whichever numeric position
# they extend, so region membership is decided by numeric position alone -- a CDR-H3
# using insertion codes for extra length must not be truncated by dropping those residues.
IMGT_REGION_BOUNDARIES: dict[str, tuple[int, int]] = {
    'FR1': (1, 26),
    'CDR1': (27, 38),
    'FR2': (39, 55),
    'CDR2': (56, 65),
    'FR3': (66, 104),
    'CDR3': (105, 117),
    'FR4': (118, 128),
}

HEAVY_CHAIN_TYPES = ('H',)
LIGHT_CHAIN_TYPES = ('K', 'L')  # K = kappa, L = lambda

# IMGT is a fixed-width scheme: every position slot exists whether or not a given
# chain has a residue there, and ANARCII emits '.'/'-' for the empty slots. Those
# placeholders are positional metadata, not sequence content, and must never reach
# a component's sequence -- they would be hashed into its checksum (splitting one
# real component into several), counted in its length, fed to the feature engine,
# and passed to an aligner whose substitution matrix has no such letter.
# They remain in NumberedComponent.imgt_numbering, which is where position
# information belongs.
IMGT_GAP_CHARACTERS = frozenset({'-', '.'})


def extract_region_sequence(numbered: NumberedComponent, region: str) -> str:
    """
    Residues of one IMGT region as actual sequence content, gap placeholders removed.

    Returns '' when the region has no residues (component too short to cover it).
    Performs no chain-type filtering -- callers that need it use extract_region.
    """
    if region not in IMGT_REGION_BOUNDARIES:
        raise ValueError(f"Unknown IMGT region: {region}")

    start, end = IMGT_REGION_BOUNDARIES[region]
    return ''.join(
        residue
        for (pos, _code), residue in numbered.imgt_numbering
        if start <= pos <= end and residue not in IMGT_GAP_CHARACTERS
    )


def extract_region(numbered: NumberedComponent, region: str, allowed_chain_types: tuple[str, ...]) -> Optional[str]:
    """
    Extract a named IMGT region (FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4) from a numbered component.

    Returns None if the component's chain type isn't in allowed_chain_types, or if no
    residues fall in the region (e.g. component too short to cover it). IMGT gap
    placeholders are excluded from the returned sequence.
    """
    if numbered.chain_type not in allowed_chain_types:
        return None

    residues = extract_region_sequence(numbered, region)
    return residues if residues else None


def extract_cdr_h3(numbered: NumberedComponent) -> Optional[str]:
    """Extract CDR-H3 (IMGT 105-117) from a numbered heavy-chain component."""
    return extract_region(numbered, 'CDR3', HEAVY_CHAIN_TYPES)


def extract_cdr_l3(numbered: NumberedComponent) -> Optional[str]:
    """Extract CDR-L3 (IMGT 105-117) from a numbered light-chain component."""
    return extract_region(numbered, 'CDR3', LIGHT_CHAIN_TYPES)
