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


def extract_variable_domain(numbered: NumberedComponent) -> str:
    """
    The numbered variable domain alone, gap placeholders removed.

    This is the component a VH/VL/VHH comparison is actually about. ANARCII
    numbers only the V-domain, so for a chain submitted with additional sequence
    (a Fab heavy chain carries CH1; some entries carry tags or fusion partners)
    `numbered.original_sequence` is strictly longer than what was numbered.
    """
    return "".join(
        residue
        for (_pos, _code), residue in numbered.imgt_numbering
        if residue not in IMGT_GAP_CHARACTERS
    )


def split_chain_at_variable_domain(
    numbered: NumberedComponent,
) -> tuple[str, str, str]:
    """
    Split a source chain into (prefix, variable_domain, suffix).

    Why this exists
    ---------------
    The source chain was being stored *whole* as the VH / VL-kappa / VL-lambda
    component. Measured on the real corpus, that made 68% of `VH` components
    Fab heavy chains rather than V-domains: median stored length 222 residues
    against a median numbered V-domain of 121, the remainder being CH1.

    That silently violated the project's comparison-integrity rule that a
    variable-domain match is never a full-chain match. Constant domains are
    near-identical between unrelated antibodies, so including them handed every
    comparison a free similarity boost -- measured at +5.6 points of identity to
    the nearest *other* antibody (0.858 whole-chain vs 0.803 V-domain-only).
    The bias is directional: it makes sequences look less novel than they are.

    Returns three strings, any of which may be empty:
      prefix  residues N-terminal to the V-domain. Unknown identity -- a tag, a
              signal peptide, a fusion partner, or another domain entirely. NOT
              a constant region, and deliberately not typed as one.
      domain  the numbered V-domain (what VH/VL comparison must use).
      suffix  residues C-terminal to the V-domain. For a Fab chain this is the
              constant domain, but it is not decomposed into CH1 / Hinge / CH2 /
              CH3 / CL here -- that needs a constant-region reference
              (SPEC S6.5) which is not bundled. It is stored undecomposed.

    The V-domain is located by substring search rather than by offset
    arithmetic: ANARCII reports positions, not input offsets, and insertion
    codes mean position count does not equal residue count. Verified contiguous
    in 11,830 / 11,830 variable-domain components of the real corpus.
    """
    domain = extract_variable_domain(numbered)
    full = numbered.original_sequence
    if not domain:
        return "", "", ""

    start = full.find(domain)
    if start < 0:
        # The numbered residues are not a contiguous substring of the input.
        # Never silently fall back to the whole chain -- that is the defect this
        # function exists to remove. Report the domain with no flanks instead,
        # so a caller indexes the V-domain and loses no correctness.
        return "", domain, ""

    return full[:start], domain, full[start + len(domain) :]


# Component types whose corpus records hold a numbered variable domain, and whose
# queries must therefore be reduced to one before comparison.
_VARIABLE_DOMAIN_TYPE_VALUES = frozenset({"VH", "VL-kappa", "VL-lambda", "VHH"})


def normalize_query_to_variable_domain(
    component_type_value: str, query: str
) -> tuple[str, Optional[str]]:
    """
    Reduce a submitted chain to its variable domain, for comparison against the corpus.

    Corpus variable-domain components hold the numbered V-domain only. A query
    submitted as a whole chain (a Fab heavy chain, or a construct carrying a tag
    or fusion partner) must be reduced the same way, or the comparison is
    domain-against-chain: exact lookup misses, and coverage is capped below the
    admissibility floor no matter how good the match is.

    Returns `(sequence_to_compare, normalization_note)`. The note is None when
    nothing was removed, and otherwise states exactly what was removed and from
    which end -- the caller must surface it, because the input is not what was
    compared. The original string is never mutated; the caller keeps it.

    Non-variable-domain types (CDRs, frameworks, constant regions) are returned
    unchanged: they are stored as extracted, not as chains.
    """
    cleaned = query.strip().upper()
    if component_type_value not in _VARIABLE_DOMAIN_TYPE_VALUES or not cleaned:
        return cleaned, None

    numbered = number_sequence(cleaned)
    if numbered is None or numbered.numbering_error is not None:
        # Unnumberable. Return the input untouched rather than guessing at a
        # domain boundary; the caller's own status handling decides what that
        # means. Silently trimming here would be the fabrication this avoids.
        return cleaned, None

    prefix, domain, suffix = split_chain_at_variable_domain(numbered)
    if not domain or (not prefix and not suffix):
        return cleaned, None

    parts = []
    if prefix:
        parts.append(f"{len(prefix)} residue(s) N-terminal to the variable domain")
    if suffix:
        parts.append(f"{len(suffix)} residue(s) C-terminal to the variable domain")
    note = (
        f"Compared the numbered variable domain ({len(domain)} residues) rather "
        f"than the full {len(cleaned)}-residue input: removed " + " and ".join(parts) +
        ". Corpus variable-domain components hold the V-domain only, so comparing "
        "the whole chain would measure constant-region and flanking similarity too."
    )
    return domain, note



def extract_cdr_h3(numbered: NumberedComponent) -> Optional[str]:
    """Extract CDR-H3 (IMGT 105-117) from a numbered heavy-chain component."""
    return extract_region(numbered, 'CDR3', HEAVY_CHAIN_TYPES)


def extract_cdr_l3(numbered: NumberedComponent) -> Optional[str]:
    """Extract CDR-L3 (IMGT 105-117) from a numbered light-chain component."""
    return extract_region(numbered, 'CDR3', LIGHT_CHAIN_TYPES)
