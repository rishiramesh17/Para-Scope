"""Approximate long-sequence index (SPEC S16.2).

Pure-Python local alignment via Bio.Align.PairwiseAligner (BLOSUM62,
gap_open=-10, gap_extend=-1). No external binary (e.g. BLAST) required.

Only long components (VH, VL, VHH, full chains, frameworks) may enter this
index. Short components (CDRs, linkers, junctions) need a length-matched
empirical null model instead (SPEC S21, CLAUDE.md) -- applying unmodified
full-length alignment scoring to them would be a defect, not a shortcut.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from Bio.Align import PairwiseAligner, substitution_matrices

from antigens.indexes.shortlist import (
    DEFAULT_K,
    DEFAULT_SHORTLIST_LIMIT,
    SHORTLIST_UNSUITABLE_TYPES,
    SHORTLIST_VERSION,
    KmerShortlistIndex,
)
from antigens.models import ComponentType

# Component types long enough for whole-sequence alignment scoring to be meaningful.
ELIGIBLE_TYPES = frozenset({
    ComponentType.VH,
    ComponentType.VL_KAPPA,
    ComponentType.VL_LAMBDA,
    ComponentType.VHH,
    ComponentType.FULL_HEAVY,
    ComponentType.FULL_LIGHT,
    ComponentType.FR_H1,
    ComponentType.FR_H2,
    ComponentType.FR_H3,
    ComponentType.FR_H4,
    ComponentType.FR_L1,
    ComponentType.FR_L2,
    ComponentType.FR_L3,
    ComponentType.FR_L4,
})

_BLOSUM62 = substitution_matrices.load("BLOSUM62")

# Letters the loaded substitution matrix can actually score. Derived from the
# matrix rather than hardcoded, so it stays correct if the matrix is ever changed.
_ALIGNABLE_ALPHABET = frozenset(_BLOSUM62.alphabet)

# Version of the scoring + admissibility rules below. Stored on evidence records
# so a score can always be traced to the rule set that produced it (SPEC S19).
#   1.0.0 -- local identity only, no coverage floor (defective: see below)
#   1.1.0 -- added global_identity + coverage admissibility floor
ALIGNMENT_SCORING_VERSION = "1.1.0"

# A local alignment reports identity over the aligned region *only*. On a real
# corpus this is actively misleading: a randomly generated 118-residue query was
# observed to align to a real VH over just 5 residues, all matching, and so
# scored sequence_identity = 1.0 at query_coverage = 0.04. Reported alone, that
# is indistinguishable from a genuine perfect match.
#
# Two independent corrections, because they answer different questions:
#   - global_identity divides identities by the longer of the two *full*
#     sequences, so a partial alignment is penalized in proportion to what it
#     failed to cover. This is the measure downstream scoring should consume.
#   - the coverage floor decides admissibility: an alignment covering a small
#     fraction of either sequence is not evidence about those sequences at all,
#     and is dropped rather than ranked low. Dropping is not the same as scoring
#     zero -- callers receive fewer hits, never a fabricated weak one.
MIN_QUERY_COVERAGE = 0.70
MIN_SUBJECT_COVERAGE = 0.70


def _make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = _BLOSUM62
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    return aligner


@dataclass
class AlignmentHit:
    record_id: str
    alignment_method: str
    sequence_identity: float   # identities / aligned region -- local, NEVER meaningful without coverage
    global_identity: float     # identities / longer full sequence -- coverage-robust, use for ranking/scoring
    query_coverage: float
    subject_coverage: float
    alignment_length: int
    insertions: int
    deletions: int
    scoring_version: str = ALIGNMENT_SCORING_VERSION
    # True when this hit came from a k-mer-shortlisted search, whose hit set is
    # ~96% complete rather than exhaustive. A consumer reporting "no closer match
    # exists" must not do so from shortlisted results -- that would turn a
    # performance heuristic into a biological claim.
    shortlisted: bool = False


def _score_alignment(record_id: str, query: str, subject: str, alignment) -> AlignmentHit:
    counts = alignment.counts()
    alignment_length = int(counts.identities + counts.mismatches + counts.gaps)
    coords = alignment.coordinates
    query_aligned_span = int(coords[0][-1] - coords[0][0])
    subject_aligned_span = int(coords[1][-1] - coords[1][0])

    # Denominator is the longer full sequence, not the aligned span: a hit that
    # matches perfectly over half the domain must not score the same as one that
    # matches perfectly over all of it.
    global_denominator = max(len(query), len(subject))

    return AlignmentHit(
        record_id=record_id,
        alignment_method="pairwise_local_blosum62",
        sequence_identity=float(counts.identities / alignment_length) if alignment_length else 0.0,
        global_identity=float(counts.identities / global_denominator) if global_denominator else 0.0,
        query_coverage=float(query_aligned_span / len(query)) if query else 0.0,
        subject_coverage=float(subject_aligned_span / len(subject)) if subject else 0.0,
        alignment_length=alignment_length,
        insertions=int(counts.insertions),
        deletions=int(counts.deletions),
    )


class ApproximateIndex:
    """Per component type, holds sequences for on-demand pairwise alignment.

    No precomputed all-vs-all index (too large to store speculatively); alignment
    is computed against every same-type sequence at query time.
    """

    def __init__(self, k: int = DEFAULT_K):
        self._sequences: dict[ComponentType, dict[str, str]] = defaultdict(dict)
        # Stage-1 shortlist, kept in lockstep with _sequences by add(). Present
        # unconditionally so a loaded snapshot always supports both search modes.
        self.shortlist_index = KmerShortlistIndex(k=k)

    def add(self, component_type: ComponentType, record_id: str, sequence: str) -> None:
        if component_type not in ELIGIBLE_TYPES:
            raise ValueError(
                f"{component_type} is not eligible for approximate long-sequence alignment; "
                f"short components require a length-matched empirical null model (SPEC S21)"
            )
        # Reject at insertion, not at query. IMGT gap placeholders were reaching
        # framework sequences and BLOSUM62 has no such letter, so every framework
        # query raised "sequence contains letters not in the alphabet" from deep
        # inside Bio.Align -- far from the ingestion that caused it. A sequence
        # this index cannot align must never enter it.
        unalignable = sorted(set(sequence) - _ALIGNABLE_ALPHABET)
        if unalignable:
            raise ValueError(
                f"Cannot index {record_id} ({component_type.value}): sequence contains "
                f"characters the substitution matrix cannot score: {unalignable}. "
                f"IMGT gap placeholders are positional metadata and must be excluded "
                f"from sequence content (see numbering.extract_region_sequence)."
            )
        self._sequences[component_type][record_id] = sequence
        self.shortlist_index.add(component_type, record_id, sequence)

    def search(
        self,
        component_type: ComponentType,
        query: str,
        top_k: int = 10,
        min_query_coverage: float = MIN_QUERY_COVERAGE,
        min_subject_coverage: float = MIN_SUBJECT_COVERAGE,
        shortlist_limit: Optional[int] = None,
    ) -> list[AlignmentHit]:
        """
        Align `query` against same-type indexed sequences and return the best hits.

        Hits failing either coverage floor are omitted entirely: a local alignment
        spanning a small fraction of the sequences says nothing about their
        similarity, and returning it with its (necessarily high) local identity
        would misrepresent a coincidence as a match. Ranking uses global_identity,
        which already accounts for uncovered sequence, so a partial alignment can
        never outrank a fuller one at equal local identity.

        `shortlist_limit` controls two-stage search (Gate 3): a k-mer shortlist
        selects at most that many candidates for full alignment. Scoring is
        identical either way -- the shortlist only decides which alignments run.

        It defaults to None (exhaustive) deliberately. Measured against
        brute-force ground truth on the real corpus, a 200-candidate shortlist
        recovers the single best hit every time (36/36 queries across VH,
        VL-kappa and VHH) but returns only ~96% of the true top-10 on average,
        and as little as 80% on individual queries. Recall does not reach 1.0
        even at limit=500, because k-mer containment and global identity are
        genuinely different rankings: a homolog with substitutions scattered
        every few residues can destroy nearly all shared 5-mers while remaining
        a legitimate match.

        For a system whose output is *comparative evidence*, quietly omitting
        real matches is a correctness problem, not a speed/accuracy tradeoff --
        so completeness is the default and speed is opted into. Hits produced via
        a shortlist are stamped `shortlisted=True` so a caller can never mistake
        an approximate result set for a complete one. Full numbers:
        reports/TWO_STAGE_SEARCH.md.
        """
        subjects = self._sequences.get(component_type, {})

        # Short, highly conserved types gain no speed from shortlisting and lose
        # substantial recall to it, so the request is downgraded to exhaustive.
        if component_type in SHORTLIST_UNSUITABLE_TYPES:
            shortlist_limit = None

        if shortlist_limit is None:
            candidate_ids = list(subjects.keys())
        else:
            candidates = self.shortlist_index.shortlist(component_type, query, limit=shortlist_limit)
            candidate_ids = [c.record_id for c in candidates]

        aligner = _make_aligner()
        hits = []
        for record_id in candidate_ids:
            subject = subjects.get(record_id)
            if subject is None:
                continue
            alignment = aligner.align(query, subject)[0]
            hit = _score_alignment(record_id, query, subject, alignment)
            if hit.query_coverage < min_query_coverage or hit.subject_coverage < min_subject_coverage:
                continue
            hit.shortlisted = shortlist_limit is not None
            hits.append(hit)

        hits.sort(key=lambda h: (h.global_identity, h.query_coverage), reverse=True)
        return hits[:top_k]

    def __len__(self) -> int:
        return sum(len(records) for records in self._sequences.values())
