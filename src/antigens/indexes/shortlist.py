"""
K-mer shortlist index for two-stage approximate search (Version 2, Gate 3).

Stage 1 (this module) cheaply narrows the corpus to plausible candidates; stage 2
(approximate.py) runs full pairwise alignment on that shortlist only.

Motivation, measured on the real corpus (reports/END_TO_END_SEARCH_VALIDATION.md):
a single 123-residue VH query aligned against all 5,473 indexed VH sequences
takes 3.93 s. A whole-construct profile touches VH, VL and eight framework
regions, so a realistic submission costs tens of seconds -- dominated entirely by
alignments against sequences that were never going to be reported.

The shortlist is an OPTIMIZATION, NOT A SCORING CHANGE. Every score a user sees
still comes from the same full alignment as before; this module only decides
which alignments are worth computing. Two consequences follow, and both are
enforced rather than assumed:

  1. It is approximate, so it can in principle miss a true match. Recall against
     brute force is measured on the real corpus and reported, never asserted.
     `shortlist_limit=None` disables stage 1 entirely, which is what the
     recall benchmark and any exactness-critical path use.
  2. A shortlist miss must never be reported as "no similarity found" -- that
     would turn a performance heuristic into a biological claim. Callers that
     surface a not_found status are responsible for saying the search was
     shortlisted; see search_meta on the returned results.

Ranking uses containment (shared k-mers / query k-mers) rather than Jaccard.
Jaccard penalizes a subject for being longer than the query, which is a length
artefact, not evidence of dissimilarity; containment asks the question actually
being posed -- how much of the query is accounted for by this subject.
"""

from dataclasses import dataclass
from typing import Optional

from antigens.models import ComponentType

SHORTLIST_VERSION = "1.0.0"

# k=5 over a 20-letter alphabet: a random 5-mer collision is ~1/3.2M, so shared
# 5-mers between two ~120-residue sequences are overwhelmingly homology rather
# than chance. Chosen against a measured k-sweep, not by assumption -- see the
# recall table in reports/TWO_STAGE_SEARCH.md.
DEFAULT_K = 5

# Default stage-1 breadth when shortlisting is explicitly requested. Measured
# recall against brute force on the real corpus is recorded in
# reports/TWO_STAGE_SEARCH.md -- do not raise or lower this without re-running
# that benchmark.
DEFAULT_SHORTLIST_LIMIT = 200

# Component types where shortlisting is measurably counterproductive, and is
# therefore ignored in favour of exhaustive search even when requested.
#
# Framework regions are short (FR-H3 is ~39 residues) and highly conserved, so
# they yield few k-mers and those k-mers are shared across most of the corpus.
# Containment then ties across thousands of candidates and the limit cuts
# essentially arbitrarily: measured FR-H3 recall@10 was 0.78 mean / 0.40 min at
# limit=200, versus 0.94-0.96 for whole domains. Crucially there is nothing to
# buy -- exhaustive FR-H3 search already runs in 0.32s, because short sequences
# align cheaply. Shortlisting here would trade real evidence for no speed.
#
# Falling back is safe in the only direction that matters: the caller gets a
# more complete result than asked for, and `shortlisted` reports False because
# the search genuinely was exhaustive.
SHORTLIST_UNSUITABLE_TYPES = frozenset({
    ComponentType.FR_H1,
    ComponentType.FR_H2,
    ComponentType.FR_H3,
    ComponentType.FR_H4,
    ComponentType.FR_L1,
    ComponentType.FR_L2,
    ComponentType.FR_L3,
    ComponentType.FR_L4,
})


def extract_kmers(sequence: str, k: int = DEFAULT_K) -> set[str]:
    """Distinct k-mers of `sequence`. Distinct, not multiset: a repeated motif
    should not let one subject dominate the shortlist by repetition alone."""
    if len(sequence) < k:
        return set()
    return {sequence[i : i + k] for i in range(len(sequence) - k + 1)}


@dataclass
class ShortlistCandidate:
    record_id: str
    shared_kmers: int
    containment: float  # shared / query k-mer count


class KmerShortlistIndex:
    """Inverted k-mer index, one posting list per component type.

    Component types never share a posting list. This is the same like-with-like
    constraint the rest of the system enforces: a shortlist that mixed CDR-H3 and
    CDR-L3 candidates would push an ineligible comparison into stage 2, where the
    only thing standing between it and an evidence record is a later check.
    """

    def __init__(self, k: int = DEFAULT_K):
        self.k = k
        # Plain nested dicts, not defaultdict -- these are pickled as part of a
        # corpus snapshot, and a defaultdict with a lambda factory is
        # unpicklable (the bug that blocked Gate 0.5).
        self._postings: dict[ComponentType, dict[str, set[str]]] = {}
        self._record_kmer_counts: dict[str, int] = {}

    def add(self, component_type: ComponentType, record_id: str, sequence: str) -> None:
        kmers = extract_kmers(sequence, self.k)
        if not kmers:
            # Sequence shorter than k. Deliberately not an error: the caller's
            # eligibility rules decide what may be searched, and a too-short
            # sequence simply cannot be shortlisted. It stays reachable via
            # brute force (shortlist_limit=None), so it is skipped here, not lost.
            return
        by_kmer = self._postings.setdefault(component_type, {})
        for kmer in kmers:
            by_kmer.setdefault(kmer, set()).add(record_id)
        self._record_kmer_counts[record_id] = len(kmers)

    def shortlist(
        self,
        component_type: ComponentType,
        query: str,
        limit: int = DEFAULT_SHORTLIST_LIMIT,
    ) -> list[ShortlistCandidate]:
        """Return up to `limit` candidate record_ids, best containment first."""
        query_kmers = extract_kmers(query, self.k)
        if not query_kmers:
            return []

        by_kmer = self._postings.get(component_type)
        if not by_kmer:
            return []

        shared: dict[str, int] = {}
        for kmer in query_kmers:
            for record_id in by_kmer.get(kmer, ()):
                shared[record_id] = shared.get(record_id, 0) + 1

        denominator = len(query_kmers)
        candidates = [
            ShortlistCandidate(
                record_id=record_id,
                shared_kmers=count,
                containment=count / denominator,
            )
            for record_id, count in shared.items()
        ]
        # Tie-break on record_id so the shortlist is deterministic: an unstable
        # shortlist would make two identical submissions return different
        # evidence, which breaks reproducibility.
        candidates.sort(key=lambda c: (-c.containment, c.record_id))
        return candidates[:limit]

    def indexed_record_ids(self, component_type: ComponentType) -> set[str]:
        ids: set[str] = set()
        for record_ids in self._postings.get(component_type, {}).values():
            ids |= record_ids
        return ids

    def __len__(self) -> int:
        return len(self._record_kmer_counts)
