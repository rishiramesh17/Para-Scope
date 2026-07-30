# Gate 3 — Two-Stage Approximate Search

Date: 2026-07-29
Corpus: `data/cache/combined_corpus_snapshot.pkl` — 46,627 components,
28,835 approximate-eligible (`corpus_schema_version` 4.0.0)
Benchmark: `scripts/benchmark_two_stage_search.py` → `reports/two_stage_search_benchmark.json`
Method: brute-force search (`shortlist_limit=None`) is ground truth. 12
deterministically-selected real queries per component type, top_k=10. A "miss"
is a record the exhaustive search reported and the shortlisted search did not.

## Why

Gate 1 measured a single 123-residue VH query against 5,473 indexed VH sequences
at **5.3 s**. A whole-construct profile touches VH, VL and eight framework
regions, so a realistic submission costs tens of seconds — spent almost entirely
on alignments against sequences that were never going to be reported.

Stage 1 is a k-mer inverted index ranked by containment (shared k-mers / query
k-mers). Containment rather than Jaccard: Jaccard penalizes a subject for being
longer than the query, which is a length artefact rather than evidence of
dissimilarity.

**The shortlist is an optimization, not a scoring change.** Every reported score
still comes from the same full pairwise alignment. Stage 1 only decides which
alignments are computed.

## Recall and speed, measured

| Type | Subjects | Brute force | limit | p50 | Speedup | recall@10 mean | recall@10 min | top-1 agree |
|---|---|---|---|---|---|---|---|---|
| VH | 5,473 | 5.301 s | 50 | 0.065 s | 82.0x | 0.850 | 0.600 | 12/12 |
| | | | 100 | 0.113 s | 47.1x | 0.925 | 0.800 | 12/12 |
| | | | **200** | 0.214 s | 24.7x | **0.958** | 0.800 | 12/12 |
| | | | 500 | 0.514 s | 10.3x | 0.967 | 0.900 | 12/12 |
| VL-kappa | 4,089 | 4.270 s | 50 | 0.069 s | 62.1x | 0.917 | 0.800 | 12/12 |
| | | | 100 | 0.122 s | 35.0x | 0.925 | 0.800 | 12/12 |
| | | | **200** | 0.240 s | 17.8x | **0.942** | 0.900 | 12/12 |
| | | | 500 | 0.597 s | 7.2x | 0.942 | 0.900 | 12/12 |
| VHH | 1,168 | 0.620 s | 50 | 0.030 s | 21.0x | 0.875 | 0.700 | 12/12 |
| | | | 100 | 0.055 s | 11.3x | 0.933 | 0.800 | 12/12 |
| | | | **200** | 0.111 s | 5.6x | **0.958** | 0.800 | 12/12 |
| | | | 500 | 0.272 s | 2.3x | 0.967 | 0.800 | 12/12 |
| FR-H3 | 4,304 | **0.324 s** | 50 | 0.008 s | 41.5x | 0.692 | **0.400** | 12/12 |
| | | | 100 | 0.011 s | 29.0x | 0.767 | **0.400** | 12/12 |
| | | | **200** | 0.018 s | 18.0x | **0.783** | **0.400** | 12/12 |
| | | | 500 | 0.039 s | 8.3x | 0.792 | **0.400** | 12/12 |

### k sweep (VH, limit=200, stage-1 recall only)

| k | recall@10 mean | recall@10 min | shortlist p50 | posting entries |
|---|---|---|---|---|
| 4 | 0.9667 | 0.800 | 17.23 ms | 1,062,096 |
| **5** | **0.9750** | 0.800 | 14.94 ms | 1,062,664 |
| 6 | 0.9583 | 0.800 | 12.59 ms | 1,058,198 |

`DEFAULT_K = 5` is chosen on this evidence. k=4 is both slower (denser postings
to walk) and slightly worse; k=6 loses recall as substitutions destroy longer
k-mers more readily.

## Findings that changed the design

### 1. Top-1 is perfectly reliable; the top-10 tail is not

Top-1 agreement was **48/48 across every type and every limit**. The shortlist
never failed to surface the single best match.

Tail recall is a different story: ~96% of the true top-10 on whole domains, and
as low as **0.80 on individual queries**. Recall does not converge to 1.0 even at
limit=500 (VL-kappa is flat at 0.942 from 200 to 500). This is not a breadth
problem that a bigger limit fixes — k-mer containment and global identity are
genuinely different rankings. A homolog carrying substitutions every few residues
can destroy nearly every shared 5-mer while remaining a legitimate match, so it
sorts below unrelated sequences that happen to share conserved motifs.

**Consequence: `shortlist_limit` defaults to `None` (exhaustive).** For a system
whose output is comparative evidence, quietly omitting real matches is a
correctness problem, not a speed/accuracy tradeoff. Speed is opted into.

Hits from a shortlisted search carry `shortlisted=True`. A consumer must not
report "no closer match exists" from such a result set — that would promote a
performance heuristic into a biological claim.

### 2. Shortlisting frameworks is strictly worse than not shortlisting them

FR-H3 recall@10 bottoms out at **0.40** — two in five true hits lost on the worst
query — and never exceeds 0.79 at any limit. Framework regions are short (~39
residues) and highly conserved, so they produce few k-mers and those k-mers are
shared across most of the corpus; containment ties across thousands of candidates
and the limit cuts arbitrarily.

There is also nothing to buy. Exhaustive FR-H3 search costs **0.324 s**, because
short sequences align cheaply. The expensive types (VH 5.3 s, VL-kappa 4.3 s) are
precisely the ones where shortlist recall holds up.

**Consequence: `SHORTLIST_UNSUITABLE_TYPES`** — all eight framework types ignore
a shortlist request and run exhaustively. The fallback errs toward completeness,
and `shortlisted` reports `False` because the search genuinely was exhaustive.

### 3. A crash that exposed a data-integrity bug

The first benchmark run died on FR-H3 with
`ValueError: sequence contains letters not in the alphabet`. Cause: IMGT gap
placeholders (`-`) were stored as region *sequence content*, and BLOSUM62 has no
such letter. 100% of FR-H1 / FR-H3 / FR-L3 / FR-L4 records were affected.

The crash was the least of it. The same characters were being hashed into
component checksums, so identical residue sequences with different IMGT gap
placement were counted as distinct components; correcting it merged 4. Regions
covering zero real residues were also being indexed as components made entirely
of gaps. Full detail in `reports/END_TO_END_SEARCH_VALIDATION.md`; fixed at the
root in `numbering.extract_region_sequence`, with an insertion-time guard in
`ApproximateIndex.add` so an unalignable sequence can never enter the index again.

## Current cost of a full profile

With exhaustive search as the default, one whole-construct submission is
approximately:

```
VH        5.3 s
VL-kappa  4.3 s
8 × FR    ~2.6 s   (0.32 s each, all exhaustive)
--------------
          ~12 s
```

Acceptable for an analysis platform, and correct. If interactive latency later
matters, the honest lever is opting whole domains into shortlisting (24.7x on VH)
while surfacing `shortlisted` in the output — not lowering the default.

## Carried forward

- Recall was measured at top_k=10. If a later gate reports more than 10 matches
  per component, recall must be re-measured at that k — it will be worse.
- 12 queries per type is enough to establish the shape of the problem, not to
  put a tight confidence interval on 0.958. Worth widening if the shortlist is
  ever promoted to a default.
- Stage-1 recall (0.975) exceeds end-to-end recall (0.958) for VH at k=5: the
  shortlist contains true hits that the coverage floor then drops. That gap is
  expected, but it means the two thresholds interact and should be revisited
  together in Gate 9 rather than independently.
