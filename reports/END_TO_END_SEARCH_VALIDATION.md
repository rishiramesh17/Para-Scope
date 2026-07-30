# Gate 1 — End-to-End Search Validation (Real Corpus)

Date: 2026-07-29
Corpus: `data/cache/sabdab2_corpus_snapshot.pkl` (SAbDab2 only, deliberately —
search behaviour is validated *before* the corpus is expanded, so a regression
cannot be blamed on new data)
Benchmark: `benchmarks/search_gold_queries.json` (13 queries, 13 categories)

## Purpose

Phase 1a proved the ingestion pipeline ran. It did not prove that a search
against the real corpus returns defensible results. This gate exercised the
built indexes against real sequences drawn from the corpus itself, plus
adversarial inputs (random, truncated, invalid).

It found one real defect. That is the entire justification for running this
gate before building anything on top of search.

## Corpus under test

| Metric | Value |
|---|---|
| Source rows | 15,641 |
| Heavy chains parsed | 14,157 |
| Light chains parsed | 10,889 |
| Unique components (deduplicated) | 39,765 |
| Observed pairs | 10,789 |
| Unique parent antibodies | 5,050 |
| Rows skipped (unsupported) | 1,375 |
| Numbering failures | 72 |

Reproducibility: rebuilding the snapshot from the numbering cache produced
byte-identical statistics (39,765 components / 10,789 pairs / 5,050 parents),
confirming ingestion is deterministic given a fixed source checksum and
ANARCII version.

## Ingestion cost, measured

The Phase 1a estimate of ~113 minutes was wrong by ~13x.

| Run | Time |
|---|---|
| First ingest (ANARCII numbering, cold) | 511.14 s |
| Snapshot save | 0.87 s |
| Re-ingest with numbering cache hit | **6.85 s** |

The numbering cache delivers a ~75x speedup on re-ingest and is the reason
every later gate can iterate on indexing logic without re-numbering. Note the
cache is keyed on (source checksum, ANARCII version) *independently* of the
corpus snapshot, so a Corpus code change invalidates the snapshot but not the
expensive numbering work.

## Query latency, measured (combined corpus, 28,840 approximate-eligible sequences)

| Operation | Latency |
|---|---|
| Exact index lookup (checksum) | p50 0.0003 ms |
| Approximate VH search (5,473 subjects, 123-residue query) | **p50 3.93 s** |

The exact path is effectively free. The approximate path is a full pairwise
alignment against every same-type sequence, and 3.93 s for a *single* component
is the binding constraint: a whole-construct profile touches VH, VL and eight
framework regions, so a realistic submission costs tens of seconds. This is the
motivation for Gate 3 (two-stage shortlist + alignment), and the number to beat.

## Defect found: fragmentary local alignment reported as a perfect match

**Severity: high.** Would have produced scientifically indefensible evidence
records, silently.

A randomly generated 118-residue sequence was submitted as a VH query. It
returned a top hit with:

```
sequence_identity = 1.0
query_coverage    = 0.042
subject_coverage  = 0.031
alignment_length  = 5
```

`sequence_identity` was computed as `identities / alignment_length` — that is,
over the *locally aligned region only*. A local aligner will always find some
short high-scoring segment between any two protein sequences, so a coincidental
5-residue match scored 1.0. Consumed downstream (Gates 9–11), that value is
indistinguishable from a genuine full-length identical match.

This inverts the constraint in CLAUDE.md — *"short components require empirical
chance correction; never apply unmodified full-length BLAST scoring to them"* —
by letting a 5-residue fragment stand in for a full-length score.

### Fix (`ALIGNMENT_SCORING_VERSION` 1.0.0 → 1.1.0)

Two independent corrections, because they answer different questions:

1. **`global_identity`** added: `identities / max(len(query), len(subject))`.
   A partial alignment is penalized in proportion to what it failed to cover.
   Ranking now uses this, so a partial alignment can never outrank a fuller one
   at equal local identity. `sequence_identity` is retained unchanged — it is a
   legitimate measurement of the aligned region — but is now explicitly
   documented as never meaningful without coverage.
2. **Coverage admissibility floor** (`MIN_QUERY_COVERAGE` /
   `MIN_SUBJECT_COVERAGE`, both 0.70, configurable per call): hits below the
   floor are *omitted*, not ranked low. A coincidental fragment is not weak
   evidence of similarity, it is no evidence — consistent with "missing evidence
   is not zero evidence."

Regression tests added at both levels (unit, in `tests/test_search.py`; real
corpus, in `tests/integration/test_real_corpus_search.py`), including an
assertion that a 10-residue fragment returns *zero* VH hits rather than a
perfect-looking one.

## Second defect found: stale corpus snapshot loaded silently

Loading the real SAbDab2 snapshot after Gate 2 added
`_therapeutic_name_to_arm_count` to `Corpus.__init__` raised an `AttributeError`
deep inside `ingest_thera_sabdab`. Cause: `pickle` restores an object's stored
`__dict__` verbatim, so a snapshot predating a new attribute yields an object
silently missing it, and the failure surfaces far from its origin.

Fix: `CORPUS_SCHEMA_VERSION` is stamped onto the instance, recorded in the
manifest, and checked on load; a mismatch raises `StaleSnapshotError` naming the
remedy. `CACHE_FORMAT_VERSION` is checked the same way. Both are covered by
tests. This matters increasingly as Gates 4–11 each add state to `Corpus`.

## Results

13 gold queries across 13 categories, all passing after the fixes:

| Category | Result |
|---|---|
| exact_match | PASS — exact checksum hit ranks first |
| single_substitution | PASS |
| several_substitution | PASS |
| similar_nonidentical | PASS — distinct records not conflated |
| long_cdr_h3 | PASS — indexed and retrievable by exact hash |
| insertion_coded_cdr | PASS — insertion codes preserved |
| heavy_only | PASS — succeeds without a light chain |
| light_only | PASS — succeeds without a heavy chain |
| paired_vh_vl | PASS — directly observed pair found |
| invalid_sequence | PASS — rejected, not scored |
| truncated_sequence | PASS — not silently misrepresented |
| duplicate_parent_record | PASS — multi-structure parent does not inflate counts |
| no_close_match | PASS *(after fix; this is the query that exposed the defect)* |

Full suite: **81 passed** (Gate 0 baseline: 36).

## Carried forward

- Approximate search at 3.93 s/component is the Gate 3 target.
- The 0.70 coverage floor is a defensible default, not a calibrated value. Gate 9
  (empirical reference distributions) should revisit it per component type —
  short frameworks (FR-H4, ~11 residues) and full domains plausibly warrant
  different floors.
- 1,375 skipped rows and 72 numbering failures are inherited Phase 1a baseline,
  unchanged by this gate and not yet individually reviewed.
