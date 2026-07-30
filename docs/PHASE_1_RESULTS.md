# Phase 1a Results: SAbDab2 Ingestion

Real run against the downloaded SAbDab2 ML-grade dataset (`data/snapshots/sabdab2/`,
Zenodo record 20083995, `ab_split.csv`, MD5 `0dbb4cc499e9eb77f14008b232f2c38c`).
Numbers below are from `scripts/run_sabdab2_ingest.py`, not estimates.

## Corpus stats

| Metric | Value |
|---|---|
| Rows in source file | 15,641 |
| Rows skipped (unsupported) | 1,375 |
| Heavy chains numbered | 14,157 |
| Light chains numbered | 10,889 |
| Numbering failures | 72 |
| **Numbering success rate** | **99.71%** (25,046 / 25,118 attempted chains) |
| Unique indexed components | 39,765 |
| Directly observed pairs (row-level) | 10,789 |
| Unique heavy+light pair combinations | 4,564 |
| Unique parent antibodies (SABDAB_ID) | 5,050 |

Acceptance target from `PHASE_1_DESIGN.md` was ">99% numbering success" and
"~5K unique sequences (order-of-magnitude)" for parent-antibody dedup -- both met.

SAbDab2 itself reports 5,301 unique SABDAB_IDs in this release; we recover 5,050
because 251 parent antibodies are represented *only* by rows this parser
intentionally skips (VNAR-only or fused-single-chain-only entries -- see below).

### Why rows are skipped

| Reason | Count |
|---|---|
| VNAR (shark IgNAR, no ComponentType mapping) | 67 |
| Engineered single-chain construct (scFv/diabody/SCDB) | 1,308 |

**Engineered single-chain constructs** were the majority finding of this run. SAbDab2
reports the *same* fused VH-linker-VL polypeptide in both `Hseq_expected` and
`Lseq_expected` for scFv, diabody, and single-chain-diabody rows -- these formats
don't have independently annotated heavy/light domains in the source data. The
parser (`parsers/sabdab2.py`) now detects `Hseq_expected == Lseq_expected` and marks
the row `unsupported_component` rather than numbering the fusion as if it were an
isolated VH or VL.

This was caught by benchmarking, not by inspection: an initial run (before this fix)
saw a 95.0% numbering success rate against a >99% target. Investigating the shortfall
showed ~1,320 failures clustered almost exactly against the count of `SCFV`/`DIABODY`/
`SCDB` rows (1,232 + 56 + 6 = 1,294, vs. 1,320 actual failures). Worse than the
failures themselves: a handful of these fused sequences *did* number successfully as
a plausible-looking VL or VH -- which would have silently indexed a wrong component
(a fusion mislabeled as an isolated chain) rather than failing loudly. Fixed by
skipping the row entirely once `Hseq_expected == Lseq_expected` is detected.

### Numbering failures (the remaining 72, post-fix)

| Category | Count |
|---|---|
| Invalid input (non-standard residue, e.g. `X`) | 59 |
| Chain-type mismatch (ANARCII returned an unexpected chain type) | 13 |
| ANARCII internal error | 0 |

Both categories are handled per CLAUDE.md's evidence-handling rules: rejected
outright, never guessed at or substituted.

## Component type breakdown

| Component type | Count |
|---|---|
| VH | 4,363 |
| VL-kappa | 3,180 |
| VL-lambda | 954 |
| VHH | 1,138 |
| CDR-H1 / CDR-H2 / CDR-H3 | 2,713 / 3,237 / 4,367 |
| CDR-L1 / CDR-L2 / CDR-L3 | 1,824 / 380 / 2,747 |
| FR-H1..FR-H4 | 2,235 / 2,421 / 3,691 / 511 |
| FR-L1..FR-L4 | 1,821 / 1,544 / 2,208 / 431 |

CDR-L2 is the smallest category (380) -- expected, since CDR-L2 is short (IMGT
56-65) and highly conserved across kappa/lambda light chains, so many collapse into
the same exact-sequence bucket.

`FULL_HEAVY` / `FULL_LIGHT` (constant-region-inclusive full chains) are **not**
populated by this corpus: SAbDab2's ML-grade release is cropped to the variable
region (IMGT 1-128) only. Full-chain evidence will show `not_applicable` /
`insufficient_evidence` against this source until a corpus with constant-region
sequences (e.g. Thera-SAbDab, Phase 1b) is ingested -- this is correct
spec-compliant behavior (CLAUDE.md: missing evidence is not zero evidence), not
a gap to paper over.

## Search latency benchmark

Single-query timings against the full ingested corpus:

| Index | Latency |
|---|---|
| Exact-hash lookup | <1 ms |
| Positional lookup (>95% overlap) | 56 ms |
| Approximate alignment (VH, top-10 of 4,363 candidates) | **6.0 s** |

**This contradicts the Phase 1 design assumption** ("<1s per query" for a ~10K-sequence
corpus). At 4,363 same-type candidates, exhaustive pairwise local alignment
(`Bio.Align.PairwiseAligner`, BLOSUM62) takes ~1.4ms per alignment, which adds up to
6 seconds for one query. This is acceptable for offline/batch evidence generation but
too slow for an interactive query path. Flagged as a Phase 2 optimization target:
a cheap pre-filter (k-mer/minhash candidate shortlist, or length-bucketing) ahead of
full alignment would cut the candidate pool before the expensive step runs, the same
way exact-sequence deduplication already happens before IMGT numbering.

## Ingestion performance

Batched ANARCII numbering (`numbering.number_sequences`, batch size 500) over ~25K
unique chain sequences took ~6,810s (~113 minutes) end-to-end for this run size and
hardware (CPU-only, Apple Silicon). This is dominated entirely by IMGT numbering, not
parsing or indexing (snapshot + parse + index steps are all sub-second by comparison).

## Known limitations (carried from PHASE_1_DESIGN.md, still accurate)

- Paired-OAS ingestion deferred to Phase 2 (storage decision pending real query patterns).
- Short-CDR empirical null models (SPEC S21) not built yet -- CDR components are
  indexed (exact + positional) but have no calibrated scoring.
- Evidence grades, calibrated percentiles, and paratope residue-set comparison all
  deferred (SPEC S19 Layers 2-4, SPEC S34).
- Study/subject-level repertoire dedup deferred (applies once OAS is added).
- Thera-SAbDab (Phase 1b) not yet ingested -- pair index and therapeutic-precedent
  evidence class are SAbDab2-only for now.

## Reproducing this run

```
.venv/bin/python scripts/run_sabdab2_ingest.py
```

Requires `data/snapshots/sabdab2/splits_final/ab_split.csv` (downloaded from
Zenodo record 20083995; not checked into git). Full report written to
`data/snapshots/sabdab2/ingest_report.json`.
