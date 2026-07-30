# Version 2 — Gate 0: Baseline Audit

Snapshot of repository state immediately before Version 2 work begins. This is
the reference point every later gate's "no regression" claims are measured
against.

## Repository structure

```
src/antigens/
  __init__.py
  cli.py                    -- CLI entry point (antigens command)
  models.py                 -- ComponentType, EvidenceClass, EvidenceRecord, compare_eligibility()
  numbering.py               -- ANARCII wrapper, IMGT region extraction, batch numbering
  corpus.py                  -- Ingestion pipeline: parse -> dedup -> number -> index
  indexes/
    __init__.py
    exact.py                 -- ExactIndex (checksum -> record_ids)
    approximate.py            -- ApproximateIndex (pairwise alignment, long sequences only)
    positional.py             -- PositionalIndex (IMGT position Jaccard overlap)
    pair.py                   -- PairIndex (directly observed heavy/light pairs)
  parsers/
    __init__.py
    sabdab2.py                -- SAbDab2 ab_split.csv parser

tests/
  test_numbering.py           -- 13 tests: numbering, CDR extraction, eligibility rules
  test_search.py              -- 23 tests: all 4 indexes, parser, corpus ingestion

scripts/
  run_sabdab2_ingest.py       -- one-off ingestion benchmark script (not in pytest suite)

docs/
  SPEC.md                    -- 51-section master specification
  PHASE_1_DESIGN.md           -- Phase 1a design doc
  PHASE_1_RESULTS.md          -- Phase 1a real ingestion results
```

No API layer, no CLI beyond single-sequence numbering, no persistence layer.
Everything downstream of `Corpus()` lives in memory only.

## Test suite baseline

```
36 passed in 6.23s
```

All Phase 1a tests pass with no failures, no skips, no warnings. Full listing
in `reports/v2_baseline_metrics.json`.

## Corpus baseline (SAbDab2 only)

| Metric | Value |
|---|---|
| Source rows | 15,641 |
| Rows skipped (unsupported) | 1,375 |
| Heavy chains numbered | 14,157 |
| Light chains numbered | 10,889 |
| Numbering success rate | 99.71% |
| Unique searchable components | 39,765 |
| Directly observed pairs | 10,789 |
| Unique parent antibodies | 5,050 |
| Full ingestion time | 6,810.86s (~113 min), ANARCII-bound |

Reproduced from `data/snapshots/sabdab2/ingest_report.json` (real run, not an
estimate). Full detail in `docs/PHASE_1_RESULTS.md`.

## Query latency baseline

| Index | Latency |
|---|---|
| Exact-hash lookup | <1 ms |
| Positional lookup (>95% overlap) | 56 ms |
| Approximate alignment (VH, top-10 of 4,363) | 6.0 s |

This is the number Gate 3 must improve without sacrificing recall.

## Software versions

| Package | Version |
|---|---|
| Python | 3.13.3 |
| anarcii | 2.0.8 |
| biopython | 1.87 |
| pytest | 9.1.1 |

## Data hygiene confirmation

- `data/` (1.1 GB on disk) is **not** tracked by git — confirmed via
  `git ls-files | grep data/` returning zero results.
- `.gitignore` excludes `data/`, `.venv/`, caches, and build artifacts.
- Every `SourceSnapshot` record carries `raw_checksum` (SHA256 of the source
  file) and `parser_version`, so the corpus is reproducible from the checksum
  plus the recorded Zenodo record ID — the raw file itself is never committed.

## Reproducibility gap identified at this gate

**The corpus cannot currently be reproduced without a ~113-minute wait.**
`Corpus.ingest_sabdab2()` builds everything in memory; there is no
serialization step. This blocks every gate from Gate 1 onward, each of which
needs to load an already-ingested corpus repeatedly (integration tests,
benchmarks, regression comparisons). Addressed immediately as Gate 0.5, before
Gate 1 begins.

## Known gaps entering Version 2 (carried forward, not new)

- Thera-SAbDab not yet ingested (Gate 2).
- `FULL_HEAVY` / `FULL_LIGHT` component types have zero records — SAbDab2's
  ML-grade release is variable-domain-only. This does **not** get resolved by
  Thera-SAbDab either (see plan revision on Thera-SAbDab scope) — it remains
  `insufficient_evidence` through all of Version 2 core, which is correct
  behavior per CLAUDE.md, not a defect.
- Approximate search is exhaustive; no shortlist stage exists yet (Gate 3).
- No calibration population exists for SPEC §19 Layers 2-4 — Gate 9 builds
  the first one.
- `SourceSnapshot.license_id` for SAbDab2 is recorded as `"unspecified"` —
  flagged for confirmation before any redistribution, not blocking for
  internal analysis use.

## Gate 0 verdict

No unexplained test failures. No data-hygiene violations. No blockers to
proceeding. Gate 0.5 (persistence layer) begins next.
