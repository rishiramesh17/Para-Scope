# Gate 2 — Thera-SAbDab Ingestion

Date: 2026-07-29
Source: `TheraSAbDab_SeqStruc_OnlineDownload.csv`
SHA256: `cc9402788297ed7b6ab3d676daa21953a6d28ebbd9cc78ec7248b19f74d49d6d`
Size: 636,357 bytes · 1,133 data rows · source last-modified 2025-03-18
Parser: `src/antigens/parsers/thera_sabdab.py` v1.0.0
Output: `data/cache/combined_corpus_snapshot.pkl`

## What the real file actually contains

The V2 plan assumed a cleaner schema than exists. Verified by direct inspection
of the downloaded file, not from documentation:

- **Variable domains only.** There is no full-length chain column. `FULL_HEAVY` /
  `FULL_LIGHT` therefore remain unpopulated by this source, exactly as with
  SAbDab2. This is a real limit on what the platform can currently claim about
  constant regions, and is recorded rather than worked around.
- **`Format` is free text, not a controlled vocabulary** — 120+ distinct values
  (`Whole mAb`: 800, `Whole mAb ADC`: 58, `Bispecific mAb`: 50, then a long tail
  of bespoke descriptions like `Trispecific Single Domains (VH-VH-VH'-VH''-VH'); Bispecific entry`).
  An exhaustive classifier over this text was **not** attempted; see decisions below.
- **`"na"` is the missing-value sentinel**, a literal string, not a blank cell.
  Treating blanks as the only missing marker would have silently ingested `"na"`
  as a sequence.
- **The clinical-stage column embeds a date**: `Highest_Clin_Trial (Feb '25)`.
  The parser requires this column by exact name so a future release that
  re-stamps it fails loudly rather than silently dropping clinical stage.

## Scope decisions

**Bispecifics are split into independent arms, never cross-paired.**
157 rows carry a genuine second specificity in `HeavySequence(ifbispec)` /
`LightSequence(ifbispec)`. Each row is parsed into up to two arms
(`arm_index` 1 = primary, 2 = bispecific). A bispecific's second heavy chain is
never paired with the primary arm's light chain: the two arms are two different
specificities, not one pair with a spare chain. This follows directly from the
constraint that observed pairing requires both chains from the *same source
record*. Covered by `test_bispecific_arms_are_not_cross_paired`.

**Single-domain (VHH) detection is a documented heuristic, not a claim.**
ANARCII's chain-type output distinguishes H from K/L but not VH from VHH, so
domain architecture has to come from metadata. Detection is a keyword match
against the `Format` text (`nanobody`, `vhh`, `single domain`) — never against
sequence content. Rows it cannot classify default to conventional VH, which
fails safe (a VHH mislabelled VH is compared VH↔VH, which is a same-class
comparison; the reverse would fabricate a missing light chain). The heuristic is
approximate by construction and labelled as such in the module docstring.

## Results

| Metric | SAbDab2 only | Combined | Thera-SAbDab delta |
|---|---|---|---|
| Source rows | 15,641 | 16,774 | +1,133 |
| Heavy chains parsed | 14,157 | 15,445 | +1,288 |
| Light chains parsed | 10,889 | 12,132 | +1,243 |
| Unique components | 39,765 | 46,632 | +6,867 |
| Observed pairs | 10,789 | 12,032 | +1,243 |
| Unique parent antibodies | 5,050 | 6,181 | +1,131 |
| Rows skipped (unsupported) | 1,375 | 1,377 | **+2** |
| Numbering failures | 72 | 72 | **+0** |

Ingestion time: 78.64 s (numbering 1,290 arms of new sequence; SAbDab2 portion
loaded from snapshot).

### Attribution check

The +1,288 heavy chains against 1,290 total arms (1,133 primary + 157
bispecific) accounts for every arm but two — and exactly **2 rows were skipped**,
matching the 2 rows found by direct inspection to have no `HeavySequence` in any
arm. Zero new numbering failures. The 1,375 pre-existing skips and 72 failures
are inherited SAbDab2 baseline, confirmed unchanged.

This attribution matters: it demonstrates the skip/failure counters are
source-attributable rather than an opaque aggregate, so a future ingestion
regression can be localised to the source that caused it.

### Component distribution (combined)

```
VH        5,473    VL-kappa  4,089    VL-lambda 1,100    VHH  1,168
CDR-H1    3,039    CDR-H2    3,734    CDR-H3    5,180
CDR-L1    2,133    CDR-L2      435    CDR-L3    3,271
FR-H1     2,423    FR-H2     2,798    FR-H3     4,305    FR-H4   526
FR-L1     2,056    FR-L2     1,800    FR-L3     2,645    FR-L4   457
```

Cross-source deduplication is active: a sequence present in both SAbDab2 and
Thera-SAbDab becomes one component record carrying both source pointers
(`test_same_sequence_from_both_sources_merges_source_pointers`), so a
therapeutic that has also been crystallised is counted once, not twice.

## Provenance

Both sources are recorded in the combined snapshot manifest with independent
checksums and parser versions:

```
SAbDab2       f0b0aa6dae7d4e06…  parser 1.0.0  15,641 records
Thera-SAbDab  cc9402788297ed7b…  parser 1.0.0   1,133 records
corpus_schema_version 2.0.0 · cache_format_version 1.1.0 · anarcii 2.0.8
```

## Open items

- **License is `unspecified`.** Thera-SAbDab's redistribution terms were not
  confirmed. This is flagged in `build_thera_sabdab_snapshot` and must be
  resolved before any public deployment or data redistribution. It does not
  affect local analysis.
- Clinical stage and development status are stored as **metadata only**. They are
  never treated as a success or failure label, and must not become one in Gate 7
  — a therapeutic reaching Phase 3 transfers no property to a similar sequence.
- `Format` free text is stored verbatim. If a later gate needs finer architecture
  classification (scFv vs. whole mAb vs. ADC), that is a separate, explicitly
  calibrated piece of work, not an extension of the VHH keyword heuristic.
