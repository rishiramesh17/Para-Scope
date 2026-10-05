# Gate 9 — Empirical reference distributions (SPEC §19 Layer 2)

Date: 2026-10-04
Artifact: `data/cache/calibration_v1.json` — 599 distributions, 8.5 MB, `calibration_version` 1.0.0
Corpus: `combined_corpus_snapshot.pkl`, schema 6.0.0, 46,727 components
Tests: 350 passing (was 316)

## What this gate changes

Before it, every gate emitted Layer 1 raw measurements and none could say whether
a number was ordinary. "This VH is 1.000 identical to its nearest therapeutic VH"
sounds extreme and *is* extreme — but "0.84" sounds unremarkable and is in fact
close to the population median of 0.829, while "0.92" is 95th percentile. None of
that was knowable, so Gate 7 and Gate 8 both carried blocked calibrated terms.

Layer 2 places a measurement in a versioned empirical distribution. Layer 1 is
untouched: calibration adds a layer, and a test asserts raw identities and match
ordering are byte-identical with and without it.

## The population, stated plainly

The reference population is **SAbDab2 structural entries plus Thera-SAbDab
therapeutic entries** — antibodies someone chose to crystallise or to develop as
a drug. It is **not** a sample of natural antibody repertoire and nothing makes
it one.

Every distribution carries that sentence, and a test asserts none of the 599 omits
it. "91st percentile of solved structures" is defensible; "91st percentile of
human antibodies" is not, and the artifact is built so a reader cannot slide from
one to the other. The population's own quantiles ship alongside every percentile
so the number can be sanity-checked rather than taken on trust.

A percentile is also not a grade. An unusual value is unusual — not worse,
riskier, or less likely to succeed. The disclaimer says so in those words, and a
test asserts it.

## What was built

| Family | Distributions | Source |
|---|---|---|
| `feature.*` | 287 | Gate 4 — length, pI, GRAVY, aromaticity, net charge, cysteine count, glycosylation sequons |
| `motif.*` | 287 | Gate 5 — per-rule counts |
| `window_coverage.coverage` | 10 | Gate 6 |
| `germline.nearest_identity` | 3 | Gate 9.0 |
| `similarity.best_global_identity_vs_therapeutics` | 12 | this gate |

### Germline identity — a biological sanity check that passes

| Class | n | p05 | median | p95 | max |
|---|---|---|---|---|---|
| VH | 698 | 0.508 | **0.680** | 0.800 | 0.845 |
| VL-kappa | 700 | 0.571 | **0.759** | 0.879 | 0.896 |
| VL-lambda | 700 | 0.550 | **0.811** | 0.883 | 0.931 |

High but nowhere near 1.0, which is what somatic hypermutation should look like.
A population pinned at 1.0 would mean assignment was returning self-matches; near
0.3 would mean it was aligning against the wrong locus. An integration test
asserts the median stays in 0.55–0.95 and the maximum stays below 1.0, so a future
regression in germline assignment fails here rather than quietly shifting every
novelty percentile.

### Similarity to nearest therapeutic — VH

| min | p05 | p25 | median | p75 | p95 | max |
|---|---|---|---|---|---|---|
| 0.630 | 0.711 | 0.792 | **0.829** | 0.866 | 0.915 | 0.992 |

Self-matches are excluded: a sampled component that *is* a therapeutic would
otherwise score 1.0 against itself, and a population piled at 1.0 would make
every genuine near-match look unremarkable. A test asserts fewer than 10% of
values reach 0.999.

This is the distribution that unblocks Gate 7. Trastuzumab's VH, which is
byte-identical to 11 therapeutic programmes, lands at the **100th percentile** —
correct, and now demonstrable rather than asserted.

## Design decisions

**Built offline, loaded read-only.** Scoring never computes a distribution
(CLAUDE.md: scoring code never queries live data). A corpus update produces a new
calibration version rather than silently changing a percentile a stored analysis
already reported. The artifact is JSON, not pickle — a calibration is a published
claim and should be inspectable and diffable.

**Full observations stored, not summary quantiles.** 8.5 MB buys exact
reproducibility and lets a later version compute a different statistic from the
same population without rebuilding it.

**Length-banded strata for CDRs.** SPEC §19 requires "CDR-H3 compared with CDR-H3
of similar length". A 22-residue CDR-H3 has no business being percentiled against
a population dominated by 12-residue ones; length drives both the measurement and
its chance distribution. Whole variable domains vary little and are not banded.
Requesting a length-stratified class without a length raises rather than silently
pooling across lengths.

**Ties sit at the centre of the tied block.** Saturated measurements — framework
window coverage, zero-valued motif counts — are mostly ties. Taking either edge
would report a confidently extreme percentile for the single most ordinary value
in the population.

**Insufficient data is a status, never a wider bin.** Below
`MINIMUM_SAMPLE_SIZE = 30` the result is `INSUFFICIENT_EVIDENCE` with the sample
size attached. It does not widen to a neighbouring band, fall back to the
unstratified distribution, or return the 50th percentile — each would answer a
question nobody asked while looking like an answer to the one that was. **128 of
599 strata (21.4%) are below the minimum** and are deliberately unreportable.

**A missing measurement is `NOT_APPLICABLE`, not a percentile of zero.**

## What this does and does not unblock

**Unblocked: Gate 7's `calibrated_component_similarity`.** Verified end-to-end on
the real corpus — the term is listed as blocked without a calibration and cleared
with one, and only where a percentile actually came back. Supplying a calibration
that holds no distribution for a stratum leaves the term blocked, with the
percentile's own reason explaining why.

**Still blocked, and why:**

- **`format_relevance`** (Gates 7 and 8) — needs the candidate's declared
  construct format, which the submission schema does not collect, plus a
  controlled vocabulary over a 120-value free-text column. Gate 10.5, not this
  gate. 0.10 of SPEC §26's weight stays missing.
- **`vh_vl_germline_pair_precedent`** (Gate 8) — the *original* blocker is gone:
  germline assignment now exists and returns a call for every variable domain
  tested. What remains is the reference table the dimension needs — which VH
  germline gene is observed paired with which VL gene, and how often, across the
  12,032 observed pairs. That means germline-assigning every unique paired
  variable domain (~9,000 components × one alignment per reference allele), a
  separate offline job. Until it exists the dimension reports nothing rather than
  inferring precedent from the candidate's own two gene calls, which would
  describe co-occurrence in one candidate rather than precedent in a database.
- **`combined_cdr_pattern_similarity`** (Gate 8) — SPEC §21 wants a length-matched
  **empirical null**: what identity would two *unrelated* CDRs of this length
  share by chance? This gate built distributions of *observed* values, which is a
  different object. The length-banded machinery is the foundation, but a null
  model needs a null population and has not been built. Not claimed as done.

## Carried forward

- `MINIMUM_SAMPLE_SIZE = 30` is a judgement, not a power calculation. Chosen so a
  single outlier moves a percentile by at most ~3.3 points.
- Germline distributions rest on a 700-per-class sample and similarity on 120 —
  both fixed-seed and recorded per distribution, but a wider sample would tighten
  the tails, which is where percentiles are least stable and most quoted.
- Window coverage yields only 10 reportable strata: most CDR and framework classes
  produce fewer than `MINIMUM_WINDOWS_FOR_FRACTION = 20` nine-mers, so coverage is
  non-`AVAILABLE` for them. That interacts with the `MINIMUM_WINDOWS` threshold,
  which is itself provisional — the two should be revisited together.
- The structure-derived population bias is **not fixed** and is not fixable with
  these sources. `promb`'s human-OAS reference is a genuine repertoire baseline
  for window coverage specifically; nothing else here has one.
