# Gate 10.5 — Construct Format Vocabulary and Format Relevance

Date: 2026-10-05
Module: `src/para_scope/construct_format.py`, vocabulary `1.2.0`
Measurement: `reports/format_vocabulary.json` (1,133 rows, 104 distinct Format values)
Tests: 500 passing (was 453)

## Why this gate existed

SPEC §26 gives format relevance **0.10** of the therapeutic-precedent composite,
and SPEC §35 needs the same term. Nothing produced it. It needs the candidate's
*declared* construct format — which the submission schema never collected — plus a
controlled vocabulary over a free-text column. Neither falls out of any other
gate, so left unowned every gate would have completed and the composite would have
shipped computed over 0.90 of its own definition with nothing flagging the gap.

**Both terms are now unblocked**, verified end-to-end on the real corpus:
`blocked_terms` for §26 goes from `['format_relevance']` to `[]` when a format is
declared, and §35's `blocked_dimensions` loses `format_relevance` likewise.

## The reframe that made it tractable

The earlier assessment was that a classifier over 104 free-text values was
infeasible. That is true of a **flat taxonomy** and false of the data.
`"Bispecific Whole mAb ADC"` is three independent facts; a flat scheme needs a
bucket per combination and drowns. Parsed as four independent axes:

| Axis | Rows resolved | Distinct values resolved |
|---|---|---|
| specificity | **1133 / 1133 (100.0%)** | 104 / 104 |
| architecture | **1112 / 1133 (98.1%)** | 84 / 104 |
| modification | **1133 / 1133 (100.0%)** | 104 / 104 |
| species | n/a — `UNSTATED` is a real answer | — |

Row coverage is the honest figure, not distinct-value coverage: 70.6% of this
column is the single string `"Whole mAb"`, so counting distinct values understates
reach badly. Both are reported.

The 21 unresolved rows are bespoke domain-topology strings
(`"Bispecific Mixed Domains (L-kappa-H-gamma1_[VH-G1(CH1-h)]_L-kappa)"`,
`"Dimerised heavy domain (VH-h-CH2-CH3-CHS dimer)"`). They return `UNKNOWN`, and
`UNKNOWN` never becomes agreement or disagreement.

## Relevance behaves like a measure of comparability

Scored against every value in the column, by declared format:

| Declared | Rows scored | Mean relevance |
|---|---|---|
| Whole mAb | 1052 / 1133 (92.9%) | **0.881** |
| Bispecific Whole mAb | 1052 / 1133 | 0.651 |
| scFv | 1052 / 1133 | 0.390 |
| Nanobody | 1052 / 1133 | 0.374 |

This is the sanity check that matters: declaring a whole mAb yields high mean
relevance because most of the corpus *is* whole mAb, and declaring an scFv yields
low mean relevance because little of it is. A metric insensitive to the
declaration would have produced four similar numbers.

## Three defects found during development, all mine

### 1. Modifier agreement masquerading as a format match

Two constructs whose architecture was `UNKNOWN` **on both sides** scored
`relevance = 1.0`, because they happened to share `"bispecific"` and an *inferred*
absence of conjugation. A perfect-format-match reading, with the axis that decides
what kind of molecule each is never assessed.

Fixed: architecture is a `REQUIRED_AXIS`. If it is uncomparable, relevance is
`None` with `INSUFFICIENT_EVIDENCE`. Specificity, modification and species *modify*
a construct; architecture *is* the construct.

### 2. The weights contradicted that argument

With architecture at 0.40 against the modifiers' 0.60, `"Whole mAb"` vs `"scFv"`
scored **0.556** while `"Whole mAb"` vs `"Bispecific Whole mAb ADC"` scored
**0.444** — disagreeing on the single most important axis scored *better* than
disagreeing on two lesser ones.

Fixed (vocabulary `1.0.0` → `1.1.0`): architecture 0.55, exceeding the other three
combined. An scFv carries no constant region and links its domains differently,
which changes what a VH match means more than adding a second specificity does.
Relevance now decreases monotonically with added format differences — asserted by
test, since this was wrong twice.

### 3. A partial unknown scored as a definite mismatch

`"Bispecific mAb"` (60 rows) scored **0.130** against declared `"Whole mAb"`,
treating `MAB_UNSPECIFIED` as an architecture *mismatch*. But that text says the
construct is an antibody and stops — it does **not** say the construct is not a
whole mAb. Scoring it as a mismatch answers a question the source never answered,
for a pair that may be architecturally identical.

Fixed (`1.1.0` → `1.2.0`): `MAB_UNSPECIFIED` is unresolved *for comparison* while
remaining a distinct enum value, because "an antibody of unstated sub-architecture"
is worth recording — it just cannot support agreement or disagreement. This is
what accounts for the 7.1% of rows that score no relevance.

## Design rules

- **Uncomparable is never 0.0.** Zero asserts the formats are maximally unrelated.
  Too few comparable axes → `None` + `INSUFFICIENT_EVIDENCE`, so §26's weight is
  reportable as missing rather than silently satisfied.
- **A missing declared format is the normal case**, not an error — the schema field
  is only now arriving. Absent it, §26's 0.10 stays explicitly unmet, which is the
  state the roadmap insisted must never be hidden by renormalising over 0.90.
- **Inference is recorded as inference.** Absence of "bispecific" reads as
  monospecific — a convention of *this column*, not a universal rule — and the
  axis is listed in `inferred_axes`. Absence of a species word is **not** read as
  human; that would be a claim, so it stays `UNSTATED`.
- **Renormalising over comparable axes is sound here and not in §26**, because
  `comparable_axes` and `comparable_weight` travel with the number. A reader sees a
  1.0 resting on 0.92 of the weight for what it is. §26's composite has no
  equivalent field, which is exactly why dropping a term there is prohibited.

## Source-data defect found

One row carries **`"Lenzilumab"`** — a therapeutic name — in the `Format` column.
It parses to `MAB_UNSPECIFIED`/`UNKNOWN` and is correctly uncomparable, so it does
no harm, but it is a curation error in Thera-SAbDab rather than an unusual format.
Not worked around: silently special-casing a drug name would hide the next one.

## Carried forward

- **An undeclared species leaves cross-species matches unpenalised.**
  `"Whole mAb"` vs `"Canine Whole mAb"` scores 1.000, because the species axis is
  uncomparable when the candidate does not declare one. The remedy is for the
  submission to declare species; the alternative — defaulting a candidate to human
  — would be a claim.
- **Axis weights are provisional**, like the grade thresholds and the coverage
  floor. They are internally consistent (architecture dominates, ordering is
  monotone) but rest on no external ground truth.
- **§35 takes the best relevance across a pair's records.** A pair backed by 37
  records of several formats scores on whichever is most comparable, since that is
  the record the evidence rests on. `FormatRelevance.matched` names which format
  produced the number, so the choice is auditable — but a consumer wanting the
  spread across records would need to compute it.
- **`DeclaredFormat` is not yet a submission-schema field.** The parameter exists
  on both `compute_therapeutic_precedent` and `compute_pairing_precedent`, and
  accepts a parsed `ConstructFormat`, but nothing collects it from a user yet.
  Gate 11 is where a submission actually arrives.
