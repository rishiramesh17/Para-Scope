# Gate 10 — Evidence Strength, Confidence, and Completeness

Date: 2026-10-04
Implements SPEC §28 (completeness must be split) and SPEC §19 Layers 3–4
Module: `src/para_scope/evidence_strength.py` (`1.0.0`), thresholds `1.1.0-provisional`
Benchmark: `reports/evidence_grades.json` — 120 real candidates across 3 classes
Tests: 453 passing (was 414)

## The four fields, kept apart

SPEC §28 is blunt: *"Do not use one ambiguous completeness score."* Four
situations get flattened by a single number, and they have different causes,
owners, and remedies:

| Field | Question | If it's low, you… |
|---|---|---|
| **28.1** Input completeness | did the submitter give us the components? | supply more input |
| **28.2** Analysis completion | did our modules run? | retry |
| **28.3** Metadata completeness | are matched source records well annotated? | distrust the annotation |
| **28.4** Evidence sufficiency | is there enough independent evidence? | distrust the conclusion |

A candidate missing a light chain, a snapshot that failed to load, a PDB entry
with no resolution recorded, and a match resting on one record are four different
problems. Averaged into `completeness: 0.62` they are indistinguishable.

So there is **no aggregate completeness field**, and a test asserts none exists —
because the moment one does, every caller uses it and the split becomes
decorative. A test also demonstrates all four moving independently: input
incomplete, analysis complete, metadata 0.0, evidence rich — a combination no
single number can express.

## The distinction that matters most

SPEC §28's closing line: *"A database outage affects analysis completion. It must
not be treated as biological evidence against the candidate."*

A failed lookup and a lookup that ran and found nothing produce the same *absence*
of evidence. Grading them alike turns an outage into a verdict. So:

- Infrastructural failure (`DATABASE_UNAVAILABLE`, `ANALYSIS_FAILED`,
  `ACCESS_RESTRICTED`) → **`Insufficient`**, with a blocked reason saying *retry
  rather than interpret*. Never `D`.
- A module returning `NOT_FOUND` / `NOT_APPLICABLE` / `UNSUPPORTED_COMPONENT`
  **completed successfully**. Counting those as failures would turn an answer into
  an outage.

`D` is a statement about a candidate. `Insufficient` is a statement about us.
Tested in both directions, including the nasty case: strong evidence plus one
outage must not come out worse than weak evidence.

## Layer 4 counts, including the ones we don't have

Reported: independent records, parent antibodies, paired records,
structure-linked records, therapeutic records, dimensions with evidence, and
counts **per evidence class** (never summed across classes — SPEC §17).

Three Layer 4 fields are **not derivable** from the current sources and say so
rather than being omitted or filled with a count that means something else:

- `study_count` — neither source records a study identifier distinct from the
  structure or therapeutic entry.
- `subject_count` — subject-level provenance exists only in repertoire sources;
  `promb`'s bundled reference is pre-aggregated at a 10%-of-subjects threshold.
- `bootstrap_interval` — computable in principle (the calibration artifact stores
  full observations) but not computed, and no interval has been validated.

A reference lookup contributes **zero** independent records. A germline allele is
one curated reference sequence, not database precedent; counting it would let a
single lookup inflate sufficiency.

## Grade benchmark, and what it found

120 candidates (40 each VH / VL-kappa / VHH), fixed seed, exhaustive search:

```
A              24 (20.0%)
B              22 (18.3%)
C              18 (15.0%)
D               0 ( 0.0%)
Insufficient   56 (46.7%)

VH         A=10  B=10  C=10  D=0  Insufficient=10
VL-kappa   A=13  B=11  C=6   D=0  Insufficient=10
VHH        A=1   B=1   C=2   D=0  Insufficient=36
```

### Finding 1 — my first D definition was wrong

The initial implementation defined D purely as low identity. It scored **0 of
120**, because `GRADE_C_MIN_IDENTITY = 0.70` coincides with the admissibility
coverage floor: anything admitted at all had already cleared C.

SPEC §19 actually defines D as *"weak, indirect, predicted, **sparse, or poorly
annotated**"*. Sparseness and annotation are two routes I had simply not
implemented. Added (thresholds `1.0.0` → `1.1.0-provisional`): evidence resting on
a single parent antibody **and** with ≤25% of expected metadata fields populated
grades D regardless of how well it aligns.

### Finding 2 — D is still unreachable, and that is a property of the data

After the fix, still 0 of 120. Measured why:

| Measure | Observed across 120 candidates | D's threshold | Reachable |
|---|---|---|---|
| `metadata_fraction` | **min 0.800**, median 1.000 | ≤ 0.25 | no |
| `best_identity` | **min 0.703**, median 0.927 | < 0.70 | no |
| single-parent candidates | 10, all well annotated | sparse **and** poor | no |

**Grade D is structurally unreachable against SAbDab2 + Thera-SAbDab.** Both
sources are too well curated to produce poorly annotated evidence, and the
coverage floor admits nothing below 0.70 identity.

This is recorded rather than fixed. Lowering thresholds until D appears would be
fitting the scale to the data instead of describing the data. The band stays
defined — it is spec-mandated and would fire on a sparser or less curated source —
and is covered by a test on constructed input so it does not become dead code that
silently stops working.

### Finding 3 — the Insufficient rate is dominated by VHH

46.7% overall looks alarming until it is split: **36 of 40 VHH** candidates have
no therapeutic precedent at all, because only 31 VHH components in the entire
corpus carry a therapeutic link. For VH and VL-kappa the rate is 10 of 40 each.
That is the correct answer — most single-domain antibodies genuinely have no
therapeutic precedent — and it is exactly why `Insufficient` must be
distinguishable from `D`.

## On "benchmarked and versioned"

SPEC §19 requires grade thresholds to be *"benchmarked and versioned"*.

**Versioned:** yes, and *separately* from the engine version
(`GRADE_THRESHOLD_VERSION`), because a threshold change alters every stored grade
without altering any measurement. A test asserts the two versions differ.

**Benchmarked:** no, and the version string says `-provisional` so a stored grade
cannot misrepresent itself. There is no external ground truth of evidence quality
to calibrate against, and inventing one would be worse than admitting its absence.
What *was* verified is the weaker, honest property: that the bands separate real
cases rather than collapsing onto one grade (4 of 5 populated, with the fifth
unreachable for a measured reason). That is a real check — thresholds that put
every candidate in one band carry no information however principled they look.

## Carried forward

- **D unreachable** with current sources (above). Revisit if a less curated source
  is ever ingested.
- **Grade A requires a structure-linked record**, which only Thera-SAbDab entries
  cross-referenced to PDB can supply. A candidate matching a structure-only
  SAbDab2 entry cannot currently reach A through the therapeutic dimension alone.
- **`parent_antibodies` is summed across dimensions**, so a candidate matching the
  same antibody via both therapeutic precedent and pairing counts it twice.
  Cross-dimension deduplication is not implemented; SPEC §15 wants one antibody
  counted once, and this module explicitly delegates that to its callers but no
  caller does it yet. This inflates sufficiency for multi-dimension candidates.
- **No bootstrap interval**, so "confidence in the score" (Layer 4) is reported as
  a grade plus counts rather than as an interval.
- Benchmark is 120 candidates across 3 classes — enough to show the distribution's
  shape, not to pin band boundaries.
