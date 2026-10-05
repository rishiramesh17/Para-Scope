# Gate 11 — Whole-Construct Profile Assembly (SPEC §37)

Date: 2026-10-05
Modules: `src/para_scope/profile.py` (`1.0.0`), `src/para_scope/short_component.py` (`1.0.0`)
Weight set: `1.0.0-provisional-igg-like` · Calibration: 617 distributions
End-to-end run: `reports/whole_construct_profile.json`
Tests: 581 passing (was 500)

## The index emits, and withholds, on real data

Four candidates derived from one real whole-mAb therapeutic, run through every
gate — search, therapeutic precedent, calibration, grading, format relevance,
assembly:

| Case | Index | Coverage | Components scored | Why |
|---|---|---|---|---|
| self | **80.8** | 95% | 16/16 | emits; only Constant/Fc absent |
| mutated ×10 | withheld | 81% | 13/16 | CDR-H3 lost its precedent |
| mutated ×30 | withheld | 59% | 9/15 | CDR-H3 + both CDR groups lost |
| novel (random) | withheld | **0%** | 0/0 | nothing matched at all |

Coverage falls monotonically as the candidate becomes more novel, and the index is
withheld — with named reasons — the moment a *required* component loses evidence.
A novel candidate does not get a low index; it gets no index and an explanation.

The self-match scores **80.8, not 100**, because structure-linked precedent is 0.0
for this therapeutic and the geometric mean refuses to let strong routine
components hide it. That is §37's stated reason for choosing a geometric mean, and
`bottleneck_slot` names the culprit.

## The property this gate exists to protect

§37's index divides by the weight that *contributed*:

```
100 × exp( Σ[wᵢ × rᵢ × ln(max(sᵢ/100, ε))] / Σ[wᵢ × rᵢ] )
```

Drop a component and numerator and denominator shrink together. Measured with
uniform component scores of 90:

| Components present | Index | Coverage |
|---|---|---|
| all 10 slots | **90.0** | 100% |
| 9 (no constant) | **90.0** | 95% |
| 3 (required only) | **90.0** | 42% |

**The index is identical at 100% and 42% coverage.** Nothing in the formula
notices. `evidence_coverage` — computed against the *full* weight set, never
against what contributed — is the only thing that distinguishes them, and
`test_index_is_identical_at_full_and_partial_coverage` pins it. If that test ever
fails because coverage stopped being reported, every index the platform emits
begins overstating itself.

## Scope is declared; the name is derived

`coverage_scope` drives the displayed name rather than the reverse:

- no constant-region evidence → `VARIABLE_DOMAIN_ONLY` →
  *"Comparative Evidence Index (variable domain only)"*
- constant-region evidence present → `FULL_CONSTRUCT` →
  *"Whole-Construct Comparative Evidence Index"*

Verified by test in both directions. The concept is never renamed — only
qualified — so when constant-region coverage improves the qualifier disappears on
its own, with no second name for one index across reports and git history.

"Critical" is per format *and* per scope: `REQUIRED_COMPONENTS` puts Constant/Fc in
the full-construct requirement set and not the variable-domain one. A VHH gets its
own weight set with no VL and no pairing slot at all, and a test asserts it is not
penalised for lacking a light chain it never had. An scFv, bispecific or fusion
format gets **no index** — §37 names reusing its weights as the mistake not to
make, and an index under the wrong weight set is wrong in a way no coverage field
reveals.

## What blocked the index, and the sub-gate it forced

The first end-to-end run scored **2 of 16 components**. Only VH and VL. Every CDR
and framework was unscored, so CDR-H3 — a required component — was missing and the
index was withheld for *every* case including the self-match.

The chain, traced rather than guessed:

1. CDRs are excluded from the alignment index, correctly — SPEC §21 requires a
   length-matched empirical null before any short-component similarity is
   meaningful.
2. So no similarity distribution exists for CDR strata.
3. So §26's composite has no `calibrated_component_similarity`.
4. So the CDR has no score.
5. So §37 withholds the index.

Gate 9 had explicitly recorded §21's null as **not done** and not claimed. It turns
out to be the load-bearing gap: without it the platform's headline number cannot
be produced at all. So §21 was built here.

### SPEC §21, implemented

`short_component.py` computes the spec's empirical p-value verbatim —
`(1 + nulls at least as strong) / (1 + total nulls)`, where the `+1` on both sides
bounds p away from zero, because a finite null cannot establish impossibility — plus
Benjamini-Hochberg q-values across the searched set, position-aware similarity, and
length compatibility reported **separately** as §21 requires.

The null is **pairs of distinct real components** of the same type and length band,
not random strings. Real CDRs share compositional biases that random sequence does
not; a null without that bias would make every real match look spectacular. 18
null strata built (`reports/` → calibration artifact), and the result justifies the
stratification immediately:

| Stratum | median null identity | p95 |
|---|---|---|
| CDR-H3 len 15–17 | **0.133** | 0.333 |
| CDR-H3 len 22–26 | 0.136 | 0.273 |
| **CDR-L1 len 12–14** | **0.583** | 0.833 |

Two unrelated CDR-H3s share 13% identity; two unrelated CDR-L1s share 58%.
Light-chain CDR1 is far more conserved, so using one null for the other would call
ordinary CDR-L1 similarity spectacular. 9 strata were skipped for too few
comparable pairs and say so.

With the null in place a real CDR-H3 exact match reaches a calibrated similarity of
**0.9992** (p = 1/1201 — it exceeds the entire background), and all 16 components
score.

## Three further defects found while assembling

### 1. `blocked_terms` contradicted the score

I cleared a term's block when *any* match supplied it, while the composite read the
best match. That produced `blocked_terms == []` beside `format_relevance == None`
and a withheld score — directly contradictory output. Blocks are now derived from
the composite's actual inputs.

### 2. One match's gaps withheld the whole composite

Taking all six terms from the single best-by-identity match was brittle: that
match's uncomparable format or absent structure flag withheld the composite even
when other matches were fully characterised. The composite now uses the best match
**among those supplying all six terms**, and records which one in
`scoring_match_id` — so the score remains attributable to a real record rather
than becoming an aggregate of several. Mixing terms across matches was rejected:
it would describe no match that exists.

### 3. Sub-region components could not see their own arm's facts

`direct_pair_or_structure_support` was `None` for every CDR, because arms name the
*variable-domain* records they contributed and no arm references a CDR record. A
single-armed therapeutic now lends its arm-level facts to its regions.

A multi-armed therapeutic deliberately does **not**. Its arms are two
specificities with different structures, and nothing records which arm a given CDR
came from — attributing arm 1's structure to a region that may have come from arm 2
is exactly the cross-arm error ingestion already refuses to make. Those regions
keep `None`, which reads as "not established" rather than "absent".

## SPEC §26's composite now emits

All six weighted terms exist for the first time (Gate 9 supplied calibrated
similarity, Gate 10.5 format relevance, §21 the short-component case). The
composite emits **only when all six are present**: dropping a term and
renormalising over the rest would produce a number that looks complete while
answering a narrower question — the same failure CLAUDE.md prohibits for the §37
index, and no more acceptable here. `SCORE_FORMULA_VERSION` is `1.0.0-provisional`.

## Carried forward

- **Constant/Fc never contributes**, so no real candidate reaches
  `FULL_CONSTRUCT` scope today. Constant-region components exist (Gate 9.1
  recovered 1,643) but are not wired into a component score, so every real index
  renders with the variable-domain qualifier. That is honest and visible, and it is
  the single biggest remaining coverage gap.
- **No index for scFv, bispecific, or fusion formats** — no versioned weight set
  exists. That is 10 of 104 format values but includes every bispecific.
- **The single-domain weight set is my own**, derived by renormalising §37's
  heavy-side weights. §37 requires a separate versioned set for VHH and this is
  one, but it is provisional in exactly the way §37's own table is and rests on no
  calibration.
- **`RELIABILITY_BY_GRADE` is a judgement.** §37 weights each component by a
  reliability without defining one; mapping A/B/C/D to 1.0/0.8/0.6/0.3 is
  defensible and uncalibrated.
- **Lunsekimig1/2/3 canonicalisation is still open.** 29 therapeutics are numbered
  variants of one name, so `therapeutics_matched` can overcount by up to 2. Not
  auto-merged: a numeric-suffix rule could merge genuinely distinct drugs.
- **`MAX_LENGTH_DIFFERENCE_FRACTION = 0.25`** and
  `MINIMUM_NULL_COMPARISONS = 200` join the provisional-threshold list alongside
  the coverage floor, grade thresholds, and `MINIMUM_WINDOWS_FOR_FRACTION`.
