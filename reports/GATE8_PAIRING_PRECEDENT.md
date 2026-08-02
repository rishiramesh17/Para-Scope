# Gate 8 — Observed Pairing Precedent (SPEC §35.1)

Built 2026-08-01. Corpus: 46,627 components, 5,686 unique pairs,
`corpus_schema_version` 5.0.0. Suite: **303 tests passing** (was 254).

Answers one question: *has this heavy chain been observed paired with this light
chain, in a single source record?* Nothing about whether the submitted chains
will pair, assemble, fold, express, or remain stable together.

---

## The rule the dimension exists to enforce

A pair is only a pair if **one source record links both chains**.

Heavy and light are searched independently — they have to be, since the candidate
is two sequences — but a result is emitted only where the pair index returns a
record for that exact `(heavy_checksum, light_checksum)` key. The index is built
exclusively from same-record co-occurrence, so no combination of independent hits
can manufacture a pair. `_paired_candidates` is the only place a pair can be
created and it is the load-bearing function in the module.

The failure mode this prevents is specific and plausible: take the best heavy hit
and the best light hit, report them together. On a corpus holding
`Alpha(heavy_A, light_A)` and `Beta(heavy_B, light_B)`, querying `heavy_A` with
`light_B` produces two perfect, genuine hits belonging to different antibodies. A
naive implementation reports a confident 1.000/1.000 pair that no experiment ever
observed.

### What the correct answer is *not*

Writing the test for this initially asserted an empty result, and that was wrong.
Alpha's own pair and Beta's own pair each half-match the query — heavy exact,
light approximate, and vice versa — and both are real observed pairs, so both are
legitimately reported. What must never appear is the `(heavy_A, light_B)`
combination.

The test now asserts the fabricated combination is absent, plus the general form:
every emitted pair must be present in the index. That is stronger, because it
fails for any route by which a pair could be manufactured rather than the one
route that was anticipated.

---

## Two counts, never one (SPEC §15)

| | |
|---|---:|
| unique `(heavy, light)` keys | **5,686** |
| backing `PairRecord` entries | **12,032** |
| records per unique pair | 2.12 |

Reporting 12,032 as "pairs" would count one antibody's dozen crystal structures
as a dozen independent pairing observations. `pairs_matched` and
`raw_observation_count` are both always emitted, so neither can be mistaken for
the other, and `independent_paired_record_breadth` counts pairs.

**Corroboration is thinner than the 2.12 average suggests.** The distribution is
heavily skewed:

| backing records | pairs |
|---:|---:|
| 1 | 3,060 |
| 2 | 1,322 |
| 3 | 564 |
| 4 | 380 |
| 5–8 | 262 |
| >8 | 98 |

**54% of observed pairs rest on a single record.** Any future weighting by
observation count is weighting a long tail, not a typical case — worth knowing
before Gate 9 treats corroboration depth as a calibration input.

---

## The two sources are almost disjoint

| source combination | unique pairs |
|---|---:|
| SAbDab2 only | 4,562 |
| Thera-SAbDab only | 1,122 |
| **both** | **2** |

Only two chain pairs in the entire corpus appear in both databases. The
structural and therapeutic sources are not overlapping views of the same
evidence; they describe nearly disjoint sets. Dropping either would lose most of
its contribution rather than a duplicated subset, and a future dedup rule that
assumes cross-source mirroring would find almost nothing to merge here.

Structure support across unique pairs: **4,852 backed, 834 stated absent, 0
unknown.**

---

## No composite is invented

Unlike SPEC §26, **§35 specifies no weights.** It lists eight "provisional
pair-match dimensions" without a formula. So this dimension emits raw
measurements and withholds nothing — there is no composite defined to withhold.
Assembling weights to produce one would be fabricating the spec, and the
resulting number would be indistinguishable downstream from a specified one.

`PairingPrecedentResult` deliberately has **no `score` field at all**, and a test
asserts its absence.

| dimension | state |
|---|---|
| heavy similarity | computable |
| light similarity | computable |
| joint coverage | computable |
| direct parent-record evidence | computable |
| independent paired-record breadth | computable |
| VH–VL germline-pair precedent | **blocked** — no germline assignment exists |
| combined CDR-pattern similarity | **blocked** — needs the §21 length-matched null |
| format relevance | **blocked** — Gate 10.5, the same term that blocks §26 |

Five of eight computable. The three blocked ones name themselves and their reason
on every result, including results where nothing matched.

### Joint coverage is the weaker chain, not the average

Defined as `min(heavy_query_coverage, light_query_coverage)`. A pair cannot be
better covered than its least-covered chain, and averaging would let a fully
covered heavy hide a light chain that barely aligned. **Provisional** — the
choice of `min` over other summaries is not calibrated against anything.

---

## Status routing

| condition | status |
|---|---|
| pairs found | `available` |
| both chains searched, no record links any hit | `not_found` |
| **no light chain supplied** | **`not_applicable`** |
| no heavy chain supplied | `unsupported_component` |
| heavy/light types swapped or invalid | `unsupported_component` |
| corpus holds no pair records | `database_unavailable` |
| alignment raised | `analysis_failed` (per chain; evidence still returned) |

The `not_applicable` case is the one worth stating plainly. **A single-domain
construct has no light chain by design.** Reporting "no pairing evidence found"
for a VHH would describe a deliberate architecture as a missing result, and
downstream that is indistinguishable from an antibody whose pair simply was never
observed. 45 heavy-only arms exist in Thera-SAbDab; every single-domain
therapeutic is among them.

Every non-`available` outcome leaves all measurements unset — never 0.0, never a
midpoint.

---

## Like with like, at pair level

A pair's relationship is the **worse** of its two chains'. A kappa/lambda
substitution downgrades the whole pair to related-class rather than being hidden
behind an exact heavy match — otherwise an exact heavy would carry a
related-class light into same-class evidence, which CLAUDE.md forbids merging.
Related-class pairs rank below same-class at equal identity; the penalty affects
ordering only and alters no reported measurement.

---

## Worked example — trastuzumab VH + VL

```
status                 available
pairs_matched          5
raw_observation_count  16
exact_pair_matches     1
best_joint_coverage    1.0

heavy id              light id              H id   L id  obs  source
VH::8b8998709aa9423e8 VL-kappa::54785c99f  1.000  1.000   11  Thera-SAbDab
VH::5cdeae182bf65f012 VL-kappa::24028a6cb  0.992  0.991    1  Thera-SAbDab
VH::d99c9fe2697371b8a VL-kappa::54785c99f  0.983  1.000    1  Thera-SAbDab
VH::8b8998709aa9423e8 VL-kappa::d15c993f5  1.000  0.982    2  SAbDab2
VH::43d1b31a3fc3f7c29 VL-kappa::54c1223c6  0.868  0.981    1  SAbDab2
```

Five pairs, sixteen backing records. The top pair alone rests on 11 records — the
deduplication rule doing exactly its job, since counting records would have
reported 16 pieces of independent pairing evidence where there are 5.

Note rows 1 and 4: the **same heavy chain** paired with **two different light
chains**, in two different sources. Both are real observations. This is why the
pair, not the chain, is the unit of evidence.

---

## What a pair *is*, not just how similar it is (engine 1.1.0)

*(added 2026-08-02 after design review)*

The measurements above say how closely an observed pair resembles the candidate.
They said nothing about **which antibody it is** — that lived in
`PairRecord.metadata`, an untyped dict, reachable only by knowing the right string
key. Promoted onto each match: therapeutic names, targets, development statuses,
construct formats, PDB entries, experimental methods, resolutions, and a
convenience `display_label`.

Same reasoning as Gate 7's `source_databases`: the information was always present,
but a wrong key returns `None` rather than raising.

### Three measured facts that constrain any display of this

- **80% of pairs have no drug name and no target.** 4,562 of 5,686 are structural
  records — PDB entry, method, resolution, no name. They are not missing an
  identity; they have a different one. A layout with fixed Name/Target slots
  renders empty for four results in five.
- **One chain pair can carry many names and many targets.** Trastuzumab's exact
  pair is backed by 11 records: 11 distinct programmes with **6 different targets**
  and mixed development statuses, because each is a bispecific sharing a
  trastuzumab arm. A single "Target: HER2" field would show one of six. This is
  why every promoted field is a list.
- **Identity is not similarity.** The reported number is `global_identity` — exact
  matches over the longer sequence. Naming it "similarity" would imply conservative
  substitutions are counted; they are not. And identity without coverage is the
  trap `approximate.py` closes, where a random sequence scored identity 1.000 at
  4% coverage. Joint coverage must accompany any identity shown.

`display_label` qualifies with the arm only when more than one exists — saying
"arm 1" on a monospecific implies an arm 2 that does not. Flattening an arm into
the whole antibody is the same conflation that produced the arms-sharing-structures
defect.

## Not implemented, deliberately (SPEC §35.2)

Predicted VH–VL interface compatibility — orientation, interface residues, steric
clashes, buried surface, model confidence — is out of scope for the first
version, and §35.2 requires that such evidence, if it ever exists, stay separate
from observed pairing rather than being converted into it. A test asserts no such
function exists in the module.

---

## Carried forward

- **Joint coverage's `min` is provisional** and needs calibration like every other
  threshold here.
- **54% single-record corroboration** bounds what observation-count weighting can
  achieve; Gate 9 should know this before designing around it.
- **Germline-pair precedent is blocked on the same germline reference set that
  Gate 9 needs** — one dependency, two gates waiting on it.
- **Format relevance is Gate 10.5**, identical to its §26 blocker. One piece of
  work unblocks a dimension in both §26 and §35.

## Reproduce

```bash
.venv/bin/python scripts/build_combined_corpus_cache.py
.venv/bin/python scripts/measure_pairing_precedent.py
.venv/bin/python -m pytest -q
```
