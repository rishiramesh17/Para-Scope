# Version 2 Roadmap — Gate Status, Decisions, and Blockers

**This file is the source of truth for project state.** `CLAUDE.md` points here and
instructs any session to read it before planning work. If it disagrees with anyone's
recollection — human or model — this file is right.

Last verified: 2026-08-29, against commit on `main`, 316 tests passing.

Work proceeds in numbered gates. A gate is a unit of work that ends in a verifiable
artifact, not a time box. Gates 0–8 are complete.

---

## Status

| Gate | Scope | Status | Artifact |
|---|---|---|---|
| 0 | Baseline audit | ✅ | `reports/V2_BASELINE_AUDIT.md` |
| 0.5 | Persistence layer | ✅ | `src/para_scope/persistence.py` |
| 1 | End-to-end search validation | ✅ | `reports/END_TO_END_SEARCH_VALIDATION.md` |
| 2 | Thera-SAbDab ingestion | ✅ | `reports/THERASABDAB_INGESTION.md` |
| 3 | Two-stage approximate search | ✅ | `reports/TWO_STAGE_SEARCH.md` |
| 4 | Deterministic feature engine | ✅ | `src/para_scope/features.py` |
| 5 | Motif observation engine | ✅ | `src/para_scope/motifs.py` |
| 6 | Human repertoire window coverage | ✅ | `reports/GATE6_WINDOW_COVERAGE.md` |
| 7 | Therapeutic-precedent profile | ✅ | `reports/GATE7_THERAPEUTIC_PRECEDENT.md` |
| 8 | Pairing precedent | ✅ | `reports/GATE8_PAIRING_PRECEDENT.md` |
| 9.0 | Germline reference data (milestone 1) | ✅ | `reports/GATE9_GERMLINE.md` |
| **9** | **Empirical reference distributions** | **ready — blocks calibrated terms in Gates 7 and 8** | SPEC §19 Layers 2–4 |
| 10 | Evidence strength / confidence / completeness | — | SPEC §28 |
| **10.5** | **Construct format vocabulary + declared format input** | **— must land before 11; see below** | SPEC §26 format relevance |
| 11 | Final profile assembly | — | SPEC §37 |

Corpus under test: 46,627 components, 28,835 approximate-eligible,
`corpus_schema_version` 5.0.0, `CACHE_FORMAT_VERSION` 1.2.0.

---

## Gate 6 — human repertoire window coverage ✅

Built 2026-07-30. Full detail in `reports/GATE6_WINDOW_COVERAGE.md`; what follows is
what a future session needs to not re-derive.

### BioPhi was rejected, and is not coming back

BioPhi requires Python 3.9, conda, and a 2020-vintage ANARCII. The project runs
Python 3.13 with `anarcii>=2.0.8`. This is an unresolvable stack collision, not an
integration difficulty. Do not re-litigate it.

**But note what shipped:** `promb` is authored by BioPhi's author (David Prihoda) and
reproduces BioPhi's OASis score. The same published metric (Prihoda et al., *mAbs*
2022), repackaged for a modern stack. MIT licensed, reference data bundled in the
wheel, zero network calls at scoring time. So the metric was never the problem — the
distribution was.

### The name: `human_window_coverage`, never "naturalness"

*(decided 2026-07-30; supersedes this gate's former title)*

The roadmap previously called this "naturalness" and flagged the llama result as the
reason that word is unsafe. The corpus-wide run makes the case far harder than four
sequences did:

| Component | n | p50 coverage |
|---|---:|---:|
| VH | 5,473 | 0.429 |
| **VHH (camelid)** | 1,168 | **0.547** |
| VL-kappa | 4,089 | 0.445 |

**Camelid VHH have a higher median coverage than VH.** Not merely overlapping —
higher, with a tighter distribution. Any name implying humanness would rank llama
single domains above the median variable domain.

The term is **Human Repertoire Window Coverage**: it states the comparison performed,
not a conclusion drawn from it, matching `Observed Repertoire Precedent` and
`Database-Relative Novelty`. The reference set's id and version are fields on the
result, so swapping reference sets later does not force a rename — the same principle
as the Gate 11 index scope decision. A test asserts no public name in the module
contains "humanness", "naturalness", or "immunogenic".

### What the metric cannot do — do not rediscover these

- **CDR-H3: uninformative.** All 28 scoreable CDR-H3s scored exactly 0.000. CDR-H3 is
  junctionally generated, so its 9-mers appear in no set built from peptides shared
  across ≥10% of subjects. Expected biology, not a defect — but a low CDR-H3 coverage
  is not a finding, and Gate 9 must exclude it rather than build a percentile on it.
- **Frameworks: saturated.** FR-L3's p50, p75 and p95 are all 1.000; FR-H3's p75 and
  p95 are 1.000. A saturated distribution cannot discriminate.
- **Only 40.3% of the corpus scores at all** — 41.2% `insufficient_evidence`, 18.4%
  `not_applicable`. Most CDRs are shorter than the 9-residue window. This is why
  status reporting is mandatory: the majority of components produce no value and
  none of them may be recorded as a zero.

### Two library behaviours deliberately not inherited

`promb` raises `ZeroDivisionError` for any sequence shorter than the window (i.e. for
most CDRs), and **silently scores sequences containing non-standard residues** — an
`X` returns 0.47 rather than an error, because windows spanning it match nothing and
deflate the fraction. The second is the dangerous one: it fails silently in the
direction that produces a confident wrong answer. Handled as `NOT_APPLICABLE` and
`UNSUPPORTED_COMPONENT` respectively, input never modified.

### Supporting change

`EvidenceStatus` now exists in `models.py` — the seven `CLAUDE.md` statuses plus
`AVAILABLE`, per SPEC §22. Before this gate they were ad-hoc string prefixes in parser
skip reasons. Every dimension added from here on should use it.

---

## Germline reference data ✅ — resolved 2026-08-29

*(was "unresolved, and larger than it looks" — resolved as Gate 9 milestone 1,
full detail in `reports/GATE9_GERMLINE.md`)*

**Source and license, verified rather than assumed:** OGRDB
(https://ogrdb.airr-community.org) publishes its germline sets under **CC0**
(public domain) — confirmed directly against
`https://ogrdb.airr-community.org/render_page/licensing_statement.html`. This
is a real difference from SAbDab2/Thera-SAbDab's `unspecified` license: the
reference set is bundled inside the package
(`src/para_scope/reference_data/germline/`, 353 V-gene alleles) rather than
gitignored and fetched separately, the same way `promb` bundles its own data.

**Option B (build on the existing alignment engine) was correct.** No new
alignment code — germline V genes are indexed into the same `ApproximateIndex`
(SPEC §16.2) Gates 7 and 8 already use, typed by locus so comparison stays
like-with-like automatically. `src/para_scope/germline.py`
(`GERMLINE_VERSION = "1.0.0"`) mirrors `therapeutic_precedent.py` /
`pairing_precedent.py`'s shape exactly.

**Scope, deliberately reduced** (matches the option originally recommended
here): human only (no camelid reference set — VHH is `UNSUPPORTED_COMPONENT`,
not scored low), V gene only (no D/J call), one reference set (SPEC §11.2's
never-silently-merge-two-germline-sets concern does not arise with only one —
revisit if IMGT/GENE-DB is ever added alongside it).

**New carried-forward item found while measuring:** ~70% of sampled VH/VL
components return `NOT_FOUND` rather than a value, traced to a real length
mismatch — a share of `ComponentType.VH`-typed components are ~210–240
residues (roughly double a typical single VH domain) against a ~98-residue
germline reference, so `query_coverage` caps near 0.44 regardless of match
quality and never clears the 0.70 floor inherited unchanged from
`indexes/approximate.py`. Not fixed here — that floor was calibrated for
whole-domain-vs-whole-domain matching, and loosening it without calibration
would repeat the exact mistake the existing "0.70 coverage floor ... not a
calibrated value" note below already warns against. Gate 9 proper should give
germline admissibility its own, separately calibrated floor.

This unblocks Database-Relative Novelty's baseline and Gate 8's blocked
`vh_vl_germline_pair_precedent` dimension (SPEC §35) — both were waiting on
this, not on each other.

---

## The cross-cutting constraint: `FULL_HEAVY` / `FULL_LIGHT` are permanently empty

**Neither SAbDab2 nor Thera-SAbDab carries constant-region sequence.** SAbDab2's
ML-grade release is variable-domain-only; Thera-SAbDab curates variable domains plus
metadata, where isotype is a *label* (`IgG1`, `IgG4`), not sequence. This is confirmed
and will not change within Version 2.

Filling these by appending germline constant regions would be fabricating sequence the
source never provided — prohibited by `CLAUDE.md` "never silently modify input."
They correctly remain `insufficient_evidence` throughout. This is correct behavior,
not a defect.

### Consequences for Gate 11 — decide before assembly, not during

1. **Make the index state its own scope — do not rename it.** *(decided 2026-07-30;
   supersedes an earlier proposal to rename it to "Variable-Domain Comparative
   Evidence Index")*

   Renaming would bake a temporary data gap into permanent vocabulary. Variable-
   domain-only is a V2 corpus limitation, not a property of the product — SPEC §6.6,
   §32, and §33 all put construct-level elements in scope. The moment constant-region
   data lands, the name would have to change back, leaving two names for one index
   across reports, formula versions, and git history.

   The overclaim does not come from the name. It comes from nothing emitting coverage.
   §37 *already* requires `Evidence coverage` and `Missing critical evidence` as
   outputs; implement those and the overclaim disappears with more precision than a
   rename would give, since a rename never says *which* slots were empty.

   Required fields on every emitted index — `coverage_scope`
   (`VARIABLE_DOMAIN_ONLY` | `FULL_CONSTRUCT`), `components_present`,
   `components_missing` (each with status and reason), `weight_set_version`. The
   displayed name is derived: today it renders "Comparative Evidence Index (variable
   domain only)"; when data improves the qualifier disappears with no rename.

   This also generalizes where a rename cannot — VHH has no light chain, scFv have
   linkers, bispecifics have two arms, and §37 already warns its weights are
   format-specific. Scope is that same parameter made explicit.

   **Never renormalize silently.** The geometric mean divides by
   `Σ(weight × reliability)`, so dropping Constant/Fc's 0.05 self-normalizes and the
   index looks complete. Coverage reporting is the only thing that makes it visible.

   Recorded in `CLAUDE.md` under "The index must state its own scope."
2. **Add isotype-level constant-region precedent.** Best value-per-effort remaining.
   From metadata already ingested: `CH1 Isotype: G1 667 · G4 200 · G2 77` (22 distinct);
   `VD LC: Kappa 847 · Lambda 117`. This supports honestly reporting "this VH's closest
   therapeutic matches are IgG1-format, kappa-paired" without a single constant-region
   residue — real evidence at class resolution.
3. **True full chains — defer past Gate 11.** IMGT publishes constant-region references
   (IGHG1, IGHG4, IGKC, IGLC); `FULL_HEAVY = user V-domain + named constant region` is
   legitimate for a candidate-assembly platform. But provenance must state "constant
   region from IMGT reference, not observed in any source record," and it can never be
   compared against observed full chains as though it were one.

**Resolved 2026-07-30 — is Constant/Fc a "critical" component?** It is not absolute;
it is per format. Define `required_components` for each construct format: Constant/Fc
is **not** required to claim `VARIABLE_DOMAIN_ONLY` scope, and **is** required to
claim `FULL_CONSTRUCT` scope. So Gate 11 does emit a number today, but may only ever
claim the scope it actually covered — which satisfies SPEC §37's "do not calculate the
index when critical evidence is missing" without declaring Fc non-critical by fiat.

---

## Known-uncalibrated items carried into Gate 9

- **`MINIMUM_WINDOWS_FOR_FRACTION = 20` is a provisional default**, not a calibrated
  threshold — same standing as the 0.70 coverage floor. It decides the status of
  41.2% of the corpus, so it is a consequential guess.
- **Window coverage needs exclusion rules, not percentiles, for CDR-H3 and the
  saturated frameworks.** See Gate 6 above.
- **Short components still lack a length-matched empirical null** (SPEC §21). Raw
  window counts are returned under every status so Gate 9 can build one.
- **Reference-population bias.** SAbDab2 is structure-derived — these are antibodies
  someone chose to crystallize. That population is not representative of natural
  repertoire. Layer 2 percentiles are percentiles *of solved structures*. Defensible
  if labeled as such; indefensible if presented as repertoire percentiles.
- **The 0.70 coverage floor** is a defensible default, not a calibrated value. Short
  frameworks (FR-H4, ~11 residues) and full domains plausibly warrant different floors.
  Confirmed to also misfire for germline matching (see "Germline reference data" above):
  ~210–240-residue VH components can never clear a query-coverage floor calibrated
  against a ~98-residue germline gene, regardless of match quality.
- **Recall was measured only at `top_k=10`.** Reporting more matches per component
  requires re-measuring; it will be worse.
- **Coverage floor and shortlist threshold interact** (stage-1 recall 0.975 vs
  end-to-end 0.958 for VH). Tune together, not independently.
- **12 queries per type** establishes the shape of the problem, not a tight confidence
  interval on 0.958.

## Gate 8 — observed pairing precedent ✅

Built 2026-08-01. Full detail in `reports/GATE8_PAIRING_PRECEDENT.md`; what
follows is what a future session must not re-derive.

### Findings that change later gates

- **54% of observed pairs rest on a single source record** (3,060 of 5,686). The
  2.12 records-per-pair average is skewed by a small deep tail — 98 pairs have
  more than 8. Any weighting by corroboration depth is weighting a long tail, not
  a typical case. **Gate 9 must know this before treating observation count as a
  calibration input.**
- **The two sources are almost disjoint: 4,562 SAbDab2-only, 1,122
  Thera-SAbDab-only, and 2 in both.** They are not overlapping views of the same
  evidence. Dropping either loses most of its contribution rather than a
  duplicated subset, and a dedup rule assuming cross-source mirroring has almost
  nothing to merge.
- Structure support across unique pairs: 4,852 backed, 834 stated absent, 0
  unknown.

### SPEC §35 defines no weights, so no composite exists

Unlike §26, §35 lists eight "provisional pair-match dimensions" with no formula.
`PairingPrecedentResult` therefore has **no `score` field at all** — nothing is
withheld because nothing is defined. Inventing weights would fabricate the spec.
Five dimensions are computable; three are blocked: germline-pair precedent (no
germline data), combined CDR-pattern similarity (needs the §21 length-matched
null), and format relevance (Gate 10.5).

**One germline reference set unblocks two gates** — Gate 9's calibration and this
dimension's germline-pair term. **One Gate 10.5 unblocks a dimension in both §26
and §35.**

### Decisions made here

- **Joint coverage is `min(heavy_coverage, light_coverage)`, not the mean.** A
  pair cannot be better covered than its weaker chain. **Provisional** — the
  choice of `min` is not calibrated.
- **A pair's relationship is the worse of its two chains'.** A kappa/lambda
  substitution downgrades the whole pair, so an exact heavy cannot carry a
  related-class light into same-class evidence.
- **No light chain returns `not_applicable`, never `not_found`.** A single-domain
  construct has no light chain by design; reporting a missing result would make a
  deliberate architecture indistinguishable from an unobserved pair.

### Presentation — decided for Gate 11, recorded here while the data is fresh

*(reviewed 2026-08-02)*

Gate 11 should present pairing evidence as **one card per observed antibody pair**,
not as a similarity table. A table hides the single thing this dimension uniquely
provides: *each row is one complete antibody whose two chains were seen together*.
Two independent search results can never establish that, so burying it wastes the
gate. Call it an "observed antibody pair", not a "pair" — the latter reads as *the
two numbers* rather than *the two chains*.

Engine 1.1.0 promoted the fields such a card needs off `PairRecord.metadata` and
onto each match: names, targets, statuses, formats, PDB entries, methods,
resolutions, plus a convenience `display_label`. Presentation stays deferred; the
data model does not need to be.

**Three constraints measured on the real corpus, all of which break the obvious
card design:**

- **80% of pairs have no drug name and no target.** 4,562 of 5,686 are structural
  records — a PDB entry, method and resolution, no name. A card with fixed
  Name/Target slots renders empty for four results in five. Those pairs are not
  missing an identity; they have a *different* one, and the card must adapt to
  which fields are populated rather than assume.
- **One chain pair can carry many names and many targets.** Trastuzumab's exact
  pair is backed by 11 records: 11 distinct programmes, **6 different targets**,
  and mixed development statuses. They are bispecifics sharing a trastuzumab arm.
  A single "Target: HER2" field would show one of six. This is why the promoted
  fields are lists.
- **Identity is not similarity, and coverage cannot be dropped.** The reported
  number is `global_identity` — exact matches over the longer sequence. Labelling
  it "similarity" would tell a scientist conservative substitutions are counted;
  they are not. And identity without coverage is the trap `approximate.py` exists
  to close: a random sequence once scored identity 1.000 at 4% coverage. Any card
  showing an identity must show joint coverage beside it.

Also: label an arm as an arm. A Thera-SAbDab pair is one arm of a construct, and
a bispecific has two with different specificities. "Observed Antibody: Evusheld"
flattens exactly the distinction that the arms-sharing-structures defect came
from. `display_label` qualifies with the arm only when more than one exists.

### The test that matters, and the mistake in writing it

`test_independent_heavy_and_light_hits_do_not_make_a_pair`. The first version
asserted an empty result and failed — correctly. Given a corpus of
`Alpha(h_A, l_A)` and `Beta(h_B, l_B)`, querying `h_A` with `l_B` legitimately
returns *both* real pairs, because each half-matches the query. The invariant is
not "nothing is found"; it is that the fabricated `(h_A, l_B)` combination never
appears. A companion test asserts the general form — every emitted pair must
exist in the index — which fails for any manufacturing route, not just the
anticipated one.

## The pair count is two different numbers

*(resolved 2026-07-31, confirmed by measurement 2026-08-01)*

This file previously said "12,032 observed pairs". Measured on the rebuilt corpus:

- **12,032** total `PairRecord` entries — what `IngestStats.observed_pairs` counts
- **5,686** unique `(heavy_checksum, light_checksum)` keys — what `len(PairIndex)` returns

Both are correct; they count different things, at ~2.1 records per unique pair
because one antibody appears across multiple structures and sources. **Gate 8 must
use 5,686**, or the same pair is counted up to a dozen times as independent
evidence — the deduplicate-before-counting rule (SPEC §15) applied directly. The
extra records are still useful as independent-observation *support*, but they are
not additional pairs.

### The four levels evidence attaches to — settle this before writing Gate 8

Gate 7 built a three-level store: **therapeutic entity** (the drug or programme)
→ **arm** (one binding specificity) → **sequence component** (VH, VL, CDR,
framework). Evidence attaches at the level it actually supports and never
propagates up or sideways — that rule is what fixed the bispecific defect where
each arm claimed the other's structures.

Gate 8 needs a fourth level those three do not express: the **source record**.
Pairing evidence is defined by co-occurrence *within a single source record*
(CLAUDE.md), so the record is the unit that licenses a pairing claim, and it is
not derivable from the entity/arm/component hierarchy. Two chains can both belong
to the same therapeutic entity and still have no source record linking them.

Two consequences, both easy to get wrong and expensive to unwind later:

- **An arm is not a pair.** 45 arms in Thera-SAbDab are heavy-only — every
  single-domain therapeutic among them. They are legitimate arms carrying real
  metadata, and they are zero evidence about pairing. A Gate 8 count that walks
  arms rather than source-linked chain pairs will silently include them.
- **Entity-level identity does not license a pairing claim.** Joining a heavy
  chain and a light chain because the same drug name appears on both is exactly
  the synthesized pair CLAUDE.md prohibits. The link must come from the record.

Built 2026-07-31. Full detail in `reports/GATE7_THERAPEUTIC_PRECEDENT.md`; what
follows is what a future session must not re-derive.

### Two defects found in the path SPEC §26 scores

- **794 of 1,131 therapeutics falsely claimed an experimental structure.** The
  sequence and metadata columns use the sentinel `na`; the three structure columns
  use `None`. `_parse_na` only knew `na`, so `"None"` survived as a truthy string
  and `structure_available` was `True` for every therapeutic. Now 337 structure-backed,
  794 not. The test fixture used `na` in those columns too, which is why 21 passing
  tests never caught it — the fixture now carries the real sentinel.
- **A bispecific's arms shared each other's structures.** Semicolon separates arms
  within a structure column, slash separates structures within an arm. The old code
  OR-ed all three columns and flagged both arms alike. Under schema drift the column
  now resolves to **unknown**, never to a guessed arm.

- **A blank structure cell asserted that no structure exists** *(found in review,
  parser 1.2.0)*. A third missing-data shape: an empty cell returned the same value
  as the explicit `None` sentinel, so `has_structure_record` answered `False` and
  SPEC §26's structure term scored 0.0 from a cell containing nothing. Defect 1
  inverted — absent evidence read as present there, unknown read as absence here.
  Blank now resolves to **unknown**. No impact on the current file (0 blanks before
  and after; the 337/794 split is unchanged), so this is defensive.

- **Structure tokens are validated by shape** *(parser 1.3.0)*. The three defects
  above were each found after the fact, and each fix taught the parser one more
  spelling of "missing" — an approach that only ever handles shapes someone has
  already met. A token surviving sentinel removal must now look like a PDB
  reference, so `N/A`, `-`, `TBD` or any future invention resolves to unknown and
  is reported through `IngestStats.warnings` instead of becoming a scored
  structure id. Measured first: 3,399 values, 679 distinct tokens, only `None`
  and one double-slash artifact are not PDB-shaped, so the rule describes the
  real data. **Use this pattern in Gate 8** — pairing has the identical failure
  mode, where a truthy sentinel becomes a fabricated pair.

`PairRecord.structure_available` is `Optional[bool]`. **`None` must never be coerced
to `False`** — parser version 1.3.0, `corpus_schema_version` 5.0.0.

Snapshot-specific counts are pinned in
`tests/integration/test_therapeutic_snapshot_counts.py`, guarded by the source
SHA-256: a Thera-SAbDab update skips the module with re-measurement instructions
rather than failing as a code regression.

### The therapeutic entity store

Metadata used to live only on `PairRecord`, so it was reachable only for arms with
*both* chains — **45 heavy-only arms, every single-domain therapeutic among them,
had no recoverable name, target, or status at all.** SPEC §26 needs those for any
supported component, so `Corpus` now holds 1,131 canonical therapeutics, 1,288 arms,
6,252 alias keys, and 10,186 component→therapeutic links (sub-regions included, so a
CDR resolves to its therapeutic).

### Aliases: 70 remain ambiguous, and that is correct

`Herceptin` → `Trastuzumab`. An alias claimed by two therapeutics is registered to
neither and counted; a poisoned set stops a third claimant re-registering it.
The surviving 70 span 22 pairs of **co-administered combination products** —
`Evusheld (combination of cilgavimab and tixagevimab)` genuinely denotes two
entities. Do not "fix" these.

Two splitting rules exist because their absence caused real collisions: semicolons
inside parentheses do not delimit (a combination name fragmented into a bare
`Bremzalerbart` that collided with the real drug), and `TBC`/`na`/`ND`/`NFD`/
`Unknown`/`None` are not names.

### The composite score is withheld — this is not an omission

`score` is `None`, `score_status = insufficient_evidence`. SPEC §26 weights
*calibrated* similarity at 0.45 and format relevance at 0.10, and neither term
exists. The four computable sub-terms are emitted as **raw proportions, never
weighted** — applying weights to a partial set yields something indistinguishable
from a composite that merely scored low.

**Gate 9 does not unblock the composite by itself.** It supplies similarity (0.45);
format relevance (0.10) stays blocked, leaving 0.90 of available weight. Emitting a
composite over that means dividing by 0.90 — silent renormalization, which this
project forbids because the result looks complete. An earlier draft of the Gate 7
report said the composite "follows" from Gate 9 alone; it does not, and that
sentence has been corrected rather than left for a future session to act on.

**Format relevance is unowned and is not Gate 9 work.** It needs the candidate's
declared construct format — which the submission schema does not collect from the
user at all, an input-contract change — plus a controlled vocabulary over a
120+-value free-text column, which is curation. Neither falls out of building
reference distributions.

**Resolved 2026-07-31 — Gate 10.5, sequenced before Gate 11.** Format relevance was
unowned: Gates 8, 9, 10 and 11 between them do not produce it, so "it will get done
eventually" was not true of it the way it is of the 0.45 similarity term. Left
unowned, all gates would complete and the composite would ship computed over 0.90
of its own definition with nothing flagging the gap.

It is scheduled before Gate 11 because Gate 11 is where the composite is actually
emitted. Landing it first means all six terms exist at assembly time and the
"publish an incomplete score?" question never arises. It has no dependency on
pairing, calibration or confidence work, so the placement is free — it is a
submission-schema field plus a curation pass, and it can move earlier if
convenient.

**If that sequencing ever slips**, the fallback is already specified and must not be
reinvented: CLAUDE.md's coverage contract (`coverage_scope`, `components_present`,
`components_missing`, `weight_set_version`) is what an incomplete composite has to
declare. What stays prohibited in every case is dropping the term and renormalizing
over the remaining 0.90 without saying so.

### Carried into later gates

- **29 therapeutics are numbered variants of one name** (`Lunsekimig1/2/3`). The
  source splits some drugs across rows, so `therapeutics_matched` can overcount by
  up to 2. Deliberately not auto-merged: a numeric-suffix rule could merge genuinely
  distinct drugs. Needs an explicit canonicalisation decision, most likely in Gate 11.
- **`format_relevance` is unimplemented, not approximated.** `Format` stays free
  text, verbatim. Finer architecture classification (scFv vs whole mAb vs ADC) is
  separate, explicitly calibrated work — not an extension of the VHH keyword heuristic.
- **Clinical stage and development status stay metadata.** 197 approved, 315
  discontinued. An antibody reaching Phase 3 transfers no property to a sequence
  resembling it; a test asserts a discontinued programme lowers no measurement.
- **Both sources carry `license_id = "unspecified"`** — Thera-SAbDab *and* SAbDab2.
  SAbDab2's was invisible until `_snapshot_to_dict` stopped discarding the field
  (fixed 2026-07-31, `CACHE_FORMAT_VERSION` 1.2.0). Blocks public deployment and
  redistribution; does not block local analysis.

### Why Gate 9 moved up the critical path

Gate 7 emits everything except its calibrated term. Gate 8 will be in the same
position, and Gate 11's index consumes both. Calibration is now the single blocker
on two finished dimensions rather than a later refinement — and it still needs the
germline baseline discussed above.

## Minor leftovers

- Exhaustive search costs ~12 s **per search**, not per profile — a deliberate
  correctness-over-speed choice. The honest lever, if interactive latency matters, is
  opting whole domains into shortlisting (24.7× on VH) while surfacing `shortlisted`
  in the output.

  **Each dimension currently runs its own searches**, so a full profile repeats them:
  Gate 7 searches the VH, Gate 8 searches the same VH again plus the VL, and Gate 11
  will invoke both. Nothing is wrong — every dimension is correct in isolation — but
  the same alignment is paid for several times.

  Deliberately *not* fixed now, and not a trap for Gate 11: the fix is additive, an
  optional "hits already computed" argument on each `compute_*` function, so no
  existing signature has to change and nothing has to be decided in advance. Worth
  doing when profile assembly makes the cost visible, not before.
- 1,375 skipped rows and 72 numbering failures inherited from Phase 1a have never been
  individually reviewed.
- **Resolved 2026-08-10:** repo/package name mismatch. Package renamed `antigens` →
  `para_scope` (distribution name `para-scope`), matching the repo, ahead of making
  it public. 107 whole-word occurrences across 32 files, mechanical except
  `pyproject.toml`'s `name`/script-entry lines and the CLI usage string, which name
  the hyphenated distribution form rather than the importable module. `reports/`
  left untouched — those are dated records of what was true when each gate was
  built, not living documentation.

---

## Definition of done for a gate

State was lost once already, between sessions, because the roadmap lived only in
conversation. A gate is not complete until **all** of the following exist:

- [ ] Implementation, with versioned constants bumped where behavior changed
- [ ] Tests passing, and the new total recorded
- [ ] A report in `reports/` when the gate produced measurements or found a defect
- [ ] **This file updated** — status row, decisions made, blockers found or cleared
- [ ] Committed to git

The fourth item is the one that gets skipped when a session runs long. It is the
reason this file exists. Update it *as part of* the gate, not afterward.
