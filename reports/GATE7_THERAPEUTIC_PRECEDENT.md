# Gate 7 — Therapeutic Precedent (SPEC §26)

Built 2026-07-31. Corpus: 46,627 components, `corpus_schema_version` 5.0.0,
`CACHE_FORMAT_VERSION` 1.2.0. Suite: **242 tests passing** (was 183 pre-gate).

*Revised 2026-07-31 after review: a third missing-data sentinel was found (below),
the "Gate 9 completes the composite" claim was wrong and is corrected, and the
snapshot-specific counts are now pinned in `tests/integration/`.*

Answers one question: *has this component, or a similar one, been observed in an
antibody that entered therapeutic development, and what does the curated record
say about it?* Nothing about how the candidate will behave.

---

## Two defects found before any new code was written

Both were in the path SPEC §26 scores directly, so they are reported first.

### 1. Every therapeutic claimed an experimental structure — 794 of them falsely

`PairRecord.structure_available` was `True` for all 1,131 therapeutics.

The sequence and metadata columns write missing data as the literal string `na`.
The three structure columns do not — they write `None`. `_parse_na` only
recognised `na`, so `"None"` survived as a truthy string and every structure
column tested as present.

Measured on the real file:

| | before | after |
|---|---:|---:|
| structure-backed | 1,131 | **337** |
| no structure record | 0 | **794** |
| unknown | 0 | 0 |

This is a fabricated evidence flag, not a display bug: SPEC §26 weights
"direct pair or structure-linked record support" at 0.10, and the flag asserted an
experimental record exists where the source says none does. The test fixture had
used `na` in those columns too, which is why 21 passing tests never caught it —
the fixture now carries the real sentinel.

### 2. A bispecific's two arms shared each other's structures

Within a structure column a semicolon separates arms (`None;8z8v:B` = nothing for
arm 1, one structure for arm 2) and a slash separates structures within an arm.
The old code OR-ed all three columns and gave both arms the same flag, so arm 1
was credited with arm 2's experimental records — two different specificities.

Alignment is exact on the real file: all 157 two-arm rows are semicolon-split,
none of the 974 single-arm rows are, 3,393/3,393 values align. Under schema drift
the column resolves to **unknown** rather than being mapped by guesswork; a wrong
arm is worse than an absent one.

`structure_available` is now `Optional[bool]`, and `None` must never be coerced to
`False`.

### 3. A blank structure cell asserted that no structure exists

Found during review, after the two above. Parser **1.2.0**.

There is a third missing-data shape and it means something different again: an
empty cell. `split_structure_column("")` returned `[[]]` — byte-identical to what
the explicit `None` sentinel returns — so `has_structure_record` answered `False`
and SPEC §26's structure term scored **0.0 off a cell containing nothing**.

This is defect 1 inverted. There, an unrecognised sentinel made absent evidence
read as present; here, an empty value made unknown evidence read as a positive
statement of absence. Both fabricate a claim the source never made, and both are
scored. A blank now resolves to *unknown*, which keeps the term out of the profile
entirely rather than contributing a zero.

An empty token anywhere in a multi-arm value makes the whole column unknown, on
the same reasoning as schema drift: a half-readable arm mapping cannot be
distinguished from a complete one, and a wrong arm is worse than an absent one.

**Impact on the current file: none.** This snapshot has no blank structure cells —
`unknown` is 0 before and after, so the 337/794 split is unchanged. The fix is
defensive. It is reported because the previous two defects were also invisible
until a fixture happened to contain the real sentinel, and the unit fixtures still
cannot see a shape the sample rows do not include.

---

## What the gate added

### Therapeutic entity store (`corpus.py`)

Therapeutic metadata previously existed only inside `PairRecord.metadata`, which
meant it was reachable only for arms that had *both* chains. **45 heavy-only arms
— every single-domain therapeutic among them — carried no recoverable name,
target, format, or development status at all.** SPEC §26 requires those as display
fields for any supported component, so they cannot depend on an arm happening to
be paired.

| | |
|---|---:|
| canonical therapeutics | 1,131 |
| arms | 1,288 |
| heavy-only arms (metadata previously unreachable) | 45 |
| components linked back to a therapeutic | 10,186 |
| resolvable alias keys | 6,252 |

`therapeutics_for_component(record_id)` is the reverse edge the dimension needs:
an alignment hit knows its record id and nothing else. Sub-regions are linked too,
so a CDR resolves to its therapeutic, not just the whole variable domain.

### Alias resolution (SPEC §15.5)

`Herceptin` → `Trastuzumab`; case-insensitive; development codes included.

An alias claimed by two therapeutics is registered to **neither** and counted,
rather than letting whichever row parsed last win. A poisoned-alias set prevents a
third claimant from re-registering it.

**70 aliases remain ambiguous across 22 therapeutic pairs, and that is correct.**
They are co-administered combination products — `Evusheld (combination of
cilgavimab and tixagevimab)` genuinely denotes two entities, so resolving it to
either would be wrong.

Fixing the splitting rule removed 3 spurious collisions of two kinds:

- **Parenthesised lists.** `REGN5713-5714-5715 (combination of Atisnolerbart;
  Bremzalerbart; and Umesolerbart)` was split on its internal semicolons, yielding
  a bare `Bremzalerbart` that collided with the real therapeutic of that name, plus
  `and Umesolerbart)`, which is not a name. Semicolons inside parentheses no longer
  delimit.
- **Placeholders.** `TBC` was listed as an alternative name by two therapeutics and
  so poisoned itself. `TBC`, `na`, `ND`, `NFD`, `Unknown`, `None` are no longer
  registered as names.

### The dimension (`therapeutic_precedent.py`)

`THERAPEUTIC_PRECEDENT_VERSION = "1.1.0"`, `CALIBRATION_VERSION = "uncalibrated"`.

Trastuzumab VH against the real corpus:

```
status                available
therapeutics_matched  15   (raw observations: 15)
exact matches         11
best_global_identity  1.0000
score                 None  [insufficient_evidence]
computable_terms      coverage 1.0 · exact_component_class_confirmation 1.0
                      direct_pair_or_structure_support 1.0 · metadata_completeness 1.0
```

All 15 are HER2-targeted. The 11 exact matches are **11 distinct therapeutic
programmes sharing a byte-identical VH** — bispecifics and derivatives built on
the trastuzumab arm — not 11 redundant records of one drug. This is why entities
and observations are reported as separate numbers.

Ties at identity 1.000 are the normal case here, not an edge case, so ordering
breaks ties by name; without that the "top match" could differ between two runs
over the same snapshot.

### Each match names its own source (engine 1.1.0)

`source_databases` and `source_pointers` are carried on every `TherapeuticMatch`.
The corpus always held this on the component record, but reaching it meant taking
`matched_record_id` back through `corpus.get(...).source_pointers` — a hop that
fails silently when fumbled, since the stored name is `Thera-SAbDab` and filtering
on a plausible `thera_sabdab` yields zero matches rather than an error. That exact
mistake produced a false "0 therapeutic-linked hits" while this gate was being
validated.

Licensing is the sharper reason. SPEC §13/§14 require a source's licence to travel
with the evidence derived from it, and licence is recorded per `SourceSnapshot`, so
a displayed match that does not name its source cannot be joined to its terms. Both
sources are `unspecified` today, so nothing is currently mis-stated; the join stops
being cosmetic the moment one source is cleared and the other is not.

A component present in both sources reports both pointers — provenance is every
source it was observed in, not the one that happened to be reached first.

---

## The composite score is withheld, deliberately

`score` is `None` with `score_status = insufficient_evidence`.

SPEC §26's composite weights **calibrated** component similarity at 0.45 and format
relevance at 0.10. Neither term exists:

| term | weight | state | owned by |
|---|---:|---|---|
| calibrated component similarity | 0.45 | **blocked** | Gate 9 (SPEC §19 layers 2–4) |
| coverage | 0.15 | computable | — |
| exact component-class confirmation | 0.15 | computable | — |
| format relevance | 0.10 | **blocked** | **no gate yet — see below** |
| direct pair / structure support | 0.10 | computable | — |
| metadata completeness | 0.05 | computable | — |

Substituting raw identity into the calibrated slot would publish an uncalibrated
number under a calibrated formula's name — reproducible, but reproducibly wrong,
and indistinguishable downstream from a finished score. The four computable
sub-terms are emitted as **raw proportions, not weighted contributions**: applying
the weights to a partial set would produce something that looks like a composite
which merely scored low.

### Gate 9 does not unblock the composite on its own

An earlier draft of this report closed by saying Gate 9 supplies the similarity
term and the composite then follows. That was wrong, and wrong in a direction the
project's rules specifically prohibit.

After Gate 9 lands, format relevance is *still* blocked. The available weight is
0.90, so emitting a composite means dividing by 0.90 — the silent renormalization
CLAUDE.md forbids, where an omitted component self-normalizes and the index looks
complete. A reader cannot tell a full-coverage score from one that quietly dropped
a tenth of its definition.

**Format relevance is not Gate 9 work and must not be filed there.** It needs two
inputs Gate 9 does not produce:

1. **The candidate's declared construct format**, which the submission schema does
   not currently collect from the user at all. That is an input-contract change.
2. **A controlled vocabulary** over Thera-SAbDab's free-text `Format` column
   (120+ distinct values). That is curation work, not calibration work.

Neither is a by-product of building empirical reference distributions. Filing this
under Gate 9 is how it becomes invisible and ships as a renormalized 0.90.

### Resolved — Gate 10.5, before Gate 11

Format relevance now has an owner: **Gate 10.5**, sequenced ahead of Gate 11.

The distinction that made this urgent: the 0.45 similarity term was always going to
be built, because Gate 9 exists and produces it. Format relevance was not — no
remaining gate produces it, so it would not have arrived on its own. Every gate
would have completed and the composite would have shipped over 0.90 of its
definition with nothing anywhere reporting the shortfall.

Gate 11 is where the composite is emitted, so landing 10.5 first means all six
terms exist at assembly and no incomplete score is ever published. The placement
costs nothing — it depends on no other gate, being a submission-schema field plus a
curation pass.

**Fallback if the sequencing slips**, already specified in CLAUDE.md and not to be
reinvented: an incomplete composite declares `coverage_scope`, `components_present`,
`components_missing` and `weight_set_version`. Dropping the term and renormalizing
over the remaining 0.90 without saying so stays prohibited regardless.

---

## Short components are not excluded from exact matching

CDRs are excluded from *alignment* — unmodified full-length scoring is invalid for
them (SPEC §21) and no length-matched null exists until Gate 9. They are **not**
excluded from *exact* matching: an identical CDR-H3 needs no chance correction
because there is no alignment to score.

Trastuzumab's IMGT CDR-H3 (`SRWGGDGFYAMDY`) returns **14 therapeutics** by exact
match, with `approximate_search_status = unsupported_component` and the reason
recorded. A naive implementation would have returned "not found" for every CDR.

Note the boundary: IMGT CDR3 opens at position 105, so it is `SRWGGDGFYAMDY`, not
the `WGGDGFYAMDY` a Kabat-shaped intuition suggests. The first probe written for
this gate used the wrong string and returned zero hits.

---

## Status routing

| condition | status |
|---|---|
| matches found | `available` |
| searched, nothing matched | `not_found` |
| corpus holds no therapeutic records | `database_unavailable` |
| empty sequence | `unsupported_component` |
| alignment raised | `analysis_failed` (per-search, evidence still returned) |

Every non-`available` outcome leaves `score` unset. Never 0.0, never a midpoint,
never a corpus mean — a failed lookup must not read as a low score.

`database_unavailable` carries "This is not evidence of absence" in its reason, so
a corpus built without Thera-SAbDab cannot read as proof that no therapeutic
precedent exists.

---

## Known-uncalibrated and carried forward

- **29 therapeutics are numbered variants of one name** — `Lunsekimig1`,
  `Lunsekimig2`, `Lunsekimig3`. The source splits some drugs across rows, so
  `therapeutics_matched` can overcount those by up to 2. Not auto-merged: a
  numeric-suffix rule could merge genuinely distinct drugs, and SPEC §15.5 wants
  one entity per therapeutic, not one entity per name that looks similar. Needs an
  explicit canonicalisation decision.
- **`format_relevance` is unimplemented, not approximated.** The `Format` column
  stays free text, verbatim, per the Gate 7 constraint in ROADMAP.
- **Both source licences remain `unspecified`** — SAbDab2 and Thera-SAbDab. Blocks
  public deployment and redistribution; does not block local analysis.
- **Clinical stage and development status are metadata.** 197 approved, 315
  discontinued. Neither is folded into any numeric term, and a test asserts a
  discontinued programme lowers no reported measurement.

## Reproduce

```bash
.venv/bin/python scripts/build_combined_corpus_cache.py
.venv/bin/python scripts/measure_therapeutic_precedent.py
.venv/bin/python -m pytest -q
```
