# Version 2 Roadmap — Gate Status, Decisions, and Blockers

**This file is the source of truth for project state.** `CLAUDE.md` points here and
instructs any session to read it before planning work. If it disagrees with anyone's
recollection — human or model — this file is right.

Last verified: 2026-07-31, against commit on `main`, 222 tests passing.

Work proceeds in numbered gates. A gate is a unit of work that ends in a verifiable
artifact, not a time box. Gates 0–7 are complete.

---

## Status

| Gate | Scope | Status | Artifact |
|---|---|---|---|
| 0 | Baseline audit | ✅ | `reports/V2_BASELINE_AUDIT.md` |
| 0.5 | Persistence layer | ✅ | `src/antigens/persistence.py` |
| 1 | End-to-end search validation | ✅ | `reports/END_TO_END_SEARCH_VALIDATION.md` |
| 2 | Thera-SAbDab ingestion | ✅ | `reports/THERASABDAB_INGESTION.md` |
| 3 | Two-stage approximate search | ✅ | `reports/TWO_STAGE_SEARCH.md` |
| 4 | Deterministic feature engine | ✅ | `src/antigens/features.py` |
| 5 | Motif observation engine | ✅ | `src/antigens/motifs.py` |
| 6 | Human repertoire window coverage | ✅ | `reports/GATE6_WINDOW_COVERAGE.md` |
| 7 | Therapeutic-precedent profile | ✅ | `reports/GATE7_THERAPEUTIC_PRECEDENT.md` |
| 8 | Pairing precedent | ready — PairIndex works; see pair-count note below | SPEC §35 |
| **9** | **Empirical reference distributions** | **ready — now blocks two gates' composites** | SPEC §19 Layers 2–4 |
| 10 | Evidence strength / confidence / completeness | — | SPEC §28 |
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

## Germline reference data — unresolved, and larger than it looks

`EvidenceClass.GERMLINE` is declared in `models.py` and **nothing populates it.**
That is its only occurrence in the codebase. The numbering layer captures no V/D/J
assignment either. ANARCII was checked directly: it has no germline capability of any
kind — no V-gene, no allele, numbering only.

**Why this outranks the `FULL_HEAVY` gap in impact:** Database-Relative Novelty is
close to meaningless without a germline baseline. "Unusual" relative to what?

### Options considered

- **A — IgBLAST.** Field-standard, authoritative, what reviewers expect. But it is an
  NCBI binary rather than a Python package: platform-specific install and separate
  germline databases to configure. Reintroduces exactly the external-tool fragility
  that rejecting BioPhi avoided.
- **B — build on the existing alignment engine. ← recommended.** Germline V genes are
  just sequences, and the project already has a pairwise aligner, a component index,
  and IMGT numbering. Load a human germline V reference set (a few hundred sequences,
  not gigabytes), index it, align the query. The best hit yields both the gene call
  and **percent identity to nearest germline** — the standard measure of drift from
  the inherited starting point, and precisely the baseline Novelty currently lacks.
  Almost entirely code already written and tested. No new heavy dependency.
- **C — a package bundling germline data** (SADIE, PyIR). Most either wrap IgBLAST
  anyway or pull large dependency trees. Only worth revisiting if B proves harder
  than expected.

### Verify before implementing

Germline source and license, which is step one and not yet done:

- **IMGT/GENE-DB** — canonical, but restrictive redistribution terms. Same category of
  problem as the `unspecified` license already flagged on Thera-SAbDab.
- **OGRDB** — the AIRR Community's open germline database. Likely the cleaner
  licensing story, which matters if this ever becomes public. **Lean toward OGRDB.**

Per SPEC §11.2: do not silently merge OGRDB and IMGT assignments. Report reference
source, germline-set version, best assignment under each set, and conflicting allele
calls.

### Open question

Whether germline becomes its own gate ahead of Gate 9. It is a Gate 9 dependency:
calibrating novelty without a germline baseline calibrates against nothing.
A reduced scope — distance-to-nearest-germline-V only, no full V(D)J gene calls —
captures most of the value for Novelty at meaningfully less work. The gene call is
useful context, an add-on rather than the core.

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
- **Recall was measured only at `top_k=10`.** Reporting more matches per component
  requires re-measuring; it will be worse.
- **Coverage floor and shortlist threshold interact** (stage-1 recall 0.975 vs
  end-to-end 0.958 for VH). Tune together, not independently.
- **12 queries per type** establishes the shape of the problem, not a tight confidence
  interval on 0.958.

## Gate 8 — the pair count is two different numbers

*(resolved 2026-07-31)*

This file previously said "12,032 observed pairs". Measured on the rebuilt corpus:

- **12,032** total `PairRecord` entries — what `IngestStats.observed_pairs` counts
- **5,686** unique `(heavy_checksum, light_checksum)` keys — what `len(PairIndex)` returns

Both are correct; they count different things, at ~2.1 records per unique pair
because one antibody appears across multiple structures and sources. **Gate 8 must
use 5,686**, or the same pair is counted up to a dozen times as independent
evidence — the deduplicate-before-counting rule (SPEC §15) applied directly. The
extra records are still useful as independent-observation *support*, but they are
not additional pairs.

## Gate 7 — therapeutic precedent ✅

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

`PairRecord.structure_available` is `Optional[bool]`. **`None` must never be coerced
to `False`** — parser version 1.1.0, `corpus_schema_version` 5.0.0.

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
from a composite that merely scored low. Gate 9 supplies the similarity term and
the composite follows with no recomputation.

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

- Exhaustive search costs ~12 s per full profile — a deliberate correctness-over-speed
  choice. The honest lever, if interactive latency matters, is opting whole domains
  into shortlisting (24.7× on VH) while surfacing `shortlisted` in the output.
- 1,375 skipped rows and 72 numbering failures inherited from Phase 1a have never been
  individually reviewed.
- Repo is `Para-Scope`; Python package is `antigens`. Renaming gets more annoying each
  gate. Undecided.

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
