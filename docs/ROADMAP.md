# Version 2 Roadmap — Gate Status, Decisions, and Blockers

**This file is the source of truth for project state.** `CLAUDE.md` points here and
instructs any session to read it before planning work. If it disagrees with anyone's
recollection — human or model — this file is right.

Last verified: 2026-07-30, against commit on `main`, 150 tests passing.

Work proceeds in numbered gates. A gate is a unit of work that ends in a verifiable
artifact, not a time box. Gates 0–5 are complete.

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
| **6** | **Human-repertoire naturalness** | **next — unblocked, not built** | — |
| 7 | Therapeutic-precedent profile | ready, no new data needed | SPEC §26 |
| 8 | Pairing precedent | ready — PairIndex works, 12,032 observed pairs | SPEC §35 |
| 9 | Empirical reference distributions | ready | SPEC §19 Layers 2–4 |
| 10 | Evidence strength / confidence / completeness | — | SPEC §28 |
| 11 | Final profile assembly | — | SPEC §37 |

Corpus under test: 46,627 components, 28,835 approximate-eligible,
`corpus_schema_version` 4.0.0.

---

## Gate 6 — human-repertoire naturalness

### BioPhi was rejected, and is not coming back

BioPhi requires Python 3.9, conda, and a 2020-vintage ANARCII. The project runs
Python 3.13 with `anarcii>=2.0.8`. This is an unresolvable stack collision, not an
integration difficulty. Do not re-litigate it.

### The replacement (validated, not yet installed)

A 9-residue sliding-window matcher against human antibody sequences. Score is the
fraction of query windows observed in real human antibodies. **207 MB, not 22 GB**,
and installs on the existing Python without downgrading anything in the pipeline.

Validated before adoption:

| Sequence | Score |
|---|---|
| Humanized (trastuzumab) | 0.71 / 0.83 |
| Mouse | 0.58 |
| **Llama (VHH)** | **0.685** |

Humanized > mouse is the correct ordering; had it inverted, the metric would be junk.

**The llama result is the load-bearing limitation.** At 0.685 it is statistically
indistinguishable from the humanized 0.705. This metric therefore **cannot measure
"how human" a sequence is.** It measures only how many of its 9-residue windows
appear in human antibodies. It must never be surfaced, named, or documented in a way
that implies immunogenicity, safety, or humanization quality — that would violate
`CLAUDE.md` prohibited claims.

Two corrections to the source document, recorded so they are not reintroduced:
a separate Python 3.11 worker process is unnecessary — skip it; and the document's
example score of 77.7% was illustrative, not a target (the real figure was 70.5%).

### State at time of writing

Never implemented. The prior session ended on `"I'm implementing it now"` and then
diverted. `promb` is **not** in the venv, `pyproject.toml` still lists only `anarcii`
and `biopython`, and no naturalness module exists. Re-verify the package name and
install before building.

### Required behavior

The naturalness slot must exist in the profile from the start, returning an explicit
`database_unavailable` status until data is wired in. It must **never** default to a
neutral or midpoint value — a neutral number reads as "average humanness" when it
actually means "we did not look." That is the `missing evidence is not zero evidence`
constraint applied directly.

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

1. **Rename the index.** "Whole-Construct Comparative Evidence Index" overclaims —
   the analysis is variable-domain-only. **Variable-Domain Comparative Evidence Index**
   costs nothing and removes the overclaim. Do this regardless of everything else.
   Note this conflicts with the term fixed in `CLAUDE.md`'s terminology table; that
   table needs updating in the same change.
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

Unresolved and blocking Gate 11's output shape: **is Constant/Fc a "critical"
component?** SPEC §37 forbids calculating the index when critical evidence is missing.
With Constant/Fc weighted 0.05 and permanently empty, this definition determines
whether Gate 11 emits a number or "Insufficient evidence for a whole-construct index."

---

## Known-uncalibrated items carried into Gate 9

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

## Gate 7 constraints

- **Thera-SAbDab license is `unspecified`** — flagged in `build_thera_sabdab_snapshot`.
  Blocks public deployment and redistribution; does not block local analysis.
- **Clinical stage and development status stay metadata.** An antibody reaching Phase 3
  transfers no property to a sequence resembling it. These must never become a success
  or failure label.
- **`Format` free text is stored verbatim.** Finer architecture classification
  (scFv vs whole mAb vs ADC) is separate, explicitly calibrated work — not an
  extension of the VHH keyword heuristic.

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
