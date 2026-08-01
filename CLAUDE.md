# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

Antigens is an **explainable antibody component comparison and candidate-assembly platform**. Users submit existing antibody sequences, or assemble a candidate construct from selected components. The platform validates the arrangement, compares every component against equivalent known components, and produces a transparent comparative evidence profile with full provenance.

It is **not** a binding predictor, docking engine, affinity model, or antibody-design AI.

Full specification: [docs/SPEC.md](docs/SPEC.md). Read the relevant section before implementing any module — the sections are numbered and the spec is the source of truth for anything not covered here.

## Non-negotiable constraints

These are correctness requirements, not style preferences. Violating one is a defect regardless of whether tests pass.

### Comparison integrity

- **Only compare like with like.** Heavy chain against heavy chain, CDR-H3 against CDR-H3. Never CDR-H3 against CDR-L3; never VH against VL. Related-class comparisons (VH↔VHH, kappa↔lambda) must be explicitly labeled and down-ranked, never silently merged into same-class evidence.
- **A variable-domain match is never a full-chain match.** Keep variable-only and full-chain records distinguishable at every layer.
- **Paired and unpaired evidence never merge.** Observed pairing requires both chains linked in the *same source record*. Never synthesize a pair by joining an independent heavy hit with an independent light hit.
- **Short components require empirical chance correction.** CDRs, linkers, and junctions need a length-matched empirical null model (spec §21). Never apply unmodified full-length BLAST scoring to them.

### Numbering

- IMGT is the canonical internal coordinate system. Kabat, Chothia, Martin, and AHo are display conversions only.
- Never combine regions extracted under different numbering definitions without explicit remapping.
- Always store: original coordinates, IMGT coordinates, insertion codes, tool name, tool version, and numbering confidence.

### Evidence handling

- **Missing evidence is not zero evidence.** Use explicit statuses (`not_found`, `not_applicable`, `insufficient_evidence`, `analysis_failed`, `database_unavailable`, `access_restricted`, `unsupported_component`). A failed database call must never lower a biological score.
- **Never silently modify input.** No residue substitution, deletion, or imputation. Preserve the original input verbatim and record every normalization action.
- **Every score must be reproducible** from stored raw evidence records plus the named formula version and calibration version.
- **Scoring code never queries a live external service.** All evidence comes from versioned, checksummed local snapshots. A database update must never silently alter a prior analysis.
- **Deduplicate before counting.** Multiple structures of one antibody, one clone expanded across a repertoire, and the same experiment mirrored in several databases each count once (spec §15).

### First-version structural boundary

Sequence-first and metadata-aware only. **Do not implement**: coordinate parsing (PDB/mmCIF), structure prediction, structural alignment, RMSD, CDR loop geometry, VH–VL orientation, molecular surfaces, coordinate-derived paratope extraction, or 3D visualization.

Structure *metadata* (PDB ID, experimental method, resolution, antigen-complex status, curated contact-residue lists) is in scope. Structure *computation* is not.

### Privacy

User sequences may be proprietary. Never transmit one to an external service without explicit per-submission consent. No sequence content in application logs.

## Required terminology

Use these exact terms in code, identifiers, schemas, comments, and user-facing copy.

| Use | Never use |
|---|---|
| Candidate antibody construct, candidate assembly | "created an antibody", "functional antibody" |
| Paratope — a numbered residue *set*, not a contiguous string | treating a paratope as a CDR-like sequence |
| Observed Repertoire Precedent | "natural precedent" |
| Database-Relative Novelty | unqualified "novelty" |
| Whole-Construct Comparative Evidence Index | binding / success / functional / safety / efficacy / clinical / manufacturability score |
| Human Repertoire Window Coverage — a count of 9-mer windows found in a reference set | "humanness", "naturalness", humanization or immunogenicity score |

## Prohibited claims

Never state or imply that a candidate will bind, has a target, will fold, express, pair, remain stable, avoid aggregation, be soluble, be manufacturable, be safe, be effective, have low immunogenicity, or succeed experimentally or clinically. Similarity to an approved therapeutic transfers none of its properties.

Any user-facing summary score must carry this text verbatim:

> This index summarizes comparative database evidence. It is not a probability of folding, expression, pairing, binding, stability, manufacturability, safety, efficacy, or experimental success.

### The index must state its own scope

The term stays "Whole-Construct Comparative Evidence Index" — it names the spec's
concept (SPEC §37) and the platform's target, not whatever subset of components
happens to be populated today. Do **not** rename it to match a temporary data gap.

Instead, every emitted index carries its coverage as required, machine-readable
fields, and the displayed name is *derived* from them:

- `coverage_scope` — `VARIABLE_DOMAIN_ONLY` | `FULL_CONSTRUCT`
- `components_present` — the slots that contributed
- `components_missing` — each with an explicit status and reason
- `weight_set_version` — §37's weights are format-specific and provisional

So an index computed today renders as "Comparative Evidence Index (variable domain
only)" because `coverage_scope` says so. When constant-region data lands, the
qualifier disappears on its own with no rename and no discontinuity in the record.

Two rules this exists to enforce:

- **Never renormalize silently.** The §37 geometric mean divides by
  `Σ(weight × reliability)`, so omitting a component self-normalizes and the index
  *looks* complete. Coverage reporting is the only thing that makes the omission
  visible. `Evidence coverage` and `Missing critical evidence` are already required
  §37 outputs — emit them.
- **"Critical" is per format, not absolute.** Define `required_components` for each
  format. Constant/Fc is not required to claim variable-domain scope; it is required
  to claim full-construct scope. The index may compute, but may only ever claim the
  scope it actually covered.

## Current state — read this before proposing work

**[docs/ROADMAP.md](docs/ROADMAP.md) is the source of truth for what is built, what
is next, and what is blocked. Read it at the start of any session before planning,
estimating, or claiming a gate is incomplete.** It is kept current because updating
it is part of a gate's definition of done — see the checklist at the end of that
file. If it disagrees with your recollection, the file is right.

Work proceeds in numbered gates. Gates 0–7 are complete; Gate 9 is next (it now
blocks the composite scores of both Gate 7 and Gate 8).

## Architecture

Python 3.13, `src/` layout, package `antigens` (repo is named Para-Scope; the
package name has not been renamed). Dependencies are deliberately few: `anarcii`
for numbering, `biopython` for alignment and residue chemistry, `promb` for the
bundled human repertoire peptide reference set (MIT; data ships inside the wheel,
so nothing is fetched at scoring time).

```
src/antigens/
  models.py            Component/evidence dataclasses, EvidenceStatus, enums
  numbering.py         ANARCII wrapper, IMGT coordinates, region extraction
  corpus.py            In-memory corpus, ingestion, CORPUS_SCHEMA_VERSION
  persistence.py       Snapshot save/load, CACHE_FORMAT_VERSION, StaleSnapshotError
  features.py          Deterministic sequence features (Gate 4)
  motifs.py            Motif observations in IMGT coordinates (Gate 5)
  window_coverage.py   Human repertoire 9-mer window coverage (Gate 6)
  therapeutic_precedent.py  Therapeutic precedent evidence, SPEC §26 (Gate 7)
  indexes/             exact, positional, pair, approximate, shortlist
  parsers/             sabdab2, thera_sabdab
```

Two invariants that cut across the code and are easy to break:

- **Versioned constants gate reproducibility.** `CORPUS_SCHEMA_VERSION`,
  `CACHE_FORMAT_VERSION`, `ALIGNMENT_SCORING_VERSION`, `FEATURE_ENGINE_VERSION`,
  `MOTIF_RULESET_VERSION`, `WINDOW_COVERAGE_VERSION`,
  `THERAPEUTIC_PRECEDENT_VERSION`.
  Adding state to `Corpus` requires bumping its schema version, or old snapshots
  deserialize silently missing attributes and fail far from the cause.
- **`data/` is gitignored and regenerable.** Corpus snapshots and the numbering
  cache never enter git. A fresh clone must rebuild via `scripts/` before the
  integration tests can run.

## Development commands

```bash
.venv/bin/python -m pytest              # full suite (~43 s, 254 tests)
.venv/bin/python -m pytest -q tests/    # unit only; tests/integration needs the corpus

# Corpus rebuild (required after a fresh clone; ~8.5 min cold, ~7 s cached)
.venv/bin/python scripts/build_v2_corpus_cache.py       # SAbDab2
.venv/bin/python scripts/build_combined_corpus_cache.py  # + Thera-SAbDab
.venv/bin/python scripts/benchmark_two_stage_search.py   # shortlist recall/latency
.venv/bin/python scripts/measure_window_coverage.py      # Gate 6 corpus distributions
.venv/bin/python scripts/measure_therapeutic_precedent.py  # Gate 7 therapeutic coverage
```
