# Para-Scope

**An explainable antibody component comparison and candidate-assembly platform.**

You submit an antibody sequence, or assemble a candidate construct from
individual components. Para-Scope compares every piece — variable domains,
CDRs, frameworks — against equivalent components in real, versioned antibody
databases, and returns what it found: identity, coverage, observed pairing,
therapeutic precedent, human-repertoire window coverage. Every number carries
its source record, its formula version, and an explicit status when evidence
is missing.

It doesn't predict whether an antibody will bind, fold, express, or succeed
clinically. It looks up what's already known and shows you, with receipts.

## What this is not

Not a binding predictor. Not a docking engine. Not an affinity model. Not an
antibody-design AI. Similarity to an approved therapeutic transfers none of
its properties, and the platform is built not to imply otherwise — see
[Design principles](#design-principles) below.

## Status

**Gates 0–8 of 11 complete.** Gate 9 (empirical calibration) is next, and is
the current blocker on two finished evidence dimensions.

| Gate | Scope | Status |
|---|---|---|
| 0–3 | Corpus ingestion, numbering, search (exact + two-stage approximate) | ✅ |
| 4–5 | Deterministic sequence features, motif observations | ✅ |
| 6 | Human repertoire window coverage | ✅ |
| 7 | Therapeutic precedent (SPEC §26) | ✅ |
| 8 | Observed heavy-light pairing precedent (SPEC §35.1) | ✅ |
| 9 | Empirical reference distributions / calibration | next |
| 10–11 | Confidence scoring, final profile assembly | planned |
| 12+ | Calibrated prediction layer ([details](#beyond-gate-11-a-calibrated-prediction-layer-planned-not-built)) | intent, unscoped |

**[docs/ROADMAP.md](docs/ROADMAP.md)** is the living source of truth — every
gate's decisions, defects found, and what a future session should not
re-derive. **[docs/SPEC.md](docs/SPEC.md)** is the full specification.
**[CLAUDE.md](CLAUDE.md)** is the concise ruleset an AI coding agent (or a
human reviewer) is held to on this codebase.

## Why this project is more than a comparison tool

Most of the work here isn't matching sequences — that part is a pairwise
aligner. The interesting part is what the platform refuses to do:

- **It withholds its own headline number.** SPEC §26 defines a weighted
  composite score for therapeutic precedent. Two of its six terms don't exist
  yet — one needs empirical calibration, the other needs a user-facing input
  the system doesn't collect. Rather than approximate them, the composite
  stays `None` with a machine-readable reason. [Full writeup](reports/GATE7_THERAPEUTIC_PRECEDENT.md).
- **It found and fixed a fabricated evidence flag.** A missing-data sentinel
  (`na` vs. the literal string `"None"`) meant every therapeutic in the
  corpus — 1,131 of them — was recorded as having an experimental structure
  behind it. The true split was 337 with structures, 794 without. Later, the
  fix's own blank-cell case was found to invert the same bug. The eventual
  solution validates data by *shape* (does this look like a PDB reference?)
  rather than by enumerating known-bad sentinels one discovery at a time.
- **A naming decision was killed by its own measurements.** The human
  repertoire coverage metric was headed for the name "humanness" until a
  corpus-wide run showed camelid VHH domains score a *higher* median than
  human VH. The metric is real; the name it almost got would have been a
  false claim baked into the schema.
- **Pairs are never synthesized.** A candidate's heavy and light chains are
  searched independently, but a "pair" is only reported when one source
  record links both chains together — never by combining the best heavy hit
  with the best light hit from different antibodies. The regression test for
  this rule initially asserted the wrong thing and caught its own author.

## Architecture

```
src/para_scope/
  models.py                 Component/evidence dataclasses, EvidenceStatus, enums
  numbering.py               ANARCII wrapper, IMGT coordinates, region extraction
  corpus.py                  In-memory corpus, ingestion, CORPUS_SCHEMA_VERSION
  persistence.py              Snapshot save/load, CACHE_FORMAT_VERSION
  features.py                 Deterministic sequence features
  motifs.py                   Motif observations in IMGT coordinates
  window_coverage.py          Human repertoire 9-mer window coverage
  therapeutic_precedent.py    Therapeutic precedent evidence (SPEC §26)
  pairing_precedent.py        Observed heavy-light pairing (SPEC §35.1)
  indexes/                    exact, positional, pair, approximate, shortlist
  parsers/                    sabdab2, thera_sabdab
```

**Python 3.13** (3.11+ required), `src/` layout. Dependencies are
deliberately few: [`anarcii`](https://pypi.org/project/anarcii/) for antibody
numbering, [`biopython`](https://biopython.org/) for alignment and residue
chemistry, [`promb`](https://pypi.org/project/promb/) for a bundled human
repertoire peptide reference set (MIT-licensed, ships inside the wheel — no
network call at scoring time). No structure prediction, no external service
calls in the scoring path, ever.

## Design principles

A short version of the rules in [CLAUDE.md](CLAUDE.md) that shape every line
of this codebase:

- **Only compare like with like.** Heavy chain against heavy chain, CDR-H3
  against CDR-H3 — never across chain types, and related-class comparisons
  (VH↔VHH, kappa↔lambda) are explicitly labeled and down-ranked, never merged
  into same-class evidence.
- **Missing evidence is not zero evidence.** A failed lookup, an inapplicable
  component, or unavailable data returns an explicit status
  (`not_found` / `not_applicable` / `insufficient_evidence` / …), never a
  default of 0.0.
- **Never silently modify input.** No residue substitution, deletion, or
  imputation — ever.
- **Every score is reproducible** from stored raw evidence plus a named
  formula version and calibration version.
- **Deduplicate before counting.** One antibody solved a dozen times in a
  structure database counts once, not twelve times (SPEC §15) — enforced
  end-to-end, including in scripts that measure corpus statistics.
- **Scoring code never queries a live external service.** All evidence comes
  from versioned, checksummed local snapshots.

## Quickstart

```bash
git clone https://github.com/rishiramesh17/Para-Scope.git
cd Para-Scope
python3.13 -m venv .venv
.venv/bin/pip install -e ".[dev]"
```

### Data (not bundled)

The corpus is built from two real antibody databases. `data/` is
`.gitignore`d — nothing here is redistributed, so you'll need your own copies:

1. **SAbDab2** (structural antibody database). Look up Zenodo record
   `20083995` and place the split CSV at
   `data/snapshots/sabdab2/splits_final/ab_split.csv`.
2. **Thera-SAbDab** (therapeutic antibody database, OPIG, Oxford). Download
   the 2025-03-18 (or later) `TheraSAbDab_SeqStruc_OnlineDownload.csv` and
   place it at `data/snapshots/thera_sabdab/TheraSAbDab_SeqStruc_OnlineDownload.csv`.

Both sources currently carry an unspecified license in their public
distributions — this blocks redistributing *the data itself*, not this
repository's code, and not local analysis. See
[docs/ROADMAP.md](docs/ROADMAP.md) for the exact status.

### Build the corpus

```bash
.venv/bin/python scripts/build_v2_corpus_cache.py        # SAbDab2 only (~8.5 min cold, ANARCII numbering)
.venv/bin/python scripts/build_combined_corpus_cache.py   # + Thera-SAbDab (~1 min)
```

Both write versioned, checksummed snapshots to `data/cache/` — regenerable
at any time, never committed. Subsequent rebuilds against the same source
files reuse the numbering cache and take seconds, not minutes.

### Run the tests

```bash
.venv/bin/python -m pytest -q             # unit suite, ~30s, no corpus required
.venv/bin/python -m pytest                # full suite incl. integration, needs the built corpus
```

### Try it

```bash
.venv/bin/python scripts/measure_therapeutic_precedent.py
.venv/bin/python scripts/measure_pairing_precedent.py
.venv/bin/python scripts/measure_window_coverage.py
```

Each prints a worked example against the real corpus — matches, raw
measurements, and exactly which evidence dimensions are blocked and why.

## Beyond Gate 11: a calibrated prediction layer (planned, not built)

Everything above is retrieval and comparison. It looks up what's already
known and reports it with provenance — it predicts nothing, by design,
through Gate 11.

The next phase adds a real prediction model on top, scoped to what
sequence-based ML can actually do today. Full binding-affinity, folding, or
clinical-outcome prediction from sequence alone remains an open, largely
unsolved problem across the field — no credible platform claims to have
solved it, and this one won't pretend to either.

What's genuinely tractable, and where the planned model is aimed:

- **Developability risk flags**, not binding prediction — aggregation
  propensity, viscosity risk, polyreactivity, the kind of in-silico
  screens (charge/hydrophobicity patch metrics in the style of TAP-type
  tools) already used industry-wide as narrow, measurable, evaluable
  endpoints. A model with a real accuracy number on a well-defined task,
  not a black-box oracle on an unsolved one.
- **The evidence this platform already computes becomes the model's
  feature set**, rather than re-deriving features from raw sequence a
  second time. Database-Relative Novelty, observed pairing precedent,
  therapeutic precedent, human repertoire window coverage — Gates 6
  through 8 exist specifically to produce calibrated, versioned inputs a
  downstream model can train on.
- Every prediction, when one ships, carries the same discipline as the
  evidence dimensions above it: a named model version, a confidence
  interval measured on held-out data, and an explicit statement of what
  it was and wasn't trained to predict. No confidence interval, no ship.
- **A plain-language layer on top, not a new source of claims.** An LLM
  can make a dense evidence profile genuinely readable — explaining what
  `qcov`, a withheld composite, or an `insufficient_evidence` status
  actually mean for a specific candidate — grounded entirely in numbers
  the platform already computed and can cite. It explains evidence; it
  does not generate new ones.

This section is intent, not a gate — no committed timeline or endpoint
ordering yet. Recorded here so it isn't miscast as a black-box AI bolted
onto a comparison tool: any predictive layer here inherits the same
evidentiary rigor as everything above it, or it doesn't ship.

## License

Code is MIT-licensed — see [LICENSE](LICENSE). This does not extend to the
third-party antibody databases the code can ingest; none of that data is
included in this repository.
