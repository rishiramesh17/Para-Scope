# Gate 6 — Human Repertoire Window Coverage

**Date:** 2026-07-30
**Engine version:** `WINDOW_COVERAGE_VERSION = 1.0.0`
**Module:** `src/antigens/window_coverage.py`
**Tests:** 30 new, suite total 150 → 180, all passing
**Corpus under test:** combined snapshot, 46,627 components, `corpus_schema_version` 4.0.0

---

## What was built

A 9-residue sliding-window matcher. For a component sequence, it chops every
overlapping 9-mer and counts how many appear verbatim in a fixed reference set of
peptides observed in human antibody repertoires. Coverage is that count over the
window total. Nothing is modelled, predicted, or inferred.

Throughput: **46,627 components in 3.0 s (15,784/s)** after a one-time 3.4 s load
of the reference set. Fast enough that no caching layer is warranted.

## Dependency resolution

The roadmap specified `promb` but flagged it as unverified. Verified:

| Property | Finding |
|---|---|
| Exists on PyPI | Yes — 1.0.2 |
| License | **MIT** — cleanly redistributable, unlike Thera-SAbDab's `unspecified` |
| Reference data | Bundled in the wheel (39.9 MB), not downloaded |
| Network calls | **None** in any module — grep of `db.py`, `cli.py`, `netmhc.py`, `parsing.py` |
| Python support | ≥3.8, runs on the project's 3.13 without downgrading anything |
| Transitive deps | scipy, pandas, tqdm |

The bundled set is `OASis_9mers_v1_10perc_subjects.txt.gz` — 5,725,303 peptides of
length 9 observed in ≥10% of subjects in the Observed Antibody Space collection.

**Not recorded in the roadmap, and worth recording:** `promb` is authored by David
Prihoda, who also authored BioPhi, and its `oasis` subcommand reproduces BioPhi's
OASis identity score. This is not a substitute metric — it is the same published
metric (Prihoda et al., *mAbs* 2022) repackaged for a modern Python stack. The
rejection of BioPhi-the-distribution stands on stack-collision grounds; the
science underneath is unchanged, which strengthens rather than weakens the choice.

Only the in-memory exact-match path is used. The scipy-backed nearest-peptide
search is not called, so scipy is carried but not exercised.

## Validation — the roadmap's figures reproduce

| Sequence | Recorded | Measured |
|---|---|---|
| trastuzumab VH (humanized) | 0.71 | **0.7054** |
| trastuzumab VL (humanized) | 0.83 | **0.8283** |
| muromonab VH (mouse) | 0.58 | **0.5225** |
| caplacizumab (VHH, llama) | 0.685 | **0.6333** |

The humanized figures reproduce exactly. The mouse and VHH figures differ from the
roadmap's because the roadmap did not record *which* mouse and VHH sequences were
used; these are muromonab and caplacizumab. The ordering — and therefore the
conclusion drawn from it — is unchanged.

These four values are pinned in `tests/test_window_coverage.py` so the module
docstring's justification cannot silently go stale.

---

## Finding 1 — the naming constraint is stronger than the roadmap recorded

The roadmap flagged the llama result as "the load-bearing limitation" on the basis
of four hand-picked sequences. The corpus-wide distribution makes the point far
more forcefully:

| Component | n | p05 | p25 | **p50** | p75 | p95 |
|---|---:|---:|---:|---:|---:|---:|
| VH | 5,473 | 0.197 | 0.333 | **0.429** | 0.655 | 0.904 |
| VHH | 1,168 | 0.351 | 0.478 | **0.547** | 0.614 | 0.739 |
| VL-kappa | 4,089 | 0.205 | 0.341 | **0.445** | 0.695 | 0.996 |
| VL-lambda | 1,100 | 0.188 | 0.396 | **0.464** | 0.704 | 0.991 |

**Camelid VHH have a higher median window coverage than VH — 0.547 against 0.429.**
Not comparable; higher. A metric presented as "humanness" would rank llama single
domains above the median antibody variable domain in a structure-derived corpus.

The VHH distribution is also markedly tighter (p05–p95 spans 0.39, against 0.71
for VH), so the overlap is not an artifact of a few outliers.

This is decisive for naming. The module is `window_coverage.py`, the field is
`human_window_coverage`, and the term is **Human Repertoire Window Coverage**,
which states the comparison performed rather than a conclusion drawn from it. The
convention matches `Observed Repertoire Precedent` and `Database-Relative
Novelty`. A test asserts that no public name in the module contains "humanness",
"naturalness", or "immunogenic", so the decision is enforced rather than merely
documented.

## Finding 2 — the metric is uninformative for CDR-H3, and near-saturated for frameworks

| Component | Available | Median coverage | Interpretation |
|---|---:|---:|---|
| CDR-H3 | 28 of 5,180 | **0.000** | Every scoreable CDR-H3 scored exactly zero |
| FR-H3 | 4,304 of 4,304 | 0.900 | p75 and p95 both 1.000 |
| FR-L3 | 2,645 of 2,645 | 1.000 | p50, p75, p95 all 1.000 |

CDR-H3 is junctionally generated and is the diversity hotspot of the repertoire,
so its 9-mers are effectively unique and appear in no reference set built from
peptides shared across ≥10% of subjects. **Coverage 0.000 for CDR-H3 is expected
biology, not a defect** — but it means the value carries no information there, and
a low CDR-H3 coverage must never be read as a finding about the candidate. Only 28
CDR-H3s were even long enough to clear the minimum-window rule.

Frameworks sit at the opposite extreme: FR-L3's median, p75 and p95 are all 1.000,
so the distribution is saturated and cannot discriminate between components. Gate 9
must not build a percentile scale on a saturated distribution.

Neither problem is fixable by tuning. Both are properties of what the metric
measures, and both belong in Gate 9's calibration as exclusions.

## Finding 3 — only 40.3% of the corpus can be scored at all

| Status | Count | Share |
|---|---:|---:|
| `available` | 18,808 | 40.3% |
| `insufficient_evidence` | 19,218 | 41.2% |
| `not_applicable` | 8,601 | 18.4% |

Most CDRs are shorter than the 9-residue window and are correctly
`not_applicable`; CDR-L2 (IMGT, ~3 residues) is 435 of 435 not applicable. Nearly
all frameworks other than FR-H3/FR-L3 fall below the 20-window minimum.

This is the honest coverage story, and it is the reason status reporting is
mandatory rather than cosmetic: **the majority of components produce no value, and
none of them may be recorded as a zero.**

One anomaly not chased down: a single FR-H2 scored `available`, meaning one
FR-H2 in the corpus is at least 28 residues long. FR-H2 is normally ~17. This is
most likely a numbering artifact and is logged here rather than silently dropped.

---

## Two library behaviours deliberately not inherited

1. **`promb` raises `ZeroDivisionError` on any sequence shorter than the window.**
   Most CDRs are shorter than nine residues, so passing a CDR to the library
   crashes it. Handled as `NOT_APPLICABLE` with `coverage=None`.

2. **`promb` silently scores sequences containing non-standard residues.**
   `EVQLVESGGGXVQPGGSLRLSCAAS` returns 0.47 rather than an error. Windows spanning
   an `X` match nothing, so an unresolvable residue *lowers* coverage and reads as
   "less human" when it means "we could not evaluate." That is absence of evidence
   rendered as a low score, which `CLAUDE.md` forbids. Handled as
   `UNSUPPORTED_COMPONENT`, with the offending residues named and the input left
   unmodified.

The second is the more serious of the two, because it fails silently and in the
direction that produces a confident wrong answer.

## Supporting change — `EvidenceStatus`

Before this gate, the seven `CLAUDE.md` statuses existed only as ad-hoc string
prefixes inside parser skip reasons. `models.py` now carries `EvidenceStatus` with
the eight values SPEC §22 names — the seven plus `AVAILABLE`. `coverage` is
populated only for `AVAILABLE`; every other status leaves it `None`, never 0.0,
never a midpoint, never a corpus mean.

A test enumerates every failure path and asserts that none returns a value and all
carry a reason.

## Known-uncalibrated, carried into Gate 9

- **`MINIMUM_WINDOWS_FOR_FRACTION = 20` is a provisional default**, not a
  calibrated threshold — the same standing as the 0.70 coverage floor. It governs
  41.2% of the corpus, so it is a consequential guess.
- **The corpus distributions above are distributions of solved structures**, not of
  natural repertoire. Same reference-population bias already recorded for Layer 2
  percentiles. Defensible if labelled; indefensible if presented as repertoire
  percentiles.
- **CDR-H3 and the saturated frameworks need explicit exclusion rules**, not
  percentile scales.
- **Short components still lack a length-matched empirical null** (SPEC §21). The
  raw window counts are returned for every status precisely so Gate 9 can build one.

## Files

| File | Change |
|---|---|
| `src/antigens/window_coverage.py` | New — Gate 6 engine |
| `src/antigens/models.py` | Added `EvidenceStatus` |
| `scripts/measure_window_coverage.py` | New — corpus-wide distributions |
| `tests/test_window_coverage.py` | New — 30 tests |
| `pyproject.toml` | Added `promb>=1.0.2` |
