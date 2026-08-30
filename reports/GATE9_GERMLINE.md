# Gate 9, milestone 1 — human germline V-gene assignment

Built 2026-08-29. Resolves the germline reference-data question `docs/ROADMAP.md`
flagged as unresolved and "verify before implementing." Full detail here; the
roadmap entry is now a short pointer to this file.

## What was verified before building anything

**License.** OGRDB (https://ogrdb.airr-community.org) publishes its germline
sets under CC0 (public domain) — confirmed directly against
`https://ogrdb.airr-community.org/render_page/licensing_statement.html`. Unlike
SAbDab2 and Thera-SAbDab (both `license_id = "unspecified"`), CC0 permits
redistribution, so the reference set is bundled inside the package
(`src/para_scope/reference_data/germline/ogrdb_human_v_genes.fasta`, 353 V-gene
alleles across IGHV/IGKV/IGLV) rather than gitignored and fetched separately —
the same pattern `promb` uses for its own reference data.

**Reduced scope**, matching the roadmap's own recommendation: human only (no
camelid/VHH reference set exists — VHH is `UNSUPPORTED_COMPONENT`, not scored
low), V gene only (no D/J call), one reference set (SPEC §11.2's
never-merge-two-germline-sets concern does not arise with only one).

## What was built

`src/para_scope/germline.py` (`GERMLINE_VERSION = "1.0.0"`,
`CALIBRATION_VERSION = "uncalibrated"`), mirroring the shape of
`therapeutic_precedent.py` and `pairing_precedent.py`: `compute_germline_assignment`
never raises, returns an explicit `EvidenceStatus` for every failure mode, and
never defaults an unset identity to 0.0. No new alignment code — germline V
genes are indexed into the same `ApproximateIndex` (SPEC §16.2) Gates 7 and 8
already use, typed by locus (`ComponentType.VH` / `VL_KAPPA` / `VL_LAMBDA`) so
comparison stays like-with-like automatically.

`scripts/build_germline_reference.py` downloads OGRDB's published Homo sapiens
IGH_VDJ / IGKappa_VJ / IGLambda_VJ sets, keeps V genes only, translates each in
frame 1 (IMGT V-REGION convention), and drops pseudogenes/ambiguous calls (14
of 367 raw V-gene entries).

## Spot check: Trastuzumab (pinned in `tests/integration/test_germline_snapshot.py`)

| Chain | Nearest germline | Identity | Query coverage | Subject coverage |
|---|---|---:|---:|---:|
| VH | IGHV3-66 | 0.667 | 0.817 | 1.000 |
| VL-kappa | IGKV1-39 | 0.766 | — | — |

Sanity-checked against public expectation before trusting the corpus-wide run:
trastuzumab's published germline is IGHV3-66/IGKV1-39 in the literature, which
this module reproduces from alignment alone.

## Corpus-wide sample (not exhaustive — see note below)

Exhaustive scoring (every unique VH/VL against all 353 reference alleles) was
too slow for a first pass — the run was stopped mid-flight and re-run against a
deterministic 1-in-8 sample per type instead. This characterizes the shape of
the distribution; it is not a final percentile table. Gate 9's actual
calibration work should re-run exhaustively, or with a validated sampling
scheme, once statistical calibration is in scope.

| Type | Sampled n (of unique) | AVAILABLE | Identity p05 / p50 / p95 |
|---|---:|---:|---|
| VH | 684 (of 5,472) | 28.8% | 0.517 / 0.718 / 0.797 |
| VL-kappa | 511 (of 4,089) | 32.3% | 0.576 / 0.764 / 0.877 |
| VL-lambda | 137 (of 1,100) | 30.7% | 0.554 / 0.819 / 0.889 |
| VHH | 145 (of 1,167) | 0% (100% UNSUPPORTED_COMPONENT) | — |

Most common nearest germline: VH → IGHV3-23, IGHV1-69, IGHV1-46; VL-kappa →
IGKV1-39, IGKV3-11, IGKV4-1 — all well-known, frequently used human germline
families, a soft sanity check that the alignment is finding real biology
rather than noise.

## Finding that changes Gate 9's calibration work: ~70% NOT_FOUND, and why

Roughly 70% of sampled VH/VL components fail to clear the 0.70 coverage floor
and return `NOT_FOUND` rather than a value. Traced directly (not guessed): a
meaningful share of `ComponentType.VH`-typed components in the corpus are
~210–240 residues — roughly double a typical single VH domain (trastuzumab's
is 120). Against a germline reference of ~98 residues (a V gene covers FR1
through FR3, not CDR3/FR4), `subject_coverage` on these still reaches ~0.98–1.0
(the germline gene fully aligns somewhere inside the longer query), but
`query_coverage` caps near 0.44 — well under the 0.70 floor — no matter how
good the local match is. These are dropped as inadmissible, correctly by the
letter of the existing rule, but the rule (`MIN_QUERY_COVERAGE = 0.70`,
`indexes/approximate.py`) was calibrated for whole-domain-vs-whole-domain
matching, not partial-gene-vs-domain matching, and was never re-examined for
this use case.

This is not fixed here. Loosening the floor without calibration would be
exactly the kind of guess CLAUDE.md and the existing "0.70 coverage floor ...
not a calibrated value" note (docs/ROADMAP.md, "Known-uncalibrated items
carried into Gate 9") already warn against. Recorded as a carried-forward item
for Gate 9 proper: germline admissibility likely needs its own floor,
calibrated against `query_coverage` capped by germline-gene-length /
component-length rather than the whole-domain floor borrowed unchanged.

## What this unblocks

- Database-Relative Novelty (Gate 9) now has a baseline to measure drift
  against, where previously `EvidenceClass.GERMLINE` was declared and
  populated by nothing.
- Gate 8's `vh_vl_germline_pair_precedent` blocked dimension (SPEC §35) can be
  computed once Gate 9 wires this module in — germline calls now exist for
  both chains of an observed pair.

## What remains out of scope, deliberately

- Camelid/VHH germline reference.
- D and J gene calls, full V(D)J assignment.
- A second reference set (IMGT/GENE-DB) and the SPEC §11.2 no-silent-merge
  reporting that would require.
- The query-coverage-floor recalibration described above.
- Turning `nearest_germline_identity` into a percentile — that is Gate 9
  proper, and needs the reference-population-bias caveat already on record
  (SAbDab2 is structure-derived, not repertoire-representative).
