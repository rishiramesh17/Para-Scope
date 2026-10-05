# Gate 9, milestone 2 — source chains split at the variable domain

Date: 2026-10-04
Corpus schema: 5.0.0 → **6.0.0** (component identity changes; snapshots are not comparable)
Tests: 316 passing

## The defect

Variable-domain components were storing the **whole source chain**, not the
numbered variable domain. Measured on the real combined corpus before the fix:

| Type | n | median stored length | median numbered V-domain |
|---|---|---|---|
| VH | 5,473 | **222** | 121 |
| VL-kappa | 4,089 | **214** | ~107 |
| VL-lambda | 1,100 | **214** | ~110 |

68.4% of `VH` components exceeded 160 residues. A V-domain is ~120. The
remainder was constant-region sequence, confirmed by inspection — this is the
C-terminal end of a component typed `VH`:

```
...GYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPR
```

That is CH1. The components were Fab chains.

### Why it mattered

This violated a non-negotiable constraint in `CLAUDE.md`:

> A variable-domain match is never a full-chain match. Keep variable-only and
> full-chain records distinguishable at every layer.

Constant domains are near-identical between unrelated antibodies, so including
them inflated every comparison. Measured, identity to the nearest *other*
antibody:

| Compared as | median identity |
|---|---|
| stored sequence (V + CH1) | 0.858 |
| V-domain only | 0.803 |
| **inflation** | **+0.056** |

The bias is systematic and directional: sequences looked **less novel than they
are**. Calibrating Gate 9's reference distributions on top of it would have
baked that into every percentile, version-stamped and authoritative-looking.

## The fix

`numbering.split_chain_at_variable_domain` splits each source chain into
`(prefix, variable_domain, suffix)`. The V-domain is located by substring search
(ANARCII reports positions, not input offsets, and insertion codes break offset
arithmetic) — verified contiguous in **11,830 / 11,830** variable-domain
components.

- **V-domain** → `VH` / `VL-kappa` / `VL-lambda` / `VHH`
- **C-terminal remainder** → new `CONSTANT_HEAVY` / `CONSTANT_LIGHT`, indexed
  when ≥ `MIN_CONSTANT_REGION_LENGTH` (30, provisional/uncalibrated)
- **N-terminal remainder** → counted, deliberately **not** typed. A segment
  before a V-domain may be a tag, signal peptide, or fusion partner; no
  component type can honestly claim it. Residues stay recoverable via
  `source_pointers`.

Queries are normalized the same way (`normalize_query_to_variable_domain`), or
the comparison would be chain-against-domain. Every normalization emits a note
onto the result's `limitations` — the input is not what was compared, and that
must be visible.

`PairIndex` keys moved to V-domain checksums via `_checksum_of_record_id`.
Registering pairs under chain checksums while lookups used domain checksums
would have silently returned zero observed pairs for every query.

## Results

| Metric | before | after | reading |
|---|---|---|---|
| VH components | 5,473 | 4,745 | 728 duplicate V-domains merged |
| VL-kappa | 4,089 | 3,540 | 549 merged |
| VL-lambda | 1,100 | 975 | 125 merged |
| VHH | 1,168 | 1,027 | 141 merged |
| `CONSTANT_HEAVY` | — | **1,091** | newly available |
| `CONSTANT_LIGHT` | — | **552** | newly available |
| VH median length | 222 | **121** | a V-domain |
| Unique pairs | 5,686 | 4,913 | same 12,032 records, better dedup |
| Pairs in **both** sources | **2** | **207** | structure ↔ therapeutic linkage |
| Germline `AVAILABLE` | ~30% | **100%** | see below |

Backing records stayed at 12,032 — no source data was lost; the same records
now collapse onto fewer deduplicated pairs (2.12 → 2.45 records per pair).

### Two roadmap conclusions were wrong

1. **"Neither SAbDab2 nor Thera-SAbDab carries constant-region sequence"** —
   false. SAbDab2 carries it; it was inside the variable-domain components. 1,643
   constant-region components now exist. `FULL_HEAVY`/`FULL_LIGHT` remain empty
   (a Fab chain is V+CH1 and carries no CH2/CH3, so claiming a full chain would
   still overstate), but constant-region evidence is no longer absent.

2. **"Germline needs its own separately calibrated coverage floor"** — not a
   calibration problem. ~70% `NOT_FOUND` was the length mismatch between
   210–240-residue "VH" components and ~98-residue germline genes. With real
   V-domains, a 180-component sample returns **100% `AVAILABLE`**, median
   nearest-germline identity **0.726** (range 0.463–0.888) — plausible for
   somatically mutated antibodies. The 0.70 floor was never the problem.

3. **"The sources barely overlap"** — an artifact. The same antibody present in
   both sources produced two different pair keys because its constant regions
   were attached, so it never matched itself across sources. Trastuzumab's exact
   pair now resolves to one deduplicated pair carrying 16 PDB structures *and* 11
   therapeutic programmes.

## Carried forward

- `MIN_CONSTANT_REGION_LENGTH = 30` is provisional, like the 0.70 coverage floor
  and `MINIMUM_WINDOWS_FOR_FRACTION`. 2,927 remainders fell below it.
- 4,152 chains carry an N-terminal flank that no component claims. Deciding
  whether any of those are typeable is separate work.
- Constant regions are **undecomposed**. Splitting into CH1 / Hinge / CH2 / CH3 /
  CL needs a constant-region reference (SPEC §6.5) that is not bundled. Length
  alone is not evidence of identity.
- Gate 7 and Gate 8 *code* is unaffected, but every measured number in their
  reports was computed on inflated identities and has not been re-measured here.
