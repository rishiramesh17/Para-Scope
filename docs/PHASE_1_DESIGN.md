# Phase 1 Design: Search Architecture and Source Ingestion

## Overview

Phase 1 implements search and indexing (SPEC §16) and source-snapshot ingestion (SPEC §13–14). The goal: take a single antibody component (e.g., VH sequence) and retrieve all matching components from versioned, deduplicated local corpora, producing raw evidence records (SPEC §18) with no scoring layers yet.

**Scope:**
- Exact-sequence indexing for all component types
- Approximate sequence search for long components (pure-Python pairwise alignment, no external binary)
- Numbered-position indexing for framework and CDR comparison (>95% position overlap threshold)
- Pair indexing for observed heavy–light chains
- Ingest and normalize two corpora: SAbDab2, Thera-SAbDab
- Produce immutable evidence records with full provenance

**Deferred to Phase 2: paired-OAS.** OAS's paired dataset is ~1.8M sequences, likely 100s of GB — well beyond what's sensible to store speculatively. SAbDab2 (~21K antibody instances) and Thera-SAbDab (~1K therapeutics) are both under 2GB combined and sufficient to validate the entire search/dedup/indexing architecture end-to-end. OAS gets added once Phase 1 proves the pipeline and real query patterns reveal how much repertoire data is actually needed — at that point we'll decide between an external drive or cloud object storage (R2/B2) based on the actual volume required, not a guess.

**Out of scope (Phase 2+):**
- Short-CDR empirical null models and FDR correction (SPEC §21)
- Calibrated percentiles (SPEC §19 Layer 2)
- Evidence grades and sufficiency (SPEC §19 Layer 3)
- Score calculations
- Paratope residue-set comparison (SPEC §34)

---

## Deduplication Strategy

From SPEC §15, deduplicate at several levels. Phase 1 implements:

### Exact Sequence (mandatory)
Hash each component sequence (SHA256). Identical sequences from different sources become one canonical record carrying all evidence pointers.

**Implementation:** Dictionary keyed by checksum; records are deduplicated on ingestion.

### Parent Antibody (mandatory for structures)
Multiple PDB structures of one antibody (e.g., trastuzumab 1FVH, 1FVG, 1N8Z) should count as one parent antibody with multiple structure observations, not three independent antibodies.

**Implementation:**
- SAbDab2 provides `sabdab_id` (unique variable region identifier), and each PDB chain links to it.
- Store: `(sabdab_id, pdb_id)` pairs separately. Dedup before counting.
- When returning evidence, report: "Trastuzumab VH (1 unique antibody, 3 PDB structures)".

### Study and Subject (deferred to Phase 2)
Repertoire data (OAS) includes study metadata. Repeated records from one subject or one clonal expansion inflate counts. Full dedup requires study-level reasoning and clonal-clustering algorithms.

**Phase 1 approach:** Flag high-count sequences in paired-OAS with metadata, but don't attempt automatic study/subject dedup. Manual review or Phase 2 heuristics.

### Therapeutic Aliases (manual, out-of-band)
Different names/codes for the same therapeutic (e.g., Herceptin = trastuzumab). Thera-SAbDab is small enough to curate by hand; Phase 1 builds the canonical mapping as a simple JSON lookup.

---

## Search Architecture

Implement four index types per SPEC §16, each optimized for a component pattern:

### 1. Exact-Sequence Index
**For:** All component types (VH, VL, CDRs, frameworks, constant regions).

**Implementation:**
- Key: SHA256(sequence)
- Value: [evidence_record_ids]
- Lookup: O(1), 100% precision.
- Source: Populate on ingestion from each corpus.

**Query flow:**
```
Input sequence → hash → exact-index lookup → hit or miss
```

**Data structure:** `dict[str, list[str]]` in memory; can be serialized to JSON for persistence.

### 2. Approximate Long-Sequence Index
**For:** VH, VL, constant regions, frameworks >60 residues.

**Implementation:**
- Use `Bio.Align.PairwiseAligner` (local mode) — `Bio.pairwise2` is deprecated in current Biopython (1.87) and slated for removal, so the aligner class is the supported pure-Python path.
- Store: Sequences + metadata as searchable records; no precomputed index (too large for laptop). Compute alignment on query.
- Top-k retrieval: Rank by sequence identity + coverage; return top 10 hits.
- Alignment details: substitution_matrix=BLOSUM62, gap_open=-10, gap_extend=-1.
- Eligibility enforced in code: only VH/VL/VHH/full-chain/framework components may enter this index (`indexes/approximate.py::ELIGIBLE_TYPES`); CDRs are rejected at insertion time, not just by convention.

**Query flow:**
```
Input sequence → align against all same-type sequences → top-k by identity + coverage → return alignments
```

**Rationale:** Pure Python, no external binaries (NCBI BLAST not guaranteed installed). Laptop-scale dataset (~10K sequences) is fast enough (<1s per query).

### 3. Numbered-Position Index
**For:** Framework and CDR position-by-position comparison (SPEC §16.3).

**Implementation:**
- Per component type (CDR-H3, FR-H1, etc.), build inverted index: IMGT position → [records at that position].
- Store: NumberedComponent tuples (from Phase 0's numbering.py) keyed by checksum.
- Lookup: Given an input component, find all records that share >80% of positions (or configurable threshold).

**Query flow:**
```
Input sequence → number under IMGT → extract positions → position-index lookup → find overlapping records
```

**Data structure:** `dict[ComponentType, dict[int, list[str]]]` (component type → position → record IDs).

### 4. Pair Index
**For:** Observed heavy–light chain pairs (SPEC §16.5).

**Implementation:**
- Store: `(heavy_checksum, light_checksum) → pair_metadata`
- Metadata: source (SAbDab2, Thera-SAbDab; OAS in Phase 2), structure availability, study.
- Constraint: Only directly-observed pairs from source records, never synthesized.

**Query flow:**
```
Input (VH, VL) → (checksum_h, checksum_l) → pair-index → hit or "no observed pairing"
```

**Data structure:** `dict[tuple[str, str], dict[str, Any]]`.

---

## Corpus Ingestion Order and Pipeline

### Phase 1a: SAbDab2 (weeks 1–2)
**Why first:** Smallest, highest quality, already deduplicated by OPIG. Parser validation is easiest here.

**Data:**
- ~11,085 PDB IDs
- ~21,237 distinct antibody instances
- ~6,540 unique variable regions (sabdab_id)
- Includes paired VH/VL from structures

**Ingestion pipeline:**
```
SAbDab2 download
  ↓
Parse: extract heavy chain, light chain, PDB metadata
  ↓
Validate: sequence, chain type (H/L/K)
  ↓
Number under IMGT using ANARCII (Phase 0 tool)
  ↓
Normalize: populate ComponentType, source_snapshot, checksum
  ↓
Deduplicate: exact sequence + parent antibody
  ↓
Index: exact-hash, approximate, positional, pair
  ↓
Snapshot metadata: version, checksum, record count, parser version
```

**Acceptance:**
- 0 parse errors on all 11K+ PDB records
- Numbering success rate >99% (some malformed structures expected)
- Exact-index has 0 hash collisions (SHA256, negligible probability)
- Deduplicated count is ~5K unique sequences (order-of-magnitude estimate)

### Phase 1b: Thera-SAbDab (weeks 2–3)
**Why second:** Small (~1K therapeutics), validates therapeutic-precedent evidence, allows spot-checking.

**Data:**
- ~1,000 therapeutic antibodies
- Known names, targets, formats
- Variable-region sequences (curated)
- Partial pairing (literature-reported, not all paired)

**Ingestion:** Same pipeline, plus:
- Build therapeutic-alias canonical map (manual JSON, if ambiguities found)
- Validate: name, target, format match literature

**Acceptance:**
- All ~1K records parse and number
- No unexpected duplicates (each therapeutic should be unique)
- Aliases resolved consistently

### Phase 2: Paired-OAS (deferred)
**Why deferred:** OAS's full paired dataset is ~1.8M sequences, likely hundreds of GB — not something to store speculatively before the pipeline is proven. SAbDab2 + Thera-SAbDab (under 2GB combined) are sufficient to validate parsing, dedup, indexing, and search end-to-end.

**Decision point for Phase 2:** Once Phase 1 is running, real query patterns will show how much repertoire breadth is actually needed. That determines the OAS subset size, which then determines storage: local disk (if small enough), an external drive (~$60–100 one-time for 1–2TB), or cloud object storage (Cloudflare R2 / Backblaze B2, ~$5–15/month, no egress fees) if the subset needs to stay off local disk.

---

## Example Data Flow: VH Query

**Input:** Trastuzumab heavy-chain variable domain

**Step 1: Exact search**
```
Input: EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQ...
Hash: 8b8998709aa9423e89f2c0b71328bf9702d979ec...
Lookup exact-index[hash] → [SAbDab2::1FVH_H, SAbDab2::1N8Z_H, ...]
Create evidence records for each hit
  component_relationship: EXACT_SAME_TYPE
  sequence_identity: 1.0
  query_coverage: 1.0
  subject_coverage: 1.0
```

**Step 2: Approximate search (if needed)**
```
Align input against all VH in index
Top hits by (identity + sqrt(coverage_harmonic_mean)):
  - VH from pembrolizumab: 92% identity, 98% coverage → score 0.90
  - VH from torcetrapib: 87% identity, 95% coverage → score 0.82
  ...
Create evidence records for top 5
  component_relationship: EXACT_SAME_TYPE
  alignment_method: "smith_waterman"
  sequence_identity: 0.92
  query_coverage: 0.98
```

**Step 3: Positional search**
```
Number input VH under IMGT → positions [1–128]
Query position-index[VH] → find all VH with >95% position overlap
Results: e.g., VH with same CDRs, different framework residues
Create evidence records
  component_relationship: EXACT_SAME_TYPE
  alignment_method: "numbered_position"
  canonical_position_agreement: 0.87
```

**Step 4: Output**
```
Emit list of EvidenceRecord (from CLAUDE.md / models.py)
Each record:
  - Raw measurements (identity, coverage, alignment_length, etc.)
  - Source provenance (source_database, snapshot_id, parser_version)
  - Component relationship (validates against compare_eligibility)
  - No scores (Layer 1 raw data only)
```

---

## Known Limitations and Phase 2 Deferral

**Not in Phase 1:**
- Paired-OAS ingestion: deferred until storage/scale decision is made (see Corpus Ingestion Order)
- Short-CDR empirical null models (SPEC §21): Will be built in Phase 2 using calibration dataset
- Evidence grades (A/B/C/D): Deferred to Phase 2, requires validation benchmark
- Study/subject deduplication for repertoire (SPEC §15.4): Simple heuristics in Phase 2, applies once OAS is added
- Calibrated percentiles (SPEC §19 Layer 2): Phase 2, requires reference distributions
- Paratope residue-set comparison (SPEC §34): Phase 2, requires curated paratope corpus

**Known risks:**
- Pair index doesn't validate whether observed pairs are actually compatible (VH–VL interface); that's Phase 2
- No handling of CDR-length variants in short-CDR; Phase 1 exact/approximate only
- Repertoire-precedent evidence class (SPEC §17.3) will show `insufficient_evidence` for all components until OAS is ingested — this is correct behavior (SPEC: missing evidence ≠ zero evidence), not a bug

---

## Deliverables

By end of Phase 1:
1. `src/antigens/parsers/` — parsers for SAbDab2, Thera-SAbDab (versioned, with checksums)
2. `src/antigens/indexes/` — exact, approximate, positional, pair index implementations
3. `src/antigens/corpus.py` — ingest pipeline and deduplication
4. `data/snapshots/` — downloaded and checksummed source data + metadata
5. `tests/test_search.py` — search correctness (e.g., trastuzumab VH finds exact matches)
6. `docs/PHASE_1_RESULTS.md` — corpus stats, dedup results, search latency benchmarks

---

## Timeline

- **Week 1:** SAbDab2 parser + exact/approximate/positional indexes
- **Week 2:** Deduplicate, validate, benchmark on SAbDab2
- **Week 3:** Thera-SAbDab parser + pair indexing
- **Week 4:** Integration tests, `PHASE_1_RESULTS.md`, decide OAS storage approach for Phase 2

Total: ~4 weeks @ 4–5 hrs/day.
