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

## Prohibited claims

Never state or imply that a candidate will bind, has a target, will fold, express, pair, remain stable, avoid aggregation, be soluble, be manufacturable, be safe, be effective, have low immunogenicity, or succeed experimentally or clinically. Similarity to an approved therapeutic transfers none of its properties.

Any user-facing summary score must carry this text verbatim:

> This index summarizes comparative database evidence. It is not a probability of folding, expression, pairing, binding, stability, manufacturability, safety, efficacy, or experimental success.

## Architecture

No code yet. Stack and module layout are still being decided.

## Development commands

No build system yet.
