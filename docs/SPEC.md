# OPTIMIZED MASTER PLANNING PROMPT

## Explainable Antibody Component Comparison and Candidate-Assembly Platform

You are acting as a multidisciplinary product, scientific, data, legal-access, and software architecture team consisting of:

- A senior antibody engineer
- A computational immunologist
- An immunoglobulin sequence-analysis specialist
- A structural bioinformatician
- A protein and antibody database architect
- A scientific data engineer
- A biostatistician specializing in evidence calibration
- A machine-learning validation specialist
- A scientific software architect
- A privacy and research-data governance specialist
- A scientific user-experience designer
- A research-software quality and validation specialist

Your assignment is to produce a complete, implementation-ready plan for an evidence-based antibody component comparison and candidate-assembly platform.

Do not begin by generating a generic website, a black-box artificial-intelligence model, an antibody-binding predictor, a docking system, or an autonomous antibody-design engine.

First define and validate:

1. The scientific scope
2. The antibody component ontology
3. Candidate-assembly rules
4. Supported input contracts
5. Sequence and evidence-record validation
6. Canonical numbering and coordinate systems
7. Required databases and legal access methods
8. Data-ingestion and versioning rules
9. Same-component comparison methods
10. Search and indexing architecture
11. Evidence classes
12. Raw measurements
13. Score calibration
14. Evidence sufficiency and uncertainty
15. Per-component scoring
16. Whole-construct summary rules
17. Provenance and reproducibility
18. Privacy and confidential-sequence protections
19. User-interface requirements
20. MVP boundaries
21. Benchmarking and acceptance criteria
22. Risks, limitations, and excluded claims

---

# 1. Non-Negotiable Governing Goal

The governing product goal is:

> **The platform allows users to submit existing antibody sequences or assemble a candidate antibody construct from user-selected components. It validates the proposed component arrangement, compares every component with equivalent known components, and creates a transparent comparative evidence profile. The assembled construct is a computational representation and not proof of a functional antibody.**

The final product must be an:

> **Explainable antibody component comparison and candidate-assembly platform—not primarily an antibody–antigen binding predictor.**

The platform must answer questions such as:

- What known heavy chains resemble this heavy chain?
- What known light chains resemble this light chain?
- What known VH, VL, VHH, or other variable domains resemble the submitted domain?
- What known CDR-H3 sequences resemble this CDR-H3?
- What known framework regions resemble this framework?
- What known constant or Fc regions resemble this region?
- What experimentally reported paratope residue sets resemble the submitted numbered residue set?
- How common or unusual is each component within the available databases?
- Has the component been observed in immune-repertoire data?
- Has it appeared in an antibody with an experimentally determined structure?
- Has it appeared in a therapeutic antibody?
- Have similar heavy and light chains been directly observed together?
- Does the proposed domain order resemble known antibody formats?
- What evidence exists for every match?
- What is missing?
- How was every measurement, percentile, grade, and score calculated?

The platform must not answer its primary question as:

> “Will this antibody bind?”

It must not convert comparative precedent into a probability of biological success.

## First-Version Structural Boundary

The first version is sequence-first and metadata-aware. It may use structural databases to identify that a matched component appears in an experimentally determined structure and to retrieve associated metadata. It must not process coordinates, predict structures, calculate three-dimensional similarity, or render molecular models. Full structural computation is a future module, not an MVP requirement.

---

# 2. Required Product Language

Use scientifically careful terminology.

## 2.1 Candidate construct

When a user combines components, call the result a:

- Candidate antibody construct
- Candidate assembly
- Digital construct representation
- Proposed antibody architecture

Do not state that the platform has “created a functional antibody.”

## 2.2 Paratope

Use the correct term **paratope**.

A paratope is normally a set of residues distributed across the folded antibody surface. It is not necessarily a contiguous sequence and must not be treated as an ordinary CDR string.

## 2.3 Observed repertoire precedent

Use **Observed Repertoire Precedent** instead of “natural precedent.”

An immune-repertoire sequence record establishes that a sequence was observed in a particular dataset. It does not necessarily prove productive expression, correct pairing, folding, binding, or biological fitness.

## 2.4 Database-relative novelty

Use **Database-Relative Novelty** instead of absolute novelty.

Absence from the supported databases does not prove biological novelty.

## 2.5 Evidence index

Any whole-construct score must be called:

> **Whole-Construct Comparative Evidence Index**

It must never be called:

- Binding score
- Success probability
- Functional probability
- Safety score
- Efficacy score
- Clinical score
- Manufacturability probability

---

# 3. Explicit Scientific Boundaries

The platform may compare and summarize evidence. It must not claim that a submitted or assembled candidate:

- Will bind an antigen
- Has a specific antigen target
- Will preserve the target of a matched antibody
- Will neutralize a pathogen
- Will activate or inhibit a biological pathway
- Will fold correctly
- Will express successfully
- Will pair correctly
- Will retain the same CDR conformations
- Will be stable
- Will avoid aggregation
- Will be soluble
- Will be manufacturable
- Will have favorable pharmacokinetics
- Will have low immunogenicity
- Will be safe
- Will be effective
- Will succeed experimentally
- Will succeed clinically

The platform may report that a component resembles one found in a known therapeutic, structure, repertoire record, or experimentally characterized antibody. It must also show the important differences and state that properties do not automatically transfer between sequences.

The platform must distinguish:

- Direct experimental evidence
- Experimental evidence about a homolog
- Curated sequence evidence
- Observed repertoire evidence
- Structure-linked experimental metadata
- Therapeutic precedent
- Literature-derived annotation
- User-supplied assertions
- Computational predictions
- Inferences made by the platform
- Missing evidence

---

# 4. Two Primary Operating Modes

## Mode A: Analyze an existing antibody

The user submits one or more existing antibody sequences and may optionally provide known database accessions as metadata.

Supported examples:

- Paired full heavy and light chains
- Heavy chain only
- Light chain only
- VH only
- VL only
- VHH or nanobody domain
- Fab
- scFv
- Full-length immunoglobulin
- Constant or Fc region
- Individual CDR or framework region with parent context
- Optional PDB, SAbDab2, or other supported structure-record accession
- User-defined paratope residue set anchored to numbered antibody positions

The system identifies, validates, decomposes, and compares the submitted material.

## Mode B: Assemble a candidate construct

The user selects or enters components and defines how they are arranged.

Examples:

- VH plus VL plus selected constant regions
- Heavy variable domain plus heavy constant domains
- Light variable domain plus kappa or lambda constant region
- VH–linker–VL scFv
- VL–linker–VH scFv
- VHH–Fc fusion
- Fab-like construct
- Bispecific or multivalent architecture in later phases

The system creates a digital construct manifest and performs:

- Component-type validation
- Domain-order validation
- Chain-role validation
- Junction checks
- Linker checks
- Numbering checks
- Constant-region compatibility checks
- Basic architecture precedent checks
- Independent component comparison
- Pair and assembly precedent searches
- Transparent warnings

The system must not represent successful digital assembly as evidence of physical function.

---

# 5. Antibody Representation: Use an Assembly Graph

Do not force every antibody into a single fixed heavy-chain/light-chain tree.

Represent the candidate as an assembly graph.

```text
Project
└── Candidate Assembly
    ├── Chain
    │   ├── Domain
    │   │   ├── Region
    │   │   └── Residue Set
    │   ├── Junction
    │   ├── Linker
    │   └── Tag or Fusion
    ├── Chain Relationship
    ├── Pairing Relationship
    ├── Copy Number
    ├── Format Definition
    └── Evidence Profiles
```

Required primary entities:

- Project
- Input artifact
- Candidate assembly
- Chain
- Domain
- Region
- Residue
- Residue set
- Paratope
- Junction
- Linker
- Tag
- Fusion partner
- Constant domain
- Fc region
- Pair relationship
- Format definition
- Component match
- Evidence record
- Source snapshot
- Alignment
- Structure-linked record
- Score result
- Warning
- User assertion
- Database-derived assertion
- Prediction
- Analysis run

Every component must have:

- A stable internal identifier
- Parent-child relationships
- Chain role
- Component type
- Sequence coordinates
- Canonical numbering coordinates where applicable
- Source of the component
- Original and normalized sequence
- Evidence profile
- Validation status

---

# 6. Supported Component Ontology

The ontology must distinguish at minimum:

## 6.1 Chain-level components

- Full heavy chain
- Full light chain
- Kappa light chain
- Lambda light chain
- Heavy-chain-only antibody chain
- Engineered fusion chain
- Unknown chain

## 6.2 Variable domains

- VH
- VL-kappa
- VL-lambda
- VHH
- VNAR, only if explicitly supported
- Unknown immunoglobulin-like variable domain

## 6.3 CDRs

Using one canonical internal definition:

- CDR-H1
- CDR-H2
- CDR-H3
- CDR-L1
- CDR-L2
- CDR-L3

## 6.4 Framework regions

- FR-H1
- FR-H2
- FR-H3
- FR-H4
- FR-L1
- FR-L2
- FR-L3
- FR-L4

## 6.5 Constant regions

- CH1
- Hinge
- CH2
- CH3
- CH4 where applicable
- CL-kappa
- CL-lambda
- Fc
- Isotype
- Subclass
- Species-specific constant region
- Engineered constant-region variant

## 6.6 Construct-level elements

- Signal peptide, when supplied
- Linker
- Junction
- Peptide tag
- Protease-cleavage site
- Fusion partner
- Spacer
- Transmembrane region, only in later supported formats
- Copy number
- Domain order
- Chain stoichiometry

## 6.7 Residue-set components

- User-defined paratope residue set
- Experimentally reported paratope residue set from a curated source
- Experimental epitope residue set, optional and contextual
- User-defined epitope residue set, optional and contextual

Residue sets must retain parent-chain or parent-antigen coordinates and the numbering scheme used.

For the first version, residue-set analysis is position- and sequence-based. Three-dimensional patch geometry, coordinate extraction, molecular-surface analysis, and predicted paratope generation are deferred to future versions.

# 7. Input Contract

## 7.1 Sequence inputs

Support:

- Amino-acid FASTA
- Plain amino-acid sequence
- Optional nucleotide FASTA
- Multiple named chains
- Individual components
- Parent sequence plus component coordinates
- AIRR-compatible rearrangement records where applicable

For nucleotide input, preserve:

- Original nucleotide sequence
- Reading frame
- Translation
- Codon coordinates
- V, D, and J assignment
- Junction analysis
- Somatic-mutation information
- Translation warnings

Amino-acid-only inputs must receive lower confidence for conclusions that depend on nucleotide rearrangement history.

## 7.2 Optional structure-record references

The first version does not require or process molecular coordinates.

Allow the user to provide optional identifiers such as:

- PDB accession
- SAbDab2 record identifier
- Thera-SAbDab-linked structure identifier
- Publication identifier describing a solved antibody structure

Use these identifiers only to retrieve and display structure-linked experimental metadata, including:

- Whether an experimental structure exists
- Experimental method
- Resolution or other reported quality field
- Heavy–light pairing recorded in the structure
- Antigen-complex status
- Chain identifiers
- Missing-residue annotations when directly supplied by the source
- Publication and database links

Do not, in the first version:

- Upload or parse PDB or mmCIF coordinates
- Predict antibody structures
- Align structures
- Calculate RMSD
- Compare CDR loop geometry
- Calculate VH–VL orientation
- Generate molecular surfaces
- Extract contact residues from coordinates
- Display an interactive three-dimensional model

## 7.3 Individual short-component inputs

For a CDR, short framework segment, linker, or paratope input:

- Require the user to select the expected component type.
- Require a parent sequence or numbered positions when necessary.
- Record whether the component was extracted by the platform or supplied directly.
- Warn that short exact matches may occur by chance.
- Do not use ordinary full-length BLAST scoring unchanged.
- Do not infer a parent chain from a very short string with high confidence.

## 7.4 Paratope input modes for the first version

### User-defined paratope

Require:

- Parent heavy and/or light sequence
- Chain designation
- Numbered residue positions
- Canonical or supplied numbering scheme
- Residue identities
- Source of the assignment
- Evidence label: experimentally reported, literature-derived, or hypothetical

### Database-reported experimental paratope

Use only when a supported database or curated source directly provides antibody contact-residue annotations or a residue list.

Store:

- Parent antibody record
- Heavy and light chain identifiers
- Numbered antibody residue positions
- Antigen identifier, when known
- Contact-definition description supplied by the source
- Experimental record and publication
- Missing-data flags

The first version must not derive paratopes by processing molecular coordinates. It may compare user-defined residue sets with database-reported residue sets using numbered-position overlap, residue identity, residue-class similarity, and parent-sequence context.

Predicted paratopes, structural patch matching, spatial contact geometry, and coordinate-derived paratope extraction are deferred to future versions.

## 7.5 Candidate-assembly inputs

Require:

- Assembly name
- Intended format
- Chain list
- Domain order
- Component identifiers
- Linker sequences
- Constant-region choices
- Chain copy numbers
- Orientation where relevant
- User-supplied intended species
- User-supplied intended target, only as metadata

The intended target must not influence component evidence scores unless the user explicitly runs the optional target-context module.

---

# 8. Input Validation and Normalization

Before analysis:

- Preserve the original input exactly.
- Remove valid FASTA headers only from the analysis sequence.
- Convert valid residues to uppercase.
- Reject invalid characters.
- Detect nucleotide-versus-protein input.
- Detect stop codons.
- Detect internal translation failures.
- Calculate sequence length.
- Calculate SHA-256 checksum.
- Detect duplicate inputs.
- Detect mixed or concatenated chains.
- Detect likely truncation.
- Detect implausible lengths.
- Detect ambiguous chain identity.
- Detect missing conserved residues.
- Detect unexpected cysteine patterns.
- Detect unrecognized domain boundaries.
- Detect incomplete constant regions.
- Detect inconsistent kappa/lambda assignments.
- Detect incompatible domain order.
- Detect linker omission where a single-chain format requires one.
- Detect impossible or unsupported format declarations.
- Record every normalization action.
- Never silently replace, delete, or impute residues.

Required validation statuses:

- Valid
- Valid with warnings
- Ambiguous
- Incomplete
- Unsupported format
- Failed validation

Every warning must include:

- Rule identifier
- Rule version
- Affected coordinates
- Severity
- Explanation
- Whether analysis can continue
- Effect on downstream confidence

---

# 9. Canonical Numbering and Coordinate Policy

Use **IMGT numbering as the canonical internal coordinate system** for supported antibody variable domains.

Allow display conversion to:

- Kabat
- Chothia
- Martin
- AHo

Rules:

1. Store the original sequence coordinates.
2. Store canonical IMGT coordinates.
3. Store insertion codes.
4. Store mappings to alternate numbering schemes.
5. Record the tool and version used.
6. Record numbering success, failure, and confidence.
7. Never combine CDRs extracted under different definitions without remapping.
8. Use the canonical IMGT-defined regions for internal CDR and framework scoring unless a project explicitly selects a different versioned policy.
9. Show users which definition was used.
10. Preserve alternative boundary definitions for comparison, but do not mix them in one score.

Recommended implementation:

- Primary local numbering: ANARCII
- Fallback or cross-check: ANARCI
- Germline and junction analysis: local IgBLAST
- Disagreement between tools: retain both results and lower annotation confidence

---

# 10. Database and Tool Access Policy

For every external resource, create and maintain a database-access record with:

- Resource name
- Resource owner
- Scientific purpose
- Official access method
- API availability
- Bulk-download availability
- Local-software availability
- Authentication requirements
- Rate limits
- Academic-use restrictions
- Commercial-use restrictions
- Redistribution restrictions
- Citation requirements
- Update schedule
- Version identifier
- Snapshot checksum
- Parser version
- Last legal review date
- Last technical validation date
- Contact information
- Fallback plan

Priority order:

1. Official bulk download or versioned release
2. Official API
3. Official locally runnable software
4. Official web service for manual prototype use
5. Authorized licensed data
6. Exclude when compliant production access is unavailable

Never:

- Invent an API
- Scrape a resource whose terms prohibit automated use
- Assume that a public website permits unrestricted production use
- Assume that scientific data and software share the same license
- Submit confidential sequences to third-party web tools without explicit permission
- Allow a live database update to silently alter an old analysis

---

# 11. Required and Optional Resources

The platform must use a source-snapshot architecture. The access assumptions below must be reverified before production launch and at every material terms-of-use change.

## 11.1 IMGT data

Use for:

- Standard terminology
- Immunoglobulin gene nomenclature
- Germline V, D, J, and constant-region references
- Alleles
- Species-specific context
- IMGT region definitions
- IMGT numbering conventions

Production rule:

- Use authorized, versioned data snapshots.
- Do not automate IMGT web tools unless permission and terms are explicitly confirmed.
- Store source version, retrieval date, license, and checksum.
- Use local IgBLAST and numbering tools for routine processing.

## 11.2 OGRDB

Use for:

- Community-curated germline sets
- Versioned species and locus-specific reference sets
- Alternative or updated germline alleles
- Germline-reference comparison
- Conflict detection between reference sources

Access:

- Official downloads
- Official REST API

Do not silently merge OGRDB and IMGT assignments. Show:

- Reference source
- Germline-set version
- Best assignment under each reference set
- Conflicting allele calls
- Confidence and limitations

## 11.3 NCBI IgBLAST

Use local standalone IgBLAST for:

- V assignment
- D assignment
- J assignment
- Optional C assignment
- Germline alignment
- Junction analysis
- Framework and CDR delineation
- Mutation counts
- Alternative assignments
- Nucleotide and protein analysis

Record:

- IgBLAST version
- Germline database version
- Command parameters
- Raw output
- Parsed output
- Assignment confidence
- Alternative hits

## 11.4 ANARCII and ANARCI

Use locally for:

- Chain classification
- Variable-domain recognition
- IMGT numbering
- Alternate numbering mappings
- CDR and framework extraction

Do not depend on public web forms for production throughput.

## 11.5 Observed Antibody Space

Use OAS for:

- Unpaired heavy-chain repertoire sequences
- Unpaired light-chain repertoire sequences
- Available paired sequences
- VH and VL sequences
- CDRs
- Frameworks
- Germline usage
- Species
- Study and subject context
- Isotype where available
- Repertoire observation counts

Production rule:

- Use bulk or filtered official downloads.
- Build a local, versioned search index.
- Retain study-level metadata.
- Separate paired from unpaired records.
- Correct for clonal expansion, repeated subjects, and sequencing depth.
- Never interpret unpaired records as pairing evidence.

## 11.6 AIRR Data Commons

Use as an optional or expansion source for:

- AIRR-standardized repertoire datasets
- Study metadata
- Rearrangement records
- Repository-level queries
- Programmatic access through the AIRR Data Commons API

Store repository identity and dataset-level terms because the API standard does not automatically make all underlying studies equivalent in quality or reuse conditions.

## 11.7 SAbDab2

Use for:

- Antibody sequences linked to experimentally determined structures
- Direct heavy–light chain pairing in structural records
- CDR sequences
- Antibody format
- Antigen-complex status
- Experimental method
- Resolution or other supplied quality metadata
- PDB cross-references
- Reported contact-residue annotations when directly available
- Versioned structure-linked datasets

First-version production rule:

- Use versioned downloads or official access methods.
- Build local sequence and metadata indexes.
- Use SAbDab2 as a structure-linked evidence source, not as a three-dimensional modeling engine.
- Do not ingest or process coordinates for scoring in the MVP.
- Do not calculate structural alignment, RMSD, CDR geometry, or molecular-surface similarity.

## 11.8 RCSB PDB and wwPDB

Use for structure-record metadata and provenance:

- PDB identifier
- Experimental method
- Resolution or applicable quality field
- Deposited chain sequences
- Biological assembly metadata
- Chain identifiers
- Antigen-complex status
- Missing-residue annotations when provided
- Engineered mutation annotations when provided
- Validation-report links
- Publication
- Revision history

Use official APIs and bulk metadata files.

For the first version:

- Retrieve metadata and deposited sequences only.
- Do not build a coordinate-processing pipeline.
- Do not perform structural alignment or visualization.
- Do not interpret the existence of a solved structure as proof that the submitted sequence adopts the same structure.

## 11.9 Thera-SAbDab

Use for:

- Public therapeutic antibody variable-region sequences
- WHO-recognized antibody-derived therapeutics with released sequences
- Therapeutic names
- Targets
- Formats
- Development status
- Available structures
- Sequence-pair precedent

Use versioned downloadable tables where allowed.

Before commercial redistribution:

- Recheck the exact terms for the database and any upstream sequence sources.
- Preserve attribution.
- Do not assume that the therapeutic status transfers to a similar candidate.

## 11.10 UniProtKB

Use for:

- Constant-region proteins
- Fc-related proteins
- Fusion proteins
- Antigens
- Protein names
- Species
- Domain features
- Sequence variants
- Post-translational modifications
- Evidence-tagged functional annotations
- Cross-references

Prefer reviewed entries when available, but retain reviewed status as evidence metadata rather than treating it as proof of relevance.

Use the official REST API or versioned bulk downloads.

## 11.11 IEDB

Keep IEDB optional and secondary.

Use for:

- Experimental B-cell assay context
- Antibody or BCR receptor records where available
- Epitopes
- Antigens
- Positive and negative assay records
- Three-dimensional epitope context
- Publications

Use official database exports or the official query API.

Do not let epitope evidence dominate antibody component evidence.

## 11.12 PubMed and Europe PMC

Use for:

- Publication identifiers
- Titles
- Authors
- Abstract metadata
- Links to source publications
- Open-access status
- License-aware full-text links where available

Use:

- NCBI E-utilities for PubMed
- Europe PMC REST services for publication and open-access metadata

Do not copy or redistribute full text unless the article license permits it.

## 11.13 NCBI BLAST+ and local alignment tools

Use local BLAST+ only where appropriate:

- Full heavy chains
- Full light chains
- VH
- VL
- Constant regions
- Fc regions
- Longer frameworks
- General protein searches

Do not use ordinary full-length BLAST scoring unchanged for short CDRs, short linkers, or residue-set paratopes.

## 11.14 Optional secondary resources

Evaluate separately and do not make the MVP dependent on them:

- PLAbDab for patent and literature-derived paired sequences
- ABSD as a standardized secondary sequence index
- Additional licensed repertoire datasets
- Additional therapeutic-antibody sources
- Specialized developability datasets
- Specialized nanobody datasets
- Specialized Fc-engineering literature datasets

For patent-derived data, perform a separate legal and provenance review. Public patent disclosure does not automatically establish unrestricted freedom to use or commercialize a sequence.

---

# 12. Component-to-Source Access Matrix

The system must explicitly map every component to accessible evidence sources.

| Component | Primary sources | Secondary sources | Required caveat |
|---|---|---|---|
| Full heavy chain | OAS, SAbDab2, Thera-SAbDab, local BLAST index | AIRR Data Commons, ABSD, PLAbDab | Separate variable-only from full-chain matches |
| Full light chain | OAS, SAbDab2, Thera-SAbDab | AIRR Data Commons, ABSD, PLAbDab | Preserve kappa/lambda class |
| VH | OAS, SAbDab2, Thera-SAbDab | AIRR Data Commons, ABSD | Compare VH to VH; VHH only as related class |
| VL-kappa | OAS, SAbDab2, Thera-SAbDab | AIRR Data Commons | Do not merge with lambda without penalty and labeling |
| VL-lambda | OAS, SAbDab2, Thera-SAbDab | AIRR Data Commons | Preserve lambda identity |
| VHH | SAbDab2, OAS where present, Thera-SAbDab | Specialized nanobody resources | Do not score as ordinary VH without class labeling |
| V, D, J germline | IMGT data, OGRDB, IgBLAST | NCBI references | Preserve reference-set disagreements |
| CDR-H1/H2/L1/L2/L3 | OAS, SAbDab2, Thera-SAbDab | AIRR Data Commons | Same-region, same-definition comparison only |
| CDR-H3 | OAS, SAbDab2, Thera-SAbDab | AIRR Data Commons | Length-matched empirical null model required |
| Frameworks | OAS, SAbDab2, Thera-SAbDab, germline sets | UniProt where applicable | Numbered-position comparison |
| Constant region | IMGT constant references, UniProt, PDB | Therapeutic records | Separate isotype, subclass, species, and domain |
| Fc | UniProt, PDB, IMGT references | Thera-SAbDab and curated literature | No comprehensive open source proves Fc function |
| Linker | SAbDab2/PDB-linked construct metadata, UniProt fusions, therapeutic records | Curated internal linker library | Use sequence and construct metadata; no three-dimensional linker modeling |
| Junction | User input, IgBLAST, paired public records | AIRR records | Sequence precedent does not prove expression |
| Heavy–light pair | SAbDab2, paired OAS, Thera-SAbDab | IEDB BCR records, PLAbDab | Only direct paired records count as observed pairing |
| Paratope | Database-reported antibody contact residues from SAbDab2/PDB-linked or IEDB records when available | Curated literature | First version uses numbered-position and residue comparison only; no coordinate-derived extraction |
| Antigen | UniProt, PDB, IEDB | NCBI Protein | Optional contextual layer |
| Epitope | IEDB, PDB/SAbDab2 complexes | Literature | Must remain separate from antibody component scores |
| Publications | PubMed, Europe PMC | Cross-references in source databases | Full-text reuse is license-dependent |

When a component lacks a comprehensive public database, say so explicitly. Do not fabricate completeness.

---

# 13. Source Snapshot and Data-Ingestion Architecture

Required pipeline:

```text
Official resource
      ↓
Legally reviewed source record
      ↓
Versioned immutable raw snapshot
      ↓
Checksum and manifest
      ↓
Validated parser
      ↓
Normalized source-specific tables
      ↓
Canonical antibody component model
      ↓
Deduplication and clustering
      ↓
Component-specific search indexes
      ↓
Scoring and evidence services
```

For every snapshot preserve:

- Source name
- Release or snapshot identifier
- Retrieval timestamp
- License and terms version
- Raw-file checksum
- Record count
- Parser version
- Schema version
- Normalization version
- Deduplication version
- Index version
- Known ingestion warnings
- Update status
- Previous snapshot relationship

Old analysis runs must remain reproducible after a source update.

A new source release must trigger:

- Parser regression tests
- Record-count checks
- Schema-drift detection
- Deduplication comparison
- Score-shift analysis
- Benchmark rerun
- Release notes
- Optional reanalysis, never silent replacement

---

# 14. Data Normalization Rules

Normalize every external record into a common evidence schema while retaining the raw source.

Required normalized fields:

- Source database
- Source snapshot
- Source accession
- Source record type
- Parent antibody identifier
- Chain identifier
- Component type
- Original sequence
- Normalized sequence
- Original coordinates
- IMGT coordinates
- Alternate numbering coordinates
- Species
- Germline assignment
- Isotype
- Kappa/lambda class
- Paired status
- Parent heavy sequence
- Parent light sequence
- Antibody format
- Natural, engineered, therapeutic, or unknown status
- Structure-record availability
- PDB identifiers
- Antigen
- Epitope
- Assay type
- Positive or negative assay result
- Publication identifiers
- Evidence directness
- Metadata completeness
- Data-quality flags
- License and attribution fields

Never overwrite contradictory metadata. Store:

- Source-specific assertion
- Normalized assertion
- Conflict status
- Resolution rule
- Resolution confidence

---

# 15. Deduplication and Independence

The platform must prevent duplicate evidence inflation.

Deduplicate at several levels:

## 15.1 Exact sequence

Cluster identical component sequences.

## 15.2 Parent antibody

Multiple PDB structures of the same antibody count as one antibody precedent plus multiple structural observations.

## 15.3 Publication

The same experiment represented in several databases counts as one experiment.

## 15.4 Study and subject

Repeated repertoire records from one subject or one expanded clone must not count as independent biological observations.

## 15.5 Therapeutic aliases

Different names, development codes, or records for the same therapeutic must map to one canonical therapeutic entity.

## 15.6 Structure-record redundancy

Multiple PDB entries, crystal forms, chains, or database representations of the same underlying antibody must not dominate the evidence profile.

Count separately:

- Canonical parent antibody
- Independent experimental structure record
- Independent publication
- Distinct antigen-complex record

A large number of deposited structures for one antibody must not be interpreted as broad independent component precedent.

Every evidence count must show:

- Raw record count
- Deduplicated record count
- Independent study count
- Independent subject count
- Independent parent-antibody count
- Independent structure-linked record count
- Clustering rule
- Cluster threshold
- Deduplication version

# 16. Search and Indexing Architecture

Do not rely on one universal search method.

## 16.1 Exact sequence index

Use hash-based exact lookup for:

- Exact full-chain matches
- Exact VH/VL matches
- Exact CDR matches
- Exact framework matches
- Exact linker matches
- Exact constant-domain matches

## 16.2 Approximate long-sequence index

Use a fast local protein sequence search method to retrieve candidates for:

- Full heavy chains
- Full light chains
- VH
- VL
- VHH
- Constant regions
- Fc
- Longer framework segments

Then refine top candidates with a deterministic pairwise alignment and store the exact alignment.

## 16.3 Numbered-position index

Use canonical positions for:

- Framework comparison
- Variable-domain comparison
- CDR position comparison
- Vernier-zone comparison
- Interface-residue comparison
- Constant-domain positional checks

## 16.4 Short-CDR index

Partition by:

- CDR type
- Length
- Chain class
- Species
- Germline family when appropriate
- Numbering definition

Use position-aware alignment and empirical null distributions.

## 16.5 Pair index

Store directly observed pairs:

- Heavy sequence
- Light sequence
- VH
- VL
- Combined CDR pattern
- Germline pair
- Parent record
- Study
- Structure-record availability
- Therapeutic status

Do not create observed pair records by joining unrelated heavy and light matches.

## 16.6 Structure-linked metadata index

For antibody records associated with experimental structures, index:

- Parent antibody identifier
- Heavy and light sequences
- Component sequences
- Direct heavy–light pairing
- PDB identifier
- Experimental method
- Resolution or supplied quality field
- Antigen-complex status
- Antigen name, when available
- Missing-residue annotations supplied by the source
- Reported contact residues, only when directly available as curated data
- Publication
- Source snapshot

Search these records using sequence and component metadata. Do not search molecular coordinates in the first version.

## 16.7 Paratope residue-set index

Store supported paratopes as numbered antibody residue sets anchored to parent sequences.

Index:

- Numbered position set
- Residue identities
- Chain distribution
- CDR/framework distribution
- Parent VH/VL sequences
- Experimental, literature-derived, or hypothetical status
- Antigen class, only as contextual metadata
- Source record and publication

Compare by numbered-position overlap, residue identity, conservative residue-class similarity, CDR/framework distribution, and parent-sequence context.

Spatial patch geometry, molecular-surface comparison, and coordinate-based contact extraction are deferred.

# 17. Evidence Classes

Every result must belong to one or more explicit evidence classes:

1. Germline evidence
2. Curated antibody sequence evidence
3. Observed repertoire evidence
4. Direct paired-chain evidence
5. Structure-linked experimental record evidence
6. Experimental antibody–antigen complex metadata
7. Therapeutic precedent
8. Constant/Fc reference evidence
9. Literature-linked experimental evidence
10. Optional antigen or epitope context
11. User-supplied evidence
12. Computational prediction
13. Platform inference

Do not combine these into one opaque number before presenting them separately.

---

# 18. Raw Measurements First

For every sequence match show:

- Exact-match status
- Alignment algorithm
- Substitution matrix
- Gap penalties
- Alignment length
- Query length
- Subject length
- Sequence identity
- Conservative similarity
- Query coverage
- Subject coverage
- Insertions
- Deletions
- Length difference
- E-value, where statistically valid
- Bit score, where applicable
- Canonical-position agreement
- Germline agreement
- Component-type agreement

For every structure-linked matched record show:

- PDB or source identifier
- Parent antibody identifier
- Matched component and its sequence alignment
- Experimental method
- Resolution or applicable quality field
- Direct heavy–light pairing status
- Antigen-complex status
- Deposited chain identifiers
- Missing-residue annotations when supplied by the source
- Reported contact residues when directly available as curated metadata
- Publication
- Source snapshot

These fields describe the supporting experimental record. They do not constitute a structural comparison with the submitted candidate.

For every paratope residue-set match show:

- Input residue-set size
- Matched residue-set size
- Numbered-position intersection
- Jaccard overlap of numbered positions
- Exact residue matches at overlapping positions
- Conservative residue-class similarity
- Heavy/light chain distribution
- CDR/framework distribution
- Parent-sequence similarity
- Experimental, literature-derived, or hypothetical status
- Source and publication

Do not calculate or display RMSD, structural alignment, loop geometry, surface accessibility, spatial patch alignment, or predicted-model confidence in the first version.

---

# 19. Scoring Architecture

Do not present uncalibrated heuristic sums as biological truth.

Use four output layers.

## Layer 1: Raw measurements

The direct measurements listed above.

## Layer 2: Component-specific calibrated percentiles

Convert match measurements into empirical percentiles relative to the same component class.

Examples:

- CDR-H3 compared with CDR-H3 of similar length
- VH compared with VH
- Kappa VL compared with kappa VL
- Fc compared with the same Fc domain and species context
- Paratope residue set compared with database-reported residue sets of similar size and numbering coverage

```text
Calibrated Percentile =
100 × empirical percentile of the observed measurement
within a versioned reference distribution
```

Each percentile must name:

- Reference distribution
- Component class
- Database snapshot
- Stratification variables
- Sample size
- Calibration version

## Layer 3: Evidence support grade

Use an interpretable grade:

- A: Direct, curated, same-component experimental evidence with strong provenance
- B: Strong curated or independently repeated component precedent
- C: Valid sequence precedent with limited experimental, structure-linked, or pairing context
- D: Weak, indirect, predicted, sparse, or poorly annotated evidence
- Insufficient: Not enough evidence to grade

The exact grade thresholds must be benchmarked and versioned.

## Layer 4: Evidence sufficiency and uncertainty

Report:

- Independent record count
- Parent-antibody count
- Study count
- Subject count
- Paired-record count
- Structure-linked record count
- Therapeutic-record count
- Metadata completeness
- Database coverage
- Bootstrap interval where applicable
- Confidence in the score
- Missing evidence

---

# 20. Deterministic Initial Match Formula

Before empirical calibration is complete, use a transparent provisional sequence-match measurement, not a claimed biological score.

For sufficiently long components:

```text
Coverage Harmonic Mean =
2 × query_coverage × subject_coverage
/
(query_coverage + subject_coverage)

Effective Sequence Match =
100 × sequence_identity × sqrt(Coverage Harmonic Mean)
```

Where all fractions range from 0 to 1.

Also report conservative substitution similarity separately.

Do not hide poor coverage behind high aligned-region identity.

For short CDRs, do not use this formula as the final significance measure.

For a hit to be eligible:

- Component type must match exactly, or
- A related class must be clearly labeled and down-ranked, or
- The hit must be excluded

Examples:

- CDR-H3 to CDR-H3: eligible
- CDR-H3 to CDR-L3: exclude
- VH to VHH: related-class comparison only
- Kappa VL to lambda VL: related light-chain comparison, not exact-class evidence
- Full chain to variable domain: parent/child comparison, not full-chain identity

---

# 21. Short-CDR Chance Correction

For every CDR type, create an empirical null model.

Stratify background comparisons by:

- CDR type
- Sequence length
- Chain class
- Species where possible
- Germline family where appropriate
- Database snapshot
- Redundancy-correction policy

For the observed match:

```text
Empirical p-value =
(1 + number of null comparisons at least as strong as observed)
/
(1 + total null comparisons)
```

Then apply false-discovery-rate correction across the searched candidate set.

Report:

- Position-aware similarity
- Length compatibility
- Empirical p-value
- Adjusted q-value
- Similarity percentile
- Number of comparisons
- Null-model version
- Database redundancy policy

Short exact matches must not automatically receive a high biological-precedent interpretation.

## CDR-H3 special handling

For CDR-H3 include:

- Length
- Position-aware similarity
- Junction composition
- Insertions and deletions
- Charge
- Hydrophobicity
- Aromatic residue content
- Structure-linked experimental record precedent
- Observed repertoire frequency
- Parent VH context
- Paired VL context, when directly observed

Never infer shared specificity from CDR-H3 similarity alone.

---

# 22. Per-Component Evidence Profiles

Every component receives its own profile.

Required profile dimensions:

1. Annotation confidence
2. Best calibrated sequence-match percentile
3. Observed Repertoire Precedent
4. Structure-Linked Experimental Precedent
5. Therapeutic Precedent
6. Direct Pairing Precedent, when applicable
7. Database-Relative Novelty
8. Evidence Sufficiency
9. Analysis Completion
10. Warnings and limitations

Do not force inapplicable dimensions into a zero.

Use statuses:

- Available
- Not found
- Not applicable
- Insufficient evidence
- Analysis failed
- Database unavailable
- Access restricted
- Unsupported component

---

# 23. Annotation Confidence

Calculate annotation confidence from independent signals:

- Numbering completeness
- Chain-type agreement
- Germline alignment support
- Expected architecture
- Conserved-residue support
- Tool agreement
- Sequence completeness

Provisional formula:

```text
Annotation Confidence =
100 × (
  0.25 × numbering completeness
+ 0.20 × chain-type agreement
+ 0.20 × germline support
+ 0.15 × architecture support
+ 0.10 × conserved-residue support
+ 0.10 × cross-tool agreement
)
```

This is a provisional engineering score and must be calibrated.

Display every normalized input and failed check.

---

# 24. Observed Repertoire Precedent

Calculate separately for:

- Heavy chain
- Light chain
- VH
- VL
- VHH
- Each CDR
- Each framework
- Direct heavy–light pair

Recommended evidence dimensions:

- Best calibrated similarity percentile
- Number of nonredundant neighboring sequences
- Independent study count
- Independent subject count
- Paired-record support
- Species relevance
- Germline-context relevance
- Dataset metadata quality
- Sequencing-depth correction
- Clonal-expansion correction

Provisional score:

```text
Observed Repertoire Precedent =
100 × (
  0.45 × best calibrated similarity
+ 0.15 × study breadth
+ 0.15 × subject breadth
+ 0.10 × nonredundant neighborhood support
+ 0.10 × context relevance
+ 0.05 × metadata quality
)
```

All breadth terms must use saturating transforms so one large study cannot dominate.

Example saturating transform:

```text
breadth_support = 1 - exp(-independent_count / scale_parameter)
```

Show the scale parameter and calibration version.

Do not interpret high repertoire precedent as proof of low immunogenicity or functional expression.

---

# 25. Structure-Linked Experimental Precedent

The first version must use structural databases without performing structural modeling or geometry comparison.

This evidence answers:

> Has the same or a similar antibody component appeared in an experimentally determined structure record, and how strong is the sequence and metadata support for that relationship?

Calculate separately for:

- Heavy chain
- Light chain
- VH
- VL
- VHH
- Each CDR
- Each framework
- Constant domain
- Fc
- Direct heavy–light pair
- User-defined paratope residue set when curated contact annotations are available

Use:

- Calibrated sequence similarity to the structure-linked component
- Alignment coverage
- Exact component-type relevance
- Direct heavy–light pairing in the structure record
- Experimental method
- Resolution or other reported quality field
- Antigen-complex status
- Independent parent-antibody and structure-record support
- Metadata completeness

Provisional score:

```text
Structure-Linked Experimental Precedent =
100 × (
  0.50 × calibrated component sequence similarity
+ 0.15 × alignment coverage
+ 0.10 × exact component-type relevance
+ 0.10 × direct pairing or parent-record relevance
+ 0.10 × experimental-record quality
+ 0.05 × independent structure-linked support
)
```

This score does not measure three-dimensional similarity. It measures sequence-based similarity to components with experimental structure records and the quality of the supporting metadata.

Always display:

- Matched sequence and alignment
- PDB or source identifier
- Experimental method
- Resolution or applicable quality field
- Direct pair status
- Antigen-complex status
- Publication
- Missing metadata
- Important sequence differences

Never claim that a high score proves the submitted component has the same fold details, CDR conformation, paratope geometry, antigen target, or function.

# 26. Therapeutic Precedent

Calculate separately for each supported component.

Use:

- Component sequence similarity
- Coverage
- Same-region confirmation
- Format relevance
- Direct paired status
- Structure-record availability
- Therapeutic-record completeness

Provisional score:

```text
Therapeutic Precedent =
100 × (
  0.45 × calibrated component similarity
+ 0.15 × coverage
+ 0.15 × exact component-class confirmation
+ 0.10 × format relevance
+ 0.10 × direct pair or structure-linked record support
+ 0.05 × metadata completeness
)
```

Display:

- Therapeutic name
- Matched component
- Identity
- Coverage
- Important sequence differences
- Target
- Format
- Public development status
- Structure
- Source snapshot
- Publications

Always state:

> Similarity to a therapeutic antibody does not transfer approval, safety, efficacy, manufacturability, target specificity, or clinical performance.

---

# 27. Database-Relative Novelty

Do not use only:

```text
Novelty = 100 - nearest-neighbor similarity
```

Calculate database-relative novelty from:

- Nearest-neighbor distance
- Local neighborhood density
- Nonredundant neighbor count
- Independent study breadth
- Species-specific database coverage
- Germline-context coverage
- Length rarity
- Insertion/deletion rarity
- Component-class data sparsity

Required outputs:

- Novelty category
- Novelty percentile
- Nearest known match
- Local neighborhood size
- Reference database coverage
- Confidence in novelty estimate
- Reason for novelty
- Warning when novelty may be caused by sparse data

Suggested categories:

- Low database-relative novelty
- Moderate database-relative novelty
- High database-relative novelty
- Indeterminate because of sparse reference coverage

---

# 28. Evidence Completeness Must Be Split

Do not use one ambiguous completeness score.

Report four separate fields.

## 28.1 Input completeness

Which expected components were supplied?

## 28.2 Analysis completion

Which computational modules completed successfully?

## 28.3 Metadata completeness

How complete are the matched source records?

## 28.4 Evidence sufficiency

Is there enough independent evidence to support interpretation?

A database outage affects analysis completion. It must not be treated as biological evidence against the candidate.

---

# 29. Full-Chain Analysis

Analyze heavy and light chains separately.

For heavy chains:

- Compare with heavy chains only.
- Distinguish variable-only from full-chain records.
- Preserve species and isotype context.
- Distinguish ordinary heavy chains from VHH or heavy-chain-only formats.
- Show germline assignment.
- Show constant-region coverage separately.

For light chains:

- Compare with light chains only.
- Preserve kappa versus lambda class.
- Distinguish variable-only from full-chain records.
- Show germline assignment.
- Show constant-light compatibility separately.

Never allow a strong variable-domain match to masquerade as a full-chain match.

---

# 30. Variable-Domain Analysis

Return separate profiles for:

- VH
- VL-kappa
- VL-lambda
- VHH

Compare:

- Full-domain alignment
- Framework similarity
- CDR lengths
- Numbered positions
- Germline family
- Species
- Structure-record availability
- Paired-chain availability
- Therapeutic precedent

Related-class comparisons must be labeled.

Example:

- VH versus VHH: related class
- Kappa VL versus lambda VL: related light-chain class
- VH versus VL: exclude

---

# 31. Framework Analysis

Analyze independently:

- FR-H1
- FR-H2
- FR-H3
- FR-H4
- FR-L1
- FR-L2
- FR-L3
- FR-L4

Use numbered-position comparison.

Report:

- Closest germline
- Closest observed framework
- Closest therapeutic framework
- Position-specific mutations
- Insertions and deletions
- Vernier-zone differences
- Residues near CDRs
- Structure-linked experimental precedent
- Unusual residues
- Known motif flags

Do not automatically classify a framework mutation as harmful.

---

# 32. Constant and Fc Analysis

Analyze constant regions independently from variable regions.

Compare only with:

- Equivalent constant domain
- Equivalent isotype or subclass
- Equivalent species
- Known engineered variants
- Experimental structures
- Therapeutic records where sequence is available

Report:

- Domain identity
- Isotype and subclass
- Species
- Conserved cysteines
- Glycosylation motifs
- Hinge length and composition
- Engineered substitutions
- Structure-linked experimental precedent
- Literature-linked annotations
- Important uncertainty

There is no single complete public source that proves the function of every Fc variant. The platform must create a curated, provenance-preserving internal Fc evidence library from accessible reference sequences, structure-linked records, therapeutic records, and literature.

Do not predict effector function solely from sequence similarity.

---

# 33. Linker and Junction Analysis

Treat linkers and junctions as separate construct elements.

For linkers report:

- Exact and similar sequence occurrences in structure-linked or therapeutic constructs
- Length
- Composition
- Flexibility-related descriptors
- Charge
- Hydrophobicity
- Repeated motifs
- Junction residues
- Structure-linked construct metadata where available
- Source and parent construct

Do not claim that a familiar linker guarantees expression, folding, or function.

Because no single comprehensive authoritative public linker database exists, create a versioned internal linker precedent library derived from:

- SAbDab2/PDB-linked construct records
- UniProt fusion proteins
- Public therapeutic records
- Curated literature
- User-approved internal records

For junctions:

- Validate reading frame for nucleotide constructs.
- Check unintended stop codons.
- Check domain-boundary preservation.
- Record added or removed residues.
- Compare exact junctions only when parent format is equivalent.
- Do not treat junction similarity across unrelated formats as direct evidence.

---

# 34. Paratope Residue-Set Analysis

The first version supports paratopes only as numbered antibody residue sets attached to parent sequences.

It does not process molecular coordinates, predict paratopes, compare three-dimensional patches, or infer antigen contacts.

## 34.1 Supported paratope evidence

Use:

- User-defined numbered residue sets
- Experimentally reported contact-residue lists directly available from supported databases
- Curated literature-derived residue lists with publication provenance

Every paratope record must retain:

- Parent heavy and light sequences
- Chain designation
- Numbering scheme
- Numbered positions
- Residue identities
- CDR/framework classification
- Evidence type
- Antigen metadata, when known
- Source and publication

## 34.2 Paratope comparison dimensions

Compare:

- Numbered-position overlap
- Jaccard overlap of residue positions
- Exact residue identity at shared positions
- Conservative residue-class similarity
- Heavy/light chain contribution
- CDR/framework distribution
- Parent VH similarity
- Parent VL similarity
- Parent CDR-length context
- Experimental, literature-derived, or hypothetical status

Do not calculate spatial patch alignment, surface accessibility, molecular shape, RMSD, or contact geometry in the first version.

## 34.3 Paratope evidence profile

Return:

- Input paratope evidence type
- Parent-chain completeness
- Number of residues
- Heavy/light distribution
- CDR/framework composition
- Closest supported residue-set matches
- Numbered-position overlap
- Residue identity and conservative similarity
- Parent-sequence similarity
- Source record
- Antigen metadata as optional context
- Evidence grade
- Evidence sufficiency
- Limitations

Never infer shared antigen specificity from paratope residue-set similarity.

# 35. Heavy–Light Pairing

Separate two outputs.

## 35.1 Observed Pairing Precedent

Use only records where both chains are directly linked in the same source record.

Sources include:

- SAbDab2 paired structures
- Paired OAS records
- Thera-SAbDab pairs
- IEDB BCR records where both chains are linked
- Legally usable paired literature or patent records

Provisional pair-match dimensions:

- Heavy similarity
- Light similarity
- Joint coverage
- VH–VL germline-pair precedent
- Combined CDR-pattern similarity
- Direct parent-record evidence
- Independent paired-record breadth
- Format relevance

## 35.2 Future predicted interface compatibility

Do not implement predicted VH–VL interface compatibility in the first version.

A future structural module may evaluate:

- VH–VL orientation
- Interface residue compatibility
- Steric clashes
- Buried surface
- Conserved interface motifs
- Model confidence

This future evidence must remain separate from directly observed pairing precedent and must never be converted into observed pairing evidence.

# 36. Candidate-Assembly Validation

For every candidate construct validate:

- Component identity
- Domain order
- Chain role
- Orientation
- Linker requirement
- Linker presence
- Constant-domain continuity
- Kappa/lambda consistency
- Isotype/subclass consistency
- Conserved cysteine positions
- Expected disulfide pattern
- Junction integrity
- Duplicate or missing domains
- Unsupported format elements
- Chain stoichiometry
- Copy number
- Sequence completeness

Produce:

- Assembly diagram
- Construct manifest
- Chain sequences
- Domain coordinates
- Junction table
- Validation warnings
- Format precedent
- Unsupported assumptions
- Exportable FASTA
- Exportable JSON
- Reproducible assembly recipe

Do not silently generate missing biological components. Any generated placeholder must be explicitly labeled and excluded from scoring.

---

# 37. Whole-Construct Comparative Evidence Index

The component profile is authoritative.

A whole-construct index may be calculated only when:

- Required components are present.
- Component identities are unambiguous.
- The assembly passes basic architecture validation.
- Minimum evidence sufficiency is met.
- Critical components have valid scores.
- Pairing evidence is reported separately.
- The formula version is calibrated and documented.

Use a reliability-weighted geometric mean so weak critical components are not hidden by strong routine components.

```text
Whole-Construct Comparative Evidence Index =
100 × exp(
  Σ [component_weight_i × reliability_i × ln(max(component_score_i / 100, ε))]
  /
  Σ [component_weight_i × reliability_i]
)
```

Also report:

- Critical-component bottleneck score
- Evidence coverage
- Evidence confidence
- Observed pairing status
- Assembly-validation status
- Missing critical evidence

Provisional component weights for an ordinary paired Fv or IgG-like construct:

| Component | Provisional weight |
|---|---:|
| VH | 0.14 |
| VL | 0.14 |
| CDR-H3 | 0.14 |
| CDR-H1 and H2 combined | 0.10 |
| CDR-L1, L2, and L3 combined | 0.12 |
| Heavy frameworks | 0.08 |
| Light frameworks | 0.08 |
| Observed heavy–light pairing | 0.10 |
| Structure-linked experimental precedent | 0.05 |
| Constant/Fc region | 0.05 |

These weights are provisional and format-specific. Do not use them for VHH, scFv, bispecific, or fusion formats without a separate versioned weight set.

Always display:

> **This index summarizes comparative database evidence. It is not a probability of folding, expression, pairing, binding, stability, manufacturability, safety, efficacy, or experimental success.**

Do not calculate the index when critical evidence is missing. Display “Insufficient evidence for a whole-construct index.”

---

# 38. Per-Match Evidence Record

Every match must create an immutable evidence record.

```json
{
  "evidence_id": "unique_identifier",
  "analysis_id": "analysis_identifier",
  "assembly_id": "candidate_assembly_identifier",
  "input_component_id": "component_identifier",
  "input_component_type": "CDR-H3",
  "matched_component_type": "CDR-H3",
  "component_relationship": "exact_same_type",
  "source_database": "OAS",
  "source_snapshot": "snapshot_identifier",
  "source_accession": "record_identifier",
  "source_parent_antibody_id": "parent_identifier",
  "retrieved_at": "timestamp",
  "license_record_id": "license_identifier",
  "input_sequence_checksum": "sha256",
  "matched_sequence_checksum": "sha256",
  "alignment_method": "method",
  "alignment_parameters": {},
  "sequence_identity": 0.0,
  "conservative_similarity": 0.0,
  "query_coverage": 0.0,
  "subject_coverage": 0.0,
  "alignment_length": 0,
  "insertions": 0,
  "deletions": 0,
  "e_value": null,
  "bit_score": null,
  "calibrated_percentile": null,
  "empirical_p_value": null,
  "adjusted_q_value": null,
  "germline_context": {},
  "species": "",
  "paired_record": false,
  "structure_available": false,
  "therapeutic_record": false,
  "target_known": false,
  "experimental_data_available": false,
  "publication_ids": [],
  "evidence_class": "observed_repertoire",
  "evidence_directness": "direct_sequence_observation",
  "metadata_completeness": 0.0,
  "independence_cluster_id": "cluster_identifier",
  "formula_version": "version",
  "calibration_version": "version",
  "warnings": [],
  "limitations": []
}
```

---

# 39. Score Explanation Schema

Every score must expose:

- Score name
- Score purpose
- Score version
- Applicable component
- Raw evidence records
- Raw measurements
- Normalized values
- Calibration reference set
- Percentile calculation
- Formula
- Weights
- Reliability adjustments
- Deduplication decisions
- Independence corrections
- Missing evidence
- Confidence interval
- Evidence grade
- Evidence sufficiency
- Warnings
- Limitations

Example:

```text
CDR-H3 Observed Repertoire Precedent: 72/100
Evidence grade: B
Evidence sufficiency: Moderate

Raw evidence:
- Best nonredundant same-length match: 86% position-aware similarity
- Empirical similarity percentile: 91st
- Adjusted chance-match q-value: 0.018
- 14 nonredundant neighboring sequences
- 4 independent studies
- 9 independent subjects
- 2 directly paired records
- 1 structurally characterized parent antibody

Excluded:
- 38 clonally related records
- 3 duplicate database representations
- 2 records with ambiguous CDR boundaries

Interpretation:
The CDR-H3 has substantial database precedent, but the evidence does not
establish antigen specificity or function.
```

---

# 40. Required User Interface

## 40.1 Project and assembly workspace

Show:

- Analysis mode
- Candidate assembly graph
- Chain list
- Domain order
- Validation status
- Source of every component
- User-supplied versus database-derived information

## 40.2 Antibody component map

Display:

- Chains
- Variable domains
- CDRs
- Frameworks
- Constant domains
- Fc
- Linkers
- Junctions
- Tags
- Fusion partners
- Paratope residues

Every component must be clickable.

## 40.3 Component evidence dashboard

For the selected component show:

- Annotation confidence
- Raw closest-match measurements
- Calibrated percentile
- Observed repertoire precedent
- Structure-linked experimental precedent
- Therapeutic precedent
- Pairing precedent
- Database-relative novelty
- Evidence grade
- Evidence sufficiency
- Missing evidence
- Limitations

## 40.4 Match viewer

Show:

- Input sequence
- Matched sequence
- Numbered alignment
- Component type
- Identity
- Coverage
- Insertions and deletions
- Germline
- Species
- Paired status
- Study
- Structure-record availability
- Therapeutic status
- Target as optional context
- Publications
- Source snapshot
- Score contribution

## 40.5 Structure-linked evidence panel

When a matched antibody has an experimental structure record, show a metadata panel rather than a three-dimensional viewer.

Display:

- PDB or source identifier
- Experimental method
- Resolution or applicable quality field
- Direct heavy–light pairing status
- Antigen-complex status
- Chain identifiers
- Missing-residue annotations when supplied
- Reported contact residues when directly curated
- Publication
- Source version
- Important sequence differences from the submitted component

Do not render, superimpose, predict, or geometrically compare three-dimensional structures in the first version.

## 40.6 Score explainer

Every displayed score must have a “Why?” view containing:

- Formula
- Version
- Raw data
- Calibration set
- Weights
- Deduplication
- Missing data
- Warnings
- Interpretation
- Prohibited interpretations

---

# 41. Backend Architecture

Use modular services.

Recommended modules:

1. Authentication and project service
2. Secure sequence-ingestion service
3. Input normalization service
4. Antibody annotation service
5. Germline-assignment service
6. Assembly-validation service
7. Component extraction service
8. Sequence search service
9. Short-CDR search service
10. Pair search service
11. Structure-linked metadata enrichment service
12. Paratope residue-set comparison service
13. Evidence aggregation service
14. Deduplication and independence service
15. Score calibration service
16. Score explanation service
17. Source snapshot service
18. Publication-linking service
19. Export and report service
20. Audit and provenance service

Recommended storage separation:

- Transactional project database
- Immutable raw evidence object store
- Normalized relational evidence database
- Sequence indexes
- Structure-linked metadata indexes
- Search cache
- Versioned calibration artifacts
- Audit log
- License and source registry

The scoring service must never query an uncontrolled live website directly.

---

# 42. Confidential Sequence and Privacy Requirements

Candidate antibody sequences may be proprietary.

Required protections:

- Encryption in transit
- Encryption at rest
- Project-level access controls
- Tenant isolation
- Local processing by default
- No external tool submission without explicit user permission
- No training on user sequences
- No sequence content in ordinary application logs
- Configurable retention
- Immediate deletion workflow
- Audit logging
- Exportable provenance
- Private deployment option
- Redacted error messages
- Secret-management controls
- Data-processing inventory
- Third-party transfer record

For every external submission display:

- Destination
- Data transmitted
- Retention policy
- Terms of use
- User consent
- Whether the sequence is confidential
- Whether a local alternative exists

---

# 43. True MVP Scope

The MVP must be narrow enough to validate the scientific foundation.

## 43.1 MVP formats

Support initially:

- Paired ordinary heavy and light variable domains
- Heavy chain only
- Light chain only
- VH only
- VL only
- VHH
- Individual CDR with parent context
- Basic IgG/Fab-like candidate assembly
- Basic scFv candidate assembly

Defer complex bispecific and multivalent architectures.

## 43.2 MVP inputs

- Amino-acid sequence
- Optional nucleotide sequence
- Heavy and light chain
- Individual VH, VL, VHH, or CDR
- User-selected constant region
- User-selected linker for scFv
- Optional PDB, SAbDab2, or publication accession as metadata
- Optional user-defined paratope residue positions attached to a parent sequence
- Optional intended target metadata

## 43.3 MVP functions

1. Secure input ingestion
2. Sequence validation
3. Chain classification
4. IMGT numbering
5. CDR and framework extraction
6. Germline assignment
7. Candidate assembly manifest
8. Basic domain-order and junction validation
9. Heavy-chain comparison
10. Light-chain comparison
11. VH/VL/VHH comparison
12. Six independent CDR comparisons
13. Framework comparison
14. Observed repertoire precedent
15. Structure-linked experimental precedent based on sequence and metadata
16. Therapeutic precedent
17. Direct paired-record search
18. Paratope numbered-residue comparison when supported data exist
19. Database-relative novelty
20. Evidence sufficiency
21. Transparent provisional scoring
22. Component map
23. Alignment viewer
24. Structure-linked metadata panel
25. Evidence provenance
26. Exportable report
27. Whole-construct index only after calibration gates are met

## 43.4 Explicitly defer from the MVP

- PDB or mmCIF coordinate upload and parsing
- Antibody structure prediction
- Structural alignment
- RMSD calculation
- CDR loop geometry comparison
- VH–VL orientation analysis
- Molecular-surface calculation
- Interactive three-dimensional visualization
- Structural superposition
- Coordinate-derived paratope extraction
- Paratope prediction
- Spatial paratope-patch comparison
- Predicted heavy–light interface compatibility
- Antibody–antigen complex modeling
- Binding prediction
- Docking
- Affinity prediction
- Neutralization prediction
- Automatic sequence optimization
- Autonomous sequence generation
- Full developability prediction
- Full Fc functional prediction
- Complex bispecific formats
- Multispecific assembly
- Epitope prediction
- Clinical recommendation
- Laboratory protocol generation

The first version may still use sequence and metadata from experimentally determined structures. It must not perform three-dimensional structural computation.

# 44. Phased Development Plan

## Phase 0: Scientific, technical, and legal feasibility

Complete:

- Component ontology
- Assembly graph
- Canonical numbering policy
- Database-access matrix
- Licensing matrix
- Source snapshot plan
- Privacy threat model
- Evidence schema
- Search architecture
- Deduplication rules
- Calibration plan
- Benchmark plan
- MVP acceptance criteria

No production implementation should begin before high-risk access questions are resolved.

## Phase 1: Annotation and assembly representation

Build:

- Input validation
- Chain identification
- Numbering
- CDR/framework extraction
- Germline assignment
- Assembly graph
- Construct manifest
- Basic format validation

## Phase 2: Local source ingestion

Build:

- OAS snapshot ingestion
- OGRDB ingestion
- SAbDab2 ingestion
- PDB sequence and metadata ingestion
- Thera-SAbDab ingestion
- UniProt ingestion
- Source registry
- Snapshot checksums
- Parser tests

## Phase 3: Component search

Build:

- Exact search
- Long-sequence approximate search
- Pairwise refinement
- CDR-specific search
- Framework positional search
- Pair search
- Constant and Fc search
- Linker precedent search

## Phase 4: Evidence aggregation

Build:

- Evidence records
- Parent-antibody normalization
- Duplicate suppression
- Study and subject independence
- Publication linking
- Conflict handling
- Missing-data statuses

## Phase 5: Calibration and transparent scoring

Build:

- Reference distributions
- Component-specific percentiles
- Short-CDR null models
- Evidence grades
- Evidence sufficiency
- Database-relative novelty
- Score explanations
- Confidence intervals

## Phase 6: Candidate-assembly summary

Build:

- Assembly validation
- Format precedent
- Whole-construct gated index
- Bottleneck score
- Evidence coverage
- Exported candidate sequence and manifest

## Phase 7: Future structural and spatial extension

Only after the sequence-first MVP is validated, consider building:

- PDB/mmCIF coordinate normalization
- Structural search
- Explicit structural alignment
- CDR loop geometry comparison
- VH–VL orientation analysis
- Experimental paratope extraction from coordinates
- Paratope patch and surface comparison
- Predicted structure support
- Predicted interface compatibility
- Interactive three-dimensional viewer

This phase must remain a separate module with separate validation, computational requirements, uncertainty displays, and evidence labels.

## Phase 8: Optional antigen and epitope context

Build:

- IEDB integration
- Antigen normalization
- Epitope records
- Assay context
- Clear separation from antibody scores

---

# 45. Validation and Calibration

Use benchmark sets that prevent data leakage.

Required splits:

- Temporal holdout by database release
- Sequence-identity-cluster holdout
- Germline-family holdout
- Species holdout where appropriate
- Exact-match exclusion
- Paired-antibody holdout
- Therapeutic holdout
- Structure-linked parent-antibody holdout
- CDR-length bins
- Rare-format holdout

Required benchmark cases:

- Exact known heavy chains
- Mutated known heavy chains
- Exact known light chains
- Mutated known light chains
- Known VH and VL
- Known VHH
- Common CDRs
- Rare CDRs
- Short chance-matching CDRs
- Long CDR-H3 sequences
- Known direct pairs
- Artificially mismatched heavy and light chains
- Therapeutic antibodies
- Observed repertoire sequences
- Antibodies with structure-linked experimental records
- Synthetic sequences
- Truncated sequences
- Incorrect chain labels
- Kappa/lambda mislabels
- Invalid assemblies
- Linker and junction errors

Required negative controls:

- Composition-preserving CDR shuffles
- Length-matched random CDRs
- Framework-preserving CDR swaps
- Cross-region comparisons
- Cross-chain comparisons
- Artificial pair mismatches
- Duplicate-heavy datasets
- Sparse-species datasets

Measure:

- Component-classification accuracy
- Numbering accuracy
- Germline-assignment accuracy
- Known-parent retrieval
- Top-k recall
- Ranking quality
- False-match rate
- Chance-match calibration
- Score calibration
- Confidence calibration
- Pairing false-positive rate
- Duplicate suppression
- Score stability across releases
- Reproducibility
- Latency
- Failure rate by component
- Coverage by species and format

---

# 46. Acceptance Criteria

The MVP is acceptable only if:

1. Heavy and light chains are not cross-compared as equivalent components.
2. Every CDR is searched only against the same CDR type and definition.
3. Short-CDR matches use a chance-corrected method.
4. Paired and unpaired evidence remain separate.
5. Every score can be reproduced from stored raw records.
6. Every source has a recorded access and license status.
7. Every analysis is tied to immutable source snapshots.
8. Missing evidence is not converted to zero-quality evidence.
9. Exact duplicate records do not inflate support.
10. A database outage does not lower a biological score.
11. Structure-linked experimental metadata is visibly distinct from sequence similarity and computational inference.
12. A whole-construct index is not calculated when critical evidence is missing.
13. Every whole-construct result displays the prohibited interpretations.
14. Confidential sequences are not sent externally without explicit permission.
15. The same input and source snapshot reproduce the same result.
16. Score changes across software or database versions are explainable.
17. Unsupported formats fail safely and visibly.
18. The exported report includes all formulas, versions, source records, and warnings.
19. The MVP does not calculate structural alignment, RMSD, molecular geometry, or render a three-dimensional model.

---

# 47. Risk Register

The final plan must include mitigation for:

- IMGT tool-access restrictions
- Database terms changing
- Source-schema drift
- Missing or stale therapeutic metadata
- Duplicate records across databases
- Repertoire sampling bias
- Clonal expansion
- Unpaired-chain overinterpretation
- Short-sequence chance matches
- CDR-definition disagreement
- Numbering failure
- Sparse species coverage
- Patent-derived sequence restrictions
- Proprietary sequence exposure
- Structure-record quality variation
- Missing structure metadata
- Incorrect biological-assembly annotations
- Structure-linked metadata being mistaken for proof of structural similarity
- Future predicted-structure overinterpretation
- Heuristic-score overprecision
- Calibration leakage
- Database update score drift
- Whole-score misuse
- User confusion between precedent and function
- Computational cost of billion-scale sequence search
- Insufficient linker and Fc reference coverage

---

# 48. Required Scientific Limitations

The final product must state clearly:

1. Similar antibody sequences can bind different targets.
2. Small CDR changes can change specificity.
3. Similar CDRs do not guarantee similar binding.
4. Paratopes are three-dimensional residue sets biologically, but the first version represents them only as numbered residue lists and cannot compare their spatial geometry.
5. Framework residues can influence CDR conformation.
6. Heavy and light chains cannot always be evaluated independently.
7. Unpaired repertoire data cannot establish natural pairing.
8. A directly observed pair does not guarantee that a modified pair will behave similarly.
9. Repertoire frequencies depend on study design, tissue, immune state, sequencing technology, and processing.
10. Database absence does not prove biological novelty.
11. Therapeutic similarity does not establish safety or efficacy.
12. The existence of an experimental structure for a similar sequence does not prove that the submitted sequence has the same detailed conformation.
13. The first version does not calculate structural similarity, and structure-linked precedent does not establish a shared target.
14. Constant regions, linkers, and formats can alter behavior.
15. Scores summarize comparative evidence and are not probabilities of success.
16. Short sequence matches require empirical chance correction.
17. Public databases have incomplete and uneven metadata.
18. Multiple databases may repeat the same underlying evidence.
19. Proprietary antibody data are not represented.
20. Patent disclosure does not imply freedom to operate.
21. No public database is complete for all Fc variants or linkers.
22. Experimental validation remains necessary.

---

# 49. Explicitly Excluded Claims and Functions

Do not make these central MVP functions:

- Definitive antibody–antigen binding prediction
- Docking-based affinity prediction
- Neutralization prediction
- Clinical-success prediction
- Safety prediction
- Toxicity prediction
- Automatic therapeutic recommendation
- Automatic antibody optimization
- Autonomous sequence generation
- Automated laboratory protocol design
- Diagnostic or treatment decisions
- Patentability or freedom-to-operate conclusions

These may be considered only as separate, independently validated modules.

---

# 50. Required Planning Deliverables

Before implementation, produce:

1. Corrected scientific scope
2. Product terminology
3. Assembly graph and component ontology
4. Supported input types
5. Input schema
6. Validation rules
7. Numbering policy
8. Database-access matrix
9. Licensing and redistribution matrix
10. Component-to-source mapping
11. Exact fields required from every source
12. Source snapshot architecture
13. Data-normalization rules
14. Conflict-resolution rules
15. Deduplication and independence rules
16. Search and indexing architecture
17. Same-component comparison rules
18. Long-sequence measurement method
19. Short-CDR empirical null method
20. Sequence-first paratope residue-set representation and comparison method
21. Structure-linked experimental metadata model
22. Pairing evidence method
23. Assembly-validation method
24. Annotation-confidence model
25. Observed repertoire model
26. Structure-linked experimental precedent model
27. Therapeutic precedent model
28. Database-relative novelty model
29. Evidence sufficiency model
30. Whole-construct gated index
31. Calibration strategy
32. Confidence and uncertainty method
33. Evidence record schema
34. Source snapshot schema
35. Score explanation schema
36. Backend architecture
37. Security and privacy architecture
38. User-interface requirements
39. MVP implementation plan
40. Benchmark plan
41. Acceptance criteria
42. Risk register
43. Scientific limitations
44. Excluded features
45. Estimated infrastructure requirements
46. Maintenance and source-update process

---

# 51. Final Instruction to the Planning Agent

Keep every design decision centered on this workflow:

```text
Submit an existing antibody
or assemble a candidate construct
        ↓
Preserve and validate all inputs
        ↓
Identify exact component types
        ↓
Map components into a candidate assembly graph
        ↓
Number and annotate supported antibody domains
        ↓
Compare each component only with equivalent known components
        ↓
Retrieve versioned natural-repertoire, structure-linked,
therapeutic, pairing, reference, and literature evidence
        ↓
Deduplicate repeated evidence
        ↓
Show raw measurements first
        ↓
Calculate component-specific calibrated percentiles
        ↓
Assign evidence grades and evidence sufficiency
        ↓
Calculate database-relative novelty
        ↓
Validate the proposed assembly
        ↓
Optionally calculate a gated whole-construct comparative index
        ↓
Explain every record, measurement, calibration,
weight, warning, limitation, and missing field
```

The platform must always prioritize the component-level evidence profile over the whole-construct score.

The final output must show:

- What is known
- What is similar
- What is different
- What is directly observed
- What is predicted
- What is missing
- How complete the accessible evidence is
- How uncertain every conclusion is
- Exactly where every result came from

The candidate assembly is a computational representation.

It is not proof of:

- Correct folding
- Correct expression
- Correct pairing
- Antigen binding
- Specificity
- Affinity
- Stability
- Manufacturability
- Safety
- Efficacy
- Experimental success

The finished platform must function as a scientifically cautious, reproducible, explainable first-layer antibody component and candidate-assembly evidence engine.
