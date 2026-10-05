"""
Human germline V-gene assignment (Version 2, Gate 9 milestone 1).

What this measures
-------------------
Which human germline V gene a submitted VH, VL-kappa, or VL-lambda sequence
aligns most closely to, and the identity of that alignment. Nothing more --
no D or J gene call, no CDR3 junction analysis, no allele-frequency claim.

This exists because Database-Relative Novelty is close to meaningless without
a baseline: "unusual" relative to what? `EvidenceClass.GERMLINE` has been
declared in models.py since early in this project and nothing populated it,
because no germline reference set or alignment path existed. This module is
that path.

Reduced scope, deliberately (docs/ROADMAP.md, "Germline reference data")
--------------------------------------------------------------------------
- **Human only.** VHH (camelid) is UNSUPPORTED_COMPONENT here, not a lower
  score -- there is no camelid germline set in this reference. Extending to
  VHH is real follow-up work, not an oversight.
- **V gene only, no D/J.** A V-gene call plus percent identity to nearest
  germline is the standard measure of drift from the inherited starting
  point, and captures most of Novelty's value; full V(D)J gene calls are a
  separate, larger effort.
- **One reference set.** SPEC S11.2 forbids silently merging assignments from
  different germline databases (e.g. OGRDB and IMGT/GENE-DB). That concern
  does not arise yet because only one set is used. If a second set is ever
  added, both assignments must be reported side by side, never merged.

Reference set and its license
------------------------------
OGRDB (https://ogrdb.airr-community.org), Homo sapiens IGH_VDJ /
IGKappa_VJ / IGLambda_VJ published germline sets, V genes only. Published
under CC0 (public domain) -- confirmed against
https://ogrdb.airr-community.org/render_page/licensing_statement.html.
Unlike SAbDab2 and Thera-SAbDab (both `license_id = "unspecified"`), CC0
permits redistribution, so this reference set is bundled inside the package
(`reference_data/germline/`) rather than gitignored and fetched separately --
the same pattern `promb` uses for its own reference data. Built by
`scripts/build_germline_reference.py`; nucleotide V-gene sequences translated
in frame 1 (the IMGT V-REGION convention: no leader, sequence starts on a
codon), pseudogenes and ambiguous calls dropped at build time.

Reuses the existing alignment engine
-------------------------------------
No new alignment code. A germline V gene is just a sequence of the same
`ComponentType` as the query it will be compared against (heavy references
indexed as VH, kappa as VL-kappa, lambda as VL-lambda), so the existing
`ApproximateIndex` (SPEC S16.2) does the alignment, coverage admissibility,
and ranking unchanged -- the same engine Gates 7 and 8 already use.

Never a neutral default
------------------------
When no germline reference clears the coverage floor, or the reference set
cannot be loaded, `nearest_germline_identity` stays None and the status says
why. It never falls back to 0.0 or to a corpus mean.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from para_scope.features import STANDARD_AMINO_ACIDS
from para_scope.indexes.approximate import (
    ALIGNMENT_SCORING_VERSION,
    MIN_QUERY_COVERAGE,
    MIN_SUBJECT_COVERAGE,
    ApproximateIndex,
)
from para_scope.models import ComponentType, EvidenceClass, EvidenceStatus
from para_scope.numbering import normalize_query_to_variable_domain

# Bump on any change that alters a computed assignment, including a change of
# reference set or translation rule. Stored on every result so a stored
# assignment is never reinterpreted under later rules (SPEC S19).
#   1.0.0 -- initial: human V-gene identity via the existing alignment engine
GERMLINE_VERSION = "1.0.0"

# No calibration exists yet (Gate 9 proper -- this module is a prerequisite
# for it). Named rather than left blank so a stored record states which
# calibration produced it, including when the answer is "none".
CALIBRATION_VERSION = "uncalibrated"

REFERENCE_SET_ID = "ogrdb-human-igh-igk-igl-v"
REFERENCE_SET_VERSION = "2026-05-27"  # OGRDB "Released" date for the sets used
REFERENCE_SET_DESCRIPTION = (
    "Human IGH_VDJ, IGKappa_VJ, and IGLambda_VJ published germline sets from "
    "OGRDB (https://ogrdb.airr-community.org), V genes only, translated in "
    "frame 1. CC0."
)
_REFERENCE_FILE = (
    Path(__file__).parent / "reference_data" / "germline" / "ogrdb_human_v_genes.fasta"
)

# Only variable-domain types have a germline V gene to align against. VHH is
# explicitly excluded (no camelid reference set), not silently skipped.
SUPPORTED_TYPES = frozenset(
    {ComponentType.VH, ComponentType.VL_KAPPA, ComponentType.VL_LAMBDA}
)
_LOCUS_BY_TYPE = {
    ComponentType.VH: "IGH",
    ComponentType.VL_KAPPA: "IGK",
    ComponentType.VL_LAMBDA: "IGL",
}

GERMLINE_DISCLAIMER = (
    "Germline assignment reports the human V gene this sequence aligns most "
    "closely to, and the identity of that alignment. It is not a measure of "
    "humanness, immunogenicity, or safety, and it is not a prediction that "
    "the candidate derives from, or will behave like, that germline gene. No "
    "D or J gene is called."
)


@dataclass
class GermlineMatch:
    """One candidate germline V gene and the alignment that produced it."""

    allele: str  # e.g. "IGHV3-23*01"
    gene_call: str  # e.g. "IGHV3-23" -- allele stripped
    locus: str  # "IGH" / "IGK" / "IGL"

    global_identity: float
    query_coverage: float
    subject_coverage: float
    alignment_length: int
    insertions: int
    deletions: int
    scoring_version: str = ALIGNMENT_SCORING_VERSION


@dataclass
class GermlineAssignmentResult:
    """Germline assignment for one component, plus everything needed to
    reproduce or discard it."""

    status: EvidenceStatus
    matches: list[GermlineMatch] = field(default_factory=list)

    # Convenience accessors for the single best match; always consistent with
    # matches[0] when status is AVAILABLE, always None otherwise.
    nearest_germline_identity: Optional[float] = None
    nearest_gene_call: Optional[str] = None

    evidence_class: EvidenceClass = EvidenceClass.GERMLINE
    engine_version: str = GERMLINE_VERSION
    scoring_version: str = ALIGNMENT_SCORING_VERSION
    calibration_version: str = CALIBRATION_VERSION

    reference_set_id: str = REFERENCE_SET_ID
    reference_set_version: str = REFERENCE_SET_VERSION
    reference_set_description: str = REFERENCE_SET_DESCRIPTION

    reason: Optional[str] = None
    limitations: list[str] = field(
        default_factory=lambda: [
            "V gene only; no D or J gene is called.",
            "Human reference only; not applicable to VHH/camelid sequences.",
            "Allele-level call is uncalibrated -- treat the gene-family call "
            "as more reliable than the specific allele.",
        ]
    )
    disclaimer: str = GERMLINE_DISCLAIMER


class ReferenceSetUnavailable(RuntimeError):
    """Raised internally when the reference set cannot be loaded; surfaced as a status."""


# Built once and reused -- ~350 sequences, cheap to build, but no reason to
# rebuild per call. Cached as a module-level singleton, mirroring
# window_coverage.py's load_reference_set().
_reference_index: Optional[ApproximateIndex] = None
_load_error: Optional[str] = None


def load_reference_index(force_reload: bool = False) -> ApproximateIndex:
    """
    Load and memoise the germline reference index.

    Raises ReferenceSetUnavailable if the bundled FASTA is missing or
    unreadable. Callers that want a status rather than an exception should
    use `compute_germline_assignment`, which converts this into
    DATABASE_UNAVAILABLE.
    """
    global _reference_index, _load_error

    if force_reload:
        _reference_index, _load_error = None, None

    if _reference_index is not None:
        return _reference_index
    if _load_error is not None:
        raise ReferenceSetUnavailable(_load_error)

    try:
        index = ApproximateIndex()
        for allele, locus, sequence in _parse_reference_fasta(_REFERENCE_FILE):
            component_type = {
                "IGH": ComponentType.VH,
                "IGK": ComponentType.VL_KAPPA,
                "IGL": ComponentType.VL_LAMBDA,
            }[locus]
            index.add(component_type, f"{allele}|{locus}", sequence)
        _reference_index = index
    except Exception as exc:
        _load_error = f"{type(exc).__name__}: {exc}"
        raise ReferenceSetUnavailable(_load_error) from exc

    return _reference_index


def _parse_reference_fasta(path: Path) -> list[tuple[str, str, str]]:
    """Parse the bundled germline FASTA into (allele, locus, sequence) tuples."""
    records: list[tuple[str, str, str]] = []
    header: Optional[str] = None
    seq_lines: list[str] = []

    def flush() -> None:
        if header is None:
            return
        allele, _, locus = header.rpartition("|")
        records.append((allele, locus, "".join(seq_lines)))

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        flush()
    return records


def compute_germline_assignment(
    component_type: ComponentType,
    sequence: str,
    top_k: int = 1,
    min_query_coverage: float = MIN_QUERY_COVERAGE,
    min_subject_coverage: float = MIN_SUBJECT_COVERAGE,
    reference_index: Optional[ApproximateIndex] = None,
) -> GermlineAssignmentResult:
    """
    Compute human germline V-gene assignment for one component sequence.

    Always returns a result. Every failure mode is expressed as an
    EvidenceStatus with a reason rather than as an exception or a substitute
    number, so an evidence profile can record what happened instead of
    losing the slot.

    Statuses returned:
      AVAILABLE              at least one reference gene cleared the coverage floor
      NOT_FOUND               no reference gene cleared the coverage floor
      UNSUPPORTED_COMPONENT   component type has no germline reference (VHH,
                               CDRs, frameworks, full chains), non-standard
                               residues present, or empty input
      DATABASE_UNAVAILABLE    reference set could not be loaded
      ANALYSIS_FAILED         alignment raised unexpectedly
    """

    def result(status: EvidenceStatus, **kwargs) -> GermlineAssignmentResult:
        return GermlineAssignmentResult(status=status, **kwargs)

    # Germline V genes are ~98 residues. A query submitted as a whole chain could
    # never clear the coverage floor against one -- measured at ~70% NOT_FOUND
    # before this normalization, which looked like missing germline evidence but
    # was a length mismatch (docs/ROADMAP.md, "Germline reference data").
    cleaned, normalization_note = normalize_query_to_variable_domain(
        component_type.value, sequence
    )

    if not cleaned:
        return result(
            EvidenceStatus.UNSUPPORTED_COMPONENT,
            reason="Empty sequence; there is nothing to align.",
        )

    if component_type not in SUPPORTED_TYPES:
        return result(
            EvidenceStatus.UNSUPPORTED_COMPONENT,
            reason=(
                f"{component_type.value} has no germline V-gene reference in "
                f"this scope. Only VH, VL-kappa, and VL-lambda are supported "
                f"-- VHH is excluded because no camelid germline set is "
                f"bundled, and CDRs/frameworks/full chains are not compared "
                f"against a whole V gene."
            ),
        )

    unknown = sorted(set(cleaned) - set(STANDARD_AMINO_ACIDS))
    if unknown:
        return result(
            EvidenceStatus.UNSUPPORTED_COMPONENT,
            reason=(
                f"Sequence contains residues outside the 20 standard amino "
                f"acids: {unknown}. The sequence is not modified to make it "
                f"alignable; the caller must decide how to handle the "
                f"residue and record that."
            ),
        )

    try:
        index = reference_index if reference_index is not None else load_reference_index()
    except ReferenceSetUnavailable as exc:
        return result(
            EvidenceStatus.DATABASE_UNAVAILABLE,
            reason=(
                f"Germline reference set could not be loaded ({exc}). "
                f"Assignment is unset rather than defaulted; this status "
                f"must not be read as an unusual or absent germline match."
            ),
        )

    try:
        hits = index.search(
            component_type,
            cleaned,
            top_k=top_k,
            min_query_coverage=min_query_coverage,
            min_subject_coverage=min_subject_coverage,
        )
    except Exception as exc:
        return result(
            EvidenceStatus.ANALYSIS_FAILED,
            reason=(
                f"Germline alignment failed ({type(exc).__name__}: {exc}). "
                f"This is a tooling failure and carries no biological meaning."
            ),
        )

    if not hits:
        return result(
            EvidenceStatus.NOT_FOUND,
            reason=(
                f"No germline V gene in the {_LOCUS_BY_TYPE[component_type]} "
                f"reference set cleared the coverage floor "
                f"(min_query_coverage={min_query_coverage}, "
                f"min_subject_coverage={min_subject_coverage})."
            ),
        )

    matches = []
    for hit in hits:
        allele, _, locus = hit.record_id.rpartition("|")
        gene_call = allele.split("*")[0]
        matches.append(
            GermlineMatch(
                allele=allele,
                gene_call=gene_call,
                locus=locus,
                global_identity=hit.global_identity,
                query_coverage=hit.query_coverage,
                subject_coverage=hit.subject_coverage,
                alignment_length=hit.alignment_length,
                insertions=hit.insertions,
                deletions=hit.deletions,
                scoring_version=hit.scoring_version,
            )
        )

    assignment = GermlineAssignmentResult(
        status=EvidenceStatus.AVAILABLE,
        matches=matches,
        nearest_germline_identity=matches[0].global_identity,
        nearest_gene_call=matches[0].gene_call,
    )
    if normalization_note is not None:
        assignment.limitations.append(normalization_note)
    return assignment
