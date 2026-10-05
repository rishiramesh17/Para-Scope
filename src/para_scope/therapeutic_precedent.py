"""
Therapeutic precedent (Version 2, Gate 7; SPEC S26).

What this answers
-----------------
Has this exact component, or a similar one, been observed in an antibody that
entered therapeutic development -- and what does the curated record say about it?

Nothing else. In particular this dimension says nothing about whether the
candidate will behave like the therapeutic it resembles. SPEC S26 requires the
disclaimer below to be stated wherever a result is surfaced, and it is a field on
every result rather than a caller's responsibility to remember.

Why the composite score is withheld
-----------------------------------
SPEC S26 defines Therapeutic Precedent as a weighted sum whose first and largest
term (0.45) is *calibrated* component similarity. Calibration is Gate 9; it does
not exist yet. Substituting raw alignment identity into that slot would publish an
uncalibrated number under a calibrated formula's name -- reproducible, but
reproducibly wrong, and indistinguishable downstream from a finished score.

So `score` stays None and `score_status` says why, exactly as an unavailable
measurement does elsewhere. What is emitted instead is every raw measurement the
formula consumes (SPEC S18: raw measurements first) plus the sub-terms that are
computable *now* without calibration, in `computable_terms`. Gate 9 then supplies
the calibrated similarity term and the composite follows, without any of this
having to be recomputed or reinterpreted.

Two of the six terms are blocked, for different reasons:

- *calibrated component similarity* -- needs Gate 9's empirical reference
  distributions. Raw `global_identity` is reported meanwhile, clearly named as raw.
- *format relevance* -- needs the candidate's declared construct format, and
  Thera-SAbDab's `Format` column is free text with 120+ distinct values and no
  controlled vocabulary. Classifying it is separate, explicitly calibrated work
  (ROADMAP, Gate 7 constraints), not something to approximate here. The matched
  therapeutic's format is reported verbatim so a reader can judge relevance
  themselves.

Deduplication
-------------
SPEC S15.5 requires one canonical entity per therapeutic. A therapeutic can match
through several arms, and several components of one therapeutic can match one
query, so matches are grouped and `therapeutics_matched` counts entities while
`raw_match_count` counts observations. Reporting the raw count as though it were
independent precedent is the same error as counting one antibody's twelve crystal
structures as twelve pieces of pairing evidence.

Short components
----------------
CDRs and other short components are excluded from *alignment* search, since
unmodified full-length alignment scoring is invalid for them (SPEC S21, CLAUDE.md)
and no length-matched empirical null exists until Gate 9. They are **not** excluded
from *exact* matching: an identical CDR-H3 in a therapeutic is real, unambiguous
evidence that needs no chance correction, because there is no alignment to score.
So a short component returns exact hits with `approximate_search_status` explaining
what was not attempted -- rather than an empty result that reads as "nothing found".
"""

from dataclasses import dataclass, field
from typing import Optional

from para_scope.corpus import Corpus, TherapeuticRecord
from para_scope.indexes.approximate import (
    ALIGNMENT_SCORING_VERSION,
    ELIGIBLE_TYPES as APPROXIMATE_ELIGIBLE_TYPES,
    MIN_QUERY_COVERAGE,
    MIN_SUBJECT_COVERAGE,
)
from para_scope.models import (
    ComponentRelationship,
    ComponentType,
    EvidenceClass,
    EvidenceStatus,
    compare_eligibility,
)
from para_scope.calibration import CalibratedPercentile, CalibrationSet
from para_scope.numbering import (
    compute_sequence_checksum,
    normalize_query_to_variable_domain,
)

# Bump on any change that alters an emitted measurement or the set of matches
# returned. Stored on every result so a stored profile is never reinterpreted
# under later rules (SPEC S19).
#   1.0.0 -- raw measurements + computable sub-terms; composite withheld pending Gate 9
#   1.1.0 -- each match carries its source database(s) and accessions, so provenance
#            and licence terms travel with displayed evidence (SPEC S13/S14) instead
#            of needing a second lookup against the corpus
THERAPEUTIC_PRECEDENT_VERSION = "1.1.0"

# No calibration exists yet (Gate 9). Named rather than left blank so a stored
# record states which calibration produced it, including when the answer is "none".
CALIBRATION_VERSION = "uncalibrated"

# The source whose records constitute therapeutic precedent. A component reaches
# this dimension only through provenance from this database -- never by resembling
# something that happens to also appear in a structure database.
THERAPEUTIC_SOURCE = "Thera-SAbDab"

# SPEC S26, verbatim and mandatory.
THERAPEUTIC_DISCLAIMER = (
    "Similarity to a therapeutic antibody does not transfer approval, safety, "
    "efficacy, manufacturability, target specificity, or clinical performance."
)

# Related-class matches (VH<->VHH, kappa<->lambda) are reported but ranked below
# same-class ones at equal identity, never merged into same-class evidence
# (CLAUDE.md). This is a rank penalty for ordering only; it does not alter any
# reported measurement.
RELATED_CLASS_RANK_PENALTY = 0.5


@dataclass
class TherapeuticMatch:
    """One therapeutic component matching the query, with its curated record."""

    therapeutic_name: str
    matched_record_id: str
    matched_component_type: ComponentType
    component_relationship: ComponentRelationship
    match_method: str  # "exact_checksum" | "pairwise_local_blosum62"

    # Raw alignment measurements (SPEC S18). An exact match reports 1.0 throughout
    # because it is identical, not because a default was applied.
    sequence_identity: float
    global_identity: float
    query_coverage: float
    subject_coverage: float
    alignment_length: int
    insertions: int = 0
    deletions: int = 0

    # SPEC S19 Layer 2: where this match's global_identity sits in the reference
    # distribution for its component class. None when no calibration was supplied;
    # a CalibratedPercentile carrying INSUFFICIENT_EVIDENCE when one was supplied
    # but holds no distribution for this stratum. Those are different facts and are
    # not collapsed.
    calibrated_similarity: Optional[CalibratedPercentile] = None

    # Curated therapeutic record (SPEC S26 display fields).
    alternative_names: list[str] = field(default_factory=list)
    target: Optional[str] = None
    format_raw: Optional[str] = None
    ch1_isotype: Optional[str] = None
    # Metadata only. An antibody reaching Phase III transfers nothing to a sequence
    # resembling it; neither field may be read as a success or failure label.
    highest_clinical_trial: Optional[str] = None
    development_status: Optional[str] = None

    # Supporting-record status. Optional[bool] throughout: None means the source did
    # not let us determine this, which is not False.
    arm_accessions: list[str] = field(default_factory=list)
    structure_ids: list[str] = field(default_factory=list)

    # Which database(s) this matched component was actually observed in, carried on
    # the match rather than left to be re-derived from matched_record_id. The
    # corpus has always held it; a consumer had to take the record id back to
    # corpus.get(...).source_pointers to reach it, which is a hop that is easy to
    # get subtly wrong -- the source name is "Thera-SAbDab", and filtering on the
    # plausible-looking "thera_sabdab" silently yields zero hits rather than an
    # error.
    #
    # It also matters for licensing, not just convenience: SPEC S13/S14 require a
    # source's licence to travel with the evidence derived from it. Licence is
    # recorded per SourceSnapshot, so a displayed match that does not name its
    # source cannot be joined to its terms. Both sources are currently
    # `unspecified`; the moment one is cleared and the other is not, that join
    # stops being cosmetic.
    source_databases: list[str] = field(default_factory=list)
    source_pointers: list[tuple[str, str]] = field(default_factory=list)
    has_structure_record: Optional[bool] = None
    directly_paired: Optional[bool] = None
    metadata_completeness: float = 0.0
    shortlisted: bool = False

    @property
    def rank_key(self) -> tuple[float, float]:
        """Ordering key: related-class matches sit below same-class at equal identity."""
        weight = (
            1.0
            if self.component_relationship is ComponentRelationship.EXACT_SAME_TYPE
            else RELATED_CLASS_RANK_PENALTY
        )
        return (self.global_identity * weight, self.query_coverage)


@dataclass
class TherapeuticPrecedentResult:
    """Therapeutic-precedent evidence for one component (SPEC S26)."""

    status: EvidenceStatus
    component_type: ComponentType
    query_length: int

    matches: list[TherapeuticMatch] = field(default_factory=list)
    # SPEC S15: entities, then observations. These are different numbers and both
    # are reported, so neither can be mistaken for the other.
    therapeutics_matched: int = 0
    raw_match_count: int = 0

    best_global_identity: Optional[float] = None
    exact_match_count: int = 0

    # The composite defined by SPEC S26. Withheld -- see the module docstring.
    score: Optional[float] = None
    score_status: EvidenceStatus = EvidenceStatus.INSUFFICIENT_EVIDENCE
    score_reason: str = ""
    computable_terms: dict[str, Optional[float]] = field(default_factory=dict)
    blocked_terms: dict[str, str] = field(default_factory=dict)

    approximate_search_status: EvidenceStatus = EvidenceStatus.AVAILABLE
    approximate_search_reason: Optional[str] = None

    evidence_class: EvidenceClass = EvidenceClass.THERAPEUTIC_PRECEDENT
    engine_version: str = THERAPEUTIC_PRECEDENT_VERSION
    scoring_version: str = ALIGNMENT_SCORING_VERSION
    calibration_version: str = CALIBRATION_VERSION
    reason: Optional[str] = None
    limitations: list[str] = field(default_factory=list)
    disclaimer: str = THERAPEUTIC_DISCLAIMER


def _arms_containing(therapeutic: TherapeuticRecord, record_id: str) -> list:
    """The therapeutic's arms that contributed this exact component."""
    return [
        arm for arm in therapeutic.arms
        if record_id in (arm.heavy_record_id, arm.light_record_id)
    ]


def _build_match(
    corpus: Corpus,
    therapeutic: TherapeuticRecord,
    record_id: str,
    matched_type: ComponentType,
    relationship: ComponentRelationship,
    match_method: str,
    sequence_identity: float,
    global_identity: float,
    query_coverage: float,
    subject_coverage: float,
    alignment_length: int,
    insertions: int = 0,
    deletions: int = 0,
    shortlisted: bool = False,
) -> TherapeuticMatch:
    arms = _arms_containing(therapeutic, record_id)

    # Pairing is per arm and requires both chains from the same source record. An
    # arm with only a heavy chain is not "unpaired evidence to be joined later";
    # it simply has no observed pair.
    if arms:
        directly_paired = any(
            arm.heavy_record_id is not None and arm.light_record_id is not None
            for arm in arms
        )
    else:
        directly_paired = None

    structure_ids: list[str] = []
    per_arm_structure = [arm.has_structure_record for arm in arms]
    for arm in arms:
        for ref in arm.structure_ids:
            if ref not in structure_ids:
                structure_ids.append(ref)
    if not per_arm_structure or all(value is None for value in per_arm_structure):
        has_structure = None
    else:
        has_structure = any(per_arm_structure)

    # Verbatim from the component record -- never inferred from the fact that a
    # therapeutic record exists. A component can be present in SAbDab2 and reached
    # here through a Thera-SAbDab therapeutic link; both pointers are real and both
    # are reported.
    matched_record = corpus.get(record_id)
    source_pointers = list(matched_record.source_pointers) if matched_record else []
    source_databases = sorted({database for database, _accession in source_pointers})

    return TherapeuticMatch(
        therapeutic_name=therapeutic.therapeutic_name,
        matched_record_id=record_id,
        matched_component_type=matched_type,
        component_relationship=relationship,
        match_method=match_method,
        sequence_identity=sequence_identity,
        global_identity=global_identity,
        query_coverage=query_coverage,
        subject_coverage=subject_coverage,
        alignment_length=alignment_length,
        insertions=insertions,
        deletions=deletions,
        alternative_names=list(therapeutic.alternative_names),
        target=therapeutic.target,
        format_raw=therapeutic.format_raw,
        ch1_isotype=therapeutic.ch1_isotype,
        highest_clinical_trial=therapeutic.highest_clinical_trial,
        development_status=therapeutic.development_status,
        arm_accessions=[arm.accession for arm in arms],
        structure_ids=structure_ids,
        source_databases=source_databases,
        source_pointers=source_pointers,
        has_structure_record=has_structure,
        directly_paired=directly_paired,
        metadata_completeness=therapeutic.metadata_completeness,
        shortlisted=shortlisted,
    )


def _exact_matches(
    corpus: Corpus, component_type: ComponentType, query: str
) -> list[TherapeuticMatch]:
    """Therapeutic components whose sequence is byte-identical to the query.

    The exact index is keyed by sequence checksum alone, so a lookup can return
    components of a *different* type that happen to share a sequence -- a short
    framework or CDR is quite capable of appearing in both a heavy and a light
    chain. Every hit is therefore re-checked with compare_eligibility, and
    excluded pairings are dropped rather than reported at identity 1.0.
    """
    matches: list[TherapeuticMatch] = []
    for record_id in corpus.exact.lookup(compute_sequence_checksum(query)):
        record = corpus.get(record_id)
        if record is None:
            continue
        relationship = compare_eligibility(component_type, record.component_type)
        if relationship is ComponentRelationship.EXCLUDED:
            continue
        for therapeutic in corpus.therapeutics_for_component(record_id):
            matches.append(
                _build_match(
                    corpus=corpus,
                    therapeutic=therapeutic,
                    record_id=record_id,
                    matched_type=record.component_type,
                    relationship=relationship,
                    match_method="exact_checksum",
                    sequence_identity=1.0,
                    global_identity=1.0,
                    query_coverage=1.0,
                    subject_coverage=1.0,
                    alignment_length=len(query),
                )
            )
    return matches


def _approximate_matches(
    corpus: Corpus,
    component_type: ComponentType,
    query: str,
    top_k: int,
    shortlist_limit: Optional[int],
    min_query_coverage: float,
    min_subject_coverage: float,
) -> list[TherapeuticMatch]:
    """Aligned therapeutic components above the coverage floors."""
    matches: list[TherapeuticMatch] = []
    hits = corpus.approximate.search(
        component_type,
        query,
        top_k=top_k,
        min_query_coverage=min_query_coverage,
        min_subject_coverage=min_subject_coverage,
        shortlist_limit=shortlist_limit,
    )
    for hit in hits:
        record = corpus.get(hit.record_id)
        if record is None:
            continue
        therapeutics = corpus.therapeutics_for_component(hit.record_id)
        if not therapeutics:
            continue  # a structure-only component; not therapeutic precedent
        relationship = compare_eligibility(component_type, record.component_type)
        if relationship is ComponentRelationship.EXCLUDED:
            continue
        for therapeutic in therapeutics:
            matches.append(
                _build_match(
                    corpus=corpus,
                    therapeutic=therapeutic,
                    record_id=hit.record_id,
                    matched_type=record.component_type,
                    relationship=relationship,
                    match_method=hit.alignment_method,
                    sequence_identity=hit.sequence_identity,
                    global_identity=hit.global_identity,
                    query_coverage=hit.query_coverage,
                    subject_coverage=hit.subject_coverage,
                    alignment_length=hit.alignment_length,
                    insertions=hit.insertions,
                    deletions=hit.deletions,
                    shortlisted=hit.shortlisted,
                )
            )
    return matches


def _computable_terms(matches: list[TherapeuticMatch]) -> dict[str, Optional[float]]:
    """The SPEC S26 sub-terms that need no calibration, taken from the best match.

    Each is a raw proportion in [0, 1], not a weighted contribution: applying the
    weights would imply the composite exists. `None` marks a term whose inputs were
    themselves unknown, and must not be read as zero.
    """
    if not matches:
        return {}
    best = max(matches, key=lambda m: m.rank_key)

    if best.directly_paired is None and best.has_structure_record is None:
        pair_or_structure: Optional[float] = None
    else:
        pair_or_structure = float(bool(best.directly_paired) or bool(best.has_structure_record))

    return {
        "coverage": best.query_coverage,
        "exact_component_class_confirmation": float(
            best.component_relationship is ComponentRelationship.EXACT_SAME_TYPE
        ),
        "direct_pair_or_structure_support": pair_or_structure,
        "metadata_completeness": best.metadata_completeness,
    }


# Measurement id of the reference distribution that calibrates a match's
# similarity. Named here so the producer (scripts/build_calibration.py) and the
# consumer cannot drift apart silently -- a typo would look like missing data.
SIMILARITY_MEASUREMENT = "similarity.best_global_identity_vs_therapeutics"

BLOCKED_TERMS = {
    "calibrated_component_similarity": (
        "requires the empirical reference distributions built in Gate 9 (SPEC S19 "
        "layers 2-4); raw global_identity is reported on each match meanwhile"
    ),
    "format_relevance": (
        "requires the candidate's declared construct format and a controlled "
        "vocabulary over Thera-SAbDab's free-text Format column (120+ distinct "
        "values); each match reports its format verbatim instead. NOT Gate 9 work "
        "and unowned: the declared format is an input the submission schema does "
        "not collect, and the vocabulary is curation, so neither falls out of "
        "building reference distributions. Gate 9 therefore does not unblock the "
        "composite by itself -- 0.10 of SPEC S26's weight remains missing"
    ),
}

SCORE_WITHHELD_REASON = (
    "SPEC S26's composite weights calibrated component similarity at 0.45 and "
    "format relevance at 0.10. Neither term exists yet, so no composite is emitted; "
    "raw measurements and the four computable sub-terms are reported instead."
)


def compute_therapeutic_precedent(
    corpus: Corpus,
    component_type: ComponentType,
    query: str,
    top_k: int = 10,
    shortlist_limit: Optional[int] = None,
    min_query_coverage: float = MIN_QUERY_COVERAGE,
    min_subject_coverage: float = MIN_SUBJECT_COVERAGE,
    calibration: Optional[CalibrationSet] = None,
) -> TherapeuticPrecedentResult:
    """
    Gather therapeutic-precedent evidence for one component (SPEC S26).

    `calibration`, when supplied, adds SPEC S19 Layer 2 percentiles to each match
    and clears the `calibrated_component_similarity` blocked term. It is optional
    rather than required because a raw measurement is still valid evidence without
    one -- what is not valid is presenting an uncalibrated number as if it were
    calibrated, which is why the blocked term stays listed until a distribution
    actually supplies the percentile.

    Never raises on an unscoreable input: an evidence dimension must be able to
    fill its slot in a profile with an explicit status rather than abort the
    surrounding analysis. `status` distinguishes "we looked and found nothing"
    (NOT_FOUND) from "this could not be evaluated" (UNSUPPORTED_COMPONENT) from
    "the corpus holds no therapeutic records at all" (DATABASE_UNAVAILABLE).
    """
    original_query = query
    query, normalization_note = normalize_query_to_variable_domain(
        component_type.value, query
    )
    result = TherapeuticPrecedentResult(
        status=EvidenceStatus.NOT_FOUND,
        component_type=component_type,
        query_length=len(query),
        blocked_terms=dict(BLOCKED_TERMS),
        score_reason=SCORE_WITHHELD_REASON,
    )

    if normalization_note is not None:
        result.limitations.append(normalization_note)

    if not query:
        result.status = EvidenceStatus.UNSUPPORTED_COMPONENT
        result.reason = "empty sequence"
        return result

    if not corpus.therapeutics:
        result.status = EvidenceStatus.DATABASE_UNAVAILABLE
        result.reason = (
            f"this corpus contains no {THERAPEUTIC_SOURCE} records; therapeutic "
            f"precedent cannot be evaluated. This is not evidence of absence."
        )
        result.approximate_search_status = EvidenceStatus.DATABASE_UNAVAILABLE
        return result

    matches = _exact_matches(corpus, component_type, query)
    result.exact_match_count = len(matches)

    if component_type in APPROXIMATE_ELIGIBLE_TYPES:
        try:
            approximate = _approximate_matches(
                corpus, component_type, query, top_k, shortlist_limit,
                min_query_coverage, min_subject_coverage,
            )
        except Exception as exc:  # index refuses unalignable content, etc.
            result.approximate_search_status = EvidenceStatus.ANALYSIS_FAILED
            result.approximate_search_reason = f"{type(exc).__name__}: {exc}"
            approximate = []
        else:
            already = {m.matched_record_id for m in matches}
            approximate = [m for m in approximate if m.matched_record_id not in already]
        matches.extend(approximate)
    else:
        result.approximate_search_status = EvidenceStatus.UNSUPPORTED_COMPONENT
        result.approximate_search_reason = (
            f"{component_type.value} is too short for unmodified alignment scoring; "
            f"it needs a length-matched empirical null model (SPEC S21), which "
            f"Gate 9 supplies. Exact matches are still reported."
        )

    if not matches:
        result.status = EvidenceStatus.NOT_FOUND
        result.reason = (
            f"no {THERAPEUTIC_SOURCE} component of type {component_type.value} matched "
            f"this sequence above the coverage floors"
        )
        return result

    # Name breaks ties ascending while identity ranks descending, so the order is
    # fully determined by the data rather than by index insertion order. Trastuzumab's
    # VH is byte-identical across 11 distinct therapeutics, so ties are the normal
    # case here, not an edge case -- without a tiebreak the "top match" could differ
    # between two runs over the same snapshot.
    matches.sort(key=lambda m: m.therapeutic_name)
    matches.sort(key=lambda m: m.rank_key, reverse=True)

    # SPEC S19 Layer 2. Attached per match rather than once per result because
    # each match may be a different component class (a related-class hit is
    # percentiled against its own population, never the query's).
    if calibration is not None:
        for match in matches:
            match.calibrated_similarity = calibration.percentile(
                SIMILARITY_MEASUREMENT,
                match.matched_component_type.value,
                match.global_identity,
                # Length of the *matched* component, since that is the population
                # it is being placed in. Length-stratified classes (the CDRs)
                # require it; for the rest it is ignored.
                length=(
                    len(matched.sequence)
                    if (matched := corpus.get(match.matched_record_id)) is not None
                    else None
                ),
            )
        # The term is only unblocked where a percentile actually came back. A
        # supplied calibration that holds no distribution for this stratum leaves
        # the term blocked, with the percentile's own reason explaining why.
        if any(
            m.calibrated_similarity is not None
            and m.calibrated_similarity.status is EvidenceStatus.AVAILABLE
            for m in matches
        ):
            result.blocked_terms = {
                k: v for k, v in result.blocked_terms.items()
                if k != "calibrated_component_similarity"
            }

    result.matches = matches
    result.raw_match_count = len(matches)
    result.therapeutics_matched = len({m.therapeutic_name for m in matches})
    result.best_global_identity = max(m.global_identity for m in matches)
    result.computable_terms = _computable_terms(matches)
    result.status = EvidenceStatus.AVAILABLE

    if any(m.shortlisted for m in matches):
        result.limitations.append(
            "Shortlisted search: the hit set is approximately 96% complete, not "
            "exhaustive. Do not report 'no closer therapeutic exists' from it."
        )
    if any(m.component_relationship is ComponentRelationship.RELATED_CLASS for m in matches):
        result.limitations.append(
            "Contains related-class matches (VH/VHH or kappa/lambda). They are "
            "labelled per match and rank below same-class evidence; they are not "
            "same-class evidence."
        )
    if result.approximate_search_status is not EvidenceStatus.AVAILABLE:
        result.limitations.append(
            f"Alignment search not performed: {result.approximate_search_reason}"
        )

    return result
