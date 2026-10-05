"""
Observed heavy-light pairing precedent (SPEC S35.1, Gate 8).

Answers: *has this heavy chain been observed paired with this light chain, in a
single source record?* Nothing about whether the candidate's chains will pair,
express, fold, or remain stable together.

The one rule this module exists to enforce
------------------------------------------
A pair is only a pair if **one source record links both chains**. The heavy hits
and the light hits are searched independently, but a result is emitted only when
the pair index -- built exclusively from same-record co-occurrence -- returns a
record for that exact `(heavy_checksum, light_checksum)` key.

That ordering matters. The obvious implementation is to take the best heavy hit
and the best light hit and report them as a pairing result; that would be a
synthesized pair, which CLAUDE.md prohibits outright. Here the index is the only
thing that can license a pair, so no combination of independent hits can produce
one. `_paired_candidates` is where this is enforced and it is the load-bearing
function in the module.

Two counts, never one (SPEC S15)
--------------------------------
`pairs_matched` counts unique `(heavy, light)` keys; `raw_observation_count`
counts the source records backing them. On the real corpus these are 5,686 and
12,032 -- a ratio of ~2.1, because one antibody appears across many structures
and both sources. Reporting the larger number as "pairs" would count a single
antibody's dozen crystal structures as a dozen independent pairing observations.
The extra records are real corroboration and are reported as such, under their
own name.

No composite is defined here
----------------------------
Unlike SPEC S26, **S35 specifies no weights** -- it lists eight "provisional
pair-match dimensions" without a formula. So this module emits raw measurements
and does not withhold a composite, because there is no composite defined to
withhold. Inventing weights to produce one would be fabricating the spec, and a
number assembled that way would be indistinguishable downstream from a specified
one. Three of the eight dimensions are blocked on work that does not exist yet
and are named in BLOCKED_DIMENSIONS.

Not implemented, deliberately (SPEC S35.2)
------------------------------------------
Predicted VH-VL interface compatibility -- orientation, interface residues,
steric clashes, buried surface -- is explicitly out of scope for the first
version, and the spec requires that such evidence, if it ever exists, stay
separate from observed pairing rather than being converted into it.
"""

from dataclasses import dataclass, field
from typing import Optional

from para_scope.corpus import Corpus
from para_scope.indexes.approximate import (
    ALIGNMENT_SCORING_VERSION,
    MIN_QUERY_COVERAGE,
    MIN_SUBJECT_COVERAGE,
)
from para_scope.indexes.pair import PairRecord
from para_scope.models import (
    ComponentRelationship,
    ComponentType,
    EvidenceClass,
    EvidenceStatus,
    compare_eligibility,
)
from para_scope.numbering import (
    compute_sequence_checksum,
    normalize_query_to_variable_domain,
)

# Bump on any change that alters an emitted measurement or the set of pairs
# returned. Stored on every result so a stored profile is never reinterpreted
# under later rules (SPEC S19).
#   1.0.0 -- observed pairing precedent; five computable dimensions, three blocked
#   1.1.0 -- curated record content (names, targets, status, PDB entry, method,
#            resolution) promoted from PairRecord.metadata onto each match, so a
#            consumer never has to look up dict keys to say what a pair *is*
PAIRING_PRECEDENT_VERSION = "1.1.0"

# No calibration exists yet (Gate 9). Named rather than left blank so a stored
# record states which calibration produced it, including when the answer is "none".
CALIBRATION_VERSION = "uncalibrated"

HEAVY_TYPES = frozenset({ComponentType.VH, ComponentType.VHH})
LIGHT_TYPES = frozenset({ComponentType.VL_KAPPA, ComponentType.VL_LAMBDA})

# Related-class pairs (VH<->VHH, kappa<->lambda) are reported but ranked below
# same-class ones, never merged into same-class evidence (CLAUDE.md). Ordering
# only -- it alters no reported measurement.
RELATED_CLASS_RANK_PENALTY = 0.5

PAIRING_DISCLAIMER = (
    "Observed pairing precedent reports heavy-light chain pairs recorded together "
    "in a single source record. It is not a prediction that the submitted chains "
    "will pair, assemble, fold, express, or remain stable together, and the "
    "absence of an observed pair is not evidence that a pair cannot form."
)

BLOCKED_DIMENSIONS = {
    "vh_vl_germline_pair_precedent": (
        "germline V-gene assignment now exists (para_scope.germline, OGRDB human "
        "reference) and returns a gene call for every variable domain tested, so "
        "the original blocker -- that nothing populated EvidenceClass.GERMLINE -- "
        "is resolved. What remains is the reference table this dimension actually "
        "needs: which VH germline gene is observed paired with which VL germline "
        "gene, and how often, across the 12,032 observed pairs. Building it means "
        "germline-assigning every unique paired variable domain, which is a "
        "separate offline job (one alignment per reference allele per component). "
        "Until that table exists this dimension reports nothing rather than "
        "guessing from the per-chain gene calls, which would describe gene "
        "co-occurrence in a candidate rather than precedent in a database"
    ),
    "combined_cdr_pattern_similarity": (
        "CDRs are short components and require a length-matched empirical null "
        "model before any similarity over them is meaningful (SPEC S21). That null "
        "is Gate 9 work. Exact CDR identity within a matched pair needs no chance "
        "correction and is reported per match instead"
    ),
    "format_relevance": (
        "requires the candidate's declared construct format, which the submission "
        "schema does not collect, plus a controlled vocabulary over a free-text "
        "column. Owned by Gate 10.5 -- the same term blocks SPEC S26"
    ),
}


@dataclass
class PairMatch:
    """One observed pair whose two chains both match the candidate's two chains.

    Every field describes a record that exists in a source database. Nothing here
    is inferred from the candidate.
    """

    heavy_record_id: str
    light_record_id: str
    heavy_component_type: ComponentType
    light_component_type: ComponentType

    # Worst of the two relationships: a pair is only as same-class as its weaker
    # chain, so a kappa/lambda substitution downgrades the whole pair rather than
    # being hidden behind an exact heavy match.
    pair_relationship: ComponentRelationship

    heavy_global_identity: float
    light_global_identity: float
    heavy_query_coverage: float
    light_query_coverage: float

    # SPEC S35's "joint coverage". Defined as the minimum of the two chains'
    # query coverage, not the mean: a pair cannot be better covered than its
    # least-covered chain, and averaging would let a fully-covered heavy hide a
    # light chain that barely aligned. PROVISIONAL -- the choice of min over other
    # summaries is not calibrated against anything.
    joint_coverage: float

    heavy_match_method: str  # "exact_checksum" | "pairwise_local_blosum62"
    light_match_method: str

    # The source records that link these two chains. Length is the corroboration
    # count for this one pair, never a count of distinct pairs.
    source_records: list[PairRecord] = field(default_factory=list)
    source_databases: list[str] = field(default_factory=list)
    source_accessions: list[str] = field(default_factory=list)

    # Three-state across the backing records: True if any record is structure
    # backed, False if all say no structure, None if none of them could tell.
    # None must never be coerced to False (CLAUDE.md).
    structure_supported: Optional[bool] = None
    structure_ids: list[str] = field(default_factory=list)

    # ---- Curated record content, promoted out of PairRecord.metadata -----------
    #
    # These are lists, not single values, and that is not defensiveness. One chain
    # pair is routinely backed by many records: trastuzumab's exact pair rests on
    # 11, which are 11 distinct therapeutic programmes built on the same chains,
    # exactly as Gate 7 found for the VH alone. Collapsing them to "the" name would
    # pick one arbitrarily and hide the other ten.
    #
    # They are typed fields rather than dict lookups for the same reason Gate 7
    # promoted source_databases: the information was always present, but reaching
    # it meant knowing the right string key, and a wrong key returns None instead
    # of raising.
    #
    # Which lists are populated depends on what kind of record backs the pair, and
    # that is deliberately left to be read from the data rather than mapped from a
    # source name. 4,562 of 5,686 pairs are structural records with a PDB entry,
    # method and resolution but no drug name; 1,124 are therapeutic arms with a
    # name and target but no experimental method. A consumer asks which fields are
    # populated; it does not need a hardcoded "SAbDab2 means structural" table that
    # would silently rot when a third source is added.
    therapeutic_names: list[str] = field(default_factory=list)
    targets: list[str] = field(default_factory=list)
    development_statuses: list[str] = field(default_factory=list)
    highest_clinical_trials: list[str] = field(default_factory=list)
    construct_formats: list[str] = field(default_factory=list)

    pdb_ids: list[str] = field(default_factory=list)
    structure_methods: list[str] = field(default_factory=list)
    structure_resolutions: list[float] = field(default_factory=list)

    # Convenience only; the lists above are authoritative. Includes the arm when
    # the record is one arm of a multi-specific construct, because an arm is not
    # the whole antibody -- flattening that at the display layer would recreate
    # the bispecific defect where one arm was credited with another's records.
    display_label: str = ""

    # True when either chain's hit came from a k-mer shortlist, whose hit set is
    # ~96% complete. A caller must not report "no closer pair exists" from these.
    shortlisted: bool = False

    @property
    def observation_count(self) -> int:
        """Source records backing this single pair. Not a count of pairs."""
        return len(self.source_records)

    @property
    def rank_key(self) -> tuple[float, float]:
        """Ordering: related-class pairs sit below same-class at equal identity."""
        weight = (
            1.0
            if self.pair_relationship is ComponentRelationship.EXACT_SAME_TYPE
            else RELATED_CLASS_RANK_PENALTY
        )
        combined = min(self.heavy_global_identity, self.light_global_identity)
        return (combined * weight, self.joint_coverage)


@dataclass
class PairingPrecedentResult:
    """Observed pairing evidence for one candidate heavy/light chain pair."""

    status: EvidenceStatus
    heavy_component_type: Optional[ComponentType] = None
    light_component_type: Optional[ComponentType] = None
    heavy_length: int = 0
    light_length: int = 0

    matches: list[PairMatch] = field(default_factory=list)

    # SPEC S15: entities then observations. Two numbers, always both, so neither
    # can be read as the other. On the real corpus the ratio is ~2.1.
    pairs_matched: int = 0
    raw_observation_count: int = 0
    exact_pair_match_count: int = 0

    best_heavy_identity: Optional[float] = None
    best_light_identity: Optional[float] = None
    best_joint_coverage: Optional[float] = None

    # SPEC S35 defines no weighted composite, so none is emitted and none is
    # withheld -- see the module docstring. These are raw measurements.
    computable_dimensions: dict[str, Optional[float]] = field(default_factory=dict)
    blocked_dimensions: dict[str, str] = field(default_factory=dict)

    heavy_search_status: EvidenceStatus = EvidenceStatus.AVAILABLE
    light_search_status: EvidenceStatus = EvidenceStatus.AVAILABLE

    evidence_class: EvidenceClass = EvidenceClass.DIRECT_PAIRING
    engine_version: str = PAIRING_PRECEDENT_VERSION
    scoring_version: str = ALIGNMENT_SCORING_VERSION
    calibration_version: str = CALIBRATION_VERSION
    reason: Optional[str] = None
    limitations: list[str] = field(default_factory=list)
    disclaimer: str = PAIRING_DISCLAIMER


def _checksum_of(record_id: str) -> str:
    """Component record ids are '<type>::<checksum>'; the pair index keys on the checksum."""
    _, _, checksum = record_id.partition("::")
    return checksum


def _collect(records: list[PairRecord], key: str) -> list[str]:
    """Deduplicated, sorted values for one metadata key across a pair's records.

    Absent keys and blank values are dropped rather than becoming empty strings:
    a record that does not state a target has no target, and "" would render as a
    populated-but-empty field.
    """
    values = {
        str(record.metadata.get(key)).strip()
        for record in records
        if record.metadata.get(key) not in (None, "")
    }
    return sorted(v for v in values if v)


def _collect_resolutions(records: list[PairRecord]) -> list[float]:
    """Resolutions as floats, skipping anything non-numeric rather than guessing."""
    resolutions = set()
    for record in records:
        raw = record.metadata.get("resolution")
        if raw in (None, ""):
            continue
        try:
            resolutions.add(float(raw))
        except (TypeError, ValueError):
            continue
    return sorted(resolutions)


def _display_label(records: list[PairRecord], therapeutic_names: list[str], pdb_ids: list[str]) -> str:
    """A single human-readable identity for the pair.

    Prefers a therapeutic name, qualified by arm where the record is one arm of a
    multi-specific construct; falls back to a PDB entry, then the accession. Never
    invents one -- a pair with no curated identity says so.
    """
    if therapeutic_names:
        label = therapeutic_names[0] if len(therapeutic_names) == 1 else (
            f"{therapeutic_names[0]} (+{len(therapeutic_names) - 1} more)"
        )
        arms = {
            record.metadata.get("arm_index")
            for record in records
            if record.metadata.get("arm_index") is not None
        }
        # Only qualify when the record really is one arm among several. Saying
        # "arm 1" on a conventional monospecific would imply an arm 2 that does
        # not exist.
        if arms and max(arms) > 1:
            return f"{label}, arm {min(arms)}"
        return label
    if pdb_ids:
        return pdb_ids[0] if len(pdb_ids) == 1 else f"{pdb_ids[0]} (+{len(pdb_ids) - 1} more)"
    return records[0].source_accession if records else ""


def _worst_relationship(
    first: ComponentRelationship, second: ComponentRelationship
) -> ComponentRelationship:
    """A pair is only as same-class as its weaker chain."""
    if ComponentRelationship.EXCLUDED in (first, second):
        return ComponentRelationship.EXCLUDED
    if ComponentRelationship.RELATED_CLASS in (first, second):
        return ComponentRelationship.RELATED_CLASS
    return ComponentRelationship.EXACT_SAME_TYPE


def _chain_hits(
    corpus: Corpus,
    declared_type: ComponentType,
    query: str,
    eligible_types: frozenset,
    top_k: int,
    shortlist_limit: Optional[int],
    min_query_coverage: float,
    min_subject_coverage: float,
) -> tuple[dict[str, dict], EvidenceStatus, Optional[str]]:
    """Exact and approximate hits for one chain, keyed by record id.

    Returns hits for *this chain alone*. Nothing here knows about pairing; the
    caller is responsible for only ever combining these through the pair index.
    """
    hits: dict[str, dict] = {}

    for record_id in corpus.exact.lookup(compute_sequence_checksum(query)):
        record = corpus.get(record_id)
        if record is None or record.component_type not in eligible_types:
            continue
        relationship = compare_eligibility(declared_type, record.component_type)
        if relationship is ComponentRelationship.EXCLUDED:
            continue
        hits[record_id] = {
            "record": record,
            "relationship": relationship,
            "global_identity": 1.0,
            "query_coverage": 1.0,
            "method": "exact_checksum",
            "shortlisted": False,
        }

    search_status = EvidenceStatus.AVAILABLE
    search_reason: Optional[str] = None
    for component_type in sorted(eligible_types, key=lambda t: t.value):
        try:
            found = corpus.approximate.search(
                component_type,
                query,
                top_k=top_k,
                min_query_coverage=min_query_coverage,
                min_subject_coverage=min_subject_coverage,
                shortlist_limit=shortlist_limit,
            )
        except Exception as exc:  # noqa: BLE001 -- must fill its slot, not abort the profile
            search_status = EvidenceStatus.ANALYSIS_FAILED
            search_reason = f"alignment search failed for {component_type.value}: {exc}"
            continue

        for hit in found:
            if hit.record_id in hits:  # an exact hit is never downgraded
                continue
            record = corpus.get(hit.record_id)
            if record is None:
                continue
            relationship = compare_eligibility(declared_type, record.component_type)
            if relationship is ComponentRelationship.EXCLUDED:
                continue
            hits[hit.record_id] = {
                "record": record,
                "relationship": relationship,
                "global_identity": hit.global_identity,
                "query_coverage": hit.query_coverage,
                "method": hit.alignment_method,
                "shortlisted": hit.shortlisted,
            }

    return hits, search_status, search_reason


def _paired_candidates(
    corpus: Corpus,
    heavy_hits: dict[str, dict],
    light_hits: dict[str, dict],
) -> list[PairMatch]:
    """Combine heavy and light hits **only** where a source record links them.

    This is the function that keeps the module honest. Every (heavy, light)
    combination is offered to the pair index, and only combinations the index
    already holds survive. A combination the index does not hold is discarded
    silently -- it is not weak evidence of pairing, it is no evidence of pairing,
    and emitting it at a low score would manufacture the synthesized pair
    CLAUDE.md prohibits.
    """
    matches: list[PairMatch] = []

    for heavy_id, heavy in heavy_hits.items():
        heavy_checksum = _checksum_of(heavy_id)
        for light_id, light in light_hits.items():
            records = corpus.pairs.lookup(heavy_checksum, _checksum_of(light_id))
            if not records:
                continue

            relationship = _worst_relationship(heavy["relationship"], light["relationship"])
            if relationship is ComponentRelationship.EXCLUDED:
                continue

            per_record_structure = [record.structure_available for record in records]
            if all(value is None for value in per_record_structure):
                structure_supported = None
            else:
                structure_supported = any(per_record_structure)

            structure_ids: list[str] = []
            for record in records:
                for ref in record.metadata.get("structure_ids", []) or []:
                    if ref not in structure_ids:
                        structure_ids.append(ref)
                pdb_id = record.metadata.get("pdb_id")
                if pdb_id and pdb_id not in structure_ids:
                    structure_ids.append(pdb_id)

            therapeutic_names = _collect(records, "therapeutic_name")
            pdb_ids = _collect(records, "pdb_id")

            matches.append(PairMatch(
                heavy_record_id=heavy_id,
                light_record_id=light_id,
                heavy_component_type=heavy["record"].component_type,
                light_component_type=light["record"].component_type,
                pair_relationship=relationship,
                heavy_global_identity=heavy["global_identity"],
                light_global_identity=light["global_identity"],
                heavy_query_coverage=heavy["query_coverage"],
                light_query_coverage=light["query_coverage"],
                joint_coverage=min(heavy["query_coverage"], light["query_coverage"]),
                heavy_match_method=heavy["method"],
                light_match_method=light["method"],
                source_records=records,
                source_databases=sorted({record.source_database for record in records}),
                source_accessions=sorted({record.source_accession for record in records}),
                structure_supported=structure_supported,
                structure_ids=structure_ids,
                therapeutic_names=therapeutic_names,
                targets=_collect(records, "target"),
                development_statuses=_collect(records, "development_status"),
                highest_clinical_trials=_collect(records, "highest_clinical_trial"),
                construct_formats=_collect(records, "format"),
                pdb_ids=pdb_ids,
                structure_methods=_collect(records, "method"),
                structure_resolutions=_collect_resolutions(records),
                display_label=_display_label(records, therapeutic_names, pdb_ids),
                shortlisted=bool(heavy["shortlisted"] or light["shortlisted"]),
            ))

    return matches


def _computable_dimensions(matches: list[PairMatch]) -> dict[str, Optional[float]]:
    """The five SPEC S35 dimensions that can be measured today, as raw values.

    Not weighted: SPEC S35 defines no weights, so there is nothing to weight them
    into. `None` marks a dimension whose inputs were themselves unknown and must
    not be read as zero.
    """
    if not matches:
        return {}
    best = max(matches, key=lambda m: m.rank_key)

    if best.structure_supported is None:
        parent_record_evidence: Optional[float] = None
    else:
        parent_record_evidence = float(best.structure_supported)

    return {
        "heavy_similarity": best.heavy_global_identity,
        "light_similarity": best.light_global_identity,
        "joint_coverage": best.joint_coverage,
        "direct_parent_record_evidence": parent_record_evidence,
        # Distinct unique pairs, not source records: breadth is how many different
        # observed pairs resemble the candidate, and one antibody solved a dozen
        # times is one pair (SPEC S15).
        "independent_paired_record_breadth": float(len(matches)),
    }


def compute_pairing_precedent(
    corpus: Corpus,
    heavy_sequence: str,
    light_sequence: Optional[str],
    heavy_component_type: ComponentType = ComponentType.VH,
    light_component_type: ComponentType = ComponentType.VL_KAPPA,
    top_k: int = 10,
    shortlist_limit: Optional[int] = None,
    min_query_coverage: float = MIN_QUERY_COVERAGE,
    min_subject_coverage: float = MIN_SUBJECT_COVERAGE,
) -> PairingPrecedentResult:
    """
    Observed pairing precedent for a candidate heavy/light chain pair (SPEC S35.1).

    Never raises on unscoreable input: an evidence dimension must fill its slot in
    a profile with an explicit status rather than abort the surrounding analysis.

    A candidate with no light chain returns `not_applicable`, not `not_found`.
    A single-domain construct has no light chain by design, and reporting "no
    pairing evidence found" for one would describe a deliberate architecture as a
    missing result.
    """
    # Both chains are reduced to their variable domains, independently, before
    # anything is compared -- corpus variable-domain components hold V-domains
    # only. Each note is surfaced because the input is then not what was compared.
    heavy_sequence, heavy_note = normalize_query_to_variable_domain(
        heavy_component_type.value, heavy_sequence or ""
    )
    light_note = None
    if light_sequence:
        light_sequence, light_note = normalize_query_to_variable_domain(
            light_component_type.value, light_sequence
        )

    result = PairingPrecedentResult(
        status=EvidenceStatus.NOT_FOUND,
        heavy_component_type=heavy_component_type,
        light_component_type=light_component_type,
        heavy_length=len(heavy_sequence or ""),
        light_length=len(light_sequence or ""),
        blocked_dimensions=dict(BLOCKED_DIMENSIONS),
    )
    if heavy_note is not None:
        result.limitations.append(f"Heavy chain: {heavy_note}")
    if light_note is not None:
        result.limitations.append(f"Light chain: {light_note}")

    if heavy_component_type not in HEAVY_TYPES or light_component_type not in LIGHT_TYPES:
        result.status = EvidenceStatus.UNSUPPORTED_COMPONENT
        result.reason = (
            f"pairing is defined for a heavy chain ({sorted(t.value for t in HEAVY_TYPES)}) "
            f"paired with a light chain ({sorted(t.value for t in LIGHT_TYPES)}); "
            f"got {heavy_component_type.value} + {light_component_type.value}"
        )
        return result

    if not heavy_sequence:
        result.status = EvidenceStatus.UNSUPPORTED_COMPONENT
        result.reason = "no heavy-chain sequence supplied; pairing requires both chains"
        return result

    if not light_sequence:
        result.status = EvidenceStatus.NOT_APPLICABLE
        result.reason = (
            "no light chain supplied. A single-domain construct has no light chain "
            "by design, so pairing precedent does not apply -- this is not an "
            "absence of evidence and must never be reported as one."
        )
        return result

    if len(corpus.pairs) == 0:
        result.status = EvidenceStatus.DATABASE_UNAVAILABLE
        result.reason = (
            "the corpus holds no observed pair records, so no pairing lookup was "
            "possible. This is not evidence that the candidate pair is unobserved."
        )
        return result

    heavy_hits, heavy_status, heavy_reason = _chain_hits(
        corpus, heavy_component_type, heavy_sequence, HEAVY_TYPES,
        top_k, shortlist_limit, min_query_coverage, min_subject_coverage,
    )
    light_hits, light_status, light_reason = _chain_hits(
        corpus, light_component_type, light_sequence, LIGHT_TYPES,
        top_k, shortlist_limit, min_query_coverage, min_subject_coverage,
    )
    result.heavy_search_status = heavy_status
    result.light_search_status = light_status
    for reason in (heavy_reason, light_reason):
        if reason:
            result.limitations.append(reason)

    matches = _paired_candidates(corpus, heavy_hits, light_hits)
    matches.sort(key=lambda m: (m.rank_key, m.heavy_record_id), reverse=True)
    result.matches = matches[:top_k]

    if not matches:
        result.status = EvidenceStatus.NOT_FOUND
        result.reason = (
            "both chains were searched, but no single source record links any "
            "heavy hit to any light hit. Independent heavy and light matches are "
            "never joined into a synthesized pair."
        )
        return result

    result.status = EvidenceStatus.AVAILABLE
    result.pairs_matched = len(matches)
    result.raw_observation_count = sum(match.observation_count for match in matches)
    result.exact_pair_match_count = sum(
        1 for match in matches
        if match.heavy_match_method == "exact_checksum"
        and match.light_match_method == "exact_checksum"
    )
    result.best_heavy_identity = max(m.heavy_global_identity for m in matches)
    result.best_light_identity = max(m.light_global_identity for m in matches)
    result.best_joint_coverage = max(m.joint_coverage for m in matches)
    result.computable_dimensions = _computable_dimensions(matches)

    if any(match.shortlisted for match in result.matches):
        result.limitations.append(
            "at least one chain hit came from a k-mer shortlist, whose hit set is "
            "~96% complete rather than exhaustive; do not report that no closer "
            "observed pair exists"
        )

    return result
