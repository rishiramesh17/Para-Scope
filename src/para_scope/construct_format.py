"""
Construct format vocabulary and format relevance (Version 2, Gate 10.5).

Why this gate exists at all
---------------------------
SPEC §26 gives format relevance 0.10 of the therapeutic-precedent composite, and
SPEC §35 needs the same term. Nothing produced it: it requires the candidate's
*declared* construct format, which the submission schema never collected, plus a
controlled vocabulary over a free-text column. Neither falls out of any other
gate, so left unowned every gate would have completed and the composite would
have shipped computed over 0.90 of its own definition with nothing flagging the
gap (docs/ROADMAP.md).

Four axes, not one taxonomy
---------------------------
The earlier read of this column was that a classifier over its 104 distinct
values was infeasible. That is true of a *flat* taxonomy and false of the data:
`"Bispecific Whole mAb ADC"` is three independent facts, and a flat scheme has to
invent a bucket for every combination. Parsed as independent axes instead,
measured against the real file (1,133 rows):

    specificity    100.0% resolved
    architecture    98.1% resolved  (21 rows unparseable)
    modification   100.0% resolved
    species        100.0% resolved

The 21 unresolved rows are bespoke domain-topology strings like
`"Bispecific Mixed Domains (L-kappa-H-gamma1_[VH-G1(CH1-h)]_L-kappa)"`. They
return `UNKNOWN`, and `UNKNOWN` never becomes agreement or disagreement -- it
makes that axis uncomparable and says so.

What format relevance means, and what it must not mean
------------------------------------------------------
It measures **how comparable a matched record is to the candidate**, not how good
either is. A scFv match for a whole-mAb candidate is weaker *evidence about this
candidate* because the constructs differ; it is not a statement that either
molecule is worse, less developable, or less likely to work (CLAUDE.md).

Two consequences that are easy to get wrong:

- **An uncomparable format is not a relevance of 0.0.** Zero would assert the
  formats are maximally unrelated. When too few axes are comparable the term is
  `None` with `INSUFFICIENT_EVIDENCE`, so §26's weight can be reported as missing
  rather than silently satisfied by a fabricated zero.
- **A missing declared format is the normal case today**, not an error. The
  submission schema is only now gaining the field. Absent it, relevance is
  uncomparable and §26's 0.10 stays explicitly unmet -- which is the honest state
  and the one the roadmap insisted must never be hidden by renormalising over the
  remaining 0.90.

Inference recorded as inference
-------------------------------
Absence of "bispecific" is read as monospecific. That is a convention of *this
column* -- Thera-SAbDab states multispecificity explicitly -- not a universal
truth about format strings, and it is recorded on the parse so a reader can see
the axis was inferred rather than stated. Absence of a species word is NOT read as
human: that would be a claim, so it stays `UNSTATED`.
"""

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from para_scope.models import EvidenceStatus

# Bump when an axis, a keyword rule, or an axis weight changes. A vocabulary
# change re-labels every stored format without altering any sequence, so it is
# versioned independently of any engine.
#   1.0.0 -- initial four axes, measured against Thera-SAbDab 2025-03-18
#   1.1.0 -- architecture weight raised above the other axes combined, so an
#            architecture mismatch cannot score better than two modifier
#            mismatches (see AXIS_WEIGHTS)
#   1.2.0 -- MAB_UNSPECIFIED treated as unresolved for comparison: it states
#            the construct is an antibody, not that it is not a whole mAb
FORMAT_VOCABULARY_VERSION = "1.2.0"

# Axis weights for the composite relevance term. PROVISIONAL, like the grade
# thresholds and the coverage floor: specificity and architecture dominate because
# they determine whether two constructs are the same *kind* of molecule, while a
# conjugate or a host species changes the molecule less than turning a whole mAb
# into an scFv does. No external ground truth exists for these.
AXIS_WEIGHTS: dict[str, float] = {
    "architecture": 0.55,
    "specificity": 0.25,
    "modification": 0.12,
    "species": 0.08,
}
# Architecture's weight deliberately exceeds the other three combined (0.55 vs
# 0.45). At a first cut it did not -- 0.40 against 0.60 -- and the consequence was
# backwards: "Whole mAb" vs "scFv" scored 0.556 while "Whole mAb" vs "Bispecific
# Whole mAb ADC" scored 0.444, so disagreeing on the one axis that decides what
# kind of molecule a construct is scored *better* than disagreeing on two
# modifiers. An scFv carries no constant region and links its domains differently,
# which changes what a VH match means more than adding a second specificity does.

# Fewest axes that must be comparable before a relevance number is reported.
# Below this the term is None with INSUFFICIENT_EVIDENCE rather than a figure
# resting on one axis.
MINIMUM_COMPARABLE_AXES = 2

# Architecture must be comparable, whatever else is.
#
# Without this, two constructs whose architecture is UNKNOWN on both sides scored
# relevance 1.0 because they happened to share "bispecific" and an *inferred*
# absence of conjugation -- a perfect-format-match reading while the axis that
# decides what kind of molecule each is had never been assessed. Specificity,
# modification and species modify a construct; architecture is the construct.
REQUIRED_AXIS = "architecture"

FORMAT_RELEVANCE_DISCLAIMER = (
    "Format relevance measures how comparable a matched record's construct format "
    "is to the declared candidate format. It is not a judgement of either "
    "construct: a format difference makes a match less informative about this "
    "candidate, not worse, less developable, or less likely to succeed."
)


class Specificity(str, Enum):
    MONOSPECIFIC = "monospecific"
    BISPECIFIC = "bispecific"
    TRISPECIFIC = "trispecific"
    TETRASPECIFIC = "tetraspecific"
    UNKNOWN = "unknown"


class Architecture(str, Enum):
    WHOLE_MAB = "whole_mab"
    # "mAb" stated without a sub-architecture. Deliberately distinct from
    # WHOLE_MAB: 60 rows say "Bispecific mAb" without saying whole, and promoting
    # those to WHOLE_MAB would assert a constant region the text never claimed.
    MAB_UNSPECIFIED = "mab_unspecified"
    FAB = "fab"
    SCFV = "scfv"
    SINGLE_DOMAIN = "single_domain"
    DUAL_VARIABLE_DOMAIN = "dual_variable_domain"
    FUSION_PROTEIN = "fusion_protein"
    # More than one architecture in one construct (e.g. "Mixed mAb and scFv").
    MIXED = "mixed"
    UNKNOWN = "unknown"


class Modification(str, Enum):
    NONE = "none"
    ADC = "antibody_drug_conjugate"
    RADIOLABELLED = "radiolabelled"
    FUSION = "fusion"
    PEGYLATED = "pegylated"
    UNKNOWN = "unknown"


class HostSpecies(str, Enum):
    # Absence of a species word is UNSTATED, never HUMAN -- inferring human from
    # silence would be a claim the column does not make.
    UNSTATED = "unstated"
    CANINE = "canine"
    FELINE = "feline"
    MOUSE = "mouse"
    RAT = "rat"
    EQUINE = "equine"


class AxisComparison(str, Enum):
    AGREES = "agrees"
    DIFFERS = "differs"
    # One or both sides unknown/unstated. Never counted as agreement or
    # disagreement -- an unobserved axis is not a matching axis.
    UNCOMPARABLE = "uncomparable"


@dataclass(frozen=True)
class ConstructFormat:
    """One construct format, parsed onto four independent axes."""

    specificity: Specificity
    architecture: Architecture
    modification: Modification
    species: HostSpecies

    # The source text, verbatim. Kept because the vocabulary is lossy by design
    # and a reader must be able to see what was actually written.
    raw: str = ""
    # Axes whose value was inferred from absence rather than stated (see module
    # docstring on the monospecific convention).
    inferred_axes: tuple[str, ...] = ()
    vocabulary_version: str = FORMAT_VOCABULARY_VERSION

    @property
    def fully_resolved(self) -> bool:
        return (
            self.specificity is not Specificity.UNKNOWN
            and self.architecture is not Architecture.UNKNOWN
            and self.modification is not Modification.UNKNOWN
        )


@dataclass
class FormatRelevance:
    """Per-axis comparison plus the SPEC §26 term, or an explicit reason there isn't one."""

    status: EvidenceStatus
    # None whenever fewer than MINIMUM_COMPARABLE_AXES axes could be compared.
    # Never 0.0 as a stand-in: zero asserts maximal format difference.
    relevance: Optional[float] = None

    by_axis: dict[str, AxisComparison] = field(default_factory=dict)
    comparable_axes: int = 0
    agreeing_axes: int = 0
    # Share of total axis weight that actually backed the number, so a consumer
    # can see a relevance of 1.0 resting on 0.55 of the weight for what it is.
    comparable_weight: float = 0.0

    declared: Optional[ConstructFormat] = None
    matched: Optional[ConstructFormat] = None

    vocabulary_version: str = FORMAT_VOCABULARY_VERSION
    reason: Optional[str] = None
    disclaimer: str = FORMAT_RELEVANCE_DISCLAIMER


def _normalise(raw: str) -> str:
    """Lowercase and collapse whitespace.

    Collapsing matters: the real file contains `"Whole  mAb"` with a double space,
    which a literal match would drop into UNKNOWN for a typing artifact.
    """
    return re.sub(r"\s+", " ", raw.strip().lower())


_SPECIFICITY_KEYWORDS: tuple[tuple[str, Specificity], ...] = (
    ("tetraspecific", Specificity.TETRASPECIFIC),
    ("trispecific", Specificity.TRISPECIFIC),
    ("bispecific", Specificity.BISPECIFIC),
    ("biparatopic", Specificity.BISPECIFIC),
    ("bifunctional", Specificity.BISPECIFIC),
)

_SPECIES_KEYWORDS: tuple[tuple[str, HostSpecies], ...] = (
    ("canine", HostSpecies.CANINE),
    ("feline", HostSpecies.FELINE),
    ("murine", HostSpecies.MOUSE),
    ("(mouse)", HostSpecies.MOUSE),
    ("(rat)", HostSpecies.RAT),
    ("equine", HostSpecies.EQUINE),
)


def _parse_architecture(text: str) -> Architecture:
    hits: list[Architecture] = []
    if "whole mab" in text or "whole igg" in text or "whole igm" in text:
        hits.append(Architecture.WHOLE_MAB)
    if re.search(r"\bfab\b", text):
        hits.append(Architecture.FAB)
    if "scfv" in text:
        hits.append(Architecture.SCFV)
    if any(k in text for k in ("single domain", "nanobody", "vhh", "vh-vh")):
        hits.append(Architecture.SINGLE_DOMAIN)
    if "dual variable domain" in text:
        hits.append(Architecture.DUAL_VARIABLE_DOMAIN)
    if "fusion protein" in text:
        hits.append(Architecture.FUSION_PROTEIN)

    if not hits:
        # "mAb" or "antibody" stated with no sub-architecture. Honest middle
        # ground between WHOLE_MAB (asserts a constant region) and UNKNOWN
        # (discards the fact that it is an antibody at all).
        if re.search(r"\bmab\b|\bantibody\b|\bigg\b|t-cell engager", text):
            return Architecture.MAB_UNSPECIFIED
        return Architecture.UNKNOWN
    if len(hits) > 1 or "mixed" in text:
        return Architecture.MIXED
    return hits[0]


def parse_format(raw: Optional[str]) -> ConstructFormat:
    """
    Parse a free-text construct format onto the four axes.

    Never raises and never guesses: an axis the text does not determine returns
    UNKNOWN (or UNSTATED for species), which downstream makes that axis
    uncomparable rather than matching or mismatching.
    """
    if not raw or not raw.strip():
        return ConstructFormat(
            specificity=Specificity.UNKNOWN,
            architecture=Architecture.UNKNOWN,
            modification=Modification.UNKNOWN,
            species=HostSpecies.UNSTATED,
            raw=raw or "",
        )

    text = _normalise(raw)
    inferred: list[str] = []

    specificity = Specificity.UNKNOWN
    for keyword, value in _SPECIFICITY_KEYWORDS:
        if keyword in text:
            specificity = value
            break
    else:
        # Convention of this column, not a universal rule -- recorded as inferred.
        specificity = Specificity.MONOSPECIFIC
        inferred.append("specificity")

    architecture = _parse_architecture(text)

    if "adc" in text or "drug conjugate" in text:
        modification = Modification.ADC
    elif "radiolab" in text:
        modification = Modification.RADIOLABELLED
    elif "pegylat" in text:
        modification = Modification.PEGYLATED
    elif "fusion" in text:
        modification = Modification.FUSION
    else:
        modification = Modification.NONE
        inferred.append("modification")

    species = HostSpecies.UNSTATED
    for keyword, value in _SPECIES_KEYWORDS:
        if keyword in text:
            species = value
            break

    return ConstructFormat(
        specificity=specificity,
        architecture=architecture,
        modification=modification,
        species=species,
        raw=raw,
        inferred_axes=tuple(inferred),
    )


# Values that mean "this axis was not determined", per axis.
#
# MAB_UNSPECIFIED counts as unresolved for *comparison* even though it is a real
# parse result. "Bispecific mAb" (60 rows) says the construct is an antibody and
# stops there -- it does not say the construct is not a whole mAb. Treating it as
# a definite mismatch against WHOLE_MAB scored 0.130 for a pair that may well be
# identical in architecture, which is a confident answer to a question the source
# text never answered. It stays a distinct enum value because the fact "an
# antibody of unstated sub-architecture" is worth recording; it just cannot
# support an agreement or a disagreement.
_UNRESOLVED = {
    "specificity": {Specificity.UNKNOWN},
    "architecture": {Architecture.UNKNOWN, Architecture.MAB_UNSPECIFIED},
    "modification": {Modification.UNKNOWN},
    "species": {HostSpecies.UNSTATED},
}


def compute_format_relevance(
    declared: Optional[ConstructFormat],
    matched: Optional[ConstructFormat],
) -> FormatRelevance:
    """
    Compare a declared candidate format against a matched record's format.

    Returns a per-axis breakdown plus SPEC §26's 0.10 term. The term is None --
    with INSUFFICIENT_EVIDENCE and a reason -- whenever fewer than
    MINIMUM_COMPARABLE_AXES axes could be compared, including the common case
    where no format was declared at all.
    """
    result = FormatRelevance(
        status=EvidenceStatus.NOT_APPLICABLE, declared=declared, matched=matched
    )

    if declared is None:
        result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
        result.reason = (
            "No construct format was declared for the candidate, so there is "
            "nothing to compare a matched record's format against. SPEC §26's "
            "0.10 format-relevance weight is unmet and must be reported as "
            "missing, not satisfied by a substituted value."
        )
        return result
    if matched is None:
        result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
        result.reason = "The matched record carries no construct format to compare."
        return result

    weighted_total = 0.0
    weighted_agreeing = 0.0
    for axis, weight in AXIS_WEIGHTS.items():
        declared_value = getattr(declared, axis)
        matched_value = getattr(matched, axis)
        unresolved = _UNRESOLVED[axis]
        if declared_value in unresolved or matched_value in unresolved:
            result.by_axis[axis] = AxisComparison.UNCOMPARABLE
            continue
        agrees = declared_value == matched_value
        result.by_axis[axis] = AxisComparison.AGREES if agrees else AxisComparison.DIFFERS
        result.comparable_axes += 1
        weighted_total += weight
        if agrees:
            result.agreeing_axes += 1
            weighted_agreeing += weight

    result.comparable_weight = weighted_total

    if result.by_axis.get(REQUIRED_AXIS) is AxisComparison.UNCOMPARABLE:
        result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
        result.reason = (
            f"Construct architecture could not be determined for "
            f"{'both' if declared.architecture is matched.architecture else 'one'} "
            f"side(s) (declared {declared.architecture.value}, matched "
            f"{matched.architecture.value}). Agreement on the remaining axes "
            f"cannot establish format relevance: specificity, modification and "
            f"species modify a construct, architecture is the construct. "
            f"Relevance is unset rather than scored from the modifiers alone."
        )
        return result

    if result.comparable_axes < MINIMUM_COMPARABLE_AXES:
        result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
        uncomparable = sorted(
            axis for axis, comparison in result.by_axis.items()
            if comparison is AxisComparison.UNCOMPARABLE
        )
        result.reason = (
            f"Only {result.comparable_axes} of {len(AXIS_WEIGHTS)} format axes "
            f"could be compared (uncomparable: {uncomparable}); at least "
            f"{MINIMUM_COMPARABLE_AXES} are required. Relevance is unset rather "
            f"than zero -- zero would assert the formats are maximally unrelated."
        )
        return result

    # Renormalised over the axes actually compared. Doing so is sound *here*, and
    # unsound for SPEC §26's own composite, because the comparable-axis count
    # travels with the number: a caller can see 2 of 4 axes backed it. The §26
    # composite has no equivalent field, which is why dropping a term there and
    # renormalising silently is prohibited.
    result.relevance = weighted_agreeing / weighted_total
    result.status = EvidenceStatus.AVAILABLE
    return result
