"""
Sequence motif observation engine (Version 2, Gate 5).

This module locates specific residue patterns in a sequence and reports where
they are, in IMGT coordinates. That is all it does.

Terminology, deliberately
-------------------------
These patterns are widely *described* in the antibody-engineering literature as
"sequence liabilities" -- sites associated with chemical degradation, charge
heterogeneity, or manufacturing variability. This module does not adopt that
framing in its output, because a motif's presence is an observation about
residues while a liability is a claim about molecular behaviour, and the second
does not follow from the first. A motif may be fully buried, structurally
constrained, or entirely inert in a given construct; whether it matters depends
on structure, expression system, and formulation, none of which is in scope for a
sequence-first system (CLAUDE.md, structural boundary).

So: `MotifObservation`, not "liability". A count, a position, and a citation for
why the pattern is of interest -- never a severity, a risk level, a penalty, or a
score. Nothing here says or implies that a construct will aggregate, degrade,
be unstable, be immunogenic, be unmanufacturable, or fail (CLAUDE.md, prohibited
claims). Downstream gates may *count* these observations and compare those counts
against a corpus distribution; they may not convert them into a risk.

Why positions matter
--------------------
A motif in a framework region and the same motif in CDR-H3 are different
observations, and the literature treats them differently -- so every observation
carries its IMGT position and region. Positions come from the numbering already
stored on the component, never from re-deriving coordinates here.
"""

from dataclasses import dataclass, field
from typing import Optional

from para_scope.numbering import (
    IMGT_GAP_CHARACTERS,
    IMGT_REGION_BOUNDARIES,
    NumberedComponent,
)

# Bump on any change to the rule set or to how observations are located. Stored on
# every result so an observation can always be traced to the rules that produced
# it, and so a stored result is never reinterpreted under a later rule set.
#   1.0.0 -- initial two-residue and single-residue motif rules
MOTIF_RULESET_VERSION = "1.0.0"


@dataclass(frozen=True)
class MotifRule:
    """One residue pattern to locate. Descriptive metadata only -- no severity, no weight."""

    rule_id: str
    name: str
    # Residue patterns, each a string of literal residues. Matching is exact and
    # case-sensitive against the ungapped sequence; no regex, so a rule cannot
    # accidentally match across a region boundary or a gap.
    patterns: tuple[str, ...]
    # What the pattern is, and why it is tracked. Phrased as an association
    # reported in the literature, never as a property of the input.
    basis: str


# Rules are ordered by pattern length then id purely for deterministic output.
MOTIF_RULES: tuple[MotifRule, ...] = (
    MotifRule(
        rule_id="deamidation_ng",
        name="Asparagine deamidation motif (NG)",
        patterns=("NG",),
        basis=(
            "Asn-Gly is the most frequently reported asparagine deamidation "
            "context; rate depends strongly on local conformation and pH, which "
            "sequence alone does not determine."
        ),
    ),
    MotifRule(
        rule_id="deamidation_ns_nt",
        name="Asparagine deamidation motif (NS/NT)",
        patterns=("NS", "NT"),
        basis=(
            "Asn followed by Ser or Thr is reported as a slower deamidation "
            "context than Asn-Gly."
        ),
    ),
    MotifRule(
        rule_id="isomerization_dg",
        name="Aspartate isomerization motif (DG)",
        patterns=("DG",),
        basis="Asp-Gly is the most commonly reported aspartate isomerization context.",
    ),
    MotifRule(
        rule_id="isomerization_ds_dt",
        name="Aspartate isomerization motif (DS/DT)",
        patterns=("DS", "DT"),
        basis="Asp followed by Ser or Thr is reported as an isomerization context.",
    ),
    MotifRule(
        rule_id="fragmentation_dp",
        name="Acid-labile bond motif (DP)",
        patterns=("DP",),
        basis=(
            "The Asp-Pro peptide bond is reported as unusually labile under "
            "acidic conditions."
        ),
    ),
    MotifRule(
        rule_id="oxidation_met",
        name="Methionine (oxidation-susceptible residue)",
        patterns=("M",),
        basis=(
            "Methionine thioether is susceptible to oxidation; solvent exposure "
            "governs whether any given residue is affected, and is not derivable "
            "from sequence."
        ),
    ),
    MotifRule(
        rule_id="oxidation_trp",
        name="Tryptophan (oxidation-susceptible residue)",
        patterns=("W",),
        basis=(
            "Tryptophan is reported as oxidation-susceptible. Note it is also a "
            "common and functionally important CDR residue, so its presence is "
            "unremarkable in isolation."
        ),
    ),
)

RULES_BY_ID = {rule.rule_id: rule for rule in MOTIF_RULES}


@dataclass(frozen=True)
class MotifObservation:
    """One located occurrence of one motif. An observation, not an assessment."""

    rule_id: str
    rule_name: str
    matched_residues: str
    # IMGT position of the first residue of the match, with its insertion code.
    # None only if the component carries no numbering for that residue.
    imgt_position: Optional[int]
    imgt_insertion_code: Optional[str]
    # IMGT region containing the first residue (FR1/CDR1/.../FR4), or None if the
    # position falls outside the defined V-domain boundaries.
    region: Optional[str]
    # 0-based offset into the ungapped sequence, so a caller can always recover
    # the match without re-deriving numbering.
    sequence_offset: int


@dataclass
class MotifReport:
    """All motif observations for one component, plus what produced them."""

    observations: list[MotifObservation] = field(default_factory=list)
    counts_by_rule: dict[str, int] = field(default_factory=dict)
    ruleset_version: str = MOTIF_RULESET_VERSION
    # Present when positions could not be resolved; the observations are still
    # valid, they just carry no IMGT coordinates.
    limitations: list[str] = field(default_factory=list)

    # Stated verbatim on any user-facing surface that presents these counts.
    disclaimer: str = (
        "These are sequence motif observations, not predictions. The presence of "
        "a motif does not establish that a construct will degrade, aggregate, be "
        "unstable, be immunogenic, or behave in any particular way."
    )


def _region_for_position(position: Optional[int]) -> Optional[str]:
    if position is None:
        return None
    for region, (start, end) in IMGT_REGION_BOUNDARIES.items():
        if start <= position <= end:
            return region
    return None


def _ungapped_with_positions(
    numbered: Optional[NumberedComponent],
) -> tuple[str, list[tuple[Optional[int], Optional[str]]]]:
    """
    Ungapped residues plus a parallel list of their IMGT (position, insertion_code).

    Built from imgt_numbering so offsets and coordinates cannot drift apart. Gap
    placeholders are excluded from both, exactly as in region extraction -- they
    are positional metadata, not residues, and a motif must never be reported as
    spanning one.
    """
    if numbered is None or not numbered.imgt_numbering:
        return "", []

    residues: list[str] = []
    coordinates: list[tuple[Optional[int], Optional[str]]] = []
    for (position, code), residue in numbered.imgt_numbering:
        if residue in IMGT_GAP_CHARACTERS:
            continue
        residues.append(residue)
        cleaned_code = code.strip() if isinstance(code, str) else None
        coordinates.append((position, cleaned_code or None))
    return "".join(residues), coordinates


def find_motifs(
    sequence: str,
    numbered: Optional[NumberedComponent] = None,
    rules: tuple[MotifRule, ...] = MOTIF_RULES,
) -> MotifReport:
    """
    Locate every motif in `sequence`, reporting IMGT coordinates when `numbered` is given.

    `sequence` is used as-is when no numbering is supplied, in which case
    observations carry sequence offsets but no IMGT coordinates -- and that
    limitation is recorded rather than left implicit. When `numbered` is supplied,
    its ungapped residues are used, so a motif can never be reported as spanning
    an IMGT gap.
    """
    report = MotifReport()

    if numbered is not None:
        residues, coordinates = _ungapped_with_positions(numbered)
        if not residues:
            residues = sequence.strip().upper()
            coordinates = [(None, None)] * len(residues)
            report.limitations.append(
                "Component carries no IMGT numbering; observations report sequence "
                "offsets only, without IMGT coordinates or region assignment."
            )
    else:
        residues = sequence.strip().upper()
        coordinates = [(None, None)] * len(residues)
        report.limitations.append(
            "No numbering supplied; observations report sequence offsets only, "
            "without IMGT coordinates or region assignment."
        )

    for rule in rules:
        count = 0
        for pattern in rule.patterns:
            span = len(pattern)
            for offset in range(len(residues) - span + 1):
                if residues[offset : offset + span] != pattern:
                    continue
                position, code = coordinates[offset] if offset < len(coordinates) else (None, None)
                report.observations.append(
                    MotifObservation(
                        rule_id=rule.rule_id,
                        rule_name=rule.name,
                        matched_residues=pattern,
                        imgt_position=position,
                        imgt_insertion_code=code,
                        region=_region_for_position(position),
                        sequence_offset=offset,
                    )
                )
                count += 1
        # Recorded for every rule, including zero. A rule that found nothing is
        # meaningfully different from a rule that was never evaluated, and a bare
        # absence from the dict could not distinguish the two.
        report.counts_by_rule[rule.rule_id] = count

    # Deterministic ordering: position in the sequence, then rule id.
    report.observations.sort(key=lambda o: (o.sequence_offset, o.rule_id))
    return report
