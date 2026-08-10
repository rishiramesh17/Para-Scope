"""
Deterministic sequence feature engine (Version 2, Gate 4).

Every value here is a **direct, reproducible computation over the input residues**
-- nothing is predicted, modelled, or inferred. Given the same sequence and the
same FEATURE_ENGINE_VERSION, the output is byte-identical, forever. That is the
whole contract: these features exist so that later gates can compare a candidate
against a distribution of known components on measurable axes, with the
comparison auditable back to arithmetic.

What these features are NOT
---------------------------
A feature is a description of the sequence, not a prediction about the molecule.
Specifically, nothing in this module says or implies that a construct will fold,
express, pair, remain soluble, avoid aggregation, be stable, be manufacturable,
be immunogenic, or bind anything (CLAUDE.md, prohibited claims). In particular:

  - `gravy` is Kyte-Doolittle average hydropathy. It is not a solubility score.
  - `isoelectric_point` is a theoretical pI from residue pKa values on the
    primary sequence alone. It ignores structure, post-translational
    modification, and local environment, and is not a formulation property.
  - `n_glycosylation_sequon_count` counts occurrences of the N-X-S/T motif. A
    sequon is a *motif observation*: whether any site is actually glycosylated
    depends on the expression system and structural accessibility, neither of
    which is in scope. It is not a liability score -- rule evaluation is Gate 5,
    and it consumes these counts rather than redefining them.
  - `free_cysteine_count` is parity-based and explicitly approximate; see below.

Extinction coefficient and molecular weight come from Biopython's
ProteinAnalysis so the constants are a maintained, citable implementation rather
than numbers transcribed here.
"""

from dataclasses import dataclass, field
from typing import Optional

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Bump on any change that alters a computed value. Stored on every FeatureSet so a
# feature can always be traced to the code that produced it (SPEC S19), and so a
# cached or persisted feature set is never silently reinterpreted under new rules.
FEATURE_ENGINE_VERSION = "1.0.0"

# The 20 standard amino acids. Composition is reported over exactly these, in a
# fixed order, so two FeatureSets are always directly comparable key-for-key.
STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolittle hydropathy. Restated here rather than taken from Biopython's
# internal table so the values used are visible at the point of use.
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Residues carrying charge at physiological pH, used for the net-charge estimate.
POSITIVE_RESIDUES = ("K", "R")
NEGATIVE_RESIDUES = ("D", "E")


class NonStandardResidueError(ValueError):
    """Raised when a sequence contains a residue the feature engine cannot score.

    Deliberately an error rather than a substitution or a skip: silently dropping
    an ambiguous residue (X, B, Z) or coercing it to a concrete one would modify
    the user's input, and every feature downstream would be computed over a
    sequence the user never submitted (CLAUDE.md: never silently modify input).
    """


@dataclass
class FeatureSet:
    """Deterministic features for one component. Descriptive only -- see module docstring."""

    sequence_length: int
    molecular_weight: float
    isoelectric_point: float
    net_charge_at_ph_7: int
    gravy: float
    aromaticity: float
    molar_extinction_coefficient_reduced: int
    molar_extinction_coefficient_cystines: int

    cysteine_count: int
    free_cysteine_count_estimate: Optional[int]
    n_glycosylation_sequon_count: int

    # Fractional composition over STANDARD_AMINO_ACIDS, summing to 1.0.
    composition: dict[str, float] = field(default_factory=dict)

    feature_engine_version: str = FEATURE_ENGINE_VERSION
    notes: list[str] = field(default_factory=list)


def _validate_sequence(sequence: str) -> str:
    if not sequence:
        raise NonStandardResidueError("Cannot compute features for an empty sequence")
    cleaned = sequence.strip().upper()
    unknown = sorted(set(cleaned) - set(STANDARD_AMINO_ACIDS))
    if unknown:
        raise NonStandardResidueError(
            f"Sequence contains residues the feature engine cannot score: {unknown}. "
            f"These are not substituted or dropped -- the caller must decide how to "
            f"handle an ambiguous or non-standard residue and record that decision."
        )
    return cleaned


def count_n_glycosylation_sequons(sequence: str) -> int:
    """
    Count N-X-S/T sequons, where X is any residue except proline.

    A motif count, not a liability assessment: actual glycosylation depends on the
    expression system and structural accessibility. Proline at X is excluded
    because it blocks the motif -- that exclusion is part of the motif definition,
    not a prediction about any particular site.
    """
    count = 0
    for i in range(len(sequence) - 2):
        if sequence[i] == "N" and sequence[i + 1] != "P" and sequence[i + 2] in ("S", "T"):
            count += 1
    return count


def compute_features(sequence: str) -> FeatureSet:
    """
    Compute deterministic features for a single component sequence.

    Raises NonStandardResidueError if the sequence contains anything outside the
    20 standard amino acids -- features are never computed over a silently
    altered sequence.
    """
    cleaned = _validate_sequence(sequence)
    analysis = ProteinAnalysis(cleaned)

    length = len(cleaned)
    counts = {aa: cleaned.count(aa) for aa in STANDARD_AMINO_ACIDS}
    composition = {aa: counts[aa] / length for aa in STANDARD_AMINO_ACIDS}

    gravy = sum(KYTE_DOOLITTLE[aa] for aa in cleaned) / length
    net_charge = sum(counts[aa] for aa in POSITIVE_RESIDUES) - sum(
        counts[aa] for aa in NEGATIVE_RESIDUES
    )

    cysteines = counts["C"]
    notes: list[str] = []

    # Parity only. Whether a cysteine is actually free depends on disulfide
    # connectivity, which is a structural property and out of scope for a
    # sequence-first system (CLAUDE.md, structural boundary). An odd count means
    # at least one cysteine cannot be paired *within this component*; an even
    # count does NOT establish that all are paired. Reported as an estimate with
    # that limitation attached, rather than as a fact.
    free_cysteine_estimate: Optional[int] = None
    if cysteines:
        free_cysteine_estimate = cysteines % 2
        notes.append(
            "free_cysteine_count_estimate is parity-based: it reflects whether "
            "cysteines can pair within this component, not observed disulfide "
            "connectivity, which requires structural data."
        )

    reduced, cystines = analysis.molar_extinction_coefficient()

    return FeatureSet(
        sequence_length=length,
        molecular_weight=float(analysis.molecular_weight()),
        isoelectric_point=float(analysis.isoelectric_point()),
        net_charge_at_ph_7=int(net_charge),
        gravy=float(gravy),
        aromaticity=float(analysis.aromaticity()),
        molar_extinction_coefficient_reduced=int(reduced),
        molar_extinction_coefficient_cystines=int(cystines),
        cysteine_count=int(cysteines),
        free_cysteine_count_estimate=free_cysteine_estimate,
        n_glycosylation_sequon_count=count_n_glycosylation_sequons(cleaned),
        composition=composition,
        notes=notes,
    )
