"""
Human repertoire window coverage (Version 2, Gate 6).

What this measures
------------------
The fraction of a sequence's overlapping 9-residue windows that appear verbatim
in a reference set of peptides drawn from human antibody repertoires. Nothing
more. Every window is an exact string lookup against a fixed, versioned local
set; the result is a count divided by a count.

Terminology, deliberately
-------------------------
This metric is widely called an *antibody humanness score*, and the underlying
reference set is BioPhi's OASis peptide set (Prihoda et al., mAbs 2022). This
module does not adopt that framing, and the reason is empirical rather than
stylistic.

Measured with this exact reference set (see `reports/`):

    trastuzumab VH (humanized)   0.7054
    trastuzumab VL (humanized)   0.8283
    muromonab  VH (mouse)        0.5225
    caplacizumab (VHH, llama)    0.6333

A camelid VHH scores far above a mouse antibody and close to a humanized one.
Whatever this quantity tracks, it is therefore *not* "how human is this
sequence" -- a metric that cannot separate llama from humanized cannot be
reported as a measure of human origin. Naming it "humanness" or "naturalness"
would license exactly the inference the numbers refuse to support, and would
imply an immunogenicity claim that CLAUDE.md prohibits outright.

So: `human_window_coverage` -- a name that states the comparison performed
(how many windows were observed, in which reference set) rather than a
conclusion drawn from it. The reference set's identity and version are fields on
the result rather than facts baked into the name, so a future reference set does
not force a rename. This follows the convention already in force for
`Observed Repertoire Precedent` and `Database-Relative Novelty`.

Nothing here says or implies that a construct will be immunogenic or
non-immunogenic, will be tolerated in humans, is safe, is humanized, or is
suitable for clinical use. High coverage is not evidence of low immunogenicity;
the reference set's own authors note that common human peptides include
low-complexity runs such as `GGGGGGGGG` that are human and undesirable at once.

Two behaviours of the underlying library that are deliberately not inherited
-----------------------------------------------------------------------------
1. `promb` raises `ZeroDivisionError` for any sequence shorter than the window.
   Most CDRs are shorter than nine residues (IMGT CDR-L2 is around three), so
   the metric is simply undefined for them. Handled here as NOT_APPLICABLE.

2. `promb` silently scores sequences containing non-standard residues. A window
   containing `X` matches nothing, so an unresolvable residue *lowers* coverage
   and reads as "less human" when it in fact means "we could not evaluate this
   window". That is missing evidence presented as a low score, which CLAUDE.md
   forbids. Handled here as UNSUPPORTED_COMPONENT, with the offending residues
   named. The input is never cleaned, substituted, or truncated to make it
   scoreable (`features.py` raises in the equivalent situation; this module
   returns a status instead, because an evidence dimension must always be able
   to fill its slot in a profile rather than abort the surrounding analysis).

Never a neutral default
-----------------------
When the reference set cannot be loaded, `coverage` stays `None` and the status
says why. It never falls back to 0.0, to 0.5, or to a corpus mean. A neutral
number reads as "average coverage" when it actually means "we did not look",
and the two must remain distinguishable at every layer.
"""

from dataclasses import dataclass, field
from typing import Optional

from antigens.features import STANDARD_AMINO_ACIDS
from antigens.models import EvidenceStatus

# Bump on any change that alters a computed coverage value, including a change of
# reference set, window length, or the minimum-window rule. Stored on every result
# so a stored value is never reinterpreted under later rules (SPEC §19).
#   1.0.0 -- initial 9-mer exact-match coverage over the OASis human-OAS set
WINDOW_COVERAGE_VERSION = "1.0.0"

# Window length. Fixed at 9 to match the reference set, which is built at that
# length -- this is not a free parameter, and changing it requires a different
# reference set rather than a different argument.
WINDOW_LENGTH = 9

# Identifier passed to promb, and the bundled file it resolves to. Recorded on
# every result: SPEC §11.2 requires the reference source and set version to travel
# with the assignment, not to be implied by the code that produced it.
REFERENCE_SET_ID = "human-oas"
REFERENCE_SET_FILE = "OASis_9mers_v1_10perc_subjects.txt.gz"
REFERENCE_SET_DESCRIPTION = (
    "Peptides of length 9 observed in at least 10% of subjects in the Observed "
    "Antibody Space (OAS) human repertoire collection, as distributed with promb "
    "(BioPhi OASis reference set)."
)

# Below this many windows, a fraction is reported as INSUFFICIENT_EVIDENCE rather
# than as a value: over three windows, "0.667" carries a precision the measurement
# does not have, and short components need a length-matched empirical null before
# any such number means anything (CLAUDE.md, short-component chance correction;
# SPEC §21). The raw counts are still returned, so a caller building that null
# model in Gate 9 has everything it needs.
#
# PROVISIONAL. 20 windows (a 28-residue component) is a defensible default, not a
# calibrated threshold -- the same standing as the 0.70 coverage floor recorded in
# ROADMAP.md. Calibrate in Gate 9 and bump WINDOW_COVERAGE_VERSION when it moves.
MINIMUM_WINDOWS_FOR_FRACTION = 20

# Stated verbatim wherever a coverage value is surfaced.
COVERAGE_DISCLAIMER = (
    "Human repertoire window coverage is the fraction of this sequence's "
    "9-residue windows found in a reference set of human antibody peptides. It "
    "is not a measure of humanness, human origin, immunogenicity, tolerability, "
    "or safety. A camelid VHH scores close to a humanized antibody under this "
    "metric, so it does not distinguish human from non-human sequences."
)


@dataclass(frozen=True)
class WindowObservation:
    """One window and whether it was found. An observation, not an assessment."""

    # 1-based window index, matching promb's CLI output so results are comparable.
    window_index: int
    # 0-based offset of the window's first residue in the submitted sequence.
    sequence_offset: int
    peptide: str
    observed: bool


@dataclass
class WindowCoverageResult:
    """
    Coverage for one component, plus everything needed to reproduce or discard it.

    `coverage` is populated only when `status` is AVAILABLE. For every other
    status it is None -- deliberately not 0.0, which would be indistinguishable
    from a sequence whose windows were all genuinely absent.
    """

    status: EvidenceStatus
    coverage: Optional[float]

    # Raw counts, always populated when they could be computed at all. These are
    # the reproducible measurement; `coverage` is their quotient.
    windows_total: int = 0
    windows_observed: int = 0

    # Provenance. Every field here is required to reproduce the number later.
    window_length: int = WINDOW_LENGTH
    reference_set_id: str = REFERENCE_SET_ID
    reference_set_file: str = REFERENCE_SET_FILE
    reference_set_description: str = REFERENCE_SET_DESCRIPTION
    reference_library_version: Optional[str] = None
    engine_version: str = WINDOW_COVERAGE_VERSION

    # Why the status is what it is. Populated for every non-AVAILABLE status;
    # a status without a reason is not actionable by whoever reads the profile.
    reason: Optional[str] = None

    # Per-window detail, so a reported fraction can always be audited back to the
    # windows that produced it. Empty when the caller opts out.
    observations: list[WindowObservation] = field(default_factory=list)

    limitations: list[str] = field(default_factory=list)
    disclaimer: str = COVERAGE_DISCLAIMER

    @property
    def windows_not_observed(self) -> int:
        return self.windows_total - self.windows_observed


class ReferenceSetUnavailable(RuntimeError):
    """Raised internally when the reference set cannot be loaded; surfaced as a status."""


# The reference set is ~5.7M peptides and takes a few seconds to load into memory,
# so it is loaded once and reused. Cached as a module-level singleton rather than
# rebuilt per call; `_load_error` is cached too, so a missing dependency does not
# cost a retry on every component in a corpus-wide run.
_reference_db = None
_load_error: Optional[str] = None


def load_reference_set(force_reload: bool = False):
    """
    Load and memoise the human repertoire reference set.

    Raises ReferenceSetUnavailable if promb is not installed or the set cannot be
    read. Callers that want a status rather than an exception should use
    `compute_window_coverage`, which converts this into DATABASE_UNAVAILABLE.
    """
    global _reference_db, _load_error

    if force_reload:
        _reference_db, _load_error = None, None

    if _reference_db is not None:
        return _reference_db
    if _load_error is not None:
        raise ReferenceSetUnavailable(_load_error)

    try:
        from promb import init_db

        _reference_db = init_db(REFERENCE_SET_ID)
    except Exception as exc:
        _load_error = f"{type(exc).__name__}: {exc}"
        raise ReferenceSetUnavailable(_load_error) from exc

    return _reference_db


def _reference_library_version() -> Optional[str]:
    try:
        from promb import __version__

        return str(__version__)
    except Exception:
        return None


def chop_windows(sequence: str, window_length: int = WINDOW_LENGTH) -> list[str]:
    """Every overlapping window of `window_length`, in order. Empty if too short."""
    if len(sequence) < window_length:
        return []
    return [sequence[i : i + window_length] for i in range(len(sequence) - window_length + 1)]


def compute_window_coverage(
    sequence: str,
    database=None,
    record_observations: bool = True,
) -> WindowCoverageResult:
    """
    Compute human repertoire window coverage for one component sequence.

    Always returns a result. Every failure mode is expressed as an EvidenceStatus
    with a reason rather than as an exception or a substitute number, so an
    evidence profile can record what happened instead of losing the slot.

    Statuses returned:
      AVAILABLE              coverage computed over enough windows to report
      NOT_APPLICABLE         sequence shorter than the 9-residue window
      INSUFFICIENT_EVIDENCE  fewer windows than MINIMUM_WINDOWS_FOR_FRACTION
      UNSUPPORTED_COMPONENT  non-standard residues present, or empty input
      DATABASE_UNAVAILABLE   reference set could not be loaded
      ANALYSIS_FAILED        lookup raised unexpectedly
    """
    library_version = _reference_library_version()

    def result(status: EvidenceStatus, **kwargs) -> WindowCoverageResult:
        kwargs.setdefault("coverage", None)
        return WindowCoverageResult(
            status=status, reference_library_version=library_version, **kwargs
        )

    cleaned = sequence.strip().upper()

    if not cleaned:
        return result(
            EvidenceStatus.UNSUPPORTED_COMPONENT,
            reason="Empty sequence; there is nothing to window.",
        )

    # Checked before length, because a sequence that is both short and ambiguous
    # should report the ambiguity -- that is the fact the caller has to act on.
    unknown = sorted(set(cleaned) - set(STANDARD_AMINO_ACIDS))
    if unknown:
        return result(
            EvidenceStatus.UNSUPPORTED_COMPONENT,
            reason=(
                f"Sequence contains residues outside the 20 standard amino acids: "
                f"{unknown}. Coverage is not computed, because a window containing "
                f"an unresolvable residue cannot match the reference set and would "
                f"lower the fraction -- reporting absence of evidence as a low "
                f"score. The sequence is not modified to make it scoreable; the "
                f"caller must decide how to handle the residue and record that."
            ),
        )

    windows = chop_windows(cleaned)
    if not windows:
        return result(
            EvidenceStatus.NOT_APPLICABLE,
            reason=(
                f"Sequence is {len(cleaned)} residues, shorter than the "
                f"{WINDOW_LENGTH}-residue window, so no window exists and the "
                f"metric is undefined. This is expected for most CDRs and is not "
                f"a failure."
            ),
        )

    try:
        db = database if database is not None else load_reference_set()
    except ReferenceSetUnavailable as exc:
        return result(
            EvidenceStatus.DATABASE_UNAVAILABLE,
            windows_total=len(windows),
            reason=(
                f"Human repertoire reference set could not be loaded ({exc}). "
                f"Coverage is unset rather than defaulted; this status must not "
                f"be read as low coverage."
            ),
        )

    try:
        observations = [
            WindowObservation(
                window_index=index + 1,
                sequence_offset=index,
                peptide=peptide,
                observed=bool(db.contains(peptide)),
            )
            for index, peptide in enumerate(windows)
        ]
    except Exception as exc:
        return result(
            EvidenceStatus.ANALYSIS_FAILED,
            windows_total=len(windows),
            reason=(
                f"Reference set lookup failed ({type(exc).__name__}: {exc}). "
                f"Coverage is unset; this is a tooling failure and carries no "
                f"biological meaning."
            ),
        )

    observed_count = sum(1 for obs in observations if obs.observed)
    kept = observations if record_observations else []

    if len(windows) < MINIMUM_WINDOWS_FOR_FRACTION:
        return result(
            EvidenceStatus.INSUFFICIENT_EVIDENCE,
            windows_total=len(windows),
            windows_observed=observed_count,
            observations=kept,
            reason=(
                f"Only {len(windows)} windows (minimum "
                f"{MINIMUM_WINDOWS_FOR_FRACTION} to report a fraction). The raw "
                f"counts are returned, but a ratio over so few windows has no "
                f"length-matched null to be interpreted against (SPEC §21). The "
                f"minimum is a provisional default, not a calibrated threshold."
            ),
        )

    return WindowCoverageResult(
        status=EvidenceStatus.AVAILABLE,
        coverage=observed_count / len(windows),
        windows_total=len(windows),
        windows_observed=observed_count,
        reference_library_version=library_version,
        observations=kept,
    )
