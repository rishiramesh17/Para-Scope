"""
Empirical reference distributions and calibrated percentiles (Version 2, Gate 9).

SPEC S19 Layer 1 is a raw measurement -- "this VH is 0.73 identical to its nearest
germline gene". Layer 2 turns that into a *relative* statement -- "0.73 sits at the
41st percentile of VH germline identity in this reference population". Gates 6, 7,
8 and the germline module all emit Layer 1 numbers today and none of them can say
whether a number is ordinary or unusual, because nothing establishes what ordinary
looks like. That is what this module is.

What a percentile here does and does not mean
---------------------------------------------
A percentile is a statement about a **population**, never about a molecule. The
reference population available to this project is structure-derived: SAbDab2 is
antibodies someone chose to crystallise, and Thera-SAbDab is antibodies someone
chose to develop as drugs. Neither is a sample of natural repertoire, and nothing
makes them one.

So every percentile this module emits names its population explicitly and is only
interpretable against it. "91st percentile of solved structures" is defensible.
"91st percentile of human antibodies" is not, and the distribution id, the
population description, and the sample size exist so a reader cannot mistake one
for the other. SPEC S19 requires exactly this: reference distribution, component
class, database snapshot, stratification variables, sample size, calibration
version -- all six travel with the number.

A percentile is also not a grade, a score, or a ranking of quality. An unusual
value is unusual. It is not worse, riskier, or less likely to work (CLAUDE.md,
prohibited claims).

Why distributions are built offline
-----------------------------------
Scoring code never queries live data (CLAUDE.md). Distributions are computed once
by `scripts/build_calibration.py` from a checksummed corpus snapshot, stored as a
versioned artifact, and loaded read-only at scoring time. A corpus update produces
a new calibration version rather than silently changing old percentiles -- which
is the whole reason a prior analysis stays reproducible.

Stratification
--------------
Like is compared with like, so a distribution is keyed by component class, and for
length-sensitive components by a length band as well (SPEC S19: "CDR-H3 compared
with CDR-H3 of similar length"). A CDR-H3 of 22 residues has no business being
percentiled against a population dominated by 12-residue CDR-H3s; length drives
both the measurement and its chance distribution.

Insufficient data is a status, not a wider bin
-----------------------------------------------
When a stratum holds fewer than MINIMUM_SAMPLE_SIZE observations, the result is
INSUFFICIENT_EVIDENCE with the sample size attached. It does not silently widen to
a neighbouring bin, fall back to the unstratified distribution, or return 50th
percentile. Each of those would answer a question nobody asked while looking like
an answer to the one that was.
"""

import json
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from para_scope.models import EvidenceStatus

# Bump on any change to how a percentile is computed, how strata are defined, or
# which population a distribution is built from. Stored on every emitted
# percentile so a stored score always names the calibration that produced it.
#   1.0.0 -- initial: empirical percentiles, length-banded strata
CALIBRATION_VERSION = "1.0.0"

# Fewest observations a stratum may hold and still yield a percentile.
#
# PROVISIONAL, like MIN_CONSTANT_REGION_LENGTH and the 0.70 coverage floor. The
# figure is chosen so a reported percentile cannot be an artifact of a handful of
# points: with 30 observations a single outlier moves a percentile by at most
# ~3.3 points. It is not derived from a power calculation, and is recorded here as
# a judgement rather than a result.
MINIMUM_SAMPLE_SIZE = 30

# Length bands for length-sensitive components, in residues. Chosen to keep CDR-H3
# -- by far the most length-variable component -- in bands that hold real
# populations rather than to be round numbers. Bands are closed on the left.
LENGTH_BANDS: tuple[tuple[int, int], ...] = (
    (0, 8), (9, 11), (12, 14), (15, 17), (18, 21), (22, 26), (27, 1_000),
)

# Component classes whose measurements are banded by length. Whole variable
# domains vary little in length and are not banded; CDRs vary enormously.
LENGTH_STRATIFIED_CLASSES = frozenset({
    "CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3",
})

POPULATION_DISCLAIMER = (
    "Percentiles are relative to the named reference population, which is "
    "structure- and therapeutic-derived (antibodies selected for crystallisation "
    "or clinical development). It is not a sample of natural antibody repertoire. "
    "An unusual value is unusual within this population; it is not evidence that a "
    "candidate is worse, riskier, or less likely to succeed."
)


def length_band(length: int) -> str:
    """The band label for a residue length, e.g. '12-14' or '27+'."""
    for low, high in LENGTH_BANDS:
        if low <= length <= high:
            return f"{low}-{high}" if high < 1_000 else f"{low}+"
    raise ValueError(f"length {length} falls in no band")


def stratum_key(component_class: str, length: Optional[int]) -> str:
    """Stratum id for a component class and (where it matters) a length."""
    if component_class in LENGTH_STRATIFIED_CLASSES:
        if length is None:
            raise ValueError(
                f"{component_class} is length-stratified; a length is required to "
                f"select its reference distribution"
            )
        return f"{component_class}|len={length_band(length)}"
    return component_class


@dataclass(frozen=True)
class ReferenceDistribution:
    """One empirical distribution, with everything SPEC S19 requires to name it."""

    measurement: str            # e.g. "nearest_germline_identity"
    component_class: str        # e.g. "CDR-H3"
    stratum: str                # stratum_key(...) -- component class plus any bands
    stratification: dict[str, str]
    sample_size: int
    # Sorted observations. Stored in full rather than as summary quantiles so a
    # percentile is reproducible exactly, and so a later version can recompute a
    # different statistic from the same population without rebuilding it.
    values: tuple[float, ...]
    population_description: str
    source_snapshots: tuple[dict, ...]
    calibration_version: str = CALIBRATION_VERSION

    def percentile_of(self, value: float) -> float:
        """
        Empirical percentile of `value`, 0-100.

        Uses the midpoint of the strictly-below and at-or-below ranks, so ties sit
        at the centre of the tied block rather than at one edge. With many
        identical values -- common for saturated measurements like framework window
        coverage -- taking either edge would report a confidently extreme
        percentile for a value that is in fact the single most ordinary one.
        """
        below = bisect_left(self.values, value)
        at_or_below = bisect_right(self.values, value)
        return 100.0 * ((below + at_or_below) / 2) / len(self.values)

    def summary(self) -> dict[str, float]:
        """Quantiles of the population itself, for display alongside a percentile."""
        n = len(self.values)
        def q(p: float) -> float:
            return self.values[min(n - 1, max(0, int(round(p * (n - 1)))))]
        return {
            "min": self.values[0], "p05": q(0.05), "p25": q(0.25),
            "median": q(0.50), "p75": q(0.75), "p95": q(0.95), "max": self.values[-1],
        }


@dataclass
class CalibratedPercentile:
    """A Layer 2 result: a percentile, or an explicit reason there isn't one."""

    status: EvidenceStatus
    measurement: str
    component_class: str
    raw_value: Optional[float] = None
    percentile: Optional[float] = None

    stratum: Optional[str] = None
    stratification: dict[str, str] = field(default_factory=dict)
    sample_size: Optional[int] = None
    population_description: Optional[str] = None
    population_summary: dict[str, float] = field(default_factory=dict)
    source_snapshots: tuple[dict, ...] = ()
    calibration_version: str = CALIBRATION_VERSION

    reason: Optional[str] = None
    disclaimer: str = POPULATION_DISCLAIMER


class CalibrationSet:
    """A loaded, read-only collection of reference distributions."""

    def __init__(
        self,
        distributions: dict[tuple[str, str], ReferenceDistribution],
        built_at: str,
        source_snapshots: tuple[dict, ...],
        calibration_version: str = CALIBRATION_VERSION,
    ):
        self._distributions = distributions
        self.built_at = built_at
        self.source_snapshots = source_snapshots
        self.calibration_version = calibration_version

    def __len__(self) -> int:
        return len(self._distributions)

    def strata(self, measurement: str) -> list[str]:
        return sorted(s for (m, s) in self._distributions if m == measurement)

    def get(self, measurement: str, stratum: str) -> Optional[ReferenceDistribution]:
        return self._distributions.get((measurement, stratum))

    def percentile(
        self,
        measurement: str,
        component_class: str,
        value: Optional[float],
        length: Optional[int] = None,
    ) -> CalibratedPercentile:
        """
        Percentile of `value` within its stratum, or a status saying why not.

        A missing `value` returns NOT_APPLICABLE rather than a percentile: there is
        nothing to locate in a distribution, and 0 is not the same as absent.
        """
        result = CalibratedPercentile(
            status=EvidenceStatus.NOT_FOUND,
            measurement=measurement,
            component_class=component_class,
            raw_value=value,
            calibration_version=self.calibration_version,
        )

        if value is None:
            result.status = EvidenceStatus.NOT_APPLICABLE
            result.reason = (
                "No measurement was supplied, so there is nothing to place in a "
                "reference distribution. This is not a percentile of zero."
            )
            return result

        try:
            stratum = stratum_key(component_class, length)
        except ValueError as exc:
            result.status = EvidenceStatus.UNSUPPORTED_COMPONENT
            result.reason = str(exc)
            return result

        result.stratum = stratum
        distribution = self.get(measurement, stratum)

        if distribution is None:
            result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
            result.reason = (
                f"No reference distribution exists for measurement {measurement!r} "
                f"in stratum {stratum!r} at calibration version "
                f"{self.calibration_version}. The stratum is not widened and no "
                f"unstratified fallback is used -- either would answer a different "
                f"question than the one asked."
            )
            return result

        if distribution.sample_size < MINIMUM_SAMPLE_SIZE:
            result.status = EvidenceStatus.INSUFFICIENT_EVIDENCE
            result.sample_size = distribution.sample_size
            result.stratification = dict(distribution.stratification)
            result.population_description = distribution.population_description
            result.reason = (
                f"Stratum {stratum!r} holds {distribution.sample_size} observations, "
                f"below the minimum of {MINIMUM_SAMPLE_SIZE} for a reportable "
                f"percentile."
            )
            return result

        result.status = EvidenceStatus.AVAILABLE
        result.percentile = distribution.percentile_of(value)
        result.sample_size = distribution.sample_size
        result.stratification = dict(distribution.stratification)
        result.population_description = distribution.population_description
        result.population_summary = distribution.summary()
        result.source_snapshots = distribution.source_snapshots
        return result


def save_calibration(
    distributions: list[ReferenceDistribution],
    built_at: str,
    source_snapshots: list[dict],
    path: Path,
) -> None:
    """Write a calibration artifact as JSON (inspectable, diffable, not pickled)."""
    payload = {
        "calibration_version": CALIBRATION_VERSION,
        "built_at": built_at,
        "minimum_sample_size": MINIMUM_SAMPLE_SIZE,
        "source_snapshots": source_snapshots,
        "population_disclaimer": POPULATION_DISCLAIMER,
        "distributions": [
            {
                "measurement": d.measurement,
                "component_class": d.component_class,
                "stratum": d.stratum,
                "stratification": d.stratification,
                "sample_size": d.sample_size,
                "population_description": d.population_description,
                "values": list(d.values),
            }
            for d in distributions
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1))


def load_calibration(path: Path) -> CalibrationSet:
    """Load a calibration artifact. Raises FileNotFoundError if absent."""
    payload = json.loads(path.read_text())
    snapshots = tuple(payload.get("source_snapshots", ()))
    distributions: dict[tuple[str, str], ReferenceDistribution] = {}
    for entry in payload["distributions"]:
        d = ReferenceDistribution(
            measurement=entry["measurement"],
            component_class=entry["component_class"],
            stratum=entry["stratum"],
            stratification=entry["stratification"],
            sample_size=entry["sample_size"],
            values=tuple(entry["values"]),
            population_description=entry["population_description"],
            source_snapshots=snapshots,
            calibration_version=payload["calibration_version"],
        )
        distributions[(d.measurement, d.stratum)] = d
    return CalibrationSet(
        distributions=distributions,
        built_at=payload["built_at"],
        source_snapshots=snapshots,
        calibration_version=payload["calibration_version"],
    )


def build_distribution(
    measurement: str,
    component_class: str,
    stratum: str,
    stratification: dict[str, str],
    values: list[float],
    population_description: str,
    source_snapshots: list[dict],
) -> ReferenceDistribution:
    """Construct a distribution from raw observations, sorting them once."""
    return ReferenceDistribution(
        measurement=measurement,
        component_class=component_class,
        stratum=stratum,
        stratification=stratification,
        sample_size=len(values),
        values=tuple(sorted(values)),
        population_description=population_description,
        source_snapshots=tuple(source_snapshots),
    )
