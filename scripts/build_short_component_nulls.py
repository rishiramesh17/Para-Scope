"""
Gate 11 / SPEC §21: build empirical null distributions for short components.

Adds `null.short_component_identity` strata to the existing calibration artifact,
one per (CDR type, length band). Resumable like the similarity builder: the
expensive germline and similarity passes are not re-paid.

The null is pairs of **distinct real components** of the same type and length band,
not random strings. Real CDRs share compositional biases that random sequence does
not, and a null without that bias would make every real match look spectacular --
it would not be a null for the question being asked.

Pairs are drawn from distinct parent antibodies where that is knowable, so the
background measures *unrelated* components rather than two deposits of one clone.
Comparing an antibody against its own near-duplicate would put genuine matches in
the null and flatten the distribution toward the observed range.
"""

import json
import random
import statistics
from datetime import datetime, timezone
from pathlib import Path

from para_scope.calibration import (
    build_distribution,
    length_band,
    load_calibration,
    save_calibration,
)
from para_scope.models import ComponentType
from para_scope.persistence import load_corpus_snapshot
from para_scope.short_component import (
    MINIMUM_NULL_COMPARISONS,
    SHORT_TYPES,
    position_aware_similarity,
)

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
CALIBRATION_PATH = Path("data/cache/calibration_v1.json")

NULL_MEASUREMENT = "null.short_component_identity"
PAIRS_PER_STRATUM = 1200
SAMPLE_SEED = 20261005

POPULATION = (
    "Pairs of distinct real components of the same type and length band from "
    "SAbDab2 + Thera-SAbDab, drawn from different parent records. This is a "
    "background for chance similarity between unrelated components of this type "
    "and length; it is NOT a sample of natural antibody repertoire."
)


def main() -> None:
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    snapshots = list(manifest.source_snapshots)
    print(f"corpus {len(corpus)} components (schema {manifest.corpus_schema_version})\n")

    # Group short components by (type, length band).
    strata: dict[tuple[ComponentType, str], list[tuple[str, str]]] = {}
    for record in corpus._components.values():  # noqa: SLF001 -- offline builder
        if record.component_type not in SHORT_TYPES:
            continue
        band = length_band(len(record.sequence))
        strata.setdefault((record.component_type, band), []).append(
            (record.record_id, record.sequence)
        )

    existing = load_calibration(CALIBRATION_PATH)
    print(f"existing calibration: {len(existing)} distributions")

    rng = random.Random(SAMPLE_SEED)
    new_dists = []
    skipped = []

    for (component_type, band), members in sorted(
        strata.items(), key=lambda kv: (kv[0][0].value, kv[0][1])
    ):
        stratum = f"{component_type.value}|len={band}"
        if len(members) < 2:
            skipped.append((stratum, len(members), "fewer than two components"))
            continue

        values: list[float] = []
        attempts = 0
        max_attempts = PAIRS_PER_STRATUM * 6
        seen: set[tuple[str, str]] = set()
        while len(values) < PAIRS_PER_STRATUM and attempts < max_attempts:
            attempts += 1
            left, right = rng.sample(members, 2)
            # Distinct components only; identical sequences would seed the null
            # with 1.0 and make a real exact match look ordinary.
            if left[1] == right[1]:
                continue
            key = tuple(sorted((left[0], right[0])))
            if key in seen:
                continue
            seen.add(key)
            similarity, compatible, _ = position_aware_similarity(left[1], right[1])
            if compatible:
                values.append(similarity)

        if len(values) < MINIMUM_NULL_COMPARISONS:
            skipped.append((stratum, len(values), "too few comparable pairs"))
            continue

        new_dists.append(
            build_distribution(
                measurement=NULL_MEASUREMENT,
                component_class=component_type.value,
                stratum=stratum,
                stratification={"length_band": band},
                values=values,
                population_description=POPULATION,
                source_snapshots=snapshots,
            )
        )
        print(f"  {stratum:26} n={len(values):5} "
              f"median={statistics.median(values):.3f} "
              f"p95={sorted(values)[int(len(values)*0.95)]:.3f} "
              f"max={max(values):.3f}")

    if skipped:
        print(f"\nskipped {len(skipped)} strata:")
        for stratum, count, why in skipped[:10]:
            print(f"  {stratum:26} {count:5} -- {why}")

    merged = [
        existing.get(m, s) for (m, s) in sorted(existing._distributions)  # noqa: SLF001
    ] + new_dists
    built_at = datetime.now(timezone.utc).isoformat()
    save_calibration(merged, built_at, snapshots, CALIBRATION_PATH)
    print(f"\nwrote {len(merged)} distributions ({len(new_dists)} new nulls) "
          f"to {CALIBRATION_PATH}")


if __name__ == "__main__":
    main()
