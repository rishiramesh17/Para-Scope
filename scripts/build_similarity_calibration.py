"""
Gate 9: add the component-similarity reference distribution to an existing calibration.

Separate from build_calibration.py and resumable so the expensive germline pass is
not re-paid: this loads the artifact that script wrote, adds only the distributions
it lacks, and writes it back under the same calibration version.

What it measures
----------------
For each sampled component, the best global_identity it achieves against the
*therapeutic-bearing* component population, excluding itself. That is the
measurement Gate 7's `calibrated_component_similarity` places a match in: "this
candidate's VH is 0.91 identical to its closest therapeutic VH" means nothing
until you know that 0.91 sits at, say, the 73rd percentile of how close an
arbitrary VH gets to the nearest therapeutic VH.

Self-exclusion matters. A corpus component that *is* a therapeutic would
otherwise score 1.0 against itself, and a population of 1.0s would make every
genuine near-match look unremarkable.

Cost
----
One alignment per therapeutic component per sampled component, so the sample is
deliberately small and its size is recorded on each distribution where
MINIMUM_SAMPLE_SIZE can see it. Exhaustive search (shortlist_limit=None) is used
within the sample: a shortlist is ~96% complete, and a *calibration* built from
systematically-incomplete best-match values would bias every percentile derived
from it.
"""

import json
import random
import time
from datetime import datetime, timezone
from pathlib import Path

from para_scope.calibration import (
    LENGTH_STRATIFIED_CLASSES,
    build_distribution,
    length_band,
    load_calibration,
    save_calibration,
    stratum_key,
)
from para_scope.indexes.approximate import ELIGIBLE_TYPES, ApproximateIndex
from para_scope.models import ComponentType
from para_scope.persistence import load_corpus_snapshot
from para_scope.therapeutic_precedent import SIMILARITY_MEASUREMENT

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
CALIBRATION_PATH = Path("data/cache/calibration_v1.json")

SAMPLE_PER_CLASS = 120
SAMPLE_SEED = 20261004

POPULATION = (
    "Best global_identity against Thera-SAbDab therapeutic components (self "
    "excluded), over a deterministic sample of SAbDab2 + Thera-SAbDab components. "
    "Antibodies selected for crystallisation or clinical development; NOT a sample "
    "of natural antibody repertoire."
)


def main() -> None:
    t0 = time.time()
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    snapshots = list(manifest.source_snapshots)
    print(f"Loaded corpus: {len(corpus)} components")

    therapeutic_ids = set(corpus._component_to_therapeutics)  # noqa: SLF001 -- offline builder
    print(f"Therapeutic-bearing components: {len(therapeutic_ids)}")

    # One index per component type holding only therapeutic components: the
    # population a match is scored against must be the therapeutic one, not the
    # whole corpus.
    by_type: dict[ComponentType, ApproximateIndex] = {}
    counts: dict[ComponentType, int] = {}
    for record_id in sorted(therapeutic_ids):
        record = corpus.get(record_id)
        if record is None or record.component_type not in ELIGIBLE_TYPES:
            continue
        by_type.setdefault(record.component_type, ApproximateIndex()).add(
            record.component_type, record.record_id, record.sequence
        )
        counts[record.component_type] = counts.get(record.component_type, 0) + 1
    print("therapeutic index sizes:", {k.value: v for k, v in sorted(counts.items(), key=lambda kv: kv[0].value)})

    existing = load_calibration(CALIBRATION_PATH)
    print(f"Existing calibration: {len(existing)} distributions, built {existing.built_at}")

    rng = random.Random(SAMPLE_SEED)
    new_dists = []

    for component_type, index in sorted(by_type.items(), key=lambda kv: kv[0].value):
        candidates = sorted(
            (r for r in corpus._components.values()  # noqa: SLF001
             if r.component_type == component_type),
            key=lambda r: r.record_id,
        )
        if not candidates:
            continue
        sample = (
            candidates if len(candidates) <= SAMPLE_PER_CLASS
            else rng.sample(candidates, SAMPLE_PER_CLASS)
        )

        t1 = time.time()
        by_stratum: dict[str, list[float]] = {}
        strat_meta: dict[str, dict[str, str]] = {}
        for record in sample:
            hits = index.search(
                component_type, record.sequence, top_k=5, shortlist_limit=None
            )
            best = None
            for hit in hits:
                if hit.record_id == record.record_id:
                    continue  # self-match: see module docstring
                best = hit.global_identity
                break
            if best is None:
                continue
            length = (
                len(record.sequence)
                if component_type.value in LENGTH_STRATIFIED_CLASSES else None
            )
            key = stratum_key(component_type.value, length)
            by_stratum.setdefault(key, []).append(float(best))
            strat_meta[key] = {"length_band": length_band(length)} if length is not None else {}

        for key, values in sorted(by_stratum.items()):
            new_dists.append(
                build_distribution(
                    measurement=SIMILARITY_MEASUREMENT,
                    component_class=component_type.value,
                    stratum=key,
                    stratification=strat_meta[key],
                    values=values,
                    population_description=POPULATION,
                    source_snapshots=snapshots,
                )
            )
        total = sum(len(v) for v in by_stratum.values())
        print(f"  {component_type.value:10} sampled {len(sample):4} -> {total:4} obs "
              f"in {len(by_stratum)} strata ({time.time()-t1:.1f}s)")

    # Merge: keep every existing distribution, add the new ones.
    merged = [
        existing.get(m, s) for (m, s) in sorted(existing._distributions)  # noqa: SLF001
    ]
    merged += new_dists
    built_at = datetime.now(timezone.utc).isoformat()
    save_calibration(merged, built_at, snapshots, CALIBRATION_PATH)
    print(f"\nWrote {len(merged)} distributions ({len(new_dists)} new) to {CALIBRATION_PATH}")
    print(f"total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
