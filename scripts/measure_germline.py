"""
Measure human germline V-gene identity across the whole corpus (Gate 9, milestone 1).

Produces the per-component-type distributions Database-Relative Novelty needs a
baseline from. Reports the status breakdown alongside the values, because VHH
components are UNSUPPORTED_COMPONENT by design (no camelid reference set) and
that must be visible, not silently excluded from a denominator that looks complete.

Exhaustive scoring (every unique VH/VL against all ~350 reference alleles) is
too slow for a first pass -- runs to a deterministic sample per type instead
(every Nth component, seeded by iteration order, not randomized). Full-corpus
percentiles are Gate 9-proper work, alongside the statistical calibration that
turns a raw identity into something interpretable; this pass exists to confirm
the shape of the distribution, not to finalize it.

Usage:
    .venv/bin/python scripts/measure_germline.py
"""

import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

SAMPLE_STRIDE = 8  # take every 8th unique component per type

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from para_scope.germline import (  # noqa: E402
    GERMLINE_VERSION,
    REFERENCE_SET_ID,
    REFERENCE_SET_VERSION,
    compute_germline_assignment,
    load_reference_index,
)
from para_scope.models import ComponentType, EvidenceStatus  # noqa: E402
from para_scope.persistence import load_corpus_snapshot  # noqa: E402

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
MEASURED_TYPES = (ComponentType.VH, ComponentType.VL_KAPPA, ComponentType.VL_LAMBDA, ComponentType.VHH)


def percentiles(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    if len(ordered) < 2:
        return {"p05": ordered[0], "p25": ordered[0], "p50": ordered[0],
                "p75": ordered[0], "p95": ordered[0]}
    cuts = statistics.quantiles(ordered, n=100, method="inclusive")
    return {"p05": cuts[4], "p25": cuts[24], "p50": cuts[49], "p75": cuts[74], "p95": cuts[94]}


def main() -> None:
    t0 = time.time()
    corpus, stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    print(f"Snapshot loaded in {time.time() - t0:.2f}s")

    index = load_reference_index()
    print(f"Reference index loaded: {REFERENCE_SET_ID} v{REFERENCE_SET_VERSION}, "
          f"{len(index)} germline V-gene alleles, engine v{GERMLINE_VERSION}")

    identities: dict[ComponentType, list[float]] = defaultdict(list)
    gene_calls: dict[ComponentType, Counter] = defaultdict(Counter)
    status_counts: dict[ComponentType, Counter] = defaultdict(Counter)

    seen: set[str] = set()
    per_type_counter: dict[ComponentType, int] = defaultdict(int)
    t1 = time.time()
    n = 0
    for record in corpus._components.values():
        if record.component_type not in MEASURED_TYPES:
            continue
        if record.sequence_checksum in seen:
            continue
        seen.add(record.sequence_checksum)

        per_type_counter[record.component_type] += 1
        if per_type_counter[record.component_type] % SAMPLE_STRIDE != 0:
            continue
        n += 1

        result = compute_germline_assignment(record.component_type, record.sequence)
        status_counts[record.component_type][result.status] += 1
        if result.status == EvidenceStatus.AVAILABLE:
            identities[record.component_type].append(result.nearest_germline_identity)
            gene_calls[record.component_type][result.nearest_gene_call] += 1

    print(f"Scored {n} sampled components (stride={SAMPLE_STRIDE}) in {time.time() - t1:.2f}s\n")

    for component_type in MEASURED_TYPES:
        total = sum(status_counts[component_type].values())
        corpus_total = per_type_counter[component_type]
        print(f"--- {component_type.value} (sampled n={total} of {corpus_total} unique) ---")
        for status, count in sorted(status_counts[component_type].items(), key=lambda kv: -kv[1]):
            print(f"  {status.value}: {count} ({count / total:.1%})")
        values = identities[component_type]
        if values:
            p = percentiles(values)
            print(f"  identity percentiles (n={len(values)}): "
                  f"p05={p['p05']:.3f} p25={p['p25']:.3f} p50={p['p50']:.3f} "
                  f"p75={p['p75']:.3f} p95={p['p95']:.3f}")
            top_genes = gene_calls[component_type].most_common(5)
            print(f"  most common nearest gene: {top_genes}")
        print()


if __name__ == "__main__":
    main()
