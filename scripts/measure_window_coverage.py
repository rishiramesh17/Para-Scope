"""
Measure human repertoire window coverage across the whole corpus (Gate 6).

Produces the per-component-type distributions that Gate 9 needs in order to turn
a raw coverage fraction into a calibrated percentile. Until those distributions
exist, a coverage value has nothing to be interpreted against.

Reports the status breakdown as prominently as the values, because the count of
components the metric could not score is itself a result -- most CDRs are shorter
than the 9-residue window and are correctly NOT_APPLICABLE, and that must be
visible rather than absorbed into a denominator.

Usage:
    .venv/bin/python scripts/measure_window_coverage.py
"""

import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from antigens.models import EvidenceStatus  # noqa: E402
from antigens.persistence import load_corpus_snapshot  # noqa: E402
from antigens.window_coverage import (  # noqa: E402
    MINIMUM_WINDOWS_FOR_FRACTION,
    REFERENCE_SET_ID,
    WINDOW_COVERAGE_VERSION,
    WINDOW_LENGTH,
    compute_window_coverage,
    load_reference_set,
)

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")


def percentiles(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    if len(ordered) < 2:
        return {"p05": ordered[0], "p25": ordered[0], "p50": ordered[0],
                "p75": ordered[0], "p95": ordered[0]}
    cuts = statistics.quantiles(ordered, n=100, method="inclusive")
    return {
        "p05": cuts[4],
        "p25": cuts[24],
        "p50": cuts[49],
        "p75": cuts[74],
        "p95": cuts[94],
    }


def main() -> int:
    if not SNAPSHOT_PATH.exists():
        print(f"Corpus snapshot not found at {SNAPSHOT_PATH}.")
        print("Rebuild it first: .venv/bin/python scripts/build_combined_corpus_cache.py")
        return 1

    print(f"engine version   {WINDOW_COVERAGE_VERSION}")
    print(f"window length    {WINDOW_LENGTH}")
    print(f"reference set    {REFERENCE_SET_ID}")
    print(f"minimum windows  {MINIMUM_WINDOWS_FOR_FRACTION} (provisional)")
    print()

    corpus, _stats, _manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    db = load_reference_set()
    print(f"corpus loaded: {len(corpus):,} components\n")

    values_by_type: dict[str, list[float]] = defaultdict(list)
    statuses_by_type: dict[str, Counter] = defaultdict(Counter)

    started = time.time()
    for record in corpus._components.values():
        result = compute_window_coverage(
            record.sequence, database=db, record_observations=False
        )
        key = record.component_type.value
        statuses_by_type[key][result.status.value] += 1
        if result.status is EvidenceStatus.AVAILABLE:
            values_by_type[key].append(result.coverage)
    elapsed = time.time() - started

    print(f"scored {len(corpus):,} components in {elapsed:.1f}s "
          f"({len(corpus) / elapsed:,.0f}/s)\n")

    print("=" * 78)
    print("COVERAGE DISTRIBUTION, by component type (AVAILABLE only)")
    print("=" * 78)
    header = f"{'component':<16}{'n':>7}{'p05':>8}{'p25':>8}{'p50':>8}{'p75':>8}{'p95':>8}"
    print(header)
    print("-" * 78)
    for key in sorted(values_by_type, key=lambda k: -len(values_by_type[k])):
        values = values_by_type[key]
        pct = percentiles(values)
        print(f"{key:<16}{len(values):>7,}"
              f"{pct['p05']:>8.3f}{pct['p25']:>8.3f}{pct['p50']:>8.3f}"
              f"{pct['p75']:>8.3f}{pct['p95']:>8.3f}")

    print()
    print("=" * 78)
    print("STATUS BREAKDOWN, by component type")
    print("=" * 78)
    all_statuses = Counter()
    for key in sorted(statuses_by_type):
        counts = statuses_by_type[key]
        all_statuses.update(counts)
        summary = "  ".join(f"{status}={count:,}" for status, count in counts.most_common())
        print(f"{key:<16}{summary}")

    print()
    print("corpus totals:")
    for status, count in all_statuses.most_common():
        print(f"  {status:<24}{count:>8,}  ({count / len(corpus):.1%})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
