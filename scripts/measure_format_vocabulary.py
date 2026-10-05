"""
Gate 10.5: measure the construct-format vocabulary against the real source column.

A controlled vocabulary over free text is only as good as its coverage, and the
number that matters is not how many *distinct* strings parse but how many *rows*
do -- 70% of this column is the single value "Whole mAb", so distinct-value
coverage would understate it badly and row coverage is the honest figure. Both are
reported.

Also reports how many rows become uncomparable against a declared format, because
a vocabulary that parses everything into buckets that cannot be compared has
bought nothing.

Writes reports/format_vocabulary.json.
"""

import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from para_scope.construct_format import (
    AXIS_WEIGHTS,
    FORMAT_VOCABULARY_VERSION,
    Architecture,
    Modification,
    Specificity,
    compute_format_relevance,
    parse_format,
)
from para_scope.models import EvidenceStatus

SOURCE = Path("data/snapshots/thera_sabdab/TheraSAbDab_SeqStruc_OnlineDownload.csv")
OUT_PATH = Path("reports/format_vocabulary.json")

# Declared formats to probe relevance against. "Whole mAb" is by far the most
# likely real submission; scFv and Nanobody probe the other end.
PROBE_FORMATS = ("Whole mAb", "scFv", "Nanobody", "Bispecific Whole mAb")


def main() -> None:
    rows = list(csv.DictReader(SOURCE.open(encoding="utf-8-sig")))
    values = Counter(r["Format"].strip() for r in rows)
    total_rows = len(rows)
    print(f"source: {total_rows} rows, {len(values)} distinct Format values")
    print(f"vocabulary version {FORMAT_VOCABULARY_VERSION}\n")

    axes = {
        "specificity": (lambda f: f.specificity, Specificity.UNKNOWN),
        "architecture": (lambda f: f.architecture, Architecture.UNKNOWN),
        "modification": (lambda f: f.modification, Modification.UNKNOWN),
    }

    per_axis = {}
    print("axis resolution (rows / distinct values):")
    for axis, (getter, unknown) in axes.items():
        unresolved_rows = 0
        unresolved_values = 0
        distribution = Counter()
        for raw, count in values.items():
            value = getter(parse_format(raw))
            distribution[value.value] += count
            if value is unknown:
                unresolved_rows += count
                unresolved_values += 1
        resolved_rows = total_rows - unresolved_rows
        per_axis[axis] = {
            "rows_resolved": resolved_rows,
            "rows_total": total_rows,
            "row_coverage": resolved_rows / total_rows,
            "distinct_values_unresolved": unresolved_values,
            "distinct_values_total": len(values),
            "distribution": dict(distribution.most_common()),
        }
        print(f"  {axis:14} {resolved_rows:5}/{total_rows} rows ({resolved_rows/total_rows:6.1%})   "
              f"{len(values) - unresolved_values:3}/{len(values)} distinct values")

    # Species has no UNKNOWN -- UNSTATED is a real answer, not a parse failure.
    species = Counter(parse_format(raw).species.value for raw, c in values.items() for _ in range(c))
    print(f"  {'species':14} UNSTATED is a real value, not a failure: {dict(species.most_common(4))}")

    print("\nrelevance status by declared format:")
    probes = {}
    for declared_raw in PROBE_FORMATS:
        declared = parse_format(declared_raw)
        statuses = Counter()
        relevances = []
        for raw, count in values.items():
            result = compute_format_relevance(declared, parse_format(raw))
            statuses[result.status.value] += count
            if result.relevance is not None:
                relevances.extend([result.relevance] * count)
        available = statuses.get(EvidenceStatus.AVAILABLE.value, 0)
        probes[declared_raw] = {
            "statuses": dict(statuses),
            "rows_scored": available,
            "row_fraction_scored": available / total_rows,
            "mean_relevance": (sum(relevances) / len(relevances)) if relevances else None,
        }
        mean = probes[declared_raw]["mean_relevance"]
        print(f"  declared {declared_raw:22} scored {available:5}/{total_rows} "
              f"({available/total_rows:6.1%})  mean relevance "
              f"{mean:.3f}" if mean is not None else "  n/a")

    unresolved_examples = sorted(
        (
            (count, raw) for raw, count in values.items()
            if parse_format(raw).architecture is Architecture.UNKNOWN
        ),
        reverse=True,
    )
    print(f"\narchitecture unresolved: {len(unresolved_examples)} distinct values, "
          f"{sum(c for c, _ in unresolved_examples)} rows")
    for count, raw in unresolved_examples[:8]:
        print(f"    {count:3}  {raw[:72]}")

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_file": str(SOURCE),
        "vocabulary_version": FORMAT_VOCABULARY_VERSION,
        "axis_weights": AXIS_WEIGHTS,
        "total_rows": total_rows,
        "distinct_values": len(values),
        "per_axis": per_axis,
        "species_distribution": dict(species.most_common()),
        "declared_format_probes": probes,
        "architecture_unresolved": [
            {"rows": c, "value": v} for c, v in unresolved_examples
        ],
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=1))
    print(f"\nwrote {OUT_PATH}")


if __name__ == "__main__":
    main()
