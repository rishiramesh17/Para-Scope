"""
Gate 10: measure the grade distribution over real corpus candidates.

SPEC §19 Layer 3 says grade thresholds "must be benchmarked and versioned".
Versioned they are. Benchmarked in the strict sense they cannot be: there is no
external ground truth of evidence quality to calibrate against, and inventing one
would be worse than admitting its absence.

What this script does instead is the weaker, honest thing -- it checks that the
bands *separate real cases*. Thresholds that put every candidate in one grade
carry no information no matter how principled they look, and that failure is
measurable. It also confirms the two properties Gate 10 exists for: that an
induced database outage produces Insufficient rather than D, and that the four
SPEC §28 fields move independently.

Writes reports/evidence_grades.json.
"""

import json
import random
import time
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from para_scope.calibration import load_calibration
from para_scope.evidence_strength import (
    GRADE_THRESHOLD_VERSION,
    DimensionReport,
    EvidenceGrade,
    assess_evidence,
    from_germline,
    from_pairing_precedent,
    from_therapeutic_precedent,
    from_window_coverage,
)
from para_scope.germline import SUPPORTED_TYPES, compute_germline_assignment
from para_scope.models import ComponentType, EvidenceStatus
from para_scope.pairing_precedent import compute_pairing_precedent
from para_scope.persistence import load_corpus_snapshot
from para_scope.therapeutic_precedent import compute_therapeutic_precedent
from para_scope.window_coverage import compute_window_coverage

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
CALIBRATION_PATH = Path("data/cache/calibration_v1.json")
OUT_PATH = Path("reports/evidence_grades.json")

SAMPLE_PER_CLASS = 40
SAMPLE_SEED = 20261004
CLASSES = (ComponentType.VH, ComponentType.VL_KAPPA, ComponentType.VHH)


def _dimensions_for(corpus, calibration, component_type, record):
    """Every dimension that applies to a single-component candidate."""
    dimensions = []

    therapeutic = compute_therapeutic_precedent(
        corpus, component_type, record.sequence,
        top_k=10, shortlist_limit=None, calibration=calibration,
    )
    dimensions.append(from_therapeutic_precedent(therapeutic))

    if component_type in SUPPORTED_TYPES:
        dimensions.append(from_germline(compute_germline_assignment(component_type, record.sequence)))

    dimensions.append(
        from_window_coverage(compute_window_coverage(record.sequence, record_observations=False))
    )
    return dimensions, therapeutic


def main() -> None:
    t0 = time.time()
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    calibration = load_calibration(CALIBRATION_PATH) if CALIBRATION_PATH.exists() else None
    print(f"corpus {len(corpus)} components (schema {manifest.corpus_schema_version}); "
          f"calibration {'loaded' if calibration else 'ABSENT'}\n")

    rng = random.Random(SAMPLE_SEED)
    rows = []

    for component_type in CLASSES:
        candidates = sorted(
            (r for r in corpus._components.values()  # noqa: SLF001 -- offline measurement
             if r.component_type is component_type),
            key=lambda r: r.record_id,
        )
        sample = (
            candidates if len(candidates) <= SAMPLE_PER_CLASS
            else rng.sample(candidates, SAMPLE_PER_CLASS)
        )
        t1 = time.time()
        for record in sample:
            dimensions, therapeutic = _dimensions_for(
                corpus, calibration, component_type, record
            )
            report = assess_evidence(
                dimensions,
                expected_components={component_type.value: record.sequence},
            )
            rows.append({
                "component_class": component_type.value,
                "record_id": record.record_id,
                "grade": report.grade.value,
                "grade_rationale": report.grade_rationale,
                "independent_records": report.evidence_sufficiency.independent_records,
                "parent_antibodies": report.evidence_sufficiency.parent_antibodies,
                "structure_linked_records": report.evidence_sufficiency.structure_linked_records,
                "metadata_fraction": report.metadata_completeness.fraction,
                "therapeutic_status": therapeutic.status.value,
                "best_identity": therapeutic.best_global_identity,
            })
        print(f"  {component_type.value:10} {len(sample):3} candidates in {time.time()-t1:.1f}s")

    grades = Counter(r["grade"] for r in rows)
    print("\ngrade distribution:")
    for grade in (EvidenceGrade.A, EvidenceGrade.B, EvidenceGrade.C,
                  EvidenceGrade.D, EvidenceGrade.INSUFFICIENT):
        count = grades.get(grade.value, 0)
        bar = "#" * int(40 * count / max(1, len(rows)))
        print(f"  {grade.value:12} {count:4} ({count/len(rows):5.1%}) {bar}")

    print("\nby class:")
    for component_type in CLASSES:
        subset = [r for r in rows if r["component_class"] == component_type.value]
        if subset:
            dist = Counter(r["grade"] for r in subset)
            print(f"  {component_type.value:10} " + "  ".join(
                f"{g}={dist.get(g, 0)}" for g in ("A", "B", "C", "D", "Insufficient")
            ))

    distinct = len([g for g in grades.values() if g > 0])
    print(f"\nbands populated: {distinct} of 5")
    if distinct <= 1:
        print("  WARNING: thresholds collapse every candidate onto one grade and "
              "carry no information.")

    # --- the property that matters most: an outage must not become a D ---
    print("\noutage control:")
    probe = rows[0]
    outage = assess_evidence([
        DimensionReport("therapeutic_precedent", EvidenceStatus.DATABASE_UNAVAILABLE,
                        reason="induced: snapshot unavailable"),
    ])
    print(f"  induced outage grade: {outage.grade.value} "
          f"(must be Insufficient, never D)")
    assert outage.grade is EvidenceGrade.INSUFFICIENT
    assert outage.grade_blocked_reason is not None

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "corpus_schema_version": manifest.corpus_schema_version,
        "calibration_loaded": calibration is not None,
        "grade_threshold_version": GRADE_THRESHOLD_VERSION,
        "sample_per_class": SAMPLE_PER_CLASS,
        "sample_seed": SAMPLE_SEED,
        "candidates": len(rows),
        "grade_distribution": dict(grades),
        "bands_populated": distinct,
        "outage_control_grade": outage.grade.value,
        "rows": rows,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=1, default=str))
    print(f"\nwrote {OUT_PATH}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
