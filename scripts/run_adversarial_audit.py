"""
Run the adversarial self-audit against the real corpus and report measured results.

Writes reports/adversarial_audit.json. Parents are drawn deterministically from
the corpus so the suite is reproducible, and from several component classes so a
class-specific failure is visible rather than averaged away.
"""

import json
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from para_scope.calibration import load_calibration
from para_scope.models import ComponentType
from para_scope.persistence import load_corpus_snapshot
from para_scope.verification.adversarial import build_case_suite
from para_scope.verification.harness import run_suite

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
CALIBRATION_PATH = Path("data/cache/calibration_v1.json")
OUT_PATH = Path("reports/adversarial_audit.json")

PARENT_CLASSES = (ComponentType.VH, ComponentType.VL_KAPPA, ComponentType.VHH)
PARENTS_PER_CLASS = 2


def main() -> None:
    t0 = time.time()
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    calibration = load_calibration(CALIBRATION_PATH) if CALIBRATION_PATH.exists() else None
    print(f"corpus: {len(corpus)} components (schema {manifest.corpus_schema_version})")
    print(f"calibration: {'loaded' if calibration else 'ABSENT'}")

    # Deterministic parents: components that actually carry therapeutic links, so
    # the precedent engine has something to find for the unmodified parent. A
    # parent with no precedent makes "worse than parent" vacuous.
    linked = set(corpus._component_to_therapeutics)  # noqa: SLF001
    parents = []
    for component_type in PARENT_CLASSES:
        candidates = sorted(
            (
                r for r in corpus._components.values()  # noqa: SLF001
                if r.component_type is component_type and r.record_id in linked
            ),
            key=lambda r: r.record_id,
        )
        for record in candidates[:PARENTS_PER_CLASS]:
            parents.append((record.record_id, record.sequence, component_type))
    print(f"parents: {len(parents)} across {len(PARENT_CLASSES)} classes")

    cases = build_case_suite(parents)
    print(f"cases: {len(cases)}\n")

    report = run_suite(corpus, cases, calibration=calibration)

    print(f"PASSED {report.passed}/{report.total}\n")
    print("by family:")
    for family, (passed, total) in report.by_family().items():
        flag = "" if passed == total else "   <-- FAILURES"
        print(f"  {family:22} {passed:3}/{total:3}{flag}")

    print("\ngrounding boundary:")
    print(f"  answers rendered            {report.answers_rendered}")
    print(f"  explicit refusals           {report.answers_refused}")
    print(f"  prohibited claims emitted   {report.prohibited_claim_count}")
    print(f"  citation coverage           {report.citation_coverage:.3f}")
    print(f"  uncited claims              {report.uncited_claim_count}")
    print(f"  crashes                     {len(report.crashes)}")

    if report.failures:
        print(f"\n{len(report.failures)} FAILURE(S):")
        for outcome in report.failures:
            print(f"  [{outcome.family}] {outcome.case_id}: {outcome.detail}")
    if report.crashes:
        print("\nCRASHES:")
        for crash in report.crashes:
            print(f"  {crash}")
    if report.prohibited_claim_examples:
        print("\nPROHIBITED-CLAIM EVENTS:")
        for example in report.prohibited_claim_examples[:10]:
            print(f"  {example}")

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "corpus_schema_version": manifest.corpus_schema_version,
        "calibration_loaded": calibration is not None,
        "adversarial_version": report.adversarial_version,
        "grounding_version": report.grounding_version,
        "harness_version": report.harness_version,
        "totals": {
            "cases": report.total,
            "passed": report.passed,
            "failed": len(report.failures),
            "crashes": len(report.crashes),
            "prohibited_claims": report.prohibited_claim_count,
            "citation_coverage": report.citation_coverage,
            "uncited_claims": report.uncited_claim_count,
            "answers_rendered": report.answers_rendered,
            "explicit_refusals": report.answers_refused,
        },
        "by_family": {k: {"passed": v[0], "total": v[1]} for k, v in report.by_family().items()},
        "outcomes": [asdict(o) for o in report.outcomes],
        "crash_detail": report.crashes,
        "prohibited_claim_detail": report.prohibited_claim_examples,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=1, default=str))
    print(f"\nwrote {OUT_PATH}  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
