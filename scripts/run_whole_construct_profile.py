"""
Gate 11: run a real candidate end to end through every gate and emit a §37 profile.

This is the integration that nothing before it exercised: Gates 1-3 produce
matches, Gate 7 scores each component, Gate 9 calibrates, Gate 10 grades, Gate 10.5
supplies format relevance, and Gate 11 assembles. A defect in the seams between
them cannot show up in any single gate's own tests.

Runs three candidates deliberately:

  self        a therapeutic's own VH/VL. Every term should be 1.0 and the index
              100 -- a trivial case, but one where any number other than 100 means
              a seam is broken.
  mutated     the same candidate with substitutions. The index must drop, and the
              per-component scores must show *where* it dropped.
  novel       a randomly generated pair. The index must be withheld, not low.

Writes reports/whole_construct_profile.json.
"""

import json
import random
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from para_scope.calibration import load_calibration
from para_scope.construct_format import parse_format
from para_scope.evidence_strength import (
    EvidenceGrade,
    assess_evidence,
    from_pairing_precedent,
    from_therapeutic_precedent,
)
from para_scope.models import ComponentType, EvidenceStatus
from para_scope.numbering import extract_region_sequence, number_sequence
from para_scope.pairing_precedent import compute_pairing_precedent
from para_scope.persistence import load_corpus_snapshot
from para_scope.profile import (
    ComponentProfileEntry,
    WholeConstructProfile,
    assemble_index,
    build_component_evidence,
    slot_for_component,
)
from para_scope.therapeutic_precedent import compute_therapeutic_precedent
from para_scope.verification.adversarial import point_mutant, random_sequence

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
CALIBRATION_PATH = Path("data/cache/calibration_v1.json")
OUT_PATH = Path("reports/whole_construct_profile.json")

REGION_TO_TYPE = {
    ("H", "FR1"): "FR-H1", ("H", "CDR1"): "CDR-H1", ("H", "FR2"): "FR-H2",
    ("H", "CDR2"): "CDR-H2", ("H", "FR3"): "FR-H3", ("H", "CDR3"): "CDR-H3",
    ("H", "FR4"): "FR-H4",
    ("L", "FR1"): "FR-L1", ("L", "CDR1"): "CDR-L1", ("L", "FR2"): "FR-L2",
    ("L", "CDR2"): "CDR-L2", ("L", "FR3"): "FR-L3", ("L", "CDR3"): "CDR-L3",
    ("L", "FR4"): "FR-L4",
}


def _entry(corpus, calibration, declared, component_type, sequence):
    """Score one component through Gates 7/9/10.5 and grade it with Gate 10."""
    precedent = compute_therapeutic_precedent(
        corpus, component_type, sequence, top_k=10, shortlist_limit=None,
        calibration=calibration, declared_format=declared,
    )
    strength = assess_evidence([from_therapeutic_precedent(precedent)])
    return ComponentProfileEntry(
        component_type=component_type.value,
        slot=(s.value if (s := slot_for_component(component_type.value)) else None),
        status=precedent.status.value,
        score=precedent.score,
        score_reason=precedent.score_reason,
        grade=strength.grade.value,
        best_identity=precedent.best_global_identity,
        therapeutics_matched=precedent.therapeutics_matched,
        blocked_terms=dict(precedent.blocked_terms),
        limitations=list(precedent.limitations),
    )


def profile_candidate(corpus, calibration, heavy, light, declared_raw="Whole mAb"):
    declared = parse_format(declared_raw)
    entries = []

    for sequence, chain in ((heavy, "H"), (light, "L")):
        if not sequence:
            continue
        numbered = number_sequence(sequence)
        if numbered is None:
            continue
        domain_type = (
            ComponentType.VH if chain == "H"
            else (ComponentType.VL_KAPPA if numbered.chain_type == "K"
                  else ComponentType.VL_LAMBDA)
        )
        entries.append(_entry(corpus, calibration, declared, domain_type, sequence))
        for region, type_name in (
            (r, t) for (c, r), t in REGION_TO_TYPE.items() if c == chain
        ):
            region_sequence = extract_region_sequence(numbered, region)
            if not region_sequence:
                continue
            entries.append(
                _entry(corpus, calibration, declared,
                       ComponentType(type_name), region_sequence)
            )

    pairing = compute_pairing_precedent(
        corpus, heavy, light, top_k=10, shortlist_limit=None,
        declared_format=declared,
    )
    pairing_dimension = from_pairing_precedent(pairing)
    pairing_strength = assess_evidence([pairing_dimension])

    # Pairing contributes its best joint coverage as a score; structure-linked
    # precedent contributes the share of matched pairs backed by a structure.
    pairing_score = (
        100.0 * pairing.best_joint_coverage
        if pairing.status is EvidenceStatus.AVAILABLE
        and pairing.best_joint_coverage is not None
        else None
    )
    structure_score = (
        100.0 * pairing_dimension.structure_linked_records / len(pairing.matches)
        if pairing.matches else None
    )

    evidence = build_component_evidence(
        entries,
        pairing_score=pairing_score,
        pairing_grade=pairing_strength.grade,
        structure_score=structure_score,
        structure_grade=pairing_strength.grade,
    )

    strength = assess_evidence(
        [from_therapeutic_precedent(
            compute_therapeutic_precedent(
                corpus, ComponentType.VH, heavy, top_k=10, shortlist_limit=None,
                calibration=calibration, declared_format=declared,
            )
        ), pairing_dimension],
        expected_components={"VH": heavy, "VL": light},
    )

    index = assemble_index(
        evidence, declared, strength=strength,
        observed_pairing_status=pairing.status,
        assembly_validated=True,
    )
    return WholeConstructProfile(
        components=entries, index=index, strength=strength,
        declared_format=declared_raw,
    )


def main() -> None:
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    calibration = load_calibration(CALIBRATION_PATH)
    print(f"corpus {len(corpus)} components (schema {manifest.corpus_schema_version})\n")

    # A whole-mAb therapeutic with both chains, chosen deterministically.
    pick = None
    for name in sorted(corpus.therapeutics):
        therapeutic = corpus.therapeutics[name]
        if therapeutic.format_raw.strip() != "Whole mAb" or not therapeutic.arms:
            continue
        arm = therapeutic.arms[0]
        if arm.heavy_record_id and arm.light_record_id:
            heavy = corpus.get(arm.heavy_record_id)
            light = corpus.get(arm.light_record_id)
            if heavy and light:
                pick = (name, heavy.sequence, light.sequence)
                break
    if pick is None:
        raise SystemExit("no paired whole-mAb therapeutic found")
    name, heavy, light = pick
    print(f"candidate source: {name} (VH {len(heavy)} aa, VL {len(light)} aa)\n")

    cases = {
        "self": (heavy, light),
        "mutated_10": (point_mutant("m10h", heavy, 10), point_mutant("m10l", light, 10)),
        "mutated_30": (point_mutant("m30h", heavy, 30), point_mutant("m30l", light, 30)),
        "novel": (random_sequence("novel_h", len(heavy)), random_sequence("novel_l", len(light))),
    }

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "corpus_schema_version": manifest.corpus_schema_version,
        "candidate_source": name,
        "cases": {},
    }

    for label, (h, l) in cases.items():
        profile = profile_candidate(corpus, calibration, h, l)
        index = profile.index
        value = f"{index.index:.1f}" if index.index is not None else "withheld"
        print(f"=== {label} ===")
        print(f"  index            {value}")
        print(f"  displayed_name   {index.displayed_name}")
        print(f"  scope            {index.coverage_scope.value}")
        print(f"  coverage         {index.evidence_coverage:.0%}")
        print(f"  bottleneck       {index.bottleneck_slot} = {index.bottleneck_score}")
        print(f"  confidence       {index.evidence_confidence}")
        print(f"  pairing          {index.observed_pairing_status}")
        print(f"  present          {len(index.components_present)} slots")
        print(f"  missing          {sorted(index.components_missing)}")
        if index.gating_failures:
            print(f"  gating           {index.gating_failures[0][:96]}")
        scored = [e for e in profile.components if e.score is not None]
        print(f"  components       {len(scored)}/{len(profile.components)} scored")
        for entry in sorted(scored, key=lambda e: e.score)[:3]:
            print(f"      weakest: {entry.component_type:10} {entry.score:6.1f} (grade {entry.grade})")
        print()
        payload["cases"][label] = asdict(profile)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(payload, indent=1, default=str))
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
