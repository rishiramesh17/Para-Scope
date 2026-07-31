"""
Gate 7 measurements: therapeutic-precedent coverage over the real corpus.

Reports what the therapeutic entity store actually contains, confirms the
structure-availability fix, and runs the SPEC S26 dimension end to end on a known
therapeutic and a known non-therapeutic input.

Run: .venv/bin/python scripts/measure_therapeutic_precedent.py
"""

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "src")

from antigens.models import ComponentType, EvidenceStatus
from antigens.persistence import load_corpus_snapshot
from antigens.therapeutic_precedent import compute_therapeutic_precedent

SNAPSHOT = Path("data/cache/combined_corpus_snapshot.pkl")

# Trastuzumab's variable domains -- a therapeutic that is certainly in the corpus.
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
# IMGT CDR3 spans positions 105-117, so it opens with SR -- not the WGGDGF... a
# Kabat-shaped intuition suggests. Region boundaries are numbering-scheme specific.
TRASTUZUMAB_CDR_H3 = "SRWGGDGFYAMDY"


def main() -> None:
    corpus, stats, manifest = load_corpus_snapshot(SNAPSHOT)
    therapeutics = corpus.therapeutics

    print("=" * 78)
    print("THERAPEUTIC ENTITY STORE")
    print("=" * 78)
    print(f"  corpus_schema_version        {corpus.schema_version}")
    print(f"  canonical therapeutics       {len(therapeutics):,}")
    print(f"  arms                         {sum(len(t.arms) for t in therapeutics.values()):,}")
    print(f"  alias keys resolvable        {len(corpus._therapeutic_aliases):,}")
    print(f"  ambiguous aliases dropped    {len(corpus._ambiguous_therapeutic_aliases):,}")
    print(f"  alias collisions (stats)     {stats.therapeutic_alias_collisions:,}")
    print(f"  components linked            {len(corpus._component_to_therapeutics):,}")

    heavy_only = [
        arm for t in therapeutics.values() for arm in t.arms
        if arm.heavy_record_id and not arm.light_record_id
    ]
    print(f"  heavy-only arms              {len(heavy_only):,} "
          f"(metadata now reachable; previously unreachable)")

    print("\n" + "=" * 78)
    print("STRUCTURE AVAILABILITY (the defect this gate fixed)")
    print("=" * 78)
    counts = Counter(t.has_structure_record for t in therapeutics.values())
    for value in (True, False, None):
        label = {True: "structure-backed", False: "no structure", None: "unknown"}[value]
        print(f"  {label:<22} {counts.get(value, 0):>6,}")
    print(f"  Before the fix every therapeutic reported structure-backed.")

    print("\n" + "=" * 78)
    print("METADATA COMPLETENESS (SPEC S26 term)")
    print("=" * 78)
    completeness = Counter(round(t.metadata_completeness, 3) for t in therapeutics.values())
    for value, count in sorted(completeness.items(), reverse=True):
        print(f"  {value:.3f}   {count:>6,}")

    print("\n" + "=" * 78)
    print("CLINICAL STAGE / DEVELOPMENT STATUS (metadata only, never a label)")
    print("=" * 78)
    for field_name in ("highest_clinical_trial", "development_status"):
        print(f"  {field_name}")
        for value, count in Counter(
            getattr(t, field_name) for t in therapeutics.values()
        ).most_common(6):
            print(f"    {str(value):<26} {count:>5,}")

    print("\n" + "=" * 78)
    print("SPEC S26 DIMENSION — trastuzumab VH")
    print("=" * 78)
    result = compute_therapeutic_precedent(corpus, ComponentType.VH, TRASTUZUMAB_VH, top_k=5)
    print(f"  status                   {result.status.value}")
    print(f"  therapeutics_matched     {result.therapeutics_matched}   "
          f"(raw observations: {result.raw_match_count})")
    print(f"  exact matches            {result.exact_match_count}")
    print(f"  best_global_identity     {result.best_global_identity:.4f}")
    print(f"  score                    {result.score}  [{result.score_status.value}]")
    print(f"  computable_terms         {result.computable_terms}")
    print(f"  blocked_terms            {sorted(result.blocked_terms)}")
    print()
    header = f"  {'therapeutic':<22}{'ident':>7}{'qcov':>7}{'pair':>6}{'struct':>8}  {'stage':<14}target"
    print(header)
    print("  " + "-" * 74)
    for match in result.matches[:5]:
        print(f"  {match.therapeutic_name[:21]:<22}{match.global_identity:>7.3f}"
              f"{match.query_coverage:>7.3f}{str(match.directly_paired):>6}"
              f"{str(match.has_structure_record):>8}  "
              f"{str(match.highest_clinical_trial):<14}{str(match.target)[:24]}")

    print("\n" + "=" * 78)
    print("ALIAS RESOLUTION (SPEC S15.5)")
    print("=" * 78)
    for probe in ("Trastuzumab", "trastuzumab", "Herceptin", "rhuMAb HER2", "not-a-drug"):
        found = corpus.get_therapeutic(probe)
        print(f"  {probe:<16} -> {found.therapeutic_name if found else '(unresolved)'}")

    print("\n" + "=" * 78)
    print("SHORT COMPONENT — exact matching runs, alignment does not")
    print("=" * 78)
    # IMGT CDR-H3 begins at position 105, so trastuzumab's is SRWGGDGFYAMDY.
    cdr = compute_therapeutic_precedent(corpus, ComponentType.CDR_H3, TRASTUZUMAB_CDR_H3)
    print(f"  CDR-H3 status                 {cdr.status.value}")
    print(f"  approximate_search_status     {cdr.approximate_search_status.value}")
    print(f"  exact matches                 {cdr.exact_match_count}")
    print(f"  therapeutics_matched          {cdr.therapeutics_matched}")
    for match in cdr.matches[:3]:
        print(f"    {match.therapeutic_name:<20} {str(match.target)[:44]}")
    for limitation in cdr.limitations:
        print(f"  limitation: {limitation}")

    print("\n" + "=" * 78)
    print("NEGATIVE CONTROL — a sequence no therapeutic contains")
    print("=" * 78)
    absent = compute_therapeutic_precedent(
        corpus, ComponentType.CDR_H3, "WWWWWWWWWWWWW"
    )
    print(f"  status                        {absent.status.value}")
    print(f"  reason                        {absent.reason}")
    print(f"  score                         {absent.score}  (never 0.0)")

    print("\n" + "=" * 78)
    print("DISCLAIMER CARRIED (SPEC S26, verbatim)")
    print("=" * 78)
    print(f"  {result.disclaimer}")


if __name__ == "__main__":
    main()
