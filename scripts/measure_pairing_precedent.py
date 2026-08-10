"""
Gate 8: measure observed pairing precedent against the real combined corpus.

Reports the two counts that must never be conflated (unique pair keys vs backing
source records), the corroboration distribution, and a worked example.
"""

from collections import Counter
from pathlib import Path

from para_scope.models import ComponentType, EvidenceStatus
from para_scope.pairing_precedent import compute_pairing_precedent
from para_scope.persistence import load_corpus_snapshot

SNAPSHOT = Path("data/cache/combined_corpus_snapshot.pkl")

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

BAR = "=" * 78


def main() -> None:
    corpus, stats, manifest = load_corpus_snapshot(SNAPSHOT)

    print(BAR)
    print("CORPUS PAIRING SHAPE")
    print(BAR)
    unique_keys = len(corpus.pairs)
    total_records = sum(len(v) for v in corpus.pairs._index.values())
    print(f"  unique (heavy, light) keys : {unique_keys:,}")
    print(f"  backing PairRecord entries : {total_records:,}")
    print(f"  IngestStats.observed_pairs : {stats.observed_pairs:,}")
    print(f"  records per unique pair    : {total_records / unique_keys:.2f}")

    depth = Counter(len(v) for v in corpus.pairs._index.values())
    print("\n  corroboration depth (how many records back one pair):")
    for n in sorted(depth)[:8]:
        print(f"    {n:>3} record(s): {depth[n]:>6,} pairs")
    deep = sum(c for n, c in depth.items() if n > 8)
    if deep:
        print(f"    >8 records: {deep:>6,} pairs")

    by_source = Counter()
    structure_state = Counter()
    for records in corpus.pairs._index.values():
        by_source[tuple(sorted({r.source_database for r in records}))] += 1
        states = [r.structure_available for r in records]
        if all(s is None for s in states):
            structure_state["unknown"] += 1
        elif any(states):
            structure_state["structure-backed"] += 1
        else:
            structure_state["stated absent"] += 1

    print("\n  unique pairs by source combination:")
    for sources, count in by_source.most_common():
        print(f"    {'+'.join(sources):<26} {count:>6,}")

    print("\n  unique pairs by structure support:")
    for state, count in structure_state.most_common():
        print(f"    {state:<26} {count:>6,}")

    print(f"\n{BAR}")
    print("WORKED EXAMPLE - trastuzumab VH + VL")
    print(BAR)
    result = compute_pairing_precedent(
        corpus, TRASTUZUMAB_VH, TRASTUZUMAB_VL,
        heavy_component_type=ComponentType.VH,
        light_component_type=ComponentType.VL_KAPPA,
        top_k=5,
        shortlist_limit=500,
    )
    print(f"  status                 {result.status.value}")
    print(f"  pairs_matched          {result.pairs_matched}")
    print(f"  raw_observation_count  {result.raw_observation_count}")
    print(f"  exact_pair_matches     {result.exact_pair_match_count}")
    print(f"  best_joint_coverage    {result.best_joint_coverage}")
    print("\n  computable dimensions (raw, unweighted):")
    for name, value in result.computable_dimensions.items():
        print(f"    {name:<38} {value}")
    print("\n  blocked dimensions:")
    for name in result.blocked_dimensions:
        print(f"    {name}")
    print("\n  top matches:")
    print(f"    {'heavy id':<22}{'light id':<22}{'H id':>7}{'L id':>7}{'obs':>5}  sources")
    for m in result.matches:
        print(f"    {m.heavy_record_id[:21]:<22}{m.light_record_id[:21]:<22}"
              f"{m.heavy_global_identity:>7.3f}{m.light_global_identity:>7.3f}"
              f"{m.observation_count:>5}  {','.join(m.source_databases)}")

    print(f"\n{BAR}")
    print("SINGLE-DOMAIN CANDIDATE (no light chain)")
    print(BAR)
    vhh = compute_pairing_precedent(
        corpus, TRASTUZUMAB_VH, None, heavy_component_type=ComponentType.VHH
    )
    print(f"  status  {vhh.status.value}   <- must be not_applicable, never not_found")

    print(f"\n{BAR}")
    print(result.disclaimer)
    print(BAR)


if __name__ == "__main__":
    main()
