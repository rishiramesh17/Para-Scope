"""
Gate 3 benchmark: does the k-mer shortlist lose true hits, and how much time does it save?

Brute-force search (shortlist_limit=None) is the ground truth. For each query we
compare the shortlisted top-k against the brute-force top-k by record_id, so a
"miss" means a record the exhaustive search would have reported and the
shortlisted search did not. Recall is measured, never assumed -- a shortlist is
only acceptable if the cost of adopting it is known.

Also sweeps k, because DEFAULT_K must be chosen against evidence rather than
intuition about alphabet size.
"""

import json
import statistics
import time
from pathlib import Path

from para_scope.indexes.shortlist import KmerShortlistIndex
from para_scope.models import ComponentType
from para_scope.persistence import load_corpus_snapshot

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
OUT_PATH = Path("reports/two_stage_search_benchmark.json")

TOP_K = 10
SHORTLIST_LIMITS = [50, 100, 200, 500]
K_VALUES = [4, 5, 6]
# Component types with enough indexed sequences for the comparison to be meaningful.
TYPES = [ComponentType.VH, ComponentType.VL_KAPPA, ComponentType.VHH, ComponentType.FR_H3]
QUERIES_PER_TYPE = 12


def _pick_queries(corpus, component_type, n):
    """Deterministic spread of real indexed sequences of this type."""
    seqs = corpus.approximate._sequences.get(component_type, {})  # noqa: SLF001 -- benchmark
    ids = sorted(seqs.keys())
    if not ids:
        return []
    step = max(1, len(ids) // n)
    return [(rid, seqs[rid]) for rid in ids[::step][:n]]


def main() -> None:
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    print(f"Loaded corpus: {len(corpus)} components, "
          f"{len(corpus.approximate)} approximate-eligible\n")

    results = {
        "corpus": {
            "components": len(corpus),
            "approximate_eligible": len(corpus.approximate),
            "corpus_schema_version": manifest.corpus_schema_version,
        },
        "top_k": TOP_K,
        "by_type": {},
        "k_sweep": {},
    }

    # ---- Recall + timing at the default k, across shortlist limits ----
    for component_type in TYPES:
        queries = _pick_queries(corpus, component_type, QUERIES_PER_TYPE)
        n_subjects = len(corpus.approximate._sequences.get(component_type, {}))  # noqa: SLF001
        if not queries:
            continue
        print(f"=== {component_type.value} ({n_subjects} subjects, {len(queries)} queries) ===")

        brute_truth = {}
        brute_times = []
        for rid, seq in queries:
            t0 = time.perf_counter()
            hits = corpus.approximate.search(component_type, seq, top_k=TOP_K, shortlist_limit=None)
            brute_times.append(time.perf_counter() - t0)
            brute_truth[rid] = [h.record_id for h in hits]

        type_entry = {
            "n_subjects": n_subjects,
            "n_queries": len(queries),
            "brute_force_p50_s": statistics.median(brute_times),
            "limits": {},
        }
        print(f"  brute force: p50 {statistics.median(brute_times):.3f}s")

        for limit in SHORTLIST_LIMITS:
            times = []
            recalls = []
            top1_agree = 0
            for rid, seq in queries:
                t0 = time.perf_counter()
                hits = corpus.approximate.search(
                    component_type, seq, top_k=TOP_K, shortlist_limit=limit
                )
                times.append(time.perf_counter() - t0)
                got = [h.record_id for h in hits]
                truth = brute_truth[rid]
                if truth:
                    recalls.append(len(set(got) & set(truth)) / len(truth))
                    if got and got[0] == truth[0]:
                        top1_agree += 1
                else:
                    recalls.append(1.0)
                    if not got:
                        top1_agree += 1

            p50 = statistics.median(times)
            mean_recall = statistics.mean(recalls)
            speedup = statistics.median(brute_times) / p50 if p50 else float("inf")
            type_entry["limits"][str(limit)] = {
                "p50_s": p50,
                "speedup_vs_brute": speedup,
                "mean_recall_at_10": mean_recall,
                "min_recall_at_10": min(recalls),
                "top1_agreement": top1_agree / len(queries),
            }
            print(f"  limit={limit:<4} p50 {p50:.4f}s  speedup {speedup:6.1f}x  "
                  f"recall@10 mean {mean_recall:.4f} min {min(recalls):.4f}  "
                  f"top1 agree {top1_agree}/{len(queries)}")
        results["by_type"][component_type.value] = type_entry
        print()

    # ---- k sweep on VH, at the default limit ----
    print("=== k sweep (VH, limit=200) ===")
    vh_seqs = corpus.approximate._sequences.get(ComponentType.VH, {})  # noqa: SLF001
    vh_queries = _pick_queries(corpus, ComponentType.VH, QUERIES_PER_TYPE)
    vh_truth = {}
    for rid, seq in vh_queries:
        hits = corpus.approximate.search(ComponentType.VH, seq, top_k=TOP_K, shortlist_limit=None)
        vh_truth[rid] = [h.record_id for h in hits]

    for k in K_VALUES:
        # Build a standalone shortlist index at this k over the same sequences.
        idx = KmerShortlistIndex(k=k)
        for rid, seq in vh_seqs.items():
            idx.add(ComponentType.VH, rid, seq)

        recalls = []
        shortlist_times = []
        for rid, seq in vh_queries:
            t0 = time.perf_counter()
            cands = idx.shortlist(ComponentType.VH, seq, limit=200)
            shortlist_times.append(time.perf_counter() - t0)
            cand_ids = {c.record_id for c in cands}
            truth = vh_truth[rid]
            # Recall of stage 1 alone: does the shortlist even contain the true hits?
            recalls.append(len(cand_ids & set(truth)) / len(truth) if truth else 1.0)

        postings = sum(len(v) for v in idx._postings.get(ComponentType.VH, {}).values())  # noqa: SLF001
        results["k_sweep"][str(k)] = {
            "stage1_mean_recall_at_10": statistics.mean(recalls),
            "stage1_min_recall_at_10": min(recalls),
            "stage1_p50_s": statistics.median(shortlist_times),
            "posting_entries": postings,
        }
        print(f"  k={k}  stage1 recall@10 mean {statistics.mean(recalls):.4f} "
              f"min {min(recalls):.4f}  shortlist p50 {statistics.median(shortlist_times)*1000:.2f}ms  "
              f"postings {postings:,}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {OUT_PATH}")


if __name__ == "__main__":
    main()
