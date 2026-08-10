"""Gate 2: extend the SAbDab2 corpus snapshot with Thera-SAbDab, save a combined snapshot.

Requires data/cache/sabdab2_corpus_snapshot.pkl to already exist (built by
scripts/build_v2_corpus_cache.py). Loads it, ingests Thera-SAbDab into the
same Corpus object (so cross-source deduplication applies, per SPEC S15),
and writes the combined result to data/cache/combined_corpus_snapshot.pkl.

Thera-SAbDab's numbering cost is small relative to SAbDab2 (1,133 rows vs.
15,641), so this step runs in minutes, not hours.
"""

import time
from pathlib import Path

from para_scope.corpus import build_thera_sabdab_snapshot
from para_scope.models import ComponentType
from para_scope.persistence import load_corpus_snapshot, save_corpus_snapshot

SABDAB2_SNAPSHOT_PATH = Path("data/cache/sabdab2_corpus_snapshot.pkl")
THERA_CSV_PATH = Path("data/snapshots/thera_sabdab/TheraSAbDab_SeqStruc_OnlineDownload.csv")
COMBINED_SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")


def main() -> None:
    t0 = time.time()
    corpus, stats, sabdab2_manifest = load_corpus_snapshot(SABDAB2_SNAPSHOT_PATH)
    t1 = time.time()
    print(f"Loaded existing SAbDab2 snapshot in {t1 - t0:.2f}s ({len(corpus)} components)")

    thera_snapshot = build_thera_sabdab_snapshot(THERA_CSV_PATH)
    corpus.ingest_thera_sabdab(THERA_CSV_PATH, stats=stats)
    t2 = time.time()
    print(f"Thera-SAbDab ingestion: {t2 - t1:.2f}s")

    component_type_counts: dict[str, int] = {}
    for record in corpus._components.values():  # noqa: SLF001 -- one-off diagnostic script
        component_type_counts[record.component_type.value] = component_type_counts.get(record.component_type.value, 0) + 1

    from para_scope.corpus import build_sabdab2_snapshot

    sabdab2_snapshot = build_sabdab2_snapshot(Path("data/snapshots/sabdab2/splits_final/ab_split.csv"))
    manifest = save_corpus_snapshot(corpus, stats, [sabdab2_snapshot, thera_snapshot], COMBINED_SNAPSHOT_PATH)
    t3 = time.time()
    print(f"Combined snapshot saved in {t3 - t2:.2f}s -> {COMBINED_SNAPSHOT_PATH}")

    print(f"\nCombined stats: rows_seen={stats.rows_seen}, heavy={stats.heavy_chains_parsed}, "
          f"light={stats.light_chains_parsed}, unique_components={stats.unique_components}, "
          f"observed_pairs={stats.observed_pairs}, unique_parent_antibodies={stats.unique_parent_antibodies}")
    print(f"Component type counts: {component_type_counts}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
