"""One-off: build the real numbering + corpus snapshot cache (Version 2, Gate 0.5).

Pays the ~113-minute ANARCII cost once against the real SAbDab2 corpus. Every
gate from here on (Gate 1 search validation, Gate 3 benchmarks, Gate 7
therapeutic precedent, etc.) loads the resulting snapshot in well under a
second instead of re-ingesting from scratch.

Not part of the pytest suite -- run manually. Writes into data/cache/, which
is derived and gitignored (regenerable from the checksummed source snapshot).
"""

import time
from pathlib import Path

from para_scope.persistence import ingest_sabdab2_cached, save_corpus_snapshot

CSV_PATH = Path("data/snapshots/sabdab2/splits_final/ab_split.csv")
CACHE_DIR = Path("data/cache")
SNAPSHOT_PATH = Path("data/cache/sabdab2_corpus_snapshot.pkl")


def main() -> None:
    t0 = time.time()
    corpus, stats, snapshot = ingest_sabdab2_cached(CSV_PATH, CACHE_DIR)
    t1 = time.time()

    manifest = save_corpus_snapshot(corpus, stats, [snapshot], SNAPSHOT_PATH)
    t2 = time.time()

    print(f"Ingest (numbering-cache build, expected slow on first run): {t1 - t0:.2f}s")
    print(f"Snapshot save: {t2 - t1:.2f}s")
    print(f"Rows seen: {stats.rows_seen}")
    print(f"Heavy chains parsed: {stats.heavy_chains_parsed}")
    print(f"Light chains parsed: {stats.light_chains_parsed}")
    print(f"Unique components: {stats.unique_components}")
    print(f"Observed pairs: {stats.observed_pairs}")
    print(f"Unique parent antibodies: {stats.unique_parent_antibodies}")
    print(f"Snapshot written to: {SNAPSHOT_PATH}")
    print(f"Manifest created_at: {manifest.created_at}")


if __name__ == "__main__":
    main()
