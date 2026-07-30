"""One-off validation/benchmark run: full SAbDab2 ingestion (Phase 1a acceptance check).

Not part of the pytest suite -- this exercises the real ~15,641-row downloaded
corpus (data/snapshots/, gitignored) and is meant to be run manually to produce
the numbers reported in docs/PHASE_1_RESULTS.md.
"""

import json
import time
from pathlib import Path

from antigens.corpus import Corpus, build_sabdab2_snapshot
from antigens.models import ComponentType

CSV_PATH = Path("data/snapshots/sabdab2/splits_final/ab_split.csv")
OUTPUT_PATH = Path("data/snapshots/sabdab2/ingest_report.json")


def main() -> None:
    t0 = time.time()
    snapshot = build_sabdab2_snapshot(CSV_PATH)
    t1 = time.time()

    corpus = Corpus()
    stats = corpus.ingest_sabdab2(CSV_PATH)
    t2 = time.time()

    component_type_counts = {}
    for record in corpus._components.values():  # noqa: SLF001 -- one-off diagnostic script
        component_type_counts[record.component_type.value] = component_type_counts.get(record.component_type.value, 0) + 1

    # Benchmark: query latency for exact + approximate + positional search, using a
    # real ingested VH sequence as the query.
    sample_record = next(
        r for r in corpus._components.values() if r.component_type == ComponentType.VH  # noqa: SLF001
    )
    query_seq = sample_record.sequence

    t3 = time.time()
    exact_hits = corpus.exact.lookup(sample_record.sequence_checksum)
    t4 = time.time()
    approx_hits = corpus.approximate.search(ComponentType.VH, query_seq, top_k=10)
    t5 = time.time()
    positions = {pos for (pos, _code), _res in sample_record.numbered.imgt_numbering}
    positional_hits = corpus.positional.lookup(ComponentType.VH, positions)
    t6 = time.time()

    report = {
        "snapshot": {
            "source_name": snapshot.source_name,
            "raw_checksum": snapshot.raw_checksum,
            "record_count": snapshot.record_count,
            "parser_version": snapshot.parser_version,
        },
        "timing_seconds": {
            "snapshot_build": round(t1 - t0, 2),
            "ingest": round(t2 - t1, 2),
        },
        "stats": {
            "rows_seen": stats.rows_seen,
            "rows_skipped_unsupported": stats.rows_skipped_unsupported,
            "heavy_chains_parsed": stats.heavy_chains_parsed,
            "light_chains_parsed": stats.light_chains_parsed,
            "numbering_failures": stats.numbering_failures,
            "unique_components": stats.unique_components,
            "observed_pairs": stats.observed_pairs,
            "unique_parent_antibodies": stats.unique_parent_antibodies,
            "warning_count": len(stats.warnings),
        },
        "index_sizes": {
            "exact": len(corpus.exact),
            "approximate": len(corpus.approximate),
            "positional": len(corpus.positional),
            "pairs": len(corpus.pairs),
        },
        "component_type_counts": component_type_counts,
        "query_benchmark_seconds": {
            "exact_lookup": round(t4 - t3, 5),
            "approximate_search_top10": round(t5 - t4, 3),
            "positional_lookup": round(t6 - t5, 5),
        },
        "sample_warnings": stats.warnings[:20],
    }

    warning_categories = {"unsupported_component": 0, "invalid_input": 0, "anarcii_numbering_failed": 0, "chain_type_mismatch": 0, "other": 0}
    for w in stats.warnings:
        if "unsupported_component" in w:
            warning_categories["unsupported_component"] += 1
        elif "invalid input" in w:
            warning_categories["invalid_input"] += 1
        elif "ANARCII numbering failed" in w:
            warning_categories["anarcii_numbering_failed"] += 1
        elif "chain-type mismatch" in w:
            warning_categories["chain_type_mismatch"] += 1
        else:
            warning_categories["other"] += 1
    report["warning_categories"] = warning_categories

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(report, indent=2))
    (OUTPUT_PATH.parent / "ingest_warnings_full.json").write_text(json.dumps(stats.warnings, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
