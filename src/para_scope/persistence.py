"""
Persistence layer for the ingested corpus and its numbering cache (Version 2, Gate 0.5).

ANARCII numbering dominates corpus ingestion cost: 511s of a 512s measured run
over the real 15,641-row SAbDab2 split (39,765 unique components), versus
sub-second parse+index. Every gate from here on needs to load an
already-ingested corpus repeatedly (integration tests, search benchmarks,
regression comparisons) -- re-paying that cost on every run is impractical.
This module adds two independent, derived caches. Neither is committed to git;
both are regenerable at any time from the checksummed source snapshot recorded
in SourceSnapshot.

1. Numbering cache: checksum -> NumberedComponent, keyed by (source file
   checksum, ANARCII version). Lets corpus.py's indexing logic change without
   forcing a renumber of sequences that haven't changed.
2. Corpus snapshot: the fully-ingested Corpus object (all four indexes plus
   components), keyed by the source file(s) it was built from. Lets tests and
   benchmarks load a ready corpus without re-running ingestion at all.
"""

import pickle
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from para_scope.corpus import CORPUS_SCHEMA_VERSION, Corpus, IngestStats
from para_scope.models import SourceSnapshot
from para_scope.numbering import NumberedComponent, get_anarcii_version

# 1.2.0 -- source_snapshots now persist all eight SourceSnapshot fields; earlier
#          formats dropped license_id, snapshot_id, retrieved_at, schema_version.
CACHE_FORMAT_VERSION = "1.2.0"


class StaleSnapshotError(RuntimeError):
    """A snapshot on disk was built under a Corpus/cache shape this code no longer matches."""


@dataclass
class CorpusManifest:
    """Human-readable summary stored alongside (not instead of) the binary snapshot."""
    cache_format_version: str
    created_at: str
    anarcii_version: str
    source_snapshots: list[dict]
    stats: dict
    corpus_schema_version: str = "unknown"


def _numbering_cache_path(cache_dir: Path, source_checksum: str) -> Path:
    # Filename embeds the raw checksum directly (not re-hashed) so a stale cache
    # is identifiable by inspection, not just by a load-time comparison.
    return cache_dir / f"numbering_{source_checksum[:16]}_{get_anarcii_version()}.pkl"


def load_numbering_cache(
    cache_dir: Path, source_checksum: str
) -> Optional[dict[str, Optional[NumberedComponent]]]:
    """Load a previously-saved numbering cache, or None if no cache exists for this source+tool version."""
    path = _numbering_cache_path(cache_dir, source_checksum)
    if not path.exists():
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def save_numbering_cache(
    cache_dir: Path, source_checksum: str, cache: dict[str, Optional[NumberedComponent]]
) -> Path:
    """Persist a numbering cache (checksum -> NumberedComponent) to disk."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = _numbering_cache_path(cache_dir, source_checksum)
    with open(path, "wb") as f:
        pickle.dump(cache, f, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def _stats_to_dict(stats: IngestStats) -> dict:
    return {
        "rows_seen": stats.rows_seen,
        "rows_skipped_unsupported": stats.rows_skipped_unsupported,
        "heavy_chains_parsed": stats.heavy_chains_parsed,
        "light_chains_parsed": stats.light_chains_parsed,
        "numbering_failures": stats.numbering_failures,
        "unique_components": stats.unique_components,
        "observed_pairs": stats.observed_pairs,
        "unique_parent_antibodies": stats.unique_parent_antibodies,
        "warning_count": len(stats.warnings),
    }


def _snapshot_to_dict(snapshot: SourceSnapshot) -> dict:
    """
    Persist every SourceSnapshot field, not a subset.

    This previously stored four fields and dropped license_id, snapshot_id,
    retrieved_at, and schema_version. The licensing omission was the serious one:
    SPEC §13/§14 require a source's licence to travel with the evidence derived
    from it, and one ingested source carries an `unspecified` licence that blocks
    redistribution -- so the single fact most needed to audit a stored corpus was
    the one being discarded at save time. A snapshot must be answerable about its
    own provenance without re-deriving it from source files that may have moved,
    changed, or gone away.
    """
    return {
        "source_name": snapshot.source_name,
        "snapshot_id": snapshot.snapshot_id,
        "retrieved_at": snapshot.retrieved_at,
        "license_id": snapshot.license_id,
        "raw_checksum": snapshot.raw_checksum,
        "record_count": snapshot.record_count,
        "parser_version": snapshot.parser_version,
        "schema_version": snapshot.schema_version,
    }


def save_corpus_snapshot(
    corpus: Corpus,
    stats: IngestStats,
    source_snapshots: list[SourceSnapshot],
    path: Path,
) -> CorpusManifest:
    """Serialize a fully-ingested Corpus (indexes + components) plus a readable manifest."""
    manifest = CorpusManifest(
        cache_format_version=CACHE_FORMAT_VERSION,
        created_at=datetime.now(timezone.utc).isoformat(),
        anarcii_version=get_anarcii_version(),
        source_snapshots=[_snapshot_to_dict(s) for s in source_snapshots],
        stats=_stats_to_dict(stats),
        corpus_schema_version=getattr(corpus, "schema_version", "unknown"),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({"manifest": manifest, "corpus": corpus, "stats": stats}, f, protocol=pickle.HIGHEST_PROTOCOL)
    return manifest


def load_corpus_snapshot(path: Path) -> tuple[Corpus, IngestStats, CorpusManifest]:
    """
    Load a previously-saved corpus snapshot.

    Raises FileNotFoundError if absent, and StaleSnapshotError if the snapshot was
    pickled under a different Corpus shape or cache format. Unpickling restores an
    object's stored __dict__ verbatim, so a snapshot predating a new Corpus
    attribute yields an object silently missing it -- which would otherwise fail
    much later, as an AttributeError inside an unrelated ingest call.
    """
    with open(path, "rb") as f:
        blob = pickle.load(f)

    manifest: CorpusManifest = blob["manifest"]
    corpus: Corpus = blob["corpus"]

    stored_schema = getattr(manifest, "corpus_schema_version", "unknown")
    if stored_schema != CORPUS_SCHEMA_VERSION:
        raise StaleSnapshotError(
            f"{path} was built under Corpus schema {stored_schema!r}, but this code "
            f"expects {CORPUS_SCHEMA_VERSION!r}. Delete it and rebuild -- the numbering "
            f"cache is keyed separately and still applies, so the rebuild skips ANARCII."
        )
    stored_cache_format = getattr(manifest, "cache_format_version", "unknown")
    if stored_cache_format != CACHE_FORMAT_VERSION:
        raise StaleSnapshotError(
            f"{path} was built under cache format {stored_cache_format!r}, but this code "
            f"expects {CACHE_FORMAT_VERSION!r}. Delete it and rebuild."
        )

    return corpus, blob["stats"], manifest


def snapshot_matches_source(manifest: CorpusManifest, source_snapshot: SourceSnapshot) -> bool:
    """Whether a loaded snapshot was built from the given (currently-hashed) source file and parser version."""
    return any(
        s["raw_checksum"] == source_snapshot.raw_checksum and s["parser_version"] == source_snapshot.parser_version
        for s in manifest.source_snapshots
    )


def ingest_sabdab2_cached(csv_path: Path, cache_dir: Path) -> tuple[Corpus, IngestStats, SourceSnapshot]:
    """
    Ingest SAbDab2, reusing a cached numbering result if one exists for this exact
    source file + ANARCII version. Always re-runs parse+index (sub-second), so
    changes to corpus.py's indexing logic take effect immediately; only ANARCII
    numbering itself is skipped when the cache hits.
    """
    from para_scope.corpus import build_sabdab2_snapshot

    snapshot = build_sabdab2_snapshot(csv_path)
    cached = load_numbering_cache(cache_dir, snapshot.raw_checksum)

    corpus = Corpus()
    stats = IngestStats()
    corpus.ingest_sabdab2(csv_path, stats=stats, numbering_cache=cached)

    if cached is None or corpus.last_numbering_cache.keys() != cached.keys():
        save_numbering_cache(cache_dir, snapshot.raw_checksum, corpus.last_numbering_cache)

    return corpus, stats, snapshot
