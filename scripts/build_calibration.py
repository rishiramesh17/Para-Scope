"""
Gate 9: build empirical reference distributions from the corpus snapshot.

Produces data/cache/calibration_v1.json, consumed read-only at scoring time.
Scoring never computes a distribution and never touches live data, so a corpus
update yields a *new* calibration version rather than silently changing the
percentile a stored analysis already reported.

Measurement cost differs by orders of magnitude, so sampling differs too and is
recorded per distribution rather than assumed:

  cheap, computed over every component -- feature values, motif counts, and
  human-repertoire window coverage are per-sequence arithmetic or dict lookups.

  expensive, computed over a deterministic sample -- germline identity needs one
  alignment per reference gene per component (~196 for heavy), so the full corpus
  would be ~1M alignments. A fixed-seed sample is used and its size is stored on
  the distribution, where MINIMUM_SAMPLE_SIZE can see it.
"""

import json
import random
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path

from para_scope.calibration import (
    CALIBRATION_VERSION,
    LENGTH_STRATIFIED_CLASSES,
    build_distribution,
    length_band,
    save_calibration,
    stratum_key,
)
from para_scope.features import NonStandardResidueError, compute_features
from para_scope.germline import SUPPORTED_TYPES, compute_germline_assignment
from para_scope.models import ComponentType, EvidenceStatus
from para_scope.motifs import MOTIF_RULES, find_motifs
from para_scope.persistence import load_corpus_snapshot
from para_scope.window_coverage import compute_window_coverage

SNAPSHOT_PATH = Path("data/cache/combined_corpus_snapshot.pkl")
OUT_PATH = Path("data/cache/calibration_v1.json")

GERMLINE_SAMPLE_PER_CLASS = 700
SAMPLE_SEED = 20261004

POPULATION = (
    "SAbDab2 structural entries plus Thera-SAbDab therapeutic entries, "
    "deduplicated by component. Antibodies selected for crystallisation or "
    "clinical development; NOT a sample of natural antibody repertoire."
)


def _grouped(corpus):
    groups: dict[ComponentType, list] = {}
    for record in corpus._components.values():  # noqa: SLF001 -- offline builder
        groups.setdefault(record.component_type, []).append(record)
    return groups


def _emit(dists, measurement, component_class, records, value_of, snapshots, log):
    """Build one distribution per stratum for a measurement over `records`."""
    by_stratum: dict[str, list[float]] = {}
    strat_meta: dict[str, dict[str, str]] = {}
    for record in records:
        value = value_of(record)
        if value is None:
            continue
        length = len(record.sequence) if component_class in LENGTH_STRATIFIED_CLASSES else None
        key = stratum_key(component_class, length)
        by_stratum.setdefault(key, []).append(float(value))
        strat_meta[key] = (
            {"length_band": length_band(length)} if length is not None else {}
        )
    for key, values in sorted(by_stratum.items()):
        dists.append(
            build_distribution(
                measurement=measurement,
                component_class=component_class,
                stratum=key,
                stratification=strat_meta[key],
                values=values,
                population_description=POPULATION,
                source_snapshots=snapshots,
            )
        )
    if by_stratum:
        total = sum(len(v) for v in by_stratum.values())
        log.append(f"  {measurement:38} {component_class:10} {len(by_stratum):2} strata, {total:6} obs")


def main() -> None:
    t0 = time.time()
    corpus, _stats, manifest = load_corpus_snapshot(SNAPSHOT_PATH)
    snapshots = list(manifest.source_snapshots)
    print(f"Loaded corpus: {len(corpus)} components (schema {manifest.corpus_schema_version})\n")

    groups = _grouped(corpus)
    dists = []
    log: list[str] = []

    # ---- cheap: deterministic features, every component ----
    feature_fields = (
        "sequence_length", "isoelectric_point", "gravy", "aromaticity",
        "net_charge_at_ph_7", "cysteine_count", "n_glycosylation_sequon_count",
    )
    for component_type, records in sorted(groups.items(), key=lambda kv: kv[0].value):
        cached: dict[str, object] = {}
        def feat(record, _cache=cached):
            key = record.record_id
            if key not in _cache:
                try:
                    _cache[key] = compute_features(record.sequence)
                except NonStandardResidueError:
                    _cache[key] = None
            return _cache[key]

        for fname in feature_fields:
            _emit(
                dists, f"feature.{fname}", component_type.value, records,
                lambda r, f=fname: (getattr(feat(r), f) if feat(r) is not None else None),
                snapshots, log,
            )

        # ---- cheap: motif counts, every component ----
        motif_cache: dict[str, object] = {}
        def motifs(record, _cache=motif_cache):
            key = record.record_id
            if key not in _cache:
                _cache[key] = find_motifs(record.sequence, numbered=record.numbered)
            return _cache[key]

        for rule in MOTIF_RULES:
            _emit(
                dists, f"motif.{rule.rule_id}", component_type.value, records,
                lambda r, rid=rule.rule_id: motifs(r).counts_by_rule.get(rid),
                snapshots, log,
            )

    # ---- cheap: human-repertoire window coverage, every eligible component ----
    for component_type, records in sorted(groups.items(), key=lambda kv: kv[0].value):
        def wc(record):
            # record_observations=False: only the fraction is calibrated here, and
            # per-window detail over ~46k components would dominate runtime.
            result = compute_window_coverage(record.sequence, record_observations=False)
            if result.status is not EvidenceStatus.AVAILABLE:
                return None
            return result.coverage
        _emit(
            dists, "window_coverage.coverage", component_type.value,
            records, wc, snapshots, log,
        )

    print("cheap measurements done in %.1fs" % (time.time() - t0))

    # ---- expensive: germline identity, deterministic sample ----
    rng = random.Random(SAMPLE_SEED)
    for component_type in sorted(SUPPORTED_TYPES, key=lambda c: c.value):
        records = groups.get(component_type, [])
        if not records:
            continue
        ordered = sorted(records, key=lambda r: r.record_id)
        sample = (
            ordered if len(ordered) <= GERMLINE_SAMPLE_PER_CLASS
            else rng.sample(ordered, GERMLINE_SAMPLE_PER_CLASS)
        )
        t1 = time.time()
        cache: dict[str, object] = {}
        for record in sample:
            res = compute_germline_assignment(component_type, record.sequence)
            cache[record.record_id] = (
                res.nearest_germline_identity
                if res.status is EvidenceStatus.AVAILABLE else None
            )
        _emit(
            dists, "germline.nearest_identity", component_type.value, sample,
            lambda r: cache.get(r.record_id), snapshots, log,
        )
        print(f"  germline {component_type.value}: {len(sample)} sampled in {time.time()-t1:.1f}s")

    for line in log:
        print(line)

    built_at = datetime.now(timezone.utc).isoformat()
    save_calibration(dists, built_at, snapshots, OUT_PATH)
    size_mb = OUT_PATH.stat().st_size / 1e6
    print(f"\nWrote {len(dists)} distributions to {OUT_PATH} ({size_mb:.1f} MB)")
    print(f"calibration_version={CALIBRATION_VERSION} built_at={built_at}")
    print(f"total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
