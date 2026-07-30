"""
Gate 1: derive a fixed gold-query benchmark set from the real, ingested SAbDab2 corpus.

Rather than hardcoding real database sequences into a committed file by hand,
this script *selects* them programmatically from an already-built corpus
snapshot (see persistence.py / scripts/build_v2_corpus_cache.py) -- so the
benchmark always reflects what is actually in the corpus, and re-running this
script against the same snapshot reproduces the same file (selection is
deterministic: sorted by record_id, first match wins).

Categories covered (Version 2 roadmap, Gate 1):
  - exact_match              -- a real VH sequence, queried verbatim
  - single_substitution      -- the same VH with exactly one residue changed
  - several_substitution     -- the same VH with five residues changed
  - similar_nonidentical     -- two distinct real VH records (not variants of each other)
  - long_cdr_h3              -- a real component with an unusually long CDR-H3
  - insertion_coded_cdr      -- a real component whose IMGT numbering uses an insertion code
  - heavy_only               -- a real VH queried with no paired light chain
  - light_only               -- a real VL queried with no paired heavy chain
  - paired_vh_vl             -- a directly observed heavy/light pair
  - invalid_sequence         -- contains non-residue characters; must be rejected, not scored
  - truncated_sequence       -- a real VH truncated to 20 residues
  - duplicate_parent_record  -- a parent antibody (SABDAB_ID) observed across >1 PDB structure
  - no_close_match           -- a randomly generated sequence unlikely to resemble any real chain
"""

import json
import random
from pathlib import Path

from antigens.models import ComponentType
from antigens.persistence import load_corpus_snapshot

SNAPSHOT_PATH = Path("data/cache/sabdab2_corpus_snapshot.pkl")
OUTPUT_PATH = Path("benchmarks/search_gold_queries.json")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _mutate(sequence: str, n_substitutions: int, seed: int) -> str:
    rng = random.Random(seed)
    positions = rng.sample(range(1, len(sequence) - 1), min(n_substitutions, len(sequence) - 2))
    chars = list(sequence)
    for pos in positions:
        original = chars[pos]
        replacement = rng.choice([c for c in AMINO_ACIDS if c != original])
        chars[pos] = replacement
    return "".join(chars)


def _has_insertion_code(numbered) -> bool:
    return any(code and code.strip() for (_pos, code), _res in numbered.imgt_numbering)


def build_gold_queries(snapshot_path: Path = SNAPSHOT_PATH) -> dict:
    corpus, stats, manifest = load_corpus_snapshot(snapshot_path)

    vh_records = sorted(
        (r for r in corpus._components.values() if r.component_type == ComponentType.VH),
        key=lambda r: r.record_id,
    )
    vl_records = sorted(
        (r for r in corpus._components.values() if r.component_type in (ComponentType.VL_KAPPA, ComponentType.VL_LAMBDA)),
        key=lambda r: r.record_id,
    )
    cdr_h3_records = sorted(
        (r for r in corpus._components.values() if r.component_type == ComponentType.CDR_H3),
        key=lambda r: r.record_id,
    )

    if len(vh_records) < 3 or len(vl_records) < 1:
        raise ValueError("Corpus snapshot too small to build a meaningful gold-query set")

    exact_record = vh_records[0]
    similar_a, similar_b = vh_records[1], vh_records[2]

    longest_cdr_h3 = max(cdr_h3_records, key=lambda r: len(r.sequence))

    insertion_coded_record = next(
        (r for r in vh_records + vl_records if _has_insertion_code(r.numbered)), None
    )

    light_record = vl_records[0]

    # First directly observed pair in the pair index, resolved back to component records.
    pair_key = next(iter(corpus.pairs._index.keys()), None)
    paired_heavy = paired_light = None
    if pair_key is not None:
        heavy_checksum, light_checksum = pair_key
        for r in vh_records:
            if r.sequence_checksum == heavy_checksum:
                paired_heavy = r
                break
        for r in vl_records:
            if r.sequence_checksum == light_checksum:
                paired_light = r
                break

    # A parent antibody (SABDAB_ID) whose structures span more than one PDB deposition.
    duplicate_sabdab_id = next(
        (sid for sid, pdb_ids in corpus._sabdab_id_to_pdb_ids.items() if len(pdb_ids) > 1), None
    )
    duplicate_record = None
    if duplicate_sabdab_id is not None:
        for r in vh_records:
            if len(r.source_pointers) > 1:
                duplicate_record = r
                break

    rng = random.Random(20260729)
    no_match_sequence = "".join(rng.choice(AMINO_ACIDS) for _ in range(118))

    queries = [
        {
            "category": "exact_match",
            "query_sequence": exact_record.sequence,
            "query_sequence_checksum": exact_record.sequence_checksum,
            "component_type_hint": exact_record.component_type.value,
            "expected_record_id": exact_record.record_id,
            "expect": {"exact_match_rank": 1, "sequence_identity": 1.0},
        },
        {
            "category": "single_substitution",
            "query_sequence": _mutate(exact_record.sequence, 1, seed=1),
            "component_type_hint": exact_record.component_type.value,
            "expected_record_id": exact_record.record_id,
            "expect": {"top_hit_rank": 1, "sequence_identity_min": 1 - 2 / len(exact_record.sequence)},
        },
        {
            "category": "several_substitution",
            "query_sequence": _mutate(exact_record.sequence, 5, seed=2),
            "component_type_hint": exact_record.component_type.value,
            "expected_record_id": exact_record.record_id,
            "expect": {"top_hit_rank": 1, "sequence_identity_min": 1 - 10 / len(exact_record.sequence)},
        },
        {
            "category": "similar_nonidentical",
            "query_sequence": similar_a.sequence,
            "component_type_hint": similar_a.component_type.value,
            "reference_record_id": similar_b.record_id,
            "expect": {"description": "two distinct real VH records; not expected to match each other"},
        },
        {
            "category": "long_cdr_h3",
            "query_sequence": longest_cdr_h3.sequence,
            "query_sequence_checksum": longest_cdr_h3.sequence_checksum,
            "component_type_hint": "CDR-H3",
            "expected_record_id": longest_cdr_h3.record_id,
            "expect": {"length": len(longest_cdr_h3.sequence)},
        },
        {
            "category": "insertion_coded_cdr",
            "query_sequence": insertion_coded_record.sequence if insertion_coded_record else None,
            "query_sequence_checksum": insertion_coded_record.sequence_checksum if insertion_coded_record else None,
            "component_type_hint": insertion_coded_record.component_type.value if insertion_coded_record else None,
            "expected_record_id": insertion_coded_record.record_id if insertion_coded_record else None,
            "expect": {"description": "component numbering must include at least one non-empty insertion code"},
        },
        {
            "category": "heavy_only",
            "query_sequence": vh_records[3].sequence,
            "query_sequence_checksum": vh_records[3].sequence_checksum,
            "component_type_hint": "VH",
            "expected_record_id": vh_records[3].record_id,
            "expect": {"description": "query with no accompanying light chain"},
        },
        {
            "category": "light_only",
            "query_sequence": light_record.sequence,
            "query_sequence_checksum": light_record.sequence_checksum,
            "component_type_hint": light_record.component_type.value,
            "expected_record_id": light_record.record_id,
            "expect": {"description": "query with no accompanying heavy chain"},
        },
        {
            "category": "paired_vh_vl",
            "heavy_sequence": paired_heavy.sequence if paired_heavy else None,
            "heavy_sequence_checksum": paired_heavy.sequence_checksum if paired_heavy else None,
            "light_sequence": paired_light.sequence if paired_light else None,
            "light_sequence_checksum": paired_light.sequence_checksum if paired_light else None,
            "expect": {"description": "must be found via PairIndex as a directly observed pair"},
        },
        {
            "category": "invalid_sequence",
            "query_sequence": "MXKZQV123INVALID",
            "component_type_hint": "VH",
            "expect": {"description": "must be rejected outright, never partially scored"},
        },
        {
            "category": "truncated_sequence",
            "query_sequence": exact_record.sequence[:20],
            "component_type_hint": "VH",
            "expect": {"description": "must fail numbering loudly or be explicitly flagged incomplete, never silently misrepresented"},
        },
        {
            "category": "duplicate_parent_record",
            "query_sequence": duplicate_record.sequence if duplicate_record else None,
            "query_sequence_checksum": duplicate_record.sequence_checksum if duplicate_record else None,
            "expected_record_id": duplicate_record.record_id if duplicate_record else None,
            "expected_sabdab_id": duplicate_sabdab_id,
            "expect": {"description": "source_pointers must list >1 structure; unique_parent_antibodies must count it once"},
        },
        {
            "category": "no_close_match",
            "query_sequence": no_match_sequence,
            "component_type_hint": "VH",
            "expect": {"description": "random sequence; must not produce a spurious high-identity top hit"},
        },
    ]

    return {
        "source_snapshot": {
            "raw_checksum": manifest.source_snapshots[0]["raw_checksum"],
            "parser_version": manifest.source_snapshots[0]["parser_version"],
            "created_at": manifest.created_at,
        },
        "corpus_stats_at_build_time": manifest.stats,
        "queries": queries,
    }


def main() -> None:
    gold = build_gold_queries()
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(gold, indent=2))
    print(f"Wrote {len(gold['queries'])} gold queries to {OUTPUT_PATH}")
    for q in gold["queries"]:
        print(f"  - {q['category']}")


if __name__ == "__main__":
    main()
