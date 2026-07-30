"""
Corpus ingestion pipeline: parse -> deduplicate -> number -> index (SPEC S13-16).

Sequences are deduplicated by checksum *before* IMGT numbering, so a sequence
observed in many source rows (e.g. the same antibody crystallized repeatedly) is
only ever numbered once. Numbering itself runs as one batched ANARCII call over
every unique sequence in the corpus, not one call per row (see numbering.number_sequences).
"""

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

from antigens.indexes import ApproximateIndex, ExactIndex, PairIndex, PairRecord, PositionalIndex
from antigens.indexes.approximate import ELIGIBLE_TYPES as APPROXIMATE_ELIGIBLE_TYPES
from antigens.models import ComponentType, SourceSnapshot
from antigens.numbering import (
    HEAVY_CHAIN_TYPES,
    IMGT_REGION_BOUNDARIES,
    LIGHT_CHAIN_TYPES,
    NumberedComponent,
    compute_sequence_checksum,
    extract_region_sequence,
    number_sequences,
)
from antigens.parsers.sabdab2 import (
    PARSER_VERSION as SABDAB2_PARSER_VERSION,
)
from antigens.parsers.sabdab2 import (
    SAbDab2Row,
    compute_file_checksum,
    parse_ab_split_csv,
)
from antigens.parsers.thera_sabdab import (
    PARSER_VERSION as THERA_SABDAB_PARSER_VERSION,
)
from antigens.parsers.thera_sabdab import (
    TheraSAbDabRow,
    compute_file_checksum as compute_thera_sabdab_checksum,
    parse_thera_sabdab_csv,
)

# Shape of the Corpus object as pickled into a snapshot (see persistence.py).
# Bump on any change to Corpus's instance attributes: a snapshot pickled under an
# older shape restores an object missing whatever attributes were added since,
# which surfaces later as an AttributeError deep inside an ingest call. Loading
# must reject the stale snapshot up front instead.
# The check covers nested index objects too, since they are pickled as part of
# the same graph: adding an attribute to ApproximateIndex invalidates a snapshot
# exactly as adding one to Corpus does.
#   2.0.0 -- added _therapeutic_name_to_arm_count (Thera-SAbDab ingestion, Gate 2)
#   3.0.0 -- ApproximateIndex gained shortlist_index (two-stage search, Gate 3)
#   4.0.0 -- IMGT gap placeholders excluded from region sequence content. This
#            changes region checksums, so component identity and deduplication
#            differ from any earlier snapshot; they are not comparable.
CORPUS_SCHEMA_VERSION = "4.0.0"

REGION_TO_COMPONENT_TYPE: dict[tuple[str, str], ComponentType] = {
    ('H', 'FR1'): ComponentType.FR_H1,
    ('H', 'CDR1'): ComponentType.CDR_H1,
    ('H', 'FR2'): ComponentType.FR_H2,
    ('H', 'CDR2'): ComponentType.CDR_H2,
    ('H', 'FR3'): ComponentType.FR_H3,
    ('H', 'CDR3'): ComponentType.CDR_H3,
    ('H', 'FR4'): ComponentType.FR_H4,
    ('L', 'FR1'): ComponentType.FR_L1,
    ('L', 'CDR1'): ComponentType.CDR_L1,
    ('L', 'FR2'): ComponentType.FR_L2,
    ('L', 'CDR2'): ComponentType.CDR_L2,
    ('L', 'FR3'): ComponentType.FR_L3,
    ('L', 'CDR3'): ComponentType.CDR_L3,
    ('L', 'FR4'): ComponentType.FR_L4,
}

ComponentTypeResolver = Callable[[str], ComponentType]


@dataclass
class CorpusComponentRecord:
    """One deduplicated, numbered component stored in the corpus."""
    record_id: str
    component_type: ComponentType
    sequence: str  # verbatim source sequence
    sequence_checksum: str
    numbered: NumberedComponent
    source_pointers: list[tuple[str, str]] = field(default_factory=list)  # (source_database, source_accession)


@dataclass
class IngestStats:
    rows_seen: int = 0
    rows_skipped_unsupported: int = 0
    heavy_chains_parsed: int = 0
    light_chains_parsed: int = 0
    numbering_failures: int = 0
    unique_components: int = 0
    observed_pairs: int = 0
    unique_parent_antibodies: int = 0
    warnings: list[str] = field(default_factory=list)


def _resolve_light_component_type(chain_type: str) -> ComponentType:
    return ComponentType.VL_KAPPA if chain_type == 'K' else ComponentType.VL_LAMBDA


def _sub_numbered(parent: NumberedComponent, region_name: str, region_seq: str) -> NumberedComponent:
    """Build a NumberedComponent for one IMGT sub-region, reusing the parent's numbering/provenance.

    `sub_numbering` retains the parent's gap placeholders: IMGT position slots that
    a given chain leaves empty are real positional information, and the positional
    index needs them. `region_seq` must NOT contain them -- see extract_region_sequence.
    """
    start, end = IMGT_REGION_BOUNDARIES[region_name]
    sub_numbering = [
        ((pos, code), res) for (pos, code), res in parent.imgt_numbering if start <= pos <= end
    ]
    checksum = compute_sequence_checksum(region_seq)
    return NumberedComponent(
        original_sequence=region_seq,
        sequence_checksum=checksum,
        imgt_numbering=sub_numbering,
        chain_type=parent.chain_type,
        scheme=parent.scheme,
        numbering_score=parent.numbering_score,
        numbering_error=None,
        tool_name=parent.tool_name,
        tool_version=parent.tool_version,
    )


class Corpus:
    """Holds deduplicated components plus all four search indexes (SPEC S16)."""

    def __init__(self):
        # Stamped onto the instance (not just the module) so it is pickled with a
        # corpus snapshot. persistence.py compares it on load and refuses a
        # snapshot built under a different shape, rather than silently restoring
        # an object missing attributes a later gate added. Bump this whenever the
        # set of instance attributes below changes.
        self.schema_version = CORPUS_SCHEMA_VERSION

        self.exact = ExactIndex()
        self.approximate = ApproximateIndex()
        self.positional = PositionalIndex()
        self.pairs = PairIndex()

        self._components: dict[str, CorpusComponentRecord] = {}
        self._checksum_to_record_id: dict[tuple[ComponentType, str], str] = {}
        self._sabdab_id_to_pdb_ids: dict[str, set[str]] = {}
        self._therapeutic_name_to_arm_count: dict[str, int] = {}
        # Numbering cache used by (or extended during) the most recent ingest_sabdab2
        # call. Exposed so callers (persistence.py) can save it to disk and skip
        # ANARCII entirely on a future ingest of the same source data -- measured
        # at 511s of a 512s full ingest, versus 6.85s once the cache hits.
        self.last_numbering_cache: dict[str, Optional[NumberedComponent]] = {}

    def get(self, record_id: str) -> Optional[CorpusComponentRecord]:
        return self._components.get(record_id)

    def __len__(self) -> int:
        return len(self._components)

    def _add_component(
        self,
        component_type: ComponentType,
        sequence: str,
        numbered: NumberedComponent,
        source_database: str,
        source_accession: str,
    ) -> str:
        """Add (or merge provenance into an existing) deduplicated component."""
        checksum = numbered.sequence_checksum
        key = (component_type, checksum)

        existing_id = self._checksum_to_record_id.get(key)
        if existing_id is not None:
            record = self._components[existing_id]
            if (source_database, source_accession) not in record.source_pointers:
                record.source_pointers.append((source_database, source_accession))
            return existing_id

        record_id = f"{component_type.value}::{checksum}"
        record = CorpusComponentRecord(
            record_id=record_id,
            component_type=component_type,
            sequence=sequence,
            sequence_checksum=checksum,
            numbered=numbered,
            source_pointers=[(source_database, source_accession)],
        )
        self._components[record_id] = record
        self._checksum_to_record_id[key] = record_id

        self.exact.add(checksum, record_id)
        if component_type in APPROXIMATE_ELIGIBLE_TYPES:
            self.approximate.add(component_type, record_id, sequence)
        positions = {pos for (pos, _code), _res in numbered.imgt_numbering}
        self.positional.add(component_type, record_id, positions)

        return record_id

    def _add_chain(
        self,
        numbered: NumberedComponent,
        resolve_component_type: ComponentTypeResolver,
        source_database: str,
        source_accession: str,
    ) -> None:
        """Index one already-numbered chain's variable domain plus all IMGT sub-regions."""
        component_type = resolve_component_type(numbered.chain_type)
        self._add_component(component_type, numbered.original_sequence, numbered, source_database, source_accession)

        chain_label = 'H' if numbered.chain_type == 'H' else 'L'
        for region_name in IMGT_REGION_BOUNDARIES:
            region_seq = extract_region_sequence(numbered, region_name)
            if not region_seq:
                continue
            region_component_type = REGION_TO_COMPONENT_TYPE[(chain_label, region_name)]
            region_numbered = _sub_numbered(numbered, region_name, region_seq)
            self._add_component(region_component_type, region_seq, region_numbered, source_database, source_accession)

    def ingest_sabdab2(
        self,
        csv_path: Path,
        stats: Optional[IngestStats] = None,
        numbering_cache: Optional[dict[str, Optional[NumberedComponent]]] = None,
    ) -> IngestStats:
        """
        Ingest SAbDab2's ab_split.csv: parse, batch-number, deduplicate, and index every row.

        Runs in two passes: (1) read all rows and collect every unique heavy/light
        sequence, (2) batch-number them once, then walk the rows again to build
        components, sub-regions, and directly-observed pairs from the cache.

        If numbering_cache is supplied (e.g. loaded from disk by persistence.py),
        only sequences missing from it are sent to ANARCII -- reused entries skip
        renumbering entirely. The resulting cache (reused + newly-numbered) is
        exposed on self.last_numbering_cache for the caller to persist.
        """
        if stats is None:
            stats = IngestStats()

        rows: list[SAbDab2Row] = []
        all_sequences: list[str] = []
        for row in parse_ab_split_csv(csv_path):
            stats.rows_seen += 1
            if not row.supported:
                stats.rows_skipped_unsupported += 1
                stats.warnings.append(f"{row.instance_id}: {row.skip_reason}")
                continue
            rows.append(row)
            if row.heavy_sequence:
                all_sequences.append(row.heavy_sequence)
            if row.light_sequence:
                all_sequences.append(row.light_sequence)

        if numbering_cache is None:
            numbering_cache = number_sequences(all_sequences)
        else:
            missing = [
                seq for seq in all_sequences
                if compute_sequence_checksum(seq) not in numbering_cache
            ]
            if missing:
                numbering_cache = {**numbering_cache, **number_sequences(missing)}

        self.last_numbering_cache = numbering_cache

        for row in rows:
            pdb_ids = self._sabdab_id_to_pdb_ids.setdefault(row.sabdab_id, set())
            pdb_ids.add(row.pdb_id)

            heavy_numbered = self._resolve_chain(row.heavy_sequence, numbering_cache, HEAVY_CHAIN_TYPES, row, stats)
            if heavy_numbered is not None:
                stats.heavy_chains_parsed += 1
                heavy_component_type = ComponentType.VHH if row.antibody_type == "SD-H" else ComponentType.VH
                self._add_chain(heavy_numbered, lambda _ct: heavy_component_type, "SAbDab2", row.instance_id)

            light_numbered = self._resolve_chain(row.light_sequence, numbering_cache, LIGHT_CHAIN_TYPES, row, stats)
            if light_numbered is not None:
                stats.light_chains_parsed += 1
                self._add_chain(light_numbered, _resolve_light_component_type, "SAbDab2", row.instance_id)

            if heavy_numbered is not None and light_numbered is not None:
                self.pairs.add(
                    heavy_numbered.sequence_checksum,
                    light_numbered.sequence_checksum,
                    PairRecord(
                        source_database="SAbDab2",
                        source_accession=row.instance_id,
                        structure_available=True,
                        metadata={
                            "pdb_id": row.pdb_id,
                            "method": row.method,
                            "resolution": row.resolution,
                            "holo": row.holo,
                            "construct": row.construct,
                        },
                    ),
                )
                stats.observed_pairs += 1

        stats.unique_components = len(self._components)
        stats.unique_parent_antibodies = len(self._sabdab_id_to_pdb_ids) + len(self._therapeutic_name_to_arm_count)
        return stats

    def ingest_thera_sabdab(
        self,
        csv_path: Path,
        stats: Optional[IngestStats] = None,
        numbering_cache: Optional[dict[str, Optional[NumberedComponent]]] = None,
    ) -> IngestStats:
        """
        Ingest Thera-SAbDab's therapeutic sequence download: parse, batch-number, and index.

        Each row may contribute up to two independent variable-domain arms (a
        primary arm, and -- for bispecifics -- a second arm from the
        HeavySequence(ifbispec)/LightSequence(ifbispec) columns). Arms are indexed
        and paired independently; a bispecific's second-arm pairing is never
        merged with the primary arm's pairing (CLAUDE.md: paired evidence requires
        both chains from the same source record).

        FULL_HEAVY / FULL_LIGHT are not populated by this source either -- Thera-SAbDab's
        download is variable-domain-only, confirmed by direct inspection of the real file
        (see parsers/thera_sabdab.py module docstring).
        """
        if stats is None:
            stats = IngestStats()

        rows: list[TheraSAbDabRow] = []
        all_sequences: list[str] = []
        for row in parse_thera_sabdab_csv(csv_path):
            stats.rows_seen += 1
            if not row.supported:
                stats.rows_skipped_unsupported += 1
                stats.warnings.append(f"{row.therapeutic_name}: {row.skip_reason}")
                continue
            rows.append(row)
            for arm in row.arms:
                if arm.heavy_sequence:
                    all_sequences.append(arm.heavy_sequence)
                if arm.light_sequence:
                    all_sequences.append(arm.light_sequence)

        if numbering_cache is None:
            numbering_cache = number_sequences(all_sequences)
        else:
            missing = [
                seq for seq in all_sequences
                if compute_sequence_checksum(seq) not in numbering_cache
            ]
            if missing:
                numbering_cache = {**numbering_cache, **number_sequences(missing)}
        self.last_numbering_cache = numbering_cache

        for row in rows:
            arm_count = self._therapeutic_name_to_arm_count.get(row.therapeutic_name, 0)

            for arm in row.arms:
                accession = f"{row.therapeutic_name}#arm{arm.arm_index}"
                warning_key = accession

                heavy_numbered = self._resolve_chain_by_key(
                    arm.heavy_sequence, numbering_cache, HEAVY_CHAIN_TYPES, warning_key, stats
                )
                heavy_component_type = ComponentType.VHH if arm.is_single_domain else ComponentType.VH
                if heavy_numbered is not None:
                    stats.heavy_chains_parsed += 1
                    self._add_chain(
                        heavy_numbered, lambda _ct: heavy_component_type, "Thera-SAbDab", accession
                    )

                light_numbered = None
                if arm.light_sequence:
                    light_numbered = self._resolve_chain_by_key(
                        arm.light_sequence, numbering_cache, LIGHT_CHAIN_TYPES, warning_key, stats
                    )
                    if light_numbered is not None:
                        stats.light_chains_parsed += 1
                        self._add_chain(light_numbered, _resolve_light_component_type, "Thera-SAbDab", accession)

                if heavy_numbered is not None and light_numbered is not None:
                    self.pairs.add(
                        heavy_numbered.sequence_checksum,
                        light_numbered.sequence_checksum,
                        PairRecord(
                            source_database="Thera-SAbDab",
                            source_accession=accession,
                            structure_available=bool(
                                row.structure_100_si or row.structure_99_si or row.structure_95_98_si
                            ),
                            metadata={
                                "therapeutic_name": row.therapeutic_name,
                                "format": row.format_raw,
                                "development_status": row.development_status,
                                "highest_clinical_trial": row.highest_clinical_trial,
                                "target": row.target,
                                "arm_index": arm.arm_index,
                            },
                        ),
                    )
                    stats.observed_pairs += 1

            arm_count += len(row.arms)
            self._therapeutic_name_to_arm_count[row.therapeutic_name] = arm_count

        stats.unique_components = len(self._components)
        stats.unique_parent_antibodies = len(self._sabdab_id_to_pdb_ids) + len(self._therapeutic_name_to_arm_count)
        return stats

    @staticmethod
    def _resolve_chain(
        sequence: Optional[str],
        numbering_cache: dict[str, Optional[NumberedComponent]],
        allowed_chain_types: tuple[str, ...],
        row: SAbDab2Row,
        stats: IngestStats,
    ) -> Optional[NumberedComponent]:
        return Corpus._resolve_chain_by_key(sequence, numbering_cache, allowed_chain_types, row.instance_id, stats)

    @staticmethod
    def _resolve_chain_by_key(
        sequence: Optional[str],
        numbering_cache: dict[str, Optional[NumberedComponent]],
        allowed_chain_types: tuple[str, ...],
        warning_key: str,
        stats: IngestStats,
    ) -> Optional[NumberedComponent]:
        """Shared by every source parser: resolve one chain sequence from the numbering
        cache, recording a warning against warning_key (a row/entry identifier) on failure."""
        if not sequence:
            return None
        checksum = compute_sequence_checksum(sequence)
        numbered = numbering_cache.get(checksum)
        if numbered is None:
            stats.numbering_failures += 1
            reason = (
                "invalid input (non-standard residue, e.g. 'X')"
                if checksum not in numbering_cache
                else "ANARCII numbering failed"
            )
            stats.warnings.append(f"{warning_key}: {reason}")
            return None
        if numbered.chain_type not in allowed_chain_types:
            stats.numbering_failures += 1
            stats.warnings.append(
                f"{warning_key}: chain-type mismatch (got '{numbered.chain_type}', "
                f"expected one of {allowed_chain_types})"
            )
            return None
        return numbered


def build_sabdab2_snapshot(csv_path: Path, retrieved_at: Optional[str] = None) -> SourceSnapshot:
    """Construct a SourceSnapshot for the SAbDab2 ab_split.csv file (SPEC S13-14)."""
    with open(csv_path, newline="", encoding="utf-8") as f:
        record_count = sum(1 for _ in f) - 1  # exclude header

    return SourceSnapshot(
        source_name="SAbDab2",
        snapshot_id="zenodo:20083995:ab_split.csv",
        retrieved_at=retrieved_at or datetime.now(timezone.utc).isoformat(),
        license_id="unspecified",  # TODO: confirm SAbDab2's license before any external redistribution
        raw_checksum=compute_file_checksum(csv_path),
        record_count=record_count,
        parser_version=SABDAB2_PARSER_VERSION,
        schema_version="sabdab2-ml-2026.06",
    )


def build_thera_sabdab_snapshot(csv_path: Path, retrieved_at: Optional[str] = None) -> SourceSnapshot:
    """Construct a SourceSnapshot for the Thera-SAbDab therapeutic download (SPEC S13-14)."""
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        record_count = sum(1 for _ in f) - 1  # exclude header

    return SourceSnapshot(
        source_name="Thera-SAbDab",
        snapshot_id="opig:therasabdab:TheraSAbDab_SeqStruc_OnlineDownload.csv:2025-03-18",
        retrieved_at=retrieved_at or datetime.now(timezone.utc).isoformat(),
        license_id="unspecified",  # TODO: confirm Thera-SAbDab's license before any external redistribution
        raw_checksum=compute_thera_sabdab_checksum(csv_path),
        record_count=record_count,
        parser_version=THERA_SABDAB_PARSER_VERSION,
        schema_version="thera-sabdab-2025-03-18",
    )
