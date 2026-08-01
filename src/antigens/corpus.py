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
#   5.0.0 -- therapeutic entity store (Gate 7): _therapeutics, _therapeutic_aliases,
#            _component_to_therapeutics replace _therapeutic_name_to_arm_count.
#            Also changes PairRecord.structure_available from bool to Optional[bool]
#            and fixes its Thera-SAbDab value, so pair metadata from an older
#            snapshot is not comparable either.
CORPUS_SCHEMA_VERSION = "5.0.0"

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
class TherapeuticArmRecord:
    """One variable-domain arm of a therapeutic, linked to the components it contributed.

    Holds record ids rather than sequences: the components themselves are already
    deduplicated in the corpus, and two therapeutics sharing an identical VH must
    resolve to the same component record while remaining two distinct therapeutics.
    """
    arm_index: int
    accession: str  # "<therapeutic name>#arm<n>", matching source_pointers
    heavy_record_id: Optional[str] = None
    light_record_id: Optional[str] = None
    is_single_domain: bool = False
    structures_100_si: Optional[list[str]] = None
    structures_99_si: Optional[list[str]] = None
    structures_95_98_si: Optional[list[str]] = None

    @property
    def has_structure_record(self) -> Optional[bool]:
        """True/False when any tier is known; None when every tier is unknown."""
        tiers = (self.structures_100_si, self.structures_99_si, self.structures_95_98_si)
        if all(tier is None for tier in tiers):
            return None
        return any(tier for tier in tiers)

    @property
    def structure_ids(self) -> list[str]:
        """Every known structure reference for this arm, highest identity tier first."""
        ordered: list[str] = []
        for tier in (self.structures_100_si, self.structures_99_si, self.structures_95_98_si):
            for ref in tier or []:
                if ref not in ordered:
                    ordered.append(ref)
        return ordered


@dataclass
class TherapeuticRecord:
    """
    One canonical therapeutic entity and its curated metadata (SPEC S15.5, S26).

    Exists because therapeutic metadata was previously reachable only through
    PairRecord.metadata, which meant a therapeutic whose arm has no light chain --
    every single-domain therapeutic among them -- carried no recoverable name,
    target, format, or development status at all. SPEC S26 requires all of those as
    display fields for any supported component, so they cannot depend on the arm
    happening to be paired.

    `highest_clinical_trial` and `development_status` are metadata and nothing more.
    An antibody reaching Phase III transfers no property whatever to a sequence that
    resembles it, and neither field may become a success or failure label.
    """
    therapeutic_name: str  # canonical
    format_raw: str        # free text, verbatim, never classified beyond is_single_domain
    alternative_names: list[str] = field(default_factory=list)
    ch1_isotype: Optional[str] = None
    light_chain_class_label: Optional[str] = None
    highest_clinical_trial: Optional[str] = None
    development_status: Optional[str] = None
    target: Optional[str] = None
    arms: list[TherapeuticArmRecord] = field(default_factory=list)

    # Fields SPEC S26 lists under "therapeutic-record completeness". Kept as an
    # explicit tuple so the denominator is auditable rather than implied by
    # whichever attributes a future reader thinks ought to count.
    COMPLETENESS_FIELDS = (
        "target", "format_raw", "ch1_isotype", "light_chain_class_label",
        "highest_clinical_trial", "development_status",
    )

    @property
    def metadata_completeness(self) -> float:
        """Fraction of SPEC S26's metadata fields this record actually carries."""
        present = sum(1 for name in self.COMPLETENESS_FIELDS if getattr(self, name))
        return present / len(self.COMPLETENESS_FIELDS)

    @property
    def has_structure_record(self) -> Optional[bool]:
        """Whether any arm is structure-backed; None only if every arm is unknown."""
        per_arm = [arm.has_structure_record for arm in self.arms]
        if all(value is None for value in per_arm):
            return None
        return any(per_arm)


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
    # Aliases claimed by more than one therapeutic. Such an alias is registered to
    # neither, since resolving it either way would merge two distinct entities.
    therapeutic_alias_collisions: int = 0
    warnings: list[str] = field(default_factory=list)


# Values that appear in the alternative-names column but are placeholders rather
# than names. Registering one as an alias would let unrelated therapeutics collide
# on it, and would let a user searching "unknown" resolve to a drug.
_ALIAS_PLACEHOLDERS = frozenset({"na", "n/a", "nd", "nfd", "tbc", "unknown", "none", "-"})


def split_alternative_names(raw: Optional[str]) -> list[str]:
    """
    Split Thera-SAbDab's semicolon-delimited alternative-names field.

    Semicolons inside parentheses do not delimit: combination products are written
    as 'REGN5713-5714-5715 (combination of Atisnolerbart; Bremzalerbart; and
    Umesolerbart)', and a naive split yields the fragments 'Bremzalerbart' and
    'and Umesolerbart)'. The first then collides with the real therapeutic of that
    name and poisons a legitimate alias; the second is not a name at all.

    The combination identifier itself is kept whole. It genuinely refers to several
    therapeutics, so it collides by design -- and being left unresolved is the
    correct outcome for a name that does not denote one entity (SPEC S15.5).
    """
    if not raw:
        return []

    names: list[str] = []
    current: list[str] = []
    depth = 0
    for character in raw:
        if character == "(":
            depth += 1
        elif character == ")":
            depth = max(0, depth - 1)
        if character == ";" and depth == 0:
            names.append("".join(current))
            current = []
        else:
            current.append(character)
    names.append("".join(current))

    return [
        name.strip() for name in names
        if name.strip() and name.strip().lower() not in _ALIAS_PLACEHOLDERS
    ]


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
        # Canonical therapeutic name -> record (SPEC S15.5: one entity per therapeutic).
        self._therapeutics: dict[str, TherapeuticRecord] = {}
        # Lowercased alias/development code -> canonical name. Ambiguous aliases are
        # absent rather than arbitrarily resolved.
        self._therapeutic_aliases: dict[str, str] = {}
        # Aliases proven ambiguous. Kept explicitly so a third claimant cannot
        # re-register an alias that a collision already removed from the map.
        self._ambiguous_therapeutic_aliases: set[str] = set()
        # Component record id -> canonical therapeutic names that contributed it. The
        # reverse edge Gate 7 needs: an alignment hit knows its record id and nothing else.
        self._component_to_therapeutics: dict[str, set[str]] = {}
        # Numbering cache used by (or extended during) the most recent ingest_sabdab2
        # call. Exposed so callers (persistence.py) can save it to disk and skip
        # ANARCII entirely on a future ingest of the same source data -- measured
        # at 511s of a 512s full ingest, versus 6.85s once the cache hits.
        self.last_numbering_cache: dict[str, Optional[NumberedComponent]] = {}

    def get(self, record_id: str) -> Optional[CorpusComponentRecord]:
        return self._components.get(record_id)

    def __len__(self) -> int:
        return len(self._components)

    def get_therapeutic(self, name: str) -> Optional[TherapeuticRecord]:
        """Resolve a therapeutic by canonical name, alias, or development code (SPEC S15.5)."""
        record = self._therapeutics.get(name)
        if record is not None:
            return record
        canonical = self._therapeutic_aliases.get(name.strip().lower())
        return self._therapeutics.get(canonical) if canonical else None

    def therapeutics_for_component(self, record_id: str) -> list[TherapeuticRecord]:
        """Every therapeutic that contributed this exact component.

        More than one is normal and is not a duplicate: distinct therapeutics can
        share an identical framework or CDR, and each is independent evidence. The
        component is deduplicated; the therapeutics behind it are not.
        """
        names = self._component_to_therapeutics.get(record_id, set())
        return [self._therapeutics[name] for name in sorted(names) if name in self._therapeutics]

    @property
    def therapeutics(self) -> dict[str, TherapeuticRecord]:
        return self._therapeutics

    def _register_therapeutic_aliases(
        self, record: TherapeuticRecord, stats: IngestStats
    ) -> None:
        """Map a therapeutic's names and development codes onto its canonical entity.

        An alias claimed by two different therapeutics is registered to neither and
        counted: silently letting the later row win would make lookups resolve to
        whichever happened to be parsed last, and quietly merge two entities that
        SPEC S15.5 requires be kept distinct.
        """
        canonical = record.therapeutic_name
        for alias in [canonical, *record.alternative_names]:
            key = alias.strip().lower()
            if not key or key in self._ambiguous_therapeutic_aliases:
                continue
            existing = self._therapeutic_aliases.get(key)
            if existing is None:
                self._therapeutic_aliases[key] = canonical
            elif existing != canonical:
                del self._therapeutic_aliases[key]
                self._ambiguous_therapeutic_aliases.add(key)
                stats.therapeutic_alias_collisions += 1
                stats.warnings.append(
                    f"therapeutic alias {alias!r} claimed by both {existing!r} and "
                    f"{canonical!r}; left unresolved rather than assigned to either"
                )

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
        therapeutic_name: Optional[str] = None,
    ) -> str:
        """Index one already-numbered chain's variable domain plus all IMGT sub-regions.

        Returns the variable domain's record id. `therapeutic_name`, when given,
        links every component this chain contributes -- the domain and each
        sub-region -- back to the therapeutic, so Gate 7 can answer "which
        therapeutics contain this CDR" from an alignment hit alone.
        """
        component_type = resolve_component_type(numbered.chain_type)
        domain_record_id = self._add_component(
            component_type, numbered.original_sequence, numbered, source_database, source_accession
        )
        if therapeutic_name is not None:
            self._component_to_therapeutics.setdefault(domain_record_id, set()).add(therapeutic_name)

        chain_label = 'H' if numbered.chain_type == 'H' else 'L'
        for region_name in IMGT_REGION_BOUNDARIES:
            region_seq = extract_region_sequence(numbered, region_name)
            if not region_seq:
                continue
            region_component_type = REGION_TO_COMPONENT_TYPE[(chain_label, region_name)]
            region_numbered = _sub_numbered(numbered, region_name, region_seq)
            region_record_id = self._add_component(
                region_component_type, region_seq, region_numbered, source_database, source_accession
            )
            if therapeutic_name is not None:
                self._component_to_therapeutics.setdefault(region_record_id, set()).add(therapeutic_name)

        return domain_record_id

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
        stats.unique_parent_antibodies = len(self._sabdab_id_to_pdb_ids) + len(self._therapeutics)
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
            if row.structure_anomalies:
                # Drift, not a skip: the row is still ingested, but the affected
                # structure columns resolved to unknown and that must be visible.
                # A sentinel this parser does not recognise is exactly how 794
                # therapeutics acquired a fabricated structure flag.
                stats.warnings.append(
                    f"{row.therapeutic_name}: unrecognised structure token(s) "
                    f"{sorted(set(row.structure_anomalies))} -- affected columns "
                    f"resolved to unknown"
                )
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
            therapeutic = self._therapeutics.get(row.therapeutic_name)
            if therapeutic is None:
                therapeutic = TherapeuticRecord(
                    therapeutic_name=row.therapeutic_name,
                    format_raw=row.format_raw,
                    alternative_names=split_alternative_names(row.alternative_names),
                    ch1_isotype=row.ch1_isotype,
                    light_chain_class_label=row.light_chain_class_label,
                    highest_clinical_trial=row.highest_clinical_trial,
                    development_status=row.development_status,
                    target=row.target,
                )
                self._therapeutics[row.therapeutic_name] = therapeutic
                self._register_therapeutic_aliases(therapeutic, stats)

            for arm in row.arms:
                accession = f"{row.therapeutic_name}#arm{arm.arm_index}"
                warning_key = accession

                arm_record = TherapeuticArmRecord(
                    arm_index=arm.arm_index,
                    accession=accession,
                    is_single_domain=arm.is_single_domain,
                    structures_100_si=arm.structures_100_si,
                    structures_99_si=arm.structures_99_si,
                    structures_95_98_si=arm.structures_95_98_si,
                )
                therapeutic.arms.append(arm_record)

                heavy_numbered = self._resolve_chain_by_key(
                    arm.heavy_sequence, numbering_cache, HEAVY_CHAIN_TYPES, warning_key, stats
                )
                heavy_component_type = ComponentType.VHH if arm.is_single_domain else ComponentType.VH
                if heavy_numbered is not None:
                    stats.heavy_chains_parsed += 1
                    arm_record.heavy_record_id = self._add_chain(
                        heavy_numbered, lambda _ct: heavy_component_type, "Thera-SAbDab", accession,
                        therapeutic_name=row.therapeutic_name,
                    )

                light_numbered = None
                if arm.light_sequence:
                    light_numbered = self._resolve_chain_by_key(
                        arm.light_sequence, numbering_cache, LIGHT_CHAIN_TYPES, warning_key, stats
                    )
                    if light_numbered is not None:
                        stats.light_chains_parsed += 1
                        arm_record.light_record_id = self._add_chain(
                            light_numbered, _resolve_light_component_type, "Thera-SAbDab", accession,
                            therapeutic_name=row.therapeutic_name,
                        )

                if heavy_numbered is not None and light_numbered is not None:
                    self.pairs.add(
                        heavy_numbered.sequence_checksum,
                        light_numbered.sequence_checksum,
                        PairRecord(
                            source_database="Thera-SAbDab",
                            source_accession=accession,
                            # This arm's own structures, not the row's. OR-ing the
                            # columns credited a bispecific's first arm with the
                            # second arm's experimental records.
                            structure_available=arm_record.has_structure_record,
                            metadata={
                                "therapeutic_name": row.therapeutic_name,
                                "format": row.format_raw,
                                "development_status": row.development_status,
                                "highest_clinical_trial": row.highest_clinical_trial,
                                "target": row.target,
                                "arm_index": arm.arm_index,
                                "structure_ids": arm_record.structure_ids,
                            },
                        ),
                    )
                    stats.observed_pairs += 1

        stats.unique_components = len(self._components)
        stats.unique_parent_antibodies = len(self._sabdab_id_to_pdb_ids) + len(self._therapeutics)
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
