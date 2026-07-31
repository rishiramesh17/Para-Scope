"""Pair index (SPEC S16.5): directly observed heavy-light chain pairs only.

CLAUDE.md: Paired and unpaired evidence never merge. Observed pairing requires
both chains linked in the *same source record*. This index must never be
populated by joining an independent heavy hit with an independent light hit.
"""

from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class PairRecord:
    source_database: str
    source_accession: str
    # Three-state. None means "the source did not let us determine this", which is
    # not False. It was previously a plain bool, and Thera-SAbDab's structure
    # columns use the sentinel 'None' rather than 'na', so every therapeutic pair
    # was recorded as structure-backed -- 794 of 1,131 of them falsely. A consumer
    # must not coerce None to False (CLAUDE.md: missing evidence is not zero
    # evidence); SPEC S26 scores this field directly.
    structure_available: Optional[bool]
    study_id: Optional[str] = None
    metadata: dict[str, Any] = field(default_factory=dict)


class PairIndex:
    """Maps (heavy_checksum, light_checksum) -> observed pair records."""

    def __init__(self):
        self._index: dict[tuple[str, str], list[PairRecord]] = {}

    def add(self, heavy_checksum: str, light_checksum: str, pair_record: PairRecord) -> None:
        key = (heavy_checksum, light_checksum)
        self._index.setdefault(key, []).append(pair_record)

    def lookup(self, heavy_checksum: str, light_checksum: str) -> list[PairRecord]:
        return list(self._index.get((heavy_checksum, light_checksum), []))

    def __len__(self) -> int:
        return len(self._index)
