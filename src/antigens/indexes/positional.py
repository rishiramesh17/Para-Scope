"""Numbered-position index (SPEC S16.3): compare components by shared IMGT positions.

Only components of the same ComponentType are ever indexed together — this index
never crosses the CDR-H3/CDR-L3 or VH/VL boundary (CLAUDE.md: only compare like with like).
"""

from collections import defaultdict

from antigens.models import ComponentType

DEFAULT_OVERLAP_THRESHOLD = 0.95


class PositionalIndex:
    """Per component type, an inverted index from IMGT position to record ids."""

    def __init__(self):
        # Plain (not defaultdict-of-defaultdict) so the index stays picklable --
        # a defaultdict with a lambda factory cannot be serialized (Gate 0.5).
        self._position_to_records: dict[ComponentType, dict[int, set[str]]] = {}
        self._record_positions: dict[str, set[int]] = {}

    def add(self, component_type: ComponentType, record_id: str, positions: set[int]) -> None:
        self._record_positions[record_id] = positions
        by_position = self._position_to_records.setdefault(component_type, {})
        for pos in positions:
            by_position.setdefault(pos, set()).add(record_id)

    def lookup(
        self,
        component_type: ComponentType,
        positions: set[int],
        threshold: float = DEFAULT_OVERLAP_THRESHOLD,
    ) -> list[tuple[str, float]]:
        """Return (record_id, overlap_fraction) for same-type records at or above threshold.

        overlap_fraction is the Jaccard similarity of the two position sets:
        |shared positions| / |union of positions|.
        """
        if not positions:
            return []

        position_map = self._position_to_records.get(component_type, {})
        shared_counts: dict[str, int] = defaultdict(int)
        for pos in positions:
            for record_id in position_map.get(pos, ()):
                shared_counts[record_id] += 1

        results = []
        for record_id, shared in shared_counts.items():
            record_positions = self._record_positions[record_id]
            union_size = len(positions | record_positions)
            overlap = shared / union_size if union_size else 0.0
            if overlap >= threshold:
                results.append((record_id, overlap))

        results.sort(key=lambda pair: pair[1], reverse=True)
        return results

    def __len__(self) -> int:
        return len(self._record_positions)
