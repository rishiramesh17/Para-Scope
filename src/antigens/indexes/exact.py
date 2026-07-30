"""Exact-sequence index (SPEC S16.1): O(1) lookup by sequence checksum."""

from collections import defaultdict


class ExactIndex:
    """Maps a sequence checksum to every record id carrying that exact sequence."""

    def __init__(self):
        self._index: dict[str, list[str]] = defaultdict(list)

    def add(self, checksum: str, record_id: str) -> None:
        if record_id not in self._index[checksum]:
            self._index[checksum].append(record_id)

    def lookup(self, checksum: str) -> list[str]:
        return list(self._index.get(checksum, []))

    def __len__(self) -> int:
        return len(self._index)
