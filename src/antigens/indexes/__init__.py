"""Search indexes (SPEC S16). Each index type enforces its own eligibility rules."""

from antigens.indexes.exact import ExactIndex
from antigens.indexes.approximate import ApproximateIndex, AlignmentHit
from antigens.indexes.positional import PositionalIndex
from antigens.indexes.pair import PairIndex, PairRecord

__all__ = [
    "ExactIndex",
    "ApproximateIndex",
    "AlignmentHit",
    "PositionalIndex",
    "PairIndex",
    "PairRecord",
]
