"""Search indexes (SPEC S16). Each index type enforces its own eligibility rules."""

from para_scope.indexes.exact import ExactIndex
from para_scope.indexes.approximate import ApproximateIndex, AlignmentHit
from para_scope.indexes.positional import PositionalIndex
from para_scope.indexes.pair import PairIndex, PairRecord

__all__ = [
    "ExactIndex",
    "ApproximateIndex",
    "AlignmentHit",
    "PositionalIndex",
    "PairIndex",
    "PairRecord",
]
