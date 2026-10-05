"""Verification layer: adversarial self-audit and grounded answering.

Two halves of one claim. `adversarial` attacks the evidence engine with inputs
that have no right to produce confident answers and measures whether confidence
actually degrades. `grounding` is the boundary every answer must cross before a
human sees it: a claim without a citation, or a claim the platform is forbidden
to make, does not get out.
"""

from para_scope.verification.adversarial import (
    ADVERSARIAL_VERSION,
    AdversarialCase,
    build_case_suite,
)
from para_scope.verification.grounding import (
    GROUNDING_VERSION,
    Claim,
    GroundedAnswer,
    ProhibitedClaimError,
    UngroundedClaimError,
    scan_for_prohibited_claims,
)

__all__ = [
    "ADVERSARIAL_VERSION",
    "AdversarialCase",
    "build_case_suite",
    "GROUNDING_VERSION",
    "Claim",
    "GroundedAnswer",
    "ProhibitedClaimError",
    "UngroundedClaimError",
    "scan_for_prohibited_claims",
]
