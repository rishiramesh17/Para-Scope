"""
Adversarial case generation for self-audit (Version 2, verification layer).

Why this exists
---------------
Gate 1 found a high-severity defect from a single adversarial query: a randomly
generated 118-residue sequence returned a top hit at `sequence_identity = 1.0`
because a local aligner found five coincidentally matching residues. The measure
was real, the alignment was real, and the answer was nonsense.

That was luck -- one hand-written query in a 13-case benchmark. This module makes
it systematic. It generates inputs that have no right to produce a confident
answer, and the harness that consumes it asserts that confidence actually
degrades. Finding nothing is a result; finding something is the point.

What an adversarial case is
---------------------------
Not a random fuzz string. Each family isolates one way a sequence can be wrong
while holding others fixed, so a failure localises:

  random          no relationship to any antibody. Nothing may match it well.
  shuffled        identical amino-acid composition, destroyed order. Separates
                  "similar because of composition" from "similar because of
                  sequence" -- a composition-driven score cannot tell them apart.
  reversed        identical composition, reversed order. A reversed domain is not
                  a domain; any method treating it as one is reading composition.
  chimera         two real halves from unrelated antibodies. Each half matches
                  something real, the whole matches nothing. Coverage rules decide
                  whether this is reported honestly.
  point_mutant    k substitutions from a real parent, k increasing. The only
                  family with a *monotone* expectation: more mutations must not
                  score closer.
  truncation      a real prefix. Tests whether partial coverage is penalised
                  rather than scored on the aligned part alone.
  duplicated      a real domain concatenated with itself. Regression guard for the
                  variable/constant split (schema 6.0.0): a ~240-residue input
                  must not be compared whole against ~120-residue domains.
  degenerate      empty, whitespace, non-residue characters, ambiguity codes.
                  Must produce an explicit status, never a number.

Determinism
-----------
Every generator is seeded from a fixed constant and derives its randomness from
the case id, so the suite is identical across runs and machines. An adversarial
result that cannot be reproduced cannot be investigated.
"""

import hashlib
import random
from dataclasses import dataclass, field
from typing import Iterable, Optional

from para_scope.models import ComponentType

# Bump when a family is added, removed, or its generation changes. Stored on the
# suite so a recorded audit names the cases it actually ran.
#   1.0.0 -- initial eight families
ADVERSARIAL_VERSION = "1.0.0"

SUITE_SEED = 20261004

STANDARD_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"

# Expectations a case asserts. Declarative so the harness, not the generator,
# owns how each is checked -- a generator that also judged its own output could
# not be wrong in a way the harness would notice.
EXPECT_NO_CONFIDENT_MATCH = "no_confident_match"
EXPECT_WORSE_THAN_PARENT = "worse_than_parent"
EXPECT_MONOTONE_IN_DISTANCE = "monotone_in_distance"
EXPECT_EXPLICIT_STATUS_NOT_A_NUMBER = "explicit_status_not_a_number"
EXPECT_NORMALIZED_TO_DOMAIN = "normalized_to_domain"


@dataclass(frozen=True)
class AdversarialCase:
    """One adversarial input and what a correct system must do with it."""

    case_id: str
    family: str
    description: str
    sequence: str
    component_type: ComponentType
    expectation: str

    # Set when the case is derived from a real corpus component, so a failure can
    # be compared against the parent's own result rather than judged in isolation.
    parent_record_id: Optional[str] = None
    parent_sequence: Optional[str] = None
    # Distance from the parent, for families where that is meaningful (mutation
    # count, residues removed). Used to check monotonicity.
    distance: Optional[int] = None
    notes: list[str] = field(default_factory=list)


def _rng(case_id: str) -> random.Random:
    """A generator seeded from the suite seed and the case id."""
    digest = hashlib.sha256(f"{SUITE_SEED}:{case_id}".encode()).hexdigest()
    return random.Random(int(digest[:16], 16))


def random_sequence(case_id: str, length: int) -> str:
    rng = _rng(case_id)
    return "".join(rng.choice(STANDARD_RESIDUES) for _ in range(length))


def shuffled_sequence(case_id: str, parent: str) -> str:
    rng = _rng(case_id)
    residues = list(parent)
    rng.shuffle(residues)
    return "".join(residues)


def reversed_sequence(parent: str) -> str:
    return parent[::-1]


def chimeric_sequence(left: str, right: str) -> str:
    """First half of one real sequence, second half of another."""
    return left[: len(left) // 2] + right[len(right) // 2 :]


def point_mutant(case_id: str, parent: str, substitutions: int) -> str:
    """`substitutions` residues replaced with a different residue, positions fixed by seed."""
    rng = _rng(case_id)
    residues = list(parent)
    positions = rng.sample(range(len(residues)), min(substitutions, len(residues)))
    for position in positions:
        original = residues[position]
        choices = [r for r in STANDARD_RESIDUES if r != original]
        residues[position] = rng.choice(choices)
    return "".join(residues)


def truncate(parent: str, keep_fraction: float) -> str:
    return parent[: max(1, int(len(parent) * keep_fraction))]


DEGENERATE_INPUTS: tuple[tuple[str, str], ...] = (
    ("empty", ""),
    ("whitespace_only", "   \t  "),
    ("digits", "1234567890123456789012345678901234567890"),
    ("ambiguity_codes", "EVQLXESGGGLVQPGGSLRLSCAASBFNIKDTYIHZVRQ"),
    ("lowercase_junk", "not a protein sequence at all"),
    ("single_residue", "E"),
    ("stop_codon_character", "EVQLVESGGG*LVQPGGSLRLSCAAS"),
    ("gap_characters", "EVQLVESGGG-LVQPGGSLRLSCAAS"),
)


def build_case_suite(
    parents: Iterable[tuple[str, str, ComponentType]],
    mutation_ladder: tuple[int, ...] = (1, 3, 10, 30),
    truncation_ladder: tuple[float, ...] = (0.75, 0.5, 0.25),
    random_lengths: tuple[int, ...] = (60, 120, 240),
) -> list[AdversarialCase]:
    """
    Build the full suite from real parent components.

    `parents` is (record_id, sequence, component_type). Real parents matter: a
    suite built only from synthetic sequences tests the aligner's behaviour on
    noise, not its behaviour on things that *look* like antibodies, which is where
    a similarity engine actually fails.
    """
    parents = list(parents)
    cases: list[AdversarialCase] = []

    for length in random_lengths:
        case_id = f"random_{length}"
        cases.append(
            AdversarialCase(
                case_id=case_id,
                family="random",
                description=f"Uniformly random {length}-residue sequence",
                sequence=random_sequence(case_id, length),
                component_type=ComponentType.VH,
                expectation=EXPECT_NO_CONFIDENT_MATCH,
                notes=[
                    "Has no evolutionary relationship to any antibody. Any "
                    "confident match is an artifact of the scoring method."
                ],
            )
        )

    for name, sequence in DEGENERATE_INPUTS:
        cases.append(
            AdversarialCase(
                case_id=f"degenerate_{name}",
                family="degenerate",
                description=f"Degenerate input: {name}",
                sequence=sequence,
                component_type=ComponentType.VH,
                expectation=EXPECT_EXPLICIT_STATUS_NOT_A_NUMBER,
                notes=["Must yield a status explaining why, never a measurement."],
            )
        )

    for index, (record_id, sequence, component_type) in enumerate(parents):
        tag = f"p{index}"

        case_id = f"shuffled_{tag}"
        cases.append(
            AdversarialCase(
                case_id=case_id, family="shuffled",
                description="Parent's exact amino-acid composition, order destroyed",
                sequence=shuffled_sequence(case_id, sequence),
                component_type=component_type,
                expectation=EXPECT_WORSE_THAN_PARENT,
                parent_record_id=record_id, parent_sequence=sequence,
                notes=[
                    "Identical composition to the parent. A score that cannot "
                    "separate this from the parent is measuring composition."
                ],
            )
        )

        cases.append(
            AdversarialCase(
                case_id=f"reversed_{tag}", family="reversed",
                description="Parent reversed end to end",
                sequence=reversed_sequence(sequence),
                component_type=component_type,
                expectation=EXPECT_WORSE_THAN_PARENT,
                parent_record_id=record_id, parent_sequence=sequence,
            )
        )

        cases.append(
            AdversarialCase(
                case_id=f"duplicated_{tag}", family="duplicated",
                description="Parent concatenated with itself",
                sequence=sequence + sequence,
                component_type=component_type,
                expectation=EXPECT_NORMALIZED_TO_DOMAIN,
                parent_record_id=record_id, parent_sequence=sequence,
                notes=[
                    "Regression guard for the variable/constant split: a "
                    "double-length input must be reduced to one domain before "
                    "comparison, not aligned whole against single domains."
                ],
            )
        )

        for substitutions in mutation_ladder:
            case_id = f"mutant{substitutions}_{tag}"
            cases.append(
                AdversarialCase(
                    case_id=case_id, family="point_mutant",
                    description=f"Parent with {substitutions} substitution(s)",
                    sequence=point_mutant(case_id, sequence, substitutions),
                    component_type=component_type,
                    expectation=EXPECT_MONOTONE_IN_DISTANCE,
                    parent_record_id=record_id, parent_sequence=sequence,
                    distance=substitutions,
                )
            )

        for fraction in truncation_ladder:
            cases.append(
                AdversarialCase(
                    case_id=f"trunc{int(fraction * 100)}_{tag}", family="truncation",
                    description=f"First {int(fraction * 100)}% of the parent",
                    sequence=truncate(sequence, fraction),
                    component_type=component_type,
                    expectation=EXPECT_MONOTONE_IN_DISTANCE,
                    parent_record_id=record_id, parent_sequence=sequence,
                    distance=len(sequence) - len(truncate(sequence, fraction)),
                )
            )

    # Chimeras need two unrelated parents of the same component type.
    for index in range(len(parents) - 1):
        left_id, left_seq, left_type = parents[index]
        right_id, right_seq, right_type = parents[index + 1]
        if left_type is not right_type:
            continue
        cases.append(
            AdversarialCase(
                case_id=f"chimera_p{index}_p{index + 1}", family="chimera",
                description="First half of one real domain, second half of another",
                sequence=chimeric_sequence(left_seq, right_seq),
                component_type=left_type,
                expectation=EXPECT_NO_CONFIDENT_MATCH,
                parent_record_id=left_id, parent_sequence=left_seq,
                notes=[
                    "Each half matches something real; the whole matches nothing. "
                    f"Right half from {right_id}."
                ],
            )
        )

    return cases
