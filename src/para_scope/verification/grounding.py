"""
Grounded claims and the output boundary (Version 2, verification layer).

The problem this solves
-----------------------
Every number this platform produces is traceable: a record id, a source
checksum, an engine version, a calibration version. That traceability is worth
nothing the moment a natural-language layer sits on top and writes a sentence,
because a sentence can assert anything and carries no provenance.

Two failure modes, and they are different:

  **Ungrounded** -- a claim no evidence record supports. "This VH appears in 12
  therapeutics" when nothing was looked up, or when 12 was the wrong field.
  Caught by requiring every claim to name the evidence it rests on.

  **Prohibited** -- a claim that *is* grounded but that this platform may not
  make at all. "This VH appears in 12 therapeutics, so it will likely express
  well." The first half is evidence; the second is a prediction about a molecule
  that no database lookup licenses. Caught by scanning the rendered text.

A prompt instructing a model to avoid both is not a guarantee. This module is a
boundary: `GroundedAnswer.render()` raises rather than returning text that fails
either check. A caller cannot forget to validate, because validation is the only
way to obtain the string.

What this is not
----------------
Not a language model, and not a wrapper around one. The generator is deliberately
pluggable and the guarantee lives in the boundary, not in the generator: a
grounded answer is one whose every claim cites a stored evidence record and whose
rendered text contains no prohibited assertion. That property is checkable
offline, deterministically, and holds whether the sentence was written by a
template, by a model, or by hand. Putting the guarantee in the boundary rather
than the prompt is the entire design.
"""

import re
from dataclasses import dataclass, field
from typing import Optional

# Bump when the prohibited-claim vocabulary or the grounding contract changes. A
# stored answer records which boundary admitted it.
#   1.0.0 -- initial: citation requirement plus prohibited-claim scan
GROUNDING_VERSION = "1.0.0"


# Assertions this platform may never make, from CLAUDE.md's prohibited-claims
# list. Each entry is (pattern, what it would be asserting).
#
# Patterns target *assertions*, not topics. "immunogenicity" must be discussable
# -- the window-coverage module's own disclaimer says it is not a measure of
# immunogenicity, and a scanner that banned the word would reject the disclaimer
# that exists to prevent the claim. So the patterns look for a subject bound to a
# predicate: "will bind", "is safe", "low immunogenicity". A negation or a
# disclaimer phrases it differently and passes.
_PROHIBITED_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"\b(will|would|should|is likely to|expected to)\s+(\w+\s+){0,2}bind\b", "binding prediction"),
    (r"\b(will|would|should|is likely to|expected to)\s+(\w+\s+){0,2}fold\b", "folding prediction"),
    (r"\b(will|would|should|is likely to|expected to)\s+(\w+\s+){0,2}express\b", "expression prediction"),
    (r"\b(will|would|should|is likely to|expected to)\s+(\w+\s+){0,2}pair\b", "pairing prediction"),
    (r"\b(will|would|should|is likely to|expected to)\s+(\w+\s+){0,2}aggregate\b", "aggregation prediction"),
    (r"\b(will|would|should|is likely to|expected to)\s+(\w+\s+){0,2}(succeed|work|fail)\b", "success prediction"),
    (r"\bis\s+(safe|effective|efficacious|manufacturable|soluble|stable|developable)\b", "property assertion"),
    (r"\b(highly|very|more)\s+(safe|effective|efficacious|manufacturable|soluble|stable|developable)\b", "property assertion"),
    (r"\b(low|high|minimal|reduced)\s+immunogenicity\b", "immunogenicity claim"),
    (r"\bimmunogenicity\s+(risk|score|is\s+low|is\s+high)\b", "immunogenicity claim"),
    (r"\b(humanness|naturalness)\s+(score|of|is)\b", "humanness claim"),
    (r"\bprobability\s+of\s+(binding|folding|expression|success|developability)\b", "probability claim"),
    (r"\b(therapeutically|clinically)\s+(viable|promising|successful)\b", "clinical claim"),
    (r"\bdrug[- ]like(ness)?\b", "drug-likeness claim"),
    (r"\bbetter\s+(candidate|antibody|binder)\b", "comparative quality claim"),
    (r"\b(good|bad|poor|excellent|promising)\s+(candidate|antibody|binder)\b", "quality judgement"),
)

_COMPILED = tuple((re.compile(p, re.IGNORECASE), label) for p, label in _PROHIBITED_PATTERNS)

# Phrases that turn a prohibited pattern into a *disclaimer about* that pattern.
#
# This distinction is not decoration. CLAUDE.md requires the window-coverage
# module to state that it is "not an immunogenicity score", and requires every
# index to carry "not a probability of folding, expression, pairing, binding...".
# A scanner without this would reject the platform's own mandatory disclaimers --
# the exact sentences that exist to prevent the claim.
#
# Only *meta* negations suppress: a statement about what a measurement is, not a
# statement about the molecule. "is not an immunogenicity score" describes the
# metric and passes. "has low immunogenicity" describes the candidate and does
# not. "is not safe" is still a safety claim and is deliberately not whitelisted.
_DISCLAIMER_CONTEXT = (
    r"not\s+(a|an)\s*$",
    r"not\s+(a|an)\s+(measure|measurement|prediction|claim|score|index|probability)\s+of\s*$",
    r"(does|do|did)\s+not\s+(imply|mean|establish|predict|indicate)\b[\w\s,'-]{0,24}$",
    r"\bnever\s+(a|an)?\s*$",
    r"\bnot\s+evidence\s+(that|of)\s*$",
    r"\bis\s+not\s*$",
    r"\bare\s+not\s*$",
    r"\bcannot\s+be\s+(read|taken)\s+as\s*$",
)
_COMPILED_DISCLAIMER = tuple(re.compile(p, re.IGNORECASE) for p in _DISCLAIMER_CONTEXT)

# How far back to look for a disclaimer phrase. Long enough to span "Window
# coverage is not an " and short enough that a negation two sentences earlier
# cannot launder a later claim.
_DISCLAIMER_LOOKBEHIND = 40


def _is_disclaimed(text: str, start: int) -> bool:
    """Whether the match at `start` sits inside a disclaimer about the metric."""
    preceding = text[max(0, start - _DISCLAIMER_LOOKBEHIND) : start]
    return any(pattern.search(preceding) for pattern in _COMPILED_DISCLAIMER)


class UngroundedClaimError(ValueError):
    """A claim was asserted with no evidence citation."""


class ProhibitedClaimError(ValueError):
    """Rendered text asserts something this platform may never assert."""


@dataclass(frozen=True)
class Citation:
    """Where a claim's evidence came from, in enough detail to re-derive it.

    These are the fields that make a claim auditable rather than plausible. A
    citation naming only a record id would not let a reader reproduce the number:
    the same record under a different engine or calibration version can yield a
    different value, and that is the whole reason both are versioned.
    """

    record_id: str
    source_database: str
    source_checksum: str
    engine_version: str
    # Present only for a calibrated (Layer 2) claim. None on a raw measurement,
    # and that distinction is itself information: an uncalibrated number is not a
    # percentile and must not read as one.
    calibration_version: Optional[str] = None

    def short(self) -> str:
        base = f"{self.record_id} ({self.source_database}@{self.source_checksum[:8]}, {self.engine_version}"
        if self.calibration_version:
            base += f", calib {self.calibration_version}"
        return base + ")"


@dataclass(frozen=True)
class Claim:
    """One assertion plus the evidence it rests on.

    `text` is what a reader sees. `citations` is why they should believe it. A
    Claim with no citations cannot be constructed -- the check is in __post_init__
    rather than in a validator a caller might not call.
    """

    text: str
    citations: tuple[Citation, ...]

    def __post_init__(self):
        if not self.text.strip():
            raise UngroundedClaimError("A claim must say something")
        if not self.citations:
            raise UngroundedClaimError(
                f"Claim has no citation and cannot be made: {self.text!r}. Every "
                f"assertion must name the evidence record supporting it; an "
                f"unsupported sentence is exactly what this boundary exists to stop."
            )


@dataclass
class Refusal:
    """An explicit statement that a question cannot be answered from evidence.

    A refusal is a result, not a failure. Returning nothing, or returning a
    hedged guess, both read as "no evidence found" -- which is a different and
    much stronger statement than "this platform cannot answer that".
    """

    question: str
    reason: str
    # One of: out_of_scope (the platform may not answer it at all),
    # no_evidence (answerable in principle, nothing found),
    # insufficient_evidence (found, but too little to report).
    category: str


def scan_for_prohibited_claims(text: str) -> list[tuple[str, str]]:
    """Every prohibited assertion in `text`, as (matched span, what it asserts).

    A match inside a disclaimer about the measurement is not an assertion and is
    skipped -- see _DISCLAIMER_CONTEXT for why that exemption has to exist.
    """
    found: list[tuple[str, str]] = []
    for pattern, label in _COMPILED:
        for match in pattern.finditer(text):
            if _is_disclaimed(text, match.start()):
                continue
            found.append((match.group(0), label))
    return found


@dataclass
class GroundedAnswer:
    """A set of cited claims and explicit refusals, renderable only if it passes.

    `render()` is the boundary. It raises on a prohibited assertion rather than
    returning text with a warning attached, because a warning beside a sentence
    does not stop the sentence being read, quoted, or pasted into a slide.
    """

    question: str
    claims: list[Claim] = field(default_factory=list)
    refusals: list[Refusal] = field(default_factory=list)
    grounding_version: str = GROUNDING_VERSION

    def citation_coverage(self) -> float:
        """Fraction of claims carrying at least one citation.

        Always 1.0 by construction -- Claim refuses to exist without one. Reported
        anyway so an audit measures the property rather than trusting the
        invariant, and so a future change that weakens Claim shows up as a number
        rather than as silence.
        """
        if not self.claims:
            return 1.0
        return sum(1 for c in self.claims if c.citations) / len(self.claims)

    def render(self) -> str:
        """The answer as text, or an exception explaining why it may not be shown."""
        if not self.claims and not self.refusals:
            raise UngroundedClaimError(
                "An answer must contain at least one cited claim or one explicit "
                "refusal. Returning empty text would read as 'nothing found', "
                "which is a claim in itself."
            )

        lines: list[str] = []
        for claim in self.claims:
            lines.append(claim.text)
            for citation in claim.citations:
                lines.append(f"    evidence: {citation.short()}")
        for refusal in self.refusals:
            lines.append(f"Cannot answer: {refusal.reason} [{refusal.category}]")

        text = "\n".join(lines)
        violations = scan_for_prohibited_claims(text)
        if violations:
            rendered = "; ".join(f"{span!r} ({label})" for span, label in violations)
            raise ProhibitedClaimError(
                f"Refusing to emit an answer asserting: {rendered}. This platform "
                f"reports comparative database evidence and may not predict "
                f"binding, folding, expression, stability, safety, efficacy, "
                f"immunogenicity, or experimental success."
            )
        return text
