# Adversarial Self-Audit and the Grounding Boundary

Date: 2026-10-04
Artifact: `reports/adversarial_audit.json`
Corpus: schema 6.0.0, 46,727 components · calibration 1.0.0
Versions: `adversarial` 1.0.0 · `grounding` 1.0.0 · `harness` 1.0.0
Tests: 414 passing (was 350)

## What this is

Two halves of one claim, built together because neither is worth much alone.

**Adversarial self-audit** attacks the evidence engine with inputs that have no
right to produce a confident answer, and measures whether confidence actually
degrades. **The grounding boundary** is the gate every answer crosses before a
human reads it: a claim without a citation, or a claim this platform may never
make, does not get out.

The motivation is specific. Gate 1 found a high-severity defect from a *single*
hand-written adversarial query — a random 118-residue sequence scoring
`sequence_identity = 1.0` off five coincidentally matching residues. That was
luck. This makes it systematic.

## Results

```
PASSED 86/86

chimera              3/3     point_mutant        24/24
degenerate           8/8     point_mutant_ladder  6/6
duplicated           6/6     random               3/3
reversed             6/6     shuffled             6/6
truncation          18/18    truncation_ladder    6/6

grounding boundary
  answers rendered           66
  explicit refusals          24
  prohibited claims emitted   0
  citation coverage        1.000
  uncited claims              0
  crashes                     0
```

Mutation ladders are cleanly monotone across all six parents, e.g.

```
d1:0.992 -> d3:0.976 -> d10:0.919 -> d30:0.758
```

More mutations never scored closer. That is the property worth having; a golden
value would break on every legitimate corpus rebuild, while "more distance must
not score closer" holds across corpora and engine versions.

## Two defects found and fixed

### 1. `NOT_FOUND` returned for inputs that are not sequences

Five of eight degenerate inputs — `"1234567890…"`, `"not a protein sequence at
all"`, a gap character, a stop codon, ambiguity codes — returned **`not_found`**
from `compute_therapeutic_precedent`.

`NOT_FOUND` asserts *"the corpus holds no therapeutic precedent for this"*. That
is a claim about the database. For a string that is not a sequence, the true
statement is that it was never a valid query. The platform's whole premise is
that these are different facts, and `germline.py` already validated the correct
way — so two modules disagreed about the same input.

Fixed in `therapeutic_precedent` and `pairing_precedent`: residues outside the
standard 20 now return `UNSUPPORTED_COMPONENT` with the offending characters
named, and no search runs. The input is not repaired to make it searchable —
stripping the character would search a sequence nobody submitted.

### 2. The prohibited-claim scanner blocked the platform's own required disclaimers

The first scanner matched on topic, so it rejected:

> "Window coverage is not an immunogenicity score."

CLAUDE.md **requires** that sentence. A guard that bans the topic rather than the
assertion rejects the exact text that exists to prevent the claim.

Fixed with disclaimer-context detection: a match inside a meta-statement about a
measurement is not an assertion. The distinction is narrow and deliberate —
*"is not an immunogenicity score"* describes the metric and passes; *"has low
immunogenicity"* describes the candidate and does not. `"is not safe"` is still a
safety claim and is deliberately **not** whitelisted. A later disclaimer cannot
launder an earlier claim (`"It will bind. It does not imply anything."` still
fails), which is asserted by test.

## Why the 86/86 is not self-congratulation

A suite that passes on first run is suspicious: a check that cannot fail reports
success on a broken engine exactly as on a correct one. Two things address that,
and one of them is a negative result.

**The detectors are validated by negative controls.** `TestDetectorsActuallyFire`
feeds each check a value it must reject — a 0.99 identity from a random sequence,
a shuffled sequence tying its parent, a mutation ladder that reverses — and
asserts failure is reported. It also asserts each detector *abstains* rather than
passing when it has too few points to judge.

**The audit did not reproduce the original Gate 1 defect, and that is recorded
rather than glossed.** Disabling the coverage floors against the therapeutic
subset produced a best local identity of 0.28, not 1.0. The original defect
needed the full 5,473-sequence VH population for a coincidental short perfect
alignment to exist; the therapeutic subset (1,118 VH) is too small. So the
engine-level reproduction is absent, and detector validity rests on the unit
negative controls instead. An audit whose teeth are only demonstrated at unit
level is weaker than one demonstrated end-to-end, and that gap is real.

## Honest limits of this run

- **The truncation ladder is vacuous.** All six truncation ladders abstained:
  only 6 of 18 truncation cases scored a match at all, because truncated
  sequences are correctly rejected by the coverage floor. "truncation 18/18" means
  the cases were recorded and the family behaved, **not** that a monotone trend
  was verified. The rejection is the right behaviour; the monotonicity claim for
  that family is untested.
- **A single residue still returns `not_found`.** `"E"` passes the alphabet check,
  so it reaches the search and finds nothing. A one-residue "VH" is not a V-domain
  and `UNSUPPORTED_COMPONENT` would be more honest. The candidate fix — reject a
  declared variable domain that ANARCII cannot number — was *not* applied, because
  it would also reject the 0.46% of legitimate chains that fail numbering, and
  that trade was not measured. Recorded rather than guessed at.
- **Six parents, three component classes.** Enough to show the shape; not enough
  to claim coverage of the corpus's diversity.
- **The grounding boundary is not a language model.** The generator is pluggable
  and the guarantee lives in the boundary, not in a prompt: a grounded answer is
  one whose every claim cites a stored evidence record and whose rendered text
  contains no prohibited assertion. That property is checkable offline and
  deterministically, and holds whether the sentence was written by a template, a
  model, or by hand. What is *not* claimed is that any particular model was tested
  — none was.
- **Citation coverage is 1.000 by construction**, since `Claim` refuses to exist
  without a citation. It is measured anyway so a future change that weakens that
  invariant shows up as a number rather than as silence.

## How to re-run

```bash
.venv/bin/python scripts/run_adversarial_audit.py     # writes reports/adversarial_audit.json
.venv/bin/python -m pytest tests/test_verification.py # detectors, generators, boundary
```
