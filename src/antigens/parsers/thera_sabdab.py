"""
Parser for Thera-SAbDab's therapeutic sequence/structure download (SPEC S13-14, Gate 2).

Source: https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv
Confirmed by direct inspection of the real downloaded file (2026-07-29,
last-modified 2025-03-18): variable-domain sequences only (VH/VL), plus
therapeutic metadata (isotype, clinical stage, target, structural links).
There is no full heavy/light chain (constant region) column in this file --
FULL_HEAVY / FULL_LIGHT remain unpopulated by this source, same as SAbDab2.

The 'Format' column is free text, not a controlled vocabulary: the real file
has 120+ distinct values, mostly bespoke bispecific/multispecific
descriptions (e.g. "Bispecific Mixed mAb and scFv-CH2-CH3", "Trispecific
Single Domains (VH-VH-VH'-VH''-VH'); Bispecific entry"). Building an
exhaustive classifier over this text is not attempted. Instead, every row is
parsed into up to two independent variable-domain arms:

  - primary arm:     HeavySequence / LightSequence
  - bispecific arm:  HeavySequence(ifbispec) / LightSequence(ifbispec)

Each arm is parsed on its own terms (present/absent heavy, present/absent
light) and tagged with an arm index. Neither arm is inferred from the other,
and a bispecific's second arm is never paired with the primary arm's light
chain (CLAUDE.md: paired evidence requires both chains from the same source
record -- a bispecific's two arms are two different specificities, not one
pair with a spare light chain).

Single-domain (VHH-like) heavy chains are detected via a keyword heuristic
against the Format text (never against sequence content), since ANARCII's
own chain-type output only distinguishes H from K/L, not VH from VHH. This
heuristic is necessarily approximate given free-text input; rows it cannot
confidently classify default to conventional VH, which fails safe. If the
Format text is later found to false-negative to VHH in a specific case,
that's a heuristic-accuracy correction, not a schema issue.

Two missing-data sentinels, not one
-----------------------------------
The sequence and metadata columns write missing data as the literal string
'na'. The three structure-link columns do not: they write the literal string
'None'. Treating only 'na' as missing made every structure column truthy, so
794 of 1,131 therapeutics were recorded as having an experimental structure
record when the source says they have none -- a fabricated evidence flag on
a field SPEC S26 scores directly. Structure columns therefore get their own
parser rather than sharing the 'na' one; the general sentinel is left alone so
a legitimate value elsewhere cannot be clobbered by widening it.

A third shape exists and means something different again: a blank cell. It is
not a sentinel -- it asserts nothing -- so it resolves to *unknown*, never to
"no structure". The metadata columns already treated blank as missing; the
structure columns did not, and read it as a positive statement of absence.

Structure links are per arm
---------------------------
Within a structure column, a semicolon separates the primary arm's links from
the bispecific arm's ('None;8z8v:B' = nothing for arm 1, one structure for
arm 2), and a slash separates multiple structures within one arm. Measured on
the real file: all 157 two-arm rows are semicolon-split and none of the 974
single-arm rows are, across all three columns (3,393/3,393 values align).
OR-ing the columns together and giving both arms the same flag would credit
arm 1 with arm 2's structures. If a future release breaks that alignment, the
affected column resolves to *unknown* for that row rather than being
misattributed -- a wrong arm is worse than an absent one.
"""

import csv
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

#   1.0.0 -- initial
#   1.1.0 -- structure columns use the 'None' sentinel, not 'na', and are split
#            per arm. Before this, structure_available was True for every
#            therapeutic (794 of them falsely).
#   1.2.0 -- a blank structure cell resolves to unknown, not to "no structure".
#            Previously indistinguishable from the 'None' sentinel.
#   1.3.0 -- structure tokens must be shaped like PDB references. Catches an
#            unrecognised sentinel by the shape of the output rather than by
#            enumerating known bad inputs, which is what the previous three
#            sentinel fixes each had to do after the fact.
PARSER_VERSION = "1.3.0"

# Column names copied verbatim from the real downloaded header, including the
# embedded date-stamp in "Highest_Clin_Trial (Feb '25)". If a future release
# renames this column (e.g. a new quarter's date-stamp), this parser must
# fail loudly rather than silently drop the clinical-stage field -- schema
# drift is a defect to surface, not paper over (CLAUDE.md).
REQUIRED_COLUMNS = frozenset({
    "Therapeutic", "Format", "CH1 Isotype", "VD LC",
    "Highest_Clin_Trial (Feb '25)", "Est. Status",
    "HeavySequence", "LightSequence",
    "HeavySequence(ifbispec)", "LightSequence(ifbispec)",
    "100% SI Structure", "99% SI Structure", "95-98% SI Structure",
    "Target", "Alternative Therapeutic Names",
})

_SINGLE_DOMAIN_KEYWORDS = ("nanobody", "vhh", "single domain")

_NA_SENTINEL = "na"

# The structure columns' own missing-data sentinel. Distinct from _NA_SENTINEL on
# purpose -- see the module docstring. Deliberately not folded into _parse_na:
# widening the general sentinel would silently drop any other column whose real
# value happened to be the word "None".
_STRUCTURE_ABSENT_SENTINEL = "none"

# Separators inside a structure column: arms first, then structures within an arm.
_ARM_SEPARATOR = ";"
_STRUCTURE_SEPARATOR = "/"

# A PDB entry id (4 characters, first is a digit) with optional :chain suffixes:
# '1n8z:BA', '6bi0:HL:IM', '7t97:HL:AB:CD:EF'.
#
# This validates the *shape of the output* rather than enumerating known bad
# inputs, and that difference is the point. Every missing-data sentinel found in
# this file so far -- 'na', then 'None', then a blank -- was handled reactively,
# after it had already produced or could produce a wrong answer. Each fix knew
# only about the shape someone had already met.
#
# A structure reference that is not shaped like a structure reference is not a
# structure reference, whatever new spelling of "missing" a future release
# invents. 'None' would have failed this check on day one, and so would 'N/A',
# '-', 'TBD' or 'unknown' without anyone having to predict them.
#
# Measured on the 2025-03-18 file: 3,399 values, 679 distinct tokens, and exactly
# two that do not match -- 'None' (the sentinel, 3,428 occurrences) and one empty
# token produced by a double slash in row 1083 (Vislarafusp). So this pattern
# describes the real vocabulary rather than an assumption about it.
_PDB_REFERENCE = re.compile(r"^[0-9][A-Za-z0-9]{3}(:[A-Za-z0-9]+)*$")


def is_structure_reference(token: str) -> bool:
    """Whether a token is shaped like a PDB entry reference."""
    return bool(_PDB_REFERENCE.match(token.strip()))


def _parse_na(value: Optional[str]) -> Optional[str]:
    """Thera-SAbDab uses the literal string 'na' (not a blank cell) for missing data."""
    if value is None:
        return None
    cleaned = value.strip()
    if not cleaned or cleaned.lower() == _NA_SENTINEL:
        return None
    return cleaned


def split_structure_column(
    value: Optional[str],
    arm_count: int,
    anomalies: Optional[list[str]] = None,
) -> Optional[list[list[str]]]:
    """
    Split one structure column into per-arm lists of structure references.

    Returns one list per arm, each holding that arm's structure identifiers with
    the 'None' sentinel removed (so an empty list means "the source states there
    is no structure", which is a finding, not missing data).

    Returns None -- meaning *unknown*, not *absent* -- in three cases, all of them
    schema drift, and under drift the honest answer is that the column cannot be
    read:

    - the number of semicolon-separated tokens does not match the number of arms
      (never happens on the file this parser was written against: 3,393/3,393
      values align). Guessing would attribute one arm's experimental structures to
      a different specificity.
    - any arm token is empty. A blank states nothing; the 'None' sentinel states
      that no structure exists. Collapsing them would score 0.0 off an empty cell.
    - any surviving token is not shaped like a PDB reference. Anything that
      outlived sentinel removal is claiming to *be* a structure record, so it has
      to look like one -- see _PDB_REFERENCE. Offending tokens are appended to
      `anomalies` when a list is supplied, so the caller can report drift instead
      of silently losing the only signal that it happened.
    """
    if value is None or arm_count <= 0:
        return None

    tokens = [token.strip() for token in value.split(_ARM_SEPARATOR)]
    if len(tokens) != arm_count:
        return None

    # An empty token is not the 'None' sentinel and must not be read as one. A
    # blank cell states nothing; the sentinel states "no structure exists". Before
    # this check they were indistinguishable -- a blank column yielded [[]], the
    # same value an explicit 'None' yields, so has_structure_record answered False
    # and SPEC S26's structure term scored 0.0 off a cell with no content. That is
    # the 794-record defect inverted: there a missing sentinel was read as evidence
    # present, here an absent value is read as evidence of absence. Unknown is the
    # only honest answer, and it keeps the term out of the score entirely.
    if any(not token for token in tokens):
        return None

    per_arm: list[list[str]] = []
    for token in tokens:
        refs = [
            ref.strip()
            for ref in token.split(_STRUCTURE_SEPARATOR)
            if ref.strip() and ref.strip().lower() != _STRUCTURE_ABSENT_SENTINEL
        ]
        # Anything surviving sentinel removal claims to be a structure record, so
        # it must look like one. An unrecognised token is either a sentinel this
        # parser has never seen or genuine corruption; both mean the column cannot
        # be trusted, and neither may become a structure id that SPEC S26 scores.
        unrecognised = [ref for ref in refs if not is_structure_reference(ref)]
        if unrecognised:
            if anomalies is not None:
                anomalies.extend(unrecognised)
            return None
        per_arm.append(refs)
    return per_arm


def is_single_domain_format(format_text: Optional[str]) -> bool:
    """Heuristic: does the free-text Format column indicate a single-domain (VHH-like) heavy chain?"""
    if not format_text:
        return False
    lowered = format_text.lower()
    return any(keyword in lowered for keyword in _SINGLE_DOMAIN_KEYWORDS)


@dataclass
class TheraSAbDabArm:
    """One independently-parseable variable-domain arm (primary, or a bispecific's second arm)."""
    arm_index: int  # 1 = primary, 2 = bispecific second arm
    heavy_sequence: Optional[str]
    light_sequence: Optional[str]
    is_single_domain: bool  # heuristic: heavy chain has no paired VL by design (VHH-like)

    # Structure references for *this arm only*, at each sequence-identity tier.
    # Three-state, and the distinction is load-bearing:
    #   None  -- unknown; the column could not be mapped to arms (schema drift)
    #   []    -- the source states this arm has no structure at this tier
    #   [ids] -- these structure records
    structures_100_si: Optional[list[str]] = None
    structures_99_si: Optional[list[str]] = None
    structures_95_98_si: Optional[list[str]] = None

    @property
    def has_structure_record(self) -> Optional[bool]:
        """True/False if any tier is known, None if every tier is unknown.

        None must never be coerced to False by a caller: "we could not determine
        whether a structure exists" is not "no structure exists" (CLAUDE.md,
        missing evidence is not zero evidence).
        """
        tiers = (self.structures_100_si, self.structures_99_si, self.structures_95_98_si)
        if all(tier is None for tier in tiers):
            return None
        return any(tier for tier in tiers)


@dataclass
class TheraSAbDabRow:
    """One parsed therapeutic entry, verbatim except type coercion and the 'na' sentinel."""
    therapeutic_name: str            # Therapeutic -- also the parent-therapeutic dedup key
    format_raw: str                  # Format -- free text, kept verbatim, never classified beyond is_single_domain
    ch1_isotype: Optional[str]
    light_chain_class_label: Optional[str]  # VD LC -- metadata only; kappa/lambda classification still comes from ANARCII
    highest_clinical_trial: Optional[str]   # metadata, never a success/failure label
    development_status: Optional[str]       # Est. Status -- metadata, never a success/failure label
    target: Optional[str]
    alternative_names: Optional[str]
    structure_100_si: Optional[str]
    structure_99_si: Optional[str]
    structure_95_98_si: Optional[str]
    arms: list[TheraSAbDabArm]
    supported: bool                  # False only if no heavy sequence exists in any arm
    skip_reason: Optional[str] = None
    # Structure tokens that survived sentinel removal but are not shaped like PDB
    # references. Non-empty means the source used a spelling of "missing" this
    # parser does not know, or the column is corrupt; the affected columns already
    # resolved to unknown, and this carries the offending tokens up so ingestion
    # can report them rather than discarding the only evidence of drift.
    structure_anomalies: list[str] = field(default_factory=list)


def parse_thera_sabdab_csv(path: Path) -> Iterator[TheraSAbDabRow]:
    """
    Stream-parse TheraSAbDab_SeqStruc_OnlineDownload.csv, yielding one TheraSAbDabRow per therapeutic.

    Raises ValueError immediately if the header is missing any column this parser
    depends on -- schema drift must fail loudly, never be guessed around.
    """
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"TheraSAbDab CSV is missing required columns: {sorted(missing)}")

        for row in reader:
            format_raw = row["Format"].strip()
            single_domain = is_single_domain_format(format_raw)

            primary_heavy = _parse_na(row["HeavySequence"])
            primary_light = _parse_na(row["LightSequence"])
            bispec_heavy = _parse_na(row["HeavySequence(ifbispec)"])
            bispec_light = _parse_na(row["LightSequence(ifbispec)"])

            arms: list[TheraSAbDabArm] = []
            if primary_heavy is not None:
                arms.append(TheraSAbDabArm(
                    arm_index=1,
                    heavy_sequence=primary_heavy,
                    light_sequence=primary_light,
                    is_single_domain=single_domain,
                ))
            if bispec_heavy is not None:
                arms.append(TheraSAbDabArm(
                    arm_index=2,
                    heavy_sequence=bispec_heavy,
                    light_sequence=bispec_light,
                    is_single_domain=single_domain,
                ))

            # Structure links resolve only once the arm count is known, since the
            # columns are split by arm. Assigned by arm_index, not list position.
            raw_structures = {
                "structures_100_si": row["100% SI Structure"],
                "structures_99_si": row["99% SI Structure"],
                "structures_95_98_si": row["95-98% SI Structure"],
            }
            structure_anomalies: list[str] = []
            for field_name, raw_value in raw_structures.items():
                per_arm = split_structure_column(raw_value, len(arms), structure_anomalies)
                for arm in arms:
                    resolved = per_arm[arm.arm_index - 1] if per_arm is not None else None
                    setattr(arm, field_name, resolved)

            supported = len(arms) > 0
            skip_reason = None if supported else (
                "unsupported_component: no heavy-chain sequence available for this therapeutic entry"
            )

            yield TheraSAbDabRow(
                therapeutic_name=row["Therapeutic"].strip(),
                format_raw=format_raw,
                ch1_isotype=_parse_na(row["CH1 Isotype"]),
                light_chain_class_label=_parse_na(row["VD LC"]),
                highest_clinical_trial=_parse_na(row["Highest_Clin_Trial (Feb '25)"]),
                development_status=_parse_na(row["Est. Status"]),
                target=_parse_na(row["Target"]),
                alternative_names=_parse_na(row["Alternative Therapeutic Names"]),
                structure_100_si=_parse_na(row["100% SI Structure"]),
                structure_99_si=_parse_na(row["99% SI Structure"]),
                structure_95_98_si=_parse_na(row["95-98% SI Structure"]),
                arms=arms,
                supported=supported,
                skip_reason=skip_reason,
                structure_anomalies=structure_anomalies,
            )


def compute_file_checksum(path: Path) -> str:
    """SHA256 of the raw source file, for SourceSnapshot provenance."""
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha256.update(chunk)
    return sha256.hexdigest()
