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
"""

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

PARSER_VERSION = "1.0.0"

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


def _parse_na(value: Optional[str]) -> Optional[str]:
    """Thera-SAbDab uses the literal string 'na' (not a blank cell) for missing data."""
    if value is None:
        return None
    cleaned = value.strip()
    if not cleaned or cleaned.lower() == _NA_SENTINEL:
        return None
    return cleaned


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
            )


def compute_file_checksum(path: Path) -> str:
    """SHA256 of the raw source file, for SourceSnapshot provenance."""
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha256.update(chunk)
    return sha256.hexdigest()
