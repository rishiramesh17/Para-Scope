"""
Parser for SAbDab2 ML-grade split metadata (SPEC S13-14).

Source: https://zenodo.org/records/20083995 (splits.tar.gz -> splits_final/ab_split.csv)
Schema documented in splits_final/README.md (Capel et al. 2026, SAbDab2).

This parser only reads the split CSV (sequences + structure metadata). It never
opens the .cif coordinate files bundled in the same archive -- coordinate parsing
is out of scope for this first version (CLAUDE.md structural boundary). Structure
*metadata* (PDB ID, method, resolution, antigen-complex status) is in scope and is
carried through on every row.
"""

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

PARSER_VERSION = "1.0.0"

# Antibody types ('type' column) this parser can place into the ComponentType
# ontology (models.py). VNAR (shark IgNAR) has no corresponding component type --
# rows of this type are read but flagged unsupported rather than silently folded
# into VHH (CLAUDE.md: never combine without explicit remapping).
SUPPORTED_TYPES = frozenset({"FV", "FAB", "FAB+FC", "SD-H", "SD-L"})
UNSUPPORTED_TYPES = frozenset({"VNAR"})

REQUIRED_COLUMNS = frozenset({
    "INSTANCE", "PDB_ID", "SABDAB_ID", "HEAVY_ID", "LIGHT_ID", "method",
    "resolution", "type", "construct", "holo", "Hchain", "Lchain",
    "Hseq_expected", "Lseq_expected", "ab_split",
})


@dataclass
class SAbDab2Row:
    """One parsed structure-instance row from ab_split.csv, verbatim except type coercion."""
    instance_id: str                     # INSTANCE -- unique structure-instance key
    pdb_id: str                           # PDB_ID
    sabdab_id: str                        # SABDAB_ID -- parent-antibody dedup key (SPEC S15)
    heavy_numbering_id: Optional[str]     # HEAVY_ID
    light_numbering_id: Optional[str]     # LIGHT_ID
    method: str                           # 'XRAY' | 'EM'
    resolution: Optional[float]
    antibody_type: str                    # 'FV' | 'FAB' | 'FAB+FC' | 'VNAR' | 'SD-H' | 'SD-L'
    construct: str                        # '' | 'SCFV' | 'DIABODY' | 'CODV' | 'SCDB' | 'DVD' | 'OTHER'
    holo: bool                            # antigen-bound structure?
    heavy_chain_author_id: Optional[str]  # Hchain
    light_chain_author_id: Optional[str]  # Lchain
    heavy_sequence: Optional[str]         # Hseq_expected: complete annotated sequence (incl. unresolved)
    light_sequence: Optional[str]         # Lseq_expected
    train_test_split: str                 # ab_split: 'train' | 'test' (carried for future calibration use)
    supported: bool                       # False for VNAR -- caller must not ingest as VHH/VL
    skip_reason: Optional[str] = None


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _parse_optional_float(value: str) -> Optional[float]:
    value = value.strip()
    return float(value) if value else None


def _parse_optional_str(value: str) -> Optional[str]:
    value = value.strip()
    return value if value else None


def parse_ab_split_csv(path: Path) -> Iterator[SAbDab2Row]:
    """
    Stream-parse ab_split.csv, yielding one SAbDab2Row per structure instance.

    Raises ValueError immediately if the header is missing any column this parser
    depends on -- a schema drift must fail loudly, never be guessed around.
    """
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"ab_split.csv is missing required columns: {sorted(missing)}")

        for row in reader:
            antibody_type = row["type"].strip()
            heavy_sequence = _parse_optional_str(row["Hseq_expected"])
            light_sequence = _parse_optional_str(row["Lseq_expected"])

            # Engineered single-chain constructs (scFv, diabody, single-chain diabody)
            # report the *same* fused VH-linker-VL polypeptide in both Hseq_expected
            # and Lseq_expected -- the source does not independently split the domains.
            # Numbering that fused sequence as if it were one clean V-domain would
            # silently misrepresent it as a plain VH or VL component (CLAUDE.md: never
            # silently modify input / treat a fusion as an isolated chain), so these
            # rows are marked unsupported rather than guessed at.
            is_fused_single_chain = (
                heavy_sequence is not None
                and light_sequence is not None
                and heavy_sequence == light_sequence
            )

            supported = antibody_type in SUPPORTED_TYPES and not is_fused_single_chain
            if antibody_type not in SUPPORTED_TYPES:
                skip_reason = f"unsupported_component: antibody type '{antibody_type}' has no ComponentType mapping"
            elif is_fused_single_chain:
                skip_reason = (
                    "unsupported_component: engineered single-chain construct "
                    "(Hseq_expected == Lseq_expected; heavy/light not independently split in source data)"
                )
            else:
                skip_reason = None

            yield SAbDab2Row(
                instance_id=row["INSTANCE"].strip(),
                pdb_id=row["PDB_ID"].strip(),
                sabdab_id=row["SABDAB_ID"].strip(),
                heavy_numbering_id=_parse_optional_str(row["HEAVY_ID"]),
                light_numbering_id=_parse_optional_str(row["LIGHT_ID"]),
                method=row["method"].strip(),
                resolution=_parse_optional_float(row["resolution"]),
                antibody_type=antibody_type,
                construct=row["construct"].strip(),
                holo=_parse_bool(row["holo"]),
                heavy_chain_author_id=_parse_optional_str(row["Hchain"]),
                light_chain_author_id=_parse_optional_str(row["Lchain"]),
                heavy_sequence=heavy_sequence,
                light_sequence=light_sequence,
                train_test_split=row["ab_split"].strip(),
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
