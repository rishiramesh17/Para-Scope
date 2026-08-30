"""One-off: build the committed human germline V-gene protein reference set.

Downloads OGRDB's published human IGH_VDJ / IGKappa_VJ / IGLambda_VJ germline
sets (nucleotide, ungapped), keeps only V genes (drops D/J), translates each
in frame 1 (IMGT V-REGION convention -- no leader, starts on a codon), and
writes a clean single-line-per-record protein FASTA.

OGRDB publishes under CC0 (public domain) -- see
https://ogrdb.airr-community.org/render_page/licensing_statement.html --
so, unlike SAbDab2/Thera-SAbDab, the output is committed to git and installed
with the package (src/para_scope/reference_data/), the same way `promb` ships
its own reference data inside its wheel, rather than gitignored under data/.
Rerun this script only when picking up a newer OGRDB release; bump
REFERENCE_SET_VERSION in para_scope/germline.py to match.
"""

import urllib.request
from pathlib import Path

from Bio.Seq import Seq

BASE_URL = "https://ogrdb.airr-community.org/download_germline_set/Homo%20sapiens"

SOURCES = {
    "IGH": (f"{BASE_URL}/IGH_VDJ/published/ungapped_ex", "V"),
    "IGK": (f"{BASE_URL}/IGKappa_VJ/published/ungapped_ex", "V"),
    "IGL": (f"{BASE_URL}/IGLambda_VJ/published/ungapped_ex", "V"),
}

OUT_PATH = Path("src/para_scope/reference_data/germline/ogrdb_human_v_genes.fasta")


def fetch_fasta(url: str) -> str:
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    seq_lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        elif line.strip():
            seq_lines.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


def translate_v_gene(nt_seq: str) -> str | None:
    usable_len = (len(nt_seq) // 3) * 3
    if usable_len < 60:  # shorter than ~20 aa is not a usable V-gene fragment
        return None
    protein = str(Seq(nt_seq[:usable_len]).translate())
    if "*" in protein or "X" in protein:
        return None  # pseudogene / ambiguous call -- not a clean reference entry
    return protein


def main() -> None:
    kept = 0
    skipped_non_v = 0
    skipped_untranslatable = 0
    out_records: list[tuple[str, str]] = []

    for locus, (url, gene_prefix) in SOURCES.items():
        text = fetch_fasta(url)
        for header, nt_seq in parse_fasta(text):
            if not header.startswith(f"IG{locus[2]}V"):
                skipped_non_v += 1
                continue
            protein = translate_v_gene(nt_seq)
            if protein is None:
                skipped_untranslatable += 1
                continue
            out_records.append((f"{header}|{locus}", protein))
            kept += 1

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w") as f:
        f.write(
            "; Source: OGRDB (https://ogrdb.airr-community.org), "
            "Homo sapiens IGH_VDJ / IGKappa_VJ / IGLambda_VJ, published sets\n"
        )
        f.write("; License: CC0 (public domain)\n")
        f.write("; V genes only, translated in frame 1, pseudogenes/ambiguous calls dropped\n")
        for header, protein in out_records:
            f.write(f">{header}\n{protein}\n")

    print(f"V genes kept: {kept}")
    print(f"Non-V headers skipped (D/J genes): {skipped_non_v}")
    print(f"Untranslatable (stop codon / ambiguous) skipped: {skipped_untranslatable}")
    print(f"Written to: {OUT_PATH}")


if __name__ == "__main__":
    main()
