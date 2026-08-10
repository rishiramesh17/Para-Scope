"""Command-line interface for Para-Scope."""

import sys
import json
from para_scope.numbering import number_sequence


def main():
    """CLI entry point: number a sequence and report raw measurements."""
    if len(sys.argv) < 2:
        print("Usage: para-scope <sequence>")
        sys.exit(1)

    sequence = sys.argv[1]
    result = number_sequence(sequence)

    if not result:
        print(json.dumps({"error": "numbering failed"}))
        sys.exit(1)

    # Report: raw IMGT numbering with provenance
    output = {
        "sequence_checksum": result.sequence_checksum,
        "chain_type": result.chain_type,
        "scheme": result.scheme,
        "numbering_score": result.numbering_score,
        "tool": result.tool_name,
        "tool_version": result.tool_version,
        "positions": [
            {"position": pos, "insertion_code": code, "residue": res}
            for (pos, code), res in result.imgt_numbering
        ]
    }

    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
