#!/usr/bin/env python3
"""Extract IS110 element sequences from JSON into FASTA for minimap2 search."""

import json
import sys
from collections import Counter

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <is110_records.json> <output.fasta> [min_length]")
        sys.exit(1)

    json_path = sys.argv[1]
    fasta_path = sys.argv[2]
    min_length = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    with open(json_path) as f:
        data = json.load(f)

    org_counts = Counter()
    written = 0
    skipped_short = 0
    skipped_empty = 0

    with open(fasta_path, "w") as out:
        for r in data:
            seq = r["is_element"].get("sequence", "")
            if not seq:
                skipped_empty += 1
                continue
            if len(seq) < min_length:
                skipped_short += 1
                continue

            organism = r.get("organism", "unknown") or "unknown"
            sample_id = r.get("sample_id", "NA")
            is_id = r.get("is_id", "NA")
            length = r["is_element"]["length"]
            contig = r["is_element"].get("contig", "NA")
            start = r["is_element"].get("start", 0)
            end = r["is_element"].get("end", 0)
            strand = r["is_element"].get("strand", ".")

            header = (
                f">{is_id}"
                f" sample={sample_id}"
                f" organism={organism.replace(' ', '_')}"
                f" loc={contig}:{start}-{end}({strand})"
                f" len={length}"
            )
            out.write(header + "\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

            org_counts[organism] += 1
            written += 1

    print(f"Total records:  {len(data)}")
    print(f"Written:        {written}")
    print(f"Skipped empty:  {skipped_empty}")
    print(f"Skipped short:  {skipped_short} (< {min_length} bp)")
    print(f"\nTop organisms:")
    for org, cnt in org_counts.most_common(10):
        print(f"  {org}: {cnt}")


if __name__ == "__main__":
    main()
