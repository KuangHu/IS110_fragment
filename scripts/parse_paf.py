"""Parse PAF alignment files from minimap2.

Provides both a CLI and importable functions for filtering and grouping hits.

Usage:
    python parse_paf.py hits.paf --min-identity 0.80 --group-by-genome
    python parse_paf.py hits.paf --target NZ_CP084874.1
"""

import argparse
import json
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass, field


@dataclass
class PafHit:
    query_name: str
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_len: int
    target_start: int
    target_end: int
    num_matches: int
    aln_block_len: int
    mapq: int
    identity: float
    divergence: float
    cigar: str | None = None
    tags: dict = field(default_factory=dict)


def _parse_tags(fields: list[str]) -> dict:
    """Parse optional PAF tags (col 13+) into a dict."""
    tags = {}
    for f in fields:
        parts = f.split(":", 2)
        if len(parts) == 3:
            key, typ, val = parts
            if typ == "i":
                tags[key] = int(val)
            elif typ == "f":
                tags[key] = float(val)
            else:
                tags[key] = val
    return tags


def _compute_identity(num_matches: int, aln_block_len: int, tags: dict) -> tuple[float, float]:
    """Compute identity and divergence from tags or columns 10/11.

    Prefers de:f: (per-base divergence) over dv:f: (gap-compressed),
    falls back to col10/col11.
    """
    # With a CIGAR: matches / aligned columns (M, =, X), the definition used by
    # every detector here (paf_utils.py). The fallbacks below are NOT the same
    # quantity: de counts each gap as one difference, col10/col11 counts every
    # gap column (a perfect copy with a 600-bp insert reads 77%).
    cg = tags.get("cg")
    if isinstance(cg, str) and cg:
        import re
        aligned = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg) if op in "M=X")
        if aligned > 0:
            identity = min(1.0, num_matches / aligned)
            return identity, 1.0 - identity
    if "de" in tags:
        div = tags["de"]
    elif "dv" in tags:
        div = tags["dv"]
    elif aln_block_len > 0:
        div = 1.0 - num_matches / aln_block_len
    else:
        div = 1.0
    identity = 1.0 - div
    return identity, div


def parse_paf(path: str) -> list[PafHit]:
    """Parse a PAF file into a list of PafHit objects."""
    hits = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue

            tags = _parse_tags(cols[12:])
            num_matches = int(cols[9])
            aln_block_len = int(cols[10])
            identity, divergence = _compute_identity(num_matches, aln_block_len, tags)

            hits.append(PafHit(
                query_name=cols[0],
                query_len=int(cols[1]),
                query_start=int(cols[2]),
                query_end=int(cols[3]),
                strand=cols[4],
                target_name=cols[5],
                target_len=int(cols[6]),
                target_start=int(cols[7]),
                target_end=int(cols[8]),
                num_matches=num_matches,
                aln_block_len=aln_block_len,
                mapq=int(cols[11]),
                identity=identity,
                divergence=divergence,
                cigar=tags.get("cg"),
                tags=tags,
            ))
    return hits


def filter_hits(
    hits: list[PafHit],
    min_identity: float = 0.0,
    min_length: int = 0,
    target: str | None = None,
) -> list[PafHit]:
    """Filter hits by identity, alignment length, and/or target genome."""
    filtered = hits
    if min_identity > 0:
        filtered = [h for h in filtered if h.identity >= min_identity]
    if min_length > 0:
        filtered = [h for h in filtered if (h.query_end - h.query_start) >= min_length]
    if target:
        filtered = [h for h in filtered if h.target_name == target]
    return filtered


def group_by_genome(hits: list[PafHit]) -> dict[tuple, list[PafHit]]:
    """Group hits by (target genome, query). Hits of different queries were
    pooled per genome and divided by the length of the first query: two full
    copies of two queries read "200% coverage"."""
    groups = defaultdict(list)
    for h in hits:
        groups[(h.target_name, h.query_name)].append(h)
    return dict(groups)


def summarize_groups(groups: dict[str, list[PafHit]]) -> list[dict]:
    """Summarize grouped hits into a table."""
    rows = []
    for (genome, query), ghits in groups.items():
        total_bp = sum(h.query_end - h.query_start for h in ghits)
        mean_ident = sum(h.identity for h in ghits) / len(ghits) if ghits else 0
        query_len = ghits[0].query_len
        # coverage = part of the query seen at least once (merged intervals);
        # total_aligned_bp / n_hits carry the copy number
        merged, end = 0, 0
        for s, e in sorted((h.query_start, h.query_end) for h in ghits):
            merged += max(0, e - max(s, end))
            end = max(end, e)
        coverage = merged / query_len * 100 if query_len > 0 else 0
        rows.append({
            "genome": genome,
            "query": query,
            "n_hits": len(ghits),
            "total_aligned_bp": total_bp,
            "mean_identity": round(mean_ident * 100, 1),
            "query_coverage_pct": round(coverage, 1),
        })
    rows.sort(key=lambda r: r["total_aligned_bp"], reverse=True)
    return rows


def print_hits_tsv(hits: list[PafHit]):
    """Print hits as TSV."""
    header = "query_start\tquery_end\tstrand\ttarget\ttarget_start\ttarget_end\tidentity\tmapq"
    print(header)
    for h in sorted(hits, key=lambda x: x.query_start):
        print(f"{h.query_start}\t{h.query_end}\t{h.strand}\t{h.target_name}\t"
              f"{h.target_start}\t{h.target_end}\t{h.identity:.3f}\t{h.mapq}")


def print_summary_tsv(rows: list[dict]):
    """Print group summary as TSV."""
    header = "genome\tquery\tn_hits\ttotal_aligned_bp\tmean_identity\tquery_coverage_pct"
    print(header)
    for r in rows:
        print(f"{r['genome']}\t{r['query']}\t{r['n_hits']}\t{r['total_aligned_bp']}\t"
              f"{r['mean_identity']}\t{r['query_coverage_pct']}")


def main():
    parser = argparse.ArgumentParser(description="Parse and filter PAF alignment files")
    parser.add_argument("paf_file", help="Input PAF file")
    parser.add_argument("--min-identity", type=float, default=0.0,
                        help="Minimum identity threshold (0-1, default: 0.0)")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Minimum alignment length in bp (default: 0)")
    parser.add_argument("--target", help="Filter to hits against this target genome")
    parser.add_argument("--group-by-genome", action="store_true",
                        help="Group and summarize hits by target genome")
    parser.add_argument("--format", choices=["tsv", "json"], default="tsv",
                        help="Output format (default: tsv)")

    args = parser.parse_args()

    hits = parse_paf(args.paf_file)
    hits = filter_hits(hits, min_identity=args.min_identity,
                       min_length=args.min_length, target=args.target)

    if not hits:
        print("No hits found matching filters.", file=sys.stderr)
        sys.exit(0)

    if args.group_by_genome:
        groups = group_by_genome(hits)
        rows = summarize_groups(groups)
        if args.format == "json":
            print(json.dumps(rows, indent=2))
        else:
            print_summary_tsv(rows)
    else:
        if args.format == "json":
            print(json.dumps([asdict(h) for h in hits], indent=2))
        else:
            print_hits_tsv(hits)


if __name__ == "__main__":
    main()
