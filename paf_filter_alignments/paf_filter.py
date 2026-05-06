#!/usr/bin/env python3
"""
Filter PAF alignments to reduce file size while preserving biological signal.

Applies three complementary filters (all optional, combinable):
  1. --min-length    : discard alignments shorter than N bp (query-side block length)
  2. --query-nested  : remove alignments whose query interval is fully contained
                       inside a longer alignment to the same (query, target) pair
  3. --target-nested : same on the target-coordinate side
  4. --target-prefix : keep only records whose target name starts with this prefix
                       (e.g. "chr_" to drop unplaced NW_/NT_ scaffolds)

Usage examples
--------------
# Aggressive: keep only long hits to named chromosomes
python paf_filter.py -i in.paf -o out.paf --min-length 10000 --target-prefix chr_

# Conservative: only remove nested redundancy
python paf_filter.py -i in.paf -o out.paf --query-nested --target-nested

# Combined
python paf_filter.py -i in.paf -o out.paf \\
    --min-length 5000 --query-nested --target-prefix chr_
"""

__version__ = "1.0.0"

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class PAFRecord:
    """Single PAF alignment record (12 mandatory columns + raw line)."""
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    num_matches: int
    alignment_length: int
    mapping_quality: int
    raw_line: str          # full original line — preserved verbatim on output

    @property
    def query_block_length(self) -> int:
        return self.query_end - self.query_start

    @property
    def target_block_length(self) -> int:
        return self.target_end - self.target_start


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter PAF alignments: remove short / nested records",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--version", "-v", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument("--input", "-i", required=True, type=Path,
                        help="Input PAF file")
    parser.add_argument("--output", "-o", required=True, type=Path,
                        help="Output (filtered) PAF file")

    # --- filter options ---
    parser.add_argument("--min-length", "-l", type=int, default=0,
                        help="Minimum alignment block length on the query (bp). "
                             "Records shorter than this are discarded. "
                             "0 = disabled.")
    parser.add_argument("--min-matches", "-m", type=int, default=0,
                        help="Minimum number of matching bases (column 10). "
                             "0 = disabled.")
    parser.add_argument("--min-mapq", "-q", type=int, default=0,
                        help="Minimum mapping quality (column 12). "
                             "0 = disabled.")
    parser.add_argument("--target-prefix", "-t", type=str, default=None,
                        help="Keep only records whose target name starts with "
                             "this prefix (e.g. 'chr_', 'chr'). "
                             "Unset = keep all targets.")
    parser.add_argument("--query-prefix", type=str, default=None,
                        help="Keep only records whose query name starts with "
                             "this prefix (e.g. 'SUPER_'). "
                             "Unset = keep all queries.")
    parser.add_argument("--query-nested", action="store_true",
                        help="Remove alignments whose query interval is fully "
                             "contained within a longer alignment to the same "
                             "(query, target) pair.")
    parser.add_argument("--target-nested", action="store_true",
                        help="Remove alignments whose target interval is fully "
                             "contained within a longer alignment to the same "
                             "(query, target) pair.")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# PAF I/O
# ---------------------------------------------------------------------------

def parse_paf(path: Path) -> List[PAFRecord]:
    """Read all PAF records from *path*."""
    records: List[PAFRecord] = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            raw = line.rstrip("\n")
            if not raw:
                continue
            fields = raw.split("\t")
            if len(fields) < 12:
                print(f"  Warning: line {lineno} has fewer than 12 fields, skipping",
                      file=sys.stderr)
                continue
            records.append(PAFRecord(
                query_name=fields[0],
                query_length=int(fields[1]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target_name=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                num_matches=int(fields[9]),
                alignment_length=int(fields[10]),
                mapping_quality=int(fields[11]),
                raw_line=raw,
            ))
    return records


def write_paf(records: List[PAFRecord], path: Path) -> None:
    """Write PAF records (raw lines) to *path*."""
    with open(path, "w") as fh:
        for rec in records:
            fh.write(rec.raw_line + "\n")


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------

def filter_by_length(records: List[PAFRecord], min_length: int) -> List[PAFRecord]:
    """Keep records whose query block length >= min_length."""
    return [r for r in records if r.query_block_length >= min_length]


def filter_by_matches(records: List[PAFRecord], min_matches: int) -> List[PAFRecord]:
    return [r for r in records if r.num_matches >= min_matches]


def filter_by_mapq(records: List[PAFRecord], min_mapq: int) -> List[PAFRecord]:
    return [r for r in records if r.mapping_quality >= min_mapq]


def filter_by_target_prefix(records: List[PAFRecord], prefix: str) -> List[PAFRecord]:
    """Keep only records whose target name starts with *prefix*."""
    return [r for r in records if r.target_name.startswith(prefix)]


def filter_by_query_prefix(records: List[PAFRecord], prefix: str) -> List[PAFRecord]:
    """Keep only records whose query name starts with *prefix*."""
    return [r for r in records if r.query_name.startswith(prefix)]


def _remove_nested_intervals(
    records: List[PAFRecord],
    get_interval,          # callable(record) -> (start, end)
    get_key,               # callable(record) -> grouping key
) -> List[PAFRecord]:
    """
    Generic helper: within each group, remove records whose interval is fully
    contained inside the interval of a longer record in the same group.

    An interval (a_s, a_e) is nested inside (b_s, b_e) if:
        b_s <= a_s  AND  a_e <= b_e  AND  (a_s, a_e) != (b_s, b_e)

    We keep the *longer* record when two intervals are identical.

    Strategy (O(n log n) per group):
      1. Sort by interval length descending.
      2. Walk through sorted list; for each record, check whether its interval
         is fully covered by any previously accepted interval.  This is done
         efficiently by maintaining a list of accepted (start, end) pairs and
         checking containment.  For large groups we use a sweep-line approach.
    """
    # Group records by key
    groups: dict = defaultdict(list)
    for rec in records:
        groups[get_key(rec)].append(rec)

    kept: List[PAFRecord] = []

    for key, grp in groups.items():
        if len(grp) == 1:
            kept.append(grp[0])
            continue

        # Sort by interval length descending so we process longer first
        grp_sorted = sorted(grp, key=lambda r: get_interval(r)[1] - get_interval(r)[0],
                            reverse=True)

        accepted_intervals: List[Tuple[int, int]] = []

        for rec in grp_sorted:
            s, e = get_interval(rec)
            is_nested = any(bs <= s and e <= be for bs, be in accepted_intervals)
            if not is_nested:
                accepted_intervals.append((s, e))
                kept.append(rec)

    return kept


def filter_query_nested(records: List[PAFRecord]) -> List[PAFRecord]:
    """
    Remove records whose query interval is fully contained in a longer record
    for the same (query_name, target_name) pair.
    """
    return _remove_nested_intervals(
        records,
        get_interval=lambda r: (r.query_start, r.query_end),
        get_key=lambda r: (r.query_name, r.target_name),
    )


def filter_target_nested(records: List[PAFRecord]) -> List[PAFRecord]:
    """
    Remove records whose target interval is fully contained in a longer record
    for the same (query_name, target_name) pair.
    """
    return _remove_nested_intervals(
        records,
        get_interval=lambda r: (r.target_start, r.target_end),
        get_key=lambda r: (r.query_name, r.target_name),
    )


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def print_stats(label: str, records: List[PAFRecord], total: int) -> None:
    if not records:
        print(f"  [{label}] 0 records (100% removed)")
        return
    n = len(records)
    lengths = [r.query_block_length for r in records]
    median = sorted(lengths)[len(lengths) // 2]
    pct_kept = 100.0 * n / total if total else 0.0
    print(f"  [{label}] {n:>9,} records kept  ({pct_kept:.1f}%)  "
          f"| query-block len — min:{min(lengths):,}  "
          f"median:{median:,}  max:{max(lengths):,}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if not args.input.exists():
        print(f"Error: input PAF not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Validate that at least one filter is active
    no_filter = (
        args.min_length == 0
        and args.min_matches == 0
        and args.min_mapq == 0
        and args.target_prefix is None
        and args.query_prefix is None
        and not args.query_nested
        and not args.target_nested
    )
    if no_filter:
        print("Warning: no filter options specified — output will be identical to input.",
              file=sys.stderr)

    print(f"Reading PAF: {args.input}")
    records = parse_paf(args.input)
    total = len(records)
    print(f"  Total records: {total:,}")
    print_stats("input", records, total)

    # --- apply filters in order ---
    print("\nApplying filters:")

    if args.query_prefix:
        records = filter_by_query_prefix(records, args.query_prefix)
        print_stats(f"--query-prefix {args.query_prefix!r}", records, total)

    if args.target_prefix:
        records = filter_by_target_prefix(records, args.target_prefix)
        print_stats(f"--target-prefix {args.target_prefix!r}", records, total)

    if args.min_mapq > 0:
        records = filter_by_mapq(records, args.min_mapq)
        print_stats(f"--min-mapq {args.min_mapq}", records, total)

    if args.min_matches > 0:
        records = filter_by_matches(records, args.min_matches)
        print_stats(f"--min-matches {args.min_matches}", records, total)

    if args.min_length > 0:
        records = filter_by_length(records, args.min_length)
        print_stats(f"--min-length {args.min_length}", records, total)

    if args.query_nested:
        print("  Removing query-nested alignments...")
        before = len(records)
        records = filter_query_nested(records)
        print_stats(f"--query-nested (removed {before - len(records):,})", records, total)

    if args.target_nested:
        print("  Removing target-nested alignments...")
        before = len(records)
        records = filter_target_nested(records)
        print_stats(f"--target-nested (removed {before - len(records):,})", records, total)

    print(f"\nFinal: {len(records):,} / {total:,} records kept "
          f"({100.0 * len(records) / total:.1f}% of input)")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_paf(records, args.output)
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
