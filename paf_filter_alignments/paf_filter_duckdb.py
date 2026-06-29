#!/usr/bin/env python3
"""
Filter PAF alignments using DuckDB for fast in-memory SQL processing.

Applies four complementary filters (all optional, combinable):
  1. --min-length    : discard alignments shorter than N bp (query-side block length)
  2. --min-matches   : minimum number of matching bases (column 10)
  3. --min-mapq      : minimum mapping quality (column 12)
  4. --target-prefix : keep only records whose target name starts with this prefix
  5. --query-prefix  : keep only records whose query name starts with this prefix
  6. --query-nested  : remove alignments whose query interval is fully contained
                       inside a longer alignment to the same (query, target) pair
  7. --target-nested : same on the target-coordinate side

Usage examples
--------------
python paf_filter_duckdb.py -i in.paf -o out.paf --min-length 10000 --target-prefix chr_
python paf_filter_duckdb.py -i in.paf -o out.paf --query-nested --target-nested
"""

__version__ = "1.0.0"

import argparse
import sys
import duckdb
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter PAF alignments using DuckDB",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", "-v", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input",  "-i", required=True, type=Path, help="Input PAF file")
    parser.add_argument("--output", "-o", required=True, type=Path, help="Output PAF file")

    parser.add_argument("--min-length",  "-l", type=int, default=0,
                        help="Minimum query block length (bp). 0 = disabled.")
    parser.add_argument("--min-matches", "-m", type=int, default=0,
                        help="Minimum matching bases (col 10). 0 = disabled.")
    parser.add_argument("--min-mapq",   "-q", type=int, default=0,
                        help="Minimum mapping quality (col 12). 0 = disabled.")
    parser.add_argument("--target-prefix", "-t", type=str, default=None,
                        help="Keep only records whose target name starts with this prefix.")
    parser.add_argument("--query-prefix", type=str, default=None,
                        help="Keep only records whose query name starts with this prefix.")
    parser.add_argument("--query-nested",  action="store_true",
                        help="Remove alignments fully nested on the query side.")
    parser.add_argument("--target-nested", action="store_true",
                        help="Remove alignments fully nested on the target side.")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def _detect_n_cols(path: Path) -> int:
    """Count tab-separated columns in the first non-empty line."""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                return len(line.split("\t"))
    return 12


def _col(i: int, n_cols: int) -> str:
    """Return DuckDB auto-generated column name (zero-padded to match column count width)."""
    width = len(str(n_cols - 1))
    return f"column{i:0{width}d}"


def load_and_filter_scalars(
    con: duckdb.DuckDBPyConnection,
    path: Path,
    args: argparse.Namespace,
) -> tuple[int, int]:
    """
    Stream PAF file through DuckDB, applying scalar filters on the fly so only
    passing rows are materialised. This keeps memory proportional to the output,
    not the (potentially multi-GB) input.

    Returns (total_input_records, filtered_record_count).
    """
    n_cols = _detect_n_cols(path)
    c = [_col(i, n_cols) for i in range(n_cols)]
    raw_line_expr = "concat_ws('\t', " + ", ".join(c) + ")"

    # Inline view alias so we can use computed aliases in WHERE
    scalar_where = _build_where(args)

    con.execute(f"""
        CREATE TABLE filtered AS
        SELECT
            {c[0]}                        AS query_name,
            {c[1]}::INT                   AS query_length,
            {c[2]}::INT                   AS query_start,
            {c[3]}::INT                   AS query_end,
            {c[4]}                        AS strand,
            {c[5]}                        AS target_name,
            {c[6]}::INT                   AS target_length,
            {c[7]}::INT                   AS target_start,
            {c[8]}::INT                   AS target_end,
            {c[9]}::INT                   AS num_matches,
            {c[10]}::INT                  AS alignment_length,
            {c[11]}::INT                  AS mapping_quality,
            {c[3]}::INT - {c[2]}::INT     AS query_block_len,
            {c[8]}::INT - {c[7]}::INT     AS target_block_len,
            {raw_line_expr}               AS raw_line
        FROM read_csv('{path}', sep='\t', header=false, all_varchar=true)
        WHERE {scalar_where}
    """)

    n_filtered = con.execute("SELECT COUNT(*) FROM filtered").fetchone()[0]

    # Count total separately with a second streaming pass (cheap — no materialisation)
    n_total = con.execute(f"""
        SELECT COUNT(*) FROM read_csv('{path}', sep='\t', header=false, all_varchar=true)
    """).fetchone()[0]

    return n_total, n_filtered


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------

def _build_where(args: argparse.Namespace) -> str:
    """Build SQL WHERE clause for scalar filters."""
    conditions = []
    if args.query_prefix:
        conditions.append(f"query_name LIKE '{args.query_prefix}%'")
    if args.target_prefix:
        conditions.append(f"target_name LIKE '{args.target_prefix}%'")
    if args.min_mapq > 0:
        conditions.append(f"mapping_quality >= {args.min_mapq}")
    if args.min_matches > 0:
        conditions.append(f"num_matches >= {args.min_matches}")
    if args.min_length > 0:
        conditions.append(f"query_block_len >= {args.min_length}")
    return " AND ".join(conditions) if conditions else "TRUE"




def apply_nested_filter(con: duckdb.DuckDBPyConnection, side: str) -> None:
    """
    Remove records whose interval on `side` ('query' or 'target') is fully
    contained within a longer alignment to the same (query_name, target_name).

    Uses NOT EXISTS: for each record, check whether any other record in the
    same (query, target) group has a strictly longer block that covers it.
    """
    start     = f"{side}_start"
    end       = f"{side}_end"
    block_len = f"{side}_block_len"

    con.execute(f"""
        CREATE TABLE filtered_next AS
        SELECT f1.* FROM filtered f1
        WHERE NOT EXISTS (
            SELECT 1 FROM filtered f2
            WHERE f2.query_name  = f1.query_name
              AND f2.target_name = f1.target_name
              AND f2.{block_len} > f1.{block_len}
              AND f2.{start}    <= f1.{start}
              AND f2.{end}      >= f1.{end}
        )
    """)
    con.execute("DROP TABLE filtered")
    con.execute("ALTER TABLE filtered_next RENAME TO filtered")


# ---------------------------------------------------------------------------
# Stats and output
# ---------------------------------------------------------------------------

def _stats(con: duckdb.DuckDBPyConnection, total: int) -> tuple:
    return con.execute("""
        SELECT COUNT(*),
               MIN(query_block_len),
               MEDIAN(query_block_len)::INT,
               MAX(query_block_len)
        FROM filtered
    """).fetchone()


def print_stats(con: duckdb.DuckDBPyConnection, label: str, total: int) -> None:
    n, mn, med, mx = _stats(con, total)
    pct = 100.0 * n / total if total else 0.0
    if n == 0:
        print(f"  [{label}]         0 records kept (0.0%) | (empty)")
    else:
        print(f"  [{label}] {n:>9,} records kept ({pct:.1f}%) "
              f"| query-block — min:{mn:,}  median:{med:,}  max:{mx:,}")


def write_output(con: duckdb.DuckDBPyConnection, path: Path, batch_size: int = 50_000) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cursor = con.execute("SELECT raw_line FROM filtered")
    with open(path, "w") as fh:
        while True:
            batch = cursor.fetchmany(batch_size)
            if not batch:
                break
            fh.writelines(line + "\n" for (line,) in batch)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    if not args.input.exists():
        print(f"Error: input PAF not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    no_filter = (
        args.min_length == 0 and args.min_matches == 0 and args.min_mapq == 0
        and args.target_prefix is None and args.query_prefix is None
        and not args.query_nested and not args.target_nested
    )
    if no_filter:
        print("Warning: no filter options specified — output will be identical to input.",
              file=sys.stderr)

    con = duckdb.connect()

    print(f"Reading PAF: {args.input}")
    total, _ = load_and_filter_scalars(con, args.input, args)
    print(f"  Total records: {total:,}")

    print("\nApplying filters:")
    where = _build_where(args)
    if where != "TRUE":
        print_stats(con, "scalar filters", total)

    if args.query_nested:
        before = con.execute("SELECT COUNT(*) FROM filtered").fetchone()[0]
        apply_nested_filter(con, "query")
        after = con.execute("SELECT COUNT(*) FROM filtered").fetchone()[0]
        print_stats(con, f"--query-nested (removed {before - after:,})", total)

    if args.target_nested:
        before = con.execute("SELECT COUNT(*) FROM filtered").fetchone()[0]
        apply_nested_filter(con, "target")
        after = con.execute("SELECT COUNT(*) FROM filtered").fetchone()[0]
        print_stats(con, f"--target-nested (removed {before - after:,})", total)

    final = con.execute("SELECT COUNT(*) FROM filtered").fetchone()[0]
    print(f"\nFinal: {final:,} / {total:,} records kept "
          f"({100.0 * final / total:.1f}% of input)")

    write_output(con, args.output)
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
