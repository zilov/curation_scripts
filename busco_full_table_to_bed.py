#!/usr/bin/env python3
"""Parse BUSCO full_table.tsv and write Complete/Fragmented/Duplicated BEDGraph files."""

import argparse
import pandas as pd


BUSCO_COLS = ["busco_id", "status", "sequence", "start", "end", "strand", "score", "length", "ortho_url", "description"]
BEDGRAPH_COLS = ["sequence", "start", "end", "score"]
STATUS_MAP   = {
    "Complete":   "complete_buscos.bedgraph",
    "Duplicated": "duplicated_buscos.bedgraph",
    "Fragmented": "fragmented_buscos.bedgraph",
}


def parse_args():
    p = argparse.ArgumentParser(description="BUSCO full_table → BEDGraph files")
    p.add_argument("full_table", help="Path to BUSCO full_table.tsv")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: .)")
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(
        args.full_table,
        sep="\t",
        comment="#",
        header=None,
        names=BUSCO_COLS,
    )

    df = df[df["status"].isin(STATUS_MAP)].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df.dropna(subset=["start", "end", "score"], inplace=True)
    df[["start", "end"]] = df[["start", "end"]].astype(int)
    df["score"] = df["score"].astype(int)

    for status, filename in STATUS_MAP.items():
        subset = df[df["status"] == status][BEDGRAPH_COLS]
        outpath = f"{args.outdir}/{filename}"
        subset.to_csv(outpath, sep="\t", index=False, header=False)
        print(f"{status}: {len(subset)} entries → {outpath}")


if __name__ == "__main__":
    main()
