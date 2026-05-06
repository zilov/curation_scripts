# PAF Filter — Nested & Short Alignment Removal

Reduce large PAF files for visualization without losing structural information.

## Problem

Whole-genome aligners like **FastGA** produce millions of small alignment records:
- Repetitive regions create thousands of short hits (<1 kb) to many targets
- Nested sub-alignments pile up inside larger alignment blocks
- Hits to unplaced scaffolds (`NW_`, `NT_`) are usually noise for visualization

A raw FastGA PAF can have 3–4 M lines while carrying only ~50k biologically
meaningful alignment blocks.

## Requirements

Python 3.7+, no external dependencies.

## Usage

```bash
# 1. Strongest reduction: keep only long hits to named chromosomes
python paf_filter.py \
    --input  assembly_vs_reference.paf \
    --output filtered.paf \
    --min-length 10000 \
    --target-prefix chr_

# 2. Also remove redundant nested alignments
python paf_filter.py \
    --input  assembly_vs_reference.paf \
    --output filtered.paf \
    --min-length 5000 \
    --target-prefix chr_ \
    --query-nested \
    --target-nested

# 3. Conservative — only deduplicate, keep everything else
python paf_filter.py \
    --input  assembly_vs_reference.paf \
    --output filtered.paf \
    --query-nested

# 4. Focus on specific query sequences
python paf_filter.py \
    --input  assembly_vs_reference.paf \
    --output filtered.paf \
    --query-prefix SUPER_ \
    --target-prefix chr_
```

## Parameters

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--input` | `-i` | required | Input PAF file |
| `--output` | `-o` | required | Output filtered PAF file |
| `--min-length` | `-l` | `0` (off) | Min query block length (bp) to keep |
| `--min-matches` | `-m` | `0` (off) | Min matching bases (PAF col 10) |
| `--min-mapq` | `-q` | `0` (off) | Min mapping quality (PAF col 12) |
| `--target-prefix` | `-t` | off | Keep only targets with this prefix |
| `--query-prefix` | | off | Keep only queries with this prefix |
| `--query-nested` | | off | Remove query-side nested alignments |
| `--target-nested` | | off | Remove target-side nested alignments |

## Filters explained

### `--min-length N`
The single most effective filter for FastGA output. Discards every record
whose query block length (`query_end - query_start`) is less than N bp.
Recommended: **5000–10000** for chromosome-level visualization.

**Risk:** losing genuine small syntenic blocks in rearranged or diverged regions.
Mitigate by lowering the threshold or combining with `--target-prefix` instead.

### `--target-prefix PREFIX`
Drops hits to unplaced scaffolds (`NW_`, `NT_`, etc.) and keeps only
named chromosomes (`chr_`, `SUPER_`, ...). Usually safe because unplaced
scaffolds represent tiny fragments already represented by the main-chromosome
hits.

**Risk:** missing assembly contigs that truly map only to unplaced scaffolds.
Usually acceptable for visualization.

### `--query-nested` / `--target-nested`
Within each **(query, target)** pair, removes an alignment if its query (or
target) coordinate interval is fully contained inside another alignment's
interval for the same pair.

Example: if SUPER_1 → chr_5 has alignments at [100, 5000] and [200, 4000],
the second is nested and gets removed.

**Risk:** losing secondary hits in repeat regions. For repeat-rich genomes this
is desirable for visualization. For structural-variant analysis you may want
to keep them.

## Recommended workflow for visualization

```
# Step 1: aggressive filter for quick visualization
python paf_filter.py -i raw.paf -o vis.paf \
    --min-length 10000 --target-prefix chr_

# Step 2: if still too large, add nested removal
python paf_filter.py -i vis.paf -o vis2.paf \
    --query-nested --target-nested

# Visualize with e.g. D-GENIES or NucPlot
```
