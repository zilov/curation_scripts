# Rename and Orient Chromosomes

Rename and orient chromosomes in a FASTA file based on alignment to a reference genome.

## Requirements

- Python 3.7+
- `matplotlib` — optional, required only for `--plot-alignments`

## Usage

```bash
# Basic usage (SUPER_* -> SUPER_*)
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf

# Rename to chr format (SUPER_* -> chr1, chr2, chrW...)
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output-chromosome-prefix chr

# From scaffold format to chr_ format
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --query-chromosome-prefix scaffold_ \
    --output-chromosome-prefix chr_

# Just numbers (1, 2, 3, W, Z...)
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output-chromosome-prefix ""

# Custom output directory and prefix
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output-dir results \
    --output-prefix my_analysis

# Generate scatter plots of alignment blocks per chromosome
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --plot-alignments
```

## Creating alignment with FastGA

```bash
# Index reference
FastGA -v -1 reference.fa

# Align query to reference
FastGA -v -P query.fa reference.fa > query_vs_reference.paf
```

## Parameters

| Parameter | Short | Default | Description |
|-----------|-------|---------|-------------|
| `--fasta` | `-f` | required | Input FASTA (gzip supported) |
| `--paf` | `-p` | required | PAF alignment file |
| `--output-dir` | `-d` | `./rename_and_orient` | Output directory for generated files |
| `--output-prefix` | `-o` | input FASTA stem | Prefix for output file names |
| `--min-coverage` | `-c` | `0.5` | Min alignment coverage for renaming |
| `--query-chromosome-prefix` | `-q` | `SUPER_` | Prefix of chromosome names in the input FASTA |
| `--output-chromosome-prefix` | `-x` | `SUPER_` | Prefix of chromosome names in output FASTA |
| `--reference-chromosome-prefix` | `-r` | auto | Prefix of chromosome names in the PAF target (auto-detected if not set) |
| `--plot-alignments` | `-P` | off | Save scatter plots of PAF alignment blocks to `<output-dir>/plots/` (requires matplotlib) |

## Output files

- `{output_dir}/{output_prefix}.fa` — renamed FASTA with corrected orientation
- `{output_dir}/{output_prefix}.chromosome.list.csv` — chromosome mapping (name, suffix, is_main)
- `{output_dir}/{output_prefix}.mapping.tsv` — detailed alignment statistics
- `{output_dir}/plots/` — scatter plots per chromosome (only with `--plot-alignments`)

## Orientation algorithm

Chromosome orientation is decided using **Pearson correlation** between query and target alignment midpoints:

- For `+`-strand blocks: $r^+ = \text{corr}(\text{query\_mid},\ \text{target\_mid})$
- For `−`-strand blocks: $r^- = \text{corr}(\text{query\_mid},\ \text{target\_length} - \text{target\_mid})$

Syntenic blocks form a clean diagonal → $r \approx 1$.  
Structural inversions scatter randomly → $r$ stays low and does not influence the decision.  
The chromosome is reverse-complemented when $r^- > r^+$.

When fewer than 2 blocks exist on either strand (correlation is uninformative), the algorithm falls back to comparing total aligned lengths.

This approach correctly handles inter-species structural inversions that previously caused the naive strand-length heuristic to make wrong orientation calls.

## Features

- Auto-detects reference chromosome prefix (`chrN`, `chr_N`, `SUPER_N`, etc.)
- Orientation using Pearson correlation (handles large structural inversions)
- Handles sex chromosomes (W, X, Y, Z, Z1, Z2…)
- Resolves conflicts when multiple query chromosomes map to the same target
- Supports unlocalized contigs (`*_unloc_*`)
- Optional scatter plots for visual QC of alignment blocks per chromosome
