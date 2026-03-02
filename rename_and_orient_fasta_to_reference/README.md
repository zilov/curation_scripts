# Rename and Orient Chromosomes

Rename and orient chromosomes in a FASTA file based on alignment to a reference genome.

## Requirements

- Python 3.7+
- [FastGA](https://github.com/thegenemyers/FASTGA) for alignment (or minimap2)

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
| `--min-coverage` | `-c` | 0.5 | Min alignment coverage for renaming |
| `--query-chromosome-prefix` | `-q` | `SUPER_` | Input chromosome prefix |
| `--output-chromosome-prefix` | `-x` | `SUPER_` | Output chromosome prefix |

## Output files

- `{output_dir}/{output_prefix}.fa` — renamed FASTA with corrected orientation
- `{output_dir}/{output_prefix}.chromosome.list.csv` — chromosome mapping (name, suffix, is_main)
- `{output_dir}/{output_prefix}.mapping.tsv` — detailed alignment statistics

## Features

- Auto-detects reference format (`chrN` or `chr_N`)
- Determines orientation based on alignment strand statistics
- Handles sex chromosomes (W, X, Y, Z, Z1, Z2...)
- Resolves conflicts when multiple chromosomes map to same target
- Supports unlocalized contigs (`*_unloc_*`)
