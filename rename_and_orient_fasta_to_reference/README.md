# Rename and Orient Chromosomes

Rename and orient chromosomes in a FASTA file based on alignment to a reference genome. Reference header format like `>chrN` or `>chr_N` is expected in PAF file. 

## Requirements

- Python 3.7+
- [FastGA](https://github.com/thegenemyers/FASTGA) for alignment (or minimap2)

## Usage

```bash
# Basic usage (SUPER_* -> SUPER_*)
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output output_prefix

# Rename to chr format (SUPER_* -> chr1, chr2, chrW...)
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output output_prefix \
    --output-prefix chr

# From scaffold format to chr_ format
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output output_prefix \
    --query-prefix scaffold_ \
    --output-prefix chr_

# Just numbers (1, 2, 3, W, Z...)
python rename_and_orient.py \
    --fasta input.fa \
    --paf alignment.paf \
    --output output_prefix \
    --output-prefix ""
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
| `--output` | `-o` | required | Output prefix |
| `--min-coverage` | `-c` | 0.5 | Min alignment coverage for renaming |
| `--query-prefix` | `-q` | `SUPER_` | Input chromosome prefix |
| `--output-prefix` | `-x` | `SUPER_` | Output chromosome prefix |

## Output files

- `{output}.fa` — renamed FASTA with corrected orientation
- `{output}.chromosome.list.csv` — chromosome mapping (name, suffix, is_main)
- `{output}.mapping.tsv` — detailed alignment statistics

## Features

- Auto-detects reference format (`chrN` or `chr_N`)
- Determines orientation based on alignment strand statistics
- Handles sex chromosomes (W, X, Y, Z, Z1, Z2...)
- Resolves conflicts when multiple chromosomes map to same target
- Supports unlocalized contigs (`*_unloc_*`)
