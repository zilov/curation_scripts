# Rename and Orient Chromosomes

Rename and orient chromosomes in a FASTA file based on alignment to a reference genome.

## Requirements

- Python 3.7+
- `matplotlib` — optional, required only for `--plot-alignments`

## Usage

### Mode 1: PAF alignment (standard, haplotype 1)

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

### Mode 2: Mapping table (haplotype 2, no re-alignment needed)

If you have already run the script on haplotype 1 and want to apply the same
renaming and orientation to haplotype 2 (which matches haplotype 1 by numbering),
pass the mapping TSV produced by the first run instead of a PAF file:

```bash
python rename_and_orient.py \
    --fasta hap2.fa \
    --mapping-table hap1_output/hap1.mapping.tsv \
    --output-dir hap2_output \
    --output-prefix hap2
```

- Chromosome renaming and RC flags are read directly from the table — no alignment required.
- Unlocalized contigs (`*_unloc_*`) present only in haplotype 2 are automatically
  mapped via their parent chromosome's new name. They are **never** reverse-complemented.
- `--min-coverage`, `--plot-alignments`, and `--reference-chromosome-prefix` are
  ignored in this mode (a warning is printed if they are supplied).

> `--paf` and `--mapping-table` are mutually exclusive — exactly one must be provided.

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
| `--paf` | `-p` | — | PAF alignment file *(mutually exclusive with `--mapping-table`)* |
| `--mapping-table` | `-mt` | — | Mapping TSV from a previous run *(mutually exclusive with `--paf`)* |
| `--output-dir` | `-d` | `./rename_and_orient` | Output directory for generated files |
| `--output-prefix` | `-o` | input FASTA stem | Prefix for output file names |
| `--min-coverage` | `-c` | `0.5` | Min alignment coverage for renaming *(PAF mode only)* |
| `--query-chromosome-prefix` | `-q` | `SUPER_` | Prefix of chromosome names in the input FASTA |
| `--output-chromosome-prefix` | `-x` | `SUPER_` | Prefix of chromosome names in output FASTA |
| `--reference-chromosome-prefix` | `-r` | auto | Prefix of chromosome names in the PAF target *(PAF mode only)* |
| `--plot-alignments` | `-P` | off | Save scatter plots of PAF alignment blocks *(PAF mode only, requires matplotlib)* |

## Output files

- `{output_dir}/{output_prefix}.fa` — renamed FASTA with corrected orientation
- `{output_dir}/{output_prefix}.chromosome.list.csv` — chromosome mapping (name, suffix, is_main)
- `{output_dir}/{output_prefix}.mapping.tsv` — detailed alignment statistics *(PAF mode only)*
- `{output_dir}/plots/` — scatter plots per chromosome *(PAF mode + `--plot-alignments` only)*

## Orientation algorithm

Chromosome orientation is decided using **Pearson correlation** between query and target alignment midpoints:

- For `+`-strand blocks: r⁺ = corr(query_mid, target_mid)
- For `−`-strand blocks: r⁻ = corr(query_mid, target_length − target_mid)

Syntenic blocks form a clean diagonal → r ≈ 1.  
Structural inversions scatter randomly → r stays low and does not influence the decision.  
The chromosome is reverse-complemented when r⁻ > r⁺.

When fewer than 2 blocks exist on either strand (correlation is uninformative), the algorithm falls back to comparing total aligned lengths.

This approach correctly handles inter-species structural inversions that previously caused the naive strand-length heuristic to make wrong orientation calls.

## Features

- Auto-detects reference chromosome prefix (`chrN`, `chr_N`, `SUPER_N`, etc.)
- Orientation using Pearson correlation (handles large structural inversions)
- Handles sex chromosomes (W, X, Y, Z, Z1, Z2…)
- Resolves conflicts when multiple query chromosomes map to the same target
- Supports unlocalized contigs (`*_unloc_*`)
- **Mapping-table mode** — apply renaming to a second haplotype without re-running alignment
- Optional scatter plots for visual QC of alignment blocks per chromosome
