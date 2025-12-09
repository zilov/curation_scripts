#!/usr/bin/env python3
"""
Rename and orient chromosomes based on reference alignment (PAF format).

This script renames chromosomes in a FASTA file and changes their orientation
based on alignment to a reference genome.
"""

import argparse
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional


@dataclass
class PAFRecord:
    """Single PAF alignment record."""
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
    
    @property
    def alignment_block_length(self) -> int:
        """Length of the alignment block on query."""
        return self.query_end - self.query_start


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Rename and orient chromosomes based on reference alignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        type=Path,
        help="Input FASTA file (can be gzipped)"
    )
    
    parser.add_argument(
        "--paf", "-p",
        required=True,
        type=Path,
        help="PAF file with alignment to reference"
    )
    
    parser.add_argument(
        "--output", "-o",
        required=True,
        type=str,
        help="Output prefix for generated files"
    )
    
    parser.add_argument(
        "--min-coverage", "-c",
        type=float,
        default=0.5,
        help="Minimum coverage threshold for renaming (0.0-1.0)"
    )
    
    parser.add_argument(
        "--query-prefix", "-q",
        type=str,
        default="SUPER_",
        help="Prefix for query chromosome names in input FASTA (e.g., SUPER_, scaffold_, contig)"
    )
    
    parser.add_argument(
        "--output-prefix", "-x",
        type=str,
        default="SUPER_",
        help="Prefix for output chromosome names (e.g., SUPER_, chr_, chr, or empty string for no prefix)"
    )
    
    return parser.parse_args()


def read_fasta(fasta_path: Path) -> Dict[str, str]:
    """
    Read FASTA file (supports gzip compression).
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    # Determine if file is gzipped
    open_func = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'
    
    with open_func(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq)
                
                # Start new sequence - take first word as name
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences


def parse_paf(paf_path: Path) -> List[PAFRecord]:
    """
    Parse PAF file.
    
    PAF format (tab-separated):
    1. Query sequence name
    2. Query sequence length
    3. Query start (0-based)
    4. Query end (0-based, open)
    5. Strand ('+' or '-')
    6. Target sequence name
    7. Target sequence length
    8. Target start (0-based)
    9. Target end (0-based, open)
    10. Number of matching bases
    11. Alignment block length
    12. Mapping quality (0-255)
    
    Args:
        paf_path: Path to PAF file
        
    Returns:
        List of PAFRecord objects
    """
    records = []
    
    with open(paf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) < 12:
                continue
            
            record = PAFRecord(
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
                mapping_quality=int(fields[11])
            )
            records.append(record)
    
    return records


def detect_reference_prefix(records: List[PAFRecord]) -> str:
    """
    Auto-detect reference chromosome prefix (chr_ or chr).
    
    Args:
        records: List of PAF records
        
    Returns:
        Detected prefix ('chr_' or 'chr')
    """
    for record in records:
        target = record.target_name
        if target.startswith('chr_'):
            return 'chr_'
        elif target.startswith('chr') and len(target) > 3:
            # Check if it's chrN format (chr1, chr2, chrW, etc.)
            suffix = target[3:]
            if suffix[0].isdigit() or suffix[0].upper() in 'WXYZ':
                return 'chr'
    return 'chr_'  # Default


def filter_paf_records(records: List[PAFRecord], 
                       query_prefix: str = "SUPER_",
                       target_prefix: str = None) -> Tuple[List[PAFRecord], str]:
    """
    Filter PAF records to keep only those matching query and target prefixes.
    Auto-detects target prefix if not specified.
    
    Args:
        records: List of PAF records
        query_prefix: Required prefix for query names
        target_prefix: Required prefix for target names (auto-detected if None)
        
    Returns:
        Tuple of (filtered records, detected target prefix)
    """
    if target_prefix is None:
        target_prefix = detect_reference_prefix(records)
    
    filtered = []
    for record in records:
        if record.query_name.startswith(query_prefix) and \
           record.target_name.startswith(target_prefix):
            filtered.append(record)
    return filtered, target_prefix


def validate_paf_fasta_consistency(
    paf_records: List[PAFRecord],
    fasta_sequences: Dict[str, str],
    query_prefix: str = "SUPER_"
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate consistency between PAF and FASTA files for chromosomes with given prefix.
    
    Checks that all chromosomes with query_prefix in PAF exist in FASTA and vice versa.
    
    Args:
        paf_records: List of PAF records
        fasta_sequences: Dictionary of FASTA sequences
        query_prefix: Prefix for chromosome names
        
    Returns:
        Tuple of:
        - is_valid: True if all chromosomes match
        - in_paf_not_fasta: List of chromosomes in PAF but not in FASTA
        - in_fasta_not_paf: List of chromosomes in FASTA but not in PAF
    """
    # Get unique chromosome names from PAF (excluding unloc)
    paf_chromosomes = set()
    for record in paf_records:
        if record.query_name.startswith(query_prefix) and "_unloc_" not in record.query_name:
            paf_chromosomes.add(record.query_name)
    
    # Get chromosome names from FASTA (excluding unloc)
    fasta_chromosomes = set()
    for name in fasta_sequences.keys():
        if name.startswith(query_prefix) and "_unloc_" not in name:
            fasta_chromosomes.add(name)
    
    # Find mismatches
    in_paf_not_fasta = sorted(paf_chromosomes - fasta_chromosomes)
    in_fasta_not_paf = sorted(fasta_chromosomes - paf_chromosomes)
    
    is_valid = len(in_paf_not_fasta) == 0 and len(in_fasta_not_paf) == 0
    
    return is_valid, in_paf_not_fasta, in_fasta_not_paf


@dataclass
class ChromosomeMapping:
    """Mapping information for a single chromosome."""
    query_name: str
    query_length: int
    target_name: str
    total_alignment_length: int
    coverage: float
    plus_strand_length: int
    minus_strand_length: int
    needs_reverse_complement: bool
    target_prefix: str = "chr_"  # Reference prefix (chr_ or chr)
    
    @property
    def target_suffix(self) -> str:
        """Extract suffix from target name (e.g., 'chr_5' -> '5', 'chrW' -> 'W')."""
        if self.target_name.startswith(self.target_prefix):
            return self.target_name[len(self.target_prefix):]
        return self.target_name


@dataclass
class UnlocMapping:
    """Mapping information for an unlocalized contig (SUPER_N_unloc_M)."""
    contig_name: str
    parent_chromosome: str  # e.g., "SUPER_5" for "SUPER_5_unloc_1"
    unloc_number: int       # e.g., 1 for "SUPER_5_unloc_1"
    needs_reverse_complement: bool  # Inherited from parent chromosome


@dataclass
class FinalChromosomeAssignment:
    """Final assignment for a chromosome after conflict resolution."""
    original_name: str           # e.g., "SUPER_1"
    new_name: str               # e.g., "SUPER_13" (after renaming)
    new_suffix: str             # e.g., "13" or "W"
    needs_reverse_complement: bool
    is_sex_chromosome: bool
    

# Sex chromosome patterns (letters or letters with numbers like Z1, Z2)
SEX_CHROMOSOME_SUFFIXES = {'W', 'X', 'Y', 'Z'}


def is_sex_chromosome_suffix(suffix: str) -> bool:
    """
    Check if a suffix indicates a sex chromosome.
    
    Sex chromosomes have non-numeric suffixes like W, Z, X, Y, Z1, Z2, etc.
    
    Args:
        suffix: The suffix to check (e.g., "5", "W", "Z1")
        
    Returns:
        True if this is a sex chromosome suffix
    """
    if not suffix:
        return False
    
    # Check if the first character is a known sex chromosome letter
    first_char = suffix[0].upper()
    if first_char in SEX_CHROMOSOME_SUFFIXES:
        # Remaining characters should be digits or empty (W, Z, Z1, Z2, etc.)
        remaining = suffix[1:]
        return remaining == '' or remaining.isdigit()
    
    return False


def extract_chromosome_suffix(name: str, prefix: str = "SUPER_") -> str:
    """
    Extract chromosome suffix from name.
    
    Args:
        name: Chromosome name (e.g., "SUPER_5", "chr_W", "chrZ1", "scaffold_10")
        prefix: Prefix to strip
        
    Returns:
        Suffix string (e.g., "5", "W", "Z1")
    """
    if name.startswith(prefix):
        return name[len(prefix):]
    return name


def is_autosome_suffix(suffix: str) -> bool:
    """
    Check if suffix represents an autosome (purely numeric).
    
    Args:
        suffix: Chromosome suffix
        
    Returns:
        True if this is an autosome suffix (numeric)
    """
    return suffix.isdigit()
    

def is_unloc_contig(name: str) -> bool:
    """
    Check if sequence name is an unlocalized contig.
    
    Unloc contigs have format: SUPER_N_unloc_M
    
    Args:
        name: Sequence name
        
    Returns:
        True if name matches unloc pattern
    """
    return "_unloc_" in name


def parse_unloc_name(name: str) -> Tuple[str, int]:
    """
    Parse unlocalized contig name to extract parent chromosome and unloc number.
    
    Args:
        name: Contig name like "SUPER_5_unloc_1"
        
    Returns:
        Tuple of (parent_chromosome, unloc_number), e.g., ("SUPER_5", 1)
    """
    # SUPER_5_unloc_1 -> ["SUPER_5", "1"]
    parts = name.split("_unloc_")
    parent = parts[0]
    unloc_num = int(parts[1]) if len(parts) > 1 else 0
    return parent, unloc_num


def group_alignments_by_query(records: List[PAFRecord]) -> Dict[str, List[PAFRecord]]:
    """
    Group PAF records by query chromosome name.
    
    Args:
        records: List of PAF records
        
    Returns:
        Dictionary mapping query names to lists of PAF records
    """
    groups = defaultdict(list)
    for record in records:
        groups[record.query_name].append(record)
    return dict(groups)


def calculate_target_alignments(records: List[PAFRecord]) -> Dict[str, Dict[str, int]]:
    """
    Calculate total alignment lengths for each query-target pair.
    
    Args:
        records: List of PAF records for a single query
        
    Returns:
        Dictionary mapping target names to dict with 'total', 'plus', 'minus' lengths
    """
    target_stats = defaultdict(lambda: {'total': 0, 'plus': 0, 'minus': 0})
    
    for record in records:
        alignment_len = record.alignment_block_length
        target_stats[record.target_name]['total'] += alignment_len
        
        if record.strand == '+':
            target_stats[record.target_name]['plus'] += alignment_len
        else:
            target_stats[record.target_name]['minus'] += alignment_len
    
    return dict(target_stats)


def determine_best_target(target_stats: Dict[str, Dict[str, int]]) -> Tuple[str, Dict[str, int]]:
    """
    Determine the best target chromosome based on total alignment length.
    
    Args:
        target_stats: Dictionary from calculate_target_alignments
        
    Returns:
        Tuple of (best_target_name, stats_dict)
    """
    best_target = None
    best_stats = None
    max_length = 0
    
    for target_name, stats in target_stats.items():
        if stats['total'] > max_length:
            max_length = stats['total']
            best_target = target_name
            best_stats = stats
    
    return best_target, best_stats


def build_chromosome_mappings(
    records: List[PAFRecord],
    min_coverage: float = 0.5,
    target_prefix: str = "chr_"
) -> List[ChromosomeMapping]:
    """
    Build chromosome mappings from PAF records.
    
    For each query chromosome:
    1. Calculate alignment lengths to each target
    2. Select best target (maximum alignment length)
    3. Check coverage threshold
    4. Determine orientation based on strand statistics
    
    Args:
        records: Filtered PAF records
        min_coverage: Minimum coverage threshold (0.0-1.0)
        target_prefix: Prefix used in reference chromosome names (chr_ or chr)
        
    Returns:
        List of ChromosomeMapping objects
    """
    # Get query lengths from records
    query_lengths = {}
    for record in records:
        query_lengths[record.query_name] = record.query_length
    
    # Group alignments by query
    groups = group_alignments_by_query(records)
    
    mappings = []
    
    for query_name, query_records in groups.items():
        query_length = query_lengths[query_name]
        
        # Calculate alignment stats for each target
        target_stats = calculate_target_alignments(query_records)
        
        # Find best target
        best_target, best_stats = determine_best_target(target_stats)
        
        if best_target is None:
            continue
        
        # Calculate coverage
        coverage = best_stats['total'] / query_length
        
        # Skip if coverage below threshold
        if coverage < min_coverage:
            print(f"  Warning: {query_name} -> {best_target} coverage {coverage:.2%} below threshold")
            continue
        
        # Determine orientation
        needs_rc = best_stats['minus'] > best_stats['plus']
        
        mapping = ChromosomeMapping(
            query_name=query_name,
            query_length=query_length,
            target_name=best_target,
            total_alignment_length=best_stats['total'],
            coverage=coverage,
            plus_strand_length=best_stats['plus'],
            minus_strand_length=best_stats['minus'],
            needs_reverse_complement=needs_rc,
            target_prefix=target_prefix
        )
        mappings.append(mapping)
    
    return mappings


def build_unloc_mappings(
    sequences: Dict[str, str],
    chromosome_mappings: List[ChromosomeMapping]
) -> List[UnlocMapping]:
    """
    Build mappings for unlocalized contigs based on parent chromosome mappings.
    
    Unloc contigs (SUPER_N_unloc_M) inherit orientation from their parent
    chromosome (SUPER_N). They don't need alignment-based mapping.
    
    Args:
        sequences: Dictionary of all sequences (name -> sequence)
        chromosome_mappings: List of chromosome mappings
        
    Returns:
        List of UnlocMapping objects
    """
    # Build lookup: parent chromosome name -> needs_reverse_complement
    parent_rc_lookup = {m.query_name: m.needs_reverse_complement for m in chromosome_mappings}
    
    unloc_mappings = []
    
    for seq_name in sequences.keys():
        if not is_unloc_contig(seq_name):
            continue
        
        parent_chr, unloc_num = parse_unloc_name(seq_name)
        
        # Get orientation from parent (default to False if parent not mapped)
        needs_rc = parent_rc_lookup.get(parent_chr, False)
        
        unloc_mapping = UnlocMapping(
            contig_name=seq_name,
            parent_chromosome=parent_chr,
            unloc_number=unloc_num,
            needs_reverse_complement=needs_rc
        )
        unloc_mappings.append(unloc_mapping)
    
    return unloc_mappings


def print_unloc_summary(unloc_mappings: List[UnlocMapping]) -> None:
    """
    Print summary of unlocalized contig mappings.
    
    Args:
        unloc_mappings: List of UnlocMapping objects
    """
    if not unloc_mappings:
        return
    
    print(f"\nUnlocalized contigs: {len(unloc_mappings)}")
    print("-" * 60)
    print(f"{'Contig':<25} {'Parent':<15} {'RC?':>5}")
    print("-" * 60)
    
    for m in sorted(unloc_mappings, key=lambda x: (x.parent_chromosome, x.unloc_number)):
        print(f"{m.contig_name:<25} {m.parent_chromosome:<15} "
              f"{'Yes' if m.needs_reverse_complement else 'No':>5}")


def natural_sort_key(name: str, prefix: str = "SUPER_") -> Tuple:
    """
    Generate a sort key for natural sorting of chromosome names.
    
    Sorts autosomes numerically (1, 2, ... 10, 11, ...) 
    and sex chromosomes alphabetically after autosomes.
    
    Args:
        name: Chromosome name (e.g., "SUPER_1", "SUPER_W")
        prefix: Prefix to strip (default "SUPER_")
        
    Returns:
        Tuple for sorting: (is_sex_chr, numeric_or_alpha_key)
    """
    suffix = extract_chromosome_suffix(name, prefix)
    
    if is_sex_chromosome_suffix(suffix):
        # Sex chromosomes come after autosomes, sorted alphabetically
        return (1, suffix)
    else:
        # Autosomes sorted numerically
        try:
            return (0, int(suffix))
        except ValueError:
            return (0, float('inf'))  # Unknown format goes to end


def print_mapping_summary(mappings: List[ChromosomeMapping]) -> None:
    """
    Print summary of chromosome mappings.
    
    Args:
        mappings: List of ChromosomeMapping objects
    """
    print("\nChromosome mapping summary:")
    print("-" * 80)
    print(f"{'Query':<15} {'Target':<10} {'Coverage':>10} {'Plus':>12} {'Minus':>12} {'RC?':>5}")
    print("-" * 80)
    
    for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
        print(f"{m.query_name:<15} {m.target_name:<10} {m.coverage:>9.1%} "
              f"{m.plus_strand_length:>12,} {m.minus_strand_length:>12,} "
              f"{'Yes' if m.needs_reverse_complement else 'No':>5}")


def save_mapping_tsv(
    mappings: List[ChromosomeMapping],
    assignments: List[FinalChromosomeAssignment],
    output_path: Path
) -> None:
    """
    Save chromosome mapping summary to TSV file.
    
    Args:
        mappings: List of ChromosomeMapping objects
        assignments: List of FinalChromosomeAssignment for renamed_to column
        output_path: Path to output TSV file
    """
    # Build lookup: original_name -> new_name
    rename_lookup = {a.original_name: a.new_name for a in assignments}
    
    with open(output_path, 'w') as f:
        # Header
        f.write("query\ttarget\trenamed_to\tquery_length\talignment_length\tcoverage\t"
                "plus_strand\tminus_strand\tneeds_reverse_complement\n")
        
        for m in sorted(mappings, key=lambda x: natural_sort_key(x.target_name, "chr_")):
            renamed_to = rename_lookup.get(m.query_name, m.query_name)
            f.write(f"{m.query_name}\t{m.target_name}\t{renamed_to}\t{m.query_length}\t"
                    f"{m.total_alignment_length}\t{m.coverage:.4f}\t"
                    f"{m.plus_strand_length}\t{m.minus_strand_length}\t"
                    f"{'yes' if m.needs_reverse_complement else 'no'}\n")
    
    print(f"Mapping summary saved to: {output_path}")


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    
    # Handle non-standard bases by keeping them as-is
    result = []
    for base in reversed(seq):
        result.append(complement.get(base, base))
    
    return ''.join(result)


def resolve_chromosome_assignments(
    mappings: List[ChromosomeMapping],
    sequences: Dict[str, str],
    query_prefix: str = "SUPER_",
    output_prefix: str = "SUPER_"
) -> Tuple[List[FinalChromosomeAssignment], Dict[str, bool]]:
    """
    Resolve chromosome assignments with conflict handling and sex chromosome logic.
    
    This function handles:
    1. Sex chromosome detection (W, Z, X, Y, Z1, Z2, etc.)
    2. Conflict resolution when multiple chromosomes map to same target
    3. Autosome -> sex chromosome mapping (reassign to max_autosome + N)
    4. Sex chromosome -> autosome mapping (skip the autosome number)
    5. Unmapped chromosomes keep original names (may cause mapped chr to shift)
    
    Args:
        mappings: List of ChromosomeMapping objects from PAF analysis
        sequences: Dictionary of all sequences (for finding unmapped ones)
        query_prefix: Input chromosome prefix (e.g., "SUPER_", "scaffold_")
        output_prefix: Output chromosome prefix (e.g., "SUPER_", "chr_", "chr", "")
        
    Returns:
        Tuple of:
        - List of FinalChromosomeAssignment for all chromosomes
        - Dictionary mapping original names to needs_reverse_complement flag
    """
    assignments = []
    rc_lookup = {}  # original_name -> needs_rc
    
    # First, identify all chromosomes (excluding unloc contigs)
    all_chr_names = {name for name in sequences.keys() 
                     if name.startswith(query_prefix) and "_unloc_" not in name}
    mapped_names = {m.query_name for m in mappings if "_unloc_" not in m.query_name}
    unmapped_names = all_chr_names - mapped_names
    
    # Separate autosomes and sex chromosomes based on QUERY names
    autosome_mappings = []
    sex_mappings = []
    
    for m in mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_prefix)
        if "_unloc_" in m.query_name:
            continue
            
        if is_sex_chromosome_suffix(query_suffix):
            sex_mappings.append(m)
        else:
            autosome_mappings.append(m)
    
    # Reserved numbers: numbers used by unmapped chromosomes
    # These chromosomes keep their original names, so we can't assign these numbers
    reserved_numbers = set()
    for name in unmapped_names:
        suffix = extract_chromosome_suffix(name, query_prefix)
        if is_autosome_suffix(suffix):
            reserved_numbers.add(int(suffix))
            print(f"  Reserved: number {suffix} (unmapped {name} keeps original name)")
    
    # Track which target numbers are reserved by sex chromosomes in reference
    # (when sex -> autosome, we skip that autosome number)
    skipped_numbers = set()
    
    # Process sex chromosomes first - they keep their original suffix
    for m in sex_mappings:
        query_suffix = extract_chromosome_suffix(m.query_name, query_prefix)
        target_suffix = m.target_suffix
        
        # If sex chromosome maps to autosome target, mark that number as skipped
        if is_autosome_suffix(target_suffix):
            skipped_numbers.add(int(target_suffix))
            print(f"  Note: {m.query_name} (sex chr) -> {m.target_name} (autosome): "
                  f"number {target_suffix} will be skipped")
        
        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{query_suffix}",  # Keep original suffix
            new_suffix=query_suffix,
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=True
        )
        assignments.append(assignment)
        rc_lookup[m.query_name] = m.needs_reverse_complement
    
    # Build mapping: target number -> list of (mapping, alignment_length)
    target_to_autosomes = defaultdict(list)
    autosomes_to_sex_target = []  # Autosomes mapping to sex chromosome targets
    
    for m in autosome_mappings:
        target_suffix = m.target_suffix
        
        if is_sex_chromosome_suffix(target_suffix):
            # Autosome maps to sex chromosome target - will get new number
            autosomes_to_sex_target.append(m)
            print(f"  Note: {m.query_name} (autosome) -> {m.target_name} (sex chr): "
                  f"will be reassigned")
        elif is_autosome_suffix(target_suffix):
            target_num = int(target_suffix)
            target_to_autosomes[target_num].append((m, m.total_alignment_length))
    
    # Resolve conflicts: for each target, pick the one with most alignment
    assigned_numbers = set()
    autosome_assignments = []  # (mapping, assigned_number)
    deferred_autosomes = []  # Mappings that lost in conflict resolution
    
    for target_num, candidates in target_to_autosomes.items():
        # Skip this number if it's reserved by an unmapped chromosome
        if target_num in reserved_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            print(f"  Number {target_num} reserved - deferring: {[m.query_name for m, _ in candidates]}")
            continue
            
        # Skip if blocked by sex chromosome mapping
        if target_num in skipped_numbers:
            deferred_autosomes.extend([m for m, _ in candidates])
            continue
            
        # Sort by alignment length, descending
        candidates.sort(key=lambda x: x[1], reverse=True)
        
        # Winner gets the number
        winner, _ = candidates[0]
        autosome_assignments.append((winner, target_num))
        assigned_numbers.add(target_num)
        
        # Losers are deferred
        for m, _ in candidates[1:]:
            deferred_autosomes.append(m)
            print(f"  Conflict: {m.query_name} lost to {winner.query_name} for number {target_num}")
    
    # Find max assigned autosome number (including reserved)
    all_used_numbers = assigned_numbers | reserved_numbers
    max_autosome = max(all_used_numbers) if all_used_numbers else 0
    
    # Assign deferred autosomes (from conflicts) to next available numbers
    next_available = max_autosome + 1
    unavailable = assigned_numbers | reserved_numbers | skipped_numbers
    
    for m in deferred_autosomes:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        print(f"  Reassigned: {m.query_name} -> {output_prefix}{next_available}")
        next_available += 1
    
    # Assign autosomes that mapped to sex chromosome targets
    for m in autosomes_to_sex_target:
        while next_available in unavailable:
            next_available += 1
        autosome_assignments.append((m, next_available))
        assigned_numbers.add(next_available)
        unavailable.add(next_available)
        print(f"  Reassigned (was sex target): {m.query_name} -> {output_prefix}{next_available}")
        next_available += 1
    
    # Create final assignments for autosomes
    for m, num in autosome_assignments:
        assignment = FinalChromosomeAssignment(
            original_name=m.query_name,
            new_name=f"{output_prefix}{num}",
            new_suffix=str(num),
            needs_reverse_complement=m.needs_reverse_complement,
            is_sex_chromosome=False
        )
        assignments.append(assignment)
        rc_lookup[m.query_name] = m.needs_reverse_complement
    
    # Handle unmapped chromosomes - they keep original suffix with new prefix
    for name in unmapped_names:
        suffix = extract_chromosome_suffix(name, query_prefix)
        assignment = FinalChromosomeAssignment(
            original_name=name,
            new_name=f"{output_prefix}{suffix}",
            new_suffix=suffix,
            needs_reverse_complement=False,
            is_sex_chromosome=is_sex_chromosome_suffix(suffix)
        )
        assignments.append(assignment)
        rc_lookup[name] = False
        print(f"  Unmapped: {name} -> {output_prefix}{suffix} (keeping orientation)")
    
    return assignments, rc_lookup


def sort_assignments_for_output(
    assignments: List[FinalChromosomeAssignment]
) -> List[FinalChromosomeAssignment]:
    """
    Sort assignments for output: autosomes by number, then sex chromosomes alphabetically.
    
    Args:
        assignments: List of FinalChromosomeAssignment
        
    Returns:
        Sorted list
    """
    autosomes = [a for a in assignments if not a.is_sex_chromosome]
    sex_chrs = [a for a in assignments if a.is_sex_chromosome]
    
    # Sort autosomes by numeric suffix
    autosomes.sort(key=lambda x: int(x.new_suffix))
    
    # Sort sex chromosomes alphabetically
    sex_chrs.sort(key=lambda x: x.new_suffix)
    
    return autosomes + sex_chrs


def write_fasta(
    sequences: Dict[str, str],
    assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    rc_lookup: Dict[str, bool],
    output_path: Path,
    output_prefix: str = "SUPER_",
    line_width: int = 60
) -> None:
    """
    Write FASTA file with renamed chromosomes and reverse complement where needed.
    
    Args:
        sequences: Original sequences dictionary
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        rc_lookup: Dictionary mapping original names to needs_reverse_complement
        output_path: Path to output FASTA file
        output_prefix: Prefix for output chromosome names
        line_width: Line width for sequence output (default 60)
    """
    # Build lookup: original_name -> new_name
    name_lookup = {a.original_name: a.new_name for a in assignments}
    
    # Build unloc parent lookup: parent_chr_original -> parent_chr_new_suffix
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}
    
    with open(output_path, 'w') as f:
        # Write chromosomes in sorted order
        for a in assignments:
            seq = sequences.get(a.original_name, '')
            if not seq:
                print(f"  Warning: No sequence found for {a.original_name}")
                continue
            
            # Apply reverse complement if needed
            if a.needs_reverse_complement:
                seq = reverse_complement(seq)
            
            # Write header and sequence
            f.write(f">{a.new_name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')
        
        # Write unloc contigs grouped by parent chromosome
        # Sort by parent, then by unloc number
        unloc_sorted = sorted(unloc_mappings, 
                              key=lambda x: (parent_suffix_lookup.get(x.parent_chromosome, x.parent_chromosome),
                                            x.unloc_number))
        
        for unloc in unloc_sorted:
            seq = sequences.get(unloc.contig_name, '')
            if not seq:
                print(f"  Warning: No sequence found for {unloc.contig_name}")
                continue
            
            # Apply reverse complement if needed (inherited from parent)
            if unloc.needs_reverse_complement:
                seq = reverse_complement(seq)
            
            # Generate new name: <output_prefix><new_parent_suffix>_unloc_<N>
            new_parent_suffix = parent_suffix_lookup.get(unloc.parent_chromosome, 
                                                         extract_chromosome_suffix(unloc.parent_chromosome, output_prefix))
            new_name = f"{output_prefix}{new_parent_suffix}_unloc_{unloc.unloc_number}"
            
            # Write header and sequence
            f.write(f">{new_name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i+line_width] + '\n')
    
    print(f"FASTA written to: {output_path}")


def write_chromosome_list(
    assignments: List[FinalChromosomeAssignment],
    unloc_mappings: List[UnlocMapping],
    output_path: Path,
    output_prefix: str = "SUPER_"
) -> None:
    """
    Write chromosome list CSV file.
    
    Format: name,suffix,yes/no (yes for main chr, no for unloc)
    
    Args:
        assignments: List of FinalChromosomeAssignment (sorted for output)
        unloc_mappings: List of UnlocMapping for unloc contigs
        output_path: Path to output CSV file
        output_prefix: Prefix for output chromosome names
    """
    # Build parent suffix lookup
    parent_suffix_lookup = {a.original_name: a.new_suffix for a in assignments}
    
    with open(output_path, 'w') as f:
        # Write main chromosomes
        for a in assignments:
            f.write(f"{a.new_name},{a.new_suffix},yes\n")
        
        # Write unloc contigs grouped by parent
        unloc_sorted = sorted(unloc_mappings,
                              key=lambda x: (parent_suffix_lookup.get(x.parent_chromosome, x.parent_chromosome),
                                            x.unloc_number))
        
        for unloc in unloc_sorted:
            new_parent_suffix = parent_suffix_lookup.get(unloc.parent_chromosome,
                                                         extract_chromosome_suffix(unloc.parent_chromosome, output_prefix))
            new_name = f"{output_prefix}{new_parent_suffix}_unloc_{unloc.unloc_number}"
            f.write(f"{new_name},{new_parent_suffix},no\n")
    
    print(f"Chromosome list written to: {output_path}")


def main():
    """Main entry point."""
    args = parse_args()
    
    query_prefix = args.query_prefix
    output_prefix = args.output_prefix
    
    # Validate inputs
    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
        
    if not args.paf.exists():
        print(f"Error: PAF file not found: {args.paf}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading FASTA file: {args.fasta}")
    sequences = read_fasta(args.fasta)
    print(f"  Found {len(sequences)} sequences")
    
    print(f"Parsing PAF file: {args.paf}")
    paf_records = parse_paf(args.paf)
    print(f"  Found {len(paf_records)} alignment records")
    
    # Filter records (auto-detect reference prefix)
    filtered_records, ref_prefix = filter_paf_records(paf_records, query_prefix)
    print(f"  Detected reference prefix: '{ref_prefix}'")
    print(f"  After filtering ({query_prefix}* -> {ref_prefix}*): {len(filtered_records)} records")
    
    # Validate PAF/FASTA consistency
    print("\nValidating PAF/FASTA consistency...")
    is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
        paf_records, sequences, query_prefix
    )
    
    if in_paf_not_fasta:
        print(f"  WARNING: Chromosomes in PAF but NOT in FASTA: {', '.join(in_paf_not_fasta)}")
        print(f"           These may indicate wrong PAF or FASTA file!")
    
    if in_fasta_not_paf:
        print(f"  WARNING: Chromosomes in FASTA but NOT in PAF: {', '.join(in_fasta_not_paf)}")
        print(f"           These will keep original suffix (no alignment data).")
    
    if is_valid:
        print(f"  OK: All {query_prefix}* chromosomes match between PAF and FASTA")
    
    # Build chromosome mappings
    print(f"\nBuilding chromosome mappings (min coverage: {args.min_coverage:.0%})...")
    mappings = build_chromosome_mappings(filtered_records, args.min_coverage, ref_prefix)
    print(f"  Successfully mapped {len(mappings)} chromosomes")
    
    # Print mapping summary
    print_mapping_summary(mappings)
    
    # Resolve chromosome assignments with conflict handling
    print("\nResolving chromosome assignments...")
    print(f"  Input prefix: '{query_prefix}' -> Output prefix: '{output_prefix}'")
    assignments, rc_lookup = resolve_chromosome_assignments(
        mappings, sequences, query_prefix, output_prefix
    )
    
    # Build unloc mappings (based on parent chromosome orientation from rc_lookup)
    unloc_mappings = build_unloc_mappings(sequences, mappings)
    if unloc_mappings:
        print(f"  Found {len(unloc_mappings)} unlocalized contigs")
    
    # Update unloc RC based on resolved assignments
    for unloc in unloc_mappings:
        unloc.needs_reverse_complement = rc_lookup.get(unloc.parent_chromosome, False)
    
    # Sort assignments for output
    sorted_assignments = sort_assignments_for_output(assignments)
    
    # Save mapping summary to TSV (after assignments are ready)
    mapping_tsv_path = Path(f"{args.output}.mapping.tsv")
    save_mapping_tsv(mappings, sorted_assignments, mapping_tsv_path)
    
    # Print final assignment summary
    print("\nFinal chromosome assignments:")
    print("-" * 70)
    print(f"{'Original':<20} {'New Name':<20} {'Suffix':<10} {'RC?':>5} {'Sex?':>5}")
    print("-" * 70)
    for a in sorted_assignments:
        print(f"{a.original_name:<20} {a.new_name:<20} {a.new_suffix:<10} "
              f"{'Yes' if a.needs_reverse_complement else 'No':>5} "
              f"{'Yes' if a.is_sex_chromosome else 'No':>5}")
    
    print_unloc_summary(unloc_mappings)
    
    # Generate output files
    fasta_output_path = Path(f"{args.output}.fa")
    csv_output_path = Path(f"{args.output}.chromosome.list.csv")
    
    print("\nWriting output files...")
    write_fasta(sequences, sorted_assignments, unloc_mappings, rc_lookup, 
                fasta_output_path, output_prefix)
    write_chromosome_list(sorted_assignments, unloc_mappings, csv_output_path, output_prefix)
    
    print(f"\nDone! Output files:")
    print(f"  - FASTA: {fasta_output_path}")
    print(f"  - Chromosome list: {csv_output_path}")
    print(f"  - Mapping summary: {mapping_tsv_path}")


if __name__ == "__main__":
    main()
