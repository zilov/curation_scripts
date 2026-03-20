"""Unit tests for individual functions in rename_and_orient.py."""
import pytest
from pathlib import Path
from rename_and_orient import (
    read_fasta,
    parse_paf,
    filter_paf_records,
    validate_paf_fasta_consistency,
    build_chromosome_mappings,
    reverse_complement,
    extract_chromosome_suffix,
    is_sex_chromosome_suffix,
    is_autosome_suffix,
    resolve_chromosome_assignments,
    PAFRecord,
    ChromosomeMapping,
    FinalChromosomeAssignment,
)


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>SUPER_1
ATCGATCG
>SUPER_2
GCTAGCTA
>SUPER_W
NNNNAAAA
"""


@pytest.fixture
def sample_fasta_dict():
    """Expected FASTA dictionary from sample_fasta_content."""
    return {
        "SUPER_1": "ATCGATCG",
        "SUPER_2": "GCTAGCTA",
        "SUPER_W": "NNNNAAAA"
    }


@pytest.fixture
def sample_paf_content():
    """Sample PAF content for testing."""
    return """SUPER_1\t8\t0\t8\t+\tchr_1\t100\t0\t8\t8\t8\t60
SUPER_2\t8\t0\t8\t-\tchr_2\t100\t0\t8\t8\t8\t60
SUPER_W\t8\t0\t8\t+\tchr_W\t100\t0\t8\t8\t8\t60
"""


@pytest.fixture
def sample_paf_records():
    """Expected PAF records from sample_paf_content."""
    return [
        PAFRecord("SUPER_1", 8, 0, 8, "+", "chr_1", 100, 0, 8, 8, 8, 60),
        PAFRecord("SUPER_2", 8, 0, 8, "-", "chr_2", 100, 0, 8, 8, 8, 60),
        PAFRecord("SUPER_W", 8, 0, 8, "+", "chr_W", 100, 0, 8, 8, 8, 60),
    ]


class TestReadFasta:
    """Test read_fasta function."""

    def test_read_fasta_basic(self, tmp_path, sample_fasta_content, sample_fasta_dict):
        """Test basic FASTA reading."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(sample_fasta_content)

        result = read_fasta(fasta_file)
        assert result == sample_fasta_dict

    def test_read_fasta_gzipped(self, tmp_path, sample_fasta_content, sample_fasta_dict):
        """Test reading gzipped FASTA."""
        import gzip
        fasta_file = tmp_path / "test.fa.gz"
        with gzip.open(fasta_file, 'wt') as f:
            f.write(sample_fasta_content)

        result = read_fasta(fasta_file)
        assert result == sample_fasta_dict


class TestParsePaf:
    """Test parse_paf function."""

    def test_parse_paf_basic(self, tmp_path, sample_paf_content, sample_paf_records):
        """Test basic PAF parsing."""
        paf_file = tmp_path / "test.paf"
        paf_file.write_text(sample_paf_content)

        result = parse_paf(paf_file)
        assert len(result) == 3
        assert result[0].query_name == "SUPER_1"
        assert result[0].query_length == 8
        assert result[0].strand == "+"
        assert result[0].target_name == "chr_1"


class TestFilterPafRecords:
    """Test filter_paf_records function."""

    def test_filter_with_auto_detect(self, sample_paf_records):
        """Test filtering with auto-detected target prefix."""
        filtered, prefix = filter_paf_records(sample_paf_records, "SUPER_")
        assert len(filtered) == 3
        assert prefix == "chr_"

    def test_filter_with_custom_prefix(self, sample_paf_records):
        """Test filtering with custom target prefix."""
        # Create records with different target prefix
        records = [
            PAFRecord("SUPER_1", 8, 0, 8, "+", "scaffold_1", 100, 0, 8, 8, 8, 60),
        ]
        filtered, prefix = filter_paf_records(records, "SUPER_", "scaffold_")
        assert len(filtered) == 1
        assert prefix == "scaffold_"


class TestValidatePafFastaConsistency:
    """Test validate_paf_fasta_consistency function."""

    def test_consistent_data(self, sample_paf_records, sample_fasta_dict):
        """Test with consistent PAF and FASTA data."""
        is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
            sample_paf_records, sample_fasta_dict, "SUPER_"
        )
        assert is_valid
        assert len(in_paf_not_fasta) == 0
        assert len(in_fasta_not_paf) == 0

    def test_inconsistent_data(self, sample_paf_records, sample_fasta_dict):
        """Test with inconsistent PAF and FASTA data."""
        # Add a chromosome in PAF but not in FASTA
        extra_record = PAFRecord("SUPER_3", 8, 0, 8, "+", "chr_3", 100, 0, 8, 8, 8, 60)
        paf_records = sample_paf_records + [extra_record]

        is_valid, in_paf_not_fasta, in_fasta_not_paf = validate_paf_fasta_consistency(
            paf_records, sample_fasta_dict, "SUPER_"
        )
        assert not is_valid
        assert "SUPER_3" in in_paf_not_fasta


class TestReverseComplement:
    """Test reverse_complement function."""

    def test_reverse_complement_basic(self):
        """Test basic reverse complement."""
        seq = "ATCG"
        result = reverse_complement(seq)
        assert result == "CGAT"

    def test_reverse_complement_lower_case(self):
        """Test reverse complement with lowercase."""
        seq = "atcg"
        result = reverse_complement(seq)
        assert result == "cgat"

    def test_reverse_complement_with_n(self):
        """Test reverse complement with N bases."""
        seq = "ATCGN"
        result = reverse_complement(seq)
        assert result == "NCGAT"

    def test_reverse_complement_empty(self):
        """Test reverse complement of empty string."""
        result = reverse_complement("")
        assert result == ""


class TestChromosomeSuffixFunctions:
    """Test chromosome suffix related functions."""

    @pytest.mark.parametrize("name,prefix,expected", [
        ("SUPER_1", "SUPER_", "1"),
        ("chr_W", "chr_", "W"),
        ("scaffold_10", "scaffold_", "10"),
        ("chrZ1", "chr", "Z1"),
    ])
    def test_extract_chromosome_suffix(self, name, prefix, expected):
        """Test extract_chromosome_suffix with various inputs."""
        result = extract_chromosome_suffix(name, prefix)
        assert result == expected

    @pytest.mark.parametrize("suffix,expected", [
        ("W", True),
        ("Z", True),
        ("X", True),
        ("Y", True),
        ("Z1", True),
        ("Z2", True),
        ("1", False),
        ("10", False),
        ("ABC", False),
    ])
    def test_is_sex_chromosome_suffix(self, suffix, expected):
        """Test is_sex_chromosome_suffix."""
        result = is_sex_chromosome_suffix(suffix)
        assert result == expected

    @pytest.mark.parametrize("suffix,expected", [
        ("1", True),
        ("10", True),
        ("123", True),
        ("W", False),
        ("Z1", False),
        ("ABC", False),
    ])
    def test_is_autosome_suffix(self, suffix, expected):
        """Test is_autosome_suffix."""
        result = is_autosome_suffix(suffix)
        assert result == expected


class TestBuildChromosomeMappings:
    """Test build_chromosome_mappings function."""

    def test_build_mappings_basic(self, sample_paf_records):
        """Test basic chromosome mapping building."""
        mappings = build_chromosome_mappings(sample_paf_records, 0.5, "chr_")

        assert len(mappings) == 3
        # Check first mapping
        mapping = mappings[0]
        assert isinstance(mapping, ChromosomeMapping)
        assert mapping.query_name == "SUPER_1"
        assert mapping.target_name == "chr_1"
        assert mapping.coverage == 1.0  # 8/8
        assert not mapping.needs_reverse_complement  # strand +

    def test_build_mappings_with_rc(self, sample_paf_records):
        """Test mapping with reverse complement needed."""
        # SUPER_2 has strand -, so should need RC
        mappings = build_chromosome_mappings(sample_paf_records, 0.5, "chr_")

        super2_mapping = next(m for m in mappings if m.query_name == "SUPER_2")
        assert super2_mapping.needs_reverse_complement

    def test_build_mappings_coverage_filter(self):
        """Test coverage filtering."""
        # Create record with low coverage
        low_coverage_record = PAFRecord("SUPER_1", 100, 0, 10, "+", "chr_1", 100, 0, 10, 10, 10, 60)
        mappings = build_chromosome_mappings([low_coverage_record], 0.5, "chr_")

        # Coverage is 10/100 = 0.1, which is below 0.5 threshold
        assert len(mappings) == 0


class TestResolveChromosomeAssignments:
    """Test resolve_chromosome_assignments function."""

    def test_resolve_basic_autosomes(self, sample_fasta_dict):
        """Test basic autosome assignment resolution."""
        mappings = [
            ChromosomeMapping("SUPER_1", 8, "chr_1", 8, 1.0, 8, 0, False, "chr_"),
            ChromosomeMapping("SUPER_2", 8, "chr_2", 8, 1.0, 0, 8, True, "chr_"),
        ]

        assignments, rc_lookup = resolve_chromosome_assignments(
            mappings, sample_fasta_dict, "SUPER_", "SUPER_"
        )

        assert len(assignments) == 3  # SUPER_1, SUPER_2, SUPER_W (unmapped)
        assert assignments[0].original_name == "SUPER_1"
        assert assignments[0].new_name == "SUPER_1"
        assert not assignments[0].needs_reverse_complement

        assert assignments[1].original_name == "SUPER_2"
        assert assignments[1].new_name == "SUPER_2"
        assert assignments[1].needs_reverse_complement

        # SUPER_W is unmapped, keeps original
        super_w = next(a for a in assignments if a.original_name == "SUPER_W")
        assert super_w.new_name == "SUPER_W"
        assert super_w.is_sex_chromosome

    def test_resolve_sex_chromosome_mapping(self, sample_fasta_dict):
        """Test sex chromosome mapping."""
        mappings = [
            ChromosomeMapping("SUPER_W", 8, "chr_W", 8, 1.0, 8, 0, False, "chr_"),
        ]

        assignments, rc_lookup = resolve_chromosome_assignments(
            mappings, sample_fasta_dict, "SUPER_", "SUPER_"
        )

        super_w = next(a for a in assignments if a.original_name == "SUPER_W")
        assert super_w.new_name == "SUPER_W"  # Keeps original suffix
        assert super_w.is_sex_chromosome