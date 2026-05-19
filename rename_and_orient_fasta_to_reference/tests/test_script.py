"""Minimal test for rename_and_orient.py script."""
import subprocess
import tempfile
import textwrap
from pathlib import Path


def test_script_runs_on_test_data(paf_file, fasta_file):
    """Test that the script runs successfully on test data."""
    # Create a temporary directory for outputs
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "output"
        output_dir.mkdir()

        # Run the script
        cmd = [
            "python", str(Path(__file__).parent.parent / "rename_and_orient.py"),
            "--fasta", str(fasta_file),
            "--paf", str(paf_file),
            "--output-dir", str(output_dir),
            "--output-prefix", "test",
            "--min-coverage", "0.5"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent.parent)

        # Check that it ran successfully
        assert result.returncode == 0, f"Script failed with stderr: {result.stderr}"

        # Check that output files were created
        assert (output_dir / "test.fa").exists()
        assert (output_dir / "test.chromosome.list.csv").exists()
        assert (output_dir / "test.mapping.tsv").exists()

        # Optionally, check that output FASTA has sequences
        with open(output_dir / "test.fa", 'r') as f:
            content = f.read()
            assert '>' in content  # Has FASTA headers
            assert len(content) > 100  # Has some content


SCRIPT = Path(__file__).parent.parent / "rename_and_orient.py"
HEADER = (
    "query\ttarget\trenamed_to\tquery_length\talignment_length\t"
    "coverage\tplus_strand\tminus_strand\tneeds_reverse_complement\n"
)


def _run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPT.parent)


def test_mapping_table_mode_produces_output(tmp_path):
    """--mapping-table mode renames FASTA using a pre-built TSV, no PAF needed."""
    fasta = tmp_path / "hap2.fa"
    fasta.write_text(">SUPER_9\nATCGATCGATCG\n>SUPER_3\nGCTAGCTAGCTA\n")

    table = tmp_path / "mapping.tsv"
    table.write_text(
        HEADER
        + "SUPER_9\tchr_1\tSUPER_1\t12\t11\t0.92\t2\t9\tyes\n"
        + "SUPER_3\tchr_2\tSUPER_2\t12\t11\t0.92\t9\t2\tno\n"
    )

    out_dir = tmp_path / "out"
    result = _run([
        "python", str(SCRIPT),
        "--fasta", str(fasta),
        "--mapping-table", str(table),
        "--output-dir", str(out_dir),
        "--output-prefix", "result",
    ])

    assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    assert (out_dir / "result.fa").exists()
    assert (out_dir / "result.chromosome.list.csv").exists()

    headers = [l.strip()[1:] for l in (out_dir / "result.fa").read_text().splitlines()
               if l.startswith(">")]
    assert "SUPER_1" in headers
    assert "SUPER_2" in headers


def test_mapping_table_mode_applies_rc(tmp_path):
    """Chromosome with needs_reverse_complement=yes must be RC'd in output FASTA."""
    seq = "AAAAACCCCC"      # RC = GGGGTTTTT
    fasta = tmp_path / "hap2.fa"
    fasta.write_text(f">SUPER_1\n{seq}\n")

    table = tmp_path / "mapping.tsv"
    table.write_text(HEADER + "SUPER_1\tchr_1\tSUPER_1\t10\t9\t0.9\t1\t8\tyes\n")

    out_dir = tmp_path / "out"
    _run(["python", str(SCRIPT), "--fasta", str(fasta),
          "--mapping-table", str(table), "--output-dir", str(out_dir),
          "--output-prefix", "result"])

    content = (out_dir / "result.fa").read_text()
    assert "GGGGTTTTT" in content, f"Expected RC sequence in output, got:\n{content}"


def test_mapping_table_mode_unloc_follow_parent(tmp_path):
    """Unloc contigs unique to haplotype 2 appear right after their parent in output."""
    fasta = tmp_path / "hap2.fa"
    fasta.write_text(
        ">SUPER_5\nAAAAA\n>SUPER_5_unloc_1\nCCCCC\n>SUPER_3\nGGGGG\n"
    )

    table = tmp_path / "mapping.tsv"
    table.write_text(
        HEADER
        + "SUPER_5\tchr_1\tSUPER_1\t5\t5\t1.0\t5\t0\tno\n"
        + "SUPER_3\tchr_2\tSUPER_2\t5\t5\t1.0\t5\t0\tno\n"
    )

    out_dir = tmp_path / "out"
    _run(["python", str(SCRIPT), "--fasta", str(fasta),
          "--mapping-table", str(table), "--output-dir", str(out_dir),
          "--output-prefix", "result"])

    headers = [l.strip()[1:] for l in (out_dir / "result.fa").read_text().splitlines()
               if l.startswith(">")]
    assert headers.index("SUPER_1_unloc_1") == headers.index("SUPER_1") + 1, (
        f"Unloc must immediately follow parent. Got order: {headers}"
    )


def test_mapping_table_and_paf_are_mutually_exclusive(tmp_path):
    """Passing both --paf and --mapping-table must produce an error."""
    fasta = tmp_path / "hap.fa"
    fasta.write_text(">SUPER_1\nAAAA\n")
    dummy = tmp_path / "dummy.tsv"
    dummy.write_text("")

    result = _run([
        "python", str(SCRIPT),
        "--fasta", str(fasta),
        "--paf", str(dummy),
        "--mapping-table", str(dummy),
        "--output-dir", str(tmp_path / "out"),
        "--output-prefix", "x",
    ])
    assert result.returncode != 0
