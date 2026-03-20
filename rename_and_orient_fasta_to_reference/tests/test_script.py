"""Minimal test for rename_and_orient.py script."""
import subprocess
import tempfile
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