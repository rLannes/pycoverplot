"""End-to-end CLI tests."""
import subprocess
import sys
import time

def test_cli_help_runs():
    """The --help flag should work and exit cleanly."""
    result = subprocess.run(
        ["pycoverplot", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0




def test_cli_produces_output(tiny_bam, tiny_gtf, tmp_path):
    """A basic CLI invocation should produce an output file."""
    out_file = tmp_path / "out.pdf"

    result = subprocess.run(
        [
            "pycoverplot",
            "--bam", str(tiny_bam),
            "--gtf", str(tiny_gtf),
            "--gene_id", str("ENSMUSG00000029177:ENSMUST00000144742"),
            "--group_name", "test",
            "--color", "PALETTE_BLUE",
            "--read_count 217",
            "--NoNormalize",
            "--out_file", str(out_file),
        ],
        capture_output=True,
        text=True,
    )

    out_file = tmp_path / "out_ENSMUSG00000029177_ENSMUST00000144742.pdf"
    assert result.returncode == 0, f"CLI failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    assert out_file.exists()
    assert out_file.stat().st_size > 0