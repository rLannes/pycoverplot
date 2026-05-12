import Rust_covpyo3
from pycoverplot import (
    Groups,
    get_intervall,
    get_file_path,
    update_group_coverage,
    color_list,
)



def test_coverage_runs_on_tiny_bam(tiny_bam, tiny_interval):
    """A single-BAM group should produce a non-empty coverage array."""
    colors = color_list(["PALETTE_BLUE"], size=1)
    group = Groups(colors=colors, bam_files=[str(tiny_bam)])
    group.group_name = "test"
    group.total_reads = [217]  # skip STAR log lookup

    update_group_coverage(
        [group],
        [tiny_interval],
        lib_scheme="frFirstStrand",
        n_thread=1,
    )
    assert len(group.cover) == 1, "Should have one coverage array per BAM"
