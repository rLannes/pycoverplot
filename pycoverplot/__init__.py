"""
pycoverplot
===========

Plot read coverage from BAM files over genomic regions.

The public API exposes the four objects needed for a typical scripting
workflow: define groups, fetch coverage, and plot.

Basic usage
-----------
python
from pathlib import Path
from pycoverplot import Groups, get_intervall, color_list, get_file_path, update_group_coverage, plot, get_reads_fromstar

# --- 1. Define groups ---

ctrl_bams  = get_file_path(["ctrl_rep1.bam", "ctrl_rep2.bam"], bam_dir="/data/bam/")
treat_bams = get_file_path(["treat_rep1.bam", "treat_rep2.bam"], bam_dir="/data/bam/")

ctrl_colors  = color_list(["PALETTE_BLUE"], size=len(ctrl_bams))
treat_colors = color_list(["PALETTE_RED"],  size=len(treat_bams))

ctrl_group  = Groups(colors=ctrl_colors,  bam_files=ctrl_bams)
treat_group = Groups(colors=treat_colors, bam_files=treat_bams)
ctrl_group.group_name  = "ctrl"
treat_group.group_name = "treatment"


groups = [ctrl_group, treat_group]

# Populate read counts from STAR logs (skip if using --NoNormalize)
get_reads_fromstar(groups)

# Optionally set read counts for normalisation (if STAR logs are not available)
# ctrl_group.total_reads  = [12_000_000, 11_500_000]
# treat_group.total_reads = [13_000_000, 12_800_000]


# --- 2. Fetch coverage ---

# Retrieve intervals from a GTF file
target_intervals = get_intervall(
    gtf="flybase.gtf",
    gene_id=["FBgn0267432"],
    inter=None
)

for target_name, target_interval in target_intervals.items():

    update_group_coverage(
        groups,
        target_interval,
        lib_scheme="frFirstStrand",
        n_thread=4,
    )

    # --- 3. Plot ---

    plot(
        groups,
        exon="intron_partial",
        intron_prop=0.3,
        normalize=True,
        norm_factor=1_000_000,
        title="Coverage — " + target_name,
        out="figure.pdf",
        color_even="gainsboro" # hilight the exon
    )

For command-line usage see plot_coverage.py or the README.
"""

__version__ = "0.1.0"

from pycoverplot.coverage_plot import (
    Groups,
    get_intervall,
    update_group_coverage,
    update_group_coverage,
    color_list,
    get_file_path,
    get_reads_fromstar,
)

from pycoverplot.plot import plot

__all__ = [
    "update_group_coverage,"
    "Groups",
    "get_intervall",
    "update_group_coverage",
    "color_list",
    "get_file_path",
    "plot",
    "get_reads_fromstar"
]