"""
pycoverplot
===========

Plot read coverage from BAM files over genomic regions.

The public API exposes the four objects needed for a typical scripting
workflow: define groups, fetch coverage, and plot.

Basic usage
-----------
>>> from pycoverplot import Groups, get_intervall, update_group_coverage, plot
>>>
>>> # 1. Build groups
>>> ctrl  = Groups(colors=ctrl_colors,  bam_files=ctrl_bams)
>>> treat = Groups(colors=treat_colors, bam_files=treat_bams)
>>> groups = [ctrl, treat]
>>>
>>> # 2. Fetch coverage
>>> intervals = get_intervall(gtf="annotation.gtf", gene_id=["ENSMUSG00000028494"], inter=None)
>>> for target_name, target_interval in intervals.items():
...     update_group_coverage(groups, target_interval)
...     plot(groups, title=target_name, out="figure.pdf")

For command-line usage see plot_coverage.py or the README.
"""

__version__ = "0.1.0"

from coverage_plot import (
    Groups,
    get_intervall,
    update_group_coverage,
    update_group_coverage,
    color_list,
    get_file_path,
)

from plot import plot

__all__ = [
    "update_group_coverage,"
    "Groups",
    "get_intervall",
    "update_group_coverage",
    "color_list",
    "get_file_path",
    "plot",
]