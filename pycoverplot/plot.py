import argparse
from copy import deepcopy
import logging
import numpy as np
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from multiprocessing import Process
from multiprocessing import Pool
from pathlib import Path
import pickle
import Rust_covpyo3
from gtf_pyparser import Interval, get_intron
from coverage import Coverage
# version april 12 2024

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42



def plot(groups, exon="exon", intron_prop=0.30, N=None, alpha=1,
            width=5, height=8, normalize=True, bg_color=None,
            norm_factor=1_000_000, linewidth=1, color_even=None,
            color_odds=None, title="NoTitle", out=None, return_fig=None):
    """
    Render per-sample read-coverage tracks for one or more groups of BAM
    files and optionally save the figure to disk.

    One coverage line is drawn per BAM file. Lines are coloured according
    to each sample's assigned color and grouped in the legend by their
    ``group_name``. Alternating genomic features (exons or introns) can be
    highlighted with background shading via ``color_even`` / ``color_odds``.

    Parameters
    ----------
    groups : list[Groups]
        Experimental groups to plot, as produced by the ``Groups`` dataclass.
        Each group is iterable and yields per-sample dictionaries containing
        ``"cover"`` (Coverage), ``"color"`` (str), ``"total_reads"`` (int),
        and ``"group_name"`` (str).
    exon : {"exon", "intron", "intron_partial"}, optional
        Coverage display mode passed to ``Coverage.get_cover``. Controls
        whether introns are collapsed, shown at full scale, or compressed.
        Default: ``"exon"``.
    intron_prop : float, optional
        Maximum fraction of the plot width that introns may occupy when
        ``exon="intron_partial"``. Must be between 0 and 1. Default: 0.30.
    N : int or None, optional
        Sliding-window size (in base pairs) for smoothing the coverage
        signal with a uniform convolution. Smoothing is disabled when
        ``None``. Default: ``None``.
    alpha : float, optional
        Opacity of the coverage lines, between 0 (transparent) and 1
        (opaque). Default: 1.
    width : float, optional
        Figure width in inches. Default: 5.
    height : float, optional
        Figure height in inches. Default: 8.
    normalize : bool, optional
        If ``True``, coverage is divided by ``total_reads`` and multiplied
        by ``norm_factor`` to produce reads-per-million (RPM) values.
        Requires ``total_reads`` to be set for every sample. Default: ``True``.
    bg_color : str or None, optional
        Background color of the figure canvas. Accepts any Matplotlib-
        compatible color string. No background color is set when ``None``.
        Default: ``None``.
    norm_factor : int, optional
        Scaling factor applied during normalisation. Default: 1 000 000.
    linewidth : float, optional
        Width of the coverage line in points. Default: 1.
    color_even : str or None, optional
        Background shading color for even-indexed genomic features (0, 2,
        4, …). Accepts any Matplotlib-compatible color string. Shading is
        disabled when ``None``. Default: ``None``.
    color_odds : str or None, optional
        Background shading color for odd-indexed genomic features (1, 3,
        5, …). Accepts any Matplotlib-compatible color string. Shading is
        disabled when ``None``. Default: ``None``.
    title : str, optional
        Title displayed above the plot. Default: ``"NoTitle"``.
    out : str or None, optional
        File path for saving the figure. The file format is inferred from
        the extension (e.g. ``.pdf``, ``.png``, ``.svg``). The figure is
        not saved when ``None``. Default: ``None``.
    return_fig : bool or None, optional
        If truthy, the ``matplotlib.figure.Figure`` object is returned.
        Default: ``None``.

    Returns
    -------
    matplotlib.figure.Figure or None
        The figure object if ``return_fig`` is truthy, otherwise ``None``.

    Notes
    -----
    The interval shading loop uses ``interval`` from the last successfully
    plotted sample. If no samples are plotted successfully (e.g. all covers
    are missing), ``interval`` will be undefined and the shading loop will
    raise a ``NameError``.
    """
    fig = plt.figure(figsize=(width, height))
    if bg_color:
        fig.patch.set_facecolor(bg_color)

    ax = fig.add_axes([0,0,1,1])

    zorder = 1
    for group in groups:

        for this in group:

            try:
                assert this["cover"].cover
            except:
                logging.ERROR("cannot plot cover for file {} no cover! skipping it".format(this["filepath"]))
                continue

            cov, interval = this["cover"].get_cover(exon, intron_prop)   

            if normalize:
                cov =  cov / this["total_reads"] * norm_factor

            if N:
                cov = np.convolve([x for x in cov], np.ones(N)/N, mode='same')

            ax.plot(range(len(cov)), cov, alpha=alpha, color=this["color"], linewidth=linewidth, zorder=zorder)
            zorder += 1
            

    for i, inter in enumerate(interval):

        if i % 2 == 0 and color_even:
            ax.axvspan(inter.start, inter.end, alpha=1, color=color_even, zorder=0)

        elif i % 2 == 1 and color_odds:
            ax.axvspan(inter.start, inter.end, alpha=1, color=color_odds, zorder=0)

    
    legend_name = []
    legend_color = []
    for group in groups:
        legend_name.append(group.group_name)
        legend_color.append(group.colors[len(group.colors) // 2])

    custom_lines = [Line2D([0], [0], color=x, lw=4) for x in legend_color]
    plt.gca().legend(handles=custom_lines, labels=legend_name,  bbox_to_anchor=(1.05, 1), loc='upper left'   )
        
    plt.margins(x=0)
    plt.margins(y=0)



    plt.title(title)
    if out:
        plt.savefig(out, bbox_inches="tight", transparent=True)

    if return_fig:
        return fig
        