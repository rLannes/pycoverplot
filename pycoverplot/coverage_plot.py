import argparse
from dataclasses import dataclass, field
import textwrap
from copy import deepcopy
from  . import coverage
from  . import plot
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import re
import matplotlib.colors as mcolors
from gtf_pyparser import Interval, parse_gtf


# TODO future NICE TO HAVE:
# 1 - support passing glob for bam name
# 2-
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42



REG_IS_HEXA = re.compile(r"^(#[A-Fa-f0-9]{6}$)|(#[A-Fa-f0-9]{8}$)")
REG_IS_INTERVALL = re.compile(r"^(\w+),([\+|\-]),(\d+)\,(\d+)")
logger = logging.getLogger(__name__)


PALETTE_BLUE = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]
PALETTE_RED    = ["#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15"]
PALETTE_GREEN  = ["#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c"]
PALETTE_ORANGE  = ["#feedde", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603"]
PALETTE_GUGN = ["#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"]
PALETTE_BUPL = ["#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"]
PALETTE_GREY =  ["#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"]

PALETTE_DICT = {
    "PALETTE_BLUE":["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"],
    "PALETTE_RED": ["#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15"],
    "PALETTE_GREEN": ["#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c"],
    "PALETTE_ORANGE": ["#feedde", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603"],
    "PALETTE_GUGN":["#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"],
    "PALETTE_BUPL": ["#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"],
    "PALETTE_GREY":  ["#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"],
}

PALETTE_LIST = [
    ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"],
    ["#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15"],
    ["#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c"],
    ["#feedde", "#fdbe85", "#fd8d3c", "#e6550d", "#a63603"],
    ["#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"],
    ["#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"],
    ["#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"],
]


def get_valid_cmap(name: str):
    """
    Check whether a string is a valid Matplotlib colormap name.

    Parameters
    ----------
    name : str
        The colormap name to validate.

    Returns
    -------
    matplotlib.colors.Colormap or None
        The colormap object if ``name`` is valid, ``None`` otherwise.

    """
    try:
        return plt.get_cmap(name)
    except ValueError:
        return None


# else use a valid cmap (non qualitative ones)
def color_list(this_color, size, PALETTE_DICT=PALETTE_DICT):
    """
    Resolve and validate a color specification into a list of hex color strings.

    Accepts three input forms, attempted in order:

    1. A built-in palette name (key in ``PALETTE_DICT``), resolved to the
       corresponding list of hex strings.
    2. A single Matplotlib colormap name, from which ``size`` evenly spaced
       colors are sampled between 0.2 and 1.0 of the colormap range.
    3. A list of hex color strings (``#RRGGBB`` or ``#RRGGBBAA``), validated
       against ``REG_IS_HEXA``.

    Parameters
    ----------
    this_color : str or list[str]
        Color specification. Either a palette name, a Matplotlib colormap
        name, or a list of hex color strings.
    size : int
        Number of colors to generate. Only used when sampling from a
        Matplotlib colormap.

    Returns
    -------
    list[str]
        List of validated hex color strings.

    Raises
    ------
    AssertionError
        If any resolved color string is not a valid hex color.
    """
    #print(this_color, PALETTE_DICT)
    if isinstance(this_color, str) and this_color in PALETTE_DICT:
        this_color = PALETTE_DICT.get(this_color)
    
    if isinstance(this_color, list) and isinstance(this_color[0], str) and this_color[0] in PALETTE_DICT:
        this_color = PALETTE_DICT.get(this_color[0])
    

    
    if (cmap := get_valid_cmap(this_color[0])):
        cpt = 0.2
        step = 0.8 / size
        this_color = []
        
        for _ in range(size):
            this_color.append(mcolors.to_hex(cmap(cpt)))
            cpt += step
            
    for e in this_color:
        try:
            assert REG_IS_HEXA.fullmatch(e)
        except:
            logger.error("{}: is not a valid hexa decimal color", e)
            raise AssertionError
    return this_color


def get_file_path(bam_files, bam_dir=None):
    """
    Resolve BAM filenames to absolute ``Path`` objects and verify they exist.

    Parameters
    ----------
    bam_files : list[str]
        BAM filenames or relative paths to resolve.
    bam_dir : str or None, optional
        Base directory to prepend to each filename. If ``None``, paths are
        resolved relative to the current working directory.

    Returns
    -------
    list[Path]
        List of resolved absolute ``Path`` objects, one per input file.

    Raises
    ------
    AssertionError
        If any resolved path does not point to an existing file.
    """
    base = Path(bam_dir) if bam_dir else Path("")
    res = []
    for bam in bam_files:
        file = base / bam
        file = file.resolve()
        try:
            assert file.is_file()
        except:
            logger.error("file: {} is not a valid path".format(file))
            raise
        res.append(file)
    return res


def get_intervall(gtf, gene_id, inter):
    """
    Build a dictionary of genomic intervals from either a GTF file or
    explicit interval strings.

    gtf parsing is slow to spead things up you can run pycoverplot_gtf --gtf <gtf file> --pkl <pickle file.pkl>
    this will read the gtf and create a pkl file wich is much faster to read. So if you give this function a pkl it will unpickle it and consider it is the gtf.

    When ``inter`` is provided, each string is parsed with ``REG_IS_INTERVALL``
    and assembled into ``Interval`` objects with ``feature_`` set to
    ``"exon"``. All intervals must be on the same chromosome.

    When ``inter`` is ``None``, a GTF file is parsed and intervals are
    extracted for each gene (or transcript) in ``gene_id``. A transcript can
    be targeted specifically using the ``GENE_ID:TRANSCRIPT_ID`` syntax.

    Parameters
    ----------
    gtf : str or None
        Path to the GTF annotation file. Required when ``inter`` is ``None``.
    gene_id : list[str] or None
        Gene IDs to look up in the GTF file. Each entry may optionally include
        a transcript ID suffix in the form ``"GENE_ID:TRANSCRIPT_ID"``.
        Required when ``inter`` is ``None``.
    inter : list[str] or None
        Explicit interval strings in the format ``"CHROM,STRAND,START,END"``.
        When provided, ``gtf`` and ``gene_id`` are ignored.

    Returns
    -------
    dict[str, list[Interval]]
        Mapping from a label string to its corresponding list of exon
        ``Interval`` objects. The label is ``"customInterval"`` for explicit
        intervals, or ``"GENE_ID:TRANSCRIPT_ID"`` for GTF-derived intervals.

    Raises
    ------
    AssertionError
        If any interval string is malformed, if START ≥ END, or if intervals
        span multiple chromosomes.
    """
    results = {}
    if inter:
        chr_p = None
        sub_res = []
        for e in inter:
            if ( m := REG_IS_INTERVALL.fullmatch(e)):
                
                chr_ = m.group(1)
                if chr_p and chr_ != chr_p:
                    logger.error("chromosome must be the same")
                    raise AssertionError
                chr_p = chr_
                strand = m.group(2)
                start = int(m.group(3))
                end = int(m.group(4))
                try:
                    assert start < end
                except:
                    logger.error("INTERVAL: {} START is GREATER then END".format("e"))
                    raise
                sub_res.append(Interval(chr_, strand, start, end, attr={"feature_": "exon"}))
            else:
                logger.error("INTERVAL: {} missformed".format(e))
                raise AssertionError 
        results["customInterval"] =  sub_res
        return results
    
    if str(gtf).split(".")[-1] == "pkl":
        logging.info("reading gtf from pickle file: {}".format(gtf))
        with open(gtf, "rb") as f:
            gtf_obj = pickle.load(f)
    else:
        gtf_obj = parse_gtf(gtf)

    for gene in gene_id:
        if ':' in gene:
            (gene_id, tr_id) = gene.strip().split(':')
        else: 
            tr_id = None
            gene_id = gene

        if (g := gtf_obj[gene_id]):
            if tr_id:
                results["{}:{}".format(g.gene_id, tr_id)] = g.transcripts[tr_id].exons
            else:
                for tr_id in g.transcripts:
                    results["{}:{}".format(g.gene_id, tr_id)] = g.transcripts[tr_id].exons
        else:
            logger.error("gene_id {} nor found".format(gene_id))
    return results
            
def gtf_to_pkl(gtf, out_):
    dict = parse_gtf(gtf)
    with open(out_, "wb") as fo:
        pickle.dump(obj=dict, file=fo, protocol=pickle.HIGHEST_PROTOCOL)




@dataclass
class Groups:
    """
    Associates a set of BAM files with their display colors and coverage data
    for a single experimental group.

    Attributes
    ----------
    colors : list[str]
        Hex color string for each BAM file, in the same order as
        ``bam_files``.
    bam_files : list[Path]
        Paths to the BAM files belonging to this group.
    group_name : str
        Label used in the plot legend. Defaults to ``"NoGroupName"``.
    cover : list[Coverage]
        Coverage objects collected via ``add_cover``, one per BAM file.
    total_reads : list[int]
        Total mapped read counts for normalisation, one per BAM file,
        collected via ``add_reads_count``.
    """
    colors: list[str]
    bam_files: list[Path]
    group_name = "NoGroupName"
    cover: list[int] = field(default_factory=list)
    total_reads: list[int] =  field(default_factory=list)
    #dico_bam = field(default_factory=dict) # will store everything as bam_file: {color, reads, cover}

    def add_cover(self, cover_dict):
        """
        Populate ``self.cover`` from a BAM-to-Coverage mapping.

        Appends the Coverage object for each BAM file in ``self.bam_files``
        to ``self.cover``, preserving the file order.

        Parameters
        ----------
        cover_dict : dict[Path, Coverage]
            Mapping from BAM file path to its corresponding Coverage object,
            as returned by ``coverage.get_Coverage_intervall``.
        """
        for b in self.bam_files:
            self.cover.append(cover_dict[b])

    def add_reads_count(self, dico_read):
        """
        Populate ``self.total_reads`` from a BAM-to-read-count mapping.

        Appends the read count for each BAM file in ``self.bam_files``
        to ``self.total_reads``, preserving the file order.

        Parameters
        ----------
        dico_read : dict[Path, int]
            Mapping from BAM file path to its total mapped read count.
        """
        for b in self.bam_files:
            self.total_reads.append(dico_read[b])

    def __iter__(self):
        """
        Iterate over per-file metadata dictionaries for this group.

        Yields
        ------
        dict
            A dictionary with the following keys for each BAM file:

            - ``"filepath"`` (Path) : path to the BAM file.
            - ``"color"`` (str) : hex color string for this file.
            - ``"group_name"`` (str) : label of the group.
            - ``"cover"`` (Coverage or None) : coverage object, or ``None``
            if ``add_cover`` has not been called.
            - ``"total_reads"`` (int or None) : read count, or ``None`` if
            ``add_reads_count`` has not been called.
        """
        for i in range(len(self.bam_files)):
            reads = None if not self.total_reads else self.total_reads[i]
            cover = None if not self.cover else self.cover[i]
            yield {"filepath": self.bam_files[i], "color": self.colors[i], "group_name": self.group_name, "cover": cover, "total_reads": reads}


def get_mapped_read(file):
    """
    Parse the number of uniquely mapped reads from a STAR ``Log.final.out``
    file.

    Parameters
    ----------
    file : str or Path
        Path to the STAR log file.

    Returns
    -------
    int or None
        Number of uniquely mapped reads if the relevant line is found,
        ``None`` otherwise.

    Notes
    -----
    Returns ``None`` silently if the ``"Uniquely mapped reads number"`` line
    is not present. The caller is responsible for checking the return value
    before use.
    """
    with open(file) as fin:
        for l in fin:
            l = l.strip()
            if l.startswith("Uniquely mapped reads number"):
                l = l.split()
                return int(l[-1])
            

def _remove_all_suffix(file):
    """
    Strip all suffixes from a ``Path`` object's filename iteratively.

    For example, ``"sample.Aligned.sortedByCoord.out.bam"`` becomes
    ``"sample"``.

    Parameters
    ----------
    file : Path
        File whose compound suffix should be removed.

    Returns
    -------
    str
        The filename stem with all suffixes removed.
    """
    name = file.name
    while Path(name).suffix:
        name = Path(name).stem
    return name  # file


def update_group_coverage(groups, target_interval, lib_scheme="frFirstStrand", n_thread=1,
                          mapq=13, flag_in=0, flag_out=256):

    """
    Fetch per-base coverage for all groups over a single genomic region
    and populate each group's coverage in place.

    Collects all BAM files across every group, computes stranded coverage
    via ``coverage.get_Coverage_intervall``, then calls ``add_cover`` on
    each group so that coverage is accessible through the group iterator.

    Parameters
    ----------
    groups : list[Groups]
        Experimental groups whose BAM files will be processed. Each group's
        ``cover`` list is populated in place after this call.
    target_interval : list[Interval]
        Exon intervals defining the genomic region of interest, as returned
        by ``get_intervall``.
    lib_scheme : str, optional
        Library strandedness scheme passed to the Rust backend.
        Default: ``"frFirstStrand"``.
    n_thread : int, optional
        Number of parallel worker processes. Default: 1.
    mapq : int, optional
        Minimum mapping quality threshold. Default: 13.
    flag_in : int, optional
        SAM flag filter for reads to include. Default: 0.
    flag_out : int, optional
        SAM flag filter for reads to exclude. Default: 256.

    Returns
    -------
    None
        Groups are modified in place.

    """
    
    files = []
    for g in groups:
        files.extend(g.bam_files)
    cover_dict = coverage.get_Coverage_intervall(files, target_interval,
                                                lib_scheme, n_thread=n_thread,
                                                mapq=mapq, flag_in=flag_in, flag_out=flag_out) 
    for g in groups:
        g.cover = []
        g.add_cover(cover_dict)


def get_reads_fromstar(groups):
    """
    Populate each group's ``total_reads`` from STAR ``Log.final.out`` files.

    For each BAM file in every group, attempts to locate the associated STAR
    log file by trying three naming conventions in order:

    1. ``{file.stem}Log.final.out`` — standard STAR output naming.
    2. ``{stem_no_suffix}Log.final.out`` — compound suffixes stripped
       (e.g. ``sample.Aligned.sortedByCoord.out.bam`` → ``sample``).
    3. ``{stem_no_suffix_no_aligned}Log.final.out`` — additionally removes
       the ``"Aligned"`` substring, for lab-specific naming conventions.

    Appends the parsed read count to ``group.total_reads`` for each file,
    in the same order as ``group.bam_files``.

    Parameters
    ----------
    groups : list[Groups]
        Experimental groups whose ``total_reads`` will be populated in place.
        Each group's ``bam_files`` must be set before calling this function.

    Returns
    -------
    None
        Groups are modified in place.

    Raises
    ------
    AssertionError
        If no uniquely mapped read count can be retrieved for any BAM file
        across all three naming conventions.

    Notes
    -----
    This function is tightly coupled to STAR aligner output conventions.
    If your aligner produces differently named log files, set
    ``group.total_reads`` manually instead.
    """
    for group in groups:
                 
        for file in group.bam_files:
            dir_ = file.parent
            tot_reads = None
            path_log = dir_ / (file.stem +  "Log.final.out")
    
            if path_log.exists():
                tot_reads = get_mapped_read(path_log)
            
            else: # custom lab requirment
                path_log = dir_ /  (_remove_all_suffix(file) +  "Log.final.out")
                
                if path_log.exists():
                    tot_reads = get_mapped_read(path_log)
                else:
                    path_log =  dir_ / (_remove_all_suffix(file).replace("Aligned", "") +  "Log.final.out")
                    if path_log.exists():
                        tot_reads = get_mapped_read(path_log)

            if not tot_reads:
                logger.error("could not collect number of uniquely mapped reads")
                raise AssertionError
            group.total_reads.append(tot_reads)

if __name__ == "__main__":

    main()




def main():

    logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M')
    
    logging.getLogger("fontTools").setLevel(logging.ERROR)


    parse = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
        Generate plots from one or more groups of sorted, indexed BAM files.
        Built with RNA-seq in mind but compatible with any aligned data in BAM format.

        Each group of BAM files is defined by a single --bam flag followed by one
        or more BAM files. Repeat --bam to define additional groups (e.g. a control
        group and a treatment group). Groups are plotted together so their coverage
        profiles can be compared side by side. Use --color to assign a color or
        palette to each group (in the same order), and --group_name to label them
        in the legend.

        The tool needs to know which portion of the genome to plot. Either:
            Use --gtf and --gene_id to plot coverage over an annotated gene.
            Note: gene_id is NOT the same as gene name. It is database-dependent —
            for Ensembl, Flybase, and RefSeq it will differ. The gene_id must
            match the value in your GTF file exactly.
        Or:
            Use --inter to plot coverage over any arbitrary genomic interval. (must be same chromosome!)


        Examples
        --------
        # Two groups with named color palettes and group labels:
        plot --bam ctrl_rep1.bam ctrl_rep2.bam --color PALETTE_BLUE \\
             --bam treat_rep1.bam treat_rep2.bam --color PALETTE_RED \\
             --group_name ctrl treatment \\
             --bam_dir /path/to/bam/files/ \\
            --gtf annotation.gtf --gene_id ENSMUSG00000028494

        # Two groups stored in separate directories:
        plot --bam ctrl_rep1.bam ctrl_rep2.bam \\
             --bam treat_rep1.bam treat_rep2.bam \\
             --bam_dir /path/to/ctrl/ /path/to/treatment/ \\
            --gtf annotation.gtf --gene_id ENSMUSG00000028494

        # Plot a specific gene from a GTF annotation:
        plot --bam ctrl.bam --bam treat.bam \\
             --gtf annotation.gtf --gene_id ENSMUSG00000028494 \\

        # Plot a custom genomic interval:
        plot --bam ctrl.bam --bam treat.bam \\
             --inter chr1,+,1000000,1050000
                                    
        # Plot a multiple genomic interval:
        plot --bam ctrl.bam --bam treat.bam \\
             --inter chr1,+,1000000,1050000 chr1,+,1060000,1070000
    """),
    )
    bam_group = parse.add_argument_group("Bam File Input Options")
    bam_group.add_argument("--bam",  action='append', nargs="+", required=True, metavar="FILE",
                help= "One or more BAM files belonging to a single group (e.g. all replicates "
             "of the same condition). Supply multiple files as a space-separated list. "
             "Repeat the flag to define additional groups. "
             "Example: --bam ctrl_rep1.bam ctrl_rep2.bam --bam treat_rep1.bam treat_rep2.bam")
    # may add CDS/UTR latter

    bam_group.add_argument("--bam_dir", "-d", nargs="+", metavar="DIR",
             help="Directory or directories containing the BAM files. When provided, "
             "--bam accepts bare filenames instead of full paths. "
             "Supply a single directory if all groups are in the same location, "
             "or one directory per group (in the same order as --bam) if groups "
             "are in different locations. "
             "Example (shared): --bam_dir /data/bam/ "
             "Example (separate): --bam_dir /data/ctrl/ /data/treatment/")
    
    bam_group.add_argument("--group_name", nargs="+", metavar="NAME",
                        help="Legend label for each BAM group, provided in the same order as --bam. "
             "Example: --group_name ctrl treatment. "
             "Defaults to numeric group indices (Group 1, Group 2, …) when omitted." )
    
    bam_group.add_argument("--read_count", nargs="+", action='append',
                        metavar="int",
                         help="Total read count(s) to use for depth normalization, one value per group "
                        "in the same order as --bam. "
                        "By default the tool automatically retrieves read counts from the "
                        "STAR-generated log files (Log.final.out) expected to be in the same "
                        "directory as each BAM file. Use this option to override that behaviour "
                        "when those files are unavailable or when you want to normalize against "
                        "a specific count (e.g. reads mapping to a particular chromosome). "
                        "Has no effect when --NoNormalize is set.")
    

    parse.add_argument("--color",  action='append', nargs="+", metavar="COLOR",
                        help="Color specification for each BAM group, in the same order as --bam. "
                    "If set for one group it must be set for all groups. "
                    "Accepted values: "
                    "(1) A built-in palette name — each palette provides 6 colors: "
                    "PALETTE_BLUE, PALETTE_RED, PALETTE_GREEN, PALETTE_ORANGE, "
                    "PALETTE_GUGN, PALETTE_BUPL, PALETTE_GREY. "
                    "If you have more than 7 groups you must use options (2) or (3). "
                    "(2) A Matplotlib colormap name (e.g. viridis, plasma, Blues). "
                    "(3) A space-separated list of explicit hex colors, one per file in "
                    "the group (e.g. '#ff0000 #00ff00'). The number of colors must be "
                    "greater than or equal to the number of BAM files in the group. "
                    "Defaults to automatic palette assignment when omitted.")
    

    
    what_to_plot_group = parse.add_argument_group("Intervall Input Options")
    what_to_plot_group.add_argument("--LibLayout", "-s", 
                        metavar="FILE", default="frFirstStrand", 
                         help="Strandedness of the RNA-seq library. This controls which reads are "
                            "counted on each strand. Default: frFirstStrand (the most common "
                            "setting for stranded paired-end RNA-seq, e.g. dUTP/TruSeq libraries). "
                            "Use 'Unstranded' if your library is unstranded. "
                            "Accepted values: Unstranded, PairedUnstranded, "
                            "frFirstStrand, frSecondStrand, fFirstStrand, fSecondStrand, "
                            "ffFirstStrand, ffSecondStrand, rfFirstStrand, rfSecondStrand, "
                            "rFirstStrand, rSecondStrand.")
    
    what_to_plot_group.add_argument("--gtf", "-g", 
                        metavar="FILE",
                         help="GTF annotation file used to define gene and transcript coordinates. "
                        "Must contain a 'gene_id' attribute. Required when using --gene_id. "
                        "Ignored when --inter is provided.")

 
    what_to_plot_group.add_argument("--gene_id", nargs="+",
                        metavar="GENE",
                        help="Gene(s) to plot, given as space-separated gene IDs matching the "
                        "'gene_id' field in the GTF file. By default all transcripts of each "
                        "gene are shown. To restrict to a specific transcript append the "
                        "transcript ID after a colon: GENE_ID:TRANSCRIPT_ID. "
                        "Multiple transcripts for the same gene can be listed separately. "
                        "Example from ensembl GAPDH gene_id is ENSG00000111640, and gene_id of TP53 ENSG00000141510:"
                        " --gene_id ENSG00000111640 ENSG00000141510:ENST00000269305 "
                        "Requires --gtf.")

    what_to_plot_group.add_argument("--inter", nargs="+",
                        metavar="CHROM,STRAND,START,END",
                        help="One or more arbitrary genomic intervals to plot, each specified as "
             "CHROM,STRAND,START,END (comma-separated, no spaces). START must be "
             "less than END (0-based coordinates). Multiple intervals are "
             "concatenated in the plot and any genomic space between them is "
             "ignored, which is useful for visualising discontinuous regions such "
             "as manually defined exons. When this option is set it overrides "
             "CHRM and strand must be the same"
             "--gtf and --gene_id, and intron proportion scaling (--intron_prop) "
             "is not available. "
             "Example: --inter chr1,+,1000000,1050000 chr1,+,1100000,1150000")
    
    plot_option_group = parse.add_argument_group("Plot Options")
    plot_option_group.add_argument("--exon",  default="exon", choices=["exon", "intron", "intron_partial"],        
                                   help="Controls how introns are rendered in the plot. "
                            "'exon': show only exonic regions, introns are collapsed (default). "
                            "'intron': show the full intron length to scale. "
                            "'intron_partial': compress introns to a fraction of their real length "
                            "set by --intron_prop, keeping the plot readable while still indicating "
                            "intron presence.")
    plot_option_group.add_argument("--intron_prop", default=0.3, type=float,
                                    help="Fraction of real intron length to display when --exon intron_partial "
             "is used. Must be between 0 and 1. Default: 0.3 (introns are shown at "
             "30%% of their actual length)." )
    
    plot_option_group.add_argument("--linewidth", default=1, type=float,
                                    help="linewidth of the coverage default 1" )
    plot_option_group.add_argument("--smooth", type=int, metavar="INT",
                help="Window size (in base pairs) for smoothing the coverage signal with a "
             "sliding average. Must be a positive integer. Useful for reducing noise "
             "Smoothing is disabled by default.")
    plot_option_group.add_argument("--alpha", type=float, metavar="FLOAT", default=1,
                                    help="Transparency of the coverage fill, between 0 (fully transparent) and "
             "1 (fully opaque). Lowering this value helps distinguish overlapping "
             "groups. Default: 1.0.")
    plot_option_group.add_argument("--height", type=float, metavar="FLOAT", default=5,
                                    help="Height of the output figure in inches. Default: 5.")
    plot_option_group.add_argument("--width", type=float, metavar="FLOAT", default=8,
                                    help="Width of the output figure in inches. Default: 8.")
    

    plot_option_group.add_argument("--bg_color", type=str, metavar="COLOR", default="whitesmoke",
                                    help="Background color of the plot. Accepts any Matplotlib-compatible color: "
                                    "a named color (e.g. 'white', 'lightgrey'), or a hex RGB value "
                                    "(e.g. '#f5f5f5'). Default: whitesmoke.")
    
    plot_option_group.add_argument("--color_even", type=str, metavar="COLOR",
                                    help="When set, alternating genomic features (exons and introns) are "
             "highlighted by shading their background in two colors. This makes it "
             "easier to visually distinguish where one feature ends and the next "
             "begins. --color_even sets the background color for even-numbered "
             "features (0, 2, 4, …). Accepts a named color or a hex RGB value "
             "(e.g. '#e8f4f8'). Use together with --color_odd. "
             "Highlighting is disabled when both options are omitted.")
    
    plot_option_group.add_argument("--color_odd", type=str, metavar="COLOR",
                                    help="Background color for odd-numbered features (1, 3, 5, …) in the "
             "alternating feature highlight. See --color_even for details. "
             "Accepts a named color or a hex RGB value (e.g. '#fdf3e3'). "
             "Use together with --color_even.")
    
    plot_option_group.add_argument("--title", type=str, metavar="STRING", default="Title",
                                    help="base title for the plot")

    plot_option_group.add_argument("--out_file", type=str, metavar="STRING",
                                    help="out file the extension will set the format myoutput.pdf")
    



    qc_group = parse.add_argument_group("QC Options")
    qc_group.add_argument("--mapq", default=13,
                    help="Minimum mapping quality (MAPQ) score required for a read to be "
                    "included. Reads with a MAPQ below this threshold are discarded. "
                    "Default: 13. For STAR alignments, MAPQ 255 indicates a uniquely "
                    "mapping read; values below 13 typically indicate multi-mappers.")
    qc_group.add_argument("--flag_in", default=0, 
                help="SAM bitwise flag: only reads whose flags match this value are "
             "included. Default: 0 (no flag required). "
             "See https://broadinstitute.github.io/picard/explain-flags.html "
             "for flag definitions.")
    qc_group.add_argument("--flag_out", default=256,
                           help="SAM bitwise flag: reads matching this value are excluded. "
             "Default: 256 (removes non-primary alignments, i.e. keeps only the "
             "primary alignment for each read). "
             "See https://broadinstitute.github.io/picard/explain-flags.html "
             "for flag definitions.")

    parse.add_argument("--NoNormalize", action='store_true', help="Disable depth normalization. Coverage values are displayed as raw "
             "read counts rather than normalized depth.")

    parse.add_argument("--thread", "-t", default=1,
                        metavar="INT",
                         help="Number of parallel threads to use when processing BAM files. "
                    "Increasing this can significantly speed up execution when working "
                    "with many or large BAM files. Default: 1.")
    
    parse.add_argument("--loglevel", "-v", default="INFO",
                        metavar="STR", choices=["INFO", "DEBUG", "WARNING", "ERROR"],
                         help="log level")

    args = parse.parse_args() 

    # validate argument
    level = logging.INFO
    match args.loglevel:
        case "INFO":
            level = logging.INFO
        case "DEBUG":
            level = logging.DEBUG
        case "WARNING":
            level = logging.WARNING
        case "ERROR":
            level = logging.ERROR

    logging.basicConfig(level=level,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M', force=True)
    

    bam_dir = args.bam_dir
    bam_files = args.bam
    colors = args.color
    
    try:
        assert not bam_dir or len(bam_dir) == 1 or len(bam_dir) == len(bam_files)
    except:
        logger.error("if bam_dir is set, then it must be either a single value or match the number of bam groups")
        raise

    try:
        assert not colors or len(colors) == 1 or len(colors) == len(bam_files)
    except:
        logger.error("if color is set, then it must be either a single value or match the number of bam groups")
        raise
    # validate args.read_count shape 

    try:
        assert (args.gtf and args.gene_id) or args.inter
    except:
        logger.error("need either (--gtf AND --gene_id) arguments set OR --inter argument set")
        raise

    #"args.group_name"
    groups_name = []
    if args.group_name: 
        groups_name = args.group_name
    else:
        groups_name = ["Group {}".format(i+1) for i in range(len(bam_files))]
    
    try:
        assert len(groups_name) == len(bam_files)
    except:
        logger.error("if using --group_name must have as many value as --bam flag ")
        raise

    logger.info("reading gtf")
    intervall = get_intervall(args.gtf, args.gene_id, args.inter)
    logger.info("reading gtf: DONE")

    groups = []
    for i, bams in enumerate(bam_files):
        if colors:
            this_color = colors[min(i, len(colors))]
        else:
            this_color = PALETTE_LIST[i]
        
        this_color = color_list(this_color, len(bams))

        this_bam_dir = None
        if bam_dir:
            this_bam_dir = bam_dir[min(i, len(bam_dir)-1)] 

        bam_path = get_file_path(bams, this_bam_dir)


        groups.append(Groups(this_color, bam_path))
        groups[-1].group_name = groups_name[i]

        if not args.NoNormalize and args.read_count:
            this_count = args.read_count[i]
            try:
                assert len(bam_files) == len(this_count)
            except:
                logging.error("if using read count must match number of bam_files: {} bam file {} reads count".format((len(bam_files)), (len(this_count))))


    if not args.NoNormalize:
        get_reads_fromstar(groups)

    files = []
    for g in groups:
        files.extend(g.bam_files)
    
    # for every target do:
    # TODO could optimze to only look for gene once by taking full genes!

    for target_name, target_interval in intervall.items():
        dico_b = {}

        logging.info("reading bam file getting value for {}.".format(target_name))
        cover_dict = coverage.get_Coverage_intervall(files, target_interval,
                                                    lib_scheme=args.LibLayout, n_thread=args.thread,
                                                    mapq=args.mapq, flag_in=args.flag_in, flag_out=args.flag_out) 
        for g in groups:
            g.add_cover(cover_dict)

        logging.info("done reading bam")



        norm = True if not args.NoNormalize else False
        logging.info("plotting")
        out_name = args.out_file.rsplit('.', maxsplit=1)[0] + "_" + target_name.replace(":", "_") + "_" + args.out_file.rsplit('.', maxsplit=1)[1]
        logging.info(f"plotting to {out_name}")
        plot.plot(groups, exon=args.exon, intron_prop=args.intron_prop,
                N=args.smooth, alpha=args.alpha,
                width=args.width, height=args.height, normalize=norm,
                bg_color=args.bg_color, norm_factor=1_000_000,
                linewidth=args.linewidth, color_even=args.color_even,
                color_odds=args.color_odd, title=args.title + " " + target_name,
                out=args.out_file, return_fig=None)



def main_pkl():

    logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M')
    
    logging.getLogger("fontTools").setLevel(logging.ERROR)


    parser = argparse.ArgumentParser()
    parser.add_argument("--file")
    parser.add_argument("--pkl", help="must not exist, delete boferore regenerating it (if you ever need to)")
    args = parser.parse_args()
    pkl = args.pkl
    try:
        assert pkl.split(".")[-1] == "pkl"
    except:
        logging.error("pkl argument must end with the .pkl extension")
        raise AssertionError
    
    try:
        assert not Path(pkl).exists()
    except:
        logging.error("pkl must nort exist erase before regenerating")
        raise AssertionError
    gtf_to_pkl(args.file, args.pkl)

