
from copy import deepcopy
import logging
import numpy as np
from multiprocessing import Pool
from pathlib import Path
import Rust_covpyo3
from gtf_pyparser import Interval, get_intron

#from Rust_covpyo3 import get_coverage, get_header
import logging


class Coverage():
    """
    Represents the read coverage of a BAM file over a genomic region.

    Intervals (exons and inferred introns) are sorted by start position and
    stored together. The class provides multiple views of the coverage signal:
    exon-only, full exon+intron, and a compressed intron version where introns
    are downscaled to occupy at most a given proportion of the total plot width.

    Attributes
    ----------
    intervalls : list[Interval]
        Sorted list of genomic intervals (exons and introns).
    leftest_position : int
        Genomic start position of the leftmost interval, used as an offset
        when indexing into the coverage array.
    cover : array-like or None (np array)
        Per-base coverage array for the full region. Must be set externally
        before calling any ``get_*`` method.
    strand : str
        Strand of the feature ('+' or '-'), inferred from the first interval.
    cover_intron_compress : any or None
        Reserved for a compressed intron coverage representation (currently unused).
    intervalls_scaled : any or None
        Reserved for scaled interval coordinates (currently unused).
    intervalls_scaled_intron_compress : any or None
        Reserved for scaled + compressed interval coordinates (currently unused).
    max_reads : int
        Maximum read depth observed; intended for normalisation across plots.
    N : any or None
        Reserved for future use (e.g. total read count).
    """
    def __init__(self, intervalls: list[Interval]) -> None:
        """
        Initialise a Coverage object.

        Sorts the provided intervals by start position, appends inferred
        intron intervals via ``gtf_pyparser.get_intron()``, and records the leftmost
        genomic position.

        Parameters
        ----------
        intervalls : list[Interval]
            Exon intervals for the feature of interest. Intron intervals are
            derived automatically from the gaps between consecutive exons.
        """
        self.intervalls = sorted(intervalls, key=lambda x: x.start)
        self.intervalls.extend(get_intron(self.intervalls))
        self.intervalls = sorted(self.intervalls, key=lambda x: x.start)
        self.leftest_position = self.intervalls[0].start # need to remove this value for position to match the list coverage

        # for plot want to return the corresponding cover and intervalls
        # i.e exons, partial intron.
        self.cover = None # is no array
        self.strand = self.intervalls[0].strand

        
        self.cover_intron_compress = None
        self.intervalls_scaled = None
        self.intervalls_scaled_intron_compress = None
        self.max_reads = 0
        self.N = None
    
    def get_cover(self, type_="exon", intron_partial=None):
        """
        Dispatch coverage retrieval based on the requested plot type.

        Parameters
        ----------
        type_ : str, optional
            Coverage mode. One of:

            - ``"exon"`` : exonic regions only (default).
            - ``"intron"`` : exonic and full intronic regions.
            - ``"intron_partial"`` : exonic regions with introns compressed
            to at most ``intron_partial`` fraction of the total plot width.
        intron_partial : float or None, optional
            Required when ``type_="intron_partial"``. Maximum allowed fraction
            (0–1) of the total plotted width that introns may occupy.

        Returns
        -------
        tuple[list, list[Interval]]
            Coverage chunks and re-indexed intervals, as returned by the
            corresponding ``get_*`` method.

        Raises
        ------
        AssertionError
            If ``type_`` is not one of the accepted values.
        """
        match type_:
            case "exon":
                return self.get_exon_cover()
            case "intron":
                return self.get_intron_cover()
            case "intron_partial":
                return self.get_intron_partial(max_prop=intron_partial)
            case _:
                logging.ERROR("wrong plot type {}, plot type must one of those: exon, intron, intron_partial".format(type_))
                raise AssertionError
    # return cover and intervall for all three
    def get_exon_cover(self):
        """
        Return coverage and intervals for exonic regions only.

        Iterates over the stored intervals, retains only those whose
        ``feature_`` attribute equals ``"exon"``, and extracts the
        corresponding slice from ``self.cover``. Interval coordinates are
        re-indexed to start from 0 in the returned list (i.e. they represent
        positions within the concatenated exon coverage, not genomic positions).

        If the feature is on the minus strand the lists are reversed.

        Returns
        -------
        tuple[list, list[Interval]]
            A 2-tuple of (coverage_chunks, intervals), where:

            - ``coverage_chunks`` is a list of array slices, one per exon.
            - ``intervals`` is the matching list of re-indexed Interval objects.

            Both lists are in 5'→3' order with respect to the transcript.
        """
        this_intervall = []
        this_cover = []

        for inter in self.intervalls:
            if inter.attribute["feature_"] != "exon":
                continue
            # position on this_cover!
            start_ = 0 if not this_cover else this_intervall[-1].end
            
            # position on cover
            cov = self.cover[inter.start - self.leftest_position: inter.end - self.leftest_position]
            
            end_ = start_ + inter.length
            this = Interval(inter.chr, start_,
                      end_, inter.strand,
                      inter.phase, inter.attribute)
            this_intervall.append(this)

            this_cover.append(cov) # append or extend?
        this_cover = np.concatenate(this_cover)
        if this.strand == '-':
            ends = this_intervall[-1].end
            news = [Interval(e.chr, ends - e.end,
                      ends - e.start, e.strand,
                      e.phase, e.attribute)
                    for e in this_intervall]
            return (this_cover[::-1], news)
        else: 
            return (this_cover, this_intervall)
        

    def get_intron_cover(self):
        """
        Return coverage and intervals for exonic and full intronic regions.

        Behaves like ``get_exon_cover`` but retains intervals whose
        ``feature_`` attribute is either ``"exon"`` or ``"intron"``, giving a
        contiguous view of the entire pre-mRNA region. Intron coverage is
        included at full (uncompressed) resolution.

        If the feature is on the minus strand the lists are reversed.

        Returns
        -------
        tuple[list, list[Interval]]
            A 2-tuple of (coverage_chunks, intervals), where:

            - ``coverage_chunks`` is a list of array slices, one per exon or intron.
            - ``intervals`` is the matching list of re-indexed Interval objects.

            Both lists are in 5'→3' order with respect to the transcript.
        """
        this_intervall = []
        this_cover = []

        for inter in self.intervalls:
            if not ((inter.attribute["feature_"] == "exon") or (inter.attribute["feature_"] == "intron")):
                continue
            
            # position on cover
            cov = self.cover[inter.start - self.leftest_position: inter.end - self.leftest_position]
            # position on this_cover!
            start_ = 0 if not this_cover else this_intervall[-1].end
            end_ = start_ + inter.length
            this = Interval(inter.chr, start_,
                      end_, inter.strand,
                      inter.phase, inter.attribute)
            this_intervall.append(this)
            this_cover.append(cov)
        this_cover = np.concatenate(this_cover)
        if this.strand == '-':
            ends = this_intervall[-1].end
            news = [Interval(e.chr, ends - e.end,
                      ends - e.start, e.strand,
                      e.phase, e.attribute)
                    for e in this_intervall]
            return (this_cover[::-1], news)
            #return (this_cover[::-1], this_intervall[::-1])
        else: 
            return (this_cover, this_intervall)
        
    def empyt(self):
        if len(self.cover) > 0:
            return False
        return True
    
    # return cover and intervall for all three
    # intron here is exon + intron (full)
    def get_intron_partial(self, max_prop):
        """
        Return coverage and intervals with introns compressed to a maximum
        proportion of the total plot width.

        Intronic sequence is typically much longer than exonic sequence, which
        makes read-coverage plots hard to read. This method downsamples the
        intron coverage signal by averaging over sliding windows so that
        introns collectively occupy at most ``max_prop`` of the combined
        exon+intron width. Exon coverage is kept at full resolution.

        The compression coefficient is derived as follows::

            final_intron_size = (exon_total * max_prop) / (1 - max_prop)
            coefficient       = final_intron_size / intron_total

        For each intron the window size is::

            window_size = max(1, int(length / (1 + coefficient * length)))

        If the observed intronic proportion is already below ``max_prop`` the
        intron intervals are included unmodified.

        If the feature is on the minus strand the lists are reversed.

        Parameters
        ----------
        max_prop : float
            Maximum allowed fraction (0–1) of the total plotted width that
            introns may occupy. For example, ``0.2`` means introns will be
            compressed until they represent at most 20 % of the plot.

        Returns
        -------
        tuple[list, list[Interval]]
            A 2-tuple of (coverage_chunks, intervals), where:

            - ``coverage_chunks`` is a list of arrays, one per exon or
              compressed intron.
            - ``intervals`` is the matching list of re-indexed Interval objects
              reflecting the compressed coordinates.

            Both lists are in 5'→3' order with respect to the transcript.
        """
        this_intervall = []
        this_cover = []

        # first determine the size exon intron.
        exon = 0
        intron = 0
        for inter in self.intervalls:
            if inter.attribute["feature_"] == "exon":
                exon += inter.length
            elif inter.attribute["feature_"] == "intron":
                intron += inter.length
        
        prop_intron = intron / (exon + intron)
        if prop_intron > max_prop:       
            final_intron_size = (exon*max_prop)/(1-max_prop)
            coefficient = final_intron_size / intron
            
        for inter in self.intervalls:
            if inter.attribute["feature_"] == "exon":

                # position on cover
                cov = self.cover[inter.start - self.leftest_position: inter.end - self.leftest_position]
                # position on this_cover!
                start_ = 0 if not this_cover else this_intervall[-1].end
                end_ = start_ + inter.length
                this = Interval(inter.chr, start_,
                        end_, inter.strand,
                        inter.phase, inter.attribute)
                this_intervall.append(this)
                this_cover.append(cov)

            elif inter.attribute["feature_"] == "intron":
                windows_s = max(1, int(inter.length / (1+(coefficient * inter.length))))
                ##
                cov = self.cover[inter.start - self.leftest_position: inter.end - self.leftest_position]
                if prop_intron > max_prop:    
                    cov = [np.mean(cov[i:i + windows_s]) for i in range(0, inter.length, windows_s)]
                start_ = 0 if not this_cover else this_intervall[-1].end
                end_ = start_ + len(cov)
                this = Interval(inter.chr, start_,
                        end_, inter.strand,
                        inter.phase, inter.attribute)
                this_intervall.append(this)
                this_cover.append(cov)
        this_cover = np.concatenate(this_cover)
        if this.strand == '-':
            ends = this_intervall[-1].end
            news = [Interval(e.chr, ends - e.end,
                      ends - e.start, e.strand,
                      e.phase, e.attribute)
                    for e in this_intervall]
            return (this_cover[::-1], news)
            #return (this_cover[::-1], this_intervall[::-1])
        else: 
            return (this_cover, this_intervall)




    def __add__(self, other):
        """
        Add the per-base coverage arrays of two Coverage objects.

        Both objects must have been built from identical interval lists and
        must have coverage arrays of the same shape.

        Parameters
        ----------
        other : Coverage
            The Coverage object to add to this one.

        Returns
        -------
        Coverage
            A new Coverage object whose ``.cover`` array is the element-wise
            sum of ``self.cover`` and ``other.cover``. Intervals are deep-copied
            from ``self``.

        Raises
        ------
        AssertionError
            If ``self.cover`` and ``other.cover`` have different shapes, or if
            the interval lists differ between the two objects.
        """
        if self.cover.shape != other.cover.shape:
            raise AssertionError("cannot add because the shape does not match {} {}".format(self.cover.shape, other.shape))
        if self.intervalls != other.intervalls:
            raise AssertionError("cannot add because the intervalls are not the same {}\n{}".format(self.intervalls, other.intervalls))

        cover = Coverage(deepcopy(self.intervalls)) 
        cover.cover = [v + other.cover[i] for i , v in enumerate(self.cover)]
        return cover

    
    def normalize(self, factor=1_000_000):
        """
        Normalise the coverage array to reads per million (RPM).

        Divides ``self.cover`` by ``self.max_reads`` and multiplies by
        ``factor``, scaling the signal so that samples with different
        sequencing depths become comparable. Modifies ``self.cover`` in place.

        Parameters
        ----------
        factor : int, optional
            Scaling factor applied after depth normalisation. Defaults to
            1 000 000 (reads per million).

        Notes
        -----
        ``self.max_reads`` must be set before calling this method, either
        directly or via ``get_max_read()``. Calling ``normalize()`` with
        ``max_reads = 0`` will raise a ``ZeroDivisionError``.
        """
        self.cover = self.cover / self.max_reads * factor

    @property
    def max_read(self):
        """
        int : Maximum read depth used for normalisation.
        """
        return self.max_read

    @max_read.setter
    def max_read(self, value):
        self.max_read = value

    def sort_intervalls(self):
        """
        Sort ``self.intervalls`` in place by the ``start`` key.

        Notes
        -----
        This method indexes intervals as dictionaries (``x["start"]``),
        which is inconsistent with the rest of the class, where intervals
        are ``Interval`` objects accessed via ``x.start``. Verify the
        expected type before calling.
        """
        self.intervalls = sorted(self.intervalls, key=lambda x: x["start"])

    def set_cover(self, cover):
        """
        Assign an external coverage array to this object.

        Parameters
        ----------
        cover : array-like
            Per-base coverage values for the full genomic region spanned by
            ``self.intervalls``. Typically a NumPy array of integers returned
            by the Rust backend.
        """
        self.cover = cover

    
def get_max_read(bam_file):
    """
    Parse the total number of mapped reads for a BAM file from an
    associated log file.

    Tries two sources in order:

    1. A STAR ``Log.final.out`` file located by replacing the
       ``"Aligned"`` suffix in ``bam_file`` with ``"Log.final.out"``.
    2. A ``samtools flagstat`` file expected alongside the BAM file,
       with the same stem and a ``.flagstat`` extension.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file. The function derives companion log file
        paths from this string.

    Returns
    -------
    int
        Number of uniquely mapped reads (STAR) or properly paired reads
        (flagstat), depending on which log file is found.

    Raises
    ------
    AssertionError
        If neither a STAR log nor a flagstat file can be located.

    Notes
    -----
    This function is tightly coupled to STAR output naming conventions.
    Users aligning with a different aligner will need to either produce
    a compatible ``flagstat`` file or supply read counts directly by
    setting ``Coverage.max_reads`` manually.
    """

    log_file = bam_file.split("Aligned")[0] + "Log.final.out"
    flagstat = Path(bam_file).parent / (str(Path(bam_file).stem) + ".flagstat")
    if Path(log_file).exists():
        with open(log_file) as f_in:
            for l in f_in:
                spt = [x.strip() for x in l.strip().split("|")]
                if spt[0] == "Uniquely mapped reads number":
                    return int(spt[1])
    elif Path(flagstat).exists():
        with open(flagstat) as f_in:
            for l in f_in:
                if "properly paired" in l:
                    return int(l.split()[0])
    else:
        raise AssertionError("no file found to retirve max reads, I tried: {} and {}".format(log_file, flagstat))



# main public function API entry point
def get_Coverage_intervall(bam_files, intervalls, lib_scheme, n_thread=6, mapq=13, flag_in=256, flag_out=0):
    """
    Compute stranded per-base coverage over a genomic region for multiple
    BAM files in parallel.

    Determines the full span of the provided intervals, then dispatches
    one BAM parsing job per file using a multiprocessing pool. Each job
    calls the Rust backend via ``parse_bam_mp``. Results are wrapped in
    ``Coverage`` objects keyed by BAM file path.

    Parameters
    ----------
    bam_files : list[str or Path]
        Paths to the BAM files to process.
    intervalls : list[Interval]
        Exon intervals defining the region of interest. The genomic span
        is taken as the union from the leftmost start to the rightmost end.
    lib_scheme : str
        Library strandedness scheme passed to the Rust backend
        (e.g. ``"fr-firststrand"``).
    n_thread : int, optional
        Number of parallel worker processes. Defaults to 6.
    mapq : int, optional
        Minimum mapping quality threshold. Defaults to 13.
    flag_in : int, optional
        SAM flag filter for reads to include. Defaults to 256.
    flag_out : int, optional
        SAM flag filter for reads to exclude. Defaults to 0.

    Returns
    -------
    dict[str, Coverage]
        Mapping from BAM file path (as passed in ``bam_files``) to the
        corresponding ``Coverage`` object with ``.cover`` set.
    """

    intervalls = sorted(intervalls, key = lambda x: x.start)
    start = intervalls[0].start
    end = intervalls[-1].end
    range_ = {
                "strand": intervalls[0].strand,
                "chr": intervalls[0].chr
             }
    if start < end:
        range_["start"] = start
        range_["end"] = end
    else:
        range_["start"] = end
        range_["end"] = start
    
    n = range_["end"] - range_["start"]
    logging.info("getting region from {} to {} on the {}. length:{}".format(range_["start"], range_["end"], range_["chr"], n))

    args = [ (str(x), lib_scheme, range_, mapq, flag_in, flag_out)  for x in bam_files ]
    logging.info("launching pool process with {} threads".format(n_thread))
    with Pool(processes=n_thread) as pool: 
        cover_list = pool.map(parse_bam_mp, args)
    

    results = {}
    for i, bam in enumerate(bam_files):
        cover = Coverage(deepcopy(intervalls))
        cover.set_cover(np.array(cover_list[i]))
        results[bam] = cover

    return results



# publicly exposed for advaned user
def get_cover_for_a_bam(intervalls, bam, lib, mapq, flag_in=0, flag_out=0):
    """
    Compute per-base coverage for a single BAM file over a list of intervals.

    Calls the Rust backend (``Rust_covpyo3.get_coverage_algo2``) once per
    interval. This is the single-file, non-parallelised counterpart to
    ``get_Coverage_intervall``, intended for advanced users who manage
    parallelisation themselves.

    Parameters
    ----------
    intervalls : list[dict]
        List of interval dictionaries, each containing the keys
        ``"start"`` (int), ``"end"`` (int), ``"chr"`` (str), and
        ``"strand"`` (str).
    bam : str
        Path to the BAM file.
    lib : str
        Library strandedness scheme passed to the Rust backend.
    mapq : int
        Minimum mapping quality threshold.
    flag_in : int, optional
        SAM flag filter for reads to include. Defaults to 0.
    flag_out : int, optional
        SAM flag filter for reads to exclude. Defaults to 0.

    Returns
    -------
    list
        List of per-base coverage arrays, one per input interval.

    Notes
    -----
    There is a scoping bug in the current implementation: ``cover`` is
    overwritten on each loop iteration but only the last interval's result
    is appended to ``coverage_rust``. Verify this matches the intended
    behaviour before use.
    """
    coverage_rust  = []
    for e in intervalls:
        cover = Rust_covpyo3.get_coverage_algo2(e["start"], e["end"],
                                           e["chr"], e["strand"],
                                           bam, lib, mapq, flag_in, flag_out)

    coverage_rust.append(cover)

    return coverage_rust

# publicly exposed for advaned user
def parse_bam_mp(args):
    """
    Compute stranded per-base coverage for a single BAM file over a
    contiguous genomic region.

    Designed to be called by ``multiprocessing.Pool.map`` inside
    ``get_Coverage_intervall``. Unpacks its arguments from a single tuple
    to comply with the ``pool.map`` single-argument contract.

    Parameters
    ----------
    args : tuple
        A 6-element tuple of ``(bam, lib, intervalls, mapq, flag_in, flag_out)``:

        - ``bam`` (str) : path to the BAM file.
        - ``lib`` (str) : library strandedness scheme.
        - ``intervalls`` (dict) : region dictionary with keys ``"start"``,
          ``"end"``, ``"chr"``, and ``"strand"``.
        - ``mapq`` (int) : minimum mapping quality threshold.
        - ``flag_in`` (int) : SAM flag filter for reads to include.
        - ``flag_out`` (int) : SAM flag filter for reads to exclude.

    Returns
    -------
    list[int]
        Per-base coverage array for the specified region, as returned by
        the Rust backend.
    """
    (bam, lib,  intervalls, mapq, flag_in, flag_out) = args
    # rust API
    # get_coverage_algo2(start:i64, end:i64, chrom: String, strand: String,
    # bam_path: String, lib: String, mapq_thr: u8, flag_in: u16, flag_exclude: u16
    results =  Rust_covpyo3.get_coverage_algo2(intervalls["start"], intervalls["end"], intervalls["chr"], intervalls["strand"],
                                                bam, lib, mapq, flag_in, flag_out)
    return results


def merge_duplicates(list):
    """
    Merge a list of Coverage objects that represent the same genomic region
    (e.g. technical replicates or sequencing lanes) into a single object.

    Coverage arrays are summed element-wise and ``max_reads`` values are
    accumulated. All objects in the list must share identical interval
    definitions and coverage array shapes.

    Parameters
    ----------
    list : list[Coverage]
        Coverage objects to merge. All must have the same ``.intervalls``
        and ``.cover.shape``.

    Returns
    -------
    Coverage
        A new Coverage object with summed coverage and accumulated
        ``max_reads``, built from the intervals of the first element.

    Raises
    ------
    AssertionError
        If any object in the list differs in coverage shape or interval
        definition from the first element.

    Notes
    -----
    Shadowing the built-in ``list`` as a parameter name will suppress
    type hints and IDE autocompletion for that name inside the function.
    Consider renaming the parameter (e.g. ``coverage_list``).
    """
    cover = Coverage(deepcopy(list[0].intervalls)) 
    cover.strand = list[0].strand
    cover.sort_intervalls()
    cover.make_intervall_scaled()

    if not all([x for x in list[1:] if x.cover.shape == list[0].cover.shape]):
        raise AssertionError("cannot add because the shape does not match {} {}".format(self.cover.shape, other.shape))
    if not all([x for x in list[1:] if x.intervalls == list[0].intervalls]):
        raise AssertionError("cannot add because the intervalls are not the same {}\n{}".format(self.intervalls, other.intervalls))
    
    for e in list:
        if cover.cover is None:
            cover.cover = np.array(e.cover)
            cover.max_reads = e.max_reads
        else:
            cover.cover = cover.cover + np.array(e.cover)
            cover.max_reads += e.max_reads

    return cover
    

def merge_duplicate(dico_bam_cover):
    """
    Merge Coverage objects from multiple BAM files that belong to the same
    biological sample, identified by parsing the BAM file stem.

    Groups BAM files by a ``genotype_sample`` key extracted via
    ``reg_get_sample_lanes`` (a module-level regex). Coverage arrays and
    ``max_reads`` values are summed within each group.

    Parameters
    ----------
    dico_bam_cover : dict[str, Coverage]
        Mapping from BAM file path to its corresponding Coverage object,
        as returned by ``get_Coverage_intervall``.

    Returns
    -------
    dict[str, Coverage]
        Mapping from ``"genotype_sample"`` identifier to the merged
        Coverage object.

    Notes
    -----
    This function depends on a module-level regex ``reg_get_sample_lanes``
    that must define three capture groups: genotype, sample, and lane.
    This regex is not defined in this file; ensure it is available in the
    module namespace before calling this function.
    """

    results = {}
    for bam, cover_reads in dico_bam_cover.items():

        m = reg_get_sample_lanes.match(Path(bam).stem)
        genotype = m.group(1)
        sample = m.group(2)
        lane = m.group(3)

        id_ = "{}_{}".format(genotype, sample)

        if id_ not in results:
            results[id_] = cover_reads
        else:
            cover = results[id_]
            cover.cover = cover.cover + np.array(cover_reads.cover)
            cover.max_reads += cover_reads.max_reads
            results[id_] = cover
        
    return results