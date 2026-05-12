# pycoverplot
![tests](https://github.com/rLannes/pycoverplot/actions/workflows/ci.yml/badge.svg)


**Fast read-coverage plots from BAM files, straight to publication-ready figures.**

![example coverage plot](asset/ENSMUSG00000039419_ENSMUST00000114641.pdf)

*12 BAM files, a 2.24 Mb gene compressed to a readable 14 kb view, replicate-averaged across 4 groups — plotted in 4.4 seconds(12 cpus HPC).*

pycoverplot reads BAM files directly through a Rust backend. No bigWig intermediate, no separate normalization step, no shell pipeline. Built for RNA-seq but works on any aligned data.

### Why pycoverplot

- **Direct BAM → plot.** Skip the `bamCoverage` → bigWig → `pyGenomeTracks` pipeline entirely.
- **Replicate-aware.** Group BAMs by condition, average automatically, render variance as a confidence band.
- **Intron compression.** Rescale long introns to a fixed fraction of the plot width so short 5′ exons stay readable in megabase-scale genes.
- **Fast by design.** Parallel BAM reading via Rust, optional GTF index caching for repeated runs.
- **Sensible defaults.** RPM-normalized from STAR logs out of the box. Strand-aware. CLI and Python API.

> Early-stage software — the API may change between versions. Pin to a commit hash if you need reproducibility.

 
#### A note on GTF parsing
 
Parsing a full Ensembl or GENCODE GTF is the slowest step in most coverage workflows. pycoverplot can builds a sidecar index (`mygtf.gtf.pbi` + `mygtf.gtf.pbi.bi`) and uses it on every subsequent run, so repeated plots against the same annotation are effectively free at the GTF stage. The index is auto-detected — no extra flag required.



> Warning: Index are pickle file You should not use pickle files you have not created yourself, as there is a security risk in running unknown pickle files. you also may want to make sure they have not be tempered with.

---

## Installation
you need a recent rust compiler (january 2026+)

```bash
# it is not on pip serveur yet
# build a wheel 
git clone https://github.com/rLannes/pycoverplot
cd pycoverplot
pyproject-build --wheel # does the heavy lifting
Successfully built pycoverplot-0.1.0-py3-none-any.whl
pip install pycoverplot-0.1.0-py3-none-any.whl
```

All dependencies are made by me, and are installed automatically. Including the Rust backend (`Rust_covpyo3`) and the GTF parser (`gtf_pyparser`).

you may need to install glibc.
---

## Requirements

- Python ≥ 3.10
- Sorted and indexed BAM files (`.bai` index required alongside each `.bam`)
- A GTF annotation file **or** explicit genomic coordinates

---

## Quick Start

### Command line

Plot coverage of exon (see --exon argument to include intron) for two groups over an annotated gene:

```bash
# index the gtf (run once, considerable speed up)
pycoverplot_gtf --file annotation.gtf


pycoverplot 
    --bam ctrl_rep1.bam ctrl_rep2.bam --color PALETTE_BLUE \
    --bam treat_rep1.bam treat_rep2.bam --color PALETTE_RED \
    --group_name ctrl treatment \
    --bam_dir /path/to/bam/files/ \
    --gtf annotation.gtf --gene_id ENSMUSG00000028494 \
    --out figure.pdf
```
--bam flag define a bam group, you can repeat it to define multiple bam group;
--color argument define the color of a given bam group( either one color or must match the nuber of bam file in a bam group)


### Some option worth knowing:
--exon [exon|intron|intron_partial]: plot only the exon, plot the exon + intron or plot the ewon + compress the intron (usefull for very large intron)
--average plot the average with enveloope (two times the standard deviation)
--smooth average windows smoothing
--thread option (multi cpu)
--gene_id: you can plot a specific transcript using geneid:transcriptid
--color_odd plot every other feature (exon/intorn in differene color) or every even intron in different color

 
## Plot coverage over a custom genomic interval instead of an annotated gene:

```bash
pycoverplot
    --bam ctrl.bam --bam treat.bam \
    --inter chr1,+,1000000,1050000 \
    --out figure.pdf
```

### real example with time:
using the data at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268126 aligned with STAR.
the genes is kl-3  




---

### Python API

The scripting API follows three steps: build your groups, fetch coverage, then plot.

```python
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

    for g in group: # reinitialise the coverage value
        g.cover = []

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
```

---

## Input

### BAM files

BAM files must be sorted and indexed. The `.bai` index file must be present in the same directory as the `.bam` file.

### Genomic region

Two options are available and are mutually exclusive:

**GTF + gene ID** — plot all transcripts of a gene, or restrict to a specific transcript using the `GENE_ID:TRANSCRIPT_ID` syntax. The `gene_id` must match the value in your GTF file exactly (it is database-dependent and differs from the gene symbol).

**Custom interval** — plot any arbitrary genomic region using `--inter CHROM,STRAND,START,END`. Multiple intervals on the same chromosome can be provided and will be concatenated in the plot. must be on same chromosome and same strand!

---

## Normalisation

By default, coverage is normalised to reads per million (RPM) using the uniquely mapped read count read from the STAR `Log.final.out` file expected alongside each BAM file. Normalisation can be disabled with `--NoNormalize`.

If STAR logs are not available, read counts can be provided manually with `--read_count` (CLI) or by setting `group.total_reads` directly (API).

---

## Color options

Colors can be specified per group in three ways:

| Format | Example |
|---|---|
| Built-in palette name | `PALETTE_BLUE`, `PALETTE_RED`, `PALETTE_GREEN`, `PALETTE_ORANGE`, `PALETTE_GUGN`, `PALETTE_BUPL`, `PALETTE_GREY` |
| Matplotlib colormap name | `viridis`, `plasma`, `Blues` |
| Explicit hex colors | `#ff0000 #00ff00` (one per file in the group) |

Each built-in palette provides 5 colors. For groups with more than 5 files, use a colormap or explicit hex colors.

---

## CLI Reference

| Argument | Description |
|---|---|
| `--bam` | One or more BAM files per group. Repeat the flag for additional groups. |
| `--bam_dir` | Base directory for BAM files. One shared directory or one per group. |
| `--group_name` | Legend label for each group, in the same order as `--bam`. |
| `--gtf` | GTF annotation file. Required with `--gene_id`. |
| `--gene_id` | Gene ID(s) to plot. Supports `GENE_ID:TRANSCRIPT_ID` syntax. |
| `--inter` | Explicit interval(s) as `CHROM,STRAND,START,END`. Overrides `--gtf`. |
| `--LibLayout` | Library strandedness. Default: `frFirstStrand`. |
| `--exon` | Intron display mode: `exon`, `intron`, or `intron_partial`. Default: `exon`. |
| `--intron_prop` | Max fraction of plot width for introns (with `intron_partial`). Default: `0.3`. |
| `--smooth` | Sliding window size in bp for coverage smoothing. |
| `--alpha` | Coverage line opacity, 0–1. Default: `1`. |
| `--color` | Color specification per group. |
| `--NoNormalize` | Disable RPM normalisation. |
| `--mapq` | Minimum mapping quality. Default: `13`. |
| `--flag_in` | SAM flag filter: reads to include. Default: `0`. |
| `--flag_out` | SAM flag filter: reads to exclude. Default: `256`. |
| `--thread` | Number of parallel threads. Default: `1`. |
| `--width` | Figure width in inches. Default: `8`. |
| `--height` | Figure height in inches. Default: `5`. |
| `--average` | plot the average for each bam group with envelope |
| `--rasterize` | rasterize the figure |
| `--out_file` | Output file path. Format inferred from extension (`.pdf`, `.png`, `.svg`). |
| `--title` | Plot title. |



### troobleshooting:

#### plot is empty or very few reads, and I am sure that should not append!
check the LibLayout, flag_in, flag_out, parameter,

#### How to include all read not just primary alignment?
use "--flag_out 0 --flag_in 0 --mapq 0" options

#### plot take a long time to open
use the rasterize option


