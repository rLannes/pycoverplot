# pycoverplot

Plot read coverage from BAM files over genomic regions. Built with RNA-seq in mind but compatible with any aligned data in BAM format.

Coverage is computed with a fast Rust backend. Reads are counted per strand and per base, then rendered as a publication-ready matplotlib figure. Multiple experimental groups can be overlaid on the same plot for direct visual comparison.


Reading BAM files is pretty fast due to the Rust backend, but reading GTF is the longest time sink in most cases. To speed it up, you can use pickle.
Run:
pycoverplot_gtf --file <gtf file> --pkl <pkl file>
This will read the GTF into a pickle file. You can then use this pickle file instead of the GTF. It is significantly faster.

⚠️ Warning: You should not use pickle files you have not created yourself, as there is a security risk in running unknown pickle files. you also may want to make sure they have not be tempered with.
---

## Installation

```bash
# it is not on pip serveur yet
pip install -e /lab/solexa_yamashita/people/Romain/Code/PackageCommonBioinfo/pycoverplot
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
pycoverplot 
    --bam ctrl_rep1.bam ctrl_rep2.bam --color PALETTE_BLUE \
    --bam treat_rep1.bam treat_rep2.bam --color PALETTE_RED \
    --group_name ctrl treatment \
    --bam_dir /path/to/bam/files/ \
    --gtf annotation.gtf --gene_id ENSMUSG00000028494 \
    --out figure.pdf
```

Plot coverage over a custom genomic interval instead of an annotated gene:

```bash
pycoverplot
    --bam ctrl.bam --bam treat.bam \
    --inter chr1,+,1000000,1050000 \
    --out figure.pdf
```

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

    for g in group:
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
| `--out_file` | Output file path. Format inferred from extension (`.pdf`, `.png`, `.svg`). |
| `--title` | Plot title. |