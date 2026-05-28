# Nanopore Analysis Pipeline - Tutorial

For installation instructions, troubleshooting and WSL/Docker setup, see the **README**.

---

## Table of Contents
 
1. [Full Pipeline Walkthrough](#1-full-pipeline-walkthrough)
2. [Basecalling](#2-basecalling)
3. [Alignment](#3-alignment)
4. [Methylation Summary](#4-methylation-summary)
5. [Deconvolution](#5-deconvolution)
6. [BAM Utilities](#6-bam-utilities)

---

## 1. Full Pipeline Walkthrough

This section describes how a complete analysis run fits together - what each step does, what it takes as input, and what it produces. Each step is covered in detail in its own section below.

The pipeline takes raw Nanopore output and produces a cell-type deconvolution result in four main steps:

```
POD5 files
    ↓  Basecalling (Dorado)
Unaligned BAM
    ↓  Alignment (minimap2 + samtools)
Aligned, sorted, indexed BAM
    ↓  Methylation Summary (modkit) — optional for UXM
BED file
    ↓  Deconvolution (UXM or NNLS)
Cell-type proportions CSV
```

> **Note on runtime:** Basecalling, alignment, and methylation summary all take significant time on a real dataset. This is normal - the pipeline is designed for whole-genome sequencing data. For learning about individual steps, the sections below each explain what you need and how to run them independently without having to run the full pipeline first.

---

### The Sample Dataset

The individual step sections use publicly available data from the Genome in a Bottle 2025.01 release from Oxford Nanopore's open data repository: the HG002 sample, sequenced on a PromethION with R10.4.1 chemistry with 5mCG methylation detection. Resources for this dataset can be downloaded with:

```bash
resource_download sample_data   # pre-basecalled BAM
resource_download genome        # hg38 reference genome
resource_download atlas         # UXM methylation atlas
```

---

### Running the Full Pipeline

Once you're familiar with the individual steps and have your data resources in place, you can run everything in sequence with:

```bash
nanopore_analysis run
```

This executes whichever steps are enabled in `config.yaml` in order. Each step's output becomes the next step's input automatically.

For your own data, the typical `config.yaml` setup is:

```yaml
pipeline_control:
  run_steps:
    basecalling: true
    align: true
    methylation_summary: true  # set false if using UXM only
    analysis: true
```

Progress for each step is logged to the terminal and written in full to a timestamped file in `logs/`.

### Expected Output Structure

After a full run, your directory will look something like this:

```
data/
├── basecalled_output/
│   └── my_experiment.bam
├── alignment/
│   ├── my_experiment.aligned.sorted.bam
│   ├── my_experiment.aligned.sorted.bam.bai
│   └── qc/
│       └── my_experiment.flagstat.txt
├── methylation/
│   └── methylation.bed
├── processed/
│   └── my_experiment.aligned.sorted.pat.gz
└── atlas/
    └── UXM_atlas.tsv
results/
└── deconvolution/
    └── deconvoluted_output.csv
logs/
└── 2024-03-15_10-30-00_pipeline_run.log
```

--- 

## 2. Basecalling

Basecalling is the process of converting the raw electrical signal recorded by the Nanopore device into DNA sequences. The tool that does this is **Dorado**, which is downloaded automatically during installation.

The output is an unaligned BAM file - a standard format for storing sequencing reads, in this case without any positional information because we haven't mapped them to a genome yet.

### What you need

- A directory of `.pod5` files (the raw output from your Nanopore device)
- Dorado installed (handled by `python install.py conda all`)

If your sequencer produced `.fast5` files instead of `.pod5`, convert them first using the `pod5` tool included in the conda environment:

```bash
pod5 convert fast5 data/fast5_input/*.fast5 --output data/pod5/
```

### Choosing a Model

Dorado uses a neural network model to interpret the raw signal. The two things to specify are:

**Model speed** - a tradeoff between accuracy and compute time:
- `hac` (high accuracy) - good accuracy, reasonable speed. Right for most runs.
- `sup` (super accuracy) - better accuracy, significantly slower. Worth it if you're not in a rush or have a good GPU.

**Modification model** - tells Dorado to detect methylation at the same time as basecalling. For CpG methylation use `5mCG_5hmCG`.

Set these in `config.yaml`:

```yaml
pipeline_steps:
  basecalling:
    params:
      method: "complex"
      complex_settings:
        model_speed: "hac"
        basecalling_modifications: "5mCG_5hmCG"
    paths:
      pod5_input_path: "data/pod5"
      basecalled_bam_name: "my_experiment.bam"
```

> **What is `method: complex`?** The pipeline has two ways to specify a model. `complex` lets you describe the model by its properties and the pipeline builds the full model name for you. `explicit` lets you name a specific model version directly. For most purposes `complex` is fine.

### Running Basecalling

```bash
nanopore_analysis basecalling run
```

Or override the input/output without modifying the config:

```bash
nanopore_analysis basecalling run \
  --input-file data/pod5/ \
  --output-dir data/basecalled_output/my_experiment.bam
```

### If your run used barcodes (demultiplexing)

Barcoding lets you sequence multiple samples in a single run by tagging each sample's DNA with a unique sequence. If your run used barcodes, the basecalled BAM will contain reads from all samples mixed together and you need to separate them. This is called **demultiplexing**.

Set `demultiplex: true` and provide your kit name:

```yaml
pipeline_steps:
  basecalling:
    params:
      demultiplex: true
      complex_settings:
        kit_name: "SQK-NBD114-24"
```

Then run as normal - basecalling and demultiplexing happen in one step, producing a separate BAM per barcode in `data/demultiplexed/`.

If you already have a basecalled BAM and just want to demultiplex it:

```bash
nanopore_analysis basecalling demux \
  --input-file data/basecalled_output/my_experiment.bam \
  --output-dir data/demultiplexed/
```

### Checking Your Output

A quick sanity check after basecalling is to look at the read length distribution:

```bash
summarise_lengths \
  data/basecalled_output/my_experiment.bam \
  --output-dir data/qc/
```

A good run shows a broad spread of lengths, typically peaking somewhere between 5-20 kb depending on the library prep. If almost all reads are below 1 kb, something likely went wrong during library preparation rather than basecalling. This of course doesn't apply if there is an experimental reason for the length distribution being different.

---

## 3. Alignment

Alignment is the process of mapping each basecalled read to its position in the reference genome. The pipeline uses **minimap2** for alignment and **samtools** for sorting and indexing the output. The result is a sorted, indexed BAM file - the standard input for all downstream steps.

--- 

### What you need

- A basecalled unaligned BAM (from the basecalling step)
- A reference genome FASTA file in your `reference_genomes/` directory

If you don't have a reference genome yet, download it first:

```bash
resource_download genome
```

> **Install the genome via wgbstools if possible.** The deconvolution step (UXM) relies on several files that wgbstools generates alongside the genome when you use `wgbstools init_genome`.
> If the genome is installed without these files (e.g. directly downloaded as a FASTA), alignment will still work but the wgbstools `bam2pat` and `uxm deconv` steps won't function correctly. The `resource_download genome` command handles this automatically, but if you're providing your own FASTA be sure to run `wgbstools init_genome` on it first.

The pipeline searches `reference_genomes/` for your genome using the `genome_id` set in `config.yaml`. It looks for the file in several common layouts, for example:

```
reference_genomes/hg38.fa
reference_genomes/hg38/genome.fa
reference_genomes/hg38/hg38.fa
```

If your FASTA is somewhere else or has an unusual name, use `custom_fasta_reference` instead of `genome_id`:

```yaml
pipeline_steps:
  align:
    paths:
      custom_fasta_reference: "reference_genomes/my_genome.fa"
```

---

### The minimap2 Index

Before alignment, minimap2 needs an index of the reference genome (a `.mmi` file). If you installed the genome via wgbstools, this index is already generated as part of `init_genome`and the pipeline will find it automatically. If it isn't present, the pipeline will automatically build it the first time you run alignment and saved alongside the FASTA:

```
reference_genomes/hg38/genome.mmi
```

Building the index takes a few minutes but only needs to be done once - subsequent runs will find it and skip straight to alignment.

---

### Running Alignment

```bash
nanopore_analysis align run
```

Or override the input and output without touching the config:

```bash
nanopore_analysis align run \
  --input-file data/basecalled_output/my_experiment.bam \
  --output-dir data/alignment/my_experiment.aligned.sorted.bam
```

This will:
1. Build the minimap2 index if it doesn't already exist
2. Align the reads to the reference genome
3. Sort the output with samtools
4. Index the sorted BAM

A flagstat QC report is also written automatically to `data/alignment/qc`.

---

### Aligning Multiple BAM files

If you demultiplexed your data, you will have a directory of BAM files rather than a single file. Pass the directory as the input and the pipeline will align every BAM it finds:

```bash
nanopore_analysis align run \
  --input-file data/demultiplexed/ \ 
  --output-dir data/alignment
```

Each file is saved as `<stem>.aligned.bam` inside the output directory.

---

### Reading the Flagstat Report

After alignment, a flagstat report is saved to `data/alignment/qc/`. This is a quick summary of alignment quality. The key line to check is the mapped percentage:

```
95.42% mapped
```

For a good whole-genome sequencing run you would typically expect >90% of reads to map. A low mapping rate could indicate a mismatch between the sequencing chemistry and the basecalling model, or the wrong reference genome (i.e. hg19 instead of hg38).

---

### Running QC on an Existing BAM

If you already have an aligned BAM and just want to regenerate the QC report:

```bash
nanopore_analysis align qc \
  --input-file data/alignment/my_experiment.aligned.sorted.bam
```

---

## 4. Methylation Summary

The methylation summary step converts the methylation tags embedded in the aligned BAM into a tabular format that the deconvolution step can read. The tool that does this is **modkit**, which produces a BED file with one row per CpG site across the genome.

> **This step is optional for deconvolution.** The UXM algorithm works directly from the aligned BAM via wgbstools `bam2pat`, bypassing the BED file entirely. The NNLS algorithm does use the BED file, but if you are only running UXM you can skip this step completely by setting `methylation_summary: false` in `config.yaml`. It is included in the tutorial for completeness and because it can be useful for quality checking your data before deconvolution.

---

### What You Need

- A sorted, indexed aligned BAM (from the alignment step)
- Optionally, a reference genome FASTA - strongly recommended

### Running the Methylation Summary

```bash
nanopore_analysis methylation run
```

Or with explicit paths

```bash
nanopore_analysis methylation run --input-file data/alignment/my_experiment.aligned.sorted.bam --output-dir data/methylation/methylation.bed --ref reference_genomes/hg38/genome.fa
```

Providing `--ref` enables CpG-only mode (`--cpg` in modkit), which restricts output to canonical CpG sites. This is significantly faster than the default and produces cleaner data for deconvolution. It is strongly recommended.

Progress is shown as a live updating status line on the terminal. This is the slowest step in the pipeline on a full BAM - expect it to run for some time.

---

### Understanding the Output

The BED file has one row per CpG site. The columns most relevant to downstream analysis are:

| Column | Content |
|---|---|
| 1 | Chromosome |
| 2 | Start position |
| 4 | Modification type (`m` = 5mCG) |
| 11 | Beta value (0–100, percent methylation) |

A beta value of 0 means fully unmethylated at that site, 100 means fully methylated. Most CpG sites in a human sample will sit at one extreme or the other - a bimodal distribution is normal and expected.

---

### Spot-checking a Region

If you have wgbstools installed you can quickly visualise the methylation pattern at any genomic region directly from the BED file:

```bash
wgbstools view -r chr1:910433-910476 data/methylation/methylation.bed
```

---

## 5. Deconvolution

Deconvolution uses the methylation data to estimate the proportions of different cell types present in your sample. The pipeline supports two algorithms - **UXM** and **NNLS** - which are fundamentally different in how they work and what inputs they require. It is important to understand this difference before running anything, since the two workflows are completely separate.

---

### Understanding the Two Algorithms

**UXM** Operates at the read level, looking at patterns of methylation across blocks of CpG sites within individual reads. It works directly from the aligned BAM via a `.pat.gz` file generated by wgbstools. Crucially, **the BED file produced by the methylation summary step is not used by UXM can can't be used as input to it** - the two are incompatible formats. UXM is more robust against lower coverage, and is recommended for most runs.

**NNLS** operates at the site level, comparing per-CpG beta values against an Illumina-based reference atlas. It requires the methylation summary BED file as its starting point, which is then converted to Illumina probe coordinates via the `deconvolution prep` step. **The PAT file used by UXM cannot be used as input to NNLS.** Because it relies on Illumina probe coordinates, it is better suited to high-coverage data where per-site estimates are reliable.

| | UXM | NNLS |
|---|---|---|
| Input | `.pat.gz` (from `wgbstools bam2pat`) | `.csv` (from `deconvolution prep`) |
| Starts from | Aligned BAM | Methylation summary BED |
| Atlas | UXM atlas | Illumina-based full atlas |
| Best for | Standard Nanopore coverage | High coverage data |

---

### UXM Deconvolution

UXM does not use the methylation summary BED file. YOu can run it directly after alignment, without running the methylation summary step at all.

**Step 1 - Generate a PAT file:**

The PAT format is a compact read-level methylation format used by wgbstools. Generate one from the aligned BAM:

```bash
wgbstools bam2pat data/alignment/my_experiment.aligned.sorted.bam --out_dir data/processed/
```

**Step 2 - Run deconvolution**

```bash
nanopore_analysis deconvolution --input-file data/processed/my_experiment.aligned.sorted.pat.gz --output-dir results/deconvolution/deconvoluted_output.csv --algorithm uxm --atlas data/atlas/UXM_atlas.tsv
```

Behind the scenes the pipeline will:
1. Filter the PAT file to only the CpG loci present in the atlas (`wgbstools view`)
2. Index the filtered file (`wgbstools index`)
3. Run `uxm deconv` against the atlas

---

### NNLS Deconvolution

First download the Illumina manifest if you haven't already:

```bash
resource_download manifest
```

Then run the deconvolution prep step to convert the BED file to Illumina coordinates:

```bash
resource_download manifest
```

Then run the deconvolution prep step to convert the BED file to Illumina coordinates:

```bash
nanopore_analysis deconvolution prep --bed-file data/methylation/methylation.bed \
  --manifest-file data/atlas/illumina_manifest.csv \
  --chunk-size 5000000
```

Then run deconvolution:

```bash
nanopore_analysis deconvolution run \ 
  --input-file data/processed/deconvolution_illumina.csv \
  --output-dir results/deconvolution/ \
  --algorithm nnls \
  --atlas data/atlas/full_atlas.csv
```

---

### Reading the Output

The result is a CSV at `results/deconvolution/deconvoluted_output.csv`. Each row is a cell type and each value is the estimated proportion for that sample, summing to approximately 1.

```
CellType,tutorial
Adipocytes,0.0000000
Bladder-Ep,0.0000000
Blood-B,0.0260316
Blood-Granul,0.0738950
Blood-Mono+Macro,0.0138440
Blood-NK,0.0317679
Blood-T,0.0000000
Bone-Osteob,0.0300103
...
```

A few things to keep in mind when interpreting results:

- **Coverage matters.** Lower coverage produces noisier estimates. A single POD5 file will give you a rough picture but not a precise one.
- **The atlas determines what cell types are reported.** If a cell type isn't in your atlas it won't appear in the output.
- **Small proportions may not be meaningful.** Values below ~2-3% are often within the noise for low-coverage samples.

--- 

## 6. BAM Utilities

These are standalone tools that operate on BAM files independeltly from the main pipeline. They are useful for exploring and quality-checking your data at any stage.

---

### Summarising Read Lengths

Before committing to a full pipeline run it is worth checking the read length distribution of your BAM. A healthy long-read run should show a broad spread of lengths, typically peaking somewhere between 5-20 kb depending on your library prep.

```bash
summarise_lengths data/alignment/my_experiment.aligned.sorted.bam --output-dir data/qc/
```

This works on both aligned and unaligned BAMs, so you can run it immediately after basecalling to catch any issues before alignment. The output is a CSV with two columns:


| Column | Content |
|---|---|
| `read_length` | Length of the read in base pairs |
| `count` | Number of reads at that length |

Saved as `<input_stem>_read_length_distribution.csv` in your output directory.

**Things to look out for:**
- Almost all reads below 1 kb - something probably went wrong during library preparation
- A sharp peak at a very specific length - may indicate adapter contamination
- A completely flat distribution - worth double-checking your basecalling model

Of course these may be expected depending on your sample, so they are only worrying if your sample isn't expected to produce these results.

---

### Filtering by Read Length

Longer reads tend to cover more CpG sites in a single read, which gives the UXM algorithm more complete methylation block information to work with. Filtering to longer reads before deconvolution can improve result quality, particularly when working with low coverage data.

```bash
filter_by_length data/alignment/my_experiment.aligned.sorted.bam --cutoff 5000 --mode above --output-dir data/filtered/
```

The `--mode` argument controls which reads are kept:

| Mode | Behaviour |
|---|---|
| `above` | Keep only reads longer than the cutoff |
| `below` | Keep only reads shorter than or equal to the cutoff |
| `both` | Write both files separately |

The output BAM is indexed automatically. You can then use it as the input to any downstream step in place of the original BAM.

**Choosing a cutoff** depends on your data. A cutoff of 5000 bp is a reasonable starting point for whole-genome methylation data, but you should check your readlength distribution first - if your run peaked at 8 kb, a 5000 bp cutoff will keep most of your reads; if it peaked at 3 kb, you may end up with very little data.

> **Note:** Filtering to long reads reduces the total number of reads available, which lowers coverage. There is a tradeoff between read quality and coverage - if your data is already low coverage then filtering aggressively may do more harm than good.
