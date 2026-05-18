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

This section walks through a complete analysis run using publicly available ONT data, starting from raw POD5 files and ending with a cell-type deconvolution result.

### The Sample Dataset

We will use the Genome in a Bottle 2025.01 release from Oxford Nanopore's open data repository: the HG002 sample, sequenced on a PromethION with R10.4.1 chemistry. The raw POD5 files are available without authentication at:

```
s3://ont-open-data/giab_2025.01/flowcells/HG002/
```

We will only download one file from the dataset to illustrate the workflow of this pipeline. We also need the **hg38 reference genome** and the **UXM methylation atlas**. The pipeline can download these and we will do so at the point in the walkthrough where they are needed..

### Step 0 - Set Up the Tutorial Directory

All tutorial files live in the `tutorial/` folder in the project root. This includes a `config.yaml` that is pre-configured for this walkthrough - all four pipeline steps enabled, hg38 as the reference genome, and the correct filenames throughout.

Start by navigating to the tutorial directory:

```bash
cd tutorial
```

All commands in this walkthrough are run from here. The pipeline picks up `config.yaml` from your current working directory automatically, so no extra flags are needed.

### Step 1 - Download and Prepare the Raw Data

The pipeline has a built-in command for downloading the tutorial dataset.

```bash
resource_download sample_data
```

This downloads a single POD5 file from HG002 flow cell data to `data/pod5/`. PromethION POD5 files are very large, so one file is enough to demonstrate the full pipeline. We won't get whole-genome coverage from a single file, but there will be sufficient reads for basecalling, alignment, and deconvolution to run meaningfully. If you want more coverage, you can increase `setup.params.num_files` in `tutorial/config.yaml`. You will also need to update `setup.downloads.data_download_url` to point one level up in the S3 bucket (i.e the flow cell directory rather than the specific file) so the pipeline has multiple files to choose from.

> **Using FAST5 data instead?**: If your dataset uses the older FAST5 format, set `input_format: "fast5"` under `basecalling.params` in `config.yaml`. The conversion to POD5 happens automatically at the start of the basecalling step.

If you want to download from a different URL or save to a different location:

```bash
resource_download sample_data \
  --url s3://my-bucket/my-data/.../ \
  --output-dir data/my_input/
```

### Step 2 - Download the Reference Resources

Download the atlas and the reference genome before starting the pipeline steps:

```bash
resource_download genome
resource_download atlas
```

The atlas (~4 MB) is saved to `data/atlas/UXM_atlas.tsv`. The reference genome (~3 GB) is downloaded and initialised via wgbstools and saved to `reference_genomes/`. Both only need to be downloaded once and can be reused across experiments.

### Step 3 - Basecalling

```bash
nanopore_analysis basecalling run
```

This runs Dorado on the POD5 file in `data/pod5` using the `hac` model with `5mCG_5hmCG` methylation detection. The output is an unaligned BAM at `data/basecalled_output/tutorial.bam`.

### Step 4 - Alignment

```bash
nanopore_analysis align run
```

This aligns the basecalled reads to hg38 with minimap2, sorts and indexes the output, and saves it to `data/alignment/tutorial.aligned.sorted.bam`. A flagstat QC report is also written to `data/alignment/qc/`.

### Step 5 - Methylation Summary

```bash
nanopore_analysis methylation run
```

This runs `modkit pileup` on the aligned BAM, producing a BED file of per-CpG methylation calls at `data/methylation/methylation.bed`.

### Step 6 - Deconvolution

Before deconvolution, generate a PAT file from the aligned BAM using wgbstools:

```bash
wgbstools bam2pat data/alignment/tutorial.aligned.sorted.bam --out-dir data/processed/
```

Then run deconvolution:

```bash
nanopore_analysis deconvolution --input-file data/processed/tutorial.aligned.sorted.pat.gz --output-dir results/deconvolution/deconvoluted_output.csv --algorithm uxm --atlas data/atlas/Atlas.U250.l4.hg38.full.tsv
```

Progress for each step is logged to the terminal and written in full to `logs/`.


### Reading the Output

The result is a CSV at `results/deconvolution/deconvoluted_output.csv`. Each row is a cell type and each value is the estimated proportion for that sample:

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

Since we're only using a single POD5 file rather than a full run, the deconvolution result will be noisy and should not be interpreted as a definitive cell-type profile. The important thing at this stage is that the pipeline completes without errors and produces a CSV with plausible non-zero values across cell types.

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

Or override the input/output without touching the config:

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

A healthy run shows a broad spread of lengths, typically peaking somewhere between 5-20 kb depending on the library prep. If almost all reads are below 1 kb, something likely went wrong during library preparation rather than basecalling.
