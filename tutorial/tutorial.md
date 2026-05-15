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

This section walks through a complete analysis run using publicly available ONT data, going from a pre-aligned BAM to a cell-type deconvolution.

### The Sample Dataset

We will use a dataset from Oxford Nanopore's open data repository: a GM24385 (human lymphoblastoid cell line) sequencing run. The raw FAST5 files are available without authentication at:
 
```
s3://ont-open-data/gm24385_mod_2021.09/flowcells/20210511_1515_X2_FAQ32637_9b683def/fast5_pass/
```

> **NOTE:** This dataset uses FAST5 format rather than the newer POD5 format, so we will need to convert it before basecalling. This is covered in step 1. We plan to replace this with a POD5 dataset in a future version of the tutorial.

A preconfigured `config.tutorial.yaml` is included in the tutorial folder. It has all 4 pipeline steps enabled, hg19 set as the reference genome, and the correct filenames for this walkthrough. You don't need to edit anything - just pass it to every command with `--user-config`

```bash
nanopore_analysis --user-config config.tutorial.yaml <command>
```

> **Why hg19?** The GM24385 dataset was sequenced against hg19, so all coordinates are hg19-based. The atlas and manifest we download must match. If you later run the pipeline on your own data aligned to hg38, use the default `config.yaml` instead.

### Step 1 - Download and Prepare the Raw Data

First, download a subset of the FAST5 files. The full dataset is very large, so we will just take a few files to keep things manageable.

```bash
aws s3 cp \
  "s3://ont-open-data/gm24385_mod_2021.09/extra_analysis/alignment/20210510_1600_X1_FAQ32172_f02f2d1c.bam" \
  "data/alignment/taken_alignment.bam" \
  --no-sign-request
```
 
Sort and index it, since modkit requires a sorted, indexed BAM:
 
```bash
samtools sort -o data/alignment/taken_alignment.sorted.bam data/alignment/taken_alignment.bam
samtools index data/alignment/taken_alignment.sorted.bam
rm data/alignment/taken_alignment.bam
```