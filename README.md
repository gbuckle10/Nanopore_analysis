# Nanopore Analysis Pipeline

A command-line pipeline for analysis of Oxford Nanopore sequencing data, from raw basecalling through to methylation-based cell type deconvolution.

## Table of Contents

- [Overview](#overview)
- [Pipeline Steps](#pipeline-steps)
- [Requirements](#requirements)
- [Installation](#installation)
  - [WSL Setup (Windows users)](#wsl-setup-windows-users)
  - [Conda Installation](#option-1-conda-recommended)
  - [Docker Installation](#option-2-docker)
- [Configuration](#configuration)
- [Usage](#usage)
- [Tutorial: Running the Full Pipeline](#tutorial-running-the-full-pipeline)
- [Working with wgbs_tools](#working-with-wgbs_tools)
- [Individual Commands](#individual-commands)
- [Troubleshooting](#troubleshooting)
- [Project Structure](#project-structure)

---

## Overview

This pipeline takes raw Nanopore sequencing output (POD5/FAST5 files) and processes it through basecalling (including native methylation detection), genome alignment and finally cell-type deconvolution based on the methylation data. 

wgbs_tools and UXM_deconv are used as a convenient deconvolution toolchain and are  not specific to bisulphite data - they operate on the methylation values extracted from the Nanopore BAM.

**Key tools used by the pipeline:**
- [Dorado](https://cdn.oxfordnanoportal.com/software/analysis/) - basecalling, native methylation calling, and alignment (downloaded from Oxford Nanopore's CDN by the installer)
- [minimap2](https://github.com/lh3/minimap2) - reference genome indexing (installed via conda)
- [samtools](http://www.htslib.org/) - BAM sorting, indexing and QC (installed via conda)
- [wgbs_tools](https://github.com/nloyfer/wgbs_tools) - PAT file generation and filtering for deconvolution (installed as a Git submodule; the pipeline uses a [slightly modified fork](https://github.com/gbuckle10/wgbs_tools) with minor import and path fixes)
- [UXM_deconv](https://github.com/nloyfer/UXM_deconv) - cell-type deconvolution (installed as a Git submodule; the pipeline uses a [slightly modified fork](https://github.com/gbuckle10/UXM_deconv) with minor import and path fixes)
- [modkit](https://github.com/nanoporetech/modkit) - optional methylation pileup to BED format (installed via conda)

--- 

## Pipeline Steps
 
```
POD5/FAST5 files
      │
      ▼
1. Basecalling        (Dorado, with methylation calling)  → unaligned BAM with methylation tags
      │
      ▼
2. Alignment          (Dorado aligner + samtools)          → aligned, sorted, indexed BAM
      │
      ▼
3. Deconvolution      (wgbs_tools + UXM or NNLS)          → cell-type proportion CSV
```

An optional demultiplexing step can be inserted between basecalling and alignment when barcoded libraries are used. Also, an optional methylation pileup step using `modkit` can be run between alignment and deconvolution to produce a BED file of per-CpG methylation levels.

---

## Requirements

- Linux or WSL (Windows Subsystem for Linux)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Miniforge](https://github.com/conda-forge/miniforge)
- Git
- AWS CLI v2 (only needed for downloading reference genomes or sample data from S3)
- Docker (optional, for the Docker-based installation)

---

## Installation

This pipeline coordinates several specialised bioinformatics tools, so the initial setup is quite involved. The installer handles as much as possible automatically, including downloading Dorado and compiling the submodules so most of the steps below are one-time only. When the initial setup is done once, running the pipeline subsequent times is very simple.

If you run into problems during installation, see the [Troubleshooting](#troubleshooting) section.

### WSL Setup (Windows users)

If you are running this on Windows, you first need to download WSL (Windows Subsystem for Linux).

**1. Install WSL**

Run the following in PowerShell, then restart your machine:
```
wsl --install
```
open WSL and create a username and password when prompted.

**2. Install Miniconda inside WSL**

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Accept the licence agreement. When asked whether to automatically initialise conda, choose **yes**.

Close and reopen WSL to finish the installation.

**3. Fix line endings (if you edit files on Windows)**

If you edit any project files in a Windows editor, you may get  `$'\r': command not found` errors in WSL due to Windows-style line endings (CRLF). If you're lucky you can change the line endings to LF directly in your editor:

- **PyCharm** - click the `CRLF` button in the bottom-right status bar and select `LF`.
- **VS Code** - click the `CRLF` button in the bottom-right status bar and select `LF`. To apply this to all new files automatically, add `"files.eol": "\n"` to your `settings.json`.
- **Vim/Neovim** - run `:set fileformat=unix`then save with `:w`. To make it permanent, add `set fileformat=unix` to your `.vimrc` or `init.vim` file.

If you can't change the setting in your editor, use `dos2unix` as a last resort:
```bash
sudo apt-get update
sudo apg-get install dos2unix
find . -type f -name "*.sh" -exec dos2unix {} \+ 
```

--- 

### Option 1: Conda (recommended)

**1. Clone the repository**

```bash
git clone https://github.com/gbuckle10/Nanopore_analysis.git
cd Nanopore_analysis
```

**2. Create the Conda environment**

```bash
conda env create -f environment.yml
conda activate nanopore_analysis
```

To update the environment later if `environment.yml` has changed:
```bash
conda env update --file environment.yml --prune
```

**3. Run the installer** (downloads Dorado and sets up submodules)

```bash
python install.py conda all
```

This will:
- Install external tools (Dorado) into `tools/`
- Initialise and compile Git submodules (wgbs_tools, UXM_deconv, meth_atlas)
- Create a `nanopore_analysis` symlink in the project root
- Copy Conda activation scripts so `PYTHONPATH` is set automatically on environment activation

**4. Install the pipeline:** (if not created automatically)

```bash
pip install -e .
```

This installs the package in editable mode and creates the `nanopore_analysis` command, which you can now call from anywhere:

```bash
nanopore_analysis --help
```

---

### Option 2: Docker

**1. Clone the repository**

```bash
git clone https://github.com/gbuckle10/Nanopore_analysis.git
cd Nanopore_analysis
```

**2. Build and install**

```bash
python install.py docker
```

This builds a Docker image tagged `nanopore-pipeline:latest` and creates a `nanopore_analysis` wrapper script that mounts `./data` and `./output` into the container automatically

```bash
./nanopore_analysis --help
```

---

## Configuration

All pipeline settings are controlled by `config.yaml` in the project root, which is intended to be self-explanatory. A few things to note if you want to edit it:

- **Use relative paths only** - do not use absolute paths for any files or directories
- **`pipeline_control.run_steps`** - toggle which steps run when using the `run` command. Steps set to `false` are skipped entirely.
- **`pipeline_steps.setup.params.dorado_version`** - controls which version of Dorado the installer downloads. Update this if you need a specific version.
- **`genome_id` vs `custom_fasta_reference`** - in the alignment section, provide one or the other but not both. `genome_id` (e.g. `hg38`) tells the pipeline to look for the reference genome in the standard `reference_genomes/` directory structure, while `custom_fasta_reference` lets you point directly to a FASTA file anywhere on disk.

### Runtime config

A second file, `runtime_config.yaml`, is generated automatically by the installer and stores paths to downloaded tools (e.g. the Dorado executable). Its values override `config.yaml` where they overlap. Do not edit this manually.

---

## Usage

Always activate the Conda environment first

```bash
conda activate nanopore_analysis
```

The general command structure is:

```
nanopore_analysis <command> [subcommand] options
```

Run `nanopore_analysis --help` or `nanopore_analysis <command> --help` for a full list of options at any given level.

---

## Tutorial: Running the Full Pipeline

This tutorial walks through a complete run from POD5 files to deconvolution output.

### Step 0: Download required resources

**Reference genome** (hg38 from iGenomes via S3):
```bash
nanopore_analysis download genome
```

The pipeline expects **UCSC-style chromosome headers** (e.g. `chr1`, `chrX`). The default download URL points to an hg38 genome with UCSC headers. Make sure any manually downloaded genome also uses this convention - using the wrong header style (e.g. `1` instead of `chr1`) will cause alignment and methylation steps to fail silently.

You can also download the genome manually using the AWS CLI:
```bash
aws s3 cp "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa" \
  reference_genomes/hg38/genome.fa --no-sign-request
```

To install AWS CLI v2:
```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```
 
**Methylation atlas** (UXM hg38 atlas):
```bash
nanopore_analysis download atlas
```

The atlas from [UXM_deconv](https://github.com/nloyfer/UXM_deconv/tree/main/supplemental) is used by default. Make sure the atlas genome assembly matches the one used for basecalling and alignment.

**Illumina manifest** (required for NNLS deconvolution):
```bash
nanopore_analysis download manifest
```

These files are saved to `data/atlas/` by default.

--- 

**Sample data for testing**

You can browse and download sample ONT data from the open data bucket. The `aws s3 ls` command is useful for exploring what is available:
```bash
aws s3 ls s3://ont-open-data/ --no-sign-request
```

To download a pre-aligned BAM for testing (bypassing basecalling entirely):
```bash
aws s3 cp "s3://ont-open-data/gm24385_mod_2021.09/extra_analysis/alignment/20210510_1600_X1_FAQ32172_f02f2d1c.bam" \
  data/alignment/taken_alignment.bam --no-sign-request
```

---

### Step 1: Basecalling

Place your POD5 files in `data/pod5/` (or update `pod5_input_path` in `config.yaml`), then run:

```bash
nanopore_analysis basecalling run
```

To override paths at the command line:

```bash
nanopore_analysis basecalling run \
  --input-file /path/to/pod5/ \
  --output-dir data/basecalled_output/basecalled.bam
```

If your library is barcoded and you want to demultiplex:

```bash
nanopore_analysis basecalling run --demultiplex
```

This will basecall first, then automatically check whether the output BAM is multiplexed and demultiplex it if so.

--- 

### Step 2: Alignment

```bash
nanopore_analysis align run
```

This will:
1. Check for a pre-build minimap2 index (`.mmi`) for your reference genome, and build one if missing
2. Align the basecalled BAM using Dorado aligner, piped into `samtools sort`
3. Index the sorted BAM with `samtools index`
4. Run `samtools flagstat` and save a QC report

To specify paths explicitly:

```bash
nanopore_analysis align run \
  --input-file data/basecalled_output/basecalled.bam \
  --output-dir data/alignment/aligned.sorted.bam \
  --ref reference_genomes/hg38/genome.fa
```
 
To run QC on an already-aligned BAM without re-aligning:
 
```bash
nanopore_analysis align qc --input-file data/alignment/aligned.sorted.bam
```

---

### Step 3: Methylation Pileup (optional)

When Dorado is run with a modification model (configured via `basecalling modifications` in `config.yaml`, e.g. `5mCG_5hmCG`), methylation tags are embedded directly in the BAM output during basecalling. In this case a separate pileup step is not needed for deconvolution; the pipeline can work directly from the BAM. Methylation calling is not on by default in Dorado and must be explicitly requested this way.

The pileup step is only necessary if you want a BED file of per-CpG methylation levels for other downstream analysis or visualisation tools. To run it:

```bash
nanopore_analysis methylation run
```

With explicit paths:

```bash
nanopore_analysis methylation run \
  --input-file data/alignment/aligned.sorted.bam \ 
  --output-dir data/methylation/methylation.bed \
  --ref reference_genomes/hg38/genome.fa
```

This runs `modkit pileup` to produce a BED file of per-CpG methylation levels. Providing `--ref` restricts output to CpG sites only, which is strongly recommended

---

### Step 4: Deconvolution

```bash
nanopore_analysis deconvolution
```

With options:

```bash
nanopore_analysis deconvolution \
  --algorithm uxm \
  --atlas data/atlas/UXM_atlas.csv
```

Two algorithms are supported:

- **`uxm`** (default) - uses wgbs_tools + UXM_deconv. Filters the PAT file against the atlas, indexes it, then runs `uxm deconv`. See [Working with wgbs_tools](#working-with-wgbs-tools) for manual PAT file preparation.
- **`nnls`** - uses the meth_atlas NNLS algorithm with a full atlas CSV.

Results are saved to `results/deconvolution/deconvoluted_output.csv` by default.

---

### Running the full pipeline in one command

If all steps are enabled in `config.yaml`:

```bash
nanopore_analysis run
```

You can override individual settings without editing the config:

```bash
nanopore_analysis run \
  --kit-name SQK-NBD114-24 \
  --algorithm uxm
```

---

## Working with wgbs_tools

The pipeline uses wgbs_tools internally during deconvolution, but you can also use it directly for exploration and QC. Full documentation is available on the [wgbs_tools Github page](https://github.com/nloyfer/wgbs_tools).

### Generating PAT files manually

`bam2pat` can process a whole directory of BAM files at once:

```bash
wgbstools bam2pat path/to/files/*.bam --out_dir path/to/output
```

### Filtering a PAT file against an atlas

Before deconvolution, filter the PAT file to retain only sites present in the atlas. The atlas used here must be the same one passed to `uxm deconv`, and both must match the genome assembly used for alignment.

```bash
wgbstools view path/to/file.pat.gz -L path/to/atlas -o path/to/output.pat
```

Then index the output:
```bash
wgbstools index path/to/output.pat.gz
```

### Running UXM deconvolution manually

```bash
uxm deconv path/*.pat.gz --atlas /path/to/atlas -o output.csv
```

### Viewing methylation data
 
View a specific genomic region from a PAT or beta file:
```bash
wgbstools view -r chr1:910433-910476 path/to/file.pat.gz
wgbstools view -r chr1:910433-910476 --genome hg19 path/to/file.beta
```
 
Visualise methylation across samples:
```bash
# Standard visualisation
wgbstools vis *.beta -r chr1:22517933-22519650
 
# Heatmap view
wgbstools vis *.beta -r chr1:22517933-22519650 --heatmap
 
# Visualise read-level methylation patterns from a PAT file
wgbstools vis path/to/file.pat.gz -r chr1:22517933-22519650
```

### Installing wgbs_tools separately
 
wgbs_tools is included as a Git submodule and compiled during installation. If you need to install it independently:
 
```bash
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
 
conda activate nanopore_analysis
python setup.py
chmod +x wgbstools 
```
 
If `python setup.py` fails with a compiler error, run `which g++` to check whether the C++ compiler is available. If it returns nothing, the symlinks are missing. Fix them by adjusting the path below to match where your Conda environments live:
 
```bash
cd /home/{USERNAME}/miniconda3/envs/wgbs_build/bin/
ln -s x86_64-conda-linux-gnu-g++ g++
ln -s x86_64-conda-linux-gnu-gcc gcc
ln -s x86_64-conda-linux-gnu-cpp cpp
ln -s x86_64-conda-linux-gnu-c++ c++
```
 
Then rerun `python setup.py` and make the file executable:
```bash
python setup.py
chmod +x wgbstools
```
 
---
 
## Individual Commands
 
### Filter a BAM by read length
 
```bash
python src/analysis/filter_bam_by_length.py input.bam \
  --cutoff 1000 \
  --mode both \
  --output-dir data/filtered/
```
 
Options for `--mode`: `above`, `below`, or `both`. Outputs separate BAM files for reads above and/or below the cutoff, and indexes them automatically.
 
### Summarise read length distribution
 
```bash
python src/analysis/summarise_lengths.py input.bam \
  --output-dir data/qc/
```
 
Produces a CSV with columns `read_length` and `count`, saved as `<input_stem>_read_length_distribution.csv`.
 
### Download tools individually
 
```bash
# Re-install Dorado at a specific version
python install.py conda tools --dorado-version 0.9.6
 
# Re-initialise submodules only
python install.py conda submodules
```
 
---
 
## Troubleshooting
 
**`PYTHONPATH` errors / module not found**
 
Run from the project root with the Conda environment active. These two exports should fix most remaining cases:
```bash
export PATH=$PATH:$PWD
export PYTHONPATH=$PYTHONPATH:$PWD/src
```
 
**`$'\r': command not found` in WSL**
 
Line endings are Windows-style (CRLF). See the [WSL setup section](#wsl-setup-windows-users) for how to fix this with `dos2unix` or your editor settings.
 
**Conda commands not found in a new WSL terminal**
 
```bash
eval "$(path/to/conda shell.bash hook)"
```
 
**Submodule directory is empty**
 
You may see an error like `The directory <Project>/externals/meth_atlas is registered as a Git root, but no Git repositories were found there`. Fix with:
```bash
git submodule sync
git submodule update --init --recursive
```
 
`git submodule sync` reads `.gitmodules` and updates your local Git configuration to match. It should print `Synchronising submodule url for 'externals/meth_atlas'` for each submodule.
 
**wgbs_tools compilation fails / `g++` not found**
 
See the [wgbs_tools standalone installation](#installing-wgbs_tools-separately) section for how to set up the missing compiler symlinks.
 
**Atlas genome assembly mismatch**
 
Make sure the atlas used for `wgbstools view` filtering and the atlas passed to `uxm deconv` are the same file, and that both match the genome assembly (hg38 or hg19) used during basecalling and alignment. Mixing assemblies will produce incorrect or empty results without an obvious error.
 
---




