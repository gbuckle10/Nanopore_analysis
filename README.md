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

Accept the licence agreement. When asked whether to automatically initialise conda, choose **no**. If you accidentally said yes you can undo it after restarting WSL with:
```bash
conda config --set auto_activate_base false
```

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

**4. Create the entry-point symlink** (if not created automatically)

```bash
ln -s src/run_pipeline.py nanopore_analysis
chmod +x src/run_pipeline.py
```

You can now run the pipeline from the project root:

```bash
./nanopore_analysis --help
```

Or add the project root to your PATH to call it from anywhere:

```bash
export PATH=$PATH:$PWD
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

