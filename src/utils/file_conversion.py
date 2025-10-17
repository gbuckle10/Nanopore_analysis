import os
from pathlib import Path
import logging
import requests
import sys
import gzip
import shutil

from src.utils.process_utils import run_command

project_root = Path(__file__).resolve().parent


def download_atlas_manifest_files(config):
    print("--- Downloading and preparing atlas files and manifests ---")
    atlas_dir = os.path.join(project_root, "data/atlas")
    os.makedirs("data/atlas", exist_ok=True)

    files_to_download = {
        "illumina_manifest.csv": "https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv",
        "full_atlas.csv": "https://github.com/nloyfer/meth_atlas/raw/refs/heads/master/full_atlas.csv.gz",
        "UXM_atlas.tsv": "https://raw.githubusercontent.com/nloyfer/UXM_deconv/refs/heads/main/supplemental/Atlas.U25.l4.hg19.tsv"
    }

    for filename, url in files_to_download.items():
        dest_path = os.path.join(atlas_dir, filename)
        if not os.path.exists(dest_path):
            download_file(url, dest_path)
        else:
            print(f"{filename} already exists. Skipping download.")

    gz_path = os.path.join(atlas_dir, "full_atlas.csv.gz")
    csv_path = os.path.join(atlas_dir, "full_atlas.csv")
    if os.path.exists(gz_path) and not os.path.exists(csv_path):
        print("Decompressing full_atlas.csv.gz...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(csv_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    # Convert UXM_atlas.tsv to .csv
    tsv_path = os.path.join(atlas_dir, "UXM_atlas.tsv")
    uxm_csv_path = os.path.join(atlas_dir, "UXM_atlas.csv")
    if os.path.exists(tsv_path) and not os.path.exists(uxm_csv_path):
        print("Converting UXM_atlas.tsv to UXM_atlas.csv")
        with open(tsv_path, 'r') as tsv_file, open(uxm_csv_path, 'w', newline='') as csv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            writer = csv.writer(uxm_csv_path, delimiter=',')
            for row in reader:
                writer.writerow(row)


def download_and_index_reference_genome_manual(config):
    """
    Downloads the reference genomes with AWS CLI and indexes it with minimap 2.
    """
    genome = config['paths']['reference_genome']
    print(f"--- Setting up reference genome {genome} ---")

    ref_dir = os.path.join(project_root, 'reference_genomes', genome)
    os.makedirs(ref_dir, exist_ok=True)
    ref_fasta = os.path.join(ref_dir, config['paths']['reference_genome_fasta_name'])
    ref_mmi = os.path.join(ref_dir, config['paths']['indexed_ref_gen_fasta_name'])

    if not os.path.exists(ref_fasta):
        # Maybe we'll hardcode the reference genome urls in this function...
        ref_url = config['paths']['reference_genome_url']
        print(f"Reference file {ref_fasta} doesn't exist. Downloading from {ref_url}")
        run_command([
            "aws", "c3", "cp", ref_url, ref_fasta, "--no-sign-request"
        ])
    else:
        print("Reference genome already exists.")

    if not os.path.exists(ref_mmi):
        print("Indexing reference genome with minimap2...")
        run_command([
            "minimap2", "-d", ref_mmi, ref_fasta
        ])
    else:
        print(f"Reference genome index already exists.")

def download_and_index_reference_genome(config):
    """
    Use wgbstools init_genome to initialise the specified genome.
    """
    genome = config['paths']['reference_genome']
    print(f"Initialising reference genome {genome}")
    wgbstools_cmd = [
        "wgbstools", "init_genome",
        genome
    ]

    run_wgbstools(wgbstools_cmd)

def download_fast5_data(config):
    print("--- Downloading sample fast5 data ---")
    fast5_dir = os.path.join(project_root, config['paths']['fast5_input_dir'])
    number_to_download = config['parameters']['setup']['num_fast5_files']
    print(f"Downloading {number_to_download} fast5 files into {fast5_dir}")
    url = config['paths']['fast5_download_url']

    os.makedirs(fast5_dir, exist_ok=True)

    cmd_file_list = [
        "aws", "s3", "ls", url, "--no-sign-request"
    ]

    try:
        result = subprocess.run(
            cmd_file_list, check=True, capture_output=True, text=True
        )

        all_files = []
        for line in result.stdout.strip().split('\n'):
            if ".fast5" in line:
                filename = line.split()[-1]
                all_files.append(filename)

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Failed to list S3 bucket contents, make sure AWS CLI is installed and working", file=sys.stderr)
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if str(number_to_download).lower() == 'all':
        files_to_download = all_files
        print(f"Preparing to download all {len(files_to_download)} files.")
    else:
        try:
            files_to_download = all_files[:int(number_to_download)]
            print(f"Preparing to download the first {len(files_to_download)} files")
        except ValueError:
            print(f"'num_fast5_files' in config is not a number or 'all'. Value is {number_to_download}",
                  file=sys.stderr)
            sys.exit(1)

    print("Starting download...")
    for filename in files_to_download:
        if not filename:
            # A blank line
            continue
        source_path = f"{url}{filename}"
        local_dest = os.path.join(fast5_dir, filename)

        if os.path.exists(local_dest):
            print(f"Skipping {filename}, already exists.")
            continue

        dl_cmd = [
            "aws", "s3", "cp", source_path, fast5_dir, "--no-sign-request", "--quiet"
        ]

        run_command(dl_cmd)

    print("--- Fast5 download complete! --- ")

def download_file(url, destination):
    print(f"Downloading from {url} to {destination}")

    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(destination, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}", file=sys.stderr)
        sys.exit(1)


def download_atlas_manifest_files(config):
    print("--- Downloading and preparing atlas files and manifests ---")
    atlas_dir = os.path.join(project_root, "data/atlas")
    os.makedirs("data/atlas", exist_ok=True)

    files_to_download = {
        "illumina_manifest.csv": "https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv",
        "full_atlas.csv": "https://github.com/nloyfer/meth_atlas/raw/refs/heads/master/full_atlas.csv.gz",
        "UXM_atlas.tsv": "https://raw.githubusercontent.com/nloyfer/UXM_deconv/refs/heads/main/supplemental/Atlas.U25.l4.hg19.tsv"
    }

    for filename, url in files_to_download.items():
        dest_path = os.path.join(atlas_dir, filename)
        if not os.path.exists(dest_path):
            download_file(url, dest_path)
        else:
            print(f"{filename} already exists. Skipping download.")

    gz_path = os.path.join(atlas_dir, "full_atlas.csv.gz")
    csv_path = os.path.join(atlas_dir, "full_atlas.csv")
    if os.path.exists(gz_path) and not os.path.exists(csv_path):
        print("Decompressing full_atlas.csv.gz...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(csv_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    # Convert UXM_atlas.tsv to .csv
    tsv_path = os.path.join(atlas_dir, "UXM_atlas.tsv")
    uxm_csv_path = os.path.join(atlas_dir, "UXM_atlas.csv")
    if os.path.exists(tsv_path) and not os.path.exists(uxm_csv_path):
        print("Converting UXM_atlas.tsv to UXM_atlas.csv")
        with open(tsv_path, 'r') as tsv_file, open(uxm_csv_path, 'w', newline='') as csv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            writer = csv.writer(uxm_csv_path, delimiter=',')
            for row in reader:
                writer.writerow(row)


def ensure_tool_symlink(link_path, target_path):
    logger = logging.getLogger('pipeline')

    logger.info(f"Ensuring '{link_path}' is a symlink to '{target_path}'")

    try:
        if link_path.exists():
            if link_path.is_symlink():
                if os.path.realpath(link_path) == os.path.realpath(target_path):
                    logger.info(f"Correct symlink for {link_path} already exists, so nothing to do.")
                else:
                    logger.info(f"Symlink exists, but points to the wrong target. Recreating.")
                    link_path.unlink()
                    link_path.symlink_to(target_path)
                    print("Symlink recreated successfully.")
            else:
                logger.info(f"'{link_path}' is a regular file, not a symlink. Removing and replacing.")
                link_path.unlink()  # remove file
                link_path.symlink_to(target_path)  # create the symlink
                logger.info(f"Symlink created successfully.")
        else:
            # Directly create the link
            logger.info(f"'{link_path}' does not exist. Creating symlink.")
            link_path.symlink_to(target_path)
            print("Symlink created successfully.")
    except Exception as e:
        logger.error(f"An error occurred: {e}")


def convert_fast5_to_pod5(config):
    print("--- Converting fast5 to pod5 ---")
    pod5_dir = os.path.join(project_root, config['paths']['pod5_dir'])
    fast5_dir = os.path.join(project_root, config['paths']['fast5_input_dir'])

    os.makedirs(config['paths']['pod5_dir'], exist_ok=True)

    pod5_cmd = [
        'pod5', 'convert', 'fast5',
        fast5_dir,
        '--output', pod5_dir,
        '--force-overwrite'
    ]

    print(f"Converting fast5 to pod5, command: {' '.join(pod5_cmd)}")
    run_command(pod5_cmd)


def apply_runtime_config(runtime_config="src/runtime_config.sh"):
    """
    Reads a shell script and applies the export PATH commands to the current
    process's environment.
    :param runtime_config:
    :return:
    """

    logger = logging.getLogger('pipeline')
    config_vars = {}
    try:
        with open(runtime_config, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if line.startswith('export '):
                    line = line[7:]  # slice off "export"

                if '=' in line:
                    key, value = line.split('=', 1)
                    config_vars[key] = value
                    logger.info(f"Loaded config: {key}={value}")
        return config_vars
    except FileNotFoundError:
        logger.error(f"Runtime config file '{runtime_config}' not found. PATH not updated.")
