import argparse
import csv
import gzip
import logging
import os
import shutil
import sys
from pathlib import Path
import requests
from tqdm import tqdm

from src.utils.cli_utils import create_io_parser
from src.utils.config_utils import resolve_param
from src.utils.file_utils import ensure_dir_exists
from src.utils.process_utils import run_command
from src.utils.tools_runner import ToolRunner

project_root = Path(__file__).resolve().parent
logger = logging.getLogger(__name__)

def _download_file_with_progress(url: str, destination: Path):
    """Downloads a file from a URL, showing a progress bar."""
    print(f"Downloading from {url} to {destination}")

    try:
        response = requests.get(url, stream=True, allow_redirects=True)
        response.raise_for_status()
        total_size = int(response.headers.get('content-length', 0))

        destination.parent.mkdir(parents=True, exist_ok=True)

        with open(destination, 'wb') as f, tqdm(
            desc=destination.name,
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                size = f.write(chunk)
                bar.update(size)
            print("Download complete.")
            return True
    except requests.exceptions.RequestException as e:
        print(f"Error: Failed to download file from {url}. Reason: {e}")
        return False

def reference_genome_handler(args, config):
    url = resolve_param(
        args, config, arg_name='url',
        config_path=['paths', 'reference_genome_url']
    )
    destination = resolve_param(
        args, config, arg_name='output_dir', construct_path=True,
        config_path=[
            ['paths', 'reference_genome_dir'],
            ['paths', 'indexed_ref_gen_fasta_name']
        ]
    )

def atlas_handler(args, config):
    url = resolve_param(
        args, config, arg_name='url',
        config_path=['paths', 'manifest_url']
    )

    destination = resolve_param(
        args, config, arg_name='output_dir', construct_path=True,
        config_path=[
            ['paths', 'atlas_dir'],
            ['paths', 'illumina_manifest.csv']
        ]
    )
    if not destination:
        sys.exit("Error: No output path specified and config doesn't contain the path.")

    # Has the user provided an output directory? If so, we're in interactive mode.
    is_interactive = args.output_dir is not None

    if is_interactive:
        logger.info("The user provided an output dir, so we are running in interactive mode.")
    else:
        logger.info("The user didn't provide an output dir, so the output is taken from config and is assumed to be fine.")

    final_destination_and_download(url, Path(destination), is_interactive)

def manifest_handler(args, config):
    url = resolve_param(
        args, config, arg_name='url',
        config_path=['paths', 'manifest_url']
    )

    destination = resolve_param(
        args, config, arg_name='output_dir', construct_path=True,
        config_path=[
            ['paths', 'atlas_dir'],
            ['paths', 'illumina_manifest']
        ]
    )

    if not destination:
        sys.exit("Error: No output path specified and config doesn't contain the path.")
    # Has the user provided an output directory? If so, we're in interactive mode.

    is_interactive = args.output_dir is not None

    if is_interactive:
        logger.info("The user provided an output dir, so we are running in interactive mode.")
    else:
        logger.info("The user didn't provide an output dir, so the output is taken from config and is assumed to be fine.")

    final_destination_and_download(url, Path(destination), is_interactive)

def final_destination_and_download(url: str, destination: Path, is_interactive: bool=False):

    logger.info(f"We will download the file from {url} to the destination {destination}")
    final_destination: Path

    if destination.suffix == '' or destination.is_dir():
        logger.info("The provided destination is a directory, so I'll take the filename from the url.")
        filename = url.split('/')[-1]
        if not filename:
            raise ValueError("Could not determine filename from URL.")
        final_destination = destination / filename
    else:
        logger.info("The provided destination is a file, so that is simply the final destination.")
        final_destination = destination

    logger.info(f"The final destination for this file is {final_destination}")

    parent_dir = final_destination.parent
    if not ensure_dir_exists(parent_dir, is_interactive):
        logger.warning("Directory creation failed or was cancelled by the user. Aborting download.")
        return False

    _download_file_with_progress(url, final_destination)



def download_atlas_manifest_files(config):
    print("--- Downloading and preparing atlas files and manifests ---")

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


def _download_and_index_reference_genome(config):
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

def download_and_index_reference_genome_wgbs(config):
    """
    Use wgbstools init_genome to initialise the specified genome.
    """
    genome = config['paths']['reference_genome']
    print(f"Initialising reference genome {genome}")
    wgbstools_cmd = [
        "wgbstools", "init_genome",
        genome
    ]

    wgbstools_exe = config['submodules']['wgbstools']
    wgbstools_runner = ToolRunner(wgbstools_exe)

    wgbstools_runner.run(wgbstools_cmd)




def setup_parsers(subparsers, parent_parser):
    io_parser = create_io_parser()
    download_parent_parser = argparse.ArgumentParser(add_help=False)
    download_parent_parser.add_argument(
        "--url",
        type=Path,
        help="URL to download from"
    )
    download_parent_parser.add_argument(
        '--force',
        action="store_true",
        help="Force redownload even if the file already exists."
    )

    download_parser = subparsers.add_parser(
        "download",
        help="Download files necessary for deconvolution",
        description="This command group contains tools for downloading and preparing files necessary for deconvolution.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    download_subparsers = download_parser.add_subparsers(
        title="Available Commands",
        description="Choose one of the following actions",
        dest='subcommand',
        metavar="<command>"
    )

    p_genome = download_subparsers.add_parser(
        'genome',
        help="Download, and prepare the specified genomes. --wgbstools will run the wgbs_tools init_genome.",
        parents=[io_parser, download_parent_parser]
    )
    p_genome.add_argument(
        '--wgbstools',
        action="store_true",
        help="Downloads and initialises the genome using wgbs_tools' init_genome function"
    )
    p_genome.set_defaults(func=reference_genome_handler)

    p_atlas = download_subparsers.add_parser(
        'atlas',
        help="Download the specified methylation atlas.",
        parents=[io_parser, download_parent_parser]
    )
    p_atlas.set_defaults(func=atlas_handler)

    p_manifest = download_subparsers.add_parser(
        'manifest',
        help="Download the specified Illumina manifest.",
        parents=[io_parser, download_parent_parser]
    )
    p_manifest.set_defaults(func=manifest_handler)