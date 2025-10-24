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
from src.utils.file_utils import ensure_dir_exists, decompress_file
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
        logger.info("Download complete.")
        return destination
    except requests.exceptions.RequestException as e:
        logger.info(f"Error: Failed to download file from {url}. Reason: {e}")
        if destination.exists():
            destination.unlink() # Delete the partially downloaded file.
        return None


def reference_genome_handler(args, config):
    url = resolve_param(
        args, config, arg_name='url',
        config_path=['paths', 'reference_genome_url']
    )

    use_wgbs=args.wgbstools
    if use_wgbs:
        logger.info("We are going to initialise the genome using wgbstools")
    else:
        logger.info("We are going to initialise the genome and index with minimap2")

    ref_fasta = resolve_param(
        args, config, arg_name="output_dir", construct_path=True,
        config_path=[
            ['paths', 'reference_genome_dir'],
            ['paths', 'indexed_ref_gen_fasta_name']
        ]
    )

    logger.info(f"Setting up reference genome {Path(ref_fasta).name}")

    ref_mmi = Path(ref_fasta).with_suffix('.mmi')

    logger.info(f"We'll index the reference genome to {ref_mmi}")
    if not os.path.exists(ref_fasta):
        # Maybe we'll hardcode the reference genome urls in this function...
        logger.info(f"Reference file {ref_fasta} doesn't exist. Downloading from {url}")
        run_command([
            "aws", "s3", "cp", url, str(ref_fasta), "--no-sign-request"
        ])
    else:
        logger.info("Reference genome already exists.")

    if not os.path.exists(ref_mmi):
        logger.info("Indexing reference genome with minimap2...")
        run_command([
            "minimap2", "-d", str(ref_mmi), str(ref_fasta)
        ])
    else:
        logger.info(f"Reference genome index already exists.")


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
        logger.info(
            "The user didn't provide an output dir, so the output is taken from config and is assumed to be fine.")

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
        logger.info(
            "The user didn't provide an output dir, so the output is taken from config and is assumed to be fine.")

    final_destination_and_download(url, Path(destination), is_interactive)


def final_destination_and_download(url: str, destination: Path, is_interactive: bool = False):
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

    downloaded_file_path = _download_file_with_progress(url, final_destination)

    if downloaded_file_path:
        if downloaded_file_path.suffix.lower() in ['.gz', '.zip']:
            # The download was successful and it's a compressed file.
            final_path = decompress_file(downloaded_file_path)
            return final_path
        else:
            return downloaded_file_path
    else:
        return None


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
        type=str,
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
