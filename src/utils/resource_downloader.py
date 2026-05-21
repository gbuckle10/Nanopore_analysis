import argparse
import fnmatch
import logging
import os
import sys
import threading
from pathlib import Path
import requests
from tqdm import tqdm
import boto3
from botocore import UNSIGNED
from botocore.config import Config

from src import PROJECT_ROOT
from src.config.models import load_and_validate_configs
from src.config.paths import build_config_paths, update_config_from_args
from src.utils.file_utils import ensure_dir_exists, decompress_file
from src.utils.logger import Logger
from src.utils.process_utils import run_command, spinner
from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def _s3_client():
    return boto3.client('s3', config=Config(signature_version=UNSIGNED))


def list_s3_files(url: str, include: str = None) -> list[str]:
    """Return S3 URIs of files at url, optionally filtered by fnmatch pattern"""
    bucket, prefix = url.replace("s3://", "").split("/", 1)
    client = _s3_client()
    paginator = client.get_paginator('list_objects_v2')
    files = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get('Contents', []):
            key = obj['Key']
            s3_uri = f"s3://{bucket}/{key}"
            if include is None or fnmatch.fnmatch(key, include):
                files.append(s3_uri)
    return files


def _download_s3(url: str, destination: Path, format: str = None, include: str = None, exclude: str = None, max_files: int = None):
    """Download using s3"""
    destination.mkdir(parents=True, exist_ok=True)
    print(f"Downloading s3 with boto3 ")
    files = list_s3_files(url, include=include)
    if max_files is not None:
        files = files[:max_files]
        logger.info(f"Downloading {max_files} file(s).")
    client = _s3_client()
    bucket = url.replace("s3://", "").split("/")[0]
    for s3_uri in files:
        key = s3_uri.replace(f"s3://{bucket}/", "")
        filename = s3_uri.split("/")[-1]
        logger.info(f"Downloading {filename}")
        client.download_file(bucket, key, str(destination / filename))


def download(url: str, destination: Path, **kwargs):
    """Download something"""
    if url.startswith("s3://"):
        _download_s3(url, destination, **kwargs)
    else:
        final_destination_and_download(url, destination)


def sample_data_handler(config):
    url = str(config.pipeline_steps.setup.downloads.data_download_url)
    destination = config.pipeline_steps.setup.paths.data_input_dir
    max_files = config.pipeline_steps.setup.params.num_files
    file_format = config.pipeline_steps.setup.params.input_format
    logger.info(f"Downloading sample data from {url} to {destination}")
    download(url, Path(destination), include=f"*.{file_format}", max_files=max_files)


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
            destination.unlink()  # Delete the partially downloaded file.
        return None


def reference_genome_handler(config):
    url = config.pipeline_steps.setup.downloads.reference_genome_url
    use_wgbs = config.pipeline_steps.setup.params.download_ref_with_wgbstools
    print(f"Use wgbstools? - {use_wgbs}")
    genome_id = config.pipeline_steps.align.paths.genome_id

    if use_wgbs:
        logger.info(f"We are going to initialise genome {genome_id} using wgbstools")
        run_command(["wgbstools", "init_genome", genome_id])
        return

    logger.info("We are going to initialise the genome and index with minimap2")

    ref_path = config.pipeline_steps.align.paths.full_ref_fasta_path
    user_path = Path(ref_path)

    # Determine final file path
    final_ref_path: Path

    if (user_path.exists() and user_path.is_dir()) or (not user_path.exists() and user_path.suffix == ''):
        # If the user path exists, it's a directory. If it doesn't exist, there is no file extension.
        logger.info(f"Reference path '{user_path}' is a directory. Appending default filename.")
        filename = str(url).split('/')[-1]
        final_ref_path = user_path / filename
    else:
        logger.info(f"Reference path is a full file path: '{user_path}'")
        final_ref_path = user_path

    logger.info(f"The final output path will be {final_ref_path}")

    if not os.path.exists(final_ref_path) or getattr(config, 'force', False):
        logger.info(f"Reference file {final_ref_path} doesn't exist. Downloading from {url}")
        s3_url = str(url).replace("s3://", "")
        bucket, key = s3_url.split("/", 1)
        _s3_client().download_file(bucket, key, str(final_ref_path))
    else:
        logger.info("Reference genome already exists.")

    ref_mmi = final_ref_path.with_suffix('.mmi')
    logger.info(f"We'll index the reference genome to {ref_mmi}")

    if not os.path.exists(ref_mmi):
        logger.info("Indexing reference genome with minimap2...")
        stop = threading.Event()
        t = threading.Thread(target=spinner, args=(stop,), kwargs={"message": "Indexing reference genome..."})
        t.start()
        try:
            run_command([
                "minimap2", "-d", str(ref_mmi), str(final_ref_path)
            ])
        finally:
            stop.set()
            t.join()
    else:
        logger.info(f"Reference genome index already exists.")


def atlas_handler(config):
    url = config.pipeline_steps.setup.downloads.uxm_atlas_url

    destination = config.pipeline_steps.analysis.paths.full_atlas_path
    if not destination:
        sys.exit("Error: No output path specified and config doesn't contain the path.")

    # The interactive tag needs to be fixed - probs will need to add the args back into this.
    '''
    # Has the user provided an output directory? If so, we're in interactive mode.
    # is_interactive = args.output_dir is not None
    if is_interactive:
        logger.info("The user provided an output dir, so we are running in interactive mode.")
    else:
        logger.info(
            "The user didn't provide an output dir, so the output is taken from config and is assumed to be fine.")
    '''
    final_destination_and_download(url, Path(destination), is_interactive=False)


def manifest_handler(config):
    url = config.pipeline_steps.setup.downloads.manifest_url

    destination = config.pipeline_steps.analysis.paths.full_manifest_path

    if not destination:
        sys.exit("Error: No output path specified and config doesn't contain the path.")
    # Has the user provided an output directory? If so, we're in interactive mode.

    '''
    is_interactive = args.output_dir is not None

    if is_interactive:
        logger.info("The user provided an output dir, so we are running in interactive mode.")
    else:
        logger.info(
            "The user didn't provide an output dir, so the output is taken from config and is assumed to be fine.")
    '''
    final_destination_and_download(url, Path(destination), is_interactive=False)


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


def _add_subcommands(subparsers, parent_parser, config):
    download_parent = argparse.ArgumentParser(add_help=False, parents=[parent_parser])
    download_parent.add_argument('--force', action="store_true",
                                 dest="force",
                                 help="Force redownload even if the file already exists.")
    p_sample = subparsers.add_parser('sample_data', help="Download sample data.", parents=[download_parent])
    p_sample.add_argument("--url", type=str,
                          default=str(config.pipeline_steps.setup.downloads.data_download_url),
                          dest="pipeline_steps.setup.downloads.data_download_url")
    p_sample.add_argument("--output-dir", type=Path, default=None,
                          dest="pipeline_steps.setup.paths.data_input_dir",
                          help="Directory to save downloaded files into.")
    sample_limit_group = p_sample.add_mutually_exclusive_group()
    sample_limit_group.add_argument("--max-files", type=int, default=config.pipeline_steps.setup.params.num_files,
                                    dest="max_files",
                                    help="Maximum number of files to download.")
    sample_limit_group.add_argument("--all-files", action="store_const", const=None, dest="max_files",
                                    help="Download all available files, ignoring any limit")
    p_sample.set_defaults(func=sample_data_handler)

    p_genome = subparsers.add_parser('genome', help="Download and prepare the reference genome.",
                                     parents=[download_parent])
    p_genome.add_argument("--url", type=str,
                          default=config.pipeline_steps.setup.downloads.reference_genome_url,
                          dest="pipeline_steps.setup.downloads.reference_genome_url")
    p_genome.add_argument("--output-dir", type=Path, default=None,
                          dest="pipeline_steps.align.paths.custom_fasta_reference")
    p_genome.set_defaults(func=reference_genome_handler)

    p_atlas = subparsers.add_parser('atlas', help="Download the methylation atlas.",
                                    parents=[download_parent])
    p_atlas.add_argument("--url", type=str,
                         default=config.pipeline_steps.setup.downloads.uxm_atlas_url,
                         dest="pipeline_steps.setup.downloads.uxm_atlas_url")
    p_atlas.add_argument("--output-dir", type=Path,
                         default=config.pipeline_steps.analysis.paths.full_atlas_path,
                         dest="pipeline_steps.analysis.paths.atlas_file_name")
    p_atlas.set_defaults(func=atlas_handler)

    p_manifest = subparsers.add_parser('manifest', help="Download the Illumina manifest.",
                                       parents=[download_parent])
    p_manifest.add_argument("--url", type=str,
                            default=config.pipeline_steps.setup.downloads.manifest_url,
                            dest="pipeline_steps.setup.downloads.manifest_url")
    p_manifest.add_argument("--output-dir", type=Path,
                            default=config.pipeline_steps.analysis.paths.full_manifest_path,
                            dest="pipeline_steps.analysis.paths.manifest_name")
    p_manifest.set_defaults(func=manifest_handler)

def main():
    Logger.setup_logger()

    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('-u', '--user-config', default='config.yaml', type=Path)
    pre_parser.add_argument('-r', '--runtime-config', default='runtime_config.yaml', type=Path)
    conf_args, _ = pre_parser.parse_known_args()

    config = load_and_validate_configs(conf_args.user_config, conf_args.runtime_config)
    build_config_paths(config)

    parser = argparse.ArgumentParser(description="Download resources for the nanopore analysis pipeline.")
    parser.add_argument('-u', '--user-config', default='config.yaml', type=Path)
    parser.add_argument('-r', '--runtime-config', default='runtime_config.yaml', type=Path)
    subparsers = parser.add_subparsers(dest='subcommand', metavar='<resource>')
    _add_subcommands(subparsers, pre_parser, config)

    args = parser.parse_args()
    update_config_from_args(config, args, parser)

    if hasattr(args, 'func'):
        args.func(config)
    else:
        parser.print_help()
