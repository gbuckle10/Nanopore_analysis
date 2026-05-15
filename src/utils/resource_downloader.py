import argparse
import logging
import os
import sys
from pathlib import Path
import requests
from tqdm import tqdm

from src.config.models import load_and_validate_configs
from src.config.paths import build_config_paths, update_config_from_args
from src.utils.file_utils import ensure_dir_exists, decompress_file
from src.utils.logger import Logger
from src.utils.process_utils import run_command
from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def _download_s3(url: str, destination: Path, include: str = None, exclude: str = None):
    """Download using s3"""
    cmd = ["aws", "s3", "cp", url, str(destination), "--no-sign-request", "--recursive", "--quiet"]
    if include:
        cmd += ["--exclude", "*"]
        cmd += ["--include", include]
    elif exclude:
        cmd += ["--exclude", exclude]
    destination.mkdir(parents=True, exist_ok=True)
    run_command(cmd)


def download(url: str, destination: Path, **kwargs):
    """Download something"""
    if url.startswith("s3://"):
        _download_s3(url, destination, **kwargs)
    else:
        final_destination_and_download(url, destination)


def sample_data_handler(config):
    url = str(config.pipeline_steps.setup.downloads.fast5_download_url)
    destination = config.pipeline_steps.setup.paths.fast5_input_dir
    logger.info(f"Downloading sample fast5 data from {url} to {destination}")
    download(url, Path(destination), include="*.fast5")


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

    if use_wgbs:
        logger.info("We are going to initialise the genome using wgbstools")
    else:
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
        run_command([
            "aws", "s3", "cp", str(url), str(final_ref_path), "--no-sign-request"
        ])
    else:
        logger.info("Reference genome already exists.")

    ref_mmi = final_ref_path.with_suffix('.mmi')
    logger.info(f"We'll index the reference genome to {ref_mmi}")

    if not os.path.exists(ref_mmi):
        logger.info("Indexing reference genome with minimap2...")
        run_command([
            "minimap2", "-d", str(ref_mmi), str(final_ref_path)
        ])
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
    download_parent = argparse.ArgumentParser(add_help=False)
    download_parent.add_argument('--force', action="store_true",
                                 dest="force",
                                 help="Force redownload even if the file already exists.")
    p_sample = subparsers.add_parser('sample_data', help="Download sample fast5 data.", parents=[download_parent])
    p_sample.add_argument("--url", type=str,
                          default=str(config.pipeline_steps.setup.downloads.fast5_download_url),
                          dest="pipeline_steps.setup.downloads.fast5_download_url")
    p_sample.add_argument("--output-dir", type=Path, default=None,
                          dest="pipeline_steps.setup.paths.fast5_input_dir")
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


def setup_parsers(subparsers, parent_parser, config):
    download_parser = subparsers.add_parser(
        "download", help="Download files necessary for the pipeline.",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[parent_parser]
    )
    download_subparsers = download_parser.add_subparsers(
        title="Available Commands", dest='subcommand', metavar="<command>"
    )
    _add_subcommands(download_subparsers, parent_parser, config)


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
