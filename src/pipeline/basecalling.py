import argparse
import logging
import os
from pathlib import Path

from src.utils.cli_utils import add_io_arguments
from src.utils.tools_runner import ToolRunner
from src.utils.file_utils import ensure_dir_exists

logger = logging.getLogger(__name__)


def full_basecalling_handler(args, config):
    logger.info("INFO: Running full basecalling step.")

    input_file = args.input_file
    if input_file is None:
        raise ValueError("Could not determine the path for the POD5 data to basecall")

    # The output should depend on whether you want to demultiplex or not. Add a tag for that in future.
    basecalled_bam = config.pipeline_steps.basecalling.paths.full_unaligned_bam_path

    demultiplexed_output_dir = args.output_dir
    kit_name = args.kit_name
    model_speed = config.pipeline_steps.basecalling.params.complex_settings.model_speed
    modifications = config.pipeline_steps.basecalling.params.complex_settings.basecalling_modifications
    batchsize = config.pipeline_steps.basecalling.params.batch_size
    dorado_exe = config.tools.dorado

    # Add validation steps in the handlers!!

    run_basecalling(dorado_exe, input_file, model_speed, modifications, kit_name, batchsize, basecalled_bam)

    run_demultiplex(dorado_exe, basecalled_bam, demultiplexed_output_dir)


def basecall_handler(args, config):
    input_file = args.input_file
    output = args.output_dir
    kit_name = args.kit_name
    model_speed = config.pipeline_steps.basecalling.params.complex_settings.model_speed
    modifications = config.basecalling.params.complex_settings.basecalling_modifications
    batchsize = config.pipeline_steps.basecalling.params.batch_size
    dorado_exe = config.tools.dorado

    run_basecalling(dorado_exe, input_file, model_speed, modifications, kit_name, batchsize, output)


def demultiplex_handler(args, config):
    input_file = args.input_file
    output_dir = args.output_dir
    dorado_exe = config.tools.dorado

    run_demultiplex(dorado_exe, input_file, output_dir)


def download_handler(args, config):
    model_name = args.model_name
    dorado_exe = config.tools.dorado

    run_model_download(dorado_exe, model_name)


def run_model_download(dorado_exe, model_name):
    # This should be better defined - give the base model name and
    # use the specified modifications to download the relevant modification models

    logger.info(f"Downloading dorado model {model_name}")

    download_cmd = [
        "dorado", "download",
        "--model", model_name,
        "--models-directory", "models/"
    ]
    dorado_runner = ToolRunner(dorado_exe)
    logger.info(f"Downloading dorado model {model_name} with \n{' '.join(download_cmd)}")

    dorado_runner.run(download_cmd)

    logger.info(f"Dorado model successfully downloaded.")


def run_demultiplex(dorado_exe, input_file, output_dir):
    dorado_runner = ToolRunner(dorado_exe, '--output-dir')

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    demux_cmd = [
        # "dorado",
        "demux",
        # "--output-dir", "analysis/demultiplexed",
        "--kit-name", "SQK-NBD114-24",
        input_file
    ]

    dorado_runner.run(demux_cmd, output_dir)


def run_basecalling(dorado_exe, pod5_input, model_speed, modifications, kit_name, batchsize, output_file=None):
    basecalling_cmd = ["basecaller",
                       f"{model_speed},{modifications}",
                       str(pod5_input),
                       "--kit-name", kit_name,
                       "--no-trim",
                       "--batchsize", str(batchsize)
                       ]

    dorado_runner = ToolRunner(dorado_exe)
    ensure_dir_exists(Path(output_file).parent)
    dorado_runner.run(basecalling_cmd, output_file)


def setup_parsers(subparsers, parent_parser, config):
    basecalling_parser = subparsers.add_parser(
        "basecalling",
        help="Run basecalling and related tasks using Dorado.",
        description="This command groups contains tools for basecalling raw Nanopore output.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['basecall'],
        parents=[parent_parser]
    )

    def show_basecalling_help(args, config):
        """Default function to show help for the basecalling command group"""
        basecalling_parser.print_help()

    basecalling_parser.set_defaults(func=show_basecalling_help)

    basecalling_subparsers = basecalling_parser.add_subparsers(
        title="Available Commands",
        description="Choose one of the following actions to perform.",
        dest='subcommand',
        metavar="<command>"
    )

    p_run = basecalling_subparsers.add_parser(
        'run', help="Basecall a POD5 file/directory and demultiplex if necessary.",
        parents=[parent_parser]
    )
    p_run.add_argument(
        "--kit-name", type=str,
        default=config.pipeline_steps.basecalling.params.complex_settings.kit_name,
        help="Specify the Nanopore kit name"
    )
    add_io_arguments(
        p_run, config,
        default_input=config.pipeline_steps.basecalling.paths.pod5_dir,
        default_output=config.pipeline_steps.basecalling.paths.demultiplexed_output_dir,
        input_file_help="Path to POD5 files",
        output_dir_help="Path to demultiplexed data directory." # Add a tag to differentiate between demultiplex and not demultiplex.
    )
    p_run.set_defaults(func=full_basecalling_handler)

    p_basecall = basecalling_subparsers.add_parser(
        'basecall', help="Basecall a POD5 file/directory.",
        parents=[parent_parser]
    )
    p_basecall.add_argument(
        "--kit-name", type=str,
        default=config.pipeline_steps.basecalling.params.complex_settings.kit_name,
        help="Specify the Nanopore kit name"
    )
    add_io_arguments(
        p_basecall, config,
        default_input=config.pipeline_steps.basecalling.paths.pod5_dir,
        default_output=config.pipeline_steps.basecalling.paths.full_unaligned_bam_path,
        input_file_help="Path to POD5 files",
        output_dir_help="Path to basecalled BAM file"
    )
    p_basecall.set_defaults(func=basecall_handler)

    p_demux = basecalling_subparsers.add_parser(
        'demux', help="Demultiplex a multiplexed BAM file.",
        parents=[parent_parser]
    )
    p_demux.add_argument("--kit-name", type=str, help="Specify the Nanopore kit name")
    add_io_arguments(
        p_demux, config,
        default_input=config.pipeline_steps.basecalling.paths.full_unaligned_bam_path,
        default_output=config.pipeline_steps.basecalling.paths.demultiplexed_output_dir,
        input_file_help="Path to basecalled BAM file",
        output_dir_help="Directory to save demultiplexed BAM files to."
    )
    p_demux.set_defaults(func=demultiplex_handler)

    p_download = basecalling_subparsers.add_parser(
        'download-model', help="Download a specific Dorado model",
        parents=[parent_parser]
    )
    p_download.add_argument(
        '--model-name', type=str,
        default=config.pipeline_steps.basecalling.params.explicit_settings.base_model_name,
        help="Name of the model to download."
    )
    add_io_arguments(
        p_download, config,
        add_input=False,
        default_output=config.pipeline_steps.basecalling.paths.dorado_model_dir,
        output_dir_help="Directory to save demultiplexed BAM files to."
    )
    p_download.set_defaults(func=download_handler)
