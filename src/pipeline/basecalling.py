import argparse
import logging
import os
from pathlib import Path

from src.utils.cli_utils import add_io_arguments
from src.utils.tools_runner import ToolRunner
from src.utils.file_utils import ensure_dir_exists, is_bam_multiplexed

logger = logging.getLogger(__name__)


def full_basecalling_handler(args, config):
    logger.info("Running full basecalling step.")

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

    run_basecalling(dorado_exe, input_file, model_speed, modifications, kit_name, batchsize, basecalled_bam)

    if is_bam_multiplexed(basecalled_bam):
        logging.info(f"The basecalled BAM {basecalled_bam} is multiplexed, so we'll demultiplex it.")
        run_demultiplex(dorado_exe, basecalled_bam, demultiplexed_output_dir)
    else:
        logging.info(f"The basecalled BAM {basecalled_bam} is not multiplexed, so we won't demultiplex it.")



def basecall_handler(args, config):
    input_file = args.input_file
    output = args.output_dir
    kit_name = args.kit_name
    model_speed = config.pipeline_steps.basecalling.params.complex_settings.model_speed
    modifications = config.pipeline_steps.basecalling.params.complex_settings.basecalling_modifications
    batchsize = config.pipeline_steps.basecalling.params.batch_size
    dorado_exe = config.tools.dorado

    run_basecalling(dorado_exe, input_file, model_speed, modifications, kit_name, batchsize, output)


def demultiplex_handler(args, config):
    input_file = args.input_file
    output_dir = args.output_dir
    dorado_exe = config.tools.dorado

    required_params = {
        "Input file (--input-file)": input_file,
        "Output directory (--output-dir)": output_dir,
        "Dorado exe (config only)": dorado_exe
    }

    errors = []
    for name, value in required_params.items():
        if value is None:
            errors.append(f" - {name}")

    if errors:
        error_report = "\n".join(errors)
        raise ValueError(f"Missing required parameters:\n{error_report}")

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


def run_demultiplex(dorado_exe, input_file, output_dir: Path):
    dorado_runner = ToolRunner(dorado_exe, '--output-dir')

    ensure_dir_exists(output_dir)

    demux_cmd = [
        "demux",
        "--kit-name", "SQK-NBD114-24",
        input_file
    ]

    dorado_runner.run(demux_cmd, str(output_dir))


def run_basecalling(dorado_exe, pod5_input, model_speed, modifications, kit_name, batchsize, output_path: Path = None):
    basecalling_cmd = ["basecaller",
                       f"{model_speed},{modifications}",
                       str(pod5_input),
                       "--kit-name", kit_name,
                       "--no-trim",
                       "--batchsize", str(batchsize)
                       ]

    final_output_arg = []
    dorado_runner = ToolRunner(dorado_exe)

    if output_path.suffix in ['.bam', '.sam']:
        # If the given output is a bam or sam file, treat it as a file.
        logging.info(f"Output path is a file. Redirecting Dorado output to: {output_path}")
        ensure_dir_exists(output_path.parent)
        # As ToolRunner output_file parameter already handles piping to an output file, don't add anything to the command.
    elif output_path.suffix == '':
        # If the user provided a directory path
        print(f"Output path is a directory. Using Dorado's -o flag: {output_path}")
        ensure_dir_exists(output_path)
        final_output_arg = ["-o", str(output_path)]
    else:
        # Invalid input
        raise ValueError(
            f"Invalid output path '{output_path}'. It must be a directory or a .bam/.sam file."
        )

    full_cmd = basecalling_cmd + final_output_arg

    if output_path.suffix in ['.bam', '.sam']:
        dorado_runner.run(full_cmd, output_path)
    else:
        dorado_runner.run(full_cmd)

def _add_kit_name_arg(parser, config):
    parser.add_argument(
        "--kit-name", type=str,
        default=config.pipeline_steps.basecalling.params.complex_settings.kit_name,
        help="Specify the Nanopore kit name"
    )

def _add_model_name_arg(parser, config):
    parser.add_argument(
        '--model-name', type=str,
        default=config.pipeline_steps.basecalling.params.explicit_settings.base_model_name,
        help="Name of the model to download."
    )

def add_all_arguments_to_parser(parser, config):
    """
    Publically available function to add arguments from the basecalling step to a given parser.
    It's used by the main run command and doesn't include the I/O args.
    :param parser:
    :param config:
    :return:
    """
    _add_kit_name_arg(parser, config)
    _add_model_name_arg(parser, config)

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
    _add_kit_name_arg(p_run, config)
    add_io_arguments(
        p_run, config,
        default_input=config.pipeline_steps.basecalling.paths.full_pod5_path,
        default_output=config.pipeline_steps.basecalling.paths.full_demultiplexed_output_dir,
        input_file_help="Path to POD5 files",
        output_dir_help="Path to demultiplexed data directory." # Add a tag to differentiate between demultiplex and not demultiplex.
    )
    p_run.set_defaults(func=full_basecalling_handler)

    p_basecall = basecalling_subparsers.add_parser(
        'basecall', help="Basecall a POD5 file/directory.",
        parents=[parent_parser]
    )
    _add_kit_name_arg(p_basecall, config)
    add_io_arguments(
        p_basecall, config,
        default_input=config.pipeline_steps.basecalling.paths.full_pod5_path,
        default_output=config.pipeline_steps.basecalling.paths.full_unaligned_bam_path,
        input_file_help="Path to POD5 files",
        output_dir_help="Path to basecalled BAM file"
    )
    p_basecall.set_defaults(func=basecall_handler)

    p_demux = basecalling_subparsers.add_parser(
        'demux', help="Demultiplex a multiplexed BAM file.",
        parents=[parent_parser]
    )
    _add_kit_name_arg(p_demux, config)
    add_io_arguments(
        p_demux, config,
        default_input=config.pipeline_steps.basecalling.paths.full_unaligned_bam_path,
        default_output=config.pipeline_steps.basecalling.paths.full_demultiplexed_output_dir,
        input_file_help="Path to basecalled BAM file",
        output_dir_help="Directory to save demultiplexed BAM files to."
    )
    p_demux.set_defaults(func=demultiplex_handler)

    p_download = basecalling_subparsers.add_parser(
        'download-model', help="Download a specific Dorado model",
        parents=[parent_parser]
    )
    _add_model_name_arg(p_download, config)
    add_io_arguments(
        p_download, config,
        add_input=False,
        default_output=config.pipeline_steps.basecalling.paths.full_dorado_model_dir,
        output_dir_help="Directory to save demultiplexed BAM files to."
    )
    p_download.set_defaults(func=download_handler)
