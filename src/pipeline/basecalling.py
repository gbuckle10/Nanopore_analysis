import argparse
import logging
import os
from pathlib import Path

from src.utils.cli_utils import add_io_arguments
from src.utils.tools_runner import ToolRunner
from src.utils.file_utils import ensure_dir_exists, is_bam_multiplexed

logger = logging.getLogger(__name__)


def full_basecalling_handler(config):
    logger.info("Running full basecalling step.")

    input_file = config.pipeline_steps.basecalling.paths.full_pod5_path
    if input_file is None:
        raise ValueError("Could not determine the path for the POD5 data to basecall")
    # The output should depend on whether you want to demultiplex or not. Add a tag for that in future.
    basecalled_bam = config.pipeline_steps.basecalling.paths.full_unaligned_bam_path
    demultiplexed_output_dir = config.pipeline_steps.basecalling.paths.full_demultiplexed_output_dir
    kit_name = config.pipeline_steps.basecalling.params.complex_settings.kit_name
    model_speed = config.pipeline_steps.basecalling.params.complex_settings.model_speed
    modifications = config.pipeline_steps.basecalling.params.complex_settings.basecalling_modifications
    batchsize = config.pipeline_steps.basecalling.params.batch_size
    should_demultiplex = config.pipeline_steps.basecalling.params.demultiplex
    dorado_exe = config.tools.dorado

    run_basecalling(dorado_exe, input_file, model_speed, modifications, kit_name, batchsize, should_demultiplex, basecalled_bam)

    if should_demultiplex:
        if is_bam_multiplexed(basecalled_bam):
            logging.info(f"The basecalled BAM {basecalled_bam} is multiplexed and you've told me to demultiplex, so we'll demultiplex it.")
            run_demultiplex(dorado_exe, basecalled_bam, demultiplexed_output_dir)
        else:
            logging.info(f"The basecalled BAM {basecalled_bam} isn't multiplexed, so I won't demultiplex")
    else:
        logging.info("The user specified not to demultiplex.")


def basecall_handler(config):
    input_file = config.pipeline_steps.basecalling.paths.full_pod5_path
    output = config.pipeline_steps.basecalling.paths.full_unaligned_bam_path
    kit_name = config.pipeline_steps.basecalling.params.complex_settings.kit_name
    model_speed = config.pipeline_steps.basecalling.params.complex_settings.model_speed
    modifications = config.pipeline_steps.basecalling.params.complex_settings.basecalling_modifications
    batchsize = config.pipeline_steps.basecalling.params.batch_size
    should_demultiplex = config.pipeline_steps.basecalling.params.demultiplex
    dorado_exe = config.tools.dorado

    run_basecalling(dorado_exe, input_file, model_speed, modifications, kit_name, batchsize, should_demultiplex, output)

    return config.pipeline_steps.basecalling.paths.full_unaligned_bam_path

def demultiplex_handler(config):
    input_file = config.pipeline_steps.basecalling.paths.full_unaligned_bam_path
    output_dir = config.pipeline_steps.basecalling.paths.full_demultiplexed_output_dir
    dorado_exe = config.tools.dorado

    run_demultiplex(dorado_exe, input_file, output_dir)


def download_handler(config):
    model_name = config.pipeline_steps.basecalling.params.explicit_settings.base_model_name
    dorado_exe = config.tools.dorado

    run_model_download(dorado_exe, model_name)


def run_model_download(dorado_exe, model_name):
    # This should be better defined - give the base model name and
    # use the specified modifications to download the relevant modification models

    logger.info(f"Downloading dorado model {model_name}")

    download_cmd = [
        "dorado", "download",
        "--model", str(model_name),
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
        str(input_file)
    ]

    dorado_runner.run(demux_cmd, str(output_dir))


def run_basecalling(dorado_exe, pod5_input, model_speed, modifications, kit_name, batchsize, should_demultiplex, output_path: Path = None):
    basecalling_cmd = ["basecaller",
                       f"{model_speed},{modifications}",
                       str(pod5_input),
                       "--batchsize", str(batchsize)
                       ]

    if should_demultiplex:
        logger.info("We're going to demultiplex, so we'll add the --no-trim tag to keep the barcodes")
        basecalling_cmd.extend(["--no-trim"])

        # We will only use the kit name if we're going to demultiplex. In future do a conditional check for kit_name and demultiplex.
        if kit_name:
            logger.info(f"Kit name '{kit_name}' provided. Adding --kit-name to command")
            basecalling_cmd.extend(["--kit-name", kit_name])
        else:
            logger.info("No kit name provided. The --kit-name argument will be omitted.")

    final_output_arg = []
    if output_path.suffix in ['.bam', '.sam']:
        ensure_dir_exists(output_path.parent)
    elif output_path.suffix == '':
        # If the user provided a directory path
        ensure_dir_exists(output_path)
        final_output_arg = ["-o", str(output_path)]
    else:
        # Invalid input
        raise ValueError(
            f"Invalid output path '{output_path}'. It must be a directory or a .bam/.sam file."
        )

    dorado_runner = ToolRunner(dorado_exe)

    full_cmd = basecalling_cmd + final_output_arg

    if output_path.suffix in ['.bam', '.sam']:
        dorado_runner.run(full_cmd, output_path)
    else:
        dorado_runner.run(full_cmd)


def _add_kit_name_arg(parser, config):
    parser.add_argument(
        "--kit-name", type=str,
        default=config.pipeline_steps.basecalling.params.complex_settings.kit_name,
        dest="pipeline_steps.basecalling.params.complex_settings.kit_name",
        help="Specify the Nanopore kit name"
    )


def _add_model_name_arg(parser, config):
    parser.add_argument(
        '--model-name', type=str,
        default=config.pipeline_steps.basecalling.params.explicit_settings.base_model_name,
        dest="pipeline_steps.basecalling.params.explicit_settings.base_model_name",
        help="Name of the model to download."
    )

def _add_demultiplex_arg(parser, config):
    parser.add_argument(
        '--demultiplex',
        action='store_true',
        default=config.pipeline_steps.basecalling.params.demultiplex,
        help='Demultiplex the basecalled bam file, if necessary.',
        dest='pipeline_steps.basecalling.params.demultiplex'
    )

def add_all_arguments_to_parser(parser, config):
    """
    Publically available function to add arguments from the basecalling step to a given parser.
    It's used by the main run command and doesn't include the I/O args.
    """
    _add_kit_name_arg(parser, config)
    _add_model_name_arg(parser, config)
    _add_demultiplex_arg(parser, config)


def setup_parsers(subparsers, parent_parser, config):
    basecalling_parser = subparsers.add_parser(
        "basecalling",
        help="Run basecalling and related tasks using Dorado.",
        description="This command groups contains tools for basecalling raw Nanopore output.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['basecall'],
        parents=[parent_parser]
    )

    def show_basecalling_help(config):
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
    _add_demultiplex_arg(p_run, config)
    add_io_arguments(
        p_run, config,
        default_input=None,
        input_file_help="Path to POD5 files",
        input_dest="pipeline_steps.basecalling.paths.user_pod5_input",
        default_output=None,
        output_dir_help="Path to basecalled data.",
        output_dest="pipeline_steps.basecalling.paths.user_basecalled_output"
    )
    p_run.set_defaults(func=full_basecalling_handler)

    p_basecall = basecalling_subparsers.add_parser(
        'basecall', help="Basecall a POD5 file/directory.",
        parents=[parent_parser]
    )
    _add_kit_name_arg(p_basecall, config)
    _add_demultiplex_arg(p_basecall, config)
    add_io_arguments(
        p_basecall, config,
        default_input=None,
        input_file_help="Path to POD5 files",
        input_dest="pipeline_steps.basecalling.paths.user_pod5_input",
        default_output=None,
        output_dir_help="Path to basecalled BAM file",
        output_dest="pipeline_steps.basecalling.paths.user_basecalled_output"
    )
    p_basecall.set_defaults(func=basecall_handler)

    p_demux = basecalling_subparsers.add_parser(
        'demux', help="Demultiplex a multiplexed BAM file.",
        parents=[parent_parser]
    )
    _add_kit_name_arg(p_demux, config)
    add_io_arguments(
        p_demux, config,
        default_input=None,
        input_file_help="Path to basecalled BAM file",
        input_dest="pipeline_steps.basecalling.paths.user_basecalled_output",
        default_output=None,
        output_dir_help="Path to demultiplexed data directory.",
        output_dest="pipeline_steps.basecalling.paths.user_demux_output"
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
        default_output=None,
        output_dir_help="Directory to save demultiplexed BAM files to.",
        output_dest="pipeline_steps.basecalling.paths.dorado_model_dir_name"
    )
    p_download.set_defaults(func=download_handler)
