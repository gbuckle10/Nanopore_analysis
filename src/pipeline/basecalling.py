import os
from pathlib import Path
import logging
import os
from pathlib import Path

from src.utils.config_utils import get_project_root
from src.utils.tools_runner import ToolRunner

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "runtime_config.yaml")
logger = logging.getLogger(__name__)


def download_dorado_model(config, model_name):
    # This should be better defined - give the base model name and
    # use the specified modifications to download the relevant modification models

    print(f"Downloading dorado model {model_name}")

    download_cmd = [
        "dorado", "download",
        "--model", model_name,
        "--models-directory", "models/"
    ]
    dorado_exe = config['tools']['dorado']
    dorado_runner = ToolRunner(dorado_exe)
    print(f"Downloading dorado model {model_name} with \n{' '.join(download_cmd)}")

    dorado_runner.run(download_cmd)

    print(f"Dorado model successfully downloaded.")


def demultiplex_bam(config, raw_bam_dir, output_dir: Path):
    raw_bam_filename = config['paths']['unaligned_bam_name']
    raw_bam_file = os.path.join(project_root, raw_bam_dir, raw_bam_filename)

    dorado_exe = config['tools']['dorado']
    dorado_runner = ToolRunner(dorado_exe)

    output_dir.mkdir(parents=True, exist_ok=True)

    demux_cmd = [
        # "dorado",
        "demux",
        "--output-dir", "analysis/demultiplexed",
        "--kit-name", "SQK-NBD114-24",
        raw_bam_file
    ]

    dorado_runner.run(demux_cmd)


def basecalling_pod5(config, pod5_input, kit_name=None, output_file=None):
    '''
    We will need to make this method modifiable, depending on whether we want to use model speed and modifications or
    specific models. This will change soon.
    :param config:
    :param kit_name:
    :return:
    '''

    model_speed = config['parameters']['basecalling']['model_speed']
    modifications = config['parameters']['basecalling']['basecalling_modifications']
    batchsize = config['parameters']['basecalling']['batch_size']

    basecalling_cmd = ["basecaller",
                       f"{model_speed},{modifications}",
                       str(pod5_input),
                       # "--kit-name", kit_name,
                       "--no-trim",
                       "--batchsize", batchsize
                       ]

    dorado_exe = Path(config['tools']['dorado'])
    dorado_runner = ToolRunner(dorado_exe)

    dorado_runner.run(basecalling_cmd, output_file)


def run_basecalling_command(args, config):
    """Handles the logic for pod5 basecalling in full-run and basecall commands"""
    print(f"We will do the full basecalling run.")
    input_pod5 = args.input_file or Path(config['paths']['pod5_dir'])
    print(f"Input pod5 is {input_pod5}")
    kit_name = args.kit_name or Path(config['parameters']['basecalling']['kit_name'])
    if args.output_dir:
        basecalled_output_dir = args.output_dir
        print(f"Basecalled output dir - {basecalled_output_dir}")
    else:
        basecalling_dir = config['paths']['basecalled_output_dir']
        basecalled_filename = config['paths']['unaligned_bam_name']
        basecalled_output_dir = os.path.join(basecalling_dir, basecalled_filename)
        print(f"Basecalled output dir - {basecalled_output_dir}")
    print(f"Basecalling pod5")
    basecalling_pod5(config, input_pod5, kit_name, basecalled_output_dir)

    if args.command == 'full-run':
        multiplexed_bam = Path(config['paths']['basecalled_output_dir'])
        demultiplex_bam(config, multiplexed_bam)


def run_demux_command(args, config):
    """Handles the logic for demux commands"""
    input_file = args.input_file or Path(config['paths']['basecalled_output_dir'])
    output_dir = args.output_dir or Path(config['paths']['demultiplexed_output_dir'])
    demultiplex_bam(config, input_file, output_dir)


def run_download_command(args, config):
    model_to_download = args.model_name or Path(config['parameters']['basecalling']['base_model_name'])

    download_dorado_model(config, model_to_download)


def add_input_file_argument(parser, help_text):
    parser.add_argument(
        "--input-file",
        type=Path,
        help=f"{help_text} (Defaults to value in config file)"
    )


def add_output_dir_argument(parser, help_text):
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        help=f"{help_text} (Defaults to value in config file)"
    )


def _setup_basecalling_parser(subparsers, name, help_text, parent_parser):
    parser = subparsers.add_parser(name, help=help_text, parents=[parent_parser])
    parser.add_argument(
        "--kit-name", type=str,
        help="Specify the Nanopore kit name"
    )
    add_input_file_argument(parser, help_text="Path to the POD5 file/directory for basecalling.")
    add_output_dir_argument(parser, help_text="Path to the basecalled BAM file ")
    parser.set_defaults(func=run_basecalling_command)

def _setup_demux_parser(subparsers, parent_parser):
    parser = subparsers.add_parser("demux", help="Demultiplex a multiplexed BAM file", parents=[parent_parser])
    add_input_file_argument(parser, help_text="Path to the multiplexed BAM file.")
    add_output_dir_argument(parser, help_text="Path to the demultiplexed BAM files")
    parser.set_defaults(func=run_demux_command)

def _setup_download_parser(subparsers, parent_parser):
    parser = subparsers.add_parser(
        'download-model',
        help="Download the dorado model specified in the config.",
        parents=[parent_parser]
    )
    parser.add_argument(
        "--model-name",
        type=str,
        help="Dorado basecalling model to download"
    )
    parser.set_defaults(func=run_download_command)

def register_subparsers(subparsers, parent_parser):
    basecalling_parser = subparsers.add_parser(
        "basecalling",
        help="Run basecalling and related tasks using Dorado."
    )

    basecalling_subparsers = basecalling_parser.add_subparsers(dest='subcommand', required=True)

    _setup_basecalling_parser(basecalling_subparsers, 'full-run', "Basecall and demultiplex", parent_parser)
    _setup_basecalling_parser(basecalling_subparsers, 'basecall', 'Run basecalling only', parent_parser)
    _setup_download_parser(basecalling_subparsers, parent_parser)
    _setup_demux_parser(basecalling_subparsers, parent_parser)

