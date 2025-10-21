import argparse
import os
import sys
from pathlib import Path
import logging
from src.utils.logger import setup_logger
from src.utils.config_utils import load_config, get_project_root, deep_merge
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


def demultiplex_bam(config, raw_bam_dir):
    raw_bam_filename = config['paths']['unaligned_bam_name']
    raw_bam_file = os.path.join(project_root, raw_bam_dir, raw_bam_filename)

    dorado_exe = config['tools']['dorado']
    dorado_runner = ToolRunner(dorado_exe)

    demux_cmd = [
        "dorado",
        "demux",
        "--output-dir", "analysis/demultiplexed"
                        "--kit-name", "SQK-NBD114-24",
        raw_bam_file
    ]

    dorado_runner.run(demux_cmd)


def basecalling_pod5(config, pod5_input, kit_name=None):
    '''
    We will need to make this method modifiable, depending on whether we want to use model speed and modifications or
    specific models. This will change soon.
    :param config:
    :param kit_name:
    :return:
    '''
    pod5_dir = config['paths']['pod5_dir']
    basecalling_dir = config['paths']['basecalled_output_dir']
    basecalled_filename = config['paths']['unaligned_bam_name']
    model_speed = config['parameters']['basecalling']['model_speed']
    modifications = config['parameters']['basecalling']['basecalling_modifications']
    batchsize = config['parameters']['basecalling']['batch_size']

    output_file = os.path.join(basecalling_dir, basecalled_filename)

    basecalling_cmd = ["basecaller",
                       f"{model_speed},{modifications}",
                       pod5_dir,
                       # "--kit-name", kit_name,
                       "--no-trim",
                       "--batchsize", batchsize
                       ]

    print(f"Executing command {' '.join(basecalling_cmd)}")

    dorado_exe = Path(config['tools']['dorado'])
    dorado_runner = ToolRunner(dorado_exe)

    dorado_runner.run(basecalling_cmd)

def add_input_file_argument(parser, help_text):
    """
    A helper function to add the --input-file argument to a parser.
    """
    parser.add_argument(
        "--input-file",
        type=Path,
        help=help_text
    )

def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '-u', '--user-config',
        type=Path,
        default=CONFIG_PATH,
        help=f"Path to the user config file. Default = {CONFIG_PATH}"
    )
    parent_parser.add_argument(
        '-r', '--runtime-config',
        type=Path,
        default=RUNTIME_CONFIG_PATH,
        help=f"Path to the runtime config file. Default = {RUNTIME_CONFIG_PATH}"
    )

    parser = argparse.ArgumentParser(
        description="Basecalling using dorado.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    parser_full = subparsers.add_parser(
        'full-run',
        help="Basecall a pod5 file/directory and demultiplex if necessary.",
        parents=[parent_parser]
    )
    parser_full.add_argument(
        "--input-file",
        type=Path,
        help="Path to the pod5 file/directory for basecalling (optional, uses config if not set."
    )
    parser_full.add_argument(
        "--kit-name",
        type=str,
        help="Specify the Nanopore analyser kit name (overrides config)."
    )

    basecall_parser = subparsers.add_parser(
        'run',
        help='Run basecalling on a POD5 directory.',
        parents=[parent_parser]
    )
    basecall_parser.add_argument(
        "--kit-name",
        type=str,
        help="Specify the Nanopore analyser kit name (overrides config)."
    )
    add_input_file_argument(basecall_parser, help_text="Path to the pod5 file/directory for basecalling.")

    download_parser = subparsers.add_parser(
        'download-model',
        help="Download the dorado model specified in the config.",
        parents=[parent_parser]
    )
    download_parser.add_argument(
        "--model-name",
        type=str,
        help="Dorado basecalling model to download"
    )

    demux_parser = subparsers.add_parser(
        'demux',
        help='Demultiplex a multiplexed BAM file',
        parents=[parent_parser]
    )
    add_input_file_argument(demux_parser, help_text="Path to the demultiplexed BAM file.")

    if not any(cmd in subparsers.choices for cmd in argv):
        print("Info: No command specified, defaulting to 'full-run'")
        argv.insert(0, 'full-run')

    print(f"Parsing arguments {argv}")
    args = parser.parse_args(argv)

    print(f"Parsed args {args}")

    print(args)

    return args


def main(argv=None):
    setup_logger()
    args = parse_args(argv)

    command_to_run = args.command

    print(f"Command to run - {command_to_run}")

    user_config = load_config(args.user_config)
    runtime_config = load_config(args.runtime_config)

    config = deep_merge(user_config, runtime_config)

    if command_to_run in ['full-run', 'run']:
        print(f"We will do the full basecalling run.")

        if args.input_file:
            input_pod5 = args.input_file
            print(f"Using --input-file from command line: {input_pod5}")
        else:
            input_pod5 = Path(config['paths']['pod5_dir'])
            print(f"No input file specified on command line, using from config: {input_pod5}")

        if args.kit_name:
            kit_name = args.kit_name
            print(f"Using --kit-name from command line: {kit_name}")
        else:
            kit_name = Path(config['parameters']['basecalling']['kit_name'])
            print(f"No kit specified on command line, using from config: {kit_name}")
        basecalling_pod5(config, input_pod5, kit_name)

        if command_to_run == 'full-run':
            demultiplex_bam(config)

    elif command_to_run == 'demux':
        demultiplex_bam(config, args.input_file)
    elif command_to_run == 'download-model':
        if args.dorado_model:
            model_to_download = args.model_name
        else:
            model_to_download = Path(config['parameters']['basecalling']['base_model_name'])
        download_dorado_model(config, model_to_download)


if __name__ == "__main__":
    main()
