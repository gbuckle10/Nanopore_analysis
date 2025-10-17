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

def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    commands = {'full-run', 'run', 'download-model', 'demux'}

    user_provided_command = False
    for arg in argv:
        if arg in commands:
            user_provided_command = True
            break

    if not user_provided_command:
        argv.insert(0, 'full-run')
        print(f"INFO: No command specified, so defaulting to 'full-run. Final args = {argv}")

    parser = argparse.ArgumentParser(
        description="Basecalling using dorado.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-u', '--user-config',
        type=Path,
        default=CONFIG_PATH,
        help=f"Path to the user config file. Default = {CONFIG_PATH}"
    )
    parser.add_argument(
        '-r', '--runtime-config',
        type=Path,
        default=RUNTIME_CONFIG_PATH,
        help=f"Path to the runtime config file. Default = {RUNTIME_CONFIG_PATH}"
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    parser_full = subparsers.add_parser(
        'full-run',
        help="Basecall a pod5 file/directory and demultiplex if necessary."
    )
    parser_full.add_argument(
        "--input-file",
        type=Path,
        help="Path to the POD5 directory/file."
    )

    run_parser = subparsers.add_parser(
        'run',
        help='Run basecalling on a POD5 directory.'
    )
    run_parser.add_argument(
        "--input-file",
        type=Path,
        help="Path to the input file/directory."
    )

    download_parser = subparsers.add_parser(
        'download-model',
        help="Download the dorado model specified in the config."
    )
    download_parser.add_argument(
        "--model-name",
        type=str,
        help="Dorado basecalling model to download"
    )

    demux_parser = subparsers.add_parser(
        'demux',
        help='Demultiplex a multiplexed BAM file'
    )
    demux_parser.add_argument(
        "--input-file",
        type=Path,
        help="Path to the input BAM file to be demultiplexed."
    )

    run_parser.add_argument(
        "--demultiplex",
        action='store_true',
        help='Demultiplex a multiplexed bam file.'
    )

    args = parser.parse_args(argv)
    return args


def main(argv=None):
    setup_logger()
    args = parse_args(argv)

    user_config = load_config(args.user_config)
    runtime_config = load_config(args.runtime_config)

    config = deep_merge(user_config, runtime_config)

    command_to_run = args.command
    if command_to_run is None:
        command_to_run = 'full-run'
        logger.info("No step specified, so I'll just do the basecalling and demux")

    print(f"Command to run - {command_to_run}")

    if command_to_run in ['full-run', 'run']:
        print(f"We will do the full basecalling run.")
        if getattr(args, 'input-file', None):
            input_pod5 = args.input_file
        else:
            input_pod5 = Path(config['paths']['pod5_dir'])

        if getattr(args, 'kit-name', None):
            kit_name = args.kit_name
        else:
            kit_name = Path(config['parameters']['basecalling']['kit_name'])

        basecalling_pod5(config, input_pod5, kit_name)

        if command_to_run == 'full-run':
            demultiplex_bam(config)

    elif command_to_run == 'demux':
        if args.input:
            input_bam_dir = args.input
        else:
            input_bam_dir = Path(config['paths']['basecalled_output_dir'])

        demultiplex_bam(config, input_bam_dir)
    elif command_to_run == 'download-model':
        if args.dorado_model:
            model_to_download = args.dorado_model
        else:
            model_to_download = Path(config['parameters']['basecalling']['base_model_name'])
        download_dorado_model(config, model_to_download)

if __name__ == "__main__":
    main()
