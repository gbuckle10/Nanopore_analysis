import argparse
import os
import sys
from pathlib import Path
from src.utils.logger import setup_logger
from utils.runner import load_config, get_project_root, run_external_command, run_dorado

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "src", "runtime_config.sh")


def download_dorado_model(config):
    # This should be better defined - give the base model name and
    # use the specified modifications to download the relevant modification models

    dorado_model_name = config['parameters']['basecalling']['base_model_name']
    print(f"Downloading dorado model {dorado_model_name}")

    download_cmd = [
        "dorado", "download",
        "--model", dorado_model_name,
        "--models-directory", "models/"
    ]

    print(f"Downloading dorado model {dorado_model_name} with \n{' '.join(download_cmd)}")
    run_dorado(download_cmd, project_root)
    # run_external_command(download_cmd)

    print(f"Dorado model successfully downloaded.")


def demultiplex_bam(config):
    raw_bam_dir = config['paths']['basecalled_output_dir']
    raw_bam_filename = config['paths']['unaligned_bam_name']
    raw_bam_file = os.path.join(project_root, raw_bam_dir, raw_bam_filename)

    demux_cmd = [
        "dorado",
        "demux",
        "--output-dir", "analysis/demultiplexed"
                        "--kit-name", "SQK-NBD114-24",
        raw_bam_file
    ]


def basecalling_pod5(config, kit_name=None):
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

    basecalling_cmd = [
        "dorado", "basecaller",
        f"{model_speed},{modifications}",
        pod5_dir,
        # "--kit-name", kit_name,
        "--no-trim",
        "--batchsize", batchsize
    ]

    print(f"Executing command {' '.join(basecalling_cmd)}")

    run_dorado(basecalling_cmd, project_root, "data/basecalled_output/calls.bam")


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    setup_logger()

    parser = argparse.ArgumentParser(
        description="Basecalling using dorado."
    )

    parser.add_argument(
        "--basecalling-pod5",
        action='store_true',
        help='Just run the basecalling method.'
    )
    parser.add_argument(
        "--download-dorado-model",
        action='store_true',
        help='Download the dorado model.'
    )
    parser.add_argument(
        "--demultiplex",
        action='store_true',
        help='Demultiplex a multiplexed bam file.'
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to the input file/directory."
    )
    parser.add_argument(
        '-c', '--config',
        type=Path,
        help="Path to the config file."
    )

    args = parser.parse_args()

    if not args.config:
        config = load_config(CONFIG_PATH)
    else:
        config = load_config(args.config)

    method_flagged = args.basecalling_pod5 or args.download_dorado_model or args.demultiplex

    print("Yes" if method_flagged else "No")

    if not method_flagged:
        print(f"No specific methods flagged. Running all of them")
        # Add an arg to decide whether to use a specific model, which model to use, and how
        # exactly we want to define the basecalling method.
        download_dorado_model(config)
        basecalling_pod5(config)
    else:
        print("Methods have been specified, so I'll run the ones you specified.")
        if args.basecalling_pod5:
            basecalling_pod5(config)
        if args.download_dorado_model:
            download_dorado_model(config)
        if args.demultiplex:
            demultiplex_bam(config)


if __name__ == "__main__":
    main()
