#!/usr/bin/env python
import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

from src import pipeline_controller
from src.pipeline import basecalling
from src.utils.config_utils import get_project_root, load_config, deep_merge

from src.utils.decorators import graceful_exit

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "runtime_config.yaml")
logger = logging.getLogger(__name__)

COMMAND_MAP = {
    'setup': 'pipeline',
    'basecalling': 'pipeline',
    'basecall': 'pipeline',
    'align': 'pipeline',
    'alignment': 'pipeline',
    'methylation_summary': 'pipeline',
    'deconvolution_prep': 'pipeline',
    'deconvolution': 'pipeline',
    'deconv': 'pipeline',
    'all': 'pipeline',

    'filter-bam-by-length': 'src/analysis/filter_bam_by_length.py',
    'summarise-lengths': 'src/analysis/summarise_lengths.py'
}


def main():
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
        description="Suite of tools for analysing nanopore data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Available command groups')

    # Register commands from modules.
    basecalling.register_subparsers(subparsers, parent_parser)

    # Parse and dispatch
    args = parser.parse_args()

    print(f"Loading user config: {args.user_config}")
    print(f"Loading runtime config: {args.runtime_config}")
    user_config = load_config(args.user_config)
    runtime_config = load_config(args.runtime_config)
    config = deep_merge(user_config, runtime_config)

    # Call the function that iss attached by set_defaults
    if hasattr(args, 'func'):
        print(f"Running function {args.func}")
        args.func(args, config)
    else:
        print("Error: You must specify a subcommand", file=sys.stderr)
        sys.exit(1)

    '''
    if len(sys.argv) < 2:
        user_command = 'all'
        args_to_send = ['all']
    else:
        user_command = sys.argv[1]
        args_to_send = sys.argv[1:]

    handler = COMMAND_MAP.get(user_command)

    if not handler:
        print(f"Error: Unknown command '{user_command}'", file=sys.stderr)
        print(f"Available commands:", ", ".join(sorted(COMMAND_MAP.keys())))
        sys.exit(1)

    if handler == 'pipeline':
        script_path = "src/pipeline_controller.py"
        command_to_run = ["python", script_path] + args_to_send

        try:
            subprocess.run(command_to_run, check=True)
        except subprocess.CalledProcessError:
            sys.exit(1)
    elif handler:
        script_path = handler
        remaining_argv = sys.argv[2:]

        command_to_run = ["python", script_path] + remaining_argv

        print(f"Running command {' '.join(command_to_run)}")
        try:
            subprocess.run(command_to_run, check=True)
        except subprocess.CalledProcessError:
            sys.exit(1)

    '''
if __name__ == '__main__':
    main()