#!/usr/bin/env python
import argparse
import logging
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from src.pipeline import basecalling, alignment
from src.utils import logger
from src.utils.config_utils import get_project_root, load_config, deep_merge

from src.utils.decorators import graceful_exit
from src.utils.logger import Logger

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "runtime_config.yaml")

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


def run_full_pipeline(args, config):
    logging.info("--- Running full pipeline from config ---")

    function_map = {
        'basecalling': basecalling.run_basecalling_command
        # Add the rest of the commands here as you go.
        # 'align': alignment.run_alignment_command
    }

    active_steps = config['pipeline_control']['run_steps']
    if not active_steps:
        logging.warning("WARNING: No active steps found in config.yaml. Nothing for me to do.")
        return

    for step_name in active_steps:
        logging.info(f">>> EXECUTING STEP: {step_name}")
        step_func = function_map.get(step_name)
        if not step_func:
            logging.warning(f"WARNING: No function found for step '{step_name}'. Skipping")
            continue


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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[parent_parser]
    )


    subparsers = parser.add_subparsers(dest='command', help='Available command groups')

    # Register commands from modules.
    run_parser = subparsers.add_parser(
        'run',
        help="Run the full pipeline using steps defined in the config file.",
        parents=[parent_parser]
    )
    run_parser.set_defaults(func=run_full_pipeline)

    # subparsers.add_parser('basecalling', aliases=['basecall'])
    basecalling.setup_parsers(subparsers, parent_parser)
    alignment.setup_parsers(subparsers, parent_parser)

    # Parse and dispatch
    args = parser.parse_args()

    # Load config
    user_config = load_config(args.user_config)
    runtime_config = load_config(args.runtime_config)
    config = deep_merge(user_config, runtime_config)

    # log_level = logging.DEBUG if args.verbose else logging.INFO
    run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_path = f"logs/{run_timestamp}_pipeline_run.log"
    Logger.setup_logger(log_level=logging.INFO, log_file=log_file_path)

    # Call the function that is attached by set_defaults
    if args.command is None:
        print("INFO: No command specified")
        run_full_pipeline(args, config)
    else:
        if hasattr(args, 'func'):
            logging.info(f"Running argument {args.func}")
            args.func(args, config)
        else:
            print(f"ERROR: You must specify a subcommand for '{args.command}'. Use -h for help.", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
