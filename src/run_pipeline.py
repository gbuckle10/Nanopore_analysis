#!/usr/bin/env python
import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from src.pipeline import basecalling, alignment, deconvolution
from src.utils import resource_downloader
from src.utils.config_utils import load_config, deep_merge, resolve_param
from src.utils.logger import Logger
from src import PROJECT_ROOT

DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config.yaml"
DEFAULT_RUNTIME_CONFIG_PATH = PROJECT_ROOT / "runtime_config.yaml"

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
    'summarise-lengths': 'src/analysis/summarise_lengths.py',
    'download': 'src/utils/resource_downloader.py'
}

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.critical("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))

    # Print the user-friendly message to the console
    print("\nFATAL ERROR: The application has crashed unexpectedly.", file=sys.stderr)
    print(f"Error Type: {exc_type.__name__}", file=sys.stderr)
    print("Please see the log file for the full technical traceback.", file=sys.stderr)
def run_full_pipeline(args, config):
    logging.info("--- Running full pipeline from config ---")

    function_map = {
        'basecalling': basecalling.full_basecalling_handler,
        'align': alignment.full_alignment_handler
    }

    active_steps = resolve_param(
        args, config, config_path=['pipeline_control', 'run_steps']
    )

    if not active_steps:
        logging.warning("WARNING: No active steps found in config.yaml. Nothing for me to do.")
        return

    for step_name in active_steps:
        logging.info(f">>> EXECUTING STEP: {step_name}")
        step_func = function_map.get(step_name)
        if not step_func:
            logging.warning(f"WARNING: No function found for step '{step_name}'. Skipping")
            continue
        step_func(args, config)


def main():
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '-u', '--user-config',
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help=f"Path to the user config file. Default = {DEFAULT_CONFIG_PATH}"
    )
    parent_parser.add_argument(
        '-r', '--runtime-config',
        type=Path,
        default=DEFAULT_RUNTIME_CONFIG_PATH,
        help=f"Path to the runtime config file. Default = {DEFAULT_RUNTIME_CONFIG_PATH}"
    )
    parent_parser.add_argument(
        '--debug', action='store_true', help="Enable debug mode. Show full tracebacks on the console."
    )
    parent_parser.add_argument(
        '--no-log', action='store_true', help='Disable logging for this run'
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
    deconvolution.setup_parsers(subparsers, parent_parser)
    resource_downloader.setup_parsers(subparsers, parent_parser)

    # Parse and dispatch
    args = parser.parse_args()

    # Load config
    user_config = load_config(args.user_config)
    runtime_config = load_config(args.runtime_config)
    config = deep_merge(user_config, runtime_config)

    log_level = logging.DEBUG if args.debug else logging.INFO
    # log_level = logging.DEBUG if args.verbose else logging.INFO
    run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_path = f"logs/{run_timestamp}_pipeline_run.log" if not args.no_log else None
    Logger.setup_logger(log_level=log_level, log_file=log_file_path)

    sys.excepthook = handle_exception
    # Call the function that is attached by set_defaults

    if args.command is None:
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
