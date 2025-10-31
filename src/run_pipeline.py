#!/usr/bin/env python
import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from src.config.models import load_and_validate_configs, AppSettings, print_config
from src.config.paths import build_config_paths
from src.pipeline import basecalling, alignment, deconvolution, methylation, full_pipeline
from src.pipeline.full_pipeline import run_full_pipeline
from src.utils import resource_downloader
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



def main():
    global_parent_parser = argparse.ArgumentParser(description="Nanopore Analysis Pipeline", add_help=False)
    global_parent_parser.add_argument(
        '-u', '--user-config',
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help=f"Path to the user config file. Default = {DEFAULT_CONFIG_PATH}"
    )
    global_parent_parser.add_argument(
        '-r', '--runtime-config',
        type=Path,
        default=DEFAULT_RUNTIME_CONFIG_PATH,
        help=f"Path to the runtime config file. Default = {DEFAULT_RUNTIME_CONFIG_PATH}"
    )
    global_parent_parser.add_argument(
        '--debug', action='store_true', help="Enable debug mode. Show full tracebacks on the console."
    )
    global_parent_parser.add_argument(
        '--no-log', action='store_true', help='Disable logging for this run'
    )

    # Use parse_known_args() to only read the arguments that the main_parser knows about - config and user config.
    conf_args, _ = global_parent_parser.parse_known_args()

    # Load config
    try:
        config = load_and_validate_configs(
            conf_args.user_config, conf_args.runtime_config
        )
        build_config_paths(config)
        #print_config(config)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error loading configuration: {e}", file=sys.stderr)
        sys.exit(1)

    # Set up the main parser
    main_parser = argparse.ArgumentParser(
        description="Suite of tools for analysing nanopore data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[global_parent_parser]
    )

    subparsers = main_parser.add_subparsers(dest='command', help='Available command groups')

    full_pipeline.setup_parsers(subparsers, global_parent_parser, config)
    basecalling.setup_parsers(subparsers, global_parent_parser, config)
    alignment.setup_parsers(subparsers, global_parent_parser, config)
    methylation.setup_parsers(subparsers, global_parent_parser, config)
    deconvolution.setup_parsers(subparsers, global_parent_parser, config)
    resource_downloader.setup_parsers(subparsers, global_parent_parser, config)

    try:
        first_arg = sys.argv[1]
    except IndexError:
        first_arg = None # No arguments were provided

    if first_arg is None:
        print("No command was given, so we'll default to run.")
        sys.argv.insert(1, 'run')

    print(f"argv - {sys.argv}")
    # Parse and dispatch
    args = main_parser.parse_args()

    log_level = logging.DEBUG if args.debug else logging.INFO
    # log_level = logging.DEBUG if args.verbose else logging.INFO
    run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_path = f"logs/{run_timestamp}_pipeline_run.log" if not args.no_log else None
    Logger.setup_logger(log_level=log_level, log_file=log_file_path)

    sys.excepthook = handle_exception
    # Call the function that is attached by set_defaults


    if hasattr(args, 'func'):
        args.func(args, config)
    else:
        print(f"ERROR: You must specify a subcommand for '{args.command}'. Use -h for help.", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
