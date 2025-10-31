#!/usr/bin/env python
import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from src.config.models import load_and_validate_configs, AppSettings, print_config
from src.config.paths import build_config_paths
from src.pipeline import basecalling, alignment, deconvolution, methylation
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


def run_full_pipeline(args, config: AppSettings):
    logging.info("--- Running full pipeline from config ---")

    function_map = {
        'basecalling': basecalling.full_basecalling_handler,
        'align': alignment.full_alignment_handler,
        'methylation': methylation.pileup_handler,
        'deconvolution': deconvolution.deconvolution_handler
    }

    steps_to_run = config.pipeline_control.run_steps
    active_steps = [step_name for step_name, should_run in steps_to_run if should_run]

    if not active_steps:
        logging.warning("WARNING: No active steps found in the final configuration. Nothing to do.")

    logging.info("STEPS TO RUN: ")
    for i, step_name in enumerate(active_steps, 1):
        logging.info(f"    Step {i}: {step_name}")
    logging.info("----------------------------------")

    for step_name in active_steps:
        logging.info(f">>> EXECUTING STEP: {step_name}")
        step_func = function_map.get(step_name)
        if not step_func:
            logging.warning(f"WARNING: No function found for step '{step_name}'. Skipping")
            continue

        step_func(args, config)


def main():
    conf_parser = argparse.ArgumentParser(description="Nanopore Analysis Pipeline", add_help=False)
    conf_parser.add_argument(
        '-u', '--user-config',
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help=f"Path to the user config file. Default = {DEFAULT_CONFIG_PATH}"
    )
    conf_parser.add_argument(
        '-r', '--runtime-config',
        type=Path,
        default=DEFAULT_RUNTIME_CONFIG_PATH,
        help=f"Path to the runtime config file. Default = {DEFAULT_RUNTIME_CONFIG_PATH}"
    )
    # Use parse_known_args() to only read the arguments that the main_parser knows about - config and user config.
    conf_args, _ = conf_parser.parse_known_args()

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

    global_opts_parser = argparse.ArgumentParser(add_help=False)

    global_opts_parser.add_argument(
        '--debug', action='store_true', help="Enable debug mode. Show full tracebacks on the console."
    )
    global_opts_parser.add_argument(
        '--no-log', action='store_true', help='Disable logging for this run'
    )

    main_parser = argparse.ArgumentParser(
        description="Suite of tools for analysing nanopore data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[conf_parser, global_opts_parser]
    )

    subparsers = main_parser.add_subparsers(dest='command', help='Available command groups')

    # Register commands from modules.
    run_parser = subparsers.add_parser(
        'run',
        help="Run the full pipeline using steps defined in the config file.",
        parents=[global_opts_parser]
    )
    # Add the arguments from each individual step to the run_parser
    basecalling.add_all_arguments_to_parser(run_parser, config)
    alignment.add_all_arguments_to_parser(run_parser, config)
    run_parser.set_defaults(func=run_full_pipeline)


    basecalling.setup_parsers(subparsers, global_opts_parser, config)
    alignment.setup_parsers(subparsers, global_opts_parser, config)
    methylation.setup_parsers(subparsers, global_opts_parser, config)
    deconvolution.setup_parsers(subparsers, global_opts_parser, config)
    resource_downloader.setup_parsers(subparsers, global_opts_parser, config)

    # Parse and dispatch
    args = main_parser.parse_args()
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
            args.func(args, config)
        else:
            print(f"ERROR: You must specify a subcommand for '{args.command}'. Use -h for help.", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
