import argparse
from pathlib import Path

from src.config.models import AppSettings


def add_input_file_argument(parser, default, help_text):
    parser.add_argument(
        "--input-file",
        type=Path,
        default=default,
        help=f"{help_text} (Defaults to value in config file)"
    )


def add_output_dir_argument(parser, default, help_text):
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=default,
        help=f"{help_text} (Defaults to value in config file)"
    )


def add_io_arguments(
        parser: argparse.ArgumentParser,
        config: AppSettings,
        default_input: Path = None,
        default_output: Path = None,
        input_file_help: str = "Path to the input file/directory",
        output_dir_help: str = "Path to the output file/directory"
):
    """
    Adds standard I/O arguments to a given parser with customisable help messages and defaults
    """
    if default_output is None:
        default_output = config.paths.root

    add_input_file_argument(parser, default_input, input_file_help)
    add_output_dir_argument(parser, default_output, output_dir_help)
