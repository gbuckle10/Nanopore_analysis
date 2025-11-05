import argparse
from pathlib import Path
from typing import Optional

from src.config.models import AppSettings


def add_input_file_argument(parser, default, input_dest, help_text):
    parser.add_argument(
        "-i", "--input-file",
        type=Path,
        default=default,
        dest=input_dest,
        metavar="<path>",
        help=f"{help_text} (Defaults to value in config file)"
    )


def add_output_dir_argument(parser, default, output_dest, help_text):
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=default,
        dest=output_dest,
        metavar="<path>",
        help=f"{help_text} (Defaults to value in config file)"
    )


def add_io_arguments(
        parser: argparse.ArgumentParser,
        config: AppSettings,
        add_input: bool = True,
        default_input: Path = None,
        input_file_help: str = "Path to the input file/directory",
        input_dest: Optional[str] = None,
        add_output: bool = True,
        default_output: Path = None,
        output_dir_help: str = "Path to the output file/directory",
        output_dest: Optional[str] = None
):
    """
    Adds standard I/O arguments to a given parser with customisable help messages and defaults
    """
    if add_input:
        if not input_dest:
            raise ValueError("'input_dest' must be provided when 'add_input' is True.")
        add_input_file_argument(
            parser, default_input, input_dest, f"{input_file_help}. (default: %(default)s)"
        )

    if add_output:
        if not output_dest:
            raise ValueError("'output_dest' must be provided when 'add_output' is True.")
        add_output_dir_argument(
            parser, default_output, output_dest, f"{output_dir_help}. (default: %(default)s)"
        )
