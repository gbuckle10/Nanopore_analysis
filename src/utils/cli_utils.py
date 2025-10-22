import argparse
from pathlib import Path


def add_input_file_argument(parser, help_text="Path to the input file/directory"):
    parser.add_argument(
        "--input-file",
        type=Path,
        help=f"{help_text} (Defaults to value in config file)"
    )


def add_output_dir_argument(parser, help_text="Path to the output file/directory"):
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        help=f"{help_text} (Defaults to value in config file)"
    )

def create_io_parser():
    """
    Creates and returns a parent parser with standard I/O arguments.
    """
    io_parser = argparse.ArgumentParser(add_help=False)
    add_input_file_argument(io_parser)
    add_output_dir_argument(io_parser)
    return io_parser
