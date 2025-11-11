import argparse
from pathlib import Path
from typing import Optional

from src.analysis.filter_bam_by_length import filter_bam_by_length
from src.analysis.summarise_lengths import summarise_lengths
from src.config.models import AppSettings


def summarise_length_handler(args: argparse.Namespace, config: AppSettings):
    """
    Handler for the 'summarise-lengths' command.
    """
    print(args)
    if args.output_dir:
        output_file = args.output_dir
    else:
        output_dir = config.paths.results_dir / "length_summaries"
        output_file = output_dir / f"{args.input_file.stem}_read_length_distribution.csv"



    summarise_lengths(
        input_file=args.input_file,
        output_file=output_file,
        must_be_aligned=args.must_be_aligned,
        min_length=args.min_length,
        max_length=args.max_length if args.max_length is not None else float('inf')
    )


def filter_by_length_handler(args: argparse.Namespace, config: AppSettings):
    """
    Handler for the 'filter_bam_by_length' command.
    """

    if args.output_dir:
        output_directory = args.output_dir
    else:
        output_directory = args.input_bam.parent

    input_file = args.input_file
    cutoff = args.cutoff
    mode = args.mode if args.mode is not None else 'both'

    filter_bam_by_length(
        input_file=input_file,
        size_cutoff=cutoff,
        output_dir=output_directory,
        side_selection=mode
    )


def _add_input_file_arg(parser: argparse.ArgumentParser, help_text: str, nargs: Optional[str] = None):
    """
    Adds a positional input file argument to a parser.
    """

    kwargs = {
        'type': Path,
        'help': help_text
    }
    if nargs is not None:
        kwargs['nargs'] = nargs

    parser.add_argument("input_file", **kwargs)


def _add_size_cutoff(parser):
    parser.add_argument(
        "-c", "--cutoff",
        type=int,
        required=True,
        help="The read length cutoff value."
    )


def _add_output_file_arg(parser, help_text):
    parser.add_argument(
        "--output-dir",
        type=Path,
        metavar="<path>",
        help=help_text
    )


def _add_cutoff_mode(parser):
    parser.add_argument(
        "--mode",
        choices=['above', 'below', 'both'],
        type=str,
        help="Specify which reads to save:\n"
             "  above: save reads with length > cutoff\n"
             "  below: save reads with length <= cutoff\n"
             "  both: save reads above and below the cutoff in separate files"
    )


def _add_must_be_aligned(parser):
    parser.add_argument(
        "-aligned", "--must-be-aligned",
        action="store_true",
        help="Add this flag if you want to only count reads that were mapped to the reference genome."
    )


def _add_minimum_length(parser):
    parser.add_argument(
        "-min", "--min-length",
        type=int,
        metavar="",
        default=0,
        help="Minimum length (in bp) to count"
    )


def _add_maximum_length(parser):
    parser.add_argument(
        "-max", "--max-length",
        type=int,
        metavar="",
        help="Maximum length (in bp) to count"
    )


def _setup_filter_length_parser(subparsers, parent_parser, config, formatter):
    p_filter_by_length = subparsers.add_parser(
        'filter-bam-by-length',
        formatter_class=formatter,
        help="Filter a bam file by a specified length and indexes the output"
    )
    _add_input_file_arg(p_filter_by_length, "Path to the input bam file.")
    _add_size_cutoff(p_filter_by_length)
    _add_output_file_arg(p_filter_by_length,"The directory to save the filtered files in. If no directory "
                                       "is specified, the files will be saved in the same directory as the input file.")
    _add_cutoff_mode(p_filter_by_length)

    p_filter_by_length.set_defaults(func=filter_by_length_handler)


def _setup_summarise_length_parser(subparsers, parent, config, formatter):
    p_summarise_lengths = subparsers.add_parser(
        'summarise-lengths',
        formatter_class=formatter,
        help='Summarise read lengths from a BAM file.'
    )
    _add_input_file_arg(p_summarise_lengths, "File to summarise the lengths of.")
    _add_output_file_arg(p_summarise_lengths, "Path to the output CSV file. If not provided, a default filename will "
                                         "be generated and put into the same directory as the input.")
    _add_must_be_aligned(p_summarise_lengths)
    _add_minimum_length(p_summarise_lengths)
    _add_maximum_length(p_summarise_lengths)

    p_summarise_lengths.set_defaults(func=summarise_length_handler)

def setup_parsers(subparsers, parent_parser, config):
    """
    Sets up the parser for the non-pipeline 'analysis' commands.
    """
    analysis_parser = subparsers.add_parser(
        "analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Run miscellaneous analysis and utility scripts",
        parents=[parent_parser]
    )
    analysis_subparsers = analysis_parser.add_subparsers(dest='subcommand', help='Available analysis commands')
    child_formatter = argparse.ArgumentDefaultsHelpFormatter

    _setup_summarise_length_parser(analysis_subparsers, parent_parser, config, child_formatter)
    _setup_filter_length_parser(analysis_subparsers, parent_parser, config, child_formatter)
