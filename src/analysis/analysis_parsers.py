import argparse
from pathlib import Path

from src.analysis.summarise_lengths import summarise_lengths
from src.config.models import AppSettings
from src.utils.cli_utils import add_io_arguments



def summarise_length_handler(args: argparse.Namespace, config: AppSettings):
    """
    Handler for the 'summarise-lengths' command.
    """
    if args.output_file:
        output_file = args.output_file
    else:
        output_dir = config.paths.results_dir / "length_summaries"
        output_file = output_dir / f"{args.input_file.stem}_read_length_distribution.csv"

    summarise_lengths(input_file=args.input_file, output_file=output_file)

def setup_parsers(subparsers, parent_parser, config):
    """
    Sets up the parser for the non-pipeline 'analysis' commands.
    """
    analysis_parser = subparsers.add_parser(
        "analysis",
        help="Run miscellaneous analysis and utility scripts",
        parents=[parent_parser]
    )
    analysis_subparsers = analysis_parser.add_subparsers(dest='subcommand', help='Available analysis commands')

    p_summarise_lengths = analysis_subparsers.add_parser(
        'summarise-lengths',
        help='Summarise read lengths from a BAM file.'
    )
    p_summarise_lengths.add_argument(
        'input-file',
        type=Path,
        help='File to summarise the lengths of'
    )
    p_summarise_lengths.add_argument(
        '--output-file', '-o',
        type=Path,
        default='.',
        help='File to save the length summary in (csv).'
    )
    p_summarise_lengths.add_argument(
        "-aligned", "--must-be-aligned",
        type=bool,
        help="Only take aligned into account"
    )
    p_summarise_lengths.add_argument(
        "-min", "--min-length",
        type=int,
        default=0,
        help="Minimum length to count"
    )
    p_summarise_lengths.add_argument(
        "-max", "--max-length",
        type=int,
        help="Maximum length to count"
    )

    p_summarise_lengths.set_defaults(func=lambda args, config: summarise_length_handler(args, config))
