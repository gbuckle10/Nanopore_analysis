import argparse
import logging

from src.utils import logger
from src.utils.cli_utils import add_io_arguments
from src.utils.file_utils import ensure_dir_exists

from src.utils.process_utils import run_command, AnsiStrippingInteractiveHandler, AnsiPassthroughHandler, \
    SelectiveAnsiInteractiveHandler

logger = logging.getLogger(__name__)

def pileup_handler(config):
    methylation_input_file = config.pipeline_steps.methylation.paths.full_aligned_input_path
    methylation_output_file = config.pipeline_steps.methylation.paths.full_bed_file_path
    reference_fasta = config.pipeline_steps.align.paths.full_ref_fasta_path

    run_methylation_pileup(methylation_input_file, methylation_output_file, reference_fasta)


def run_methylation_pileup(aligned_sorted_file, output_bed, reference_fasta_path):
    logger.debug(f"Ensuring output directory exists: {output_bed.parent}.")
    ensure_dir_exists(output_bed.parent)
    pileup_cmd = [
        'modkit', 'pileup',
        str(aligned_sorted_file),
        str(output_bed)
    ]
    if reference_fasta_path is not None:
        pileup_cmd.extend(['--cpg', '-r', str(reference_fasta_path)])

    run_command(pileup_cmd, output_handler_class=SelectiveAnsiInteractiveHandler)


def setup_parsers(subparsers, parent_parser, config):
    methylation_parser = subparsers.add_parser(
        "methylation-summary",
        help="Run tasks related to summarising the methylation pattern of an aligned BAM file.",
        description="This command groups contains tools for summarising and analysing the methylation in an aligned BAM file.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['methylation'],
        parents=[parent_parser]
    )

    def show_methylation_help(config):
        """Default function to show help for the basecalling command group"""
        methylation_parser.print_help()

    methylation_parser.set_defaults(func=show_methylation_help)

    methylation_subparsers = methylation_parser.add_subparsers(
        title="Available Commands",
        description="Choose one of the following actions to perform.",
        dest='subcommand',
        metavar="<command>"
    )

    p_pileup = methylation_subparsers.add_parser(
        'run', help="Summarise methylation in an aligned BAM file using modkit pileup",
        parents=[parent_parser]
    )
    add_io_arguments(
        p_pileup, config,
        default_input=None,
        input_file_help="Path to the aligned and sorted BAM file.",
        input_dest="pipeline_steps.methylation.paths.user_methylation_input",
        default_output=None,
        output_dir_help="Path to the BED file",
        output_dest="pipeline_steps.methylation.paths.user_methylation_ouput"
    )
    p_pileup.add_argument(
        '--ref',
        help="Path to reference fasta for pileup"
    )
    p_pileup.set_defaults(func=pileup_handler)
