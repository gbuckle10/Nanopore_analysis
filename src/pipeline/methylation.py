import argparse
from src.utils.cli_utils import add_io_arguments

from src.utils.process_utils import run_command


def pileup_handler(config):
    aligned_file = config.pipeline_steps.align.paths.full_aligned_bam_path
    output_bed = config.pipeline_steps.methylation.paths.final_bed_file

    run_methylation_pileup(aligned_file, output_bed)


def run_methylation_pileup(aligned_sorted_file, output_bed):
    pileup_cmd = [
        'modkit', 'pileup',
        aligned_sorted_file,
        output_bed
    ]

    run_command(pileup_cmd)


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
        input_file_help="Path to the aligned BAM path",
        input_dest="pipeline_steps.align.paths.aligned_bam_name",
        default_output=None,
        output_dir_help="Path to the BED file",
        output_dest="pipeline_steps.methylation.paths.methylation_bed_name"
    )
    p_pileup.set_defaults(func=pileup_handler)
