import argparse
import os
import subprocess
import sys
from pathlib import Path

from src.utils.cli_utils import create_io_parser
from src.utils.config_utils import resolve_param

from src.utils.process_utils import run_command


def pileup_handler(args, config):
    aligned_file = resolve_param(
        args, config, arg_name="input_file", construct_path=True,
        config_path=[
            ['paths', 'alignment_output_dir'],
            ['paths', 'aligned_bam_name']
        ]
    )

    output_bed = resolve_param(
        args, config, arg_name="output_dir", construct_path=True,
        config_path=[
            ['paths', 'methylation_dir'],
            ['paths', 'methylation_bed_name']
        ]
    )

    run_methylation_pileup(aligned_file, output_bed)


def run_methylation_pileup(aligned_sorted_file, output_bed):
    pileup_cmd = [
        'modkit', 'pileup',
        aligned_sorted_file,
        output_bed
    ]

    run_command(pileup_cmd)


def setup_parsers(subparsers, parent_parser):
    io_parser = create_io_parser()
    methylation_parser = subparsers.add_parser(
        "methylation-summary",
        help="Run tasks related to summarising the methylation pattern of an aligned BAM file.",
        description="This command groups contains tools for summarising and analysing the methylation in an aligned BAM file.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['methylation'],
        parents=[parent_parser, io_parser]
    )

    def show_methylation_help(args, config):
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
        parents=[parent_parser, io_parser]
    )
    p_pileup.set_defaults(func=pileup_handler)
