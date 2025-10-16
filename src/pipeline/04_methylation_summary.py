import argparse
import os
import subprocess
import sys
from pathlib import Path

from src.utils.logger import setup_logger
from utils.runner import load_config, get_project_root, run_external_command

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "src", "runtime_config.sh")


def methylation_pileup(config):
    alignment_output_dir = config['paths']['alignment_output_dir']
    aligned_sorted_file = os.path.join(
        project_root,
        alignment_output_dir,
        config['paths']['aligned_bam_name']
    )
    methylation_dir = config['paths']['methylation_dir']
    methylation_bed_name = config['paths']['methylation_bed_name']
    output_bed = os.path.join(
        project_root,
        methylation_dir,
        methylation_bed_name
    )

    pileup_cmd = [
        'modkit', 'pileup',
        aligned_sorted_file,
        output_bed
    ]


def add_args(parser):
    """Add setup-specific args to the parser"""

    parser.add_argument(
        "-c", "--config",
        type=Path,
        help="Path to the config file."
    )

    return parser


def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Methylation summary of sequenced and aligned data using modkit."
    )

    parser = add_args(parser)
    args = parser.parse_args(argv)
    return args


def main(argv=None):

    setup_logger()
    args = parse_args(argv)

    if not args.config:
        config = load_config(CONFIG_PATH)
    else:
        config = load_config(args.config)

    methylation_pileup(config)


if __name__ == '__main__':
    main()
