#!/usr/bin/env python

import logging
import argparse
from pipeline_controller import PipelineController


def main():
    """
    Main entry point for the pipeline
    """
    parser = argparse.ArgumentParser(description="Nanopore analysis pipeline controller.")

    # --- Sub-parser setup ---
    subparsers = parser.add_subparsers(dest='command', help='Pipeline step to run')
    subparsers.required = False # Subcommand is optional, if not given it'll default to "all".

    subparsers.add_parser('setup', help="Run the initial setup step (download tools, data etc).")
    subparsers.add_parser('basecalling', aliases=['basecall'], help="Run the basecalling step.")
    subparsers.add_parser('align', aliases=['alignment'], help="Run the alignment step.")
    subparsers.add_parser('methylation_summary', help="Run the methylation summary step.")
    subparsers.add_parser('deconvolution_prep', help="Run the deconvolution prep.")
    subparsers.add_parser('deconvolution', aliases=['deconv'], help="Run the deconvolution step.")
    subparsers.add_parser('all', help="Run all steps enabled in config.yaml")

    args = parser.parse_args()

    controller = PipelineController()

    if args.command == 'setup':
        controller.run_setup()
    elif args.command == 'basecalling':
        controller.run_basecalling()
    elif args.command == 'align':
        controller.run_alignment()
        controller.run_alignment_qc()
    elif args.command == "methylation_summary":
        controller.run_methylation_summary()
    elif args.command == "all":
        controller.run_active_steps()
    else:
        print(f"Unknown command: {args.command}")
        parser.print_help()


if __name__ == '__main__':
    main()
