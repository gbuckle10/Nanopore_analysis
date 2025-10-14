#!/usr/bin/env python

import logging
import argparse
from pipeline_controller import PipelineController
import subprocess
import sys


def main():
    """
    Main entry point for the pipeline
    """

    command_map = {
        'setup': 'pipeline',
        'basecalling': 'pipeline',
        'align': 'pipeline',
        'methylation_summary': 'pipeline',
        'deconvolution_prep': 'pipeline',
        'deconvolution': 'pipeline',
        'all': 'pipeline',

        'filter_bam_by_size': 'src/deconvolution_analysis.py'
    }

    alias_map = {
        'basecall': 'basecalling',
        'alignment': 'align',
        'deconv': 'deconvolution'
    }

    if len(sys.argv) < 2:
        user_command = 'all'
    else:
        user_command = sys.argv[1]

    handler = command_map.get(user_command)

    if not handler:
        sys.exit(1)

    if handler == 'pipeline':
        controller = PipelineController()

        controller.main(sys.argv[1:])
    elif handler:
        script_path = handler
        remaining_argv = sys.argv[2:]
        command_to_run = ["python", script_path, user_command] + remaining_argv

        try:
            subprocess.run(command_to_run, check=True),
        except subprocess.CalledProcessError:
            sys.exit(1)


if __name__ == '__main__':
    main()