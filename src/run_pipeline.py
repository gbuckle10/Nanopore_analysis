#!/usr/bin/env python

import logging
import argparse
from pipeline_controller import PipelineController
import subprocess
import sys

COMMAND_MAP = {
    'setup': 'pipeline',
    'basecalling': 'pipeline',
    'basecall': 'pipeline',
    'align': 'pipeline',
    'alignment': 'pipeline',
    'methylation_summary': 'pipeline',
    'deconvolution_prep': 'pipeline',
    'deconvolution': 'pipeline',
    'deconv': 'pipeline',
    'all': 'pipeline',

    'filter-bam-by-length': 'src/analysis/filter_bam_by_length.py',
    'summarise-lengths': 'src/analysis/summarise_lengths.py'
}

def main():

    if len(sys.argv) < 2:
        user_command = 'all'
        args_to_send = ['all']
    else:
        user_command = sys.argv[1]
        args_to_send = sys.argv[1:]

    handler = COMMAND_MAP.get(user_command)

    if not handler:
        print(f"Error: Unknown command '{user_command}'", file=sys.stderr)
        print(f"Available commands:", ", ".join(sorted(COMMAND_MAP.keys())))
        sys.exit(1)

    if handler == 'pipeline':
        script_path = "src/pipeline_controller.py"
        command_to_run = ["python", script_path] + args_to_send

        print(f"Running command {' '.join(command_to_run)}")
        try:
            subprocess.run(command_to_run, check=True)
        except subprocess.CalledProcessError:
            print(f"There was an error...")
            sys.exit(1)
    elif handler:
        script_path = handler
        remaining_argv = sys.argv[2:]

        command_to_run = ["python", script_path] + remaining_argv

        print(f"Running command {' '.join(command_to_run)}")
        try:
            subprocess.run(command_to_run, check=True)
        except subprocess.CalledProcessError:
            sys.exit(1)


if __name__ == '__main__':
    main()