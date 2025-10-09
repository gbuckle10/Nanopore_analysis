#!/usr/bin/env python

import subprocess
import sys
import yaml
import logging
import argparse
from collections import deque
from datetime import datetime
from pathlib import Path
from src.utils.runner import run_command, run_wgbstools, run_uxm
from src.utils.logger import setup_logger
from src.utils.file_conversion import apply_runtime_config, ensure_tool_symlink
from src.deconvolution import Deconvolution
from src.deconvolution_prep import generate_deconvolution_files
import os
import re


class Pipeline:
    def __init__(self, config_path):
        """Initialise pipeline, set up logging and load config"""
        run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file_path = f"logs/{run_timestamp}_pipeline_run.log"
        self.logger = setup_logger(log_file=log_file_path)
        self.logger.info(f"PIPELINE RUN STARTED AT: {run_timestamp}")
        self.logger.info(f"The log for this run will be saved to: {log_file_path}")
        self.logger.info("=" * 80)
        self.tool_paths = {}
        self.config = self.load_config(config_path)
        self.steps_to_run = self.config['pipeline_control']['run_steps']
        self.project_root = Path(__file__).resolve().parent

    def load_config(self, config_file="config.yaml"):
        """ Loads the pipeline config from a YAML file """

        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            self.logger.info(f"Configuration loaded from {config_file}")
            return config
        except FileNotFoundError:
            self.logger.critical(f"Configuration file {config_file} not found.")
            raise

    def run_setup(self):
        """ Executes the 00_setup.sh script """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting Step 0: Setup")

        '''
        script_path = "src/00_setup.sh"
        config_file = "config.yaml"
        command = ["bash", script_path, config_file]
        '''

        script_path = "00_setup.py"
        command = [
            sys.executable,
            script_path
        ]

        run_command(command)

        self.tool_paths = apply_runtime_config()

        wgbstool_path = self.project_root / "externals" / "wgbs_tools"
        wgbstools_sl_path = self.project_root / "externals" / "wgbs_tools" / "wgbstools"
        wgbstools_py_path = self.project_root / "externals" / "wgbs_tools" / "src" / "python" / "wgbs_tools.py"
        ensure_tool_symlink(wgbstools_sl_path, wgbstools_py_path)

        self.logger.info("Symlink for wgbstools good.")

        run_wgbstools(['wgbstools', '--version'], self.project_root)

        uxm_sl_path = self.project_root / "externals" / "UXM_deconv" / "uxm"
        uxm_py_path = self.project_root / "externals" / "UXM_deconv" / "src" / "uxm.py"
        ensure_tool_symlink(uxm_sl_path, uxm_py_path)

        self.logger.info("Symlink for uxm good.")

        run_uxm(['uxm', '--help'], self.project_root)

    def run_basecalling(self):
        """ Executes the basecalling script using the 01_basecalling.sh script """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting step 1: Basecalling")

        '''
        script_path = "src/01_basecalling.sh"
        config_file = "config.yaml"
        command = ["bash", script_path, config_file]
        '''

        script_path = "01_basecalling.py"
        command = [
            sys.executable,
            script_path
        ]
        run_command(command)


    def run_alignment(self):
        """ Executes the alignment script: 02_alignment.sh """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting step 2: Alignment and Indexing")

        '''
        script_path = "src/02_alignment.sh"
        config_file = "config.yaml"
        command = ["bash", script_path, config_file]
        '''
        script_path = "02_alignment.py"
        command = [
            sys.executable,
            script_path
        ]

        run_command(command)

    def run_alignment_qc(self):
        self.logger.info("=" * 80)
        self.logger.info(">>> Running QC on aligned and indexed data")

        '''
        script_path = "src/03_alignment_qc.sh"
        config_file = "config.yaml"
        command = ["bash", script_path, config_file]
        '''

        script_path = "03_alignment_qc.py"
        command = [
            sys.executable,
            script_path
        ]

        run_command(command)

    def run_methylation_summary(self):
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting step 4: Methylation Summary")
        '''
        script_path = "src/04_methylation_summary.sh"
        config_file = "config.yaml"
        command = ["bash", script_path, config_file]
        '''

        script_path = "04_methylation_summary.py"
        command = [
            sys.executable,
            script_path
        ]

        run_command(command)

    def run_deconvolution_prep(self):
        """
        Executes the full deconvolution analysis in two steps. The first step is to
        arrange the data in the bed file in such a way that it can be compared to the
        meth_atlas.
        """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting preparation for deconvolution.")

        atlas_file = self.project_root / "data" / "atlas" / "full_atlas_geco.csv"
        file_to_deconv = self.project_root / "data" / "processed" / "deconvolution_geco.csv"
        output_file = self.project_root / "data" / "processed"

        deconv_handler = Deconvolution(self.config, self.tool_paths, atlas_file, file_to_deconv, output_file)

        if "deconvolution_prep" in self.steps_to_run:
            deconv_handler.prepare()
        deconv_handler.run()




    def ___run_deconvolution_submodule(self):
        self.logger.info(">>> Starting the deconvolution process using the meth_atlas submodule")

        atlas_file = f"data/atlas/{self.config['paths']['atlas_file_uxm']}"
        file_to_deconvolve = f"data/processed/{self.config['paths']['file_for_deconvolution_uxm']}"

        self.logger.info(f"Deconvolving file {file_to_deconvolve} using atlas file {atlas_file}")

        deconvolve_script = "externals/meth_atlas/deconvolve_genome_coordinates.py"
        # deconvolve_script = "externals/meth_atlas/deconvolve.py"

        command = [
            "python",
            deconvolve_script,
            "-a",
            atlas_file,
            file_to_deconvolve,
            "--out_dir", self.config['paths']['analysis_dir']
        ]

        self.logger.info(f"--- Running : {' '.join(command)} --- ")

        try:
            with subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    universal_newlines=True
            ) as process:
                if process.stdout:
                    for line in iter(process.stdout.readline, ''):
                        self.logger.info(f"[Subprocess] {line.strip()}")

            if process.returncode != 0:
                self.logger.error(f"Deconvolution script failed with exit code {process.returncode}.")
                raise subprocess.CalledProcessError(process.returncode, process.args)

            self.logger.info("Deconvolution script completed successfully.")
        except FileNotFoundError:
            self.logger.error(f"Command not found: {command[0]}. Make sure python is in PATH")
            raise
        except subprocess.CalledProcessError as e:
            self.logger.critical("Halting process due to deconvolution script failure.")
            raise e

    def run(self):
        self.logger.info(" ---------------- Starting main run ----------------")
        try:
            steps_to_run = self.config['pipeline_control']['run_steps']
        except (FileNotFoundError, KeyError) as e:
            self.logger.error(f"FATAL ERROR: Could not load required configuration.")
            self.logger.error(f"Details: {e}")
            sys.exit(1)

        if not steps_to_run:
            self.logger.error("Error, there aren't any steps for me to run. Check config.yaml.")
            sys.exit(1)

        self.logger.info('--- Pipeline Started ---')
        self.logger.info(f"--- Pipeline will execute the following steps: {steps_to_run} ---")

        if 'setup' in steps_to_run:
            self.run_setup()
        if 'basecalling' in steps_to_run:
            self.run_basecalling()
        if 'align' in steps_to_run:
            self.run_alignment()
        if 'align_qc' in steps_to_run:
            self.run_alignment_qc()
        if 'methylation_summary' in steps_to_run:
            self.run_methylation_summary()
        if 'deconvolution_prep' in steps_to_run:
            self.run_deconvolution_prep()
        if 'deconvolution' in steps_to_run:
            self.run_deconvolution_submodule()

def main():
    CONFIG_FILE = "../config.yaml"

    parser = argparse.ArgumentParser(description="Nanopore analysis pipeline controller.")

    subparsers = parser.add_subparsers(dest='command', help='Pipeline step to run')
    subparsers.required = True
    subparsers.add_parser('setup', help="Run the initial setup step (download tools, data etc).")
    subparsers.add_parser('basecall', help="Run the basecalling step.")
    subparsers.add_parser('align', help="Run the alignment step.")
    subparsers.add_parser('methylation_summary', help="Run the methylation summary step.")
    subparsers.add_parser('deconvolution_prep', help="Run the deconvolution prep.")
    subparsers.add_parser('deconvolution', help="Run the deconvolution step.")

    subparsers.add_parser('all', help="Run all steps enabled in config.yaml")

    args = parser.parse_args()

    controller = Pipeline(CONFIG_FILE)

    if args.command == 'setup':
        controller.run_setup()
    elif args.command == 'basecall':
        controller.run_basecalling()
    elif args.command == 'setup':
        controller.run_setup()
    else:
        print(f"Unknown command: {args.command}")
        parser.print_help()

if __name__ == '__main__':
    CONFIG_FILE = "../config.yaml"
    main()
    #pipeline = Pipeline(CONFIG_FILE)

    try:
        pipeline.run()
    except Exception as e:
        logging.getLogger('pipeline').critical(f"--- PIPELINE HALTED DUE TO UNHANDLED ERROR: {e} ---", exc_info=True)
        sys.exit(1)
