import subprocess
import sys
import yaml
import logging
from collections import deque
from datetime import datetime
from scripts.utils.runner import load_config, run_command
from scripts.utils.logger import setup_logger
from scripts.deconvolution_prep import generate_deconvolution_files
import os
import re

def load_tool_paths(env_file="runtime.env"):
    logger = logging.getLogger('pipeline')
    paths = {}

    try:
        with open(env_file, 'r') as f:
            for line in f:
                if '=' in line and not line.strip().startswith('#'):
                    key, value = line.strip().split('=', 1)
                    paths[key] = value

            logger.info("Successfully loaded tool paths from runtime.env")
            return paths
    except FileNotFoundError:
        logger.critical(f"Runtime environment file not found.")
        raise

def apply_runtime_config(runtime_config="scripts/runtime_config.sh"):
    """
    Reads a shell script and applies the export PATH commands to the current
    process's environment.
    :param runtime_config:
    :return:
    """

    logger = logging.getLogger('pipeline')
    config_vars = {}
    try:
        with open(runtime_config, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if line.startswith('export '):
                    line = line[7:] # slice off "export"

                if '=' in line:
                    key, value = line.split('=', 1)
                    config_vars[key] = value
                    logger.info(f"Loaded config: {key}={value}")
        return config_vars
    except FileNotFoundError:
        logger.error(f"Runtime config file '{runtime_config}' not found. PATH not updated.")
def run_setup(config):
    """ Executes the 00_setup.sh script """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting Step 0: Setup")

    script_path = "scripts/00_setup.sh"
    config_file = "config.yaml"

    command = ["bash", script_path, config_file]

    run_command(command, config)

    tool_paths = apply_runtime_config()
    wgbs_tools_path = tool_paths.get("WGBSTOOLS_EXE")
    uxm_tools_path = tool_paths.get("UXM_EXE")
    wgbs_command = [wgbs_tools_path, "--version"]
    run_command(wgbs_command, config)
    uxm_command = [uxm_tools_path, '--help']
    run_command(uxm_command, config)

def run_basecalling(config):
    """ Executes the basecalling script using the 01_basecalling.sh script """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting step 1: Basecalling")

    script_path = "scripts/01_basecalling.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_command(command, config)


def run_alignment(config):
    """ Executes the alignment script: 02_alignment.sh """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting step 2: Alignment and Indexing")

    script_path = "scripts/02_alignment.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_command(command, config)


def run_alignment_qc(config):
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Running QC on aligned and indexed data")

    script_path = "scripts/03_alignment_qc.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_command(command, config)


def run_methylation_summary(config):
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting step 4: Methylation Summary")

    script_path = "scripts/04_methylation_summary.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_command(command, config)

def run_deconvolution_prep(config):
    """
    Executes the full deconvolution analysis in two steps. The first step is to
    arrange the data in the bed file in such a way that it can be compared to the
    meth_atlas.
    """

    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting preparation for deconvolution.")


    try:
        # One day change this so that the yaml structure mirrors the file structure. config['paths']['methylation_dir']['methylation_bed_name']
        bed_file_path = config['paths']['methylation_dir'] + config['paths']['methylation_bed_name']
        manifest_file_path = config['paths']['atlas_dir'] + config['paths']['illumina_manifest']
        uxm_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['uxm_atlas_name']
        chunk_size = int(config['parameters']['analysis']['methylation_aggregation_chunksize'])

        command = [
            sys.executable,
            "-u",
            "scripts/deconvolution_prep.py",
            "--bed-file", bed_file_path,
            "--manifest-file", manifest_file_path,
            "--uxm-atlas-file", uxm_atlas_file_path,
            "--chunk-size", str(chunk_size)
        ]

        run_command(command, config, use_conda=False)


    except Exception as e:
        print(f"--- ERROR during deconvolution prep: {e} ---")
        raise

def run_deconvolution_submodule(config):
    logger = logging.getLogger('pipeline')

    logger.info(">>> Starting the deconvolution process using the meth_atlas submodule")

    atlas_file = f"data/atlas/{config['paths']['atlas_file_uxm']}"
    file_to_deconvolve = f"data/processed/{config['paths']['file_for_deconvolution_uxm']}"

    logger.info(f"Deconvolving file {file_to_deconvolve} using atlas file {atlas_file}")

    deconvolve_script = "externals/meth_atlas/deconvolve_genome_coordinates.py"
    # deconvolve_script = "externals/meth_atlas/deconvolve.py"

    command = [
        "python",
        deconvolve_script,
        "-a",
        atlas_file,
        file_to_deconvolve,
        "--out_dir", config['paths']['analysis_dir']
    ]

    logger.info(f"--- Running : {' '.join(command)} --- ")

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
                    logger.info(f"[Subprocess] {line.strip()}")

        if process.returncode != 0:
            logger.error(f"Deconvolution script failed with exit code {process.returncode}.")
            raise subprocess.CalledProcessError(process.returncode, process.args)

        logger.info("Deconvolution script completed successfully.")
    except FileNotFoundError:
        logger.error(f"Command not found: {command[0]}. Make sure python is in PATH")
        raise
    except subprocess.CalledProcessError as e:
        logger.critical("Halting process due to deconvolution script failure.")
        raise e





def main():
    """ Main entry point for the pipeline controller """

    run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file_path = f"logs/{run_timestamp}_pipeline_run.log"
    logger = setup_logger(log_file=log_file_path)

    logger.info(f"PIPELINE RUN STARTED AT: {run_timestamp}")
    logger.info(f"The log for this run will be saved to: {log_file_path}")
    logger.info("=" * 80)

    try:
        config = load_config()
        steps_to_run = config['pipeline_control']['run_steps']
    except (FileNotFoundError, KeyError) as e:
        logger.error(f"FATAL ERROR: Could not load required configuration.")
        logger.error(f"Details: {e}")
        sys.exit(1)

    if not steps_to_run:
        logger.error("Error, there aren't any steps for me to run. Check config.yaml.")
        sys.exit(1)

    logger.info('--- Pipeline Started ---')
    logger.info(f"--- Pipeline will execute the following steps: {steps_to_run} ---")

    if 'setup' in steps_to_run:
        run_setup(config)
    if 'basecalling' in steps_to_run:
        run_basecalling(config)
    if 'align' in steps_to_run:
        run_alignment(config)
    if 'align_qc' in steps_to_run:
        run_alignment_qc(config)
    if 'methylation_summary' in steps_to_run:
        run_methylation_summary(config)
    if 'deconvolution_prep' in steps_to_run:
        run_deconvolution_prep(config)
    if 'deconvolution' in steps_to_run:
        run_deconvolution_submodule(config)


if __name__ == '__main__':
    main()

