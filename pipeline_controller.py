import subprocess
import sys
import yaml
import logging
from collections import deque
from datetime import datetime
from scripts.deconvolution_prep import *
from scripts.utils.logger import setup_logger


def load_config(config_file="config.yaml"):
    """ Loads the pipeline config from a YAML file """
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def run_and_log(config, command):
    """
    Runs a command inside the conda environment, captures the output in real time and logs it using the
    'pipeline' logger
    """
    # Load the logger
    logger = logging.getLogger("pipeline")

    # Build command with the conda wrapper
    conda_env = config['conda_env_name']
    full_command = ["conda", "run", "-n", conda_env] + command
    logger.info(f"Executing command: {' '.join(full_command)}")

    last_n_lines = deque(maxlen=20)

    try:
        with subprocess.Popen(
                full_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
        ) as process:
            if process.stdout:
                for line in iter(process.stdout.readline, ''):
                    clean_line = line.strip()
                    logger.info(clean_line)
                    last_n_lines.append(clean_line)
        if process.returncode != 0:
            logger.error(f"Bash script failed with exit code {process.returncode}.")
            logger.error("--- Start of error output ---")
            for line in last_n_lines:
                logger.error(line)
            logger.error("--- End of error output ---")

            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            logger.info("Command completed successfully")
    except subprocess.CalledProcessError as e:
        raise e
    except FileNotFoundError:
        logger.critical(f"Command not found: {full_command[0]}. Ensure conda is installed and in PATH")
        raise
    except Exception as e:
        logger.critical(f"An unexpected error occurred while running subprocess: {e}")
        raise


def run_deconvolution_submodule(config):
    logger = logging.getLogger('pipeline')

    logger.info(">>> Starting the deconvolution process using the meth_atlas submodule")

    atlas_file = f"data/atlas/{config['paths']['atlas_file_gc']}"
    file_to_deconvolve = f"data/processed/{config['paths']['file_for_deconvolution_ilmn']}"
    output_file = f"{config['paths']['analysis_dir']}{config['paths']['deconvolution_results']}"

    logger.info(f"Deconvolving file {file_to_deconvolve} using atlas file {atlas_file}")
    logger.info(f"Deconvolution results will be found in {output_file}")

    # deconvolve_script = "externals/meth_atlas/deconvolve.py"
    deconvolve_script = "externals/meth_atlas/deconvolve_genome_coordinates.py"

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
        result = subprocess.run(command, check=True, capture_output=False, text=True)
        if result.stderr:
            logger.error(result.stderr.strip())
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        logger.error(f"--- ERROR IN COMMAND ---")
        logger.error(f"Exit code: {e.returncode}")
        logger.error(f"STDOUT:\n{e.stdout}")
        logger.error(f"STDERR:\n{e.stderr}")
        raise
    except Exception as error:
        logger.error(error)


def run_setup(config):
    """ Executes the 00_setup.sh script """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting Step 0: Setup")

    script_path = "scripts/00_setup.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_and_log(config, command)

    # if not dorado_path or not os.path.exists(dorado_path):
    #    raise FileNotFoundError(f"Setup script failed to return a valid path.")

    # print(f"--- Successfully found Dorado executable at: {dorado_path} ---\n ")

    # return dorado_path


def run_basecalling(config):
    """ Executes the basecalling script using the 01_basecalling.sh script """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting step 1: Basecalling")

    script_path = "scripts/01_basecalling.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_and_log(config, command)


def run_alignment(config):
    """ Executes the alignment script: 02_alignment.sh """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting step 2: Alignment and Indexing")

    script_path = "scripts/02_alignment.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_and_log(config, command)


def run_alignment_qc(config):
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Running QC on aligned and indexed data")

    script_path = "scripts/03_alignment_qc.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_and_log(config, command)


def run_methylation_summary(config):
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting step 4: Methylation Summary")

    script_path = "scripts/04_methylation_summary.sh"
    config_file = "config.yaml"
    command = ["bash", script_path, config_file]

    run_and_log(config, command)


def run_analysis(config):
    """
    Executes the full deconvolution analysis in two steps. The first step is to
    arrange the data in the bed file in such a way that it can be compared to the
    meth_atlas.
    """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 80)
    logger.info(">>> Starting analysis workflow")

    try:

        # One day change this so that the yaml structure mirrors the file structure. config['paths']['methylation_dir']['methylation_bed_name']
        bed_file_path = config['paths']['methylation_dir'] + config['paths']['methylation_bed_name']
        manifest_file_path = config['paths']['atlas_dir'] + config['paths']['illumina_manifest']
        illumina_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['atlas_file_ilmn']
        geco_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['atlas_file_gc']
        uxm_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['uxm_atlas_name']
        deconvolution_path = config['paths']['deconvolution_dir']
        chunk_size = int(config['parameters']['analysis']['methylation_aggregation_chunksize'])

        generate_deconvolution_files(
            bed_file=bed_file_path,
            manifest_file=manifest_file_path,
            uxm_atlas_file=uxm_atlas_file_path,
            chunk_size=chunk_size
        )

        '''
        convert_atlas_to_genome_coordinates(
            output_file=geco_atlas_file_path,
            atlas_file=illumina_atlas_file_path,
            manifest_file=manifest_file_path
        )
        '''

        # format_atlas_file(atlas_file=uxm_atlas_file_path)

        # run_deconvolution_submodule(config)

    except Exception as e:
        print(f"--- ERROR during deconvolution prep: {e} ---")
        raise


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
        run_analysis(config)
    if 'deconvolution' in steps_to_run:
        run_deconvolution_submodule(config)


if __name__ == '__main__':
    main()

