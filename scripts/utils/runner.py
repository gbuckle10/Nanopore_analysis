import subprocess
import logging
import yaml


def load_config(config_file="config.yaml"):
    """ Loads the pipeline config from a YAML file """
    logger = logging.getLogger('pipeline')

    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Configuration loaded from {config_file}")
        return config
    except FileNotFoundError:
        logger.critical(f"Configuration file {config_file} not found.")
        raise


def run_command(command: list, config: dict, use_conda: bool=True):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    logger = logging.getLogger('pipeline')

    if use_conda:
        conda_env = config['conda_env_name']
        full_command = ["conda", "run", "-n", conda_env] + command
    else:
        full_command = command

    logger.info(f"Executing command: {' '.join(full_command)}")

    try:
        with subprocess.Popen(
                full_command,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
        ) as process:
            if process.stdout:
                for line in iter(process.stdout.readline, ''):
                    logger.info(line.strip())
        if process.returncode != 0:
            logger.error(f"Bash script failed with exit code {process.returncode}.")
            logger.error("See the output above for the full error details.")

            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            logger.info("Command completed successfully")
    except subprocess.CalledProcessError as e:
        logger.critical(f"Pipeline halting due to command failure: {' '.join(e.cmd)}")
        raise e
    except FileNotFoundError:
        logger.critical(f"Command not found: {full_command[0]}. Ensure conda is installed and in PATH")
        raise
    except Exception as e:
        logger.critical(f"An unexpected error occurred while running subprocess: {e}")
        raise e



