import pty
import subprocess
import logging
import yaml
import os
import signal

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

    print(command)
    logger.info(f"Executing command: {' '.join(command)}")

    process = None

    try:
        # create pseudoterminal
        parent_fd, child_fd = pty.openpty()

        # start subprocess, connecting the output to the pseudoterminal.
        process = subprocess.Popen(
            command,
            stdout=child_fd,
            stderr=child_fd,
            text=True,
            preexec_fn=os.setpgrp
        )

        # We don't need the child part of the terminal anymore.
        os.close(child_fd)

        # Read the clean, live output from the parent part of the terminal
        with os.fdopen(parent_fd) as master_file:
            try:
                for line in iter(master_file.readline, ''):
                    logger.info(line.strip())
            except OSError:
                pass
        process.wait()

        if process.returncode != 0:
            logger.error(f"Command failed with exit code {process.returncode}")
            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            logger.info("Command completed successfully")

    except KeyboardInterrupt:
        logger.warning(">>> Keyboard interrupt detected. Terminating subprocess...")
        if process:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        raise
    except Exception as e:
        logger.critical(f"An unexpected error occurred: {e}")
        if process:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        raise e







