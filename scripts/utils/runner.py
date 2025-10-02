import pty
import subprocess
import logging
import yaml
import os
import signal
import sys
from pathlib import Path

def get_project_root() -> Path:
    """
    Finds the projects root directory by navigating up from the file's location.
    """

    current_dir = Path(__file__).resolve().parent

    while current_dir != current_dir.parent:
        if (current_dir / "config.yaml").is_file():
            return current_dir
        # Move up one level
        current_dir = current_dir.parent

    raise FileNotFoundError("Could not find project root containing 'config.yaml'.")


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

def run_external_command(command: list, cwd=None):
    """
    A simpler version of the runner.run_command method which worker scripts will use to execute
    external tools.
    Prints output to console and exists on failure.
    """

    print(f"Executing: {' '.join(command)}")

    try:
        subprocess.run(
            command,
            cwd=cwd,
            check=True
        )
    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. Make sure it's in your PATH", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}. Aborting.", file=sys.stderr)
        sys.exit(1)

def run_command(command: list, env=None):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    logger = logging.getLogger('pipeline')

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
            preexec_fn=os.setpgrp,
            env=env
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







