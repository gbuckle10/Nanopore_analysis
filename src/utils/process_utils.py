import logging
import os
import pty
import signal
import subprocess

from src.utils.runner import logger


def kill_process_group(pgid):
    """
    Safely terminates a process group.
    """

    if pgid is None:
        return
    try:
        logger.warning(f"Sending SIGTERM to process group {pgid}")
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        logger.info(f"Process group {pgid} already terminated.")
    except Exception as e:
        logger.error(f"Error during process group termination: {e}")


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
        pgid = os.getpgid(process.pid)
        # Read the clean, live output from the parent part of the terminal
        with os.fdopen(parent_fd) as master_file:
            parent_fd = None  # fd is now managed by master_file
            try:
                for line in iter(master_file.readline, ''):
                    clean_line = line.strip()
                    if not clean_line:
                        continue
                    if '\r' in line:
                        logger.debug(clean_line)
                    else:
                        logger.info(clean_line)
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
        kill_process_group(pgid)
        raise
    except Exception as e:
        logger.critical(f"An unexpected error occurred: {e}")
        kill_process_group(pgid)
        raise e

    finally:
        if parent_fd is not None:
            os.close(parent_fd)
