import logging
import os
import pty
import re
import signal
import subprocess
import sys
from typing import Callable


logger = logging.getLogger('pipeline')

class LiveLoggingHandler:
    def __init__(self, log_level=logging.INFO):
        self.buffer = ""
        self.log_level = log_level

        self.ansi_escape_pattern = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

    def __call__(self, byte_string: bytes):
        text = byte_string.decode(sys.stdout.encoding, errors='replace')

        for char in text:
            if char == '\r':
                self._flush_buffer(log_level=logging.DEBUG)
                self.buffer = ""
            elif char == '\n':
                self._flush_buffer()
                self.buffer = ""
            else:
                self.buffer += char

    def _flush_buffer(self, log_level=None):
        if not self.buffer:
            return

        clean_line = self.ansi_escape_pattern.sub('', self.buffer).strip()

        if clean_line:
            level_to_use = log_level if log_level is not None else self.log_level
            logger.log(level_to_use, clean_line)

    def flush(self):
        self._flush_buffer()
def kill_process_group(pgid):
    """
    Safely terminates a process group.
    """
    logger = logging.getLogger('pipeline')
    if pgid is None:
        return
    try:
        logger.warning(f"Sending SIGTERM to process group {pgid}")
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        logger.info(f"Process group {pgid} already terminated.")
    except Exception as e:
        logger.error(f"Error during process group termination: {e}")


def raw_print_handler(line: str):
    print(line)
    sys.stdout.write(line)


def log_info_handler(line: str):
    clean_line = line.strip()
    if clean_line:
        logger.info(clean_line)


def run_command(command: list, output_handler: Callable[[str], None] = log_info_handler, output_handler_class=LiveLoggingHandler, env=None):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    print(f"Executing command: {' '.join(command)}")
    logger.info(f"Executing command: {' '.join(command)}")
    process = None
    pgid = None

    handler = output_handler_class()

    try:
        # create pseudoterminal
        parent_fd, child_fd = pty.openpty()

        # start subprocess, connecting the output to the pseudoterminal.
        process = subprocess.Popen(
            command,
            stdout=child_fd,
            stderr=child_fd,
            text=False,
            preexec_fn=os.setpgrp,
            env=env
        )

        # We don't need the child part of the terminal anymore.
        os.close(child_fd)
        pgid = os.getpgid(process.pid)

        # Read the clean, live output from the parent part of the terminal
        with os.fdopen(parent_fd, 'rb', 0) as master_file:
            while True:
                try:
                    chunk = master_file.read(1024)
                    if not chunk:
                        break

                    handler(chunk)
                    '''
            parent_fd = None  # fd is now managed by master_file
            try:
                for line in iter(master_file.readline, ''):
                    handler(line)
                    #output_handler(line)
                    
                    clean_line = line.strip()
                    if not clean_line:
                        continue
                    if '\r' in line:
                        logger.debug(clean_line)
                    else:
                        logger.info(clean_line)
                    
                handler.flush()
                    '''
                except OSError:
                    break

        handler.flush()
        process.wait()

        if process.returncode != 0:
            logger.error(f"Command failed with exit code {process.returncode}")
            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            logger.info("Command completed successfully")

    except KeyboardInterrupt:
        logger.warning(">>> Keyboard interrupt detected. Terminating subprocess...")
        kill_process_group(pgid)
        logger.warning("Process terminated by user. Exiting")
        raise
    except Exception as e:
        logger.critical(f"An unexpected error occurred: {e}")
        kill_process_group(pgid)
        raise
