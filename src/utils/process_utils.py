import logging
import os
import pty
import re
import shutil
import signal
import subprocess
import sys
from typing import Callable


logger = logging.getLogger('pipeline')

class LiveDisplayHandler:
    '''
    A class that processes a stream of bytes from a subprocess. It displays progress bars on the terminal by
    overwriting the current line, and prints a clean format to the logging file.
    '''

    def __init__(self, log_level=logging.INFO):
        self.buffer = ""
        self.is_interactive = sys.stdout.isatty()
        self.log_level = log_level


    def __call__(self, byte_string: bytes):
        text = byte_string.decode(sys.stdout.encoding, errors='replace')
        self.buffer += text

        while '\n' in self.buffer or '\r' in self.buffer:
            idx_n = self.buffer.find('\n')
            idx_r = self.buffer.find('\r')

            # Which terminator comes first?
            if idx_n == -1: idx_n = float('inf')
            if idx_r == -1: idx_r = float('inf')

            split_idx = min(idx_n, idx_r)

            line = self.buffer[:split_idx]
            terminator = self.buffer[split_idx]
            self.buffer = self.buffer[split_idx + 1:]

            if terminator == '\n':
                self._handle_newline(line)
            elif terminator == '\r':
                self._handle_carriage_return(line)
    def _handle_newline(self, line: str):
        clean_line = line.strip()

        if self.is_interactive:
            # Clear the previous line before printing a new, normal log line
            terminal_width = shutil.get_terminal_size((80, 20)).columns
            sys.stdout.write('\r' + ' ' * terminal_width + '\r')
            sys.stdout.flush()

        if clean_line:
            logger.info(clean_line)

    def _handle_carriage_return(self, line: str):
        clean_line = line.strip()

        if clean_line:
            if self.is_interactive:
                terminal_width = shutil.get_terminal_size((80, 20)).columns
                sys.stdout.write(clean_line.ljust(terminal_width) + '\r')
                sys.stdout.flush()
        else:
            logger.debug(clean_line)

    def flush(self):
        if self.buffer:
            self._handle_newline(self.buffer)

        if self.is_interactive:
            #sys.stdout.write('\n')
            sys.stdout.flush()

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


def run_command(command: list, output_handler: Callable[[str], None] = log_info_handler, output_handler_class=LiveDisplayHandler, env=None, cwd=None):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

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
            env=env,
            cwd=cwd
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

                    #print(f"DEBUG CHUNK: {repr(chunk)}")

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


