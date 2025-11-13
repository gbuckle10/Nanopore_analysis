import logging
import os
import pty
import re
import shutil
import signal
import subprocess
import sys
from pathlib import Path
from typing import Callable, List, Union

logger = logging.getLogger(__name__)

class TunableHandler:
    """
    A tunable handler which suppresses the main data stream but actively watches for specific, important messages
    like file skipping warnings and records their occurrences in its own internal state.
    """
    def __init__(self, log_prefixes: List[str] = None, skip_signal_regex: str = None, *args, **kwargs):
        self.buffer = b""
        self.skip_signal_detected = False
        if log_prefixes is None:
            self.log_prefixes = [b'[', b'Error', b'Warning', b'info']
        else:
            self.log_prefixes = [p.encode() for p in log_prefixes]

        if skip_signal_regex:
            self.skip_signal_regex = re.compile(skip_signal_regex.encode(), re.IGNORECASE)
        else:
            self.skip_signal_regex = None

    def __call__(self, byte_chunk: bytes):
        self.buffer += byte_chunk
        lines = self.buffer.splitlines(keepends=True)
        if lines and not lines[-1].endswith((b'\n', b'\r')):
            self.buffer = lines.pop()
        else:
            self.buffer = b""
        for line_bytes in lines:
            self._process_line(line_bytes)

    def _process_line(self, line_bytes: bytes):
        stripped_line = line_bytes.strip()
        if not stripped_line:
            return

        # Check for skip signal, if configured
        if self.skip_signal_regex and self.skip_signal_regex.search(stripped_line):
            self.skip_signal_detected = True
            line_text = stripped_line.decode(sys.stdout.encoding, 'replace')
            logging.warning(f"Stateful signal detected: {line_text}")
            return

        # Check for other generic log messages
        if any(stripped_line.startswith(p) for p in self.log_prefixes):
            line_text = stripped_line.decode(sys.stdout.encoding, 'replace')
            logger.info(line_text)
            return

        logging.debug("DATA_STREAM (suppressed): %s", stripped_line[:200])

    def flush(self):
        if self.buffer:
            self._process_line(self.buffer)
            self.buffer = b""

class SilentHandler:
    """
    A handler class that does nothing with the output from a subprocess. It silences the stderr/stdout scream while
    still allowing the main execution wrapper to monitor the process for exit codes.
    """
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, byte_chunk: bytes):
        # Don't do anything with the incoming data
        pass

    def flush(self):
        pass

class LineBufferingHandler:
    """
    A simple handler for tools that produce standard text output. It does no
    interactive terminal manipulation.
    """
    def __init__(self, *args, **kwargs):
        self.buffer = b''

    def __call__(self, byte_chunk: bytes):
        """Buffers bytes and processes complete lines."""
        self.buffer += byte_chunk
        # Use splitlines to handle \n \r and \r\n
        lines = self.buffer.splitlines(keepends=True)
        if lines and not lines[-1].endswith((b'\n', b'\r')):
            self.buffer = lines.pop()
        else:
            self.buffer = b""

        for line_bytes in lines:
            self._process_line(line_bytes)

    def _process_line(self, line_bytes: bytes):
        """Decodes and logs a single, complete line."""
        line_text = line_bytes.decode(sys.stdout.encoding, errors='replace').strip()
        if line_text:
            logger.info(line_text)

    def flush(self):
        if self.buffer:
            self._process_line(self.buffer)
            self.buffer = b""

class LiveDisplayHandler:
    '''
    A handler class that processes a stream of bytes from a subprocess. It displays progress bars on the terminal by
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
        pass
    except Exception as e:
        logger.error(f"Error during process group termination: {e}")

def run_command(command: list, output_handler_class=LiveDisplayHandler,
                stdout_redirect_path: Union[str, Path, None] = None,
                env=None, cwd=None, **handler_kwargs):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    str_command = list(map(str, command))
    logger.info(f"Executing command: {' '.join(command)}")
    if cwd:
        logger.debug(f"Running command in directory: {cwd}")
    if stdout_redirect_path:
        stdout_redirect_path = str(stdout_redirect_path)
        logger.debug(f"Redirecting command stdout to: {stdout_redirect_path}")

    process = None
    pgid = None
    f_out = None


    try:
        if stdout_redirect_path:
            f_out = open(stdout_redirect_path, 'wb')

        # create pseudoterminal
        parent_fd, child_fd = pty.openpty()
        stdout_target = f_out if f_out else child_fd

        # start subprocess, connecting the output to the pseudoterminal.
        process = subprocess.Popen(
            command,
            stdout=stdout_target,
            stderr=child_fd,
            text=False,
            preexec_fn=os.setpgrp,
            env=env,
            cwd=cwd
        )

        # We don't need the child part of the terminal anymore.
        os.close(child_fd)
        pgid = os.getpgid(process.pid)

        handler = output_handler_class()

        # Read the clean, live output from the parent part of the terminal
        with os.fdopen(parent_fd, 'rb', 0) as master_file:
            while True:
                try:
                    chunk = master_file.read(1024)
                    if not chunk:
                        break

                    #logging.debug(f"RAW_CHUNK: {repr(chunk)}")

                    handler(chunk)

                except OSError:
                    break

        parent_fd = -1
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


