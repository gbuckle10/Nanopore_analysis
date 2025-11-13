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
FULL_ANSI_REGEX = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
LAYOUT_CONTROL_REGEX = re.compile(r'\x1B\[[0-9;]*[A-DJK]')

class MultiLineAnsiHandler:
    """
    Configurable handler for any tool with a multi-line display. It can be configured to
    handle UIs with stable headers and footers, redrawing only the middle section.
    """
    def __init__(self, header_lines_count: int = 2, footer_lines_count: int = 2, **kwargs):
        self.buffer = b""
        self.is_interactive = sys.stdout.isatty()

        # Set the number of stable header and footer lines.
        self.header_lines_count = header_lines_count
        self.footer_lines_count = footer_lines_count

        # When the header is written, this flag is set to true and we don't write the header anymore
        self.header_captured = False
        # How many lines the UI block was the last time we drew it, so we know how many lines to erase
        self.last_drawn_height = 0
        self.final_log_message = ""

    def __call__(self, byte_chunk: bytes):

        if not self.is_interactive:
            # If we're in an interactive environment, strip all of the ansi and log the updates at the debug level.
            clean_text = FULL_ANSI_REGEX.sub('', byte_chunk.decode(sys.stdout.encoding, 'replace')).strip()
            if clean_text:
                logger.debug(clean_text)
            return # No interactive display needed

        # Add the new data to the buffer
        self.buffer += byte_chunk

        # Only proceed if the buffer contains a newline, signalling the end of a multi-line update from the tool
        if b'\n' in self.buffer:
            # Get the last complete update with rsplit. data contains the full update and self.buffer is left with anything else
            # to be used in the next chunk.
            data, self.buffer = self.buffer.rsplit(b'\n', 1)
            # Pass the complete, stable update to main processing engine.
            self._process_update(data)

    def _process_update(self, full_update_bytes):
        # Prepare the data.
        # Decode from raw bytes to a python string
        full_text = full_update_bytes.decode(sys.stdout.encoding, 'replace')
        # Strip all ANSI to get pure text
        clean_text = FULL_ANSI_REGEX.sub('', full_text).strip()

        lines = clean_text.splitlines()
        if not lines:
            return

        # Decide whether to do the initial draw or a redraw
        if not self.header_captured:
            # First time we draw. Slice the list of lines into header, footer and middle
            header = lines[:self.header_lines_count]
            footer = lines[-self.footer_lines_count:] if self.footer_lines_count > 0 else []
            middle = lines[self.header_lines_count : len(lines) - self.footer_lines_count]

            # Print everything to the terminal the first time.
            if header:
                sys.stdout.write('\n'.join(header) + '\n')
            if middle:
                sys.stdout.write('\n'.join(middle) + '\n')
            if footer:
                sys.stdout.write('\n'.join(footer) + '\n')
            # Flag that we've already done this.
            self.header_captured = True
        else:
            # Only care about dynamic middle and footer sections
            middle = lines[self.header_lines_count : len(lines) - self.footer_lines_count]
            # Calculatehow far up the cursor needs to go to get to the start of the previous dynamic block
            prev_dynamic_height = self.last_drawn_height - self.header_lines_count - self.footer_lines_count
            lines_to_move_up = prev_dynamic_height + self.footer_lines_count

            if lines_to_move_up > 0:
                # ANSI command to move the cursor up N lines.
                sys.stdout.write(f'\x1b[{lines_to_move_up}A')

            # Erase everything from the cursor's current position to the end of the screen
            sys.stdout.write('\x1b[J')

            # Draw the new middle and footer to the terminal
            if middle:
                sys.stdout.write('\n'.join(middle) + '\n')
            if self.footer_lines_count > 0:
                footer = lines[-self.footer_lines_count:]
                sys.stdout.write('\n'.join(footer) + '\n')

        sys.stdout.flush()
        self.last_drawn_height = len(lines)
        self.final_log_message = clean_text

    def flush(self):
        if self.buffer:
            self._process_update(self.buffer)

        if self.is_interactive and self.last_drawn_height > 0:
            sys.stdout.write(f'\x1b[{self.last_drawn_height}A')
            sys.stdout.write('\x1b[J')
            sys.stdout.flush()

        if self.final_log_message:
            logger.info(f"Final status: \n{self.final_log_message}")

class SelectiveAnsiInteractiveHandler:
    def __init__(self, *args, **kwargs):
        self.buffer = b""
        self.is_interactive = sys.stdout.isatty()
        self.last_line_with_colour = ""

    def __call__(self, byte_chunk: bytes):
        line_text = byte_chunk.decode(sys.stdout.encoding, errors='replace')

        partially_cleaned = LAYOUT_CONTROL_REGEX.sub('', line_text)

        if '\r' in partially_cleaned:
            final_message = partially_cleaned.rsplit('\r', 1)[-1]
        else:
            final_message = partially_cleaned

        final_message = final_message.strip()

        if not final_message:
            return

        self.last_line_with_colour = final_message

        if self.is_interactive:
            terminal_width = shutil.get_terminal_size((80, 20)).columns

            reset_code = '\x1b[0m'
            display_line = (final_message + reset_code).ljust(terminal_width + len(reset_code))

            sys.stdout.write(display_line + '\r')
            sys.stdout.flush()
    def flush(self):

        if self.is_interactive:
            terminal_width = shutil.get_terminal_size((80, 20)).columns
            sys.stdout.write(' ' * terminal_width + '\r')
            sys.stdout.flush()

        if self.last_line_with_colour:
            final_clean_log_message = FULL_ANSI_REGEX.sub('', self.last_line_with_colour).strip()
            if final_clean_log_message:
                logger.info(f"Final status: {final_clean_log_message}")


class AnsiPassthroughHandler:
    def __init__(self, *args, **kwargs):
        self.log_buffer = b""

    def __call__(self, byte_chunk: bytes):
        try:
            sys.stdout.buffer.write(byte_chunk)
            sys.stdout.buffer.flush()
        except (IOError, ValueError):
            pass

        self.log_buffer += byte_chunk
        lines = self.log_buffer.splitlines(keepends=True)
        if lines and not lines[-1].endswith((b'\n', b'\r')):
            self.buffer = lines.pop()
        else:
            self.buffer = b""

        for line_bytes in lines:
            self._process_line_for_logging(line_bytes)

    def _process_line_for_logging(self, line_bytes: bytes):
        line_text = line_bytes.decode(sys.stdout.encoding, errors='replace')
        clean_line = FULL_ANSI_REGEX.sub('', line_text).strip()
        if clean_line:
            logger.debug(f"STREAM: {clean_line}")

    def flush(self):
        if self.log_buffer:
            self._process_line_for_logging(self.log_buffer)

class AnsiStrippingInteractiveHandler:
    def __init__(self, *args, **kwargs):
        self.buffer = b""
        self.is_interactive = sys.stdout.isatty()
        self.last_cleaned_line = ""

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
        line_text = line_bytes.decode(sys.stdout.encoding, errors='replace')
        clean_line = FULL_ANSI_REGEX.sub('', line_text).strip()
        if not clean_line:
            return

        self.last_cleaned_line = clean_line

        logger.debug(clean_line)

        if self.is_interactive:
            # Clear the previous line before printing a new, normal log line
            terminal_width = shutil.get_terminal_size((80, 20)).columns
            sys.stdout.write(clean_line.ljust(terminal_width) + '\r')
            sys.stdout.flush()

    def flush(self):
        if self.buffer:
            self._process_line(self.buffer)

        if self.is_interactive:
            terminal_width = shutil.get_terminal_size((80, 20)).columns
            sys.stdout.write(' ' * terminal_width + '\r')
            sys.stdout.flush()

        if self.last_cleaned_line:
            logger.info(f"Final modkit statis: {self.last_cleaned_line}")

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

                    logger.debug(f"DEBUG CHUNK: {repr(chunk)}")

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


