import os
import subprocess
import sys
from typing import List, Optional, Union
from pathlib import Path
from src.utils.logger import logging
from src.utils.process_utils import run_command, LiveDisplayHandler

logger = logging.getLogger(__name__)


class ToolRunner:
    """
    A runner for external command-line tools. Object handles finding the executable and running commands.
    """

    def __init__(self, executable_path: Union[Path, str], output_flag: Optional[str] = None,
                 handler_class=LiveDisplayHandler):
        """
        Initialise the runner.
        Args:
            executable_path: the path to the tool's executable path.
            output_flag: The flag that the tool uses for file output (e.g. '-o', '--output')
        """

        self.handler_class = handler_class
        self.executable_path = Path(executable_path)
        if not self.executable_path.is_file():
            raise FileNotFoundError(f"Tool executable not found at: {self.executable_path}")
        if not os.access(self.executable_path, os.X_OK):
            raise PermissionError(f"File exists but is not executable: {self.executable_path}")
        self.tool_name = self.executable_path.name
        self.output_flag = output_flag

    def _run_and_capture_stdout(self, full_command: List[str], output_path: str):
        logger.info(f"Executing: {' '.join(full_command)} > {output_path}")
        try:
            with open(output_path, 'wb') as f_out:
                subprocess.run(
                    full_command,
                    check=True,
                    stdout=f_out,
                    stderr=sys.stderr.buffer
                )
            logger.info(f"Successfully wrote output to {output_path}")
        except subprocess.CalledProcessError as e:
            logger.critical(f"{self.tool_name} command failed with code {e.returncode}.")
            if os.path.exists(output_path):
                os.remove(output_path)
            raise

    def run(self, args: List[str], output_path: Union[str, Path, None] = None):  # output_path is optional, but must be a string.
        """
        Constructs and runs the tool command.
        """
        full_command = [str(self.executable_path)] + args
        if output_path:
            if self.output_flag:
                full_command.extend([self.output_flag, str(output_path)])
            else:
                logging.debug(f"An output path is given but the ToolRunner has no output flag. The output will be "
                              f"redirected to {output_path}")
                # The tool prints to stdout by default. Redirect it.
                self._run_and_capture_stdout(full_command, output_path)

        run_command(full_command, output_handler_class=self.handler_class)

    def start(self, args: List[str], **kwargs) -> Optional[subprocess.Popen]:
        """
        Starts the command using subprocess.Popen and returns the process object.
        """
        full_command = [str(self.executable_path)] + args

        logger.info(f"Starting process: {' '.join(full_command)} as part of a pipe.")

        try:
            process = subprocess.Popen(full_command, **kwargs)
            return process
        except FileNotFoundError:
            logger.critical(f"Executable not found for command: {full_command[0]}")
            raise
        except Exception as e:
            logger.critical(f"Failed to start process for command:: {' '.join(full_command)}. Error: {e}")
            raise
