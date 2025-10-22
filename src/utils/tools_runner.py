import os
import subprocess
import sys
from typing import List, Optional, Union
from pathlib import Path
from src.utils.logger import logging
from src.utils.process_utils import run_command, raw_print_handler

logger = logging.getLogger(__name__)

class ToolRunner:
    """
    A runner for external command-line tools. Object handles finding the executable and running commands.
    """

    def __init__(self, executable_path: Union[Path, str]):
        self.executable_path = Path(executable_path)
        if not self.executable_path.is_file():
            raise FileNotFoundError(f"Tool executable not found at: {self.executable_path}")
        if not os.access(self.executable_path, os.X_OK):
            raise PermissionError(f"File exists but is not executable: {self.executable_path}")
        self.tool_name = self.executable_path.name

    def run(self, args: List[str], stdout_file: Optional[str] = None):  # stdout_file is optional, but must be a string.
        """
        Constructs and runs the tool command.
        """
        full_command = [str(self.executable_path)] + args

        if stdout_file:
            logger.info(f"Executing: {' '.join(full_command)} > {stdout_file}")
            try:
                with open(stdout_file, 'wb') as f_out:
                    subprocess.run(
                        full_command,
                        check=True,
                        stdout=f_out,
                        stderr=sys.stderr.buffer
                    )
                logger.info(f"Successfully wrote output to {stdout_file}")
            except subprocess.CalledProcessError as e:
                logger.critical(f"{self.tool_name} command failed with code {e.returncode}.")
                if os.path.exists(stdout_file):
                    os.remove(stdout_file)
                raise
        else:
            run_command(full_command, output_handler=raw_print_handler)

    def start(self, args: List[str], **kwargs) -> Optional[subprocess.Popen]:
        """
        Starts the command using subprocess.Popen and returns the process object.
        """
        full_command = [str(self.executable_path)] + args

        logger.info(f"Starting process: {' '.join(full_command)}")

        try:
            process = subprocess.Popen(full_command, **kwargs)
            return process
        except FileNotFoundError:
            logger.critical(f"Executable not found for command: {full_command[0]}")
            raise
        except Exception as e:
            logger.critical(f"Failed to start process for command:: {' '.join(full_command)}. Error: {e}")
            raise