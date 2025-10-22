import logging
import colorlog
from colorlog import ColoredFormatter
from pathlib import Path
import re
import os

class AnsiStrippingFormatter(logging.Formatter):
    """
    A logging formatter that strips ANSI escape codes
    """

    ANSI_ESCAPE_REGEX = re.compile(r'\x1b\[[0-9;]*[mK]')

    def format(self, record):
        # Initial formatting
        formatted_record = super().format(record)

        # Strip ANSI codes from the result.
        return self.ANSI_ESCAPE_REGEX.sub('', formatted_record)

class Logger:

    @staticmethod
    def setup_logger(log_level=logging.INFO, log_file=None):
        """
        Sets up a standardised logger with coloured console output and file output.
        - If the logger is already configured, it does nothing.
        - If log_file is provided, it logs to a file in addition to the console.
        - If log_file is not provided, it logs only to the console.
        """
        # Create logger
        root_logger = logging.getLogger()
        root_logger.setLevel(log_level)

        # If the logger already has handlers, it means that it's already been configured.
        # Therefore, don't add more handlers.
        if root_logger.hasHandlers():
            root_logger.handlers.clear()

        console_formatter = ColoredFormatter(
            '%(log_color)s%(levelname)-8s - %(message)s%(reset)s',
            log_colors={
                'DEBUG': 'cyan',
                'INFO': 'green',
                'WARNING': 'yellow',
                'ERROR': 'red',
                'CRITICAL': 'red,bg_white',
            }
        )

        console_handler = logging.StreamHandler()
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)

        if log_file:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)

            file_handler = logging.FileHandler(log_file, mode='w')

            file_formatter = AnsiStrippingFormatter(
                '%(asctime)s - %(levelname)-8s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            file_handler.setFormatter(file_formatter)
            root_logger.addHandler(file_handler)
