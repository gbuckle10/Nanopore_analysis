import logging
import colorlog
from colorlog import ColoredFormatter
from pathlib import Path
import re

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
def setup_logger(name='pipeline', log_file='logs/pipeline.log'):
    """
    Sets up a standardised logger with coloured console output and file output
    """
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    console_formatter = ColoredFormatter(
        '%(log_color)s%(levelname)-8s%(reset)s %(blue)s%(message)s',
        log_colors={
            'DEBUG': 'cyan',
            'INFO': 'green',
            'WARNING': 'yellow',
            'ERROR': 'red',
            'CRITICAL': 'red,bg_white',
        }
    )
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file, mode=w)

    file_formatter = AnsiStrippingFormatter(
        '%(asctime)s - %(levelname)-8s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    '''
    # Prevent from propagating to the root logger
    logger.propagate = False

    # Create console handler with colours
    color_handler = colorlog.StreamHandler()
    color_format = colorlog.ColoredFormatter(
        '%(log_color)s%(asctime)s - %(levelname)-8s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        log_colors={
            'DEBUG':    'cyan',
            'INFO':     'green',
            'WARNING':  'yellow',
            'ERROR':    'red',
            'CRITICAL': 'bold_red'
        }
    )

    color_handler.setFormatter(color_format)

    # Create file handler
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_file, mode='w')
    file_format = logging.Formatter(
        '%(asctime)s - %(levelname)-8s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_format)

    if not logger.handlers:
        logger.addHandler(color_handler)
        logger.addHandler(file_handler)
    '''
    return logger