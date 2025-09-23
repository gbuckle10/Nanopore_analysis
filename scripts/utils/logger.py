import logging
import colorlog

def setup_logger(name='pipeline', log_file='logs/pipeline.log'):
    """
    Sets up a standardised logger with coloured console output and file output
    """
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

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
    file_handler = logging.FileHandler(log_file, mode='w')
    file_format = logging.Formatter(
        '%(asctime)s - %(levelname)-8s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_format)

    if not logger.handlers:
        logger.addHandler(color_handler)
        logger.addHandler(file_handler)

    return logger