import pandas as pd
import os
import pybedtools
import logging

def apply_runtime_config(runtime_config="scripts/runtime_config.sh"):
    """
    Reads a shell script and applies the export PATH commands to the current
    process's environment.
    :param runtime_config:
    :return:
    """

    logger = logging.getLogger('pipeline')
    config_vars = {}
    try:
        with open(runtime_config, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if line.startswith('export '):
                    line = line[7:]  # slice off "export"

                if '=' in line:
                    key, value = line.split('=', 1)
                    config_vars[key] = value
                    logger.info(f"Loaded config: {key}={value}")
        return config_vars
    except FileNotFoundError:
        logger.error(f"Runtime config file '{runtime_config}' not found. PATH not updated.")
