import pandas as pd
import os
from pathlib import Path
import logging

def ensure_tool_symlink(link_path, target_path):
    logger = logging.getLogger('pipeline')

    logger.info(f"Ensuring '{link_path}' is a symlink to '{target_path}'")

    try:
        if link_path.exists():
            if link_path.is_symlink():
                if os.path.realpath(link_path) == os.path.realpath(target_path):
                    logger.info(f"Correct symlink for {link_path} already exists, so nothing to do.")
                else:
                    logger.info(f"Symlink exists, but points to the wrong target. Recreating.")
                    link_path.unlink()
                    link_path.symlink_to(target_path)
                    print("Symlink recreated successfully.")
            else:
                logger.info(f"'{link_path}' is a regular file, not a symlink. Removing and replacing.")
                link_path.unlink() # remove file
                link_path.symlink_to(target_path) # create the symlink
                logger.info(f"Symlink created successfully.")
        else:
            # Directly create the link
            logger.info(f"'{link_path}' does not exist. Creating symlink.")
            link_path.symlink_to(target_path)
            print("Symlink created successfully.")
    except Exception as e:
        logger.error(f"An error occurred: {e}")

def apply_runtime_config(runtime_config="src/runtime_config.sh"):
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
