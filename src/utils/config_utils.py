import argparse
import collections
import logging
import os
from functools import reduce
from pathlib import Path
from .config import AppSettings
import yaml

def resolve_combined_path(args: argparse.Namespace, config: AppSettings, arg_name: str, config_path_components: List[str]):
    """
    Resolved a path, prioritising a single CLI override argument.
    If no CLI override is present, it constructs the path from config components.
    """
    cli_path = resolve_param(args, config, arg_name=arg_name)

    # If there is a cli override, just return that value
    if cli_path:
        return Path(cli_path)

    path_parts = []
    for component_path in config_path_components:
        # Get each part of the config ONLY
        part = resolve_param(args, config, arg_name=None, config_path=component_path)
        if part is None:
            # If any part is missing, we can't construct the path.
            return None
        path_parts.append(str(part))

    # Join the parts into a single path.
    if path_parts:
        return Path(*path_parts)

    # If there are no path_parts, return None
    return None

def resolve_param(args: argparse.Namespace, config: AppSettings, arg_name=None, config_path=None):
    """
    Get a parameter from the command line, if provided. Otherwise, find it at the provided
    config location. Values are taken from Pydantic AppSettings objects.
    """

    # If the CLI value is given, then we just return that.
    if arg_name:
        cli_value = getattr(args, arg_name, None)
        if cli_value is not None:
            return cli_value

    # If no CLI value, or no arg_name was given, check the config.
    if config_path:
        try:
            keys = config_path.split('.')
            config_value = reduce(getattr, keys, config)
            return config_value
        except AttributeError:
            # The path didn't exist in the config, but it won't cause an error.
            pass

    # Default to None
    return None

def deep_merge(d1, d2):
    """
    Recursively merges 2 dictionaries. Values from
    d1 will overwrite d2.

    This is used to combine config.yaml and runtime_config.yaml. runtime_config.yaml
    overwrites config.yaml.
    """
    for k, v in d2.items():
        if k in d1 and isinstance(d1[k], dict) and isinstance(v, collections.abc.Mapping):
            deep_merge(d1[k], v)
        else:
            d1[k] = v

    return d1

def load_config(config_file) -> dict:
    """ Loads the pipeline config from a YAML file """
    logger = logging.getLogger('pipeline')
    if not Path(config_file).is_file():
        raise FileNotFoundError(f"Config file not found at: {config_file}")
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)
