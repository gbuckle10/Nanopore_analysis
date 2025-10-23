import collections
import configparser
import logging
import os
from pathlib import Path

import yaml

def resolve_param(args, config, arg_name, config_path=None, construct_path=False):
    """
    Get a parameter from the command line, if provided. Otherwise, find it at the provided
    config location.
    """

    # If the CLI value is given, then we just return that.
    cli_value = getattr(args, arg_name, None)
    if cli_value is not None:
        return cli_value

    if not config_path:
        return None # No config path to check

    # Look in the config file
    try:
        if construct_path:
            path_components = []
            for path_keys in config_path:
                current_level = config
                for key in path_keys:
                    current_level = current_level[key]
                path_components.append(str(current_level))
            return os.path.join(*path_components)
        else:
            current_level = config
            for key in config_path:
                current_level = current_level[key]
            return current_level

    except (KeyError, TypeError) as e:
        print(f"Warning: Could not resolve '{arg_name}' from config. Missing key: {e}")
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


def get_project_root() -> Path:
    """
    Finds the projects root directory by navigating up from the file's location.
    """

    current_dir = Path(__file__).resolve().parent

    while current_dir != current_dir.parent:
        if (current_dir / "config.yaml").is_file():
            return current_dir
        # Move up one level
        current_dir = current_dir.parent

    raise FileNotFoundError("Could not find project root containing 'config.yaml'.")
