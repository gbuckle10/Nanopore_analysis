from pathlib import Path


def validate_path(path: Path, must_exist: bool = True, must_be_file: bool = False, must_be_dir: bool = False, param_name: str = "Path"):
    """
    Validates a given path and raises specific errors on failure.
    """

    if path is None:
        raise ValueError(f"Missing required parameter: {param_name}")

    if must_exist and not path.exists():
        raise FileNotFoundError(f"Path for {param_name} doesn't exist: {path}")

    if must_be_file and not path.is_file():
        raise ValueError(f"Path for {param_name} must be a file, but a directory was provided: {path}")

    if must_be_dir and not path.is_dir():
        raise ValueError(f"Path for {param_name} must be a directory, but a file was provided: {path}")
