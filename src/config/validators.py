from pathlib import Path
from typing import Optional



def validate_path(path: Path, must_exist: bool = True, must_be_file: bool = False, must_be_dir: bool = False, param_name: str = "Path"):
    """
    Validates a given path and raises specific errors on failure.
    """

    if path is None:
        raise ValueError(f"Missing required parameter: {param_name}")

    if must_exist and not path.exists():
        raise FileNotFoundError(f"Path for {param_name} doesn't exist: {path}")

    if path.exists():
        if must_be_file and not path.is_file():
            raise ValueError(f"Path for {param_name} must be a file, but a directory was provided: {path}")

        if must_be_dir and not path.is_dir():
            raise ValueError(f"Path for {param_name} must be a directory, but a file was provided: {path}")

def validate_pod5(path_to_check: Optional[Path], param_name: str):
    """
    Checks a pod5 input path added to config.yaml:
    - Checks that the path is not None
    - If it's a directory, checks that it contains at least 1 .pod5 file.
    - If it's a file, checks that it's a .pod5 file.
    """
    if path_to_check is None:
        raise ValueError(f"Configuration Error for {param_name}: Path wasn't provided or couldn't be built.")
    if not path_to_check.exists():
        raise FileNotFoundError(f"Configuration Error for {param_name}: The specified path doesn't exist: {path_to_check}")
    if path_to_check.is_dir():
        # Check for .pod5 files in the given directory.
        contains_pod5 = any(path_to_check.rglob('*.pod5'))

        if not contains_pod5:
            raise FileNotFoundError(
                f"Configuration Error for {param_name}: The directory exists but doesn't contain any .pod5 files. Path: {path_to_check}"
            )
    elif path_to_check.is_file():
        if path_to_check.suffix != '.pod5':
            raise ValueError(
                f"Configuration Error for {param_name}: The provided input file does not have a '.pod5' extension. Path: {path_to_check}"
            )

    return True


def resolve_path(root_dir: Path, path_from_config: str):
    """
    Checks whether the path specified in config is absolute or relative (an absolute path starts with / and a
    relative path doesn't).
    - If the path is absolute, it returns the path as-provided.
    - If the path is relative, it adds the relative path to the root directory and returns the full path.
    :param root_dir:
    :param path_from_config:
    :return:
    """

    user_path = Path(path_from_config)

    if user_path.is_absolute():
        return user_path
    else:
        return root_dir / user_path
