from pathlib import Path


def path_must_exist(value: Path) -> Path:
    """
    Pydantic validator to ensure that a path exists
    """
    if not value.exists():
        raise ValueError(f"Path does not exist: {value}")
    return value

def path_must_be_file(value: Path) -> Path:
    """
    Pydantic validator to ensure that a path is a file.
    """
    path_must_exist(value) # Ensure that it exists, first of all.
    if not value.is_file():
        raise ValueError(f"Path must be a file, but it's not: {value}")
    return value

def path_must_be_dir(value: Path) -> Path:
    """
    Pydantic validator to ensure that a path is a directory.
    """
    path_must_exist(value)
    if not value.is_dir():
        raise ValueError(f"Path must be a directory, but it's not: {value}")
    return value

def path_is_file_or_dir(value: Path) -> Path:
    """
    Pydantic validator for paths that can be either a file or a directory.
    """
    path_must_exist(value)
    if not (value.is_file() or value.is_dir()):
        raise ValueError(f"Path is neither a file nor a directory: {value}")
    return value