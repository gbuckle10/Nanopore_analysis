import os
import zipfile
from pathlib import Path
import logging
import gzip
import shutil

logger = logging.getLogger(__name__)


def decompress_file(file_path: Path, delete_original: bool = True):
    """
    Decompresses a .gz or .zip file.
    """
    suffix = file_path.suffix.lower()

    if suffix == '.gz':
        output_path = file_path.with_suffix('')
        logger.info(f"'{file_path.name}' is a .gz file, so we'll unzip it")
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Decompressed to {output_path}")
    elif suffix == '.zip':
        logger.info(f"'{file_path.name}' is a .zip file, so we'll unzip it.")
        output_path = file_path.parent  # .zip files extract into a directory
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(output_path)
        logger.info(f"Extracted contacts of {file_path.name} to {output_path}")
    else:
        logger.warning(f"File {file_path.name} is not a .gz or a .zip file. I haven't done anything.")
        return file_path

    if delete_original:
        file_path.unlink()
        logger.info(f"Deleted original file: {file_path.name}")

    return output_path


def ensure_dir_exists(dir_path: Path, interactive: bool = False) -> bool:
    """
    Checks whether a given directory exists. If in interactive mode, it will prompt the user to create the directory
    if it doesn't already. In non-interactive mode it will create the directory automatically.
    """

    logger.info(f"We are going to check whether the directory {dir_path} exists. Interactive mode? {interactive}")
    if dir_path.is_dir():
        logger.info(f"Path '{dir_path}' exists and is a directory, so we will continue")
        return True
    if dir_path.exists():
        logger.info(f"Error: Path '{dir_path}' exists, but it's not a directory.")
        return False

    if interactive:
        try:
            choice = input(f"Directory '{dir_path}' does not exist. Create it? [y/N] ")
            if choice.lower() != 'y':
                logger.info("Operation cancelled by user.")
                return False

        except (EOFError, KeyboardInterrupt):
            logger.info("\nOperation cancelled.")
            return False

    # Either we're in non-interactive mode or the user said yes to making the directory.
    try:
        logger.info(f"Creating directory: {dir_path}")
        dir_path.mkdir(parents=True, exist_ok=True)
        return True
    except OSError as e:
        logger.error(f"Error: Couldn't create directory {dir_path}. Error - {e}")
        return False
