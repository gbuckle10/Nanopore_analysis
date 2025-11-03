import importlib.metadata
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
import socket

from src.config.models import AppSettings, RunMetadata


def _parse_dorado_version(version_string: str) -> str:
    """
    Parses the version output format from Dorado
    Split on '+' and take the first part
    """

    version = version_string.split('\n')[1]

    return version


def _parse_samtools_version(version_string: str) -> str:
    """
    Parses the version output format from samtools by taking the second word
    """
    for line in version_string.split('\n'):
        if line.lower().startswith('samtools'):
            return line.split()[1]
    return "unknown"


def _parse_default_version(version_string: str) -> str:
    """
    A generic parser that takes the first word of the first line.
    """
    return version_string.split('\n')[0].split()[0]

# Parser dispatch table, mapping tool names to their parsers
VERSION_PARSERS = {
    'dorado': _parse_dorado_version,
    'samtools': _parse_samtools_version
}   # Add more if you need to

def get_tool_version(executable_path: str) -> str:

    if not executable_path:
        return "Not configured"


    try:

        tool_name = Path(executable_path).name
        print(f"Tool name for exe path {executable_path} is {tool_name}")

        cmd = [executable_path, '--version']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        stdout_clean = result.stdout.strip()
        stderr_clean = result.stderr.strip()

        output_string = stdout_clean if stdout_clean else stderr_clean
        parser_func = VERSION_PARSERS.get(tool_name, _parse_default_version)
        return parser_func(output_string)

    except Exception:
        return "unknown"


def get_git_hash() -> str:
    """
    Gets the current Git commit hash of the repository. If an analysis run is to be replicated exactly,
    this is how to find the correct version to run.
    Returns "unknown" if not in a git repository or if git command fails.
    """
    try:
        repo_root = Path(__file__).parent
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            cwd=repo_root
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        # If we're not in a git repository or if git isn't installed.
        return "unknown"


def build_metadata(config: AppSettings) -> None:
    """
    Populates the metadata section of the object in-place

    """
    print("--- Gathering and building run metadata ---")

    meta = RunMetadata()

    try:
        meta.pipeline_version = importlib.metadata.version("nanopore_analysis")
    except importlib.metadata.PackageNotFoundError:
        meta.pipeline_version = "unknown"

    meta.run_timestamp = datetime.now().isoformat()
    meta.launch_command = ' '.join(sys.argv)
    meta.python_version = sys.version
    meta.hostname = socket.gethostname()
    meta.git_hash = get_git_hash()
    meta.working_directory = os.getcwd()
    meta.dorado_version = get_tool_version(config.tools.dorado)
    meta.samtools_version = get_tool_version("samtools")
    meta.minimap2_version = get_tool_version("minimap2")

    config.metadata = meta
