import argparse
import subprocess
from pathlib import Path

from src import PROJECT_ROOT
from src.utils.process_utils import run_command
import os
import sys
import requests
import tarfile
import yaml

def link_tools_to_conda_bin(executable_paths):
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if not conda_prefix:
        print("WARNING: CONDA_PREFIX not set, skipping tool linking.")
        return

    bin_dir = os.path.join(conda_prefix, 'bin')
    tool_names = {'uxm_exe': 'uxm', 'wgbstools_exe': 'wgbstools'}

    for key, name in tool_names.items():
        exe_path = executable_paths.get(key)
        if not exe_path:
            continue
        os.chmod(exe_path, 0o755)
        link_path = os.path.join(bin_dir, name)
        if os.path.islink(link_path) or os.path.exists(link_path):
            os.remove(link_path)
        os.symlink(exe_path, link_path)
        print(f"Linked '{name}' -> {exe_path}")

def download_file(url, destination):
    print(f"Downloading from {url} to {destination}")

    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(destination, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}", file=sys.stderr)
        sys.exit(1)

def setup_submodules(config):
    print(" --- Setting up git submodules. ---")
    if not os.path.exists(f"{PROJECT_ROOT}/.gitmodules"):
        print("No gitmodules file found. Skipping submodule setup.")
        return

    repository_path = PROJECT_ROOT
    executable_paths = {}

    print(f">>> Configuring Git to trust {repository_path}")
    subprocess.run(
        ['git', 'config', '--global', '--add', 'safe.directory', '*'], check=True
    )

    sync_command = ["git", "submodule", "sync", "--recursive"]
    update_command = ["git", "submodule", "update", "--init", "--recursive", "--force"]

    print(f"Running sync command {sync_command}")
    run_command(sync_command)
    run_command(update_command)

    print("Compiling wgbstools")

    comp_command = ["python", "install.py"]
    # run_external_command(comp_command, cwd=wgbstools_dir)

    submodule_paths_config = config.get('pipeline_steps', {}).get('analysis', {}).get('tools') or {}

    uxm_rel = submodule_paths_config.get('uxm_dir', 'externals/UXM_deconv')
    wgbstools_rel = submodule_paths_config.get('wgbstools_dir', 'externals/wgbs_tools')
    methatlas_rel = submodule_paths_config.get('meth_atlas_dir', 'externals/meth_atlas')

    uxm_abs = PROJECT_ROOT / uxm_rel / "uxm"
    wgbstools_abs = PROJECT_ROOT / wgbstools_rel / "wgbstools"
    methatlas_abs = PROJECT_ROOT / methatlas_rel

    executable_paths['uxm_exe'] = str(uxm_abs)
    executable_paths['wgbstools_exe'] = str(wgbstools_abs)
    executable_paths['methatlas_exe'] = str(methatlas_abs)

    return executable_paths


def install_dorado(config):
    """
    Downloads and extracts the correct version of Dorado.
    """
    print(" --- Setting up Dorado ---")
    version = config['pipeline_steps']['setup']['params']['dorado_version']
    archive_filename = f"dorado-{version}-linux-x64.tar.gz"
    download_url = f"https://cdn.oxfordnanoportal.com/software/analysis/dorado-{version}-linux-x64.tar.gz"

    dorado_dir = f"{PROJECT_ROOT}/tools/dorado-{version}-linux-x64"

    print(f"Checking for dorado at {dorado_dir}")

    if os.path.isdir(dorado_dir):
        print(f"Dorado already found at {dorado_dir}. Writing to config.yaml and skipping download.")
        dorado_exe_path = os.path.abspath(os.path.join(dorado_dir, "bin", "dorado"))
        return {'dorado': dorado_exe_path}
    else:
        print(f"Downloading dorado version {version} from {download_url}")
        archive_path = os.path.join(PROJECT_ROOT, "tools", archive_filename)
        os.makedirs(os.path.dirname(archive_path), exist_ok=True)
        download_file(download_url, archive_path)

        print(f"Extracting {archive_path}...")
        with tarfile.open(archive_path, "r:gz") as tar:
            tar.extractall(path="tools")

        os.remove(archive_path)
        print("Extraction complete")

    dorado_exe_path = os.path.abspath(os.path.join(dorado_dir, "bin", "dorado"))
    print(f"The dorado executable path is {dorado_exe_path}")

    return {'dorado': dorado_exe_path}


def add_args(parser):
    """Add setup-specific args to the parser"""
    parser.add_argument(
        "--dorado-version",
        type=str,
        help="Dorado version to download."
    )
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        type=Path,
        help="Path to the config file."
    )
    subparsers = parser.add_subparsers(dest='command', help="Installation task to run")
    parser.set_defaults(command='all') # If no command is given, default to all.

    all_parser = subparsers.add_parser('all', help="Run all installation steps (default)")
    all_parser.add_argument('--dorado-version', help='Override Dorado version from config file.')

    tool_parser = subparsers.add_parser('tools', help="Install/update command-line tools (e.g. Dorado)")
    tool_parser.add_argument('--dorado-version', help='Override Dorado version from config file.')

    subparsers.add_parser('submodules', help="Initialise/update Git submodules.")

    return parser


def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description="Setup for pipeline run"
    )
    parser = add_args(parser)
    args = parser.parse_args(argv)
    return args


def main(argv=None):
    args = parse_args(argv)

    print(f">>> Loading configurations")

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    # If dorado version is in args, override config
    if args.dorado_version is not None:
        print(f"Overriding config with user-defined dorado version '{args.dorado_version}'")
        config['pipeline_steps']['setup']['params']['dorado_version'] = args.dorado_version

    # Load existing runtime_config, otherwise make a new one.
    try:
        with open('runtime_config.yaml', 'r') as f:
            runtime_config = yaml.safe_load(f)
    except FileNotFoundError:
        runtime_config = {}

    if args.command in ['all', 'tools']:
        tool_paths = install_dorado(config)
        runtime_config.setdefault('tools', {}).update(tool_paths)

    if args.command in ['all', 'submodules']:
        submodule_paths = setup_submodules(config)
        link_tools_to_conda_bin(submodule_paths)
        runtime_config.setdefault('pipeline_steps', {}).setdefault('analysis', {}).setdefault('tools', {}).update(submodule_paths)
    print(">>> Writing updated runtime_config.yaml...")
    with open('runtime_config.yaml', 'w') as f:
        yaml.dump(runtime_config, f, sort_keys=False)


if __name__ == "__main__":
    main()
