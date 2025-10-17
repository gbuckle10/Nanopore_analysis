#!/usr/bin/env python
import os
import sys
import subprocess
import argparse

DOCKER_WRAPPER_SCRIPT = """#!/bin/bash
# This script runs the nanopore pipeline inside its Docker container.
IMAGE_NAME="nanopore-analysis:latest"
docker run --rm -it \\
  -v "$(pwd)/data:/app/data" \\
  -v "$(pwd)/output:/app/output" \\
  "$IMAGE_NAME" "$@"
"""


def run_command(command, description):
    print(f">>> {description}...")
    try:
        subprocess.run(command, check=True)
        print(">>> Success!")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"ERROR during '{description}': {e}", file=sys.stderr)
        sys.exit(1)


def install_docker():
    """
    Builds the docker image and creates the wrapper script.
    """
    run_command(["docker", "build", "-t", "nanopore-pipeline:latest", "."], "Building Docker image")

    print(">>> Creating Docker wrapper script: ./nanopore_analysis")
    with open("nanopore_analysis", "w") as f:
        f.write(DOCKER_WRAPPER_SCRIPT)
    os.chmod("nanopore_analysis", 0o755)  # Give it rwxr-xr-x permissions
    print(">>> Docker setup complete")


def install_conda(task):
    """
    Creates the Conda environment and the local symlink.
    """
    env_name = "nanopore_analysis"

    if task == 'all':
        run_command([
            "mamba", "env", "create", "-f", "environment.yml"
        ], "Creating conda environment and running internal setup step.")

    print(">>> Installing Conda activation scripts")
    conda_info = subprocess.check_output(["conda", "info", "--json"], text=True)
    import json
    conda_root = json.loads(conda_info)['conda_prefix']
    env_path = os.path.join(conda_root, 'envs', env_name)
    activation_dir = os.path.join(env_path, 'etc', 'conda', 'activate.d')
    os.makedirs(activation_dir, exist_ok=True)
    import shutil
    shutil.copy(
        './etc/conda/activate.d/env_vars.sh',
        os.path.join(activation_dir, 'env_vars.sh')
    )

    command_to_run = [
        "mamba", "run", "-n", env_name, "python", "-m", "src.setup_logic", task
    ]
    run_command(command_to_run, f"Running step logic for '{task}'")

    print(">>> Creating local symlink: ./nanopore_analysis")
    if os.path.exists("nanopore_analysis"):
        os.remove("nanopore_analysis")
    os.symlink("src/run_pipeline.py", "nanopore_analysis")
    os.chmod("src/run_pipeline.py", 0o755)
    print(">>> Local setup complete.")


def add_args(parser):
    """
    Define command-line arguments for the installer
    """
    subparsers = parser.add_subparsers(dest='mode', required=True, help='Installation_mode')

    conda_parser = subparsers.add_parser('conda', help="Install for a local Conda environment.")
    conda_parser.add_argument(
        'task', nargs='?', default='all', choices=['all', 'tools', 'submodules'], help='Specific task to run (default: all)'
    )

    docker_parser = subparsers.add_parser('docker', help="Install by building a Docker image")

    return parser


def parse_args(argv=None):
    if argv is None:
        # If None is passed to this method, get the list of strings from the command line, excluding the script name.
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description="Installer for the nanopore analysis pipeline.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser = add_args(parser)

    args = parser.parse_args(argv)

    return args


def main(argv=None):
    args = parse_args(argv)

    if args.mode == 'docker':
        install_docker()
    elif args.mode == 'conda':
        install_conda(args.task)


if __name__ == "__main__":
    main()
