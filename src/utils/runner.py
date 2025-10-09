import pty
import subprocess
import logging
import yaml
import os
import signal
import sys
from pathlib import Path


def run_dorado(dorado_command: list, project_root, output_file: str = None):
    '''
    Temporarily add dorado executable to PATH and runs the specified command.
    This function, along with the run_wgbstools and run_uxm, are temporary and will be combined in future.
    '''

    # This is hard coded for now, but won't be in future.
    dorado_exe_path = project_root / "tools" / "dorado-0.9.6-linux-x64" / "bin" / "dorado"
    if not dorado_exe_path.exists():
        raise FileNotFoundError(f"dorado not found at {dorado_exe_path}")
    dorado_dir = str(dorado_exe_path.parent)
    env = os.environ.copy()
    env['PATH'] = f"{dorado_dir}{os.pathsep}{env['PATH']}"
    # Check whether the wgbstools is the first item in the command list, and if not prepend wgbstools to the list.
    if dorado_command[0] == 'dorado':
        full_command = dorado_command
    else:
        full_command = ['dorado'] + dorado_command

    print(f"Running dorado command: {' '.join(full_command)}")

    try:
        if output_file:
            with open(output_file, 'wb') as f_out:
                print(f"Writing to output file {output_file}")
                subprocess.run(full_command, check=True, env=env, stdout=f_out)
        else:
            subprocess.run(full_command, check=True, env=env)
    except Exception as e:
        print(f"CRITICAL: wgbstools command failed.", file=sys.stderr)
        raise e


def run_wgbstools(wgbstools_args: list, project_root):
    """
    Temporarily adds the wgbstools executable to the PATH, and runs the specified command.
    :param command:
    :return:
    """

    # THE RUNNER SHOULD BE ABLE TO FIND THE PROJECT_ROOT ON ITS OWN, THIS IS TEMPORARY
    wgbstools_exe_path = project_root / "externals" / "wgbs_tools" / "wgbstools"
    if not wgbstools_exe_path.exists():
        raise FileNotFoundError(f"wgbstools not found at {wgbstools_exe_path}")
    wgbs_dir = str(wgbstools_exe_path.parent)
    env = os.environ.copy()
    env['PATH'] = f"{wgbs_dir}{os.pathsep}{env['PATH']}"
    # Check whether the wgbstools is the first item in the command list, and if not prepend wgbstools to the list.
    if wgbstools_args[0] == 'wgbstools':
        full_command = wgbstools_args
    else:
        full_command = ['wgbstools'] + wgbstools_args

    print(f"Running wgbstools command: {' '.join(full_command)}")

    try:
        subprocess.run(full_command, check=True, env=env)
    except Exception as e:
        print(f"CRITICAL: wgbstools command failed.", file=sys.stderr)
        raise e

def run_uxm(uxm_args: list, project_root):
    """
        Temporarily adds the uxm executable to the PATH, and runs the specified command.
        :param command:
        :return:
        """

    # THE RUNNER SHOULD BE ABLE TO FIND THE PROJECT_ROOT ON ITS OWN, THIS IS TEMPORARY
    uxm_exe_path = project_root / "externals" / "UXM_deconv" / "uxm"
    if not uxm_exe_path.exists():
        raise FileNotFoundError(f"uxm not found at {uxm_exe_path}")
    uxm_dir = str(uxm_exe_path.parent)
    env = os.environ.copy()
    env['PATH'] = f"{uxm_dir}{os.pathsep}{env['PATH']}"

    # Check whether the wgbstools is the first item in the command list, and if not prepend wgbstools to the list.
    if uxm_args[0] == 'uxm':
        full_command = uxm_args
    else:
        full_command = ['uxm'] + uxm_args

    print(f"Running uxm command: {' '.join(full_command)}")

    try:
        subprocess.run(full_command, check=True, env=env)
    except Exception as e:
        print(f"CRITICAL: uxm command failed.", file=sys.stderr)
        raise e

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


def load_config(config_file="config.yaml"):
    """ Loads the pipeline config from a YAML file """
    logger = logging.getLogger('pipeline')
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Configuration loaded from {config_file}")
        return config
    except FileNotFoundError:
        logger.critical(f"Configuration file {config_file} not found.")
        raise


def run_external_command(command: list, cwd=None):
    """
    A simpler version of the runner.run_command method which worker src will use to execute
    external tools.
    Prints output to console and exists on failure.

    args:
        cwd: current working directory - runs the command from the specified directory.
    """
    print(f"Executing: {' '.join(command)}")

    try:
        subprocess.run(
            command,
            cwd=cwd,
            check=True
        )
    except FileNotFoundError:
        print(f"Error: Command '{command[0]}' not found. Make sure it's in your PATH", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}. Aborting.", file=sys.stderr)
        sys.exit(1)

def run_command(command: list, env=None):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    logger = logging.getLogger('pipeline')

    logger.info(f"Executing command: {' '.join(command)}")

    process = None

    try:
        # create pseudoterminal
        parent_fd, child_fd = pty.openpty()

        # start subprocess, connecting the output to the pseudoterminal.
        process = subprocess.Popen(
            command,
            stdout=child_fd,
            stderr=child_fd,
            text=True,
            preexec_fn=os.setpgrp,
            env=env
        )
        pgid = os.getpgid(process.pid)

        # We don't need the child part of the terminal anymore.
        os.close(child_fd)

        # Read the clean, live output from the parent part of the terminal
        with os.fdopen(parent_fd) as master_file:
            parent_fd = None # fd is now managed by master_file
            try:
                for line in iter(master_file.readline, ''):
                    clean_line = line.strip()
                    if not clean_line:
                        continue
                    if '\r' in line:
                        logger.debug(clean_line)
                    else:
                        logger.info(clean_line)
            except OSError:
                pass

        process.wait()

        if process.returncode != 0:
            logger.error(f"Command failed with exit code {process.returncode}")
            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            logger.info("Command completed successfully")

    except KeyboardInterrupt:
        logger.warning(">>> Keyboard interrupt detected. Terminating subprocess...")
        if process:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        raise
    except Exception as e:
        logger.critical(f"An unexpected error occurred: {e}")
        if process:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        raise e

    finally:
        if parent_fd is not None:
            os.close(parent_fd)
        if pgid is not None:
            try:
                logger.info(f"Sending SIGTERM ot process group {pgid}")
                os.killpg(pgid, signal.SIGTERM) # Sends a (polite) termination request to the entire process group.
                logger.info("Signal sent. Child process should now terminate.")
            except ProcessLookupError:
                logger.info("Child process already terminated.")
            except Exception as e:
                logger.error(f"Error during cleanup: {e}")






