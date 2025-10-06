import pty
import subprocess
import logging
import yaml
import os
import signal
import sys



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

def run_command(command: list):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    logger = logging.getLogger('pipeline')

    print(command)
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
            preexec_fn=os.setpgrp
        )

        # We don't need the child part of the terminal anymore.
        os.close(child_fd)

        # Read the clean, live output from the parent part of the terminal
        with os.fdopen(parent_fd) as master_file:
            try:
                for line in iter(master_file.readline, ''):
                    logger.info(line.strip())
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







