import signal
import subprocess
import logging
import yaml
import shlex
import os

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


def run_command(command: list, config: dict, use_conda: bool=True):
    '''
    Runs the commands for each step, logs the outputs in real time and handles errors.
    '''

    logger = logging.getLogger('pipeline')

    command_str = shlex.join(command)

    full_command = [
        "script",
        "-q",
        "--flush",
        "-c", command_str,
        "/dev/null"
    ]

    logger.info(f"Executing command: {' '.join(full_command)}")

    process = None

    try:
        process = subprocess.Popen(
            full_command,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
            preexec_fn=os.setsid
        )
        if process.stdout:
            for line in iter(process.stdout.readline, ''):
                logger.info(line.strip())

        process.wait()

        if process.returncode != 0:
            logger.error(f"{process.returncode}")
            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            logger.info("Command completed successfully")

    except KeyboardInterrupt:
        logger.warning(">>> Keyboard interrupt")
        if process:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        raise
    except subprocess.CalledProcessError as e:
        logger.critical(' '.join(e.cmd))
        raise e
    except FileNotFoundError:
        logger.critical(f"{full_command[0]}")
        raise
    except Exception as e:
        logger.critical(e)
        if process:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        raise e







