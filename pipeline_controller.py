import subprocess
import os
import sys
import yaml

# --- Configuration ---
def load_config(config_file="config.yaml"):
    """ Loads the pipeline config from a YAML file """
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def run_and_capture(config, command):
    """ Runs a command inside the conda environment and stores the output as a variable"""
    conda_env = config['conda_env_name']

    full_command = ["conda", "run", "-n", conda_env] + command

    print(f"--- Running : {' '.join(full_command)} --- ")

    try:
        result = subprocess.run(full_command, check=True, capture_output=True, text=True)
        if result.stderr:
            print(result.stderr.strip())
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"--- ERROR IN COMMAND ---")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise

    print(result)

def run_and_stream(config, command):
    """ Runs a command inside the conda environment and streams the output """
    conda_env = config['conda_env_name']

    full_command = ["conda", "run", "-n", conda_env] + command

    print(f"--- Running : {' '.join(full_command)} --- ")

    try:
        result = subprocess.run(full_command, check=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"--- ERROR IN COMMAND ---")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise

    print(result)

def run_setup(config):
    """ Executes the 00_setup.sh script """
    print(">>> Starting step 0: Setup")
    script_path = "scripts/00_setup.sh"

    tasks = config.get('run_setup_tasks', {})

    download_fast5 = str(tasks.get('download_fast5_data', False)).lower()
    convert_fast5_to_pod5 = str(tasks.get('convert_fast5_to_pod5', False)).lower()

    command = [
        "bash", script_path,
        config['dorado_version'],
        config['fast5_download_url'],
        config['fast5_input_dir'],
        config['num_fast5_files'],
        config['pod5_dir'],
        download_fast5,
        convert_fast5_to_pod5
    ]

    dorado_path = run_and_capture(config, command)

    #if not dorado_path or not os.path.exists(dorado_path):
    #    raise FileNotFoundError(f"Setup script failed to return a valid path.")

    #print(f"--- Successfully found Dorado executable at: {dorado_path} ---\n ")

    #return dorado_path

def run_basecalling(config):
    """ Executes the basecalling script using the 01_basecalling.sh script """
    print(">>> Starting step 1: Basecalling")
    script_path = "scripts/01_basecalling.sh"

    command = [
        "bash", script_path,
        config['dorado_model_name'],
        config['model_speed'],
        config['basecalling_modifications'],
        config['basecalling_batch_size']
    ]

    run_and_stream(config, command)
def main():
    """ Main entry point for the pipeline controller """
    config = load_config()
    steps_to_run = config.get('run_steps', [])

    if not steps_to_run:
        print("Error, there aren't any steps for me to run. Check config.yaml.")
        sys.exit(1)

    print(f"--- Pipeline will execute the following steps: {steps_to_run} ---")

    if 'setup' in steps_to_run:
        run_setup(config)
    if 'basecalling' in steps_to_run:
        run_basecalling(config)

if __name__ == '__main__':
    main()



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
