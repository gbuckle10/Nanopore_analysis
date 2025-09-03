import subprocess
import argparse
import sys
import yaml

# --- Configuration ---
def load_config(config_file="config.yaml"):
    """ Loads the pipeline config from a YAML file """
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def run_in_conda(config, command):
    """ Runs a command inside the conda environment """
    conda_env = config['conda_env_name']

    full_command = ["conda", "run", "-n", conda_env] + command

    print(f"--- Running : {' '.join(full_command)} --- ")

    try:
        result = subprocess.run(full_command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"--- ERROR IN COMMAND ---")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise

def run_setup(config):
    """ Executes the 00_setup.sh script """
    print(">>> Starting step 0: Setup")
    script_path = "scripts/00_setup.sh"
    command = [
        "bash", script_path,
        config['reference_genome_url'],
        config['reference_genome_dir'],
        config['reference_fasta'],
        config['reference_index'],
        config['dorado_version']
    ]

    run_in_conda(config, command)

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


if __name__ == '__main__':
    main()



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
