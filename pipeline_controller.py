import subprocess
import sys
import yaml
from scripts.analysis_logic import *


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


def run_and_log(config, command, log_path="logs/wowow.txt"):
    """ Runs a command inside the conda environment and streams the output """
    conda_env = config['conda_env_name']
    full_command = ["conda", "run", "-n", conda_env] + command

    print(f"--- Running : {' '.join(full_command)} --- ")

    with open(log_path, 'w') as log_file:
        process = subprocess.Popen(
            full_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )

        for line in process.stdout:
            log_file.write(line)
            print(line, end='')
            sys.stdout.flush()

        process.wait()
        if process.returncode != 0:
            print(f"\n--- ERROR in command. Check log for details: {log_path} ---")
            raise subprocess.CalledProcessError(process.returncode, process.args)
        else:
            print("--- Command successful ---\n")

    '''
    try:
        subprocess.run(full_command, check=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"--- ERROR IN COMMAND ---")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise
    '''


def run_deconvolution_submodule(config):
    print(">>> Starting the deconvolution process using the meth_atlas submodule")

    atlas_file = f"{config['paths']['atlas_dir']}{config['paths']['atlas_file_illumina']}"
    file_to_deconvolve = f"{config['paths']['analysis_dir']}{config['paths']['file_for_deconvolution']}"
    output_file = f"{config['paths']['analysis_dir']}{config['paths']['deconvolution_results']}"
    deconvolve_script = "externals/meth_atlas/deconvolve.py"
    # deconvolve_script = "externals/meth_atlas/deconvolve_genome_coordinates.py"

    command = [
        "python",
        deconvolve_script,
        "-a",
        atlas_file,
        file_to_deconvolve,
        "--out_dir", config['paths']['analysis_dir']
    ]

    print(f"--- Running : {' '.join(command)} --- ")

    try:
        result = subprocess.run(command, check=True, capture_output=False, text=True)
        if result.stderr:
            print(result.stderr.strip())
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"--- ERROR IN COMMAND ---")
        print(f"Exit code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise
    except Exception as error:
        print(error)


def run_setup(config):
    """ Executes the 00_setup.sh script """
    print(">>> Starting step 0: Setup")
    script_path = "scripts/00_setup.sh"
    config_file = "config.yaml"

    command = ["bash", script_path, config_file]

    run_and_log(config, command)

    # if not dorado_path or not os.path.exists(dorado_path):
    #    raise FileNotFoundError(f"Setup script failed to return a valid path.")

    # print(f"--- Successfully found Dorado executable at: {dorado_path} ---\n ")

    # return dorado_path


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

    run_and_log(config, command)


def run_alignment(config):
    """ Executes the alignment script: 02_alignment.sh """
    print(">>> Starting step 2: Alignment and Indexing")
    script_path = "scripts/02_alignment.sh"

    command = [
        "bash", script_path,
        config['reference_genome_dir'],
        config['reference_genome_url'],
        config['reference_genome_name'],
        config['indexed_ref_gen_name'],
        config['basecalled_output_dir'],
        config['unaligned_bam_name'],
        config['alignment_output_dir'],
        config['aligned_bam_name'],
        config['threads'],
        config['sort_memory_limit']
    ]

    run_and_log(config, command)


def run_alignment_qc(config):
    print(">>> Running QC on aligned and indexed data")
    script_path = "scripts/03_alignment_qc.sh"

    command = [
        "bash", script_path,
        config['alignment_output_dir'],
        config['aligned_bam_name'],
        config['qc_dir'],
        config['alignment_flagstat_name'],
        config['alignment_stats_name']
    ]

    run_and_log(config, command)


def run_methylation_summary(config):
    print(">>> Running methylation calling")
    script_path = "scripts/04_methylation_summary.sh"

    command = [
        "bash", script_path,
        config['aligned_bam_name'],
        config['alignment_output_dir'],
        config['methylation_bed_name'],
        config['reference_fasta'],
        config['threads'],
        config['methylation_dir'],
        config['methylation_log_file'],
        config['reference_fasta']
    ]

    run_and_log(config, command)


def run_analysis(config):
    """
    Executes the full deconvolution analysis in two steps. The first step is to
    arrange the data in the bed file in such a way that it can be compared to the
    meth_atlas.
    """

    print(">>> Starting analysis workflow")

    try:

        # One day change this so that the yaml structure mirrors the file structure. config['paths']['methylation_dir']['methylation_bed_name']
        bed_file_path = config['paths']['methylation_dir'] + config['paths']['methylation_bed_name']
        manifest_file_path = config['paths']['atlas_dir'] + config['paths']['illumina_manifest']
        file_for_decon_path = config['paths']['analysis_dir'] + config['paths']['file_for_deconvolution']
        illumina_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['atlas_file_illumina']
        geco_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['atlas_file_genome_coordinate']
        uxm_atlas_file_path = config['paths']['atlas_dir'] + config['paths']['uxm_atlas_name']
        chunk_size = int(config['parameters']['analysis']['methylation_aggregation_chunksize'])

        generate_deconvolution_files(
            bed_file=bed_file_path,
            manifest_file=manifest_file_path,
            output_file=file_for_decon_path,
            range_atlas_file=uxm_atlas_file_path,
            chunk_size=chunk_size
        )


        '''
        convert_atlas_to_genome_coordinates(
            output_file=geco_atlas_file_path,
            atlas_file=illumina_atlas_file_path,
            manifest_file=manifest_file_path
        )
        '''

        #format_atlas_file(atlas_file=uxm_atlas_file_path)

        #generate_aggregated_deconvolution_file(file_for_decon_path, uxm_atlas_file_path)

        #run_deconvolution_submodule(config)

    except Exception as e:
        print(f"--- ERROR during deconvolution prep: {e} ---")
        raise


def main():
    """ Main entry point for the pipeline controller """
    config = load_config()

    try:
        steps_to_run = config['pipeline_control']['run_steps']
    except (FileNotFoundError, KeyError) as e:
        print(f"FATAL ERROR: Could not load required configuration.")
        print(f"Details: {e}")
        sys.exit(1)

    if not steps_to_run:
        print("Error, there aren't any steps for me to run. Check config.yaml.")
        sys.exit(1)

    print(f"--- Pipeline will execute the following steps: {steps_to_run} ---")

    if 'setup' in steps_to_run:
        run_setup(config)
    if 'basecalling' in steps_to_run:
        run_basecalling(config)
    if 'alignment' in steps_to_run:
        run_alignment(config)
    if 'methylation_summary' in steps_to_run:
        run_methylation_summary(config)
    if 'analysis' in steps_to_run:
        run_analysis(config)


if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
