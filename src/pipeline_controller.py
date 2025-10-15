import argparse
import subprocess
import sys
import yaml
from datetime import datetime
from utils.runner import run_command, run_wgbstools, run_uxm, get_project_root
from utils.logger import setup_logger
from utils.file_conversion import apply_runtime_config, ensure_tool_symlink
from deconvolution import Deconvolution


class PipelineController:
    def __init__(self, config_name="config.yaml"):
        """Initialise pipeline, set up logging and load config"""
        run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file_path = f"logs/{run_timestamp}_pipeline_run.log"
        self.logger = setup_logger(log_file=log_file_path)
        self.logger.info(f"PIPELINE RUN STARTED AT: {run_timestamp}")
        self.logger.info(f"The log for this run will be saved to: {log_file_path}")
        self.logger.info("=" * 80)
        self.tool_paths = {}
        self.project_root = get_project_root()

        config_path = self.project_root / config_name
        self.config = self.load_config(str(config_path))

        try:
            self.steps_to_run = self.config['pipeline_control']['run_steps']
        except (FileNotFoundError, KeyError) as e:
            self.logger.error(f"FATAL ERROR: Could not load required configuration.")
            self.logger.error(f"Details: {e}")
            sys.exit(1)
        if not self.steps_to_run:
            self.logger.error("Error, there aren't any steps for me to run. Check config.yaml.")
            sys.exit(1)
        self.logger.info(" ---------------- Starting main run ----------------")

    def load_config(self, config_file="config.yaml"):
        """ Loads the pipeline config from a YAML file """

        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            self.logger.info(f"Configuration loaded from {config_file}")
            return config
        except FileNotFoundError:
            self.logger.critical(f"Configuration file {config_file} not found.")
            raise

    def run_active_steps(self, script_args):
        """
        Executes the scripts specified in the config.yaml file.
        """
        self.logger.info("=" * 80)
        self.logger.info("------------ Starting whole pipline ------------")

        try:
            steps_to_run = self.config['pipeline_control']['run_steps']
        except (FileNotFoundError, KeyError) as e:
            self.logger.error(f"FATAL ERROR: Could not load required configuration.")
            self.logger.error(f"Details: {e}")
            sys.exit(1)

        if not steps_to_run:
            self.logger.error("Error, there aren't any steps for me to run. Check config.yaml.")
            sys.exit(1)

        self.logger.info('--- Pipeline Started ---')
        self.logger.info(f"--- Pipeline will execute the following steps: {steps_to_run} ---")

        if 'setup' in steps_to_run:
            self.run_setup()
        if 'basecalling' in steps_to_run:
            self.run_basecalling()
        if 'align' in steps_to_run:
            self.run_alignment()
        if 'align_qc' in steps_to_run:
            self.run_alignment_qc()
        if 'methylation_summary' in steps_to_run:
            self.run_methylation_summary()
        if 'deconvolution' in steps_to_run:
            self.run_deconvolution()

    def run_setup(self, script_args):
        """ Executes the 00_setup.sh script """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting Step 0: Setup")

        script_path = self.project_root / "src" / "00_setup.py"
        command = [
            sys.executable,
            str(script_path)
        ]

        run_command(command)

        self.tool_paths = apply_runtime_config()

        wgbstool_path = self.project_root / "externals" / "wgbs_tools"
        wgbstools_sl_path = self.project_root / "externals" / "wgbs_tools" / "wgbstools"
        wgbstools_py_path = self.project_root / "externals" / "wgbs_tools" / "src" / "python" / "wgbs_tools.py"
        ensure_tool_symlink(wgbstools_sl_path, wgbstools_py_path)

        self.logger.info("Symlink for wgbstools good.")

        run_wgbstools(['wgbstools', '--version'], self.project_root)

        uxm_sl_path = self.project_root / "externals" / "UXM_deconv" / "uxm"
        uxm_py_path = self.project_root / "externals" / "UXM_deconv" / "src" / "uxm.py"
        ensure_tool_symlink(uxm_sl_path, uxm_py_path)

        self.logger.info("Symlink for uxm good.")

        run_uxm(['uxm', '--help'], self.project_root)

    def run_basecalling(self, script_args):
        """ Executes the basecalling script using the 01_basecalling.sh script """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting step 1: Basecalling")

        script_path = self.project_root / "src" / "01_basecalling.py"
        command = [
                      sys.executable,
                      str(script_path)
                  ] + script_args

        run_command(command)

    def run_alignment(self, script_args):
        """ Executes the alignment script: 02_alignment.sh """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting step 2: Alignment and Indexing")

        script_path = self.project_root / "src" / "02_alignment.py"
        command = [
            sys.executable,
            str(script_path)
        ]

        run_command(command)

    def run_alignment_qc(self, script_args):
        self.logger.info("=" * 80)
        self.logger.info(">>> Running QC on aligned and indexed data")

        script_path = "03_alignment_qc.py"
        command = [
            sys.executable,
            str(script_path)
        ]

        run_command(command)

    def run_methylation_summary(self, script_args):
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting step 4: Methylation Summary")

        script_path = self.project_root / "src" / "04_methylation_summary.py"
        command = [
            sys.executable,
            str(script_path)
        ]

        run_command(command)

    def run_deconvolution(self, script_args):
        """
        Executes the full deconvolution analysis in two steps. The first step is to
        arrange the data in the bed file in such a way that it can be compared to the
        meth_atlas.
        """
        self.logger.info("=" * 80)
        self.logger.info(">>> Starting preparation for deconvolution.")

        atlas_file = self.project_root / "data" / "atlas" / "full_atlas_geco.csv"
        file_to_deconv = self.project_root / "data" / "processed" / "deconvolution_geco.csv"
        output_file = self.project_root / "data" / "processed"

        deconv_handler = Deconvolution(self.config, self.tool_paths, atlas_file, file_to_deconv, output_file)

        if "deconvolution_prep" in self.steps_to_run:
            deconv_handler.prepare()
        deconv_handler.run()

    def ___run_deconvolution_submodule(self, script_args):
        self.logger.info(">>> Starting the deconvolution process using the meth_atlas submodule")

        atlas_file = f"data/atlas/{self.config['paths']['atlas_file_uxm']}"
        file_to_deconvolve = f"data/processed/{self.config['paths']['file_for_deconvolution_uxm']}"

        self.logger.info(f"Deconvolving file {file_to_deconvolve} using atlas file {atlas_file}")

        deconvolve_script = "externals/meth_atlas/deconvolve_genome_coordinates.py"
        # deconvolve_script = "externals/meth_atlas/deconvolve.py"

        command = [
            "python",
            deconvolve_script,
            "-a",
            atlas_file,
            file_to_deconvolve,
            "--out_dir", self.config['paths']['analysis_dir']
        ]

        self.logger.info(f"--- Running : {' '.join(command)} --- ")

        try:
            with subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    universal_newlines=True
            ) as process:
                if process.stdout:
                    for line in iter(process.stdout.readline, ''):
                        self.logger.info(f"[Subprocess] {line.strip()}")

            if process.returncode != 0:
                self.logger.error(f"Deconvolution script failed with exit code {process.returncode}.")
                raise subprocess.CalledProcessError(process.returncode, process.args)

            self.logger.info("Deconvolution script completed successfully.")
        except FileNotFoundError:
            self.logger.error(f"Command not found: {command[0]}. Make sure python is in PATH")
            raise
        except subprocess.CalledProcessError as e:
            self.logger.critical("Halting process due to deconvolution script failure.")
            raise e

    def main(self, argv=None):
        parser = argparse.ArgumentParser(description='Internal pipeline steps.')
        parser.add_argument('command', help='The pipeline step to run')

        args, remaining_argv = parser.parse_known_args(argv)
        user_command = args.command if args.command is not None else 'all'

        command_map = {
            'setup': self.run_setup,
            'basecalling': self.run_basecalling,
            'align': self.run_alignment,
            'methylation_summary': self.run_methylation_summary,
            'deconvolution': self.run_deconvolution,
            'all': self.run_active_steps
        }

        alias_map = {
            'basecall': 'basecalling',
            'alignment': 'align',
            'deconv': 'deconvolution'
        }
        function_to_run = command_map.get(alias_map.get(user_command, user_command))

        if function_to_run:
            function_to_run(remaining_argv)
        else:
            print(f"Error: Unknown command '{args.command}'")
            parser.print_help()


if __name__ == '__main__':
    controller = PipelineController()
    controller.main()
