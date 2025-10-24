import argparse
import subprocess
import sys
from datetime import datetime
from src.utils.process_utils import run_command
from src.utils.config_utils import get_project_root, load_config, deep_merge
from deconvolution import Deconvolution


class PipelineController:
    def __init__(self, args, config_name="config.yaml"):
        """Initialise pipeline, set up logging and load config"""
        run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file_path = f"logs/{run_timestamp}_pipeline_run.log"
        self.args = args
        self.logger = setup_logger(log_file=log_file_path)
        self.logger.info(f"PIPELINE RUN STARTED AT: {run_timestamp}")
        self.logger.info(f"The log for this run will be saved to: {log_file_path}")
        self.logger.info("=" * 80)
        self.tool_paths = {}
        self.project_root = get_project_root()
        self.scripts_path = self.project_root / "src" / "pipeline"
        config_path = self.project_root / config_name
        self.user_config = load_config(str(config_path))
        self.runtime_config = load_config("runtime_config.yaml")
        self.config = deep_merge(self.user_config, self.runtime_config)

        self.pipeline_steps = {
            'basecalling': {
                'step_number': 1,
                'script_name': '01_basecalling.py',
                'description': 'Basecalling'
            },
            'alignment': {
                'step_number': 2,
                'script_name': '02_alignment.py',
                'description': 'Alignment'
            },
            'align_qc': {
                'step_number': 3,
                'script_name': '03_alignment_qc.py',
                'description': 'Alignment QC'
            },
            'methylation_summary': {
                'step_number': 4,
                'script_name': '04_methylation_summary.py',
                'description': 'Methylation Summary'
            }
        }

        self.command_map = {
            'basecalling': self.run_basecalling,
            'align': self.run_alignment,
            'methylation_summary': self.run_methylation_summary,
            'deconvolution': self.run_deconvolution,
            'all': self.run_active_steps
        }

        self.alias_map = {
            'basecall': 'basecalling',
            'alignment': 'align',
            'deconv': 'deconvolution'
        }

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

        self.logger.info('--- pipeline Started ---')
        self.logger.info(f"--- pipeline will execute the following steps: {steps_to_run} ---")

        if 'basecalling' in steps_to_run:
            self._run_step('basecalling', script_args)
        if 'align' in steps_to_run:
            self._run_step('align', script_args)
        if 'align_qc' in steps_to_run:
            self._run_step('align_qc', script_args)
        if 'methylation_summary' in steps_to_run:
            self._run_step('methylation_summary', script_args)
        if 'deconvolution' in steps_to_run:
            self._run_step('deconvolution', script_args)

    def _run_step(self, step_name, script_args):
        """
        A generic method to run a pipeline step
        """
        step_info = self.pipeline_steps.get(step_name)
        if not step_info:
            self.logger.error(f"Undefined pipeline step: {step_name}")
            sys.exit(1)

        step_num = step_info['step_number']
        description = step_info['description']
        script_path = self.scripts_path / step_info['script_name']

        self.logger.info("=" * 80)
        self.logger.info(f">>> Starting step {step_num}: {description}")

        command = [
            sys.executable,
            str(script_path)
        ] + script_args

        run_command(command)

    def run_basecalling(self, script_args):
        self._run_step('basecalling', script_args)

    def run_alignment(self, script_args):
        """ Executes the alignment script: 02_alignment.sh """
        self._run_step('alignment', script_args)

    def run_alignment_qc(self, script_args):
        self._run_step('alignment_qc', script_args)

    def run_methylation_summary(self, script_args):
        self._run_step('methylation_summary', script_args)

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


    def main(self, argv=None):
        try:
            parser = argparse.ArgumentParser(description='Internal pipeline steps.')
            parser.add_argument('command', help='The pipeline step to run')

            args, remaining_argv = parser.parse_known_args(argv)
            user_command = args.command if args.command is not None else 'all'

            command_map = {
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
        except KeyboardInterrupt:
            print("Process terminated by user.", file=sys.stderr)
            sys.exit(130)

if __name__ == '__main__':
    controller = PipelineController()
    controller.main()
