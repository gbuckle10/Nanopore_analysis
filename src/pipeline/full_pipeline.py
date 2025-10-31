from src.config.models import AppSettings
from src.pipeline import basecalling, alignment, deconvolution, methylation
import logging
from src.utils.logger import Logger

def full_pipeline_handler(args, config: AppSettings):
    logging.info("--- Running full pipeline from config ---")

    function_map = {
        'basecalling': basecalling.full_basecalling_handler,
        'align': alignment.full_alignment_handler,
        'methylation': methylation.pileup_handler,
        'deconvolution': deconvolution.deconvolution_handler
    }

    io_map = {
        'basecalling': {
            'input_file': 'config.pipeline_steps.setup.paths.full_pod5_path',
            'output_file': 'config.pipeline_steps.basecalling.paths.full_demultiplexed_output_dir'
        },
        'alignment': {
            'input_file': 'config.pipeline_steps.basecalling.paths.full_unaligned_bam_path',
            'output_file': 'config.pipeline_steps.alignment.paths.full_aligned_bam_path'
        },
        'methylation': {
            'input_file': 'config.pipeline_steps.alignment.paths.full_aligned_bam_path',
            'output_file': 'config.pipeline_steps.methylation.paths.final_bed_file'
        },
        'deconvolution': {
            'input_file': 'config.pipeline_steps.analysis.paths.full_deconv_input_path',
            'output_file': 'config.pipeline_steps.analysis.paths.full_deconv_results_path'
        }
    }

    steps_to_run = config.pipeline_control.run_steps
    active_steps = [step_name for step_name, should_run in steps_to_run if should_run]

    if not active_steps:
        logging.warning("WARNING: No active steps found in the final configuration. Nothing to do.")

    logging.info("STEPS TO RUN: ")
    for i, step_name in enumerate(active_steps, 1):
        logging.info(f"    Step {i}: {step_name}")
    logging.info("----------------------------------")

    for step_name in active_steps:
        logging.info(f">>> EXECUTING STEP: {step_name}")
        step_func = function_map.get(step_name)
        if not step_func:
            logging.warning(f"WARNING: No function found for step '{step_name}'. Skipping")
            continue

        # Create a deep copy of the main 'args' object


        step_func(args, config)
def setup_parsers(subparsers, parent_parser, config):

    run_parser=subparsers.add_parser(
        'run',
        help="Run the full pipeline using steps defined in the config file.",
        parents=[parent_parser]
    )

    print("Setting up run parser")

    # Register commands from modules.
    # Add the arguments from each individual step to the run_parser
    basecalling.add_all_arguments_to_parser(run_parser, config)
    alignment.add_all_arguments_to_parser(run_parser, config)
    deconvolution.add_all_arguments_to_parser(run_parser, config)
    run_parser.set_defaults(func=full_pipeline_handler)
