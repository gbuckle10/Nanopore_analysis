from functools import reduce

from src.config.models import AppSettings
from src.pipeline import basecalling, alignment, deconvolution, methylation
import logging
from src.utils.logger import Logger
import copy


def get_nested_attr(obj, attr_string: str):
    """
    Gets a nested attribute from an object using a dot-separated string.
    For example, get_nested_attr(config, 'globals.threads') is the same as config.globals.threads. We can't just use dot
    notation if the attribute we're trying to get is stored in a string variable.
    :param obj:
    :param attr_string:
    :return:
    """

    # Split the string into a list of attribute names - 'a.b.c' becomes ['a', 'b', 'c']
    attributes = attr_string.split('.')

    # Use 'reduce' to apply getattr cumulatively.
    return reduce(getattr, attributes, obj)


def full_pipeline_handler(args, config: AppSettings):
    logging.info("--- Running full pipeline from config ---")

    function_map = {
        'basecalling': basecalling.full_basecalling_handler,
        'align': alignment.full_alignment_handler,
        'methylation': methylation.pileup_handler,
        'analysis': deconvolution.deconvolution_handler
    }

    # Make an I/O map giving the default input and output of each step.
    io_map = {
        'basecalling': {
            'input_file': 'pipeline_steps.setup.paths.full_pod5_path',
            'output_dir': 'pipeline_steps.basecalling.paths.full_demultiplexed_output_dir'
        },
        'align': {
            'input_file': 'pipeline_steps.basecalling.paths.full_unaligned_bam_path',
            'output_dir': 'pipeline_steps.align.paths.full_aligned_bam_path'
        },
        'methylation': {
            'input_file': 'pipeline_steps.align.paths.full_aligned_bam_path',
            'output_dir': 'pipeline_steps.methylation.paths.final_bed_file'
        },
        'deconvolution': {
            'input_file': 'pipeline_steps.analysis.paths.full_deconv_input_path',
            'output_dir': 'pipeline_steps.analysis.paths.full_deconv_results_path'
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

        # Create a deep copy of the main 'args' object to avoid messing up the variables for later steps.
        step_args = copy.deepcopy(args)
        step_io_args = io_map.get(step_name, {})
        # Add the io args from the io args map.
        for key, config_path_str in step_io_args.items():
            print(f"Setting attribute {key} using {config_path_str}")
            setattr(
                step_args,
                key,
                get_nested_attr(config, config_path_str)
            )
        print(f"After adding the io args, step args is:\n{step_args}")

        step_func(step_args, config)


def setup_parsers(subparsers, parent_parser, config):
    run_parser = subparsers.add_parser(
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
