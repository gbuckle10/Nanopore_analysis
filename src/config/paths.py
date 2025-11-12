import argparse
import logging
from pathlib import Path
import inspect
from src import PROJECT_ROOT
from .models import AppSettings


def build_config_paths(config: AppSettings) -> None:
    """
    Builds all absolute Path objects from the experiment_root and filenames.
    This function modifies the config object in-place.
    """

    # Make the central common paths.
    root = PROJECT_ROOT  # One day we'll be able to change the root, let the user select it.
    common_paths = config.paths
    common_paths.root = root
    common_paths.data_dir = root / "data"
    common_paths.log_dir = root / "logs"
    common_paths.results_dir = root / "results"
    common_paths.reference_genome_dir = root / "reference_genomes"
    common_paths.externals_dir = root / "externals"

    # Run build_paths method for each specific step
    config.pipeline_steps.setup.paths.build_and_validate(common_paths)
    config.pipeline_steps.basecalling.paths.build_and_validate(common_paths)

    if config.pipeline_steps.basecalling.params.demultiplex:
        unaligned_bam_source = config.pipeline_steps.basecalling.paths.full_demultiplexed_output_dir
    else:
        unaligned_bam_source = config.pipeline_steps.basecalling.paths.full_unaligned_bam_path
    config.pipeline_steps.align.paths.build_and_validate(common_paths, unaligned_bam_source)
    config.pipeline_steps.methylation.paths.build_and_validate(common_paths)
    config.pipeline_steps.analysis.paths.build_and_validate(common_paths)


def update_config_from_args(config: AppSettings, args: argparse.Namespace, parser: argparse.ArgumentParser):
    """
    Updates the Pydantic config object in-place with values from argparse
    """
    args_dict = vars(args)

    for dest, value in args_dict.items():
        default_value = parser.get_default(dest)
        if '.' in dest and value is not None and value != default_value:
            # print(f"Setting {dest} to {value}, because it's different to {default_value}")
            _set_config_attribute(config, dest, value)


def _set_config_attribute(obj, path, value):
    keys = path.split('.')
    current_obj = obj
    for key in keys[:-1]:
        current_obj = getattr(current_obj, key)

    final_key = keys[-1]

    setattr(current_obj, final_key, value)


def validate_active_steps(config: AppSettings):
    """
    Validates that if a step is active, all of its required inputs are present.
    :return:
    """
    steps = config.pipeline_steps

    first_active_step = True

    for step_name, should_run in config.pipeline_control.run_steps:
        if not should_run:
            continue  # Skip validation for inactive steps

        step_model = getattr(steps, step_name, None)
        if not step_model or not hasattr(step_model, 'paths'):
            logging.warning(f"No validation path model found for active step: {step_name}")
            continue
        validate_method = getattr(step_model.paths, '_validate', None)

        if callable(validate_method):
            try:
                is_standalone_check = first_active_step
                sig = inspect.signature(validate_method)

                if 'is_standaline_run' in sig.parameters:
                    validate_method(is_standalone_check)
                else:
                    validate_method()
                first_active_step = False
            except (ValueError, FileNotFoundError) as e:
                raise ValueError(f"Validation failed for step '{step_name}': {e}") from e
