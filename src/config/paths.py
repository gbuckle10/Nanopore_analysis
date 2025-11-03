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
    config.pipeline_steps.align.paths.build_and_validate(common_paths)
    config.pipeline_steps.methylation.paths.build_and_validate(common_paths)
    config.pipeline_steps.analysis.paths.build_and_validate(common_paths)


def validate_active_steps(config: AppSettings):
    """
    Validates that if a step is active, all of its required inputs are present.
    :return:
    """
    steps = config.pipeline_steps

    for step_name, should_run in config.pipeline_control.run_steps:

        if not should_run:
            continue  # Skip validation for inactive steps
        step_model = getattr(steps, step_name, None)

        if not step_model:
            continue
        validate_method = getattr(step_model.paths, '_validate', None)

        if callable(validate_method):
            try:
                validate_method()
            except (ValueError, FileNotFoundError) as e:
                raise ValueError(f"Validation failed for step '{step_name}': {e}") from e