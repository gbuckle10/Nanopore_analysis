from src import PROJECT_ROOT
from src.config.models import AppSettings


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
    config.pipeline_steps.setup.paths.build_paths(common_paths)
    config.pipeline_steps.basecalling.paths.build_paths(common_paths)
    config.pipeline_steps.alignment.paths.build_paths(common_paths)
    config.pipeline_steps.methylation.paths.build_paths(common_paths)
    config.pipeline_steps.analysis.paths.build_paths(common_paths)




