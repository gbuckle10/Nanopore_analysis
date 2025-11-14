import argparse
import logging
import os
import sys

from src.config.models import AppSettings
from src.utils.cli_utils import add_io_arguments
from pathlib import Path

from src.utils.file_utils import ensure_dir_exists
from src.utils.process_utils import LiveDisplayHandler, SilentHandler, TunableHandler, LineBufferingHandler
from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def deconvolution_handler(config: AppSettings):
    logger.info("Running deconvolution step.")
    analysis_settings = config.pipeline_steps.analysis
    algorithm = analysis_settings.params.deconv_algorithm



    if algorithm == "uxm":
        wgbstools_runner = ToolRunner(analysis_settings.tools.wgbstools_exe, '-o', handler_class=SilentHandler)
        uxm_runner = ToolRunner(analysis_settings.tools.uxm_exe, '-o', handler_class=LiveDisplayHandler)
        _run_uxm_algorithm(config, wgbstools_runner, uxm_runner)
    elif algorithm == "nnls":
        nnls_runner = ToolRunner(analysis_settings.tools.methatlas_exe, '--out_dir', handler_class=LiveDisplayHandler)
        _run_nnls_algorithm()
    else:
        logger.error("You have chosen an algorithm that doesn't exist. Unfortunately I can't do this.")


# ===============================
# UXM Algorithm
# ===============================

def _run_uxm_algorithm(config: AppSettings, wgbstools_runner: ToolRunner, uxm_runner: ToolRunner):
    """
    Orchestrates the entire UXM deconvolution pipeline
    """

    logger.info("Running UXM Deconvolution...")

    deconv_tasks = _generate_uxm_deconvolution_tasks(config)

    if not deconv_tasks:
        logger.warning("No valid UXM deconvolution jobs found. Exiting.")
        return

    logger.info(f"There are {len(deconv_tasks)} deconvolution task(s) to do.")

    atlas_path = config.pipeline_steps.analysis.paths.full_atlas_path
    failed_samples = []
    for i, task in enumerate(deconv_tasks):
        logger.info(f"--- Running deconvolution task {i + 1}/{len(deconv_tasks)}")
        try:
            _run_single_uxm_deconv(task, atlas_path, wgbstools_runner, uxm_runner)
        except Exception as e:
            logger.error(f"Failed to process sample {task['base_name']}. Error: {e}")
            failed_samples.append(task['base_name'])

    if failed_samples:
        logger.error(f"UXM pipeline finished with {len(failed_samples)} failed sample(s): {', '.join(failed_samples)}")
    else:
        logger.info("All UXM deconvolution jobs finished successfully.")


def _generate_uxm_deconvolution_tasks(config: AppSettings) -> list[dict]:
    """
    Generates a list of deconvolution jobs by getting all the paths from the config object and ensuring a consistent
    architecture for the output.
    """
    analysis_paths = config.pipeline_steps.analysis.paths
    input_path = analysis_paths.full_deconv_input_path
    logger.info(f"Deconvoluting the files in {input_path} using the uxm algorithm")

    pat_dir = analysis_paths.full_pat_dir
    filtered_pat_dir = analysis_paths.full_pat_filtered_dir
    deconvolution_dir = analysis_paths.full_deconv_output_path

    logger.debug(f"Pat files from bam2pat will be saved in {pat_dir}")
    logger.debug(f"Filtered pat files will be saved in {filtered_pat_dir}")
    logger.debug(f"Deconvoluted files will be saved in {deconvolution_dir}")

    tasks = []
    if not input_path.exists():
        logger.error(f"Input path does not exist: {input_path}")
        return []

    bam_files = []
    # The input is a directory
    if input_path.is_dir():
        logger.info(f"The input is a directory. I am searching for .bam files in {input_path}")
        bam_files = sorted(list(input_path.glob('*.bam')))
    elif input_path.is_file() and input_path.suffix == '.bam':
        logger.info(f"The input is a single bam file: {input_path}")
        bam_files = [input_path]

    if not bam_files:
        logger.warning(f"No .bam files found for deconvolution at: {input_path}")
        return []

    for bam_file in bam_files:
        base_name = bam_file.stem # Maybe we need to find a way to clean up the file names if there are multiple extensions
        tasks.append({
            'input_bam': bam_file,
            'base_name': base_name,
            'pat_dir': pat_dir,
            'pat_file_path': pat_dir / f"{base_name}.pat.gz",  # We also need to know what file to look for
            'filtered_pat_path': filtered_pat_dir / f"{base_name}.pat",
            'filtered_index_pat_path': filtered_pat_dir / f"{base_name}.pat.gz",
            'deconvoluted_file': deconvolution_dir / f"{base_name}_uxm_results.csv"
        })

    return tasks

def _run_single_uxm_deconv(task: dict, atlas_path: Path, wgbstools_runner: ToolRunner, uxm_runner: ToolRunner):
    """Executes the complete UXM pipeline for a single sample"""
    base_name = task['base_name']
    logger.info(f"--- Starting UXM deconvolution for sample: {base_name} ---")

    _generate_pats(task['input_bam'], task['pat_dir'], wgbstools_runner)
    _wgbstools_pat_filter(task['pat_file_path'], task['filtered_pat_path'], atlas_path, wgbstools_runner)
    _uxm_deconvolution(uxm_runner, task['filtered_index_pat_path'], atlas_path, task['deconvoluted_file'])

    logger.info(f"--- Finished UXM deconvolution for sample: {base_name} ---")


def _generate_pats(input_bam: Path, pat_output_dir: Path, wgbstools_runner: ToolRunner):
    """
    Generates a single pat file from a single bam file
    """
    logger.info(f"Generating .pat file for {input_bam.name}")
    bam2pat_command = ['bam2pat', str(input_bam)]
    ensure_dir_exists(pat_output_dir)

    wgbstools_runner.run(bam2pat_command, pat_output_dir)


def _wgbstools_pat_filter(input_pat: Path, output_pat: Path, atlas_path: Path, wgbstools_runner: ToolRunner):
    """
    Filters the input pat file based on the atlas provided, and then indexes it.
    """

    logger.debug(f"Filtering {input_pat.name} against the atlas.")
    ensure_dir_exists(output_pat.parent)

    pat_filter_command = [
        'view',
        str(input_pat),
        '-L', str(atlas_path)
    ]
    wgbstools_runner.run(pat_filter_command, output_pat)
    logger.info(f"Pat file {input_pat} has been filtered based on {atlas_path.name}")

    index_command = ['index', str(output_pat)]
    wgbstools_runner.run(index_command)

    logger.info(f"Finished filtering and indexing: {output_pat}")


def _uxm_deconvolution(uxm_runner: ToolRunner, deconv_input: Path, atlas_path: Path, deconv_output: Path):
    """Runs the final deconvolution step with uxm."""
    logger.info(f"Running UXM deconvolution on {deconv_input.name}")
    ensure_dir_exists(deconv_output.parent)

    deconvolution_command = [
        'deconv',
        str(deconv_input),
        '--atlas', str(atlas_path)
    ]
    uxm_runner.run(deconvolution_command, deconv_output)


def _run_nnls_algorithm(nnls_exe, input_data_path, output_dir, atlas_path):
    '''
    Runs the NNLS deconvolution algorithm
    '''

    logger.info("Running NNLS deconvolution...")

    deconv_runner = ToolRunner(nnls_exe, "--out_dir")
    command = [
        "-a",
        str(atlas_path),
        str(input_data_path)
        # "--out_dir", str(output_dir)
    ]

def _run_single_nnls_deconv(task: dict, atlas_path: Path, nnls_runner: ToolRunner):
    base_name = task['base_name']
    logger.info(f"--- Starting NNLS deconvolution for sample: {base_name} ---")

    output_dir = task['deconvoluted_file']



def _generate_nnls_deconvolution_tasks(input_path: Path) -> list[dict]:
    """
    Looks at the input path and creates a list of deconvolution jobs (from bed file to deconvolution results in a
    csv file).
    """
    logger.info(f"Deconvoluting the files in {input_path} using the nnls algorithm")

    parent_dir = Path("analysis")

    tasks = []
    if not input_path.exists():
        logger.error(f"Input path does not exist: {input_path}")
        return []

    if input_path.is_dir():
        logger.info(f"Input is a directory. Searching for .bed files in {input_path}")
        bed_files = sorted(list(input_path.glob('*.bed')))
        if not bed_files:
            logger.warning(f"No .bed files found in directory: {input_path}")
            return []

        for bed_file in bed_files:
            base_name = bed_file.stem
            tasks.append({
                'input_bed': bed_file,
                'base_name': base_name,
                'sample_output_dir': parent_dir / "deconvoluted_output" / f"{base_name}_nnls.csv"
            })

    elif input_path.is_file():
        logger.info(f"Input is a single file: {input_path}")
        if input_path.suffix != '.bed':
            logger.warning(f"Input file is not a .bed file. Skipping: {input_path}")
            return []

        base_name = input_path.stem
        tasks.append({
            'input_bed': input_path,
            'base_name': base_name,
            'sample_output_dir': parent_dir / "deconvoluted_output" / f"{base_name}_nnls.csv"
        })

    return tasks



def _add_atlas_arg(parser, config):
    parser.add_argument(
        '--atlas', type=Path,
        default=None,
        metavar="<path>",
        dest="pipeline_steps.analysis.paths.atlas_dir_name",
        help="Path to the atlas used for deconvolution"
    )


def _add_algorithm_arg(parser, config):
    parser.add_argument(
        '-a', '--algorithm', type=str,
        default=None,
        dest="pipeline_steps.analysis.params.deconv_algorithm",
        choices=['uxm', 'nnls'],
        help="Algorithm to use for deconvolution."
    )


def add_all_arguments_to_parser(parser, config):
    """
    Publically available function to add deconvolution arguments to a given parser
    :param parser:
    :param config:
    :return:
    """
    _add_algorithm_arg(parser, config)
    _add_atlas_arg(parser, config)


def setup_parsers(subparsers, parent_parser, config):
    deconv_parser = subparsers.add_parser(
        'deconvolution',
        help="Deconvolute sequenced and aligned data using methylation information",
        description="This command group contains tools for deconvolution of sequenced and aligned data.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['deconv', 'analysis'],
        parents=[parent_parser]
    )
    validation_func = lambda: config.pipeline_steps.analysis.paths._validate(True)
    deconv_parser.set_defaults(validation_func=validation_func)

    _add_atlas_arg(deconv_parser, config)
    _add_algorithm_arg(deconv_parser, config)

    add_io_arguments(
        deconv_parser, config,
        default_input=None,
        input_file_help="Path to file for deconvolution.",
        input_dest="pipeline_steps.analysis.paths.user_deconv_input",
        default_output=None,
        output_dir_help="File to save the deconvolution results in.",
        output_dest="pipeline_steps.analysis.paths.user_deconv_output"
    )
    deconv_parser.set_defaults(func=deconvolution_handler)
