import argparse
import logging
import os
import sys

from src.utils.cli_utils import add_io_arguments
from pathlib import Path

from src.utils.file_utils import ensure_dir_exists
from src.utils.process_utils import LiveDisplayHandler, SilentHandler, TunableHandler, LineBufferingHandler
from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def deconvolution_handler(config):
    logger.info("Running deconvolution step.")
    input_file = config.pipeline_steps.analysis.paths.full_deconv_input_path
    atlas_path = config.pipeline_steps.analysis.paths.full_atlas_path
    output_dir = config.pipeline_steps.analysis.paths.full_deconv_output_path
    algorithm = config.pipeline_steps.analysis.params.deconv_algorithm

    print(f"Input file - {input_file}")
    if algorithm == "uxm":
        wgbstools_exe = config.pipeline_steps.analysis.tools.wgbstools_exe
        uxm_exe = config.pipeline_steps.analysis.tools.uxm_exe
        _run_uxm_algorithm(wgbstools_exe, uxm_exe, input_file, atlas_path, output_dir)
    elif algorithm == "nnls":
        logger.info("Running deconvolution with NNLS algorithm")
        nnls_exe = config.pipeline_steps.analysis.tools.methatlas_exe
        _run_nnls_algorithm(nnls_exe, input_file, output_dir, atlas_path)
    else:
        logger.error("You have chosen an algorithm that doesn't exist. Unfortunately I can't do this.")

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

def _generate_uxm_deconvolution_tasks(input_path: Path) -> list[dict]:
    """
    Looks at the input path and creates a list of deconvolution jobs (from aligned sorted bam file to
    deconvolution results in a csv file).
    """

    logger.info(f"Deconvoluting the files in {input_path} using the uxm algorithm")

    parent_dir = Path("analysis")
    pat_dir = parent_dir / "pat_files"
    filtered_pat_dir = parent_dir / "pat_filtered"
    deconvolution_dir = parent_dir / "deconvoluted_output"

    tasks = []
    if not input_path.exists():
        logger.error(f"Input path does not exist: {input_path}")
        return []

    # The input is a directory
    if input_path.is_dir():
        logger.info(f"The input is a directory. I am searching for .bam files in {input_path}")
        bam_files = sorted(list(input_path.glob('*.bam')))
        if not bam_files:
            logger.warning(f"No .bam files found in directory: {input_path}")
            return []

        for bam_file in bam_files:
            base_name = bam_file.stem
            tasks.append({
                'input_bam': bam_file,
                'base_name': base_name,
                'pat_dir': pat_dir,
                'pat_file_path': pat_dir / f"{base_name}.pat.gz", # We also need to know what file to look for
                'filtered_pat_path': filtered_pat_dir / f"{base_name}.pat",
                'filtered_index_pat_path': filtered_pat_dir / f"{base_name}.pat.gz",
                'deconvoluted_file': deconvolution_dir / f"{base_name}.csv"
            })

    # The input is a single file
    elif input_path.is_file():
        logger.info(f"The input is a single file: {input_path}")
        if input_path.suffix != '.bam':
            logger.warning(f"Input file is not a .bam file. Skipping: {input_path}")
            return []

        base_name = input_path.stem
        tasks.append({
                'input_bam': input_path,
                'base_name': base_name,
                'pat_dir': pat_dir,
                'pat_file_path': pat_dir / f"{base_name}.pat.gz",
                'filtered_pat_path': filtered_pat_dir / f"{base_name}.pat",
                'filtered_index_pat_path': filtered_pat_dir / f"{base_name}.pat.gz",
                'deconvoluted_file': deconvolution_dir / f"{base_name}.csv"
            })

    return tasks

def _run_single_nnls_deconv(task: dict, atlas_path: Path, nnls_runner: ToolRunner):
    base_name = task['base_name']
    logger.info(f"--- Starting NNLS deconvolution for sample: {base_name} ---")

    output_dir = task['deconvoluted_file']

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

    bam2pat_command = ['bam2pat', str(input_bam)]
    ensure_dir_exists(pat_output_dir)

    wgbstools_runner.run(bam2pat_command, pat_output_dir)
    logger.info("Finished generating pat files.")

def _wgbstools_pat_filter(input_pat: Path, output_pat: Path, atlas_path: Path, wgbstools_runner: ToolRunner):
    """
    Goes through all of the pat files in the pat_dir and filters them based on the atlas
    """

    logger.debug(f"The file {input_pat} will be filtered and saved to {output_pat}")
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
    ensure_dir_exists(deconv_output.parent)

    deconvolution_command = [
        'deconv',
        str(deconv_input),
        '--atlas', str(atlas_path)
    ]
    uxm_runner.run(deconvolution_command, deconv_output)

def _run_uxm_algorithm(wgbstools_exe, uxm_exe, input_data_path, atlas_path, output_dir):
    '''
    Runs the UXM deconvolution algorithm
    '''
    logger.info("Running UXM Deconvolution...")

    deconv_tasks = _generate_uxm_deconvolution_tasks(input_data_path)

    if not deconv_tasks:
        logger.warning("No valid deconvolution jobs found. Exiting.")
        return

    logger.info(f"There are {len(deconv_tasks)} deconvolution task(s) to do.")

    wgbstools_runner = ToolRunner(wgbstools_exe, '-o', handler_class=SilentHandler)
    uxm_runner = ToolRunner(uxm_exe, '-o', handler_class=LiveDisplayHandler)

    failed_samples = []
    for i, job in enumerate(deconv_tasks):
        logger.info(f"--- Running deconvolution task {i+1}/{len(deconv_tasks)}")
        try:
            _run_single_uxm_deconv(job, atlas_path, wgbstools_runner, uxm_runner)
        except Exception as e:
            logger.error(f"Failed to process sample {job['base_name']}. Error: {e}")
            failed_samples.append(job['base_name'])

    if failed_samples:
        logger.error(f"UXM pipeline finished with {len(failed_samples)} failed sample(s): {', '.join(failed_samples)}")
    else:
        logger.info("All UXM deconvolution jobs finished successfully.")




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


