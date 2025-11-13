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



def _generate_pats(input_data_path: Path, pat_dir: Path, wgbstools_runner: ToolRunner):
    """
    Generates pat files using wgbstools
    """

    bam2pat_command = [
        'bam2pat',
        str(input_data_path)
    ]
    ensure_dir_exists(pat_dir)
    wgbstools_runner.run(bam2pat_command, pat_dir)
    logger.info("Finished generating pat files.")

def _wgbstools_pat_filter(filter_input: Path, filter_output: Path, atlas_path: Path, wgbstools_runner: ToolRunner):
    """
    Goes through all of the pat files in the pat_dir and filters them based on the atlas
    """

    logger.debug(f"The file {filter_input} will be filtered and saved to {filter_output}")
    ensure_dir_exists(filter_output.parent)


    pat_filter_command = [
        # 'wgbstools',
        'view',
        str(filter_input),
        '-L', str(atlas_path)
    ]
    wgbstools_runner.run(pat_filter_command, filter_output)

    index_command = [
        'index', str(filter_output)
    ]
    wgbstools_runner.run(index_command)

def _uxm_deconvolution(uxm_runner: ToolRunner, deconv_input: Path, atlas_path: Path, deconv_output: Path):
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
    logger.info(f"Deconvolving {input_data_path.name} using atlas {atlas_path.name}")

    wgbstools_runner = ToolRunner(wgbstools_exe, '-o', handler_class=SilentHandler)

    #========================================
    # If input_data_path is a directory, it should loop through the pat files in that folder and do the bam2pat command on each one.
    # The subsequent filtering step should also loop through all of the files in pat_dir
    #========================================
    parent_dir = Path("analysis")
    pat_dir = parent_dir / "pat_files"
    filtered_pat_dir = parent_dir / "pat_filtered"
    deconvolution_dir = parent_dir / "deconvoluted_output"

    base_name = input_data_path.stem
    logger.debug(f"File base name - {base_name}")
    # The first step is to convert the bam file to a pat file.
    _generate_pats(input_data_path, pat_dir, wgbstools_runner)

    logger.debug(f"Pat files saved to {pat_dir}")

    filter_input = pat_dir / f"{base_name}.pat.gz"
    filter_output = filtered_pat_dir / f"{base_name}.pat"
    filter_index_output = filter_output.with_suffix(".pat.gz")
    logger.debug(f"Filtered index output - {filter_index_output}")
    # Then we need to filter the pat file(s) using the atlas.
    _wgbstools_pat_filter(filter_input, filter_output, atlas_path, wgbstools_runner)

    logger.info(f"Finished filtering the pat file {filter_input}")

    deconvoluted_file = deconvolution_dir / f"{base_name}.csv"
    ensure_dir_exists(deconvolution_dir)

    uxm_runner = ToolRunner(uxm_exe, '-o', handler_class=LiveDisplayHandler)
    _uxm_deconvolution(uxm_runner, filter_index_output, atlas_path, deconvoluted_file)



def _run_nnls_algorithm(nnls_exe, input_data_path, output_dir, atlas_path):
    '''
    This will contain the basic nnls algorithm for a generic atlas/bed pair.
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


