import argparse
import logging
import sys

from src.utils.cli_utils import add_io_arguments
from pathlib import Path

from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def deconvolution_handler(args, config):
    input_data_path = args.input_file
    atlas_path = args.atlas
    output_dir = args.output_dir
    algorithm = args.algorithm


    if algorithm == "uxm":
        wgbstools_exe = config.pipeline_steps.analysis.tools.wgbstools_exe
        uxm_exe = config.pipeline_steps.analysis.tools.uxm_exe
        _run_uxm_algorithm(wgbstools_exe, uxm_exe, input_data_path, atlas_path, output_dir)
    elif algorithm == "nnls":
        logger.info("Running deconvolution with NNLS algorithm")
        nnls_exe = config.pipeline_steps.analysis.tools.methatlas_exe
        _run_nnls_algorithm(nnls_exe, input_data_path, output_dir, atlas_path)
    else:
        logger.error("You have chosen an algorithm that doesn't exist. Unfortunately I can't do this.")


def _run_uxm_algorithm(wgbstools_exe, uxm_exe, input_data_path, atlas_path, output_dir):
    '''
    Runs the UXM deconvolution algorithm
    '''
    logger.info("Running UXM Deconvolution...")
    logger.info(f"Deconvolving {input_data_path.name} using atlas {atlas_path.name}")

    pat_index_suffix = '.pat.gz'
    pat_suffix = '.pat'
    input_base_name = input_data_path.name.removesuffix(pat_index_suffix)
    filtered_pat_name = f"{input_base_name}_tempnew{pat_suffix}"
    indexed_pat_name = f"{input_base_name}_tempnew{pat_index_suffix}"
    filtered_pat_dir = input_data_path.with_name(filtered_pat_name)
    indexed_pat_dir = input_data_path.with_name(indexed_pat_name)
    logger.info(f"Filtering input based on the atlas, filtered file {filtered_pat_dir}")

    wgbstools_runner = ToolRunner(wgbstools_exe, '-o')

    pat_filter_command = [
        # 'wgbstools',
        'view',
        str(input_data_path),
        '-L', str(atlas_path)
    ]

    logger.info(f"Filtering the pat file based on the atlas using command: {' '.join(pat_filter_command)}")

    wgbstools_runner.run(pat_filter_command, filtered_pat_dir)

    logger.info(f"Finished filtering the pat based on the atlas. Indexing...")
    pat_index_command = [
        # 'wgbstools',
        'index',
        str(filtered_pat_dir)
    ]
    wgbstools_runner.run(pat_index_command)

    logger.info(f"Finished indexing the pat. Deconvoluting...")

    uxm_runner = ToolRunner(uxm_exe, '--output')

    deconvolution_command = [
        # 'uxm',
        'deconv',
        str(indexed_pat_dir),
        '--atlas', str(atlas_path)
        # '--output', str(output_dir),
    ]
    uxm_runner.run(deconvolution_command, output_dir)


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
        default=config.pipeline_steps.analysis.paths.full_atlas_path,
        help="Path to the atlas used for deconvolution"
    )
def _add_algorithm_arg(parser, config):
    parser.add_argument(
        '-a', '--algorithm', type=str,
        default=config.pipeline_steps.analysis.params.deconv_algorithm, choices=['uxm', 'nnls']
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
        aliases=['deconv'],
        parents=[parent_parser]
    )
    _add_atlas_arg(deconv_parser, config)
    _add_algorithm_arg(deconv_parser, config)

    add_io_arguments(
        deconv_parser, config,
        default_input=config.pipeline_steps.analysis.paths.full_deconv_input_path,
        default_output=config.pipeline_steps.analysis.paths.full_deconv_results_path,
        input_file_help="Path to file for deconvolution.",
        output_dir_help="File to save the deconvolution results in."
    )
    deconv_parser.set_defaults(func=deconvolution_handler)
