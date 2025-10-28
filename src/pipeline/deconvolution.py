import argparse
import logging
import sys

from src.utils.cli_utils import create_io_parser
from src.utils.config_utils import resolve_param, resolve_combined_path
from src.utils.process_utils import run_command
from pathlib import Path

from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def deconvolution_handler(args, config):
    input_data_path = resolve_combined_path(
        args, config, arg_name="input_file",
        config_path_components=[
            'paths.deconvolution_dir',
            'paths.file_to_deconvolute'
        ]
    )
    atlas_path = resolve_combined_path(
        args, config, arg_name="atlas",
        config_path_components=[
            'paths.atlas_dir',
            'paths.atlas_file_name'
        ]
    )
    output_dir = resolve_param(
        args, config, arg_name='output_dir',
        config_path='paths.deconvolution_dir'
    )
    algorithm = resolve_param(
        args, config, arg_name='algorithm',
        config_path='parameters.analysis.deconv_algorithm'
    )

    if algorithm == "uxm":
        wgbstools_exe = resolve_param(args, config, config_path='paths.submodules.wgbstools_exe')
        uxm_exe = resolve_param(args, config, config_path='paths.submodules.uxm_exe')
        _run_uxm_algorithm(wgbstools_exe, uxm_exe, input_data_path, atlas_path, output_dir)
    elif algorithm == "nnls":
        logger.info("Running deconvolution with NNLS algorithm")
        nnls_exe = resolve_param(args, config, config_path='paths.submodules.methatlas_exe')
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


def setup_parsers(subparsers, parent_parser):
    io_parser = create_io_parser()

    deconv_parser = subparsers.add_parser(
        'deconvolution',
        help="Deconvolute sequenced and aligned data using methylation information",
        description="This command group contains tools for deconvolution of sequenced and aligned data.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['deconv'],
        parents=[parent_parser, io_parser]
    )

    deconv_parser.add_argument(
        '--atlas', type=Path, default=None, help="Path to the atlas used for deconvolution"
    )
    deconv_parser.add_argument(
        '-a', '--algorithm',
        type=str, default='uxm', choices=['uxm', 'nnls']
    )
    deconv_parser.set_defaults(func=deconvolution_handler)

    def _prepare_for_nnls(self):
        """
        Prepare the data output from previous steps for basic nnls algorithm.
        """
        self.logger.info("Running preparation for meth_atlas-based nnls deconvolution algorithm.")

        bed_file_path = self.config['paths']['methylation_dir'] + self.config['paths']['methylation_bed_name']
        manifest_file_path = self.config['paths']['atlas_dir'] + self.config['paths']['illumina_manifest']
        uxm_atlas_file_path = self.config['paths']['atlas_dir'] + self.config['paths']['uxm_atlas_name']
        chunk_size = int(self.config['parameters']['analysis']['methylation_aggregation_chunksize'])

        command = [
            sys.executable,
            "-u",
            "src/deconvolution_prep.py",
            "--bed-file", bed_file_path,
            "--manifest-file", manifest_file_path,
            "--chunk-size", str(chunk_size)
        ]

        run_command(command)

    def run(self):
        """
        Main method to execute deconvolution step.
        """

        algorithm = self.config['parameters']['deconvolution']['algorithm']
        self.logger.info(f"Selected deconvolution algorithm: {algorithm}")

        if algorithm == "uxm":
            self._run_uxm_algorithm()
        elif algorithm == "nnls":
            self._run_nnls_algorithm()
        else:
            self.logger.error(f"Unknown deconvolution algorithm: {algorithm}")
            raise ValueError("Invalid deconvolution algorithm specified in config.")
