import logging
import sys

from scripts.utils.runner import run_command
from pathlib import Path

class Deconvolution:
    def __init__(self, config, tool_paths, atlas_path, input_data_path, output_dir):
        """
        Initialise deconvolution handler.
        """
        self.config = config
        self.tool_paths = tool_paths
        self.logger = logging.getLogger('pipeline')
        self.atlas_path = atlas_path
        self.input_data_path = input_data_path
        self.output_dir = output_dir


    def _run_uxm_algorithm(self):
        '''
        Runs the UXM deconvolution algorithm
        '''
        self.logger.info("Running UXM Deconvolution...")
        self.logger.info(f"Deconvolving {self.input_data_path.name} using atlas {self.atlas_path.name}")

        pat_index_suffix = '.pat.gz'
        pat_suffix = '.pat'
        input_base_name = self.input_data_path.name.removesuffix(pat_index_suffix)
        filtered_pat_name = f"{input_base_name}_tempnew{pat_suffix}"
        indexed_pat_name = f"{input_base_name}_tempnew{pat_index_suffix}"
        filtered_pat_dir = self.input_data_path.with_name(filtered_pat_name)
        indexed_pat_dir = self.input_data_path.with_name(indexed_pat_name)
        self.logger.info(f"Filtering input based on the atlas, filtered file {filtered_pat_dir}")

        pat_filter_command = [
            'wgbstools', 'view',
            str(self.input_data_path),
            '-L', str(self.atlas_path),
            '-o', str(filtered_pat_dir)
        ]
        #run_command(pat_filter_command)

        self.logger.info(f"Finished filtering the pat based on the atlas. Indexing...")
        pat_index_command = [
            'wgbstools', 'index',
            str(filtered_pat_dir)
        ]
        #run_command(pat_index_command)

        self.logger.info(f"Finished indexing the pat. Deconvoluting...")

        deconvolution_command = [
            'uxm', 'deconv',
            str(indexed_pat_dir),
            '--atlas', str(self.atlas_path),
            '--output', str(self.output_dir),
        ]
        run_command(deconvolution_command)

    def _run_nnls_algorithm(self):
        '''
        This will contain the basic nnls algorithm for a generic atlas/bed pair.
        '''
        self.logger.info("Running NNLS deconvolution...")
        deconvolve_script = "externals/meth_atlas/deconvolve_genome_coordinates.py"

        command = [
            sys.executable,
            deconvolve_script,
            "-a",
            str(self.atlas_path),
            str(self.input_data_path),
            "--out_dir", str(self.output_dir)
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