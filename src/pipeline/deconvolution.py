import argparse
import logging
import sys

import pandas as pd

from src.utils.cli_utils import add_io_arguments
from pathlib import Path

from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def process_chunk(chunk_df, manifest_df):
    '''
    Takes a DataFrame chunk and processes it for the deconvolution outputs.
    Returns a dictionary of results

    :param chunk_df:
    :param manifest_df:
    :return:
    '''
    logger = logging.getLogger('pipeline')
    # Prepare chunk
    chunk_df = chunk_df[chunk_df['mod'].str.contains('m')].copy()
    if chunk_df.empty:
        logger.warning("The current chunk is empty.")
        return None

    # Beta value should be between 0 and 1.
    chunk_df['beta_value'] = chunk_df['beta_value'] / 100.0
    chunk_df['site_id'] = chunk_df['chr'].astype(str) + ':' + chunk_df['start'].astype(str)

    logger.info(f"Methylation dataframe chunk loaded and filtered. {len(chunk_df)} rows.")
    logger.info("Here are the first 5 rows of this chunk:")
    logger.info(chunk_df.head())

    # --- Illumina merges ---
    logger.info("Merging the chunk dataframe with the Illumina manifest for the meth_atlas deconvolution")
    merged_chunk = pd.merge(chunk_df, manifest_df, on='site_id', how='inner')

    logger.info(f"The chunk dataframe has been merged with the illumina manifest. Final size: {len(merged_chunk)}")
    logger.info("Here are the first 5 rows of the merged chunk:")
    logger.info(merged_chunk.head())

    # Return as a dictionary containing the 4 deconvolution-ready chunks. If the merged chunk is empty we explicity
    # return an empty dataframe to avoid an error.
    return {
        "illumina": merged_chunk[['IlmnID', 'beta_value']] if not merged_chunk.empty else pd.DataFrame(),
        "geco": merged_chunk[['site_id', 'beta_value']] if not merged_chunk.empty else pd.DataFrame(),
        "point_methylations": chunk_df[['site_id', 'beta_value']]
    }


def generate_deconvolution_files(bed_file, manifest_file, chunk_size=1000000):
    """
    Filters the Nanopore methylation data, provided in a .bed file, to retain only the
    methylation sites we can use to deconvolute. This is done by retaining only those
    methylation sites which are present in the downloaded illumina manifest.
    """
    logger = logging.getLogger('pipeline')
    processed_manifest_dir = "data/processed/illumina_manifest.parquet"

    logger.info("--- Generating files for deconvolution from the BED file ---")

    output_dir = Path("data/processed/")
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = {"illumina": [], "geco": []}

    logger.info("Creating chunk iterator")
    chunk_iterator = pd.read_csv(
        bed_file, sep='\t', header=None, usecols=[0, 1, 3, 10], names=['chr', 'start', 'mod', 'beta_value'],
        chunksize=chunk_size
    )

    logger.info(f"Getting manifest file from {manifest_file}")
    manifest_df = get_or_create_manifest_df(manifest_file, processed_manifest_dir)

    logger.info(f"Processing {bed_file} in chunks...")

    for i, chunk_df in enumerate(chunk_iterator):
        logger.info(f"  - Processing chunk {i + 1}...")

        chunk_results = process_chunk(chunk_df, manifest_df)

        if chunk_results:
            all_results["illumina"].append(chunk_results["illumina"])
            all_results["geco"].append(chunk_results["geco"])
            chunk_results["point_methylations"].to_csv("data/processed/single_site_methylations.csv", mode='a',
                                                       header=False, index=False)

    finalise_outputs(all_results, output_dir)
    # ----- ADD A CHECK HERE TO SEE WHETHER THE FILES EXIST

    # format_atlas_file("data/atlas/UXM_atlas.tsv", "data/atlas/UXM_atlas_deconvolution.csv")
    convert_atlas_to_genome_coordinates("data/atlas/full_atlas_gc.csv", "data/atlas/full_atlas.csv",
                                        "data/atlas/illumina_manifest.csv")


def finalise_outputs(results, output_dir):
    '''
    Takes the collected results and puts them into their final files.
    '''

    logger = logging.getLogger('pipeline')
    logger.info("--- Chunk loop finished. Finalising all output files.")

    for key in ["illumina", "geco"]:
        if results[key]:
            output_path = f"data/processed/deconvolution_{key}.csv"
            pd.concat(results[key], ignore_index=True).to_csv(output_path, index=False)
            logger.info(f"Saved {key} deconvolution file to {output_path}")


def get_or_create_manifest_df(manifest_csv_dir: str, processed_manifest_parquet: str) -> pd.DataFrame:
    '''
    Loads the processed manifest dataframe from a parquet file if it exists.
    Otherwise The resulting dataframe has the following columns:
        - IlmnID - Illumina ID
        - site_id - location of the CpG on the chromosome with format "chrX:YYYY"

    :param manifest_csv_dir:
    :return:
    '''

    logger = logging.getLogger('pipeline')
    processed_path = Path(processed_manifest_parquet)

    # Check whether the parquet exists already
    if processed_path.is_file():
        logger.info(f"Found an existing processed manifest. Loading from: {processed_path}")
        manifest_df = pd.read_parquet(processed_path)
    else:
        logger.info(f"Processed manifest not found. Creating from raw csv: {manifest_csv_dir}")

        # No need to chunk because it's not a huge file.
        manifest_df = pd.read_csv(
            manifest_csv_dir, sep=',', skiprows=7, usecols=['IlmnID', 'CHR', 'MAPINFO'],
            dtype={'CHR': str}
        )

        # Remove the decimal point from MAPINFO points.
        manifest_df = manifest_df.dropna(subset=['MAPINFO'])
        manifest_df['MAPINFO'] = manifest_df['MAPINFO'].astype(int)

        manifest_df['site_id'] = 'chr' + manifest_df['CHR'].astype(str) + ':' + manifest_df['MAPINFO'].astype(str)
        manifest_df.drop(['CHR', 'MAPINFO'], axis=1, inplace=True)

        logger.info("I have processed the Illumina manifest. Here are the first 5 rows of the manifest dataframe.")
        logger.info(manifest_df.head())

        # Make sure the output directory exists before saving
        processed_path.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Saving processed manifest to: {processed_path}")
        manifest_df.to_parquet(processed_path, index=False)

    logger.info(f"Manifest DataFrame ready with {len(manifest_df)} sites.")

    return manifest_df


def format_atlas_file(atlas_file, new_atlas_file, sep='\t'):
    """
    Takes an atlas file and formats the first column to work with the meth_atlas algorithm.

    At the moment this only works for the UXM atlas, eventually should allow the user to specify which
    columns are relevant.

    :param atlas_file: File path to the atlas file.
    :return:
    """

    logger = logging.getLogger('pipeline')
    logger.info("Formatting the genome coordinate atlas for deconvolution")
    logger.info(f"Atlas saved to {atlas_file}")

    atlas_df = pd.read_csv(
        atlas_file,
        sep=sep
    )

    all_columns = atlas_df.columns.tolist()
    cols_to_keep = [col for col in all_columns if
                    col not in ['chr', 'start', 'end', 'startCpG', 'endCpG', 'target', 'direction']]

    atlas_df = atlas_df[cols_to_keep]

    atlas_df.to_csv(new_atlas_file, index=False)

    logger.info(f"Atlas for genome coordinate deconvoluted saved to {atlas_file}.")


def convert_atlas_to_genome_coordinates(output_file, atlas_file, manifest_file):
    """
    Given a methylation atlas indexed by Illumina IDs, output a methylation atlas indexed by
    genome coordinates.
    :param output_file: The methylation atlas file indexed by genome coordinates.

    """

    logger = logging.getLogger('pipeline')

    logger.info(
        f"Converting Illumina atlas {atlas_file} to genome coordinate atlas {output_file} based on manifest file {manifest_file}")

    # Load the atlas and Illumina manifest
    atlas_df = pd.read_csv(
        atlas_file, sep=','
    )
    atlas_df.rename(columns={'Unnamed: 0': 'IlmnID'}, inplace=True)

    manifest_df = pd.read_csv(
        manifest_file, sep=',', skiprows=7, usecols=['IlmnID', 'CHR', 'MAPINFO'],
        dtype={'CHR': str}
    )

    # Remove the decimal point from MAPINFO points.
    manifest_df = manifest_df.dropna(subset=['MAPINFO'])
    manifest_df['MAPINFO'] = manifest_df['MAPINFO'].astype(int)

    manifest_df['site_id'] = 'chr' + manifest_df['CHR'].astype(str) + ':' + manifest_df['MAPINFO'].astype(str)
    manifest_df.drop(['CHR', 'MAPINFO'], axis=1, inplace=True)

    geco_atlas = pd.merge(atlas_df, manifest_df, on='IlmnID', how='inner')
    all_columns = geco_atlas.columns.tolist()
    cell_type_columns = [col for col in all_columns if col not in ['IlmnID', 'site_id']]
    new_column_order = ['site_id'] + cell_type_columns
    geco_atlas = geco_atlas[new_column_order]

    logger.info("The Illumina IDs in the atlas have been converted to genome locations")

    geco_atlas.to_csv(output_file, index=False)


def deconvolution_prep_handler(config):
    generate_deconvolution_files(
        bed_file=config.pipeline_steps.methylation.paths.full_bed_file_path,
        manifest_file=config.pipeline_steps.analysis.paths.full_manifest_path,
        chunk_size=config.pipeline_steps.analysis.params.methylation_aggregation_chunksize
    )

def deconvolution_handler(config):
    input_data_path = config.pipeline_steps.analysis.paths.full_deconv_input_path

    if input_data_path.suffix == ".bed" or input_data_path.name.endswith(".bed.gz"):
        logger.error(
            f"BED file input({input_data_path.name}) is not supported by 'deconv run'. "
            f"Run 'nanopore_analysis deconv prep' first to convert it, then pass the "
            f"resulting CSV as --input."
        )
        raise sys.exit(1)

    atlas_path = config.pipeline_steps.analysis.paths.full_atlas_path
    output_dir = config.pipeline_steps.analysis.paths.full_deconv_results_path
    algorithm = config.pipeline_steps.analysis.params.deconv_algorithm

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
        default=None,
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
        help="Deconvove sequenced and aligned data using methylation information",
        description="This command group contains tools for deconvolution of sequenced and aligned data.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['deconv', 'analysis'],
        parents=[parent_parser]
    )
    def show_deconv_help(config):
        """Default function to show help for the deconvolution command group."""
        deconv_parser.print_help()

    deconv_parser.set_defaults(func=show_deconv_help)

    deconv_subparsers = deconv_parser.add_subparsers(
        title="Available Commands",
        description="Choose one of the following actions to perform.",
        dest='subcommand',
        metavar="<command>"
    )

    p_run = deconv_subparsers.add_parser(
        'run', help="Deconvolve aligned sequencing data with methylation information.",
        parents=[parent_parser]
    )
    _add_atlas_arg(p_run, config)
    _add_algorithm_arg(p_run, config)
    add_io_arguments(
        p_run, config,
        default_input=None,
        input_file_help="Path to file for deconvolution.",
        input_dest="pipeline_steps.analysis.paths.deconvolution_input_name",
        default_output=None,
        output_dir_help="File to save the deconvolution results in.",
        output_dest="pipeline_steps.analysis.paths.deconvolution_results_name"
    )
    p_run.set_defaults(func=deconvolution_handler)

    p_prep = deconv_subparsers.add_parser(
        'prep', help="Prepare methylation data for deconvolution",
        parents=[parent_parser]
    )

    p_prep.add_argument(
        '-b', '--bed-file',
        type=Path,
        default=config.pipeline_steps.methylation.paths.full_bed_file_path,
        dest="pipeline_steps.methylation.paths.methylation_bed_name",
        metavar="PATH",
        help="Path to modkit BED file."
    )
    p_prep.add_argument(
        '-m', '--manifest-file',
        type=Path,
        default=config.pipeline_steps.analysis.paths.full_manifest_path,
        dest="pipeline_steps.analysis.paths.manifest_name",
        metavar="PATH",
        help="Path to the Illumina manifest file."
    )
    p_prep.add_argument(
        '-c', '--chunk-size',
        type=int,
        default=config.pipeline_steps.analysis.params.methylation_aggregation_chunksize,
        dest="pipeline_steps.analysis.params.methylation_aggregation_chunksize",
        metavar="INT",
        help="Number of rows to process at a time."
    )
    p_prep.set_defaults(func=deconvolution_prep_handler)
