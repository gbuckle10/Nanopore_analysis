from pathlib import Path

import pandas as pd
import logging
import sys
import argparse

from src.utils.atlas_tools import convert_atlas_to_genome_coordinates
from src.utils.manifest_tools import get_or_create_manifest_df


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


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stdout)

    parser = argparse.ArgumentParser(
        description="Deconvolution Preparation"
    )
    parser.add_argument(
        "-b", "--bed-file",
        required=True
    )
    parser.add_argument(
        "-m", "--manifest-file",
        required=True
    )
    parser.add_argument(
        "-c", "--chunk-size",
        required=True,
        type=int)

    args = parser.parse_args()

    try:
        generate_deconvolution_files(
            bed_file=args.bed_file,
            manifest_file=args.manifest_file,
            chunk_size=args.chunk_size
        )
    except Exception as e:
        logging.error(f"An error occurred during deconvolution prep: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()


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

    #format_atlas_file("data/atlas/UXM_atlas.tsv", "data/atlas/UXM_atlas_deconvolution.csv")
    convert_atlas_to_genome_coordinates("data/atlas/full_atlas_gc.csv", "data/atlas/full_atlas.csv",
                                        "data/atlas/illumina_manifest.csv")
