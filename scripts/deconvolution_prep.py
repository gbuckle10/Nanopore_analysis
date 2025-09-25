import pandas as pd
import os
import pybedtools
import logging
from pathlib import Path
from scripts.utils.file_conversion import *


def process_chunk(chunk_df, manifest_df, uxm_atlas_bed):
    '''
    Takes a DataFrame chunk and processes it for the deconvolution outputs.
    Returns a dictionary of results

    :param chunk_df:
    :param manifest_df:
    :param uxm_atlas_bed:
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

    # UXM Aggregation
    chunk_df['end'] = chunk_df['start'] + 1
    chunk_bed = pybedtools.BedTool.from_dataframe(chunk_df[['chr', 'start', 'end', 'beta_value']])
    intersections = chunk_bed.intersect(uxm_atlas_bed, wa=True, wb=True)

    # Return as a dictionary containing the 4 deconvolution-ready chunks. If the merged chunk is empty we explicity
    # return an empty dataframe to avoid an error.
    return {
        "illumina": merged_chunk[['IlmnID', 'beta_value']] if not merged_chunk.empty else pd.DataFrame(),
        "geco": merged_chunk[['site_id', 'beta_value']] if not merged_chunk.empty else pd.DataFrame(),
        "point_methylations": chunk_df[['site_id', 'beta_value']],
        "uxm_intersections": intersections
    }


def generate_deconvolution_files(bed_file, manifest_file, uxm_atlas_file, chunk_size=1000000):
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

    illumina_deconvolution_path = output_dir / "data_to_deconvolute_illumina.csv"
    gc_illumina_deconvolution_path = output_dir / "data_to_deconvolute_geco.csv"
    uxm_deconvolution_path = output_dir / "data_to_deconvolute_uxm.csv"
    uxm_aggregated_path = output_dir / "data_to_deconvolute_uxm_aggregated.csv"
    uxm_aggregated_intermediate = output_dir / "temp_uxm_intersections.tsv"

    all_results = {"illumina": [], "geco": []}

    # Prepare UXM atlas for pybedtools
    logger.info(f"Putting UXM atlas at {uxm_atlas_file} into a BedTool object")
    uxm_atlas_bed = pybedtools.BedTool.from_dataframe(
        pd.read_csv(uxm_atlas_file, sep='\t', usecols=['chr', 'start', 'end', 'name'])
    )

    logger.info("Creating chunk iterator")
    chunk_iterator = pd.read_csv(
        bed_file, sep='\t', header=None, usecols=[0, 1, 3, 10], names=['chr', 'start', 'mod', 'beta_value'],
        chunksize=chunk_size
    )

    logger.info(f"Getting manifest file from {manifest_file}")
    manifest_df = get_or_create_manifest_df(manifest_file, processed_manifest_dir)

    uxm_aggregated_header_written = False
    logger.info(f"Processing {bed_file} in chunks...")

    for i, chunk_df in enumerate(chunk_iterator):
        logger.info(f"  - Processing chunk {i + 1}...")

        chunk_results = process_chunk(chunk_df, manifest_df, uxm_atlas_bed)

        if chunk_results:
            all_results["illumina"].append(chunk_results["illumina"])
            all_results["geco"].append(chunk_results["geco"])
            chunk_results["point_methylations"].to_csv("data/processed/single_site_methylations.csv", mode='a',
                                                       header=False, index=False)

            intersections = chunk_results["uxm_intersections"]
            if len(intersections) > 0:
                map_cols = ['meth_chr', 'meth_start', 'meth_end', 'beta_value',
                            'atlas_chr', 'atlas_start', 'atlas_end', 'atlas_name']
                mapped_df = intersections.to_dataframe(names=map_cols)
                if not uxm_aggregated_header_written:
                    mapped_df.to_csv(uxm_aggregated_intermediate, sep='\t', index=False, header=True)
                    uxm_aggregated_header_written = True
                else:
                    mapped_df.to_csv(uxm_aggregated_intermediate, sep='\t', index=False, header=False, mode='a')

    finalise_outputs(all_results, output_dir, uxm_aggregated_intermediate)

    format_atlas_file("data/atlas/UXM_atlas.tsv")


def finalise_outputs(results, output_dir, agg_intermediate_file):
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

    if agg_intermediate_file.exists():
        logger.info("Aggregating UXM intersection data.")

        overview_path = "data/processed/uxm_overview.csv"
        deconvolution_path = "data/processed/uxm_deconvolution.csv"

        mapped_df = pd.read_csv(agg_intermediate_file, sep='\t')

        aggregated_df = mapped_df.groupby(
            ['atlas_chr', 'atlas_start', 'atlas_end', 'atlas_name']
        ).agg(
            mean_methylation=('beta_value', 'mean'),
            cpg_count=('beta_value', 'count')
        ).reset_index(
        )

        aggregated_df.rename(columns={'beta_value': 'mean_methylation'}, inplace=True)

        logger.info("--- Formatting final output and cleaning up ---")
        aggregated_df['site_id'] = aggregated_df['atlas_chr'].astype(str) + ':' + aggregated_df[
            'atlas_start'].astype(str) + '-' + aggregated_df['atlas_end'].astype(str)

        logger.info("Aggregated df head:")
        logger.info(aggregated_df.head())

        final_df = aggregated_df[
            ['site_id', 'atlas_name', 'mean_methylation', 'cpg_count', 'atlas_chr', 'atlas_start', 'atlas_end']]
        final_deconvolution_df = aggregated_df[['atlas_name', 'mean_methylation']]

        logger.info(f"Top 5 rows of final aggregated data (out of {len(final_df)}:")
        logger.info(final_df.head())

        # Save final aggregrated result
        final_df.to_csv(overview_path, sep=',', index=False)
        logger.info(f"Final aggregated data saved to {overview_path}")

        final_deconvolution_df.to_csv(deconvolution_path, sep=',', index=False)
        logger.info(f"Final aggregated deconvolution file saved to {deconvolution_path}")


def load_meth_df_from_bed(bed_dir):
    '''
        Loads a bed or bedmethyl file into a pandas dataframe for further processing.
        :param bed_dir: filepath of bed(methyl) file to be loaded
        :return:
    '''
    logger = logging.getLogger('pipeline')
    logger.info(f"Loading and filtering Nanopore methylation data from: {bed_dir}...")
    # Load into a dataframe. Col 3 is taken so that we can filter out everything apart from "m" modification.
    # This needs to be more selective, it only works if we have a bedmethyl file but if it's a different type of bed
    # file it won't work.

    bed_cols = ['chr', 'start', 'mod', 'beta_value']
    chunk_size = 5_000_000

    chunk_iterator = pd.read_csv(
        bed_dir,
        sep='\t',
        header=None,
        usecols=[0, 1, 3, 10],  # read only chr, start, mod and beta value
        names=bed_cols,
        chunksize=chunk_size
    )

    filtered_chunks = []
    for i, chunk in enumerate(chunk_iterator):
        logger.info(f" - Processing chunk {i + 1}...")
        filtered_chunk = chunk[chunk['mod'].str.contains('m')]
        filtered_chunks.append(filtered_chunk)

    logger.info("Combining filtered chunks...")
    meth_df = pd.concat(filtered_chunks, ignore_index=True)

    # Beta value should be between 0 and 1
    meth_df['beta_value'] = meth_df['beta_value'] / 100.0
    meth_df['site_id'] = meth_df['chr'].astype(str) + ':' + meth_df['start'].astype(str)

    logger.info(f"Methylation dataframe loaded and filtered. Final size: {len(meth_df)} rows.")
    logger.info("Here are the first 5 rows of our methylation dataframe:")
    logger.info(meth_df.head())

    logger.info("Saving meth_df to csv and parquet...")

    meth_df.to_parquet("data/methylation/methylation_dataframe.parquet")
    meth_df.to_csv("data/methylation/methylation_dataframe.csv", index=False)

    return meth_df


def load_meth_df_from_bed__(bed_dir):
    '''
    Loads a bed or bedmethyl file into a pandas dataframe for further processing.
    :param bed_dir: filepath of bed(methyl) file to be loaded
    :return:
    '''

    print(f'Loading Nanopore methylation data from: {bed_dir} ... ')
    # Load into a dataframe. Col 3 is taken so that we can filter out everything apart from "m" modification.
    # This needs to be more selective, it only works if we have a bedmethyl file but if it's a different type of bed
    # file it won't work.

    meth_df = pd.read_csv(
        bed_dir, sep='\t', header=None, usecols=[0, 1, 3, 10], names=['chr', 'start', 'mod', 'beta_value']
    )

    # Only keep
    meth_df = meth_df[meth_df['mod'] == 'm'].copy()
    meth_df.drop(columns=['mod'], inplace=True
                 )
    meth_df['beta_value'] = meth_df['beta_value'] / 100.0

    meth_df['site_id'] = meth_df['chr'].astype(str) + ':' + meth_df['start'].astype(str)

    print("Methylation dataframe loaded from the bed file.")
    print("Here are the first 5 rows of our methylation dataframe:")
    print(meth_df.head())

    meth_df.to_parquet("data/methylation/methylation_dataframe.parquet")
    meth_df.to_csv("data/methylation/methylation_dataframe.csv", index=False)

    return meth_df


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


def generate_illumina_deconvolution_file(methylation_df, manifest_df):
    '''
    Generate a file to be deconvoluted with the meth_atlas deconvolution algorithm, in which the location is based on
    Illumina ID.
    The output file will contain 2 columns:
        - IlmnID - The Illumina ID of the CpG location
        - beta_value - The beta value for that specific CpG
    '''

    logger = logging.getLogger('pipeline')

    illumina_deconvolution_df = pd.merge(methylation_df, manifest_df, on="site_id", how='inner')[
        ['IlmnID', 'beta_value']]

    # illumina_deconvolution_df.to_parquet("analysis/data_to_deconvolute.parquet", index=False)
    illumina_deconvolution_df.to_csv("analysis/data_to_deconvolute_illumina.csv", index=False)


def generate_gc_illumina_deconvolution_file(methylation_df, manifest_df):
    '''
    Using an Illumina manifest and genome coordinate methylation dataframe, make a deconvolution-ready file
    matching the Illumina probe IDs to genome locations.
    '''

    # Only keep the rows for which the location matches an Illumina probe
    illumina_deconvolution_df_gl = pd.merge(methylation_df, manifest_df, on="site_id", how='inner')[
        ['site_id', 'beta_value']]
    # illumina_deconvolution_df_gl.to_parquet("analysis/data_to_deconvolute_geco.parquet", index=False)
    illumina_deconvolution_df_gl.to_csv("analysis/data_to_deconvolute_geco.csv", index=False)


def generate_uxm_deconvolution_file(methylation_df):
    print("Generating uxm deconvolution file from dataframe.")
    print("Top 5 rows of methylation dataframe:")
    print(methylation_df.head())
    uxm_deconvolution_df = methylation_df[['site_id', 'beta_value']]

    # uxm_deconvolution_df.to_parquet("analysis/data_to_deconvolute_uxm.parquet", index=False)
    uxm_deconvolution_df.to_csv("analysis/data_to_deconvolute_big.csv", index=False)


def generate_aggregated_deconvolution_file(methylation_file, range_atlas_file, chunk_size):
    '''
    Assigns each CpG site from the methylation data to a region from the
    atlas.

    :param deconvolution_file: Methylation dataframe containing columns 'chr', 'loc' and 'beta_value'
    :param range_atlas_file:  Atlas file with location ranges
    :return:
    '''

    intermediate_file = "analysis/temp_mapped_methylation.tsv"
    output_data_file = "analysis/aggregated_deconvolution_overview.csv"
    output_deconvolution_file = "analysis/aggregated_deconvolution.csv"
    print("--- Starting to make the aggregated deconvolution file ---")
    atlas_df = pd.read_csv(range_atlas_file, sep='\t')

    all_columns = atlas_df.columns.tolist()
    unneeded_cols = ['startCpG', 'endCpG', 'target', 'direction']
    cols_to_keep = [col for col in all_columns if
                    col not in unneeded_cols]
    atlas_df = atlas_df[cols_to_keep]
    print("Top 5 rows atlas file:")
    print(atlas_df.head())

    atlas_regions = atlas_df[['chr', 'start', 'end', 'name']]
    atlas_bed = pybedtools.BedTool.from_dataframe(atlas_regions)

    # Read in chunks
    print(f"Reading methylation file in chunks of size {chunk_size}")

    chunk_iterator = pd.read_csv(
        methylation_file,
        chunksize=chunk_size
    )

    header_written = False
    for i, chunk_df in enumerate(chunk_iterator):
        print(f" - Processing chunk {i + 1} - row {i * chunk_size} to {(i + 1) * chunk_size}...")
        chunk_df[['chr', 'start']] = chunk_df['site_id'].str.split(':', expand=True)
        chunk_df['end'] = chunk_df['start'].astype(int) + 1
        chunk_df_pbt = chunk_df[['chr', 'start', 'end', 'beta_value']]
        chunk_bed = pybedtools.BedTool.from_dataframe(chunk_df_pbt)

        mapped_chunk_bed = chunk_bed.intersect(atlas_bed, wa=True, wb=True)

        if len(mapped_chunk_bed) > 0:
            map_cols = ['meth_chr', 'meth_start', 'meth_end', 'beta_value',
                        'atlas_chr', 'atlas_start', 'atlas_end', 'atlas_name']
            mapped_df = mapped_chunk_bed.to_dataframe(names=map_cols)

            if not header_written:
                mapped_df.to_csv(intermediate_file, sep='\t', index=False, header=True)
                header_written = True
            else:
                mapped_df.to_csv(intermediate_file, sep='\t', index=False, header=False, mode='a')

    print("--- All chunks processed ---")
    if not os.path.exists(intermediate_file):
        print("No overlaps found. Exiting...")
        open(output_data_file, 'w').close()
    else:
        final_mapped_df = pd.read_csv(intermediate_file, sep='\t')
        print(f"we did it! First few rows of the final mapped df (out of {len(final_mapped_df)}:")
        print(final_mapped_df.head())

        aggregated_df = final_mapped_df.groupby(
            ['atlas_chr', 'atlas_start', 'atlas_end', 'atlas_name']
        ).agg(
            mean_methylation=('beta_value', 'mean'),
            cpg_count=('beta_value', 'count')
        ).reset_index()

        aggregated_df.rename(columns={'beta_value': 'mean_methylation'}, inplace=True)

        print("--- Formatting final output and cleaning up ---")
        aggregated_df['site_id'] = aggregated_df['atlas_chr'].astype(str) + ':' + aggregated_df[
            'atlas_start'].astype(str) + '-' + aggregated_df['atlas_end'].astype(str)

        print("Aggregated df head:")
        print(aggregated_df.head())

        final_df = aggregated_df[
            ['site_id', 'atlas_name', 'mean_methylation', 'cpg_count', 'atlas_chr', 'atlas_start', 'atlas_end']]
        final_deconvolution_df = aggregated_df[['atlas_name', 'mean_methylation']]

        print(f"Top 5 rows of final aggregated data (out of {len(final_df)}:")
        print(final_df.head())

        # Save final aggregrated result
        final_df.to_csv(output_data_file, sep=',', index=False)
        print(f"Final aggregated data saved to {output_data_file}")

        final_deconvolution_df.to_csv(output_deconvolution_file, sep=',', index=False)
        print(f"Final aggregated deconvolution file saved to {output_deconvolution_file}")


def format_atlas_file(atlas_file, sep='\t'):
    """
    Takes an atlas file and formats the first column to work with the meth_atlas algorithm.

    At the moment this only works for the UXM atlas, eventually should allow the user to specify which
    columns are relevant.

    :param atlas_file: File path to the atlas file.
    :return:
    """

    logger = logging.getLogger('pipeline')
    logger.log_info("Formatting the genome coordinate atlas for deconvolution")
    logger.log_info(f"Atlas saved to {atlas_file}")

    atlas_df = pd.read_csv(
        atlas_file,
        sep=sep
    )

    all_columns = atlas_df.columns.tolist()
    cols_to_keep = [col for col in all_columns if
                    col not in ['chr', 'start', 'end', 'startCpG', 'endCpG', 'target', 'direction']]

    atlas_df = atlas_df[cols_to_keep]

    atlas_df.to_csv(atlas_file, index=False)

    logger.log_info(f"Atlas for genome coordinate deconvoluted saved to {atlas_file}.")


def convert_atlas_to_genome_coordinates(output_file, atlas_file, manifest_file):
    """
    Given a methylation atlas indexed by Illumina IDs, output a methylation atlas indexed by
    genome coordinates.
    :param output_file: The methylation atlas file indexed by genome coordinates.

    """

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

    print("The top 5 rows of the atlas file:")
    print(atlas_df.head())

    print("The top 5 rows of the manifest file:")
    print(manifest_df.head())

    print(f"The atlas currently contains {len(atlas_df)} rows")
    print(f"The manifest df currently contains {len(manifest_df)} rows")

    geco_atlas = pd.merge(atlas_df, manifest_df, on='IlmnID', how='inner')
    all_columns = geco_atlas.columns.tolist()
    cell_type_columns = [col for col in all_columns if col not in ['IlmnID', 'site_id']]
    new_column_order = ['site_id'] + cell_type_columns
    geco_atlas = geco_atlas[new_column_order]

    print("The top 5 rows of the genomic location atlas:")
    print(geco_atlas.head())

    print(f"The genome location manifest currently contains {len(geco_atlas)} rows")

    geco_atlas.to_csv(output_file, index=False)

