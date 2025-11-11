import logging
import pandas as pd

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
