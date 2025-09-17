import pandas as pd
import os
import pybedtools

def calculate_methylation_range(deconvolution_file, range_atlas_file):
    '''
    Assigns each CpG site from the methylation data to a region from the
    atlas.

    :param deconvolution_file: File for deconvolution
    :param range_atlas_file:  Atlas file with ranges
    :return:
    '''

    atlas_df = pd.read_csv(range_atlas_file, sep='\t')

    all_columns = atlas_df.columns.tolist()
    unneeded_cols = ['startCpG', 'endCpG', 'target', 'direction']
    cols_to_keep = [col for col in all_columns if
                    col not in unneeded_cols]
    atlas_df = atlas_df[cols_to_keep]

    decon_df = pd.read_csv("analysis/data_to_deconvolute_geco.csv", sep=',')
    decon_df[['chr', 'start']] = decon_df['site_id'].str.split(':', expand=True)

    print("Top 5 rows deconvolution file:")
    print(decon_df.head())

    print("Top 5 rows atlas file:")
    print(atlas_df.head)

    print("--- Converting dataframes to BedTool objects")

    decon_points_df = decon_df.copy()
    decon_points_df['end'] = decon_points_df['start'].astype(int) + 1
    decon_points_df = decon_points_df[['chr', 'start', 'end', 'percentage']]
    decon_bed = pybedtools.BedTool.from_dataframe(decon_points_df)

    atlas_regions = atlas_df[['chr', 'start', 'end']]
    atlas_bed = pybedtools.BedTool.from_dataframe(atlas_regions)

    print("--- Aggregating single-point methylation data into atlas regions ---")
    aggregated_bed = atlas_bed.map(decon_bed, c=4, o='mean', null='NaN')

    mapped_bed = decon_bed.intersect(atlas_bed, wa=True, wb=True)

    map_cols = ['decon_chr', 'decon_start', 'decon_end', 'percentage', 'atlas_chr', 'atlas_start', 'atlas_end']
    mapped_df = mapped_bed.to_dataframe(names=map_cols)

    print(f"Top 5 rows of mapped dataframe (there are {len(mapped_df)} samples mapped")
    print(mapped_df.head())

    aggregated_df = aggregated_bed.to_dataframe(
        names=['chr', 'start', 'end', 'mean_methylation']
    )

    aggregated_df.dropna(subset=['mean_methylation'], inplace=True)

    print(f"Top 5 rows of aggregated dataframe (where there are {len(aggregated_df)} total rows")
    print(aggregated_df.head())


def format_atlas_file(atlas_file, sep='\t'):
    """
    Takes an atlas file and formats the first column to work with the meth_atlas algorithm.

    At the moment this only works for the UXM atlas, eventually should allow the user to specify which
    columns are relevant.

    :param atlas_file: File path to the atlas file.
    :return:
    """
    print("Formatting the atlas file to include site_id as the first column.")

    atlas_df = pd.read_csv(
        atlas_file,
        sep=sep
    )

    print(atlas_df.head())

    all_columns = atlas_df.columns.tolist()
    cols_to_keep = [col for col in all_columns if col not in ['chr', 'start', 'end', 'startCpG', 'endCpG', 'target', 'direction']]

    atlas_df = atlas_df[cols_to_keep]

    atlas_df.to_csv("data/atlas/UXM_atlas_decon_ready.csv", index=False)

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

def generate_deconvolution_file(bed_file, manifest_file, output_file):
    """
    Filters the Nanopore methylation data, provided in a .bed file, to retain only the
    methylation sites we can use to deconvolute. This is done by retaining only those
    methylation sites which are present in the downloaded illumina manifest.

    """
    '''
    # Load the relevant information from the bed file.
    print(f'--- Loading Nanopore methylation data from: {bed_file} ---')
    meth_df = pd.read_csv(
        bed_file, sep='\t', header=None, usecols=[0, 1, 10], names=['chr', 'start', 'percentage']
    )
    
    meth_df['percentage'] = meth_df['percentage']/100.0
    meth_df['site_id'] = meth_df['chr'].astype(str) + ':' + meth_df['start'].astype(str)
    '''

    # --- MAYBE WE SHOULDN'T REMOVE THE COLUMNS WE DON'T NEED SO WE CAN MANUALLY QC ---
    meth_df = pd.read_parquet("data/methylation/methylation_dataframe.parquet")

    print("This is the first 5 rows of our methylation dataframe.")
    print(meth_df.head())


    # Load the manifest to get the list of relevant sites
    print(f"--- Loading manifest from: {manifest_file} ---")
    manifest_df = pd.read_csv(
        manifest_file, sep=',', skiprows=7, usecols=['IlmnID', 'CHR', 'MAPINFO'],
        dtype={'CHR': str}
    )

    # Remove the decimal point from MAPINFO points.
    manifest_df = manifest_df.dropna(subset=['MAPINFO'])
    manifest_df['MAPINFO'] = manifest_df['MAPINFO'].astype(int)

    manifest_df['site_id'] = 'chr' + manifest_df['CHR'].astype(str) + ':' + manifest_df['MAPINFO'].astype(str)
    manifest_df.drop(['CHR', 'MAPINFO'], axis=1, inplace=True)

    print("This is the first 5 rows of the manifest dataframe.")
    print(manifest_df.head())


    # Find the sites that are in the manifest and in our data.
    print(f"--- Finding common sites between your data and the manifest ---")

    # This filter and print is not necessary, it's only here for QC purposes atm.
    relevant_cpgs = meth_df[meth_df.site_id.isin(manifest_df.site_id.unique().tolist())]
    print(f"Your sample contains {len(relevant_cpgs)} CpGs which were also found in the manifest.")

    deconvolution_df = pd.merge(meth_df, manifest_df, on="site_id", how='inner')[['IlmnID', 'percentage']]
    deconvolution_df_gl = pd.merge(meth_df, manifest_df, on="site_id", how='inner')[['site_id', 'percentage']]

    deconvolution_df.to_parquet("analysis/data_to_deconvolute.parquet", index=False)
    deconvolution_df.to_csv("analysis/data_to_deconvolute.csv", index=False)
    deconvolution_df_gl.to_parquet("analysis/data_to_deconvolute_geco.parquet", index=False)
    deconvolution_df_gl.to_csv("analysis/data_to_deconvolute_geco.csv", index=False)

    #deconvolution_df = pd.read_parquet("data/methylation/data_to_deconvolute.parquet")

    print(f"This is the first 5 rows (out of {len(deconvolution_df)} total) of the deconvolution dataframe:")
    print(deconvolution_df.head())

    print("We are done preparing the data for deconvolution.")