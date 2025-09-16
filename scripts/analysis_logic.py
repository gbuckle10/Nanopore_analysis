import pandas as pd
import os

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
    meth_df = pd.read_parquet("analysis/methylation_dataframe.parquet")

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