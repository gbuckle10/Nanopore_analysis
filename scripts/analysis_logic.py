import pandas as pd
import os
import pybedtools

def load_meth_df_from_bed(bed_dir):
    print(f'Loading Nanopore methylation data from: {bed_dir} ... ')

    meth_df = pd.read_csv(
        bed_dir, sep='\t', header=None, usecols=[0, 1, 10], names=['chr', 'start', 'percentage']
    )
    meth_df['percentage'] = meth_df['percentage'] / 100.0

    meth_df['site_id'] = meth_df['chr'].astype(str) + ':' + meth_df['start'].astype(str)

    print("Methylation dataframe loaded from the bed file.")
    print("Here are the first 5 rows of our methylation dataframe:")
    print(meth_df.head())

    return meth_df

def manifest_df_from_csv(manifest_dir):
    # Load the manifest to get the list of relevant sites
    print(f"--- Loading manifest from: {manifest_dir} ---")
    manifest_df = pd.read_csv(
        manifest_dir, sep=',', skiprows=7, usecols=['IlmnID', 'CHR', 'MAPINFO'],
        dtype={'CHR': str}
    )

    # Remove the decimal point from MAPINFO points.
    manifest_df = manifest_df.dropna(subset=['MAPINFO'])
    manifest_df['MAPINFO'] = manifest_df['MAPINFO'].astype(int)

    manifest_df['site_id'] = 'chr' + manifest_df['CHR'].astype(str) + ':' + manifest_df['MAPINFO'].astype(str)
    manifest_df.drop(['CHR', 'MAPINFO'], axis=1, inplace=True)

    print("This is the first 5 rows of the manifest dataframe.")
    print(manifest_df.head())

    return manifest_df

def generate_illumina_deconvolution_file(methylation_df, manifest_df):
    '''
    Generate a file to be deconvoluted with the meth_atlas deconvolution algorithm.
    '''

    illumina_deconvolution_df = pd.merge(methylation_df, manifest_df, on="site_id", how='inner')[
        ['IlmnID', 'percentage']]

    #illumina_deconvolution_df.to_parquet("analysis/data_to_deconvolute.parquet", index=False)
    illumina_deconvolution_df.to_csv("analysis/data_to_deconvolute_illumina.csv", index=False)

def generate_gc_illumina_deconvolution_file(methylation_df, manifest_df):
    '''
    Using an Illumina manifest and genome coordinate methylation dataframe, make a deconvolution-ready file
    matching the Illumina probe IDs to genome locations.
    '''

    # Only keep the rows for which the location matches an Illumina probe
    illumina_deconvolution_df_gl = pd.merge(methylation_df, manifest_df, on="site_id", how='inner')[
        ['site_id', 'percentage']]
    # illumina_deconvolution_df_gl.to_parquet("analysis/data_to_deconvolute_geco.parquet", index=False)
    illumina_deconvolution_df_gl.to_csv("analysis/data_to_deconvolute_geco.csv", index=False)

def generate_uxm_deconvolution_file(methylation_df):

    uxm_deconvolution_df = methylation_df[['site_id', 'percentage']]

    #uxm_deconvolution_df.to_parquet("analysis/data_to_deconvolute_uxm.parquet", index=False)
    uxm_deconvolution_df.to_csv("analysis/data_to_deconvolute_uxm.csv", index=False)

def generate_deconvolution_files(bed_file, manifest_file, output_file, range_atlas_file, chunk_size=1000000):
    """
    Filters the Nanopore methylation data, provided in a .bed file, to retain only the
    methylation sites we can use to deconvolute. This is done by retaining only those
    methylation sites which are present in the downloaded illumina manifest.

    """
    print("--- Generating files for deconvolution from the BED file ---")
    # Load the relevant information from the bed file.

    methylation_df = load_meth_df_from_bed(bed_file)
    #methylation_df = pd.read_parquet("data/methylation/methylation_dataframe.parquet")

    manifest_df = manifest_df_from_csv(manifest_file)

    # Find the sites that are in the manifest and in our data.
    print(f"--- Prepping files for deconvolution ---")

    generate_illumina_deconvolution_file(methylation_df, manifest_df)

    generate_uxm_deconvolution_file(methylation_df)


    generate_aggregated_deconvolution_file("analysis/data_to_deconvolute_uxm.csv", range_atlas_file, chunk_size)

    #print(f"This is the first 5 rows (out of {len(illumina_deconvolution_df)} total) of the deconvolution dataframe:")
    #print(illumina_deconvolution_df.head())

    print("We are done preparing the data for deconvolution.")



def generate_aggregated_deconvolution_file(methylation_file, range_atlas_file, chunk_size):
        '''
        Assigns each CpG site from the methylation data to a region from the
        atlas.

        :param deconvolution_file: Methylation dataframe containing columns 'chr', 'loc' and 'percentage'
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
        print(atlas_df.head)

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
            chunk_df_pbt = chunk_df[['chr', 'start', 'end', 'percentage']]
            chunk_bed = pybedtools.BedTool.from_dataframe(chunk_df_pbt)

            mapped_chunk_bed = chunk_bed.intersect(atlas_bed, wa=True, wb=True)

            if len(mapped_chunk_bed) > 0:
                map_cols = ['meth_chr', 'meth_start', 'meth_end', 'percentage',
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
                mean_methylation=('percentage', 'mean'),
                cpg_count=('percentage', 'count')
            ).reset_index()

            aggregated_df.rename(columns={'percentage': 'mean_methylation'}, inplace=True)

            print("--- Formatting final output and cleaning up ---")
            aggregated_df['site_id'] = aggregated_df['atlas_chr'].astype(str) + ':' + aggregated_df[
                'atlas_start'].astype(str) + '-' + aggregated_df['atlas_end'].astype(str)

            print("Aggregated df head:")
            print(aggregated_df.head())

            final_df = aggregated_df[['site_id', 'atlas_name', 'mean_methylation', 'cpg_count', 'atlas_chr', 'atlas_start', 'atlas_end']]
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

