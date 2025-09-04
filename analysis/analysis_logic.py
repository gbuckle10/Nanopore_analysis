import pandas as pd

def make_deconvolution_file(bed_file, manifest_file, output_file):
    """
    Filters the Nanopore methylation data, provided in a .bed file, to retain only the
    methylation sites we can use to deconvolute. This is done by retaining only those
    methylation sites which are present in the downloaded illumina manifest.

    """

    # Load the relevant information from the bed file.
    print(f"--- Loading Nanopore methylation data from: {bed_file} ---")
    meth_df = pd.read_csv(
        bed_file, sep='\t', header=None, usecols=[0, 1, 4], names=['chr', 'start', 'percentage']
    )
    meth_df['site_id'] = meth_df['chr'].astype(str) + ':' + meth_df['start'].astype(str)
    meth_series = meth_df.set_index('site_id')['percentage']/100.00

    #print(meth_series)

    # Load the manifest to get the list of relevant sites
    print(f"--- Loading manifest from: {manifest_file} ---")
    manifest_df = pd.read_csv(
        manifest_file, sep='\t', skip_rows=7, usecols=['IlmID', 'CHR', 'MAPINFO']
    )

    manifest_df['site_id'] = 'chr' + manifest_df['CHR'] + ':' + manifest_df['MAPINFO'].astype(str)

    manifest_sites = set(manifest_df['site_id'])

    # Find the sites that are in the manifest and in our data.
    print(f"--- Finding common sites between your data and the manifest ---")
