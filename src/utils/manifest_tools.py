import logging
from pathlib import Path


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
