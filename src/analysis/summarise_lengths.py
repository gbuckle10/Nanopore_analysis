import argparse
import sys

import pysam
import logging
from pathlib import Path
import csv
from collections import Counter
from src.utils.logger import setup_logger


def summarise_lengths(input_file, output_dir=None, must_be_aligned=False, min_length=0, max_length=float('inf')):
    '''
    Outputs a csv file with the number of reads of each length.

    Default behaviour is to count all reads, whether they're aligned or not.
    :param input_file:
    :param output_dir:
    :return:
    '''

    read_lengths = Counter()
    logger = logging.getLogger('pipeline')
    logger.info(f"Summarising sequence lengths in file: {input_file}")
    prefix = Path(input_file).stem
    output_file = output_dir / f"{prefix}_read_length_distribution.csv"

    try:
        with pysam.AlignmentFile(input_file, "rb", check_sq=False) as infile:
            # check_sq=False tells pysam not to crash if there's no header, so it'll be more robust for unaligned bams.
            for read in infile:
                if read.is_secondary or read.is_supplementary:
                    continue
                read_lengths[read.query_length] += 1
        logger.info(f"Process complete.")

        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow(['read_length', 'count'])
            for length, count in sorted(read_lengths.items()):
                writer.writerow([length, count])

        print(f"Distribution saved to {output_file}")

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error processing BAM file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    setup_logger()

    parser = argparse.ArgumentParser(
        description="Summarise the read length distribution in a bam file."
    )

    parser.add_argument(
        "input_bam",
        type=Path,
        help="Path to the input bam file."
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        help="Path to save output"
    )
    parser.add_argument(
        "-aligned", "--must-be-aligned",
        type=bool,
        help="Only take aligned into account"
    )
    parser.add_argument(
        "-min", "--min-length",
        type=int,
        help="Minimum length to count"
    )
    parser.add_argument(
        "-max", "--max-length",
        type=int,
        help="Maximum length to count"
    )


    args = parser.parse_args()

    if args.output_dir:
        output_directory = args.output_dir
    else:
        output_directory = args.input_bam.parent

    summarise_lengths(args.input_bam, output_directory, args.must_be_aligned, args.min_length, args.max_length)


