import argparse
import sys

import pysam
import logging
from pathlib import Path
import csv
from collections import Counter


def summarise_lengths(input_file, output_file, must_be_aligned=False, min_length=0, max_length=float('inf')):
    '''
    Outputs a csv file with the number of reads of each length.
    '''

    read_lengths = Counter()
    logger = logging.getLogger('pipeline')
    logger.info(f"Summarising sequence lengths in file: {input_file}")
    prefix = Path(input_file).stem

    try:
        with pysam.AlignmentFile(input_file, "rb", check_sq=False) as infile:
            # check_sq=False tells pysam not to crash if there's no header, so it'll be more robust for unaligned bams.
            for i, read in enumerate(infile):
                if (i+1) % 1_000_000 == 0:
                    logger.info(f"Processed {i+1:,} reads...")

                if read.is_secondary or read.is_supplementary:
                    continue
                if must_be_aligned and read.is_unmapped:
                    # Skip unaligned reads.
                    continue
                if not (min_length <= (read.query_length or 0) <= max_length):
                    continue
                read_lengths[read.query_length] += 1
        logger.info(f"Process complete.")

        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow(['read_length', 'count'])
            for length, count in sorted(read_lengths.items()):
                writer.writerow([length, count])

        logger.info(f"Read length distribution saved to {output_file}")

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error processing BAM file: {e}", file=sys.stderr)
        sys.exit(1)

def setup_parsers(subparsers, parent_parser, config):

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
        default='.',
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


