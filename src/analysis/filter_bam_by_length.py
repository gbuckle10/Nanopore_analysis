import argparse
import pysam
import logging
from pathlib import Path
import sys
import subprocess



def filter_bam_by_length(input_file, size_cutoff, output_dir, side_selection):
    '''
    Takes a bam file and filters based on the size cutoff. A new bam file is made based on the filtering,
    with the option to save the reads above the cutoff, below the cutoff, or both.

    The default behaviour is to output in the same directory as the input file and to save both above and below.
    :param input_file:
    :param size_cutoff:
    :param output_dir:
    :return:
    '''

    logger = logging.getLogger('pipeline')
    logger.info(f"Splitting file: {input_file}")
    logger.info(f"Read size cutoff - {size_cutoff}")
    logger.info(f"Outputting file(s) to {output_dir}")
    logger.info(f"Selection mode: {side_selection}")
    prefix = Path(input_file).stem

    output_above = None
    output_below = None
    count_above = 0
    count_below = 0
    count_total = 0

    file_above = output_dir / f"{prefix}_above_{size_cutoff}bp.bam"
    file_below = output_dir / f"{prefix}_below_{size_cutoff}bp.bam"

    try:
        with pysam.AlignmentFile(input_file, "rb") as infile:
            if side_selection in ["both", "above"]:
                logger.info(f"Outputting the reads above {size_cutoff} bp to {file_above}")
                output_above = pysam.AlignmentFile(file_above, "wb", template=infile)
            if side_selection in ["both", "below"]:
                logger.info(f"Outputting the reads below {size_cutoff} bp to {file_below}")
                output_below = pysam.AlignmentFile(file_below, "wb", template=infile)

            for read in infile:
                count_total += 1
                read_len = read.query_length

                if read_len > size_cutoff:
                    if output_above:
                        output_above.write(read)
                        count_above += 1
                else:
                    if output_below:
                        output_below.write(read)
                        count_below += 1

            logger.info("--- Filtering Complete ---")
            logger.info(f"Total reads processed - {count_total}")
            if output_above:
                print(f"Reads > {size_cutoff} bp: {count_above}")
            if output_below:
                print(f"Reads < {size_cutoff} bp: {count_below}")

        logger.info("--- Indexing output files ---")
        if side_selection in ['above', 'both']:
            logger.info("Indexing {output_above}")
            index_above_cmd = ['samtools', 'index', file_above]
            subprocess.run(index_above_cmd, check=True)
            logger.info("Indexing complete")
        if side_selection in ['below', 'both']:
            logger.info("Indexing {output_below}")
            index_below_cmd = ['samtools', 'index', file_below]
            subprocess.run(index_below_cmd, check=True)
            logger.info("Indexing complete")

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error processing BAM file: {e}", file=sys.stderr)
    finally:
        if output_above:
            output_above.close()
        if output_below:
            output_below.close()
