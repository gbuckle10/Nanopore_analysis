import argparse
import os
import subprocess
import sys
from pathlib import Path

from src.utils.logger import setup_logger
from utils.runner import load_config, get_project_root, run_external_command

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "src", "runtime_config.sh")


def align_and_index(config):
    unaligned_bam_dir = os.path.join(
        project_root,
        config['paths']['basecalled_output_dir']
    )
    unaligned_bam = os.path.join(
        unaligned_bam_dir,
        config['paths']['unaligned_bam_name']
    )
    threads = config['parameters']['general']['threads']
    reference_index = os.path.join(
        project_root,
        'reference_genomes',
        config['paths']['indexed_ref_gen_fasta_name']
    )
    sort_memory_limit = config['parameters']['general']['sort_memory_limit']
    aligned_bam_dir = os.path.join(
        project_root,
        config['paths']['alignment_output_dir']
    )
    aligned_bam_file = os.path.join(
        aligned_bam_dir,
        config['paths']['aligned_bam_name']
    )

    # Add a check - allow the user to decide whether to align a single bam or a directory of bams.
    # If the user wants to a directory, you need to specify an --output-dir
    # Sorting and indexing is automatic if a directory is given instead of a specific file.
    alignment_cmd = [
        'dorado', 'aligner',
        '-t', threads,
        reference_index,
        unaligned_bam
    ]

    # The sort shouldn't happen if the alignment_cmd is for an entire directory.
    sort_cmd = [
        "samtools", "sort",
        "-@", threads,
        "-m", sort_memory_limit,
        "-o", aligned_bam_file
    ]

    print(f"Executing pipe: {' '.join(alignment_cmd)} | {' '.join(sort_cmd)}")

    try:
        # Start aligner process
        # stdout will be captured for the pipe
        align_process = subprocess.Popen(alignment_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Start sort process
        # stdin is connected to stdout of the first process.
        sort_process = subprocess.Popen(sort_cmd, stdin=align_process.stdout, stderr=subprocess.PIPE)

        # Close the pipe in the alignment process so that it can receive errors from sort process.
        align_process.stdout.close()

        # Wait for the processes to finish and capture their stderr output.
        align_stderr = align_process.communicate()[1]
        sort_stderr = sort_process.communicate()[1]

        if align_process.returncode != 0:
            print("---ERROR in alignment step ---", file=sys.stderr)
            print(align_stderr.decode(), file=sys.stderr)
            sys.exit(1)
        if sort_process.returncode != 0:
            print("--- ERROR in sorting step ---", file=sys.stderr)
            print(sort_stderr.decode(), file=sys.stderr)
            sys.exit(1)
        print(f"--- Alignment and sorting complete. Output sent to {aligned_bam_file}")

    except FileNotFoundError as e:
        print(f"CRITICAL: Command '{e.filename}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    print("Now it's time to index.")
    index_cmd = [
        "samtools", "index",
        aligned_bam_file
    ]

    run_external_command(index_cmd)

    print("Indexing complete")

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    setup_logger()

    parser = argparse.ArgumentParser(
        description="Align a basecalled bam file to a reference genome."
    )

    parser.add_argument(
        '-c', '--config',
        type=Path,
        help="Path to the config file."
    )
    parser.add_argument(
        "--index-genome",
        action='store_true',
        help='Index a genome with samtools'
    )

    args = parser.parse_args()

    if not args.config:
        config = load_config(CONFIG_PATH)
    else:
        config = load_config(args.config)

    method_flagged = args.index_genome

    align_and_index(config)

if __name__ == '__main__':
    main()