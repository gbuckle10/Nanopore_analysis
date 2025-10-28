import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

from src.utils.cli_utils import create_io_parser
from src.utils.process_utils import run_command
from src.utils.config_utils import resolve_param
from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)

def full_alignment_handler(args,config):

    unaligned_bam = resolve_param(
        args,
        config,
        arg_name='input_file',
        construct_path=True,
        config_path=[
            ['paths', 'basecalled_output_dir'],
            ['paths', 'unaligned_bam_name']
        ]
    )

    threads = resolve_param(
        args, config, arg_name='threads', config_path=['parameters', 'general', 'threads']
    )

    reference_index = resolve_param(
        args, config, arg_name='ref', config_path=['paths', 'indexed_ref_gen_fasta_name']
    )

    sort_memory_limit = resolve_param(
        args, config, config_path=['parameters', 'general', 'sort_memory_limit']
    )

    # Intermediate aligned sorted bam
    aligned_bam_file = resolve_param(
        args, config, construct_path=True,
        config_path=[
            ['paths', 'alignment_output_dir'],
            ['paths', 'aligned_bam_name']
        ]
    )

    dorado_exe = resolve_param(
        args, config, config_path=['tools', 'dorado']
    )

    run_alignment_command(dorado_exe, unaligned_bam, aligned_bam_file, reference_index, sort_memory_limit, threads)

    flagstat_report = resolve_param(
        args, config, arg_name='output_dir', construct_path=True,
        config_path=[
            ['paths', 'flagstat_report'],
            ['paths', 'alignment_flagstat_name']
        ]
    )

    run_qc_command(aligned_bam_file, flagstat_report)


def alignment_handler(args, config):
    unaligned_bam = resolve_param(
        args,
        config,
        arg_name='input_file',
        construct_path=True,
        config_path=[
            ['paths', 'basecalled_output_dir'],
            ['paths', 'unaligned_bam_name']
        ]
    )

    aligned_bam_file = resolve_param(
        args, config, arg_name='output_dir', construct_path=True,
        config_path=[
            ['paths', 'alignment_output_dir'],
            ['paths', 'aligned_bam_name']
        ]
    )

    threads = resolve_param(
        args, config, arg_name='threads', config_path=['parameters', 'general', 'threads']
    )

    reference_index = resolve_param(
        args, config, arg_name='ref', config_path=['paths', 'indexed_ref_gen_fasta_name']
    )

    sort_memory_limit = resolve_param(
        args, config, config_path=['parameters', 'general', 'sort_memory_limit']
    )

    dorado_exe = resolve_param(
        args, config, config_path=['tools', 'dorado']
    )

    run_alignment_command(dorado_exe, unaligned_bam, aligned_bam_file, reference_index, sort_memory_limit, threads)

def qc_handler(args, config):
    aligned_sorted_file = resolve_param(
        args, config, arg_name="input_file", construct_path=True,
        config_path=[
            ['paths', 'alignment_output_dir'],
            ['paths', 'aligned_bam_name']
        ]
    )

    flagstat_report = resolve_param(
        args, config, arg_name='output_dir', construct_path=True,
        config_path=[
            ['paths', 'flagstat_report'],
            ['paths', 'alignment_flagstat_name']
        ]
    )

    run_qc_command(aligned_sorted_file, flagstat_report)


def run_alignment_command(dorado_exe, unaligned_bam, aligned_bam_file, reference_index, sort_memory_limit, threads):
    # Add a check - allow the user to decide whether to align a single bam or a directory of bams.
    # If the user wants to a directory, you need to specify an --output-dir
    # Sorting and indexing is automatic if a directory is given instead of a specific file.

    # ADD A METHOD TO PRODUCE MMI REFERENCE HERE IF IT DOESN'T ALREADY EXIST

    alignment_cmd = [
        'aligner',
        '-t', str(threads),
        str(reference_index),
        str(unaligned_bam)
    ]

    dorado_runner = ToolRunner(dorado_exe)

    # The sort shouldn't happen if the alignment_cmd is for an entire directory.
    sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-m", sort_memory_limit,
        "-o", str(aligned_bam_file)
    ]

    print(f"Executing pipe: {' '.join(alignment_cmd)} | {' '.join(sort_cmd)}")

    align_process = dorado_runner.start(
        alignment_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    if align_process is None:
        return

    try:
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

    run_command(index_cmd)

    print("Indexing complete")


def run_qc_command(aligned_sorted_file, flagstat_report):


    flagstat_cmd = [
        "samtools", "flagstat",
        aligned_sorted_file
    ]

    print(f"Running command: {' '.join(flagstat_cmd)}")

    try:
        result = subprocess.run(
            flagstat_cmd,
            check=True,
            capture_output=True,  # Capture stdout/stderr
            text=True  # Decodes output as strings
        )

        with open(flagstat_report, 'w') as f:
            f.write(result.stdout)

        print(f"Flagstat report saved to {flagstat_report}")

    except FileNotFoundError:
        print(f"CRITICAL: 'samtools' not found.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"CRITICAL: samtools flagstat failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(1)


def setup_parsers(subparsers, parent_parser):
    # Make new parents
    io_parser = create_io_parser()
    alignment_parent_parser = argparse.ArgumentParser(add_help=False)
    alignment_parent_parser.add_argument(
        "--ref",
        type=Path,
        help="Path to the reference genome"
    )

    alignment_parser = subparsers.add_parser(
        "align",
        help="Align basecalled BAM file to a specified genome, and QC the alignment. Use subcommands to run "
             "individual steps.",
        description="This command group contains tools for aligning reads to a reference and outputting quality "
                    "analysis on the resulting alignments.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['alignment'],
        epilog="""
Example Usage:
    nanopore_analysis align run --ref genome.fa --input-file reads.bam --output-dir analysis
    nanopore_analysis align qc --input-file aligned_reads.bam
"""
    )

    def show_align_help(args, config):
        """Default function to show help for the align command group"""
        alignment_parser.print_help()

    alignment_parser.set_defaults(func=show_align_help)
    alignment_subparsers = alignment_parser.add_subparsers(
        title="Available Commands",
        description="Choose one of the following actions to perform.",
        dest='subcommand',
        metavar="<command>"
    )

    p_run = alignment_subparsers.add_parser(
        'run', help="Align a BAM file to a specified genome and QC the alignment.",
        parents=[parent_parser, io_parser, alignment_parent_parser]
    )
    p_run.add_argument(
        '--threads', type=int, help="Number of threads for alignment and samtools."
    )
    p_run.set_defaults(func=full_alignment_handler)

    p_qc_only = alignment_subparsers.add_parser(
        'qc', help="Generate alignment statistics for a provided BAM file.",
        parents=[parent_parser, io_parser, alignment_parent_parser]
    )
    p_qc_only.set_defaults(func=qc_handler)

    p_align_only = alignment_subparsers.add_parser(
        'align', help="Align a BAM file to a specified genome and index.",
        parents=[parent_parser, io_parser, alignment_parent_parser]
    )
    p_align_only.set_defaults(func=alignment_handler)