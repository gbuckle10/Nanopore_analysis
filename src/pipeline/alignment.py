import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import List

from src.config.models import AppSettings
from src.utils.cli_utils import add_io_arguments
from src.utils.process_utils import run_command
from src.utils.tools_runner import ToolRunner

logger = logging.getLogger(__name__)


def _get_file_basename(reference_path: Path):
    """
    Get the basename of a genome file, handling one and two extensions. e.g. genome.fa and genome.fa.gz
    """
    fasta_suffixes = ['.fasta', '.fa', '.fna']
    compression_suffixes = ['.gz', '.bz2', '.xz']

    name = reference_path.name

    for suffix in compression_suffixes:
        if name.endswith(suffix):
            name = name.removesuffix(suffix) # Might slice for compatibility with python <3.9
            break
    for suffix in fasta_suffixes:
        if name.endswith(suffix):
            name = name.removesuffix(suffix)
            break

    return name

def _ensure_mmi_exists(ref_fasta_path: Path, threads: int) -> Path:
    """
    Given a path to a FASTA file, ensures the corresponding .mmi index exists, and creates it if necessary. It can
    handle .fa and .fa.gz files.

    """
    basename = _get_file_basename(ref_fasta_path)
    index_path = ref_fasta_path.parent / f"{basename}.mmi"

    if index_path.is_file():
        logging.info(f"Found existing minimap2 index: {index_path}")
        return index_path

    # If the index is missing, we need to find the source FASTA.
    source_fasta_to_index = None

    if ref_fasta_path.is_file():
        # If the provided fasta exists, just index that.
        source_fasta_to_index = ref_fasta_path
    else:
        # If the provided fasta doesn't exist, check for the gzipped version (.fa.gz)
        gz_fasta_path = ref_fasta_path.with_suffix('.fa.gz')
        if gz_fasta_path.is_file():
            source_fasta_to_index = gz_fasta_path

    # Give an error if we can't find any source fasta
    if source_fasta_to_index is None:
        raise FileNotFoundError(f"Couldn't find minimap2 index ({index_path}) or a source fasta file to index at "
                                f"'{ref_fasta_path}' or {ref_fasta_path.with_suffix('.fa.gz')}")

    logging.info(f"Minimap2 index not found. Building index from {source_fasta_to_index}")

    index_cmd = [
        "minimap2",
        "-d", str(index_path),
        "-t", str(threads),
        str(source_fasta_to_index)
    ]

    run_command(index_cmd)

    if not index_path.is_file():
        raise RuntimeError(f"Minimap2 indexing failed. Index was not created at {index_path}")

    logging.info(f"Successfully built index: {index_path}")
    return index_path

def _resolve_alignment_inputs(input_path: Path) -> List[Path]:
    """
    Finds all BAM files to be aligned from a given input path
    """

    logging.info(f"Resolving alignment input {input_path}")
    if not input_path.exists():
        raise FileNotFoundError(f"Input for alignment doesn't exist: {input_path}")

    if input_path.is_file():
        logging.info(f"The input {input_path} is a file, and will be aligned.")
        if input_path.suffix != '.bam':
            logging.warning(f"Input file {input_path} doesn't have a .bam extension")
        return [input_path]
    if input_path.is_dir():
        logging.info(f"The input {input_path} is a directory, so we'll look for .bam files inside it.")
        found_files = list(input_path.rglob('*.bam'))
        if not found_files:
            raise FileNotFoundError(f"No BAM files found in directory {input_path}")
        return found_files
    raise ValueError(f"Input path is neither a file nor a directory: {input_path}")

def _prepare_alignment_io(input_files: List[Path], output_path: Path):

    is_single_file_mapping = (len(input_files) == 1 and output_path.suffix)

    if is_single_file_mapping:
        # A single input and output file was specified
        input_file = input_files[0]
        output_file = output_path
        output_dir = output_file.parent
        logging.info(f"Single alignment: Aligning '{input_file.name}' directly to '{output_file}'.")

        output_dir.mkdir(parents=True, exist_ok=True)

        return output_dir, [(input_file, output_file)]

    output_dir: Path
    if output_path.suffix:
        # The specified output is a file, but the input contained multiple files
        output_dir = output_path.parent
        logging.warning(f"Input is a directory, so the output path '{output_path}' will be treated as a directory. "
                        f"Using its parent: '{output_dir}")
    else:
        # The specified output is a directory, so we'll just save the file into the directory
        output_dir = output_path
        logging.info(f"The alignment output will be saved in {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"We'll be aligning a batch of files: All alignments will be saved in: {output_dir}.")

    io_pairs = []
    for input_file in input_files:
        # Get the input file stem - can't just use stem in case there are more than one extension
        stem = input_file.name.replace(".bam", "").replace(".unaligned", "")
        output_bam = output_dir / f"{stem}.aligned.bam"
        io_pairs.append((input_file, output_bam))

    return output_dir, io_pairs

def full_alignment_handler(config: AppSettings):
    unaligned_bam = config.pipeline_steps.align.paths.full_unaligned_input_path
    aligned_bam_file = config.pipeline_steps.align.paths.full_aligned_bam_path
    threads = config.globals.threads
    reference_fasta_path = config.pipeline_steps.align.paths.full_ref_fasta_path
    sort_memory_limit = config.globals.sort_memory_limit
    dorado_exe = config.tools.dorado

    if reference_fasta_path is None:
        raise ValueError("Missing required argument: --ref")

    try:
        logging.info(f"Checking whether reference index exists at {reference_fasta_path.parent}")
        reference_index_path = _ensure_mmi_exists(reference_fasta_path, threads)
    except (FileNotFoundError, RuntimeError) as e:
        logging.info(f"Error preparing reference genome: {e}", file=sys.stderr)
        sys.exit(1)

    run_alignment_command(dorado_exe, unaligned_bam, aligned_bam_file, reference_index_path, sort_memory_limit, threads)

    flagstat_report = config.pipeline_steps.align.paths.full_flagstat_path

    run_qc_command(aligned_bam_file, flagstat_report)


def alignment_handler(config):
    unaligned_bam = config.pipeline_steps.align.paths.full_unaligned_input_path
    aligned_bam_file = config.pipeline_steps.align.paths.full_aligned_bam_path
    threads = config.globals.threads
    reference_fasta_path = config.pipeline_steps.align.paths.full_ref_fasta_path

    sort_memory_limit = config.globals.sort_memory_limit
    dorado_exe = config.tools.dorado

    if reference_fasta_path is None:
        raise ValueError("Missing required argument: --ref")

    try:
        logging.info(f"Checking whether reference index exists at {reference_fasta_path.parent}")
        reference_index_path = _ensure_mmi_exists(reference_fasta_path, threads)
    except (FileNotFoundError, RuntimeError) as e:
        logging.info(f"Error preparing reference genome: {e}", file=sys.stderr)
        sys.exit(1)

    run_alignment_command(dorado_exe, unaligned_bam, aligned_bam_file, reference_index_path, sort_memory_limit, threads)


def qc_handler(config):
    aligned_bam_file = config.pipeline_steps.align.paths.full_aligned_bam_path
    flagstat_report = config.pipeline_steps.align.paths.full_flagstat_path

    run_qc_command(aligned_bam_file, flagstat_report)


def run_alignment_command(dorado_exe, input_path, output_path, reference_index, sort_memory_limit, threads):
    # Add a check - allow the user to decide whether to align a single bam or a directory of bams.
    # If the user wants to a directory, you need to specify an --output-dir
    # Sorting and indexing is automatic if a directory is given instead of a specific file.

    try:
        input_files_to_align = _resolve_alignment_inputs(input_path)

        # Validate output path
        output_dir, io_pairs = _prepare_alignment_io(input_files_to_align, output_path)

        logging.info(f"The output dir is {output_dir}")
        for pair in io_pairs:
            logging.info(f"    {pair[0]} -> {pair[1]}")
    except (FileNotFoundError, ValueError) as e:
        logging.error(f"Error preparing for alignment: {e}")

    for input_file, output_bam_file in io_pairs:
        logging.info(f"Processing input file: {input_file.name} -> {output_bam_file.name}")

        alignment_cmd = [
            'aligner',
            '-t', str(threads),
            str(reference_index),
            str(input_file)
        ]

        sort_cmd = [
            "samtools", "sort",
            "-@", str(threads),
            "-m", str(sort_memory_limit),
            "-o", str(output_bam_file)
        ]

        try:
            dorado_runner = ToolRunner(dorado_exe)
            logger.info(f"Executing pipe: dorado {' '.join(alignment_cmd)} | {' '.join(sort_cmd)}")
            align_process = dorado_runner.start(
                alignment_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            if align_process is None:
                return
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
            print(f"--- Alignment and sorting complete. Output sent to {output_bam_file}")

        except FileNotFoundError as e:
            print(f"CRITICAL: Command '{e.filename}' not found.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}", file=sys.stderr)
            sys.exit(1)

        print("Now it's time to index.")
        index_cmd = [
            "samtools", "index",
            str(output_bam_file)
        ]

        run_command(index_cmd)

    print("Indexing complete")


def run_qc_command(aligned_sorted_file, flagstat_report):
    flagstat_cmd = [
        "samtools", "flagstat",
        str(aligned_sorted_file)
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


def _make_alignment_parent_parser(config):
    """
    Creates a parent parser with the shared alignment arguments.
    :param config:
    :return:
    """
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        "--ref",
        default=None,
        dest="pipeline_steps.align.paths.custom_fasta_reference",
        type=Path,
        help="Path to the reference genome"
    )
    parent.add_argument(
        "--threads",
        default=None,
        dest="globals.threads",
        help="Number of threads for alignment and samtools."
    )
    return parent


def add_all_arguments_to_parser(parser, config):
    """
    Publically available function to add alignment arguments to a given parser
    :param parser:
    :param config:
    :return:
    """
    temp_parent = _make_alignment_parent_parser(config)
    for action in temp_parent._actions:
        parser._add_action(action)

    # If there are any individual arguments, add those underneath.


def setup_parsers(subparsers, parent_parser, config):
    # Make new parents
    alignment_parent_parser = _make_alignment_parent_parser(config)

    alignment_parser = subparsers.add_parser(
        "align",
        help="Align basecalled BAM file to a specified genome, and QC the alignment. Use subcommands to run "
             "individual steps.",
        description="This command group contains tools for aligning reads to a reference and outputting quality "
                    "analysis on the resulting alignments.",
        formatter_class=argparse.RawTextHelpFormatter,
        aliases=['alignment'],
        parents=[parent_parser],
        epilog="""
Example Usage:
    nanopore_analysis align run --ref genome.fa --input-file reads.bam --output-dir analysis
    nanopore_analysis align qc --input-file aligned_reads.bam
"""
    )
    validation_func = lambda: config.pipeline_steps.align.paths._validate()
    alignment_parser.set_defaults(validation_func=validation_func)

    def show_align_help(config):
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
        parents=[parent_parser, alignment_parent_parser]
    )

    add_io_arguments(
        p_run, config,
        default_input=None,
        input_file_help="Path to full unaligned BAM file",
        input_dest="pipeline_steps.align.paths.user_alignment_input",
        default_output=None,
        output_dir_help="Path to aligned, sorted and indexed BAM file",
        output_dest="pipeline_steps.align.paths.user_alignment_output"
    )
    p_run.set_defaults(func=full_alignment_handler)

    p_qc_only = alignment_subparsers.add_parser(
        'qc', help="Generate alignment statistics for a provided BAM file.",
        parents=[parent_parser, alignment_parent_parser]
    )
    add_io_arguments(
        p_qc_only, config,
        default_input=None,
        input_file_help="Path to aligned, sorted and indexed BAM file",
        input_dest="pipeline_steps.align.paths.user.alignment_output",
        default_output=None,
        output_dir_help="Filepath of saved output",
        output_dest="pipeline_steps.align.paths.alignment_flagstat_name"
    )
    p_qc_only.set_defaults(func=qc_handler)

    p_align_only = alignment_subparsers.add_parser(
        'align', help="Align a BAM file to a specified genome and index.",
        parents=[parent_parser, alignment_parent_parser]
    )
    add_io_arguments(
        p_align_only, config,
        default_input=None,
        input_file_help="Path to full unaligned BAM file",
        input_dest="pipeline_steps.align.paths.user_alignment_input",
        default_output=None,
        output_dir_help="Path to aligned, sorted and indexed BAM file",
        output_dest="pipeline_steps.align.paths.user_alignment_output"
    )

    p_align_only.set_defaults(
        func=alignment_handler
    )
