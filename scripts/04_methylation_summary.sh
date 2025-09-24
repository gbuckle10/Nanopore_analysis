#!/bin/bash
# Exit on error
set -e

# Check to make sure that the runtime_config is there.
if [ ! -f "scripts/runtime_config.sh" ]; then
  log_error "ERROR: runtime_config.sh not found. Please run 'setup' step first..."
  exit 1
fi

# Source the logging utility and runtime_config
# Assumes that logging.sh is in the utils folder of the current folder.
source "$(dirname "$0")/utils/logging.sh"
source "scripts/runtime_config.sh"

CONFIG_FILE=$1
ALIGNED_BAM_NAME=$(yq e '.paths.aligned_bam_name' "${CONFIG_FILE}")
ALIGNED_BAM_DIR=$(yq e '.paths.alignment_output_dir' "${CONFIG_FILE}")
METHYLATION_BED=$(yq e '.paths.methylation_bed_name' "${CONFIG_FILE}")
REFERENCE_FASTA=$(yq e '.paths.reference_genome_name' "${CONFIG_FILE}")
THREADS=$(yq e '.parameters.general.threads' "${CONFIG_FILE}")
METHYLATION_DIR=$(yq e '.paths.methylation_dir' "${CONFIG_FILE}")
REFERENCE_GENOME_DIR="reference_genomes/"
LOG_FILE=$(yq e '.paths.methylation_log_file' "${CONFIG_FILE}")

check_vars \
  "ALIGNED_BAM_NAME" \
  "ALIGNED_BAM_DIR" \
  "METHYLATION_BED" \
  "REFERENCE_FASTA" \
  "THREADS" \
  "METHYLATION_DIR" \
  "REFERENCE_GENOME_DIR" \
  "LOG_FILE"

if [[ -z "${ALIGNED_BAM_NAME}" || -z "$METHYLATION_BED" || -z "$REFERENCE_FASTA" || -z "$THREADS" ]]; then
  log_error "Missing one or more required arguments." >&2
  log_error "Usage: $0 <aligned_bam> <output_bed> <ref_fasta> <threads> [log_file]" >&2
  exit 1
fi

mkdir -p "$(dirname "${METHYLATION_BED}")"

index_genome() {

  log_info "Looking for reference file at ${REFERENCE_GENOME_DIR}${REFERENCE_FASTA}.fai".
  if [ -f "${REFERENCE_GENOME_DIR}${REFERENCE_FASTA}.fai" ]; then
    log_info "Reference file ${REFERENCE_FASTA}.fai already exists, so we will skip the download"
    return 0
  fi
  log_info "Indexing genome to ${REFERENCE_GENOME_DIR}${REFERENCE_FASTA}.fai."


  local index_cmd=(
    "samtools"
    "faidx"
    "${REFERENCE_GENOME_DIR}${REFERENCE_FASTA}"
  )

  log_info "Executing indexing command: ${index_cmd[*]}"
  "${index_cmd[@]}"

  log_info "Index created at ${REFERENCE_GENOME_DIR}${REFERENCE_FASTA}.fai"
}

methylation_pileup() {

  PATH_TO_ALIGNED_BAM="${ALIGNED_BAM_DIR}${ALIGNED_BAM_NAME}"
  PATH_TO_OUTPUT_BED="${METHYLATION_DIR}${METHYLATION_BED}"


  log_info "Running methylation pileup on ${PATH_TO_ALIGNED_BAM}"
  log_info "Sending data to ${PATH_TO_OUTPUT_BED}"

  # Define command
  local pileup_cmd=(
    "modkit"
    "pileup"
    "${PATH_TO_ALIGNED_BAM}"
    "${PATH_TO_OUTPUT_BED}"
    --filter-threshold 0.667
    --cpg
    --force-allow-implicit
    --ref "$REFERENCE_FASTA"
    --log-filepath "${LOG_FILE}"
  )
  log_info "Generating bedMethyl file with command: ${pileup_cmd[*]}"

  # Run command
  "${pileup_cmd[@]}"

  log_info "Pileup completed"
}

index_genome
methylation_pileup