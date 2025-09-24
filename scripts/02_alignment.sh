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
REFERENCE_GENOME_URL=$(yq e '.paths.reference_genome_url' "${CONFIG_FILE}")
REFERENCE_GENOME_NAME=$(yq e '.paths.reference_genome_name' "${CONFIG_FILE}")
REFERENCE_GENOME_INDEX=$(yq e '.paths.indexed_ref_gen_name' "${CONFIG_FILE}")
BASECALLED_OUTPUT_DIR=$(yq e '.paths.basecalled_output_dir' "${CONFIG_FILE}")
ALIGNED_OUTPUT_DIR=$(yq e '.paths.alignment_output_dir' "${CONFIG_FILE}")
UNALIGNED_BAM_NAME=$(yq e '.paths.unaligned_bam_name' "${CONFIG_FILE}")
ALIGNED_BAM_NAME=$(yq e '.paths.aligned_bam_name' "${CONFIG_FILE}")
THREADS=$(yq e '.parameters.general.threads' "${CONFIG_FILE}")
REFERENCE_GENOME_DIR='reference_genomes/'
SORT_MEMORY_LIMIT=$(yq e '.parameters.general.sort_memory_limit' "${CONFIG_FILE}")

check_vars \
  "REFERENCE_GENOME_URL" \
  "REFERENCE_GENOME_NAME" \
  "REFERENCE_GENOME_INDEX" \
  "BASECALLED_OUTPUT_DIR" \
  "ALIGNED_OUTPUT_DIR" \
  "UNALIGNED_BAM_NAME" \
  "ALIGNED_BAM_NAME" \
  "THREADS" \
  "REFERENCE_GENOME_DIR" \
  "SORT_MEMORY_LIMIT"

if [ ! -f "scripts/runtime_config.sh" ]; then
  echo "ERROR: runtime_config.sh not found. Please run 'setup' step first..."
  exit 1
fi
source "scripts/runtime_config.sh"

# --- Main script ---
log_info "Using Dorado executable at ${DORADO_EXECUTABLE}"


download_reference_genome() {

  log_info "Checking for reference genome"

  REF_FASTA="${REFERENCE_GENOME_DIR}${REFERENCE_GENOME_NAME}"


  # Check to see whether the file already exists.
  if [ -f "${REF_FASTA}" ]; then
    log_info "Reference file ${REF_FASTA} already exists, so we will skip the download"
    return 0
  fi

  log_info "--- Reference genome not found at ${REF_FASTA}. Downloading from ${REFERENCE_GENOME_URL} to ${REFERENCE_GENOME_DIR} ---"

  aws s3 cp "${REFERENCE_GENOME_URL}" "${REF_FASTA}" --no-sign-request

  log_info "Genome has been downloaded and indexed."

  # The output of this should be in the format >chr1 >chr2 etc.
  # grep ">" "${REF_FASTA}" | head -n 5
}

index_reference_genome() {
  REF_FASTA="${REFERENCE_GENOME_DIR}${REFERENCE_GENOME_NAME}"
  REF_MMI="${REFERENCE_GENOME_DIR}${REFERENCE_GENOME_INDEX}"

  log_info "Checking for reference index at ${REF_MMI}"
  if [ ! -f "${REF_MMI}" ]; then
    log_info "Index not found, creating minimap2 index"
    minimap2 -d "${REF_MMI}" "${REF_FASTA}"
    log_info "Genome index created at ${REF_MMI}"
  else
    log_info "Index already exists, so we'll skip indexing."
  fi
}

align_and_index() {
  # In the end it'd be good to have the option to align and basecall at once, or optionally do the two separately.
  local UNALIGNED_BAM="${BASECALLED_OUTPUT_DIR}/${UNALIGNED_BAM_NAME}"
  local REFERENCE_INDEX="reference_genomes/${REFERENCE_GENOME_INDEX}"
  local ALIGNED_BAM="${ALIGNED_OUTPUT_DIR}/${ALIGNED_BAM_NAME}"


  log_info "--- Starting alignment ---"
  mkdir -p "data/alignment_output"

  index_reference_genome

  local alignment_cmd=(
    "${DORADO_EXECUTABLE}"
    aligner
    -t "${THREADS}"
    "${REFERENCE_INDEX}"
    "${UNALIGNED_BAM}"
  )

  local sort_cmd=(
    "samtools"
    "sort"
    -@ "${THREADS}"
    -m "${SORT_MEMORY_LIMIT}"
    -o "${ALIGNED_BAM}"
  )

  log_info "Executing command: ${alignment_cmd[*]} | ${sort_cmd[*]}" # [*] expands the array into a single string for printing

  # Execute the pipeline
  "${alignment_cmd[@]}" | "${sort_cmd[@]}" # [@] expands the array into a list of individually quoted words for running the command

  log_info "Alignment and sorting complete"

  # Define and run the indexing command
  local index_cmd=(
    "samtools"
    "index"
    "${ALIGNED_BAM}"
  )

  log_info "Executing indexing command: ${index_cmd[*]}"

  "${index_cmd[@]}"

  log_info "Indexing complete."

  # To make sure the alignment worked well, you can use something like
  # samtools view data/alignment_output/aligned.sorted.bam | grep 'MM:Z:' | head -n 5

}


download_reference_genome
align_and_index