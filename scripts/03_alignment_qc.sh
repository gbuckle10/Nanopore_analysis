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

ALIGNED_OUTPUT_DIR=$(yq e '.paths.alignment_output_dir' "${CONFIG_FILE}")
ALIGNED_BAM_NAME=$(yq e '.paths.aligned_bam_name' "${CONFIG_FILE}")
QC_DIR=$(yq e '.paths.qc_dir' "${CONFIG_FILE}")
FLAGSTAT_NAME=$(yq e '.paths.alignment_flagstat_name' "${CONFIG_FILE}")
ALIGNMENT_STATS_NAME=$(yq e '.paths.alignment_stats_name' "${CONFIG_FILE}")

alignment_qc() {
  local ALIGNED_BAM="${ALIGNED_OUTPUT_DIR}${ALIGNED_BAM_NAME}"

  # Make QC output directory
  mkdir -p "${QC_DIR}"

  # --- Summary with flagstat ---
  local flagstat_report="${QC_DIR}${FLAGSTAT_NAME}"
  log_info "Alignment summary will be saved to ${flagstat_report}"

  # Define command
  local flagstat_cmd=(
    "samtools"
    "flagstat"
    "${ALIGNED_BAM}"
  )

  log_info "Generating summary with flagstat, command: ${flagstat_cmd[*]} > ${flagstat_report}"

  # Run command
  "${flagstat_cmd[@]}" > "${flagstat_report}"

  log_info " Flagstat report saved to ${flagstat_report}"

  # --- Detailed report with samtools ---
  local stats_report="${QC_DIR}${ALIGNMENT_STATS_NAME}"
  log_info "Detailed report will be saved to ${stats_report}"

  local stats_cmd=(
    "samtools"
    "stats"
    "${ALIGNED_BAM}"
  )

  log_info "Generating alignment report with samtools, command: ${stats_cmd[*]} > ${stats_report}"

  "${stats_cmd[@]}" > "${stats_report}"

  log_info "Stats report saved to ${stats_report}"

  #You can print the top 5 lines of the aligned.sorted.bam file using samtools
  # samtools view aligned.bam | head -n 5
}

alignment_qc