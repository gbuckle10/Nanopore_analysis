#!/bin/bash

set -e

CONFIG_FILE=$1
echo ${CONFIG_FILE}
yq eval "${CONFIG_FILE}"

ALIGNED_OUTPUT_DIR=$(yq e '.paths.alignment_output_dir' "${CONFIG_FILE}")
ALIGNED_BAM_NAME=$(yq e '.paths.aligned_bam_name' "${CONFIG_FILE}")
QC_DIR=$(yq e '.paths.qc_dir' "${CONFIG_FILE}")
FLAGSTAT_NAME=$(yq e '.paths.alignment_flagstat_name' "${CONFIG_FILE}")
ALIGNMENT_STATS_NAME=$(yq e '.paths.alignment_stats_name' "${CONFIG_FILE}")

alignment_qc() {
  local ALIGNED_BAM="${ALIGNED_OUTPUT_DIR}${ALIGNED_BAM_NAME}"


  echo "--- Starting alignment quality control ---"

  # Make QC output directory
  mkdir -p "${QC_DIR}"

  # Summary with flagstat
  local flagstat_report="${QC_DIR}${FLAGSTAT_NAME}"
  echo " Generating summary with flagstat"
  samtools flagstat "${ALIGNED_BAM}" > "${flagstat_report}"
  echo " Flagstat report saved to ${flagstat_report}"

  # Detailed report with stats
  local stats_report="${QC_DIR}${ALIGNMENT_STATS_NAME}"
  echo " Generating detailed report"
  samtools stats "${ALIGNED_BAM}" > "${stats_report}"
  echo " Stats report saved to ${stats_report}"

  #You can print the top 5 lines of the aligned.sorted.bam file using samtools
  # samtools view aligned.bam | head -n 5

}
