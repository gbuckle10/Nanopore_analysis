#!/bin/bash

ALIGNED_OUTPUT_DIR=$1
ALIGNED_BAM_NAME=$2
QC_DIR=$3
FLAGSTAT_NAME=$4
ALIGNMENT_STATS_NAME=$5

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
