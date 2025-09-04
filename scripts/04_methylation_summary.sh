#!/bin/bash

ALIGNED_BAM_NAME=$1
ALIGNED_BAM_DIR=$2
METHYLATION_BED=$3
REFERENCE_FASTA=$4
THREADS=$5
METHYLATION_DIR=$6
LOG_FILE=$7
REFERENCE_FASTA=$8

if [[ -z "${ALIGNED_BAM_NAME}" || -z "$METHYLATION_BED" || -z "$REFERENCE_FASTA" || -z "$THREADS" ]]; then
  echo "ERROR: Missing one or more required arguments." >&2
  echo "Usage: $0 <aligned_bam> <output_bed> <ref_fasta> <threads> [log_file]" >&2
  exit 1
fi

echo "INFO: Starting modkit pileup to summarise methylation frequencies..." >&2

mkdir -p "$(dirname "${METHYLATION_BED}")"

index_genome() {
  if [ -f "${REFERENCE_FASTA}.fai" ]; then
    echo "Reference file ${REFERENCE_FASTA}.fai already exists, so we will skip the download"
    return 0
  fi
  echo "Indexing genome to an fai file."
  samtools faidx "$REFERENCE_FASTA"
}

methylation_pileup() {

  PATH_TO_ALIGNED_BAM="${ALIGNED_BAM_DIR}${ALIGNED_BAM_NAME}"
  PATH_TO_OUTPUT_BED="${METHYLATION_DIR}${METHYLATION_BED}"


  echo "Running methylation pileup on ${PATH_TO_ALIGNED_BAM}"
  echo "Sending data to ${PATH_TO_OUTPUT_BED}"

  # --filter-threshold 0.667 is the same as -a 0.333 -b 0.667 in modbam2bed
  modkit pileup "${PATH_TO_ALIGNED_BAM}" "${PATH_TO_OUTPUT_BED}" --filter-threshold 0.667 --cpg --force-allow-implicit --ref "$REFERENCE_FASTA" --log-filepath "${LOG_FILE}"

  echo "modkit done with methylation summary"
}

index_genome
methylation_pileup