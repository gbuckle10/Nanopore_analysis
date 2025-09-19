#!/bin/bash

CONFIG_FILE=$1

ALIGNED_BAM_NAME=$(yq e '.paths.aligned_bam_name' "${CONFIG_FILE}")
ALIGNED_BAM_DIR=$(yq e '.paths.alignment_output_dir' "${CONFIG_FILE}")
METHYLATION_BED=$(yq e '.paths.methylation_bed_name' "${CONFIG_FILE}")
REFERENCE_FASTA=$(yq e '.paths.reference_genome_name' "${CONFIG_FILE}")
THREADS=$(yq e '.parameters.general.threads' "${CONFIG_FILE}")
METHYLATION_DIR=$(yq e '.paths.methylation_dir' "${CONFIG_FILE}")
LOG_FILE=$(yq e '.paths.methylation_log_file' "${CONFIG_FILE}")

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