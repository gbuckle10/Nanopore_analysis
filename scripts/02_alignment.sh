#!/bin/bash


CONFIG_FILE=$1

REFERENCE_GENOME_URL=$(yq e '.paths.reference_genome_url' "${CONFIG_FILE}")
REFERENCE_GENOME_NAME=$(yq e '.paths.reference_genome_name' "${CONFIG_FILE}")
REFERENCE_GENOME_INDEX=$(yq e '.paths.indexed_ref_gen_name' "${CONFIG_FILE}")
BASECALLED_OUTPUT_DIR=$(yq e '.paths.basecalled_output_dir' "${CONFIG_FILE}")
ALIGNED_OUTPUT_DIR=$(yq e '.paths.alignment_output_dir' "${CONFIG_FILE}")
UNALIGNED_BAM_NAME=$(yq e '.paths.unaligned_bam_name' "${CONFIG_FILE}")
ALIGNED_BAM_NAME=$(yq e '.paths.aligned_bam_name' "${CONFIG_FILE}")
THREADS=$(yq e '.parameters.general.threads' "${CONFIG_FILE}")
SORT_MEMORY_LIMIT=$(yq e '.parameters.general.sort_memory_limit' "${CONFIG_FILE}")

if [ ! -f "scripts/runtime_config.sh" ]; then
  echo "ERROR: runtime_config.sh not found. Please run 'setup' step first..."
  exit 1
fi
source "scripts/runtime_config.sh"

# --- Main script ---
echo "INFO: Using Dorado executable at ${DORADO_EXECUTABLE}" >&2


download_reference_genome() {

  echo "Checking for reference genome"

  REF_FASTA="reference_genomes/${REFERENCE_GENOME_NAME}"


  # Check to see whether the file already exists.
  if [ -f "${REF_FASTA}" ]; then
    echo "Reference file ${REF_FASTA} already exists, so we will skip the download"
    return 0
  fi

  echo "--- Reference genome not found. Downloading from ${REFERENCE_GENOME_URL} to ${REFERENCE_GENOME_DIR} ---"

  aws s3 cp "${REFERENCE_GENOME_URL}" "${REF_FASTA}" --no-sign-request

  # To save on RAM during the alignment, it's much better to index the reference genome before doing the alignment.


  echo "Genome has been downloaded and indexed."

  # The output of this should be in the format >chr1 >chr2 etc.
  # grep ">" "${REF_FASTA}" | head -n 5
}

index_reference_genome() {
  REF_FASTA="reference_genomes/${REFERENCE_GENOME_NAME}"
  REF_MMI="reference_genomes/${REFERENCE_GENOME_INDEX}"

  echo "Checking for reference index at ${REF_MMI}"
  if [ ! -f "${REF_MMI}" ]; then
    echo "INFO: Index not found, creating minimap2 index"
    minimap2 -d "${REF_MMI}" "${REF_FASTA}"
    echo "We did it! index created"
  else
    echo "INFO: Index already exists, so we'll skip indexing."
  fi
}

align_and_index() {
  # In the end it'd be good to have the option to align and basecall at once, or optionally do the two separately.
  local UNALIGNED_BAM="${BASECALLED_OUTPUT_DIR}/${UNALIGNED_BAM_NAME}"
  local REFERENCE_INDEX="reference_genomes/${REFERENCE_GENOME_INDEX}"
  local ALIGNED_BAM="${ALIGNED_OUTPUT_DIR}/${ALIGNED_BAM_NAME}"


  echo "--- Starting alignment ---"
  mkdir -p "data/alignment_output"

  index_reference_genome

  "${DORADO_EXECUTABLE}" aligner -t 1 "${REFERENCE_INDEX}" "${UNALIGNED_BAM}" \
  | samtools sort -@ ${THREADS} -m ${SORT_MEMORY_LIMIT} -o "${ALIGNED_BAM}"

  samtools index "${ALIGNED_BAM}"

  echo "INFO: Alignment and indexing finished successfully."

  # To make sure the alignment worked well, you can use something like
  # samtools view data/alignment_output/aligned.sorted.bam | grep 'MM:Z:' | head -n 5


}


download_reference_genome
align_and_index