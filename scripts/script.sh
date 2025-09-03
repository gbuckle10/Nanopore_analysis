#!/bin/bash

download_reference_genome() {
  local url=$1

  REF_DIR="reference_genomes"
  REF_FASTA="${REF_DIR}/hg38.fa"

  echo "--- Downloading reference genome ---"
  echo "URL: ${url}"
  echo "Final destination: ${REF_DIR}"

  # Create the directory
  mkdir -p "${REF_DIR}"

  # Check to see whether the file already exists.
  if [ -f "${REF_FASTA}" ]; then
    echo "Reference file ${REF_FASTA} already exists, so we will skip the download"
    return 0
  fi

  echo "Reference genome not found, so let's download it."
  echo "We'll download from ${url}"

  aws s3 cp "${url}" "${REF_FASTA}" --no-sign-request

  echo "Genome has been downloaded."

  # The output of this should be in the format >chr1 >chr2 etc.
  grep ">" "${REF_FASTA}" | head -n 5
}

alignment_qc() {
  local ALIGNED_BAM="data/alignment_output/aligned.sorted.bam"
  local QC_DIR="data/qc_output/"

  echo "--- Starting alignment quality control ---"

  # Make QC output directory
  mkdir -p "${QC_DIR}"

  # Summary with flagstat
  local flagstat_report="${QC_DIR}/alignment.flagstat.txt"
  echo " Generating summary with flagstat"
  samtools flagstat "${ALIGNED_BAM}" > "${flagstat_report}"
  echo " Flagstat report saved to ${flagstat_report}"

  # Detailed report with stats
  local stats_report="${QC_DIR}/alignment.stats.txt"
  echo " Generating detailed report"
  samtools stats "${ALIGNED_BAM}" > "${stats_report}"
  echo " Stats report saved to ${stats_report}"

  #You can print the top 5 lines of the aligned.sorted.bam file using samtools
  # samtools view aligned.bam | head -n 5

}

index_reference_genome() {
  REF_FASTA="reference_genomes/hg38.fa"
  REF_MMI="reference_genomes/hg38.mmi"

  echo "Checking for reference index"
  if [ ! -f "${REF_MMI}" ]; then
    echo "INFO: Index not found, creating minimap2 index"
    dorado index "${REF_MMI}" "${REF_FASTA}"
    echo "We did it! index created"
  else
    echo "INFO: Index already exists, so we'll skip indexing."
  fi

}

align_and_index() {
  # In the end it'd be good to have the option to align and basecall at once, or optionally do the two separately.
  local UNALIGNED_BAM="data/basecalled_output/taken_bam.bam"
  local REFERENCE_INDEX="reference_genomes/$1"
  local ALIGNED_BAM="data/alignment_output/aligned.sorted.bam"

  local THREADS=1
  local SORT_MEMORY_LIMIT="1G"

  echo "--- Starting alignment ---"
  mkdir -p "data/alignment_output"

  # To save on RAM during the alignment, it's much better to index the reference genome before doing the alignment.
  index_reference_genome


  dorado aligner -t 1 "${REFERENCE_INDEX}" "${UNALIGNED_BAM}" \
  | samtools sort -@ ${THREADS} -m ${SORT_MEMORY_LIMIT} -o "${ALIGNED_BAM}"

  samtools index "${ALIGNED_BAM}"

  echo "INFO: Alignment and indexing finished successfully."

  # To make sure the alignment worked well, you can use something like
  # samtools view data/alignment_output/aligned.sorted.bam | grep 'Mm:Z:' | head -n 5


}

run_script(){

  DORADO_VERSION="1.0.0"
  install_dorado "$DORADO_VERSION"
  #download_fast5_data 1
  #convert_fast5_to_pod5
  #basecalling_pod5 64
  download_reference_genome s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
  align_and_index "hg38.mmi"

  alignment_qc

  ## MAKE A GENERIC DOWNLOADER ONE DAY
}

run_script