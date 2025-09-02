#!/bin/bash

check_device_for_dorado() {

  if command -v nvidia-smi &> /dev/null; then
    echo "NVIDIA GPU found, so we'll use that for basecalling"

    # Tell dorado to tun with --device flag cuda:0
    export DORADO_DEVICE_FLAG="--device cuda:0"
    # Default batchsize is 128
    export DORADO_BATCHSIZE_FLAG="--batchsize 128"

  else
    echo "No NVIDIA GPU found, so we'll run with the CPU only"

    # Tell dorado to run with --device flag cpu
    export DORADO_DEVICE_FLAG="--device cpu"
    # Default batchsize is 128
    export DORADO_BATCHSIZE_FLAG=""
  fi

  echo "Dorado device flags set: ${DORADO_DEVICE_FLAG} ${DORADO_BATCHSIZE_FLAG}"

}

# --- Dorado installation ---

install_dorado() {
  local VERSION=$1
  echo "--- Downloading and setting up Dorado ---"

  echo "-- Detecting OS --"
  case "$(uname -s)" in
    Linux*)
      os_type="linux-x64"
      archive_filename="dorado.tar.gz"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VERSION-${os_type}.tar.gz"
      ;;
    CYGWIN*|MINGW*|MSYS*)
      os_type="win64"
      archive_filename="dorado.zip"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-$VERSION-${os_type}.zip"
      ;;
    *)
      echo "I don't have a clue what operating system you're trying to use..."
      return 1
      ;;
  esac
  echo "The OS has been found, it's $os_type"

  dorado_folder=tools/dorado-${VERSION}-${os_type}
  echo "--- Checking for existing Dorado installation in ${dorado_folder}---"
  if [ -d "${dorado_folder}" ]; then
    echo "Dorado folder already exists"
    LOCAL_DORADO_BIN="${dorado_folder}/bin"
    export PATH="$LOCAL_DORADO_BIN:$PATH"
    echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"
    return 0
  fi

  echo "Dorado folder ${dorado_folder} not found, so I'll install it now."

  echo "Creating installation directory in tools"
  mkdir -p "tools"

  # Include a check here for whether the archive_filename is there to save more time.
  echo "Downloading dorado version $VERSION for ${os_type} from ${download_url}"
  curl -o tools/${archive_filename}  "$download_url"

  echo "Download complete. Unzipping and setting up."

  case ${os_type} in
    "linux-x64")
      tar -xzvf "tools/${archive_filename}" -C "tools/"
      ;;
    "win64")
      unzip -q "tools/${archive_filename}" -d "tools/"
      ;;
  esac

  folder_to_remove=tools/${archive_filename}
  echo "Removing ${folder_to_remove}"

  rm "${folder_to_remove}"

  echo "-------- Dorado was successfully installed ---------------"

  echo "-------- Setting Dorado exe to path --------------"
  LOCAL_DORADO_BIN="${dorado_folder}/bin"
  export PATH="$LOCAL_DORADO_BIN:$PATH"
  echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"

}

download_fast5_data() {
  local NUM_FILES=$1
  alias mkdir='/usr/bin/mkdir'
  set -e

  echo "----Starting Data Download Script----"

  #DOWNLOAD_FILE_PATH="s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/"
  DOWNLOAD_FILE_PATH="s3://ont-open-data/gm24385_mod_2021.09/flowcells/20210511_1515_X2_FAQ32637_9b683def/fast5_pass/"
  DESTINATION_DIR="data/fast5_input/"


  echo "--- Making fast5 input and basecalled output folders ---"
  mkdir -p data
  mkdir -p data/fast5_input
  mkdir -p data/basecalled_output

  echo "--- Downloading ${NUM_FILES} fast5 files ---"
  #aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request | head -n 5 | awk '{print $4}' | xargs -I {} aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/{} data/fast5_input/ --no-sign-request
  DATA_DOWNLOAD_DIR="s3://ont-open-data/gm24385_mod_2021.09/flowcells/20210511_1515_X2_FAQ32637_9b683def/fast5_pass/"


  aws s3 ls "${DATA_DOWNLOAD_DIR}" --no-sign-request \
  | head -n "$NUM_FILES" \
  | awk '{print $4}' \
  | while read -r filename; do
      aws s3 cp "${DATA_DOWNLOAD_DIR}${filename}" data/fast5_input/"$filename" --no-sign-request || true
  done

}

convert_fast5_to_pod5() {
  echo "--- Converting fast5 to pod5 ---"

  POD5_DIR="data/pod5_output"
  mkdir -p "${POD5_DIR}"

  pod5 convert fast5 data/fast5_input --output "${POD5_DIR}/all_reads.pod5" --force-overwrite
}

download_dorado_model() {
  MODEL_ALIAS="hac"
  MODEL_NAME="dna_r9.4.1_e8_sup@v3.3"
  MODEL_DIR="models"

  mkdir -p "${MODEL_DIR}"

  dorado download --model "${MODEL_NAME}" --models-directory "${MODEL_DIR}"


}

basecalling_pod5() {
  BATCHSIZE=$1
  echo "Adding dorado to the path"
  LOCAL_DORADO_BIN="$(pwd)/tools/dorado_v0.9.6/bin/"
  export PATH="$LOCAL_DORADO_BIN:$PATH"
  local MODEL_ALIAS="hac"
  local MODEL_MODIFICATION="5mCG_5hmCG"
  local POD5_INPUT="data/pod5_output/all_reads.pod5"
  local BASECALLED_OUTPUT="data/basecalled_output/calls.bam"

  # ---- You need to add the download_dorado_model here
  #  download_dorado_model
  # ---- You could also make it so that the alignment and basecalling is done at once. That should be the next thing to test.

  which dorado
  dorado basecaller --modified-bases 5mCG_5hmCG --batchsize "${BATCHSIZE}" hac data/pod5_output/all_reads.pod5 > data/basecalled_output/calls.bam

}

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

align_and_index() {
  # In the end it'd be good to have the option to align and basecall at once, or optionally do the two separately.
  local UNALIGNED_BAM="data/basecalled_output/taken_bam.bam"
  local REFERENCE_GENOME="reference_genomes/$1"
  local ALIGNED_BAM="data/alignment_output/aligned.sorted.bam"
  local THREADS=1
  local SORT_MEMORY_LIMIT="1G"

  echo "--- Starting alignment ---"
  mkdir -p "data/alignment_output"

  # Step 1: Alignment and sorting. It's a 3-stage pipe:
  # 1. samtools fastq converts the input BAM to FASTQ format
  # 2. minimap2 aligns the FASTQ reads to the reference
  # 3. samtools sort sorts the aligned output and saves it as a bam.bai file.

  echo "Aligning and sorting reads"

  minimap2 -ax map-ont -t "${THREADS}" -K5M "${REFERENCE_GENOME}" "${UNALIGNED_BAM}" \
  | samtools sort -@ "${THREADS}" -m "${SORT_MEMORY_LIMIT}" -o "${ALIGNED_BAM}" -

  echo "Indexing the sorted BAM file"
  # Index the sorted BAM file
  samtools index "${ALIGNED_BAM}"


}

run_script(){

  DORADO_VERSION="1.0.0"
  install_dorado "$DORADO_VERSION"
  #download_fast5_data 1
  #convert_fast5_to_pod5
  #basecalling_pod5 64
  download_reference_genome s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
  align_and_index "hg38.fa"

  alignment_qc

  ## MAKE A GENERIC DOWNLOADER ONE DAY
}

run_script