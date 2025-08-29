#!/bin/bash

# --- Config ---

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


  S3_BASE_PATH="s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/"

  echo "--- Making fast5 input and basecalled output folders ---"
  mkdir -p data
  mkdir -p data/fast5_input
  mkdir -p data/basecalled_output

  echo "--- Downloading ${NUM_FILES} fast5 files ---"
  #aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request | head -n 5 | awk '{print $4}' | xargs -I {} aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/{} data/fast5_input/ --no-sign-request


  aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request \
  | head -n "$NUM_FILES" \
  | awk '{print $4}' \
  | while read -r filename; do
      aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/"$filename" data/fast5_input/"$filename" --no-sign-request || true
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
  MODEL_DIR="models"

  mkdir -p "${MODEL_DIR}"

  echo "--- Downloading model '${MODEL_ALIAS}' to local directory: ${MODEL_DIR} ---"

  # Use dorado download with --path flag. If the model already exists it'll just do nothing.
  # Additionally, it outputs the full path to the downloaded model
  FULL_MODEL_PATH=$(dorado download --model "${MODEL_ALIAS}" --data "${POD5_DIR}/all_reads.pod5" --models-directory "./${MODEL_DIR}/")

  echo "Model has been downloaded and put in ${FULL_MODEL_PATH}"

}

basecalling_pod5() {
  BATCHSIZE=$1

  # ---- You need to add the download_dorado_model here
  #  download_dorado_model

  # ---- You could also make it so that the alignment and basecalling is done at once. That should be the next thing to test.

  dorado basecaller --batchsize "${BATCHSIZE}" hac data/pod5_output/all_reads.pod5 > data/basecalled_output/calls.bam
}

download_reference_genome() {
  local url=$1
  local compressed_filename

  REF_DIR="reference_genomes"
  compressed_filename=$(basename "${url}")
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

  COMPRESSED_FILE="${REF_DIR}/${compressed_filename}"

  echo "And put it in $COMPRESSED_FILE"

  curl -L -f -o "${REF_DIR}/${compressed_filename}" "${url}"

  echo "Decompressing genome from ${COMPRESSED_FILE} and putting it into ${REF_FASTA}"

  gunzip -c "${COMPRESSED_FILE}" > "$REF_FASTA"

  echo "Cleaning up compressed file"
  rm "${COMPRESSED_FILE}"

}

run_alignment() {
  # In the end it'd be good to have the option to align and basecall at once, or optionally do the two separately.
  local UNALIGNED_BAM="data/basecalled_output/calls.bam"
  local REFERENCE_GENOME="reference_genome/$1"
  local ALIGNED_BAM="data/alignment_output/aligned.sorted.bam"
  local THREADS=8

  echo "--- Starting alignment ---"
  mkdir -p "data/alignment_output"

  # A 3-stage pipe
  samtools fastq -T


}

run_script(){
  #DORADO_VERSION="0.9.6"
  #install_dorado "$DORADO_VERSION"
  #download_fast5_data 1
  #convert_fast5_to_pod5
  #basecalling_pod5 64
  download_reference_genome "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz"
  run_alignment "hg38.fa"

  ## MAKE A GENERIC DOWNLOADER ONE DAY
}

run_script