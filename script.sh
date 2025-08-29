#!/bin/bash

# --- Config ---

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

  echo "Dorado device flags set: ${DORADO_DEVICE_FLAG} ${DORAGO_BATCHSIZE_FLAG}"

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

  expected_folder=tools/dorado-${VERSION}-${os_type}
  echo "--- Checking for existing Dorado installation in ${expected_folder}---"
  if [ -d "${expected_folder}" ]; then
    echo "Dorado folder already exists"
    LOCAL_DORADO_BIN="${expected_folder}/bin"
    export PATH="$LOCAL_DORADO_BIN:$PATH"
    echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"
    return 0
  fi

  echo "Dorado folder ${expected_folder} not found, so I'll install it now."

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
  LOCAL_DORADO_BIN="${expected_folder}/bin"
  export PATH="$LOCAL_DORADO_BIN:$PATH"
  echo "The local dorado binary is in ${LOCAL_DORADO_BIN}"

  echo "Installing pod5 package"
  pip install pod5
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
  #aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request | head -n 20 | awk '{print $4}' | xargs -I {} aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/{} data/fast5_input/ --no-sign-request

  aws s3 ls s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/ --no-sign-request \
  | head -n $NUM_FILES \
  | awk '{print $4}' \
  | while read -r filename; do
      aws s3 cp s3://ont-open-data/rrms_2022.07/flowcells/Benchmarking_ASmethylation_COLO829_1-5/COLO829_1/20211102_1709_X1_FAR52193_a64b5c94/fast5_pass/"$filename" data/fast5_input/"$filename" --no-sign-request
  done
}

convert_fast5_to_pod5() {
  echo "Converting fast5 to pod5"

  mkdir -p data/pod5_output/

  pod5 convert fast5 data/fast5_input --output data/pod5_output/all_reads.pod5 --force-overwrite
}

download_dorado_model() {
  MODEL_ALIAS="hac"
  MODEL_DIR="models"

  mkdir -p "${MODEL_DIR}"
}

basecalling_pod5() {
  BATCHSIZE=$1
  echo "Adding dorado to the path"
  LOCAL_DORADO_BIN="$(pwd)/tools/dorado_v0.9.6/bin/"
  export PATH="$LOCAL_DORADO_BIN:$PATH"

  local POD5_INPUT="data/pod5_output/all_reads.pod5"
  local BASECALLED_OUTPUT="data/basecalled_output/calls.bam"

  check_device_for_dorado

  dorado basecaller --batchsize "${BATCHSIZE}" hac ${POD5_INPUT} > ${BASECALLED_OUTPUT}
}

run_script(){

  DORADO_VERSION="0.9.6"
  install_dorado "$DORADO_VERSION"
  download_fast5_data 10
  convert_fast5_to_pod5
  basecalling_pod5 128


}

run_script