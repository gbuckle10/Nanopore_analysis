#!/bin/bash

# Exit on error
set -e

# Source the logging utility
# Assumes that logging.sh is in the utils folder of the current folder.
source "$(dirname "$0")/utils/logging.sh"

# Read config values
CONFIG_FILE="config.yaml"
DORADO_VERSION=$(yq e '.parameters.setup.dorado_version' "${CONFIG_FILE}")
FAST5_DOWNLOAD_URL=$(yq e '.paths.fast5_download_url' "${CONFIG_FILE}")
FAST5_DESTINATION_DIR=$(yq e '.paths.fast5_input_dir' "${CONFIG_FILE}")
NUM_FAST5_FILES=$(yq e '.parameters.setup.num_fast5_files' "${CONFIG_FILE}")
POD5_DIR=$(yq e '.paths.pod5_dir' "${CONFIG_FILE}")
SHOULD_DOWNLOAD_FAST5_FILES=$(yq e '.pipeline_control.run_setup_tasks.download_fast5_data' "${CONFIG_FILE}")
SHOULD_CONVERT_FAST5_TO_POD5=$(yq e '.pipeline_control.run_setup_tasks.convert_fast5_to_pod5' "${CONFIG_FILE}")
RUNTIME_CONFIG="scripts/runtime_config.sh"

make_directories() {
  log_info "--- Creating directories ---"
  DIRECTORIES_TO_CREATE=$(yq '.paths.core_directories[]' "$CONFIG_FILE")

  for dir in ${DIRECTORIES_TO_CREATE}; do
    mkdir -p "${dir}"
  done

  log_info "Created top level directories"
}

install_dorado() {
  log_info "--- Downloading and setting up Dorado ---"

  case "$(uname -s)" in
    Linux*)
      os_type="linux-x64"
      archive_filename="dorado.tar.gz"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-${os_type}.tar.gz"
      ;;
    CYGWIN*|MINGW*|MSYS*)
      os_type="win64"
      archive_filename="dorado.zip"
      download_url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-${os_type}.zip"
      ;;
    *)
      log_error "I don't have a clue what operating system you're trying to use..."
      return 1
      ;;
  esac

  DORADO_DIR=tools/dorado-${DORADO_VERSION}-${os_type}
  log_info "--- Checking for existing Dorado installation in ${DORADO_DIR}---"
  if [ -d "${DORADO_DIR}" ]; then
    log_info "Dorado folder already exists, so I won't download it again."

  else
    log_info "Downloading dorado version ${DORADO_VERSION} for ${os_type} from ${download_url}" >&2
    curl -o tools/${archive_filename}  "$download_url"
  
    log_info "Download complete. Unzipping and setting up." >&2
  
    case ${os_type} in
      "linux-x64")
        tar -xzvf "tools/${archive_filename}" -C "tools/"
        ;;
      "win64")
        unzip -q "tools/${archive_filename}" -d "tools/"
        ;;
    esac

    folder_to_remove=tools/${archive_filename}
    log_info "Unzipping complete, removing ${folder_to_remove}"
    rm "${folder_to_remove}"
  fi
  
  LOCAL_DORADO_BIN="${DORADO_DIR}/bin"
  export PATH="$LOCAL_DORADO_BIN:$PATH"
  log_info "The local dorado binary is in ${LOCAL_DORADO_BIN}"

  DORADO_EXECUTABLE_PATH=$(realpath "${LOCAL_DORADO_BIN}/dorado")

  echo "export DORADO_EXECUTABLE=\"${DORADO_EXECUTABLE_PATH}\"" > "${RUNTIME_CONFIG}"
  log_info "Dorado was successfully installed. Dorado path was saved to ${RUNTIME_CONFIG}" >&2
}

download_fast5_data() {
  mkdir -p "${FAST5_DESTINATION_DIR}"

  log_info "--- Downloading ${NUM_FAST5_FILES} fast5 files into '${FAST5_DESTINATION_DIR}' ---"
  log_info "--- Downloading from '${FAST5_DOWNLOAD_URL}'"

  ALL_S3_FILENAMES=$(aws s3 ls "${FAST5_DOWNLOAD_URL}" --no-sign-request | grep '\.fast5' | awk '{print $4}')

  if [[ "${NUM_FAST5_FILES}" == "all" ]]; then
    log_info "Preparing to download all available files."
    filenames_to_download="${ALL_S3_FILES}"
  else
    log_info "INFO: Preparing to download the first ${NUM_FAST5_FILES} files."
    filenames_to_download=$(echo "${ALL_S3_FILENAMES}" | head -n "${NUM_FAST5_FILES}")
  fi

  log_info "Starting download..." >&2
  echo "${filenames_to_download}" | while read -r filename; do
    # skip any empty lines
    if [[ -z "$filename" ]]; then
      continue
    fi

    local destination="${FAST5_DESTINATION_DIR}"
    if [ ! -f "$destination" ]; then
      aws s3 cp "${FAST5_DOWNLOAD_URL}${filename}" "${destination}" --no-sign-request --quiet
    fi
  done


  #aws s3 ls "${FAST5_DOWNLOAD_URL}" --no-sign-request \
  #| head -n "${NUM_FAST5_FILES}" \
  #| awk '{print $4}' \
  #| while read -r filename; do
  #    aws s3 cp "${FAST5_DOWNLOAD_URL}${filename}" "${FAST5_DESTINATION_DIR}" --no-sign-request --quiet || true
  #done

  log_info "--- Fast5 input downloaded ---"

}

convert_fast5_to_pod5() {
  log_info "--- Converting fast5 to pod5 ---"

  mkdir -p "${POD5_DIR}"

  pod5 convert fast5 "${FAST5_DESTINATION_DIR}" --output "${POD5_DIR}/all_pod.pod5" --force-overwrite
}

download_methylation_atlas_and_illumina_manifest(){
  mkdir -p "data/atlas/"

  if ! [ -f "data/atlas/illumina_manifest.csv" ]; then
    log_info "The illumina manifest doesn't already exist, so we will download it."
    wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv -q -O "data/atlas/illumina_manifest.csv"
  else
    log_info "The illumina manifest is already there so won't be re-downloaded."
  fi
  if ! [ -f "data/atlas/full_atlas.csv" ]; then
    log_info "The methylation atlas doesn't already exist, so we will download it."
    wget https://github.com/nloyfer/meth_atlas/raw/refs/heads/master/full_atlas.csv.gz -q -O "data/atlas/full_atlas.csv.gz"
    gunzip -v "data/atlas/full_atlas.csv.gz"
    #rm "data/atlas/full_atlas.csv.gz"
  else
    log_info "The methylation atlas is already there so won't be re-downloaded."
  fi
  if ! [ -f "data/atlas/UXM_atlas.tsv" ]; then
    log_info "The UXM atlas doesn't already exist, so we'll download it."
    wget https://raw.githubusercontent.com/nloyfer/UXM_deconv/refs/heads/main/supplemental/Atlas.U25.l4.hg19.tsv -q -O "data/atlas/UXM_atlas.tsv"
  else
    log_info "The UXM atlas is already there so won't be re-downloaded."
  fi

}

make_directories
install_dorado
download_methylation_atlas_and_illumina_manifest
if [[ "$SHOULD_DOWNLOAD_FAST5_FILES" == "true" ]]; then
  download_fast5_data
fi
if [[ "$SHOULD_CONVERT_FAST5_TO_POD5" == "true" ]]; then
  convert_fast5_to_pod5
fi